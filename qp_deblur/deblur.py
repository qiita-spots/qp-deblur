# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os import mkdir, environ
from os.path import join, exists

from future.utils import viewitems
from functools import partial
from collections import OrderedDict
import json
from skbio import TreeNode
import pandas as pd

from biom import Table, load_table
from biom.util import biom_open
from qiita_client import ArtifactInfo
from qiita_client.util import system_call

from qiita_files.demux import to_per_sample_files
import qp_deblur


DEBLUR_PARAMS = {
    'Positive filtering database': 'pos-ref-fp',
    'Negative filtering database': 'neg-ref-fp',
    'Indexed positive filtering database': 'pos-ref-db-fp',
    'Indexed negative filtering database': 'neg-ref-db-fp',
    'Mean per nucleotide error rate': 'mean-error',
    'Error probabilities for each Hamming distance': 'error-dist',
    'Insertion/deletion (indel) probability': 'indel-prob',
    'Maximum number of insertion/deletion (indel)': 'indel-max',
    'Sequence trim length (-1 for no trimming)': 'trim-length',
    'Minimum dataset-wide read threshold': 'min-reads',
    'Minimum per-sample read threshold': 'min-size',
    'Threads per sample': 'threads-per-sample',
    'Jobs to start': 'jobs-to-start',
    'Reference phylogeny for SEPP': 'Greengenes_13.8'}


def generate_deblur_workflow_commands(preprocessed_fp, out_dir, parameters):
    """Generates the deblur commands

    Parameters
    ----------
    preprocessed_fp : list of str
        A list of one element with the input fastq or per-sample folder
    out_dir : str
        The job output directory
    parameters : dict
        The command's parameters, keyed by parameter name

    Returns
    -------
    list of str
        The deblur commands

    Raises
    ------
    ValueError
        If there is more than 1 file passed as preprocessed_fastq
    """
    if len(preprocessed_fp) != 1:
        raise ValueError("deblur doesn't accept more than one filepath: "
                         "%s" % ', '.join(preprocessed_fp))

    translated_params = {DEBLUR_PARAMS[k]: v
                         for k, v
                         in parameters.items()
                         if k != 'Reference phylogeny for SEPP'}
    params = OrderedDict(sorted(translated_params.items(), key=lambda t: t[0]))
    params = ['--%s "%s"' % (k, v) if v is not True else '--%s' % k
              for k, v in viewitems(params) if v != 'default']
    cmd = 'deblur workflow --seqs-fp "%s" --output-dir "%s" %s' % (
        preprocessed_fp[0], out_dir, ' '.join(params))

    return cmd


def generate_sepp_placements(seqs, out_dir, threads, reference_phylogeny=None,
                             reference_alignment=None):
    """Generates the SEPP commands

    Parameters
    ----------
    seqs : list of str
        A list of seqs to generate placements
    out_dir : str
        The job output directory
    threads : int
        Number if CPU cores to use
    reference_phylogeny : str, optional
        A filepath to an alternative reference phylogeny for SEPP.
        If None, default phylogeny (Greengenes 13.8 99% id) is used.
    reference_alignment : str, optional
        A filepath to an alternative reference alignment for SEPP.
        If None, default alignment (Greengenes 13.8 99% id) is used.

    Returns
    -------
    dict of strings
        keys are the seqs, values are the new placements as JSON strings

    Raises
    ------
    ValueError
        If run-sepp.sh does not produce expected file placements.json which is
        an indicator that something failed.
    """
    # return an empty dict if no sequences have been passed to the function
    if len(seqs) < 1:
        return {}

    # Create a multiple fasta file for all input seqs
    file_input = "%s/input.fasta" % out_dir
    with open(file_input, 'w') as fh_input:
        for seq in seqs:
            fh_input.write(">%s\n%s\n" % (seq, seq))

    # execute SEPP
    run_name = 'qiita'
    param_phylogeny = ''
    if reference_phylogeny is not None:
        param_phylogeny = ' -t %s ' % reference_phylogeny
    param_alignment = ''
    if reference_alignment is not None:
        param_alignment = ' -a %s ' % reference_alignment
    # SEPP writes output into the current working directory (cwd), therefore
    # we here first need to store the cwd, then move into the output directory,
    # perform SEPP and move back to the stored cwd for a clean state
    curr_pwd = environ['PWD']
    std_out, std_err, return_value = system_call(
        'cd %s && run-sepp.sh %s %s -x %s %s %s; cd %s' %
        (out_dir, file_input, run_name, threads,
         param_phylogeny, param_alignment, curr_pwd))

    # parse placements from SEPP results
    file_placements = '%s/%s_placement.json' % (out_dir, run_name)
    if exists(file_placements):
        with open(file_placements, 'r') as fh_placements:
            plcmnts = json.loads(fh_placements.read())
            return {p['nm'][0][0]: p['p'] for p in plcmnts['placements']}
    else:
        # due to the wrapper style of run-sepp.sh the actual exit code is never
        # returned and we have no way of finding out which sub-command failed
        # Therefore, we can only assume that something went wrong by not
        # observing the expected output file.
        # If the main SEPP program fails, it reports some information in two
        # files, the content of which we can read and report
        file_stderr = '%s/sepp-%s-err.log' % (out_dir, run_name)
        if exists(file_stderr):
            with open(file_stderr, 'r') as fh_stderr:
                std_err = fh_stderr.readlines()
        file_stdout = '%s/sepp-%s-out.log' % (out_dir, run_name)
        if exists(file_stdout):
            with open(file_stdout, 'r') as fh_stdout:
                std_out = fh_stdout.readlines()
        error_msg = ("Error running run-sepp.sh:\nStd out: %s\nStd err: %s"
                     % (std_out, std_err))
        raise ValueError(error_msg)


def _generate_template_rename(file_reference_phylogeny,
                              file_reference_alignment,
                              out_dir):
    """Produces placement template and rename script for reference phylogeny.

    Parameters
    ----------
    file_reference_phylogeny : str
        A filepath to an alternative reference phylogeny for SEPP.
    file_reference_alignment : str
        A filepath to an alternative reference alignment for SEPP.
    out_dir : str
        The job output directory

    Returns
    -------
    (str, str) : Filepaths of reference_template json file and
    reference_rename python script.

    Raises
    ------
    ValueError
        If a) the given out_dir directory does not exist.
        b) the given reference phylogeny or alignment does not exist.
        c) the run-sepp.sh wrapper script fails for any reason.

    Notes
    -----
    This function only needs to be called once per reference phylogeny/
    alignment, i.e. if we update Greengenes or extend SEPP for Silva or other
    reference phylogenies. I am including this function for easier maintenance
    in the future.
    """
    if not exists(out_dir):
        raise ValueError("Output directory '%s' does not exist!" % out_dir)
    if not exists(file_reference_phylogeny):
        raise ValueError("Reference phylogeny file '%s' does not exits!" %
                         file_reference_phylogeny)
    if not exists(file_reference_alignment):
        raise ValueError("Reference alignment file '%s' does not exits!" %
                         file_reference_alignment)

    # create a dummy sequence input file
    file_input = '%s/input.fasta' % out_dir
    with open(file_input, 'w') as f:
        f.write('>dummySeq\n')
        f.write('TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGATGGA'
                'CAAGTCTGATGTGAAAGGCTGGGGCCCAACCCCGGGACTGCATTGGAAACTGCCCGTCTT'
                'GAGTG\n')
    std_out, std_err, return_value = system_call(
        'cd %s; run-sepp.sh %s dummy -x 1 -a %s -t %s' %
        (out_dir, file_input, file_reference_alignment,
         file_reference_phylogeny))
    if return_value != 0:
        error_msg = ("Error running SEPP:\nStd out: %s\nStd err: %s"
                     % (std_out, std_err))
        raise ValueError(error_msg)

    # take resulting placement.json and turn it into the template by
    # clearing the list of placements
    file_template = '%s/tmpl_dummy_placement.json' % out_dir
    with open('%s/dummy_placement.json' % out_dir, 'r') as f:
        placements = json.loads(f.read())
        placements['placements'] = []
        with open(file_template, 'w') as fw:
            json.dump(placements, fw)

    # Another file produced by SEPP is xxx_rename-json.py, where xxx is the
    # name of the run, here "dummy". SEPP needs to escape node names before the
    # reference tree is given to guppy which can only handle a limited name
    # format. Thus, after guppy, the result needs to be back translated to
    # original names with the rename-json.py script that is generated by SEPP.
    return (file_template, '%s/dummy_rename-json.py' % out_dir)


def generate_insertion_trees(placements, out_dir,
                             reference_template=None,
                             reference_rename=None):
    """Generates phylogenetic trees by inserting placements into a reference

    Parameters
    ----------
    placements : dict of strings
        keys are the seqs, values are the new placements as JSON strings
    out_dir : str
        The job output directory
    reference_template : str, optional
        Filepath to the reference placement json file.
        This file can be produced via _generate_template_rename() and should be
        stored in the plugin package, because it can re used.
        If None, it falls back to the Greengenes 13.8 99% reference.
    reference_rename : str, optional
        Similar to reference_template, but a filepath to the generated python
        renaming script to undo the name scaping post guppy.
        If None, it falls back to the Greengenes 13.8 99% reference.

    Returns
    -------
    str
        The filepath of the phylogenetic insertion tree in Newick format.

    Raises
    ------
    ValueError
        If a) the given reference_template or reference_rename files do not
        exist
        b) or the guppy binary exits with non-zero return code
        c) or the given rename script exists with non-zero return code.
    """
    # test if reference file for rename script actually exists.
    file_ref_rename = qp_deblur.get_data(
        join('sepp', 'tmpl_gg13.8-99_rename-json.py'))
    if reference_rename is not None:
        file_ref_rename = reference_rename
    if not exists(file_ref_rename):
        raise ValueError("Reference rename script '%s' does not exits!" %
                         file_ref_rename)

    # create a valid placement.json file as input for guppy
    file_ref_template = qp_deblur.get_data(
        join('sepp', 'tmpl_gg13.8-99_placement.json'))
    if reference_template is not None:
        file_ref_template = reference_template
    if not exists(file_ref_template):
        raise ValueError("Reference template '%s' does not exits!" %
                         file_ref_template)
    with open(file_ref_template, 'r') as f:
        plcmnts = json.loads(f.read())

    plcmnts['placements'].extend(
        [{'p': placement, 'nm': [[sequence, 1]]}
         for sequence, placement
         in placements.items()])

    file_placements = '%s/placements.json' % out_dir
    with open(file_placements, 'w') as f:
        json.dump(plcmnts, f)

    # execute guppy
    file_tree_escaped = join(out_dir, 'insertion_tree.tre')
    std_out, std_err, return_value = system_call(
        'guppy tog %s -o %s' % (file_placements, file_tree_escaped))
    if return_value != 0:
        error_msg = ("Error running guppy:\nStd out: %s\nStd err: %s"
                     % (std_out, std_err))
        raise ValueError(error_msg)

    # execute node name re-labeling (to revert the escaping of names necessary
    # for guppy)
    file_tree = join(out_dir, 'insertion_tree.relabelled.tre')
    std_out, std_err, return_value = system_call(
        'cat %s | python %s > %s' %
        (file_tree_escaped, file_ref_rename, file_tree))
    if return_value != 0:
        error_msg = (("Error running %s:\n"
                      "Std out: %s\nStd err: %s")
                     % (file_ref_rename, std_out, std_err))
        raise ValueError(error_msg)

    # making sure that all branches in the generated tree have branch lenghts
    tree = TreeNode.read(file_tree)
    for node in tree.preorder(include_self=False):
        if node.length is None:
            node.length = 0.0
    tree.write(file_tree)

    return file_tree


def deblur(qclient, job_id, parameters, out_dir):
    """Run deblur with the given parameters

    Parameters
    ----------
    qclient : qiita_client.QiitaClient
        The Qiita server client
    job_id : str
        The job id
    parameters : dict
        The parameter values to run deblur
    out_dir : str
        The path to the job's output directory

    Returns
    -------
    boolean, list, str
        The results of the job

    Notes
    -----
    The code will check if the artifact has a preprocessed_demux element, if
    not it will use the preprocessed_fastq. We prefer to work with the
    preprocessed_demux as running time will be greatly improved
    """
    out_dir = join(out_dir, 'deblur_out')
    # Step 1 get the rest of the information need to run deblur
    qclient.update_job_step(job_id, "Step 1 of 4: Collecting information")
    artifact_id = parameters['Demultiplexed sequences']
    # removing input from parameters so it's not part of the final command
    del parameters['Demultiplexed sequences']

    # Get the artifact filepath information
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    fps = artifact_info['files']

    # Getting preparation information
    prep_info = qclient.get(
        '/qiita_db/prep_template/%s/' % artifact_info['prep_information'][0])
    df = pd.read_csv(prep_info['prep-file'], sep='\t')
    if 'platform' not in [x.lower() for x in df.columns]:
        error_msg = ('Preparation Information File does not have a platform '
                     'column, which is requiered')
        return False, None, error_msg
    if [x.lower() for x in df.platform.unique()] != ['illumina']:
        error_msg = ('deblur is only valid for Illumina `platform`, current '
                     'values in the Preparation Information File: %s' %
                     ', '.join(sorted(df.platform.unique())))
        return False, None, error_msg

    # Step 2 generating command deblur
    if 'preprocessed_demux' in fps:
        qclient.update_job_step(job_id, "Step 2 of 4: Generating per sample "
                                "from demux (1/2)")

        if not exists(out_dir):
            mkdir(out_dir)
        split_out_dir = join(out_dir, 'split')
        if not exists(split_out_dir):
            mkdir(split_out_dir)

        # using the same number of parallel jobs as defined by the command
        n_jobs = int(parameters['Jobs to start'])
        # [0] cause there should be only 1 file
        to_per_sample_files(fps['preprocessed_demux'][0],
                            out_dir=split_out_dir, n_jobs=n_jobs)

        qclient.update_job_step(job_id, "Step 2 of 4: Generating per sample "
                                "from demux (2/2)")
        out_dir = join(out_dir, 'deblured')
        cmd = generate_deblur_workflow_commands([split_out_dir],
                                                out_dir, parameters)
    else:
        qclient.update_job_step(job_id, "Step 2 of 4: Generating deblur "
                                "command")
        cmd = generate_deblur_workflow_commands(fps['preprocessed_fastq'],
                                                out_dir, parameters)

    # Step 3 execute deblur
    qclient.update_job_step(job_id, "Step 3 of 4: Executing deblur job")
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error running deblur:\nStd out: %s\nStd err: %s"
                     % (std_out, std_err))
        return False, None, error_msg

    # Generating artifact
    pb = partial(join, out_dir)

    # Generate the filepaths
    final_biom = pb('all.biom')
    final_seqs = pb('all.seqs.fa')
    final_biom_hit = pb('reference-hit.biom')
    final_seqs_hit = pb('reference-hit.seqs.fa')

    if not exists(final_biom_hit):
        # Create an empty table. We need to send something to Qiita that is
        # a valid BIOM, so we are going to create an empty table
        t = Table([], [], [])
        with biom_open(final_biom_hit, 'w') as f:
            t.to_hdf5(f, 'qp-deblur generated')

    if not exists(final_seqs_hit):
        # Same as before, create an empty sequence file so we can send it
        with open(final_seqs_hit, 'w') as f:
            f.write("")

    # Step 4, communicate with archive to check and generate placements
    qclient.update_job_step(job_id, "Step 4 of 4 (1/4): Retrieving "
                            "observations information")
    features = list(load_table(final_biom_hit).ids(axis='observation'))

    fp_phylogeny = None
    if features:
        observations = qclient.post(
            "/qiita_db/archive/observations/", data={'job_id': job_id,
                                                     'features': features})
        novel_fragments = list(set(features) - set(observations.keys()))

        qclient.update_job_step(job_id, "Step 4 of 4 (2/4): Generating %d new "
                                "placements" % len(novel_fragments))

        # Once we support alternative reference phylogenies for SEPP in the
        # future, we need to translate the reference name here into
        # filepaths pointing to the correct reference alignment and
        # reference tree. If left 'None' the Greengenes 13.8 reference
        # shipped with the fragment-insertion conda package will be used.
        fp_reference_alignment = None
        fp_reference_phylogeny = None
        fp_reference_template = None
        fp_reference_rename = None
        if 'Reference phylogeny for SEPP' in parameters:
            if parameters['Reference phylogeny for SEPP'] == 'tiny':
                fp_reference_alignment = qp_deblur.get_data(join(
                    'sepp', 'reference_alignment_tiny.fasta'))
                fp_reference_phylogeny = qp_deblur.get_data(join(
                    'sepp', 'reference_phylogeny_tiny.nwk'))
                fp_reference_template = qp_deblur.get_data(join(
                    'sepp', 'tmpl_tiny_placement.json'))
                fp_reference_rename = qp_deblur.get_data(join(
                    'sepp', 'tmpl_tiny_rename-json.py'))
        try:
            new_placements = generate_sepp_placements(
                novel_fragments, out_dir, parameters['Threads per sample'],
                reference_alignment=fp_reference_alignment,
                reference_phylogeny=fp_reference_phylogeny)
        except ValueError as e:
            return False, None, str(e)

        qclient.update_job_step(job_id, "Step 4 of 4 (3/4): Archiving %d "
                                "new placements" % len(novel_fragments))
        # values needs to be json strings as well
        for fragment in new_placements.keys():
            new_placements[fragment] = json.dumps(new_placements[fragment])

        # fragments that get rejected by a SEPP run don't show up in
        # the placement file, however being rejected is a valuable
        # information and should be stored in the archive as well.
        # Thus, we avoid re-computation for rejected fragments in the
        # future.
        for fragment in novel_fragments:
            if fragment not in new_placements:
                new_placements[fragment] = ""
        if len(new_placements.keys()) > 0:
            qclient.patch(url="/qiita_db/archive/observations/", op="add",
                          path=job_id, value=json.dumps(new_placements))

        # retrieve all fragments and create actuall tree
        qclient.update_job_step(job_id, "Step 4 of 4 (4/4): Composing "
                                "phylogenetic insertion tree")
        placements = qclient.post(
            "/qiita_db/archive/observations/", data={'job_id': job_id,
                                                     'features': features})
        # remove fragments that have been rejected by SEPP, i.e. whoes
        # placement is the empty string and
        # convert all other placements from string to json
        placements = {frag: json.loads(placements[frag])
                      for frag, plc
                      in placements.items()
                      if plc != ''}
        try:
            fp_phylogeny = generate_insertion_trees(
                placements, out_dir,
                reference_template=fp_reference_template,
                reference_rename=fp_reference_rename)
        except ValueError as e:
            return False, None, str(e)
    else:
        new_placements = None

    ainfo = [ArtifactInfo('deblur final table', 'BIOM',
                          [(final_biom, 'biom'),
                           (final_seqs, 'preprocessed_fasta')])]
    if fp_phylogeny is not None:
        ainfo.append(ArtifactInfo('deblur reference hit table', 'BIOM',
                     [(final_biom_hit, 'biom'),
                      (final_seqs_hit, 'preprocessed_fasta'),
                      (fp_phylogeny, 'plain_text')], new_placements))

    return True, ainfo, ""


def generate_tree_from_fragments(fp_placements,
                                 fp_biom,
                                 out_dir,
                                 fp_reference_template=None,
                                 fp_reference_rename=None):
    """Generates a phylogenetic tree by inserting placements into a reference,
       and trims observations in BIOMs to those successfully matched to the
       tree.

    Parameters
    ----------
    fp_placements : str
        The path to a file containing a JSON dump of fragments
        (json dump is of dictionary, where keys are the seqs,
         values are the new placements as JSON strings)
    fp_biom : str
        The path to a BIOM file.
    out_dir : str
        The path to the job's output directory
    fp_reference_template : str, optional
        The path to a file containing the reference placement file.
        If None, it falls back to the Greengenes 13.8 99% reference.
    fp_reference_rename : str, optional
        Similar to fp_reference_template, but a filepath to the generated
        python renaming script to undo the name scaping post guppy.
        If None, it falls back to the Greengenes 13.8 99% reference.

    Returns
    -------
        The filepath of the phylogenetic insertion tree in Newick
        The filepath of the final, trimmed biom file.

    Raises
    ------
    ValueError
        If the given reference_template or reference_rename files do not exist
        If the guppy binary exits with non-zero return code
        If the given rename script exists with non-zero return code.
    """
    with open(fp_placements) as placements_file:
        placements = json.load(placements_file)

        try:
            fp_phylogeny = generate_insertion_trees(
                                    placements,
                                    out_dir,
                                    reference_template=fp_reference_template,
                                    reference_rename=fp_reference_rename)
        except Exception:
            # we can get an exception if the tree can't be build; there are
            # many reasons for this but perhaps the most important is a
            # different target region.
            fp_phylogeny = None
            fp_biom_out = None

        if fp_biom is not None and fp_phylogeny is not None:
            # read tree
            tree = TreeNode.read(str(fp_phylogeny))
            fragments_tree = {str(tip.name) for tip in tree.tips()
                              if tip.name is not None}

            # load biom w/fragments not found in fp_phylogeny
            biom_table = load_table(fp_biom)
            fragments_table = set(map(str, biom_table.ids(axis='observation')))

            # filter biom file
            tbl_matched = biom_table.filter(fragments_table & fragments_tree,
                                            axis='observation',
                                            inplace=False)

            fp_biom_out = '%s_insertion_filter.biom' % fp_biom[:-len('.biom')]
            with biom_open(fp_biom_out, 'w') as f:
                tbl_matched.to_hdf5(
                    f, "Generated by Qiita, qp-deblur")
        else:
            fp_biom_out = None

        return fp_phylogeny, fp_biom_out
