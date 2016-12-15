# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os import mkdir
from os.path import join, exists

from future.utils import viewitems
from functools import partial
from collections import OrderedDict

from biom import Table
from biom.util import biom_open
from qiita_client import ArtifactInfo
from qiita_client.util import system_call

from qiita_files.demux import to_per_sample_files


def generate_deblur_workflow_commands(preprocessed_fp, out_dir, parameters):
    """Generates the deblur commands

    Parameters
    ----------
    preprocessed_fp : list of str
        A list of 1 element with the input fastq or per sample folder
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

    params = OrderedDict(sorted(parameters.items(), key=lambda t: t[0]))
    params = ['--%s "%s"' % (k, v) if v is not True else '--%s' % k
              for k, v in viewitems(params) if v != 'default']
    cmd = 'deblur workflow --seqs-fp "%s" --output-dir "%s" %s' % (
        preprocessed_fp[0], out_dir, ' '.join(params))

    return cmd


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
    qclient.update_job_step(job_id, "Step 1 of 3: Collecting information")
    artifact_id = parameters['seqs-fp']
    # removing input from parameters so it's not part of the final command
    del parameters['seqs-fp']

    # Get the artifact filepath information
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    fps = artifact_info['files']

    # Step 2 generating command deblur
    if 'preprocessed_demux' in fps:
        qclient.update_job_step(job_id, "Step 2 of 3: Generating per sample "
                                "from demux (1/2)")

        if not exists(out_dir):
            mkdir(out_dir)
        split_out_dir = join(out_dir, 'split')
        if not exists(split_out_dir):
            mkdir(split_out_dir)

        # using the same number of parallel jobs as defined by the command
        n_jobs = parameters['jobs-to-start']
        # [0] cause there should be only 1 file
        to_per_sample_files(fps['preprocessed_demux'][0],
                            out_dir=split_out_dir, n_jobs=n_jobs)

        qclient.update_job_step(job_id, "Step 2 of 3: Generating per sample "
                                "from demux (2/2)")
        out_dir = join(out_dir, 'deblured')
        cmd = generate_deblur_workflow_commands([split_out_dir],
                                                out_dir, parameters)
    else:
        qclient.update_job_step(job_id, "Step 2 of 3: Generating deblur "
                                "command")
        cmd = generate_deblur_workflow_commands(fps['preprocessed_fastq'],
                                                out_dir, parameters)

    # Step 3 execute deblur
    qclient.update_job_step(job_id, "Step 3 of 3: Executing deblur job")
    std_out, std_err, return_value = system_call(cmd)
    if return_value != 0:
        error_msg = ("Error running deblur:\nStd out: %s\nStd err: %s"
                     % (std_out, std_err))
        return False, None, error_msg

    # Generating artifact
    pb = partial(join, out_dir)

    # Generate the filepaths
    final_biom = pb('final.biom')
    final_seqs = pb('final.seqs.fa')
    final_biom_16s = pb('final.only-16s.biom')
    final_seqs_na = pb('final.seqs.fa.no_artifacts')

    if not exists(final_biom_16s):
        # Create an empty table. We need to send something to Qiita that is
        # a valid BIOM, so we are going to create an empty table
        t = Table([], [], [])
        with biom_open(final_biom_16s, 'w') as f:
            t.to_hdf5(f, 'qp-deblur generated')

    if not exists(final_seqs_na):
        # Same as before, create an empty sequence file so we can send it
        with open(final_seqs_na, 'w') as f:
            f.write("")

    ainfo = [ArtifactInfo('deblur final table', 'BIOM',
                          [(final_biom, 'biom')]),
             ArtifactInfo('deblur final seqs', 'FASTA',
                          [(final_seqs, 'raw_fasta')]),
             ArtifactInfo('deblur 16S only table', 'BIOM',
                          [(final_biom_16s, 'biom')]),
             ArtifactInfo('deblur 16S only seqs', 'FASTA',
                          [(final_seqs_na, 'raw_fasta')])]

    return True, ainfo, ""
