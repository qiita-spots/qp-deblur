# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from os.path import join

from future.utils import viewitems
from functools import partial
from collections import OrderedDict

from qiita_client import ArtifactInfo
from qiita_client.util import system_call


def generate_deblur_workflow_commands(preprocessed_fastq, out_dir, parameters,
                                      assure_od=True):
    """Generates the deblur commands

    Parameters
    ----------
    preprocessed_fastq : str
        The preprocessed_fastq filepaths
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
    if len(preprocessed_fastq) != 1:
        raise ValueError("deblur doesn't accept more than one FASTQ: "
                         "%s" % ', '.join(preprocessed_fastq))

    params = OrderedDict(sorted(parameters.items(), key=lambda t: t[0]))
    params = ['--%s "%s"' % (k, v) if v is not True else '--%s' % k
              for k, v in viewitems(params) if v]
    cmd = 'deblur workflow --seqs-fp "%s" --output-dir "%s" %s' % (
        preprocessed_fastq[0], out_dir, ' '.join(params))

    return cmd


def _run_commands(qclient, job_id, commands, msg):
    for i, cmd in enumerate(commands):
        qclient.update_job_step(job_id, msg % i)
        std_out, std_err, return_value = system_call(cmd)
        if return_value != 0:
            error_msg = ("Error running deblur:\nStd out: %s\nStd err: %s"
                         % (std_out, std_err))
            return False, error_msg

    return True, ""


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
    """
    # Step 1 get the rest of the information need to run deblur
    qclient.update_job_step(job_id, "Step 1 of 3: Collecting information")
    artifact_id = parameters['seqs-fp']
    # removing input from parameters so it's not part of the final command
    del parameters['seqs-fp']

    # Get the artifact filepath information
    artifact_info = qclient.get("/qiita_db/artifacts/%s/" % artifact_id)
    fps = artifact_info['files']

    # Step 2 generating command deblur
    qclient.update_job_step(job_id, "Step 2 of 3: Generating deblur command")
    cmd = generate_deblur_workflow_commands(fps['preprocessed_fastq'],
                                            out_dir, parameters)

    # Step 3 execute deblur
    msg = "Step 3 of 3: Executing deblur job (%s/1)"
    success, msg = _run_commands(qclient, job_id, [cmd], msg)
    if not success:
        return False, None, msg

    # Generating artifact
    pb = partial(join, out_dir)
    filepaths = [(pb('final.biom'), 'biom'),
                 (pb('final.seqs.fa'), 'preprocessed_fasta')]
    ainfo = [ArtifactInfo('OTU table', 'BIOM', filepaths)]

    return True, ainfo, ""
