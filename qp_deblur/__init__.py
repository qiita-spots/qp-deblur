# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client import QiitaPlugin, QiitaCommand

from .deblur import deblur

__all__ = ['deblur']


# Initialize the plugin
plugin = QiitaPlugin(
    'deblur', '1.0.3', 'A greedy deconvolution algorithm based on Illumina '
    'Miseq/Hiseq error profiles')

# Define the deblur-workflow command
req_params = {'Demultiplexed sequences': ('artifact', ['Demultiplexed'])}
opt_params = {
    # parameters not being passed
    # output-dir
    # keep-tmp-files
    # log-level
    # log-file
    # overwrite
    # is-worker-thread
    'Positive filtering database': ['choice:["default"]', 'default'],
    'Negative filtering database': ['choice:["default"]', 'default'],
    'Indexed positive filtering database': ['choice:["default"]', 'default'],
    'Indexed negative filtering database': ['choice:["default"]', 'default'],
    'Mean per nucleotide error rate': ['float', '0.005'],
    'Error probabilities for each Hamming distance': [
        'string', ('1, 0.06, 0.02, 0.02, 0.01, 0.005, 0.005, '
                   '0.005, 0.001, 0.001, 0.001, 0.0005')],
    'Insertion/deletion (indel) probability': ['float', '0.01'],
    'Maximum number of insertion/deletion (indel)': ['integer', '3'],
    'Sequence trim length (-1 for no trimming)': ['integer', '-1'],
    'Minimum dataset-wide read threshold': ['integer', '0'],
    'Minimum per-sample read threshold': ['integer', '2'],
    'Threads per sample': ['integer', '1'],
    'Jobs to start': ['integer', '1'],
    'Reference phylogeny for SEPP': ['choice:["Greengenes_13.8"]',
                                     'Greengenes_13.8']
}
outputs = {'deblur final table': 'BIOM',
           'deblur reference hit table': 'BIOM'}
dflt_param_set = {
    'Defaults': {'Positive filtering database': 'default',
                 'Negative filtering database': 'default',
                 'Indexed positive filtering database': 'default',
                 'Indexed negative filtering database': 'default',
                 'Mean per nucleotide error rate': 0.005,
                 'Error probabilities for each Hamming distance': (
                    '1, 0.06, 0.02, 0.02, 0.01, 0.005, 0.005, '
                    '0.005, 0.001, 0.001, 0.001, 0.0005'),
                 'Insertion/deletion (indel) probability': 0.01,
                 'Maximum number of insertion/deletion (indel)': 3,
                 'Sequence trim length (-1 for no trimming)': -1,
                 # the default min-reads is 10, however to ensure that deblur
                 # is actually per-sample, min-reads must be set to 0 otherwise
                 # filtering is applied over the samples included in a single
                 # run
                 'Minimum dataset-wide read threshold': 0,
                 'Minimum per-sample read threshold': 2,
                 'Threads per sample': 1, 'Jobs to start': 1,
                 'Reference phylogeny for SEPP': 'Greengenes_13.8'}
}
deblur_cmd = QiitaCommand(
    "Deblur", "deblurring workflow", deblur, req_params, opt_params,
    outputs, dflt_param_set)
plugin.register_command(deblur_cmd)
