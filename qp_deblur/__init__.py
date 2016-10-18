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
    'deblur', '0.1.0', 'A greedy deconvolution algorithm based on Illumina '
    'Miseq/Hiseq error profiles')

# Define the deblur-workflow command
req_params = {'seqs-fp': ('artifact', ['Demultiplexed'])}
opt_params = {
    # parameters not being passed
    # output-dir
    # keep-tmp-files
    # log-level
    # log-file
    # overwrite
    # is-worker-thread
    'pos-ref-fp': ['choice:["default"]', 'default'],
    'neg-ref-fp': ['choice:["default"]', 'default'],
    'pos-ref-db-fp': ['choice:["default"]', 'default'],
    'neg-ref-db-fp': ['choice:["default"]', 'default'],
    'mean-error': ['float', '0.005'],
    'error-dist': ['string', ('1, 0.06, 0.02, 0.02, 0.01, 0.005, 0.005, '
                              '0.005, 0.001, 0.001, 0.001, 0.0005')],
    'indel-prob': ['float', '0.01'],
    'indel-max': ['integer', '3'],
    'trim-length': ['integer', '100'],
    'min-reads': ['integer', '0'],
    'min-size': ['integer', '2'],
    'negate': ['boolean', 'True'],
    'threads-per-sample': ['integer', '1'],
    'jobs-to-start': ['integer', '1']
}
outputs = {'deblur table': 'BIOM', 'deblur seqs': 'FASTA'}
dflt_param_set = {
    'Defaults': {'pos-ref-fp': 'default', 'neg-ref-fp': 'default',
                 'pos-ref-db-fp': 'default', 'neg-ref-db-fp': 'default',
                 'mean-error': 0.005,
                 'error-dist': ('1, 0.06, 0.02, 0.02, 0.01, 0.005, 0.005, '
                                '0.005, 0.001, 0.001, 0.001, 0.0005'),
                 'indel-prob': 0.01, 'indel-max': 3, 'trim-length': 100,
                 'min-reads': 10, 'min-size': 2, 'negate': True,
                 'threads-per-sample': 1, 'jobs-to-start': 1}
}
deblur_cmd = QiitaCommand(
    "deblur-workflow", "deblurring workflow", deblur, req_params, opt_params,
    outputs, dflt_param_set)
plugin.register_command(deblur_cmd)
