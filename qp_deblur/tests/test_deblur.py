# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os import close, remove
from shutil import copyfile, rmtree
from tempfile import mkstemp, mkdtemp
from json import dumps
from os.path import exists, isdir, join

from qiita_client.testing import PluginTestCase


from qp_deblur import plugin
from qp_deblur.deblur import (
    deblur, generate_deblur_workflow_commands)


class deblurTests(PluginTestCase):
    def setUp(self):
        plugin("https://localhost:21174", 'register', 'ignored')
        self.params = {
            'ref-fp': '', 'ref-db-fp': '', 'mean-error': 0.005,
            'error-dist': ('1, 0.06, 0.02, 0.02, 0.01, 0.005, 0.005, '
                           '0.005, 0.001, 0.001, 0.001, 0.0005'),
            'indel-prob': 0.01, 'indel-max': 3, 'trim-length': 100,
            'min-reads': 0, 'min-size': 2, 'negate': True,
            'threads-per-sample': 1, 'jobs-to-start': 1, 'overwrite': True}
        self._clean_up_files = []

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_generate_deblur_workflow_commands_error(self):
        with self.assertRaises(ValueError):
            generate_deblur_workflow_commands(
                ['fastq/s1.fastq', 'fastq/s1.fastq'],
                'output', self.params)

    def test_generate_deblur_workflow_commands(self):
        exp = ('deblur workflow --seqs-fp "fastq/s1.fastq" '
               '--output-dir "output" --error-dist "1, 0.06, 0.02, 0.02, '
               '0.01, 0.005, 0.005, 0.005, 0.001, 0.001, 0.001, 0.0005" '
               '--indel-max "3" --indel-prob "0.01" --jobs-to-start "1" '
               '--mean-error "0.005" --min-size "2" --negate --overwrite '
               '--threads-per-sample "1" --trim-length "100"')
        obs = generate_deblur_workflow_commands(
            ['fastq/s1.fastq'], 'output', self.params)

        self.assertEqual(obs, exp)

    def test_deblur(self):
        # generating filepaths
        fd, fp = mkstemp(suffix='_seqs.fastq')
        close(fd)
        self._clean_up_files.append(fp)
        copyfile('support_files/filtered_5_seqs.fastq', fp)

        # inserting new prep template
        prep_info_dict = {
            'SKB7.640196': {'description': 'SKB7'},
            'SKB8.640193': {'description': 'SKB8'}
        }
        data = {'prep_info': dumps(prep_info_dict),
                # magic #1 = testing study
                'study': 1,
                'data_type': '16S'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        # inserting artifacts
        data = {
            'filepaths': dumps([(fp, 'preprocessed_fastq')]),
            'type': "Demultiplexed",
            'name': "New demultiplexed artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['seqs-fp'] = aid

        data = {'user': 'demo@microbio.me',
                'command': dumps(['deblur', '0.1.0', 'deblur-workflow']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        success, ainfo, msg = deblur(self.qclient, jid, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        # there is only one artifact, thus [0]
        ainfo = ainfo[0]
        self.assertEqual("BIOM", ainfo.artifact_type)
        exp = [(join(out_dir, 'final.biom'), 'biom'),
               (join(out_dir, 'final.seqs.fa'), 'preprocessed_fasta')]
        self.assertEqual(exp, ainfo.files)


if __name__ == '__main__':
    main()
