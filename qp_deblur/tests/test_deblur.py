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
        # this will allow us to see the full errors
        self.maxDiff = None

        plugin("https://localhost:21174", 'register', 'ignored')
        self.params = {
            'pos-ref-fp': 'default', 'neg-ref-fp': 'default',
            'mean-error': 0.005,
            'error-dist': ('1, 0.06, 0.02, 0.02, 0.01, 0.005, 0.005, '
                           '0.005, 0.001, 0.001, 0.001, 0.0005'),
            'indel-prob': 0.01, 'indel-max': 3, 'trim-length': 100,
            'min-reads': 0, 'min-size': 2, 'negate': True,
            'threads-per-sample': 1, 'jobs-to-start': 1, 'skip-trimming': True}
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
        exp = ('deblur workflow --seqs-fp "fastq/s1.fastq" --output-dir '
               '"output" --error-dist "1, 0.06, 0.02, 0.02, 0.01, '
               '0.005, 0.005, 0.005, 0.001, 0.001, 0.001, 0.0005" '
               '--indel-max "3" --indel-prob "0.01" --jobs-to-start "1" '
               '--mean-error "0.005" --min-reads "0" --min-size "2" '
               '--negate --skip-trimming --threads-per-sample "1" '
               '--trim-length "100"')
        obs = generate_deblur_workflow_commands(
            ['fastq/s1.fastq'], 'output', self.params)

        self.assertEqual(obs, exp)

    def test_deblur_empty(self):
        # generating filepaths
        fd, fp = mkstemp(suffix='_seqs.fastq')
        close(fd)
        self._clean_up_files.append(fp)
        with open('support_files/filtered_5_seqs.fastq', 'r') as f:
            with open(fp, 'w') as outf:
                for _ in range(4):
                    outf.write(f.readline())
                    outf.write(f.readline())
                    outf.write(f.readline())
                    outf.write(f.readline())

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

        self.assertEqual(msg, "")
        self.assertTrue(success)

        self.assertEqual(ainfo[0].artifact_type, "BIOM")
        self.assertEqual(ainfo[1].artifact_type, "FASTA")
        self.assertEqual(ainfo[2].artifact_type, "BIOM")
        self.assertEqual(ainfo[3].artifact_type, "FASTA")

        exp_final_biom = join(out_dir, 'deblur_out', 'final.biom')
        exp_final_seqs = join(out_dir, 'deblur_out', 'final.seqs.fa')
        exp_final_biom_16s = join(out_dir, 'deblur_out', 'final.only-16s.biom')
        exp_final_seqs_na = join(out_dir, 'deblur_out',
                                 'final.seqs.fa.no_artifacts')

        self.assertEqual(ainfo[0].files, [(exp_final_biom, 'biom')])
        self.assertEqual(ainfo[1].files,
                         [(exp_final_seqs, 'preprocessed_fasta')])
        self.assertEqual(ainfo[2].files, [(exp_final_biom_16s, 'biom')])
        self.assertEqual(ainfo[3].files,
                         [(exp_final_seqs_na, 'preprocessed_fasta')])

        self.assertTrue(exists(exp_final_biom))
        self.assertTrue(exists(exp_final_seqs))
        self.assertTrue(exists(exp_final_biom_16s))
        self.assertTrue(exists(exp_final_seqs_na))


    def test_deblur(self):
        # generating filepaths
        fd, fp = mkstemp(suffix='_seqs.demux')
        close(fd)
        self._clean_up_files.append(fp)
        copyfile('support_files/filtered_5_seqs.demux', fp)

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

        self.assertEqual(msg, "")
        self.assertTrue(success)

        self.assertEqual(ainfo[0].artifact_type, "BIOM")
        self.assertEqual(ainfo[1].artifact_type, "FASTA")
        self.assertEqual(ainfo[2].artifact_type, "BIOM")
        self.assertEqual(ainfo[3].artifact_type, "FASTA")

        exp_final_biom = join(out_dir, 'deblur_out', 'final.biom')
        exp_final_seqs = join(out_dir, 'deblur_out', 'final.seqs.fa')
        exp_final_biom_16s = join(out_dir, 'deblur_out', 'final.only-16s.biom')
        exp_final_seqs_na = join(out_dir, 'deblur_out',
                                 'final.seqs.fa.no_artifacts')

        self.assertEqual(ainfo[0].files, [(exp_final_biom, 'biom')])
        self.assertEqual(ainfo[1].files,
                         [(exp_final_seqs, 'preprocessed_fasta')])
        self.assertEqual(ainfo[2].files, [(exp_final_biom_16s, 'biom')])
        self.assertEqual(ainfo[3].files,
                         [(exp_final_seqs_na, 'preprocessed_fasta')])

        self.assertTrue(exists(exp_final_biom))
        self.assertTrue(exists(exp_final_seqs))
        self.assertTrue(exists(exp_final_biom_16s))
        self.assertTrue(exists(exp_final_seqs_na))

    def test_deblur_demux(self):
        # generating filepaths
        fd, fp = mkstemp(suffix='_seqs.demux')
        close(fd)
        self._clean_up_files.append(fp)
        copyfile('support_files/filtered_5_seqs.demux', fp)

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
            'filepaths': dumps([(fp, 'preprocessed_demux')]),
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

        self.assertEqual("BIOM", ainfo[0].artifact_type)
        self.assertEqual("FASTA", ainfo[1].artifact_type)
        self.assertEqual("BIOM", ainfo[2].artifact_type)
        self.assertEqual("FASTA", ainfo[3].artifact_type)

        self.assertEqual([(join(out_dir, 'deblur_out', 'deblured',
                                'final.biom'),
                           'biom')], ainfo[0].files)
        self.assertEqual([(join(out_dir, 'deblur_out', 'deblured',
                                'final.seqs.fa'),
                           'preprocessed_fasta')], ainfo[1].files)
        self.assertEqual([(join(out_dir, 'deblur_out', 'deblured',
                                'final.only-16s.biom'),
                          'biom')], ainfo[2].files)
        self.assertEqual([(join(out_dir, 'deblur_out', 'deblured',
                                'final.seqs.fa.no_artifacts'),
                          'preprocessed_fasta')], ainfo[3].files)

if __name__ == '__main__':
    main()
