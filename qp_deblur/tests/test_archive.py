# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

import json
from unittest import main
from os import close, remove
from shutil import copyfile, rmtree
from tempfile import mkstemp, mkdtemp
from json import dumps
from os.path import exists, isdir

from qiita_client.testing import PluginTestCase


from qp_deblur import plugin
from qp_deblur.deblur import deblur


class deblurTests(PluginTestCase):
    def setUp(self):
        # this will allow us to see the full errors
        self.maxDiff = None

        plugin("https://localhost:8383", 'register', 'ignored')
        self.params = {
            'Positive filtering database': 'default',
            'Negative filtering database': 'default',
            'Mean per nucleotide error rate': 0.005,
            'Error probabilities for each Hamming distance': (
                '1, 0.06, 0.02, 0.02, 0.01, 0.005, 0.005, '
                '0.005, 0.001, 0.001, 0.001, 0.0005'),
            'Insertion/deletion (indel) probability': 0.01,
            'Maximum number of insertion/deletion (indel)': 3,
            'Sequence trim length (-1 for no trimming)': 100,
            'Minimum dataset-wide read threshold': 0,
            'Minimum per-sample read threshold': 2,
            'Threads per sample': 1, 'Jobs to start': 1,
            'Reference phylogeny for SEPP': 'Greengenes_13.8'}
        self._clean_up_files = []
        self.novel_seqs = [
            ('TACGAAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCGTTAAG'
             'TTGGATGTGAAAGCCCCGGGCTCAACCTGGGAACTG'),
            ('TACGTAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTCGGAAAG'
             'AAAGATGTGAAATCCCAGGGCTTAACCTTGGAACTG'),
            ('TACGAAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTCAGCAAG'
             'TTGAAGGTGAAATCCCCGGGCTCAACCTGGGAACTG'),
            ('TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGTTAAG'
             'TCAGGGGTGAAATACGGTGGCTCAACCATCGCAGTG'),
            ('TACGAAGGGTGCAAGCGTTACTCGGAATTACTGGGCGTAAAGCGTGCGTAGGTGGTTGTTTAAG'
             'TCTGTTGTGAAAGCCCTGGGCTCAACCTGGGAATTG'),
            ('TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTTTGTAAG'
             'TCAGATGTGAAAGCCCCGGGCTCAACCTGGGAACTG'),
            ('TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAG'
             'TCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTG'),
            ('TACGAAGGGGGCTAGCGTTGTTCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGATATTTAAG'
             'TCAGGGGTGAAATCCCAGAGCTCAACTCTGGAACTG'),
            ('TACGTAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTTGTGCAAG'
             'ACAGATGTGAAATCCCCGGGCTCAACCTGGGAATTG'),
            ('TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGTAGGCGGTTCGTTAAG'
             'CCAGCTGTGAAATCCCCGGGCTCAACCTGGGAATTG')
        ]

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_fragment_archiving(self):
        # generating filepaths
        fd, fp = mkstemp(suffix='_seqs.demux')
        close(fd)
        self._clean_up_files.append(fp)
        copyfile('support_files/filtered_5_seqs.demux', fp)

        # inserting new prep template
        prep_info_dict = {
            'SKB7.640196': {'description_prep': 'SKB7'},
            'SKB8.640193': {'description_prep': 'SKB8'}
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
        self.params['Demultiplexed sequences'] = aid

        data = {'user': 'demo@microbio.me',
                'command': dumps(['deblur', '1.1.0', 'Deblur']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        # populate Qiita archive with some precomputed placements
        # placements.json is output from a SEPP run for the resulting Deblur
        # table, but with "tree" value removed for the sake of space
        features = dict()
        with open('support_files/sepp/placements.json', 'r') as f:
            for placement in json.load(f)['placements']:
                fragment = placement['nm'][0][0]
                # exclude 10 sequences to trigger SEPP computation later on
                if fragment not in self.novel_seqs:
                    features[fragment] = json.dumps(placement['p'])
        # add in a feature which should be rejected by SEPP
        features['A' * len(self.novel_seqs[0])] = ""

        # 1) check that archive is currently empty:
        observations = self.qclient.post(
            "/qiita_db/archive/observations/",
            data={'job_id': jid, 'features': list(features.keys())})
        self.assertTrue(len(observations.keys()) == 0)

        # 2) insert placements into archive ...
        self.qclient.patch(url="/qiita_db/archive/observations/",
                           op="add", path=jid,
                           value=json.dumps(features))
        # ... and check that archive does hold those placements now:
        observations = self.qclient.post(
            "/qiita_db/archive/observations/",
            data={'job_id': jid, 'features': list(features.keys())})
        self.assertTrue(len(observations.keys()) == len(features.keys()))

        # 3) execute deblur job with subsequent SEPP run and tiny reference
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)
        self.params['Reference phylogeny for SEPP'] = 'tiny'
        success, ainfo, msg = deblur(self.qclient, jid, self.params, out_dir)
        self.assertEqual("", msg)
        self.assertTrue(success)

        # ensure number of stored placements did grow or at least did not
        observations_2 = self.qclient.post(
            "/qiita_db/archive/observations/",
            data={'job_id': jid,
                  'features': list(features.keys()) + [self.novel_seqs[0]]})
        self.assertTrue(len(observations.keys()) <= len(observations_2.keys()))
        self.maxDiff = None
        # test specific placement values for one fragment that has
        # been pre-populated ...
        exp_placement = (
            '[[226990, -15902.052, 0.14311954, 9.856619e-06, 6.113515e-06], '
            '[226989, -15902.052, 0.14311936, 7.0000096e-06, 6.113515e-06], '
            '[226993, -15902.052, 0.14311917, 8.61664e-06, 6.113515e-06], '
            '[226991, -15902.052, 0.14311911, 6.3553584e-06, 6.113515e-06], '
            '[227443, -15902.052, 0.14311688, 6.7868327e-06, 6.113515e-06], '
            '[226994, -15902.052, 0.14311177, 5.000002e-07, 6.113515e-06], '
            '[227452, -15902.064, 0.14129417, 0.00160019, 6.113515e-06]]'
        )
        self.assertEqual(observations_2[
            ('TACGTAGGGCGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTCA'
             'CGTCGGATGTGAAAGCCCGGGGCTTAACCCCGGGTCTG')], exp_placement)

        # ... and one fragment that was recompted via SEPP during this test
        exp_placement = (
            '[[78, -18489.055, 0.8486466, 0.015792055, 6.113515e-06], '
            '[74, -18491.146, 0.10484001, 0.017408343, 0.010122812], '
            '[77, -18491.959, 0.046513416, 0.015838308, 0.010945947]]')
        self.assertEqual(observations_2[self.novel_seqs[0]], exp_placement)


if __name__ == '__main__':
    main()
