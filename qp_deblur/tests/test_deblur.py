# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from os import close, remove, chmod
from shutil import copyfile, rmtree
from tempfile import mkstemp, mkdtemp
from json import dumps, load
from os.path import exists, isdir, join
from os import environ

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

        # to prevent test from timing out, we need to pre-populate the archive
        # with placements for Deblur fragments to avoid long running SEPP.
        # Therefore, we need a job-id to infer merging scheme which can only be
        # done after demultiplexing job is created. Thus, actuall population
        # needs to be done within the test.
        self.features = dict()
        with open('support_files/sepp/placements.json', 'r') as f:
            for placement in load(f)['placements']:
                self.features[placement['nm'][0][0]] = \
                    dumps(placement['p'])

        # saving current value of PATH
        self.oldpath = environ['PATH']

    def tearDown(self):
        # restore eventually changed PATH env var
        environ['PATH'] = self.oldpath
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
               '--threads-per-sample "1" '
               '--trim-length "100"')
        obs = generate_deblur_workflow_commands(
            ['fastq/s1.fastq'], 'output', self.params)

        self.assertEqual(obs, exp)

    def test_deblur_no_tree(self):
        # generating filepaths
        fd, fp = mkstemp(suffix='_seqs.demux')
        close(fd)
        self._clean_up_files.append(fp)
        copyfile('support_files/no_sepp_seqs.demux', fp)
        prep_info_dict = {
            'SKB7.640196': {'description_prep': 'SKB7'},
            'SKB8.640193': {'description_prep': 'SKB8'}
        }
        data = {'prep_info': dumps(prep_info_dict),
                # magic #1 = testing study
                'study': 1,
                'data_type': 'ITS'}
        pid = self.qclient.post('/apitest/prep_template/', data=data)['prep']

        # inserting artifacts
        data = {
            'filepaths': dumps([(fp, 'preprocessed_fastq')]),
            'type': "Demultiplexed",
            'name': "New demultiplexed artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['Demultiplexed sequences'] = aid

        data = {'user': 'demo@microbio.me',
                'command': dumps(['deblur', '1.0.4', 'Deblur']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # pre-populate archive with fragment placements
        self.qclient.patch(url="/qiita_db/archive/observations/",
                           op="add", path=jid,
                           value=dumps(self.features))
        success, ainfo, msg = deblur(self.qclient, jid, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        self.assertEqual("BIOM", ainfo[0].artifact_type)
        self.assertEqual(1, len(ainfo))

        self.assertEqual(
            [(join(out_dir, 'deblur_out', 'all.biom'), 'biom'),
             (join(out_dir, 'deblur_out', 'all.seqs.fa'),
              'preprocessed_fasta')], ainfo[0].files)

    def test_deblur(self):
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
            'filepaths': dumps([(fp, 'preprocessed_fastq')]),
            'type': "Demultiplexed",
            'name': "New demultiplexed artifact",
            'prep': pid}
        aid = self.qclient.post('/apitest/artifact/', data=data)['artifact']

        self.params['Demultiplexed sequences'] = aid

        data = {'user': 'demo@microbio.me',
                'command': dumps(['deblur', '1.0.4', 'Deblur']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # pre-populate archive with fragment placements
        self.qclient.patch(url="/qiita_db/archive/observations/",
                           op="add", path=jid,
                           value=dumps(self.features))
        success, ainfo, msg = deblur(self.qclient, jid, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        self.assertEqual("BIOM", ainfo[0].artifact_type)

        self.assertEqual(
            [(join(out_dir, 'deblur_out', 'all.biom'), 'biom'),
             (join(out_dir, 'deblur_out', 'all.seqs.fa'),
              'preprocessed_fasta')], ainfo[0].files)

    def test_deblur_demux(self):
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
                'command': dumps(['deblur', '1.0.4', 'Deblur']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # pre-populate archive with fragment placements
        self.qclient.patch(url="/qiita_db/archive/observations/",
                           op="add", path=jid,
                           value=dumps(self.features))
        success, ainfo, msg = deblur(self.qclient, jid, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)

        self.assertEqual("BIOM", ainfo[0].artifact_type)
        self.assertEqual("BIOM", ainfo[1].artifact_type)

        self.assertEqual(
            [(join(out_dir, 'deblur_out', 'deblured', 'all.biom'), 'biom'),
             (join(out_dir, 'deblur_out', 'deblured', 'all.seqs.fa'),
              'preprocessed_fasta')], ainfo[0].files)

        file_tree = join(out_dir, 'deblur_out',
                         'deblured', 'insertion_tree.relabelled.tre')
        self.assertEqual(
            [(join(out_dir, 'deblur_out', 'deblured', 'reference-hit.biom'),
              'biom'),
             (join(out_dir, 'deblur_out', 'deblured',
                   'reference-hit.seqs.fa'), 'preprocessed_fasta'),
             (file_tree, 'plain_text')
             ], ainfo[1].files)

        # finally we are gonna test that the patch is added to the tree
        with open(file_tree, 'r') as tree_fp:
            tree = tree_fp.read()
            self.assertTrue(tree.endswith("'k__Bacteria':0.0);\n"))

    def test_deblur_failingbin(self):
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
                'command': dumps(['deblur', '1.0.4', 'Deblur']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # create a fake deblur binary that will always fail
        fp_fake_deblur = join(out_dir, 'deblur')
        with open(fp_fake_deblur, 'w') as f:
            f.write('#!/bin/bash\nexit 123\n')
        chmod(fp_fake_deblur, 0o775)
        environ['PATH'] = '%s:%s' % (out_dir, self.oldpath)
        success, ainfo, msg = deblur(self.qclient, jid, self.params, out_dir)
        self.assertFalse(success)
        self.assertEqual(ainfo, None)
        self.assertIn('Error running deblur', msg)

    def test_deblur_failing_guppy(self):
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
                'command': dumps(['deblur', '1.0.4', 'Deblur']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # pre-populate archive with fragment placements
        self.qclient.patch(url="/qiita_db/archive/observations/",
                           op="add", path=jid,
                           value=dumps(self.features))
        # create a fake guppy binary that will always fail
        fp_fake_guppy = join(out_dir, 'guppy')
        with open(fp_fake_guppy, 'w') as f:
            f.write('#!/bin/bash\nexit 123\n')
        chmod(fp_fake_guppy, 0o775)
        environ['PATH'] = '%s:%s' % (out_dir, self.oldpath)
        success, ainfo, msg = deblur(self.qclient, jid, self.params, out_dir)
        self.assertFalse(success)
        self.assertEqual(ainfo, None)
        self.assertIn('Error running guppy', msg)

    def test_deblur_keyerror(self):
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
                'command': dumps(['deblur', '1.0.4', 'Deblur']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # pre-populate archive with fragment placements
        # make sure that at least one sequence got no placements via SEPP
        self.features[('TACGGAGGGTGCAAGCGTTATCCGGATTCACTGGGTTTAAAGGGTGCGTAGGT'
                       'GGGTTGGTAAGTCAGTGGTGAAATCTCCGGGCTTAACTCGGAAACTG')] = ''
        self.qclient.patch(url="/qiita_db/archive/observations/",
                           op="add", path=jid,
                           value=dumps(self.features))
        success, ainfo, msg = deblur(self.qclient, jid, self.params, out_dir)

        self.assertEqual("", msg)
        self.assertTrue(success)


class deblurSilva(PluginTestCase):
    def setUp(self):
        # this will allow us to see the full errors
        self.maxDiff = None

        plugin("https://localhost:21174", 'register', 'ignored')
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
            'Reference phylogeny for SEPP': 'Silva_12.8'}
        self._clean_up_files = []

        # saving current value of PATH
        self.oldpath = environ['PATH']

    def tearDown(self):
        # restore eventually changed PATH env var
        environ['PATH'] = self.oldpath
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_deblur_silva(self):
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
                'command': dumps(['deblur', '1.0.4', 'Deblur']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # create a fake run-sepp.sh binary.
        # Basically, I want to test here if the right filepaths are passed to
        # run-sepp.sh when requesting silva instead of the default reference.
        # Unfortunately, a real computation would take too long for Travis.
        # Thus, I replace the actual binary script here with something that
        # takes the parameters and echos them to stderr and fails. This way,
        # we can expect to find file paths in the msg error message.
        fp_fake_sepp = join(out_dir, 'run-sepp.sh')
        with open(fp_fake_sepp, 'w') as f:
            f.write('#!/bin/bash\n\n'
                    'while [[ $# -gt 0 ]]\n'
                    'do\n'
                    '	key="$1"\n'
                    '	case $key in\n'
                    '		-a)\n'
                    '			a="$2"\n'
                    '			shift # past argument\n'
                    '			;;\n'
                    '		-t)\n'
                    '			t="$2"\n'
                    '			shift # past argument\n'
                    '			;;\n'
                    '		-r)\n'
                    '			r="$2"\n'
                    '			shift # past argument\n'
                    '			;;\n'
                    '		*)\n'
                    '			opts="$opts"" ""$key"" ""$2"\n'
                    '			shift # past argument\n'
                    '			;;\n'
                    '	esac\n'
                    '	shift # past argument or value\n'
                    'done\n'
                    '>&2 echo "info file \'$r\'"\n'
                    '>&2 echo "tree file \'$t\'"\n'
                    '>&2 echo "alignment file \'$a\'"\n'
                    'exit 1\n')
        chmod(fp_fake_sepp, 0o775)
        environ['PATH'] = '%s:%s' % (out_dir, self.oldpath)
        success, ainfo, msg = deblur(self.qclient, jid, self.params, out_dir)
        self.assertFalse(success)
        self.assertEqual(ainfo, None)
        self.assertIn(('share/fragment-insertion/ref/'
                       'silva12.8_99otus_aligned_masked1977.tre'), msg)
        self.assertIn(('share/fragment-insertion/ref/'
                       'silva12.8_99otus_aligned_masked1977.fasta'), msg)
        self.assertIn(('share/fragment-insertion/ref/'
                       'silva12.8_99otus_aligned_masked1977.info'), msg)


class deblurTests_binaryfail(PluginTestCase):
    def setUp(self):
        # this will allow us to see the full errors
        self.maxDiff = None

        plugin("https://localhost:21174", 'register', 'ignored')
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

        # saving current value of PATH
        self.oldpath = environ['PATH']

    def tearDown(self):
        # restore eventually changed PATH env var
        environ['PATH'] = self.oldpath
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_deblur_failing_sepp(self):
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
                'command': dumps(['deblur', '1.0.4', 'Deblur']),
                'status': 'running',
                'parameters': dumps(self.params)}
        jid = self.qclient.post('/apitest/processing_job/', data=data)['job']

        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        # create a fake sepp binary that will always fail
        fp_fake_sepp = join(out_dir, 'run-sepp.sh')
        with open(fp_fake_sepp, 'w') as f:
            f.write('#!/bin/bash\nexit 123\n')
        chmod(fp_fake_sepp, 0o775)
        environ['PATH'] = '%s:%s' % (out_dir, self.oldpath)
        success, ainfo, msg = deblur(self.qclient, jid, self.params, out_dir)
        self.assertFalse(success)
        self.assertEqual(ainfo, None)
        self.assertIn('Error running run-sepp.sh', msg)


if __name__ == '__main__':
    main()
