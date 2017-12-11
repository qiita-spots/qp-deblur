# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from subprocess import Popen, PIPE

from os import remove
from os.path import exists, isdir, join
from shutil import copyfile, rmtree
from tempfile import mkstemp, mkdtemp

from qiita_client.testing import PluginTestCase
from qp_deblur import plugin
from qp_deblur.deblur import (generate_sepp_placements)


TESTPREFIX = 'foo'


class seppNativeTests(TestCase):
    def setUp(self):
        # this will allow us to see the full errors
        self.maxDiff = None

    def tearDown(self):
        for file in ["%s_placement.json" % TESTPREFIX,
                     "%s_placement.tog.relabelled.tre" % TESTPREFIX,
                     "%s_placement.tog.relabelled.xml" % TESTPREFIX,
                     "%s_placement.tog.tre" % TESTPREFIX,
                     "%s_placement.tog.xml" % TESTPREFIX,
                     "%s_rename-json.py" % TESTPREFIX,
                     "sepp-%s-err.log" % TESTPREFIX,
                     "sepp-%s-out.log" % TESTPREFIX]:
            try:
                remove(file)
            except OSError:
                pass

    def test_execution(self):
        fp_ref_alignment = join('support_files', 'sepp',
                                'reference_alignment_tiny.fasta')
        fp_ref_pytholgeny = join('support_files', 'sepp',
                                 'reference_phylogeny_tiny.nwk')
        fp_input = join('support_files', 'sepp', 'input_fragments.fasta')

        cmd = 'run-sepp.sh %s %s -a %s -t %s -x 1' % (
            fp_input,
            TESTPREFIX,
            fp_ref_alignment,
            fp_ref_pytholgeny)
        process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        pr_out, pr_err = process.communicate()

        self.assertIn('INFO: All checkpointed executions Finished in',
                      pr_out.decode())
        self.assertEqual(b'', pr_err)

        with open("%s_placement.tog.relabelled.tre" % TESTPREFIX, 'r') as f:
            tree = "\n".join(f.readlines())
            self.assertIn('f__Halomonadaceae', tree)
            self.assertIn('testseqd', tree)


class seppTests(TestCase):
    def setUp(self):
        # this will allow us to see the full errors
        self.maxDiff = None
        self.seqs = [
            ('TACGTAGGATGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGCAGGCGGGACTGTAAG'
             'TTGGATGTGAAATACCGTGGCTTAACCACGGAACTGCATCCAAAACTGTAGTTCTTGAGTG'),
            ('TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGTGGTTTAATAAG'
             'TCTGATGTGAAAGCCCACGGCTCAACCGTGGAGGGTCATTGGAAACTGTCAAACTTGAGTG'),
            ('TACGTAGGTCCCGAGCGTTATCCGGATTTATTGGGCGTAAAGCGAGCGCAGGCGGTTAGATAAG'
             'TCTGAAGTTAAAGGCTGTGGCTTAACCATAGTACGCTTTGGAAACTGTTTTACTTGAGTGC'),
            ('TACGTAGGGAGCGAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGTAGGCGGAATCGCAAG'
             'TCAGATGTGAAAACTATGGGCTTAACCCATAAACTGCATTTGAAACTGTGGTTCTTGAGTG'),
            ('TACGTATGGATCGAGCGTTGTCCGGAATCATTGGGCGTAAAGGGTACGTAGGCGGCCTAGTAAG'
             'TTAGAAGTGAAAGAATATAGCTCAACTATATAAAGCTTTTAAAACTGTTAGGCTTGAGAGA'),
            ('TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGATGGACAAG'
             'TCTGATGTGAAAGGCTGGGGCCCAACCCCGGGACTGCATTGGAAACTGCCCGTCTTGAGTG'),
            ('TACGTAGGGGGCGAGCGTTATCCGGAATGATTGGGCGTAAAGCGCGCGCAGGCGGCCGCTCAAG'
             'CGGGACCTCTAACCCCGGGGCTCAACCTCGGGCCGGGTCCCGAACTGGGCGGCTCGAGTGC'),
            ('TACGGAGGATCCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTTTAGTAAG'
             'TCAGCGGTGAAATTTTGGTGCTTAACACCAAGCGTGCCGTTGATACTGCTGGGCTAGAGAG'),
            ('TACGTAGGGAGCAAGCGTTATCCGGATTTATTGGGTGTAAAGGGTGCGTAGACGGGAATACAAG'
             'TTAGTTGTGAAATACCTCGGCTTAACTGAGGAACTGCAACTAAAACTATATTTCTTGAGTA'),
            ('TACGGAGGGTGCAAGCGTTAATCGGAATCACTGGGCGTAAAGCGCACGTAGGCGGCTTGGTAAG'
             'TCAGGGGTGAGATCCCACAGCCCAACTGTGGAACTGCCTTTGATACTGCCAGGCTTGAGTA')]

    def tearDown(self):
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_generate_sepp_placements(self):
        out_dir = mkdtemp()
        placements = generate_sepp_placements(self.seqs, out_dir, 5)
        #self._clean_up_files.append(out_dir)


if __name__ == '__main__':
    main()
