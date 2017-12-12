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
from os.path import join, abspath
from shutil import rmtree
from tempfile import mkdtemp

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
        fp_ref_alignment = abspath(join('support_files', 'sepp',
                                        'reference_alignment_tiny.fasta'))
        fp_ref_pytholgeny = abspath(join('support_files', 'sepp',
                                         'reference_phylogeny_tiny.nwk'))
        fp_input = abspath(join('support_files', 'sepp',
                                'input_fragments.fasta'))

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
        self.exp = {
            self.seqs[6]:
            [[957, -25552.186, 1, 0.06003037, 0.032140907]],
            self.seqs[5]:
            [[896, -25279.383, 1, 8.754862e-06, 0.008639097]],
            self.seqs[3]:
            [[795, -25287.547, 0.76562905, 0.0036086612, 0.02935211],
             [796, -25289.215, 0.14446084, 0.04416296, 0.02090377],
             [797, -25289.69, 0.089910105, 7.837499e-06, 0.02483988]],
            self.seqs[9]:
            [[458, -25315.588, 0.91818565, 0.07924193, 0.0701725],
             [462, -25318.62, 0.04435069, 0.0167353, 0.092968106],
             [461, -25318.787, 0.03746365, 0.022151958, 0.097728774]],
            self.seqs[7]:
            [[414, -25326.35, 0.9521992, 0.12123759, 0.11371933],
             [379, -25329.64, 0.035452154, 0.09661898, 0.15217614],
             [377, -25330.695, 0.012348609, 0.13954082, 0.1661857]],
            self.seqs[2]:
            [[756, -25374.973, 0.9637653, 5.5463693e-06, 0.008664635],
             [761, -25379.564, 0.009759569, 0.016187644, 0.01309992],
             [757, -25379.879, 0.0071310354, 5.0176413e-06, 0.017295655],
             [758, -25379.967, 0.0065307072, 5.019576e-06, 0.017373659],
             [760, -25379.975, 0.0064835623, 7.6980505e-06, 0.017384922],
             [759, -25379.998, 0.006329803, 2.5000013e-06, 0.01740473]],
            self.seqs[4]:
            [[781, -24898.662, 0.40177065, 0.022466354, 0.096478306],
             [780, -24898.957, 0.2992006, 0.028304012, 0.11486658],
             [779, -24898.957, 0.29902875, 0.038573176, 0.11486842]],
            self.seqs[8]:
            [[967, -25336.166, 0.81720763, 8.235129e-06, 0.22352734],
             [804, -25338.336, 0.09330978, 0.018047476, 0.2054955],
             [794, -25339.684, 0.024266021, 0.087570414, 0.19646433],
             [965, -25339.695, 0.023959905, 0.071342826, 0.2269426],
             [802, -25340.07, 0.016466064, 0.07623818, 0.21268782],
             [962, -25340.195, 0.014542027, 0.15002388, 0.21349123],
             [975, -25340.545, 0.010248571, 0.004942888, 0.22604735]],
            self.seqs[1]:
            [[765, -25280.664, 1, 0.038505234, 6.113515e-06]],
            self.seqs[0]:
            [[964, -25330.373, 0.28328982, 0.052007545, 0.15232158],
             [805, -25330.59, 0.22800335, 6.190029e-06, 0.19001018],
             [804, -25330.59, 0.22790368, 0.030604076, 0.1900029],
             [955, -25330.838, 0.177996, 6.588367e-06, 0.19621426],
             [823, -25332.365, 0.038650226, 0.044555224, 0.18255654],
             [962, -25332.752, 0.02623141, 0.14301622, 0.16091308],
             [932, -25333.133, 0.017925516, 0.17963046, 0.16220896]]}

    def test_generate_sepp_placements(self):
        fp_ref_alignment = abspath(join('support_files', 'sepp',
                                        'reference_alignment_tiny.fasta'))
        fp_ref_phylogeny = abspath(join('support_files', 'sepp',
                                        'reference_phylogeny_tiny.nwk'))
        out_dir = mkdtemp()
        # testing default references might take up to 10 minutes of compute
        placements = generate_sepp_placements(
            self.seqs, out_dir, 1, reference_alignment=fp_ref_alignment,
            reference_phylogeny=fp_ref_phylogeny)

        # test that every sequence has a placement
        for seq in self.seqs:
            self.assertIn(seq, placements)

        # test that placement for seq nr 7 has only one placement ...
        self.assertEqual(len(placements[self.seqs[6]]), 1)
        # ... and the sequence gets placed into node 957
        self.assertEqual(placements[self.seqs[6]][0][0], 957)

        # test the full output
        self.assertCountEqual(placements, self.exp)

        # clean up working directory
        rmtree(out_dir)

    def test_generate_sepp_placements_noseqs(self):
        self.assertEqual(generate_sepp_placements([], None, 1), {})

    def test_generate_sepp_placements_nonzero(self):
        out_dir = mkdtemp()
        self.assertRaisesRegex(
            ValueError, "Error running run-sepp.sh", generate_sepp_placements,
            self.seqs, out_dir, 1, reference_phylogeny='/dev/null')

        rmtree(out_dir)


if __name__ == '__main__':
    main()
