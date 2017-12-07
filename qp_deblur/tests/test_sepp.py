# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from pkg_resources import Requirement, resource_filename
import subprocess

from os.path import join
from os import remove

TESTPREFIX = 'foo'


class seppNativeTests(TestCase):
    def test_execution(self):
        fp_sepp_binary = resource_filename(Requirement.parse('qp-deblur'),
                                           'assets/sepp-package/run-sepp.sh')
        fp_ref_alignment = join('support_files', 'sepp',
                                'reference_alignment_tiny.fasta')
        fp_ref_pytholgeny = join('support_files', 'sepp',
                                 'reference_phylogeny_tiny.nwk')
        fp_input = join('support_files', 'sepp', 'input_fragments.fasta')

        cmd = '%s %s %s -a %s -t %s' % (fp_sepp_binary,
                                        fp_input,
                                        TESTPREFIX,
                                        fp_ref_alignment,
                                        fp_ref_pytholgeny)
        with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                              executable="bash") as call_x:
            self.assertTrue(call_x.wait() == 0)
        with open("%s_placement.tog.relabelled.tre" % TESTPREFIX, 'r') as f:
            tree = "\n".join(f.readlines())
            self.assertIn('f__Halomonadaceae', tree)
            self.assertIn('testseqd', tree)

    def tearDown(self):
        remove("%s_placement.json" % TESTPREFIX)
        remove("%s_placement.tog.relabelled.tre" % TESTPREFIX)
        remove("%s_placement.tog.relabelled.xml" % TESTPREFIX)
        remove("%s_placement.tog.tre" % TESTPREFIX)
        remove("%s_placement.tog.xml" % TESTPREFIX)
        remove("%s_rename-json.py" % TESTPREFIX)
        remove("sepp-%s-err.log" % TESTPREFIX)
        remove("sepp-%s-out.log" % TESTPREFIX)


if __name__ == '__main__':
    main()
