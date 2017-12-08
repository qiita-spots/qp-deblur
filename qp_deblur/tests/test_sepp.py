# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
from pkg_resources import Requirement, resource_filename
from subprocess import Popen, PIPE

from os.path import join
from os import remove
import subprocess

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
        assets_dir = resource_filename(Requirement.parse('qp-deblur'), 'assets/')
        print("ASSERTS_DIR >%s<" % assets_dir)
        subprocess.run(['ls', '-laR', assets_dir], check=True, cwd=assets_dir)
        subprocess.run(['cat', assets_dir+'sepp-package/sepp/.sepp/main.config'], check=True, cwd=assets_dir)


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
        process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        pr_out, pr_err = process.communicate()
        print("==========OUT========", pr_out.decode())
        print("==========ERR========", pr_err.decode())

        self.assertIn('INFO: All checkpointed executions Finished in',
                      pr_out.decode())
        self.assertEqual(b'', pr_err)

        with open("%s_placement.tog.relabelled.tre" % TESTPREFIX, 'r') as f:
            tree = "\n".join(f.readlines())
            self.assertIn('f__Halomonadaceae', tree)
            self.assertIn('testseqd', tree)


if __name__ == '__main__':
    main()
