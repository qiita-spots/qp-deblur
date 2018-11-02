# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import main
from subprocess import Popen, PIPE

from os import remove
from shutil import rmtree
from tempfile import mkdtemp
from os.path import exists, isdir
from os import environ

from qiita_client.testing import PluginTestCase


class TestCmdGenTree(PluginTestCase):
    def setUp(self):
        self._clean_up_files = []
        self.oldpath = environ['PATH']

    def tearDown(self):
        environ['PATH'] = self.oldpath
        for fp in self._clean_up_files:
            if exists(fp):
                if isdir(fp):
                    rmtree(fp)
                else:
                    remove(fp)

    def test_cmd_gen_tree(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        archive_file = 'support_files/test_archive_file.json'
        ref_template_file = 'support_files/sepp/tmpl_tiny_placement.json'

        p = Popen(["./scripts/generate_tree_from_fragments \
                    --fp_archive %s \
                    --output_dir %s \
                    --fp_ref_template=%s" % (archive_file,
                                             out_dir,
                                             ref_template_file)],
                  shell=True,
                  stdout=PIPE)

        p_out, p_err = p.communicate()
        p_out = p_out.decode("utf-8").rstrip()

        # make sure cmd returns successfully (0)
        self.assertEqual(p.returncode, 0)
        # make sure path to output file matches temporary directory
        self.assertEqual(out_dir, p_out[:len(out_dir)])
        # make sure file is at least named with a .tre extension
        self.assertEqual(p_out[-4:], '.tre')


if __name__ == '__main__':
    main()
