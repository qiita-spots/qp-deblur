# -----------------------------------------------------------------------------
# Copyright (c) 2014--, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from qiita_client.testing import PluginTestCase

from unittest import main
from subprocess import Popen, PIPE

from os import remove
from shutil import rmtree
from tempfile import mkdtemp
from os.path import exists, isdir, join
from os import environ
import shutil
from hashlib import md5
from json import loads


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
        fp_biom = 'support_files/otu_table.biom'

        # since this cmd alters the target biom, create a copy to out_dir
        fp_output_biom = join(out_dir, 'output_table.biom')
        shutil.copyfile(fp_biom, fp_output_biom)

        p = Popen(["./scripts/generate_tree_from_fragments \
                    --fp_archive %s \
                    --fp_biom %s \
                    --output_dir %s \
                    --fp_ref_template=%s" % (archive_file,
                                             fp_output_biom,
                                             out_dir,
                                             ref_template_file)],
                  shell=True,
                  stdout=PIPE)

        p_out, p_err = p.communicate()
        p_out = p_out.decode("utf-8").rstrip()

        # make sure cmd returns successfully (0)
        self.assertEqual(p.returncode, 0)

        # convert JSON-based output back into a dictionary
        p_out = loads(p_out)

        # make sure path to output file matches temporary directory
        self.assertEqual(out_dir, p_out['archive'][:len(out_dir)])

        # make sure file is at least named with a .tre extension
        self.assertEqual(p_out['archive'][-4:], '.tre')

        # given the known inputs, fp_output_biom should be heavily stripped,
        # and no longer matching fp_biom.
        with open(fp_biom, 'rb') as original_data:
            checksum_original = md5(original_data.read()).hexdigest()

        with open(p_out['biom'], 'rb') as output_data:
            checksum_output = md5(output_data.read()).hexdigest()

        self.assertNotEqual(checksum_original, checksum_output)


if __name__ == '__main__':
    main()
