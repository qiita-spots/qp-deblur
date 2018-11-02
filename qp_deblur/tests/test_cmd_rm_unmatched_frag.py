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
from os.path import exists, isdir, join
from os import environ
import shutil
from hashlib import md5

from qiita_client.testing import PluginTestCase


class TestCmdRemoveUnmatchedFragments(PluginTestCase):
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

    def test_cmd_rm_unmatched_frag(self):
        out_dir = mkdtemp()
        self._clean_up_files.append(out_dir)

        fp_phylogeny = 'support_files/insertion_tree.relabelled.tre'
        fp_biom = 'support_files/otu_table.biom'

        # since this cmd alters the target biom, create a copy to out_dir
        fp_output_biom = join(out_dir, 'output_table.biom')
        shutil.copyfile(fp_biom, fp_output_biom)

        p = Popen(["./scripts/remove_unmatched_fragments_from_biom \
                    --fp_phylogeny %s \
                    --fp_biom %s" % (fp_phylogeny,
                                     fp_output_biom)],
                  shell=True,
                  stdout=PIPE)

        # wait for cmd to complete, collect stderr, stdout
        p_out, p_err = p.communicate()
        p_out = p_out.decode("utf-8").rstrip()

        # assert cmd returned successfully (0)
        self.assertEqual(p.returncode, 0)

        # given the known inputs, fp_output_biom should be heavily stripped,
        # and no longer matching fp_biom.
        with open(fp_biom, 'rb') as original_data:
            checksum_original = md5(original_data.read()).hexdigest()
            print(checksum_original)

        with open(fp_output_biom, 'rb') as output_data:
            checksum_output = md5(output_data.read()).hexdigest()
            print(checksum_output)

        self.assertNotEqual(checksum_original, checksum_output)


if __name__ == '__main__':
    main()
