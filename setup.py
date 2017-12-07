#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2013, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from setuptools import setup
from setuptools.command.install import install
import tarfile
import subprocess

__version__ = "1.0.3"

classes = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering :: Bio-Informatics
    Topic :: Software Development :: Libraries :: Application Frameworks
    Topic :: Software Development :: Libraries :: Python Modules
    Programming Language :: Python
    Programming Language :: Python :: 2.7
    Programming Language :: Python :: Implementation :: CPython
    Operating System :: POSIX :: Linux
    Operating System :: MacOS :: MacOS X
"""

# configuration taken from: https://github.com/biocore/q2-fragment-insertion


def _config_sepp(assets_dir):
    subprocess.run(['python', 'setup.py', 'config', '-c'], check=True,
                   cwd=assets_dir + '/sepp-package/sepp')


class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        import urllib.request
        import shutil
        import os

        # using a tagged version from Siavash's repo
        git_tag = '4.3.4b'
        src_url = ('https://github.com/smirarab/sepp-refs/archive/%s.tar.gz' %
                   git_tag)

        assets_dir = os.path.join(self.install_libbase,
                                  'assets/')

        if not os.path.exists(assets_dir):
            os.mkdir(assets_dir)

        out_f = 'tagged-sepp-package.tar.gz'
        # 1/3: download git tagged version sources ...
        with urllib.request.urlopen(
                src_url) as response, open(out_f, 'wb') as out:
            shutil.copyfileobj(response, out)

        # 2/3: ... which come as one tar archive that needs to be extracted ...
        opened = tarfile.open(out_f, "r:gz")
        opened.extractall(path=self.install_libbase)
        opened.close()

        # 3/3: ... and contains another tar archive which is extracted here.
        opened = tarfile.open(os.path.join(self.install_libbase,
                                           'sepp-refs-%s' % git_tag,
                                           'gg', 'sepp-package.tar.bz'),
                              "r:bz2")
        opened.extractall(path=assets_dir)
        opened.close()

        # copy default taxonomy Greengenes 99%: OTU-ID to lineage
        # LEAVING AS WE MIGHT NEED IT IN THE FUTURE
        # shutil.copy('taxonomy_gg99.qza', assets_dir)

        self.execute(_config_sepp, [assets_dir], 'Configuring SEPP')


with open('README.rst') as f:
    long_description = f.read()

classifiers = [s.strip() for s in classes.split('\n') if s]

setup(name='qp-deblur',
      version=__version__,
      long_description=long_description,
      license="BSD",
      description='Qiita Plugin: Deblur',
      author="Qiita development team",
      author_email="qiita.help@gmail.com",
      url='https://github.com/qiita-spots/qp-deblur',
      test_suite='nose.collector',
      packages=['qp_deblur'],
      scripts=['scripts/configure_deblur', 'scripts/start_deblur'],
      extras_require={'test': ["nose >= 0.10.1", "pep8"]},
      install_requires=['click >= 3.3', 'future', 'deblur'],
      dependency_links=[
        ('https://github.com/qiita-spots/qiita-files/archive/master.zip#'
         'egg=qiita-files-0.1.0-dev')],
      classifiers=classifiers,
      cmdclass={'install': PostInstallCommand})
