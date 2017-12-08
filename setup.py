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
from setuptools.command.develop import develop
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


def _initial():
    import shutil

    if shutil.which('java') is None:
        raise ValueError('java not found')

    result = subprocess.run(['java', '-version'], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    err_content = result.stderr.decode('ascii')
    if ('java version' not in err_content) and \
       ('jdk version' not in err_content):
        raise ValueError(('Please verify that java is installed and working. '
                          'As a first test, please execute "java -version" '
                          'and make sure the output shows there is an actual '
                          'version installed. OSX lies. If Java needs to be '
                          'installed, please obtain the 1.8 or greater Java '
                          'Runtime Environment (JRE) from Oracle.com. A '
                          'google search for "download jre" is likely to be '
                          'sufficient.'))


def _post(obj):
    import urllib.request
    import shutil
    import os

    # using a tagged version from Siavash's repo
    git_tag = '4.3.4b'
    src_url = ('https://github.com/smirarab/sepp-refs/archive/%s.tar.gz' %
               git_tag)

    assets_dir = os.path.join(obj.install_libbase,
                              'qp_deblur/assets/')

    if not os.path.exists(assets_dir):
        os.mkdir(assets_dir)

    out_f = 'tagged-sepp-package.tar.gz'
    # 1/3: download git tagged version sources ...
    with urllib.request.urlopen(
            src_url) as response, open(out_f, 'wb') as out:
        shutil.copyfileobj(response, out)

    # 2/3: ... which come as one tar archive that needs to be extracted ...
    opened = tarfile.open(out_f, "r:gz")
    opened.extractall(path=obj.install_libbase)
    opened.close()

    # 3/3: ... and contains another tar archive which is extracted here.
    opened = tarfile.open(os.path.join(obj.install_libbase,
                                       'sepp-refs-%s' % git_tag,
                                       'gg', 'sepp-package.tar.bz'),
                          "r:bz2")
    opened.extractall(path=assets_dir)
    opened.close()

    # copy default taxonomy Greengenes 99%: OTU-ID to lineage
    # LEAVING AS WE MIGHT NEED IT IN THE FUTURE
    # shutil.copy('taxonomy_gg99.qza', assets_dir)

    # copy patch file
    name_patch = 'onlyplacements.patch'
    shutil.copy(os.path.join('support_files', 'sepp', name_patch),
                assets_dir)
    name_patch2 = 'debug.patch'
    shutil.copy(os.path.join('support_files', 'sepp', name_patch2),
                assets_dir)

    obj.execute(_patch_sepp, [assets_dir, name_patch],
                 'Patch run-sepp.sh')
    obj.execute(_patch_sepp, [assets_dir, name_patch2],
                 'Patch run-sepp.sh debug')
    obj.execute(_config_sepp, [assets_dir], 'Configuring SEPP')


def _patch_sepp(assets_dir, name_patch):
    subprocess.run(['patch', 'sepp-package/run-sepp.sh', name_patch],
                   check=True, cwd=assets_dir)


def _config_sepp(assets_dir):
    subprocess.run(['python', 'setup.py', 'config', '-c'], check=True,
                   cwd=assets_dir + '/sepp-package/sepp')


class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        _initial()
        install.run(self)
        _post(self)


class PostDevelopCommand(develop):
    """Post-installation for development mode."""
    def run(self):
        _initial()
        develop.run(self)
        _post(self)

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
      cmdclass={'install': PostInstallCommand,
                'develop': PostDevelopCommand})
