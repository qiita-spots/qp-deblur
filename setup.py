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


class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        import shutil
        import os

        install.run(self)

        assets_dir = os.path.join(self.install_libbase, 'qp_deblur/assets/')
        if not os.path.exists(assets_dir):
            os.mkdir(assets_dir)

        shutil.copy(os.path.join('support_files', 'sepp',
                                 'tmpl_gg13.8-99_placement.json'), assets_dir)
        shutil.copy(os.path.join('support_files', 'sepp',
                                 'tmpl_gg13.8-99_rename-json.py'), assets_dir)


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
