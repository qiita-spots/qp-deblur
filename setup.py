#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2013, The Qiita Development Team.
#
# Distributed under the terms of the BSD 3-clause License.
#
# The full license is in the file LICENSE, distributed with this software.
# -----------------------------------------------------------------------------

from setuptools import setup

__version__ = "1.1.0"

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
      package_data={'qp_deblur': [
          '../support_files/sepp/*.json',
          '../support_files/sepp/*.py',
          '../support_files/sepp/reference_alignment_tiny.fasta',
          '../support_files/sepp/reference_phylogeny_tiny.nwk']},
      scripts=['scripts/configure_deblur', 'scripts/start_deblur',
               'scripts/generate_tree_from_fragments'],
      extras_require={'test': ["nose >= 0.10.1", "pep8"]},
      install_requires=['click >= 3.3', 'future', 'deblur>=1.1.0'],
      dependency_links=[
          ('https://github.com/qiita-spots/qiita-files/archive/master.zip#'
           'egg=qiita-files-0.1.0-dev')],
      classifiers=classifiers,
      )
