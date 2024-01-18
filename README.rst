Deblur Qiita Plugin
====================

|Build Status| |Coverage Status|

Qiita (canonically pronounced *cheetah*) is an analysis environment for microbiome (and other "comparative -omics") datasets.

This package includes the deblur functionality for Qiita.

`Deblur <https://github.com/biocore/deblur>`__ is a greedy deconvolution algorithm based on Illumina Miseq/Hiseq error profiles.

Note that deblur has vsearch, mafft, SortMeRNA==2.0 and fragment-insertion as requirements but these are not installed as part of the install.

.. |Build Status| image:: https://travis-ci.org/qiita-spots/qp-deblur.svg?branch=master
   :target: https://travis-ci.org/qiita-spots/qp-deblur
.. |Coverage Status| image:: https://coveralls.io/repos/github/qiita-spots/qp-deblur/badge.svg?branch=master
   :target: https://coveralls.io/github/qiita-spots/qp-deblur?branch=master
