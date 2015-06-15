.. PyGauss documentation master file, created by
   sphinx-quickstart on Sun Jun 14 01:13:38 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyGauss!
===================================

+--------------------------+-----------------------------------------------+
|**Author**                | Chris Sewell                                  |
+--------------------------+-----------------------------------------------+
|**Project Page**          | https://github.com/chrisjsewell/PyGauss       |
+--------------------------+-----------------------------------------------+
|**Conda Distro**          | https://conda.binstar.org/cjs14               |
+--------------------------+-----------------------------------------------+
|**PyPi Distro**           | https://pypi.python.org/pypi/pygauss          |
+--------------------------+-----------------------------------------------+
|**Communication**         | pygauss@googlegroups.com                      |
+--------------------------+-----------------------------------------------+

.. include:: mission_statement.txt

Intitally PyGauss has been designed for the purpose of examining one or more
`Gaussian <http://www.gaussian.com/>`__ quantum chemical
computations, both **geometrically** and **electronically**. It is built on top of the
`cclib <http://cclib.github.io/>`__/`chemview <http://chemview.readthedocs.org/en/latest/>`__/`chemlab <http://chemlab.readthedocs.org/en/latest/index.html>`__
suite of packages and python scientific stack though, and so should be extensible 
to other types of computational chemical analysis. 
PyGauss is primarily designed to be used interactively in the 
`IPython Notebook <http://ipython.org/notebook.html>`__. 

As shown in the examples, a molecular optimisation can be assesed
individually (much like in
`gaussview <http://www.gaussian.com/g_prod/gv5b.htm>`__), but also as
part of a group. The advantages of this package are then:

-  Faster, more efficient analysis
-  Extensible analysis
-  Reproducible analysis

Contents
--------

.. toctree::
   :maxdepth: 2

   install
   example
   status
   history
   enhancements
   user_api


License
-------

Pygauss is released under the `GNU GPLv3
<http://www.gnu.org/licenses/gpl.html>`_ and its main developer is
Chris Sewell.
