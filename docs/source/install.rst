Installation
------------

The Easy Way (OSX and Linux)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The recommended was to use pygauss is to download the
`Anaconda <http://continuum.io/downloads>`__ Scientific Python
Distribution (64-bit). Once downloaded a new environment can be created
in terminal and pygauss installed:

::

    conda create -n pg_env python=2.7
    conda install -c https://conda.binstar.org/cjs14 -n pg_env pygauss


The Hard Way (Windows)
~~~~~~~~~~~~~~~~~~~~~~

There is currently no pygauss conda distributable for Windows or for
chemlab, which has C-extensions that need to be built using a compiler.
Therefore it will need to be cloned from GitHub. the extensions built,
dependancies installed and finally installed.

::

    conda create -n pg_env python=2.7
    conda install -n pg_env -c https://conda.binstar.org/cjs14 cclib
    conda install -n pg_env -c https://conda.binstar.org/cjs14 chemview
    conda install -n pg_env -c https://conda.binstar.org/cjs14 pyopengl     
    git clone --recursive https://github.com/chemlab/chemlab.git
    cd chemlab
    python setup.py build_ext --inplace
    conda install -n pg_env <pil, pandas, matplotlib, scikit-learn, ...> 
    activate pg_env
    pip install . # or add to PYTHONPATH
    pip install pygauss

If you encounter difficulties it may be useful for you to look in
`working\_conda\_environments <https://github.com/chrisjsewell/PyGauss/tree/master/working_conda_environments>`__
at conda environments known to work.

