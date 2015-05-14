
Python Gaussian Analysis Tool (PyGauss)
=======================================

PyGauss is designed to be an API for parsing one or more input/output
files from a `Gaussian <http://www.gaussian.com/>`__ quantum chemical
computation and provide functionality to assess **molecular geometry**
and **electronic distribution** both visually and quantitatively.

It is built on top of the
`cclib <http://cclib.github.io/>`__/`chemview <http://chemview.readthedocs.org/en/latest/>`__/`chemlab <http://chemlab.readthedocs.org/en/latest/index.html>`__
suite of packages and python scientific stack and is primarily designed
to be used interactively in the `IPython
Notebook <http://ipython.org/notebook.html>`__ (within which this readme
was created). As shown below, a molecular optimisation can be assesed
individually (much like in
`gaussview <http://www.gaussian.com/g_prod/gv5b.htm>`__), but also as
part of a group. The advantages of this package are then:

-  Faster, more efficient analysis
-  Reproducible analysis
-  Trend analysis

Instillation
------------

-  The source code is hosted on GitHub;
   https://github.com/chrisjsewell/PyGauss
-  A PyPi distribution is available at;
   https://pypi.python.org/pypi/pygauss
-  A Conda distribution is available at; https://conda.binstar.org/cjs14

The Easy Way (OSX)
~~~~~~~~~~~~~~~~~~

The recommended was to use pygauss is to download the
`Anaconda <http://continuum.io/downloads>`__ Scientific Python
Distribution (64-bit). Once downloaded a new environment can be created
in terminal and pygauss installed:

::

    conda create -n pg_env python=2.7
    conda install -c https://conda.binstar.org/cjs14 -n pg_env pygauss

The Middle Road (Linux)
~~~~~~~~~~~~~~~~~~~~~~~

There is currently no pygauss conda distributable for Linux, but there
is for chemlab. So chemlab can be installed, then install a few
dependancies that pip finds difficult / doesn't have, and finally
install pygauss using pip (make sure to activate the required
environment)

::

    conda create -n pg_env python=2.7
    conda install -n pg_env -c https://conda.binstar.org/cjs14 chemlab  
    conda install -n pg_env <pil, pandas, matplotlib, scikit-learn> 
    activate pg_env
    pip install pygauss

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

Example Assessment
------------------

You should then be able to open an assessment in IPython Notebook
starting with the following:

.. code:: python

    from IPython.display import display
    %matplotlib inline
    import pygauss as pg
    folder = pg.get_test_folder()
    pg.__version__




.. parsed-literal::

    '0.3.0'



Single Molecule Analysis
~~~~~~~~~~~~~~~~~~~~~~~~

A *molecule* can be created containg data about the inital geometry,
optimisation process and analysis of the final configuration. Molecules
can be viewed statically or interactively (not currently supported by
Firefox).

.. code:: python

    mol = pg.molecule.Molecule(folder,
                    init_fname='CJS1_emim-cl_B_init.com', 
                    opt_fname=['CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_opt-modredundant_difrz.log',
                               'CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_opt-modredundant_difrz_err.log',
                               'CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_opt-modredundant_unfrz.log'],
                    freq_fname='CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_freq_unfrz.log',
                    nbo_fname='CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_pop-nbo-full-_unfrz.log', 
                    alignto=[3,2,1])
    
    #mol.show_initial(active=True)
    display(mol.show_initial(zoom=0.5, rotations=[[0,0,90], [-90, 90, 0]]))
    display(mol.show_optimisation(ball_stick=True, rotations=[[0,0,90], [-90, 90, 0]]))



.. image:: output_8_0.png



.. image:: output_8_1.png


Basic analysis of optimisation...

.. code:: python

    print('Optimised? {0}, Conformer? {1}, Energy = {2} a.u.'.format(
        mol.is_optimised(), mol.is_conformer(), round(mol.get_optimisation_E(units='hartree'),3)))
    ax = mol.plot_optimisation_E(units='hartree')
    ax.get_figure().set_size_inches(3, 2)


.. parsed-literal::

    Optimised? True, Conformer? True, Energy = -805.105 a.u.
    


.. image:: output_10_1.png


Geometric analysis...

.. code:: python

    print 'Cl optimised polar coords from aromatic ring : ({0}, {1},{2})'.format(
        *[round(i, 2) for i in mol.calc_polar_coords_from_plane(20,3,2,1)])
    ax = mol.plot_opt_trajectory(20, [3,2,1])
    ax.set_title('Cl optimisation path')
    ax.get_figure().set_size_inches(4, 3)


.. parsed-literal::

    Cl optimised polar coords from aromatic ring : (0.11, -116.42,-170.06)
    


.. image:: output_12_1.png


Potential Energy Scan analysis of geometric conformers...

.. code:: python

    mol2 = pg.molecule.Molecule(folder, alignto=[3,2,1],
                pes_fname=['CJS_emim_6311_plus_d3_scan.log', 
                           'CJS_emim_6311_plus_d3_scan_bck.log'])   
    ax = mol2.plot_pes_scans([1,4,9,10], rotation=[0,0,90], img_pos='local_maxs', zoom=0.5)
    ax.set_title('Ethyl chain rotational conformer analysis')
    ax.get_figure().set_size_inches(7, 3)



.. image:: output_14_0.png


Natural Bond Orbital and Second Order Perturbation Theory analysis...

.. code:: python

    print '+ve charge centre polar coords from aromatic ring: ({0} {1},{2})'.format(
        *[round(i, 2) for i in mol.calc_nbo_charge_center(3, 2, 1)])
    display(mol.show_nbo_charges(ball_stick=True, axis_length=0.4, 
                                  rotations=[[0,0,90], [-90, 90, 0]]))
    display(mol.show_SOPT_bonds(min_energy=15., rotations=[[0, 0, 90]]))


.. parsed-literal::

    +ve charge centre polar coords from aromatic ring: (0.02 -51.77,-33.15)
    


.. image:: output_16_1.png



.. image:: output_16_2.png


Multiple Computations Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Multiple computations, for instance of different starting conformations,
can be grouped into an *Analysis* class.

.. code:: python

    analysis = pg.Analysis(folder)
    errors = analysis.add_runs(headers=['Cation', 'Anion', 'Initial'], 
                                   values=[['emim'], ['cl'],
                                           ['B', 'BE', 'BM', 'F', 'FE', 'FM']],
                init_pattern='*{0}-{1}_{2}_init.com',
                opt_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_opt*unfrz.log',
                freq_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_freq*.log',
                nbo_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_pop-nbo-full-*.log')
    print 'Read Errors:'
    errors.File


.. parsed-literal::

    Read Errors:
    



.. parsed-literal::

    0                                 *emim-cl_FM_init.com
    1         *emim-cl_FM_6-311+g-d-p-_gd3bj_opt*unfrz.log
    2             *emim-cl_FM_6-311+g-d-p-_gd3bj_freq*.log
    3    *emim-cl_FM_6-311+g-d-p-_gd3bj_pop-nbo-full-*.log
    Name: File, dtype: object



**New Feature:** you can now access files on a server over ssh in the
following manner:

::

    analysis = pg.Analysis( '/path/to/folder', 
                    ssh_server='login.server.com',
                    ssh_username='username')

The methods mentioned for indivdiual molecules can then be applied to
all or a subset of these computations.

.. code:: python

    analysis.add_mol_property_subset('Opt', 'is_optimised', rows=[2,3])
    analysis.add_mol_property('Energy (au)', 'get_optimisation_E', units='hartree')
    analysis.add_mol_property('Cation chain, $\\psi$', 'calc_dihedral_angle', [1, 4, 9, 10])
    analysis.add_mol_property('Cation Charge', 'calc_nbo_charge', range(1, 20))
    analysis.add_mol_property('Anion Charge', 'calc_nbo_charge', [20])
    analysis.add_mol_property(['Anion-Cation, $r$', 'Anion-Cation, $\\theta$', 'Anion-Cation, $\\phi$'], 
                                   'calc_polar_coords_from_plane', 3, 2, 1, 20)
    analysis.get_table(row_index=['Anion', 'Cation', 'Initial'], 
                       column_index=['Cation', 'Anion', 'Anion-Cation'])
    analysis




.. parsed-literal::

      Anion Cation Initial   Opt  Energy (au)  Cation chain, $\psi$  Cation Charge  Anion Charge  Anion-Cation, $r$  Anion-Cation, $\theta$  Anion-Cation, $\phi$
    0    cl   emim       B   NaN     -805.105                80.794          0.888        -0.888              0.420                -123.392               172.515
    1    cl   emim      BE   NaN     -805.105                80.622          0.887        -0.887              0.420                -123.449               172.806
    2    cl   emim      BM  True     -805.104                73.103          0.874        -0.874              0.420                 124.121              -166.774
    3    cl   emim       F  True     -805.118               147.026          0.840        -0.840              0.420                  10.393                 0.728
    4    cl   emim      FE   NaN     -805.117                85.310          0.851        -0.851              0.417                 -13.254                -4.873



**NEW FEATURE:** there is now an option (requiring
`pdflatex <http://www.tug.org/applications/pdftex/>`__ and
`ghostscript <http://www.ghostscript.com/download/gsdnld.html>`__\ +\ `imagemagik <http://www.imagemagick.org/script/binary-releases.php>`__)
to output the tables as a latex formatted image.

.. code:: python

    analysis.get_table(row_index=['Anion', 'Cation', 'Initial'],
                       column_index=['Cation', 'Anion', 'Anion-Cation'],
                       as_image=True, font_size=12)




.. image:: output_24_0.png



RadViz is a way of visualizing multi-variate data.

.. code:: python

    ax = analysis.plot_radviz_comparison('Anion', columns=range(4, 10))



.. image:: output_26_0.png


The KMeans algorithm clusters data by trying to separate samples into n
groups of equal variance.

.. code:: python

    kwargs = {'mtype':'optimised', 'align_to':[3,2,1], 
                'rotations':[[0, 0, 90], [-90, 90, 0]],
                'axis_length':0.3}
    pg.utils.iprint_kmean_groups(analysis, 'Anion', 'cl', 4, 
                                 range(4, 10), output=['Initial'],
                                 **kwargs)


.. parsed-literal::

    -------------
    Category 0:
    -------------
    Initial: B
    


.. image:: output_28_1.png


.. parsed-literal::

    Initial: BE
    


.. image:: output_28_3.png


.. parsed-literal::

    -------------
    Category 1:
    -------------
    Initial: BM
    


.. image:: output_28_5.png


.. parsed-literal::

    -------------
    Category 2:
    -------------
    Initial: FE
    


.. image:: output_28_7.png


.. parsed-literal::

    -------------
    Category 3:
    -------------
    Initial: F
    


.. image:: output_28_9.png


MORE TO COME!!
