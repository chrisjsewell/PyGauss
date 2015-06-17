
Example Assessment
------------------

After installing PyGauss you should be able to open this IPython
Notebook from;
https://github.com/chrisjsewell/PyGauss/blob/master/Example\_Assessment.ipynb,
and run the following...

.. code:: python

    from IPython.display import display, Image
    %matplotlib inline
    import pygauss as pg
    print 'pygauss version: {}'.format(pg.__version__)


.. parsed-literal::

    pygauss version: 0.4.3
    

The test folder has a number of example Gaussian outputs to play around
with.

.. code:: python

    folder = pg.get_test_folder()
    len(folder.list_files())




.. parsed-literal::

    37



**Note:** the *folder* object will act identical whether using a local
path or one on a server over ssh (using
`paramiko <http://www.paramiko.org/>`__):

::

    folder = pg.Folder('/path/to/folder', 
                    ssh_server='login.server.com',
                    ssh_username='username')

Single Molecule Analysis
~~~~~~~~~~~~~~~~~~~~~~~~

A *molecule* can be created containg data about the inital geometry,
optimisation process and analysis of the final configuration. Molecules
can be viewed statically or interactively.

.. code:: python

    mol = pg.molecule.Molecule(folder_obj=folder,
                    init_fname='CJS1_emim-cl_B_init.com', 
                    opt_fname=['CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_opt-modredundant_difrz.log',
                               'CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_opt-modredundant_difrz_err.log',
                               'CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_opt-modredundant_unfrz.log'],
                    freq_fname='CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_freq_unfrz.log',
                    nbo_fname='CJS1_emim-cl_B_6-311+g-d-p-_gd3bj_pop-nbo-full-_unfrz.log', 
                    atom_groups={'emim':range(20), 'cl':[20]},
                    alignto=[3,2,1])
    
    #mol.show_initial(active=True)
    vdw = mol.show_initial(represent='vdw', rotations=[[0,0,90], [-90, 90, 0]])
    ball_stick = mol.show_optimisation(represent='ball_stick', rotations=[[0,0,90], [-90, 90, 0]])
    display(vdw, ball_stick)



.. image::  images/output_7_0.png



.. image::  images/output_7_1.png


Basic analysis of optimisation...

.. code:: python

    print('Optimised? {0}, Conformer? {1}, Energy = {2} a.u.'.format(
        mol.is_optimised(), mol.is_conformer(), 
        round(mol.get_optimisation_E(units='hartree'),3)))
    ax = mol.plot_optimisation_E(units='hartree')
    ax.get_figure().set_size_inches(3, 2)
    ax = mol.plot_freq_analysis()
    ax.get_figure().set_size_inches(4, 2)


.. parsed-literal::

    Optimised? True, Conformer? True, Energy = -805.105 a.u.
    


.. image::  images/output_9_1.png



.. image::  images/output_9_2.png


Geometric analysis...

.. code:: python

    print 'Cl optimised polar coords from aromatic ring : ({0}, {1},{2})'.format(
        *[round(i, 2) for i in mol.calc_polar_coords_from_plane(20,3,2,1)])
    ax = mol.plot_opt_trajectory(20, [3,2,1])
    ax.set_title('Cl optimisation path')
    ax.get_figure().set_size_inches(4, 3)


.. parsed-literal::

    Cl optimised polar coords from aromatic ring : (0.11, -116.42,-170.06)
    


.. image::  images/output_11_1.png


Potential Energy Scan analysis of geometric conformers...

.. code:: python

    mol2 = pg.molecule.Molecule(folder_obj=folder, alignto=[3,2,1],
                pes_fname=['CJS_emim_6311_plus_d3_scan.log', 
                           'CJS_emim_6311_plus_d3_scan_bck.log'])   
    ax = mol2.plot_pes_scans([1,4,9,10], rotation=[0,0,90], img_pos='local_maxs', zoom=0.5)
    ax.set_title('Ethyl chain rotational conformer analysis')
    ax.get_figure().set_size_inches(7, 3)



.. image::  images/output_13_0.png


Natural Bond Orbital and Second Order Perturbation Theory analysis...

.. code:: python

    print '+ve charge centre polar coords from aromatic ring: ({0} {1},{2})'.format(
        *[round(i, 2) for i in mol.calc_nbo_charge_center(3, 2, 1)])
    display(mol.show_nbo_charges(represent='ball_stick', axis_length=0.4, 
                                  rotations=[[0,0,90], [-90, 90, 0]]))


.. parsed-literal::

    +ve charge centre polar coords from aromatic ring: (0.02 -51.77,-33.15)
    


.. image::  images/output_15_1.png


.. code:: python

    print 'H inter-bond energy = {} kJmol-1'.format(
            mol.calc_hbond_energy(eunits='kJmol-1', atom_groups=['emim', 'cl']))
    print 'Other inter-bond energy = {} kJmol-1'.format(
        mol.calc_sopt_energy(eunits='kJmol-1', no_hbonds=True, atom_groups=['emim', 'cl']))
    display(mol.show_sopt_bonds(min_energy=1, eunits='kJmol-1',
                                atom_groups=['emim', 'cl'],
                                no_hbonds=True,
                                rotations=[[0, 0, 90]]))
    display(mol.show_hbond_analysis(cutoff_energy=5.,alpha=0.6, 
                                    atom_groups=['emim', 'cl'],
                                    rotations=[[0, 0, 90], [90, 0, 0]]))


.. parsed-literal::

    H inter-bond energy = 111.7128 kJmol-1
    Other inter-bond energy = 11.00392 kJmol-1
    


.. image::  images/output_16_1.png



.. image::  images/output_16_2.png


Multiple Computations Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Multiple computations, for instance of different starting conformations,
can be grouped into an *Analysis* class.

.. code:: python

    analysis = pg.Analysis(folder_obj=folder)
    errors = analysis.add_runs(headers=['Cation', 'Anion', 'Initial'], 
                                   values=[['emim'], ['cl'],
                                           ['B', 'BE', 'BM', 'F', 'FE']],
                init_pattern='*{0}-{1}_{2}_init.com',
                opt_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_opt*unfrz.log',
                freq_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_freq*.log',
                nbo_pattern='*{0}-{1}_{2}_6-311+g-d-p-_gd3bj_pop-nbo-full-*.log',
                alignto=[3,2,1], atom_groups={'emim':range(1,20), 'cl':[20]})
    
    fig, caption = analysis.plot_mol_images(mtype='initial', max_cols=3,
                            info_columns=['Cation', 'Anion', 'Initial'],
                            rotations=[[0,0,90]])
    print caption


.. parsed-literal::

    Figure: (A) emim, cl, B, (B) emim, cl, BE, (C) emim, cl, BM, (D) emim, cl, F, (E) emim, cl, FE
    


.. image::  images/output_19_1.png


The methods mentioned for indivdiual molecules can then be applied to
all or a subset of these computations.

.. code:: python

    analysis.add_mol_property_subset('Opt', 'is_optimised', rows=[2,3])
    analysis.add_mol_property('Energy (au)', 'get_optimisation_E', units='hartree')
    analysis.add_mol_property('Cation chain, $\\psi$', 'calc_dihedral_angle', [1, 4, 9, 10])
    analysis.add_mol_property('Cation Charge', 'calc_nbo_charge', 'emim')
    analysis.add_mol_property('Anion Charge', 'calc_nbo_charge', 'cl')
    analysis.add_mol_property(['Anion-Cation, $r$', 'Anion-Cation, $\\theta$', 'Anion-Cation, $\\phi$'], 
                                   'calc_polar_coords_from_plane', 3, 2, 1, 20)
    analysis.add_mol_property('Anion-Cation h-bond', 'calc_hbond_energy', 
                              eunits='kJmol-1', atom_groups=['emim', 'cl'])
    analysis.get_table(row_index=['Anion', 'Cation', 'Initial'], 
                       column_index=['Cation', 'Anion', 'Anion-Cation'])



.. role:: raw-latex(raw)
    :format: latex html

.. raw:: html

    <script type="text/javascript" src="http://localhost/mathjax/MathJax.js?config=TeX-AMS_HTML"></script>

.. raw:: latex html

    <div style="max-height:1000px;max-width:1500px;overflow:auto;">
    <table border="1" class="dataframe">
      <thead>
        <tr>
          <th></th>
          <th></th>
          <th></th>
          <th colspan="2" halign="left"></th>
          <th colspan="2" halign="left">Cation</th>
          <th>Anion</th>
          <th colspan="4" halign="left">Anion-Cation</th>
        </tr>
        <tr>
          <th></th>
          <th></th>
          <th></th>
          <th>Opt</th>
          <th>Energy (au)</th>
          <th>chain, $\psi$</th>
          <th>Charge</th>
          <th>Charge</th>
          <th>$r$</th>
          <th>$\theta$</th>
          <th>$\phi$</th>
          <th>h-bond</th>
        </tr>
        <tr>
          <th>Anion</th>
          <th>Cation</th>
          <th>Initial</th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
          <th></th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th rowspan="5" valign="top">cl</th>
          <th rowspan="5" valign="top">emim</th>
          <th>B</th>
          <td>NaN</td>
          <td>-805.105</td>
          <td>80.794</td>
          <td>0.888</td>
          <td>-0.888</td>
          <td>0.420</td>
          <td>-123.392</td>
          <td>172.515</td>
          <td>111.713</td>
        </tr>
        <tr>
          <th>BE</th>
          <td>NaN</td>
          <td>-805.105</td>
          <td>80.622</td>
          <td>0.887</td>
          <td>-0.887</td>
          <td>0.420</td>
          <td>-123.449</td>
          <td>172.806</td>
          <td>112.382</td>
        </tr>
        <tr>
          <th>BM</th>
          <td>True</td>
          <td>-805.104</td>
          <td>73.103</td>
          <td>0.874</td>
          <td>-0.874</td>
          <td>0.420</td>
          <td>124.121</td>
          <td>-166.774</td>
          <td>130.624</td>
        </tr>
        <tr>
          <th>F</th>
          <td>True</td>
          <td>-805.118</td>
          <td>147.026</td>
          <td>0.840</td>
          <td>-0.840</td>
          <td>0.420</td>
          <td>10.393</td>
          <td>0.728</td>
          <td>202.004</td>
        </tr>
        <tr>
          <th>FE</th>
          <td>NaN</td>
          <td>-805.117</td>
          <td>85.310</td>
          <td>0.851</td>
          <td>-0.851</td>
          <td>0.417</td>
          <td>-13.254</td>
          <td>-4.873</td>
          <td>177.360</td>
        </tr>
      </tbody>
    </table>
    </div>



**NEW FEATURE:** there is now an option (requiring
`pdflatex <http://www.tug.org/applications/pdftex/>`__ and
`ghostscript <http://www.ghostscript.com/download/gsdnld.html>`__\ +\ `imagemagik <http://www.imagemagick.org/script/binary-releases.php>`__)
to output the tables as a latex formatted image.

.. code:: python

    analysis.get_table(row_index=['Anion', 'Cation', 'Initial'],
                       column_index=['Cation', 'Anion', 'Anion-Cation'],
                       as_image=True, font_size=12)




.. image::  images/output_23_0.png



RadViz is a way of visualizing multi-variate data.

.. code:: python

    ax = analysis.plot_radviz_comparison('Anion', columns=range(4, 10))



.. image::  images/output_25_0.png


The KMeans algorithm clusters data by trying to separate samples into n
groups of equal variance.

.. code:: python

    pg.utils.imgplot_kmean_groups(
        analysis, 'Anion', 'cl', 4, range(4, 10), 
        output=['Initial'], mtype='optimised', 
        rotations=[[0, 0, 90], [-90, 90, 0]],
        axis_length=0.3)



.. image::  images/output_27_0.png


.. parsed-literal::

    Figure: (A) BM
    


.. image::  images/output_27_2.png


.. parsed-literal::

    Figure: (A) FE
    


.. image::  images/output_27_4.png


.. parsed-literal::

    Figure: (A) B, (B) BE
    


.. image::  images/output_27_6.png


.. parsed-literal::

    Figure: (A) F
    

Documentation (MS Word)
~~~~~~~~~~~~~~~~~~~~~~~

After analysing the computations, it would be reasonable to want to
document some of our findings. This can be achieved by outputting
individual figure or table images via the folder object.

.. code:: python

    file_path = folder.save_ipyimg(vdw, 'image_of_molecule')
    Image(file_path)




.. image::  images/output_30_0.png



But you may also want to produce a more full record of your analysis,
and this is where `python-docx <https://python-docx.readthedocs.org>`__
steps in. Building on this package the pygauss MSDocument class can
produce a full document of your analysis.

.. code:: python

    d = pg.MSDocument()
    d.add_heading('A Pygauss Example Assessment', level=1)
    
    d.add_paragraph('We have looked at the following aspects;')
    d.add_list(['geometric conformers', 'electronic structure'])
    
    d.add_heading('Geometric Conformers', level=2)
    fig, caption = analysis.plot_mol_images(max_cols=2, 
                    rotations=[[90,0,0], [0,0,90]], 
                    info_columns=['Anion', 'Cation', 'Initial'])
    d.add_mpl(fig, dpi=96, height=9)
    fig.clear()
    d.add_markdown(caption.replace('Figure:', '**Figure:**'))
    d.add_paragraph()
    df = analysis.get_table(columns=['Anion Charge', 'Cation Charge', 
                                     'Energy (au)'],
                       row_index=['Anion', 'Cation', 'Initial'])
    d.add_dataframe(df, incl_indx=True, style='Medium Shading 1 Accent 1')
    d.add_markdown('**Table:** Analysis of Conformer Charge')
    
    d.save('exmpl_assess.docx')



Which gives us the following:

.. image:: /images/example_docx.png
   :alt: DocX Image

MORE TO COME!!

