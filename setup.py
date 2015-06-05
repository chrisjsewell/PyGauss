# -*- coding: utf-8 -*-
"""
Created on Thu May 07 14:34:45 2015

@author: chris
"""
#!/usr/bin/env python

from distutils.core import setup

import os
def readme(file, git_path, img_folder):
    if os.path.exists(file):
        descript = open(file).read()
        descript = descript.replace('image:: ', 
                                    'image:: {0}/raw/master/{1}/'.format(git_path, img_folder))
        return descript
    return ''

import re
def version(path):
    verstrline = open(path, "rt").read()
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    mo = re.search(VSRE, verstrline, re.M)
    if mo:
        verstr = mo.group(1)
    else:
        raise RuntimeError("Unable to find version string in %s." % (path,))
    return verstr

setup(name='pygauss',
      version=version("pygauss/_version.py"),
      description='Python API for Analysis of Gaussian Quantum Chemical Compuations',
      keywords = "chemistry gaussian dft",
      long_description=readme('setup_README.rst',
                              'https://github.com/chrisjsewell/PyGauss',
                              'readme_images'),
      author='Chris Sewell',
      author_email='chrisj_sewell@hotmail.com',
      url='https://github.com/chrisjsewell/PyGauss/wiki',
      license = "GPL3",
      platforms = ["Any."],
      packages=['pygauss', 
                'pygauss.cclib_patch','pygauss.cclib_patch.parser', 
                'pygauss.chemlab_patch',
                'pygauss.chemlab_patch.graphics', 
                'pygauss.chemlab_patch.graphics.renderers',
                'pygauss.chemlab_patch.io','pygauss.chemlab_patch.io.handlers',
                'pygauss.chemview_patch',
                'pygauss.test_data',
                'pygauss.tests'],
      package_data={'': ['*.rst', '*.txt'],
                    'pygauss.test_data': ['*.log', '*.com']},
      install_requires=[
                          "numpy>=1.9", #.2",
                          "scipy>=0.15", #.1",
                          "matplotlib>=1.4", #.3",
                          "pandas>=0.15", #.2,
                          "ipython>=3",
                          "scikit-learn>=0.15",
                          "paramiko",
                        #"scikit-image",
                          "pillow",
                        #"numexpr",
                          "cclib",
                          "chemview",
                          "pyopengl==3.0.2",
                          "chemlab"
                       ],               
     )
