# -*- coding: utf-8 -*-
"""
Created on Thu May 07 14:34:45 2015

@author: chris
"""
#!/usr/bin/env python

from distutils.core import setup

import os
def readme(file, git_img_path):
    if os.path.exists(file):
        descript = open(file).read()
        descript = descript.replace('image:: ', 
                                    'image:: {0}raw/master/'.format(git_img_path))
        return descript
    return ''

setup(name='pygauss',
      version='0.1.6',
      description='PYthon GAUSSian DFT output analysis',
      long_description=readme('setup_README.rst',
                              'https://github.com/chrisjsewell/PyGauss/readme'),
      author='Chris Sewell',
      author_email='chrisj_sewell@hotmail.com',
      url='https://github.com/chrisjsewell/PyGauss',
      platforms = ["Any."],
      packages=['pygauss', 
                'pygauss.cclib_patch','pygauss.cclib_patch.parser', 
                'pygauss.chemlab_patch',
                'pygauss.chemlab_patch.graphics', 
                'pygauss.chemlab_patch.graphics.renderers',
                'pygauss.chemlab_patch.io','pygauss.chemlab_patch.io.handlers',
                'pygauss.chemview_patch',
                'pygauss.test_data'],
      package_data={'': ['*.rst', '*.txt'],
                    'pygauss.test_data': ['*.log', '*.com']},
      install_requires=[
                          "numpy>=1.9.2",
                          "scipy>=0.15.1",
                          "pil>=1.1.7",
                          "matplotlib>=1.4.3",
                          "pandas>=0.15.2",
                          "ipython>=3.0.0",
                          "scikit-learn>=0.15.2, <0.16"
                       ],               
     )