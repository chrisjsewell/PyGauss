# -*- coding: utf-8 -*-
"""
Created on Thu May 07 14:34:45 2015

@author: chris
"""
#!/usr/bin/env python

from distutils.core import setup

setup(name='PyGauss',
      version='0.1',
      description='PYthon GAUSSian DFT output analysis',
      author='Chris Sewell',
      author_email='chrisj_sewell@hotmail.com',
      url='https://github.com/chrisjsewell/PyGauss',
      platforms = ["Any."],
      packages=['pygauss', 
                'pygauss.cclib_patch', 
                'pygauss.chemlab_patch', 
                'pygauss.chemview_patch',
                'pygauss.test_data'],
                
     )