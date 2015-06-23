# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 04:57:18 2015

@author: chris
"""
from nose import tools as that
from nose_parameterized import parameterized

import pygauss as pg
import matplotlib.pyplot as plt
import pandas as pd


class Test_MSDocuments(object):
    def setUp(self):
        self.doc = pg.MSDocument()
    @parameterized(['no markdown', 
    '**just bold**', '*just italic*', 
    '**bold** first', 'last **bold**', 'middle **bold** middle',
    '*italic* first', 'last *italic*', 'middle *italic* middle',
    '*italic* **bold**', '**bold** *italic*',
    'a *full* phrase of **mixed bold** and *italic*'])
    def test_add_markdown(self, phrase):
        para = self.doc.add_markdown(phrase)
        that.assert_equal(para.text, phrase.replace('*', ''))

    @that.raises(Exception)
    def test_add_markdown_fails(self):
        self.doc.add_markdown('**open bold')
    
    def test_add_docstring(self):
        paras = self.doc.add_docstring("""
        # Heading 1
        ## Heading 2
        
        **paragraph** 1
        
        *paragraph* 2
        
        some_{subscript}
        
        some^{superscript}
        
        - bullet 1
        - bullet 2
        
        1. numbered 1
        2. numbered 2
        """)

        that.assert_equal([p.text for p in paras], 
                          ['Heading 1', 'Heading 2',
                          'paragraph 1', 'paragraph 2',
                          'somesubscript', 'somesuperscript',
                          'bullet 1', 'bullet 2',
                          'numbered 1', 'numbered 2'])        
        
    def test_add_list(self):
        plist = self.doc.add_list(['one', 'two', 'three'])
        that.assert_equal([p.text for p in plist], ['one', 'two', 'three'])

    def test_add_mpl(self):
        fig, ax = plt.subplots()
        self.doc.add_mpl(fig, width=1, height=1)
        
    def test_add_dataframe(self):

        df = pd.DataFrame({('Test1', 'Energy'): {('B', 'cl'): 'emim',
          ('BE', 'cl'): 'emim',
          ('BM', 'cl'): 'emim',
          ('F', 'cl'): 'emim',
          ('FE', 'cl'): 'emim'},
         ('Test2', r'math $\theta$'): {('B', 'cl'): 123456.12345,
          ('BE', 'cl'): -12345.12345,
          ('BM', 'cl'): 1,
          ('F', 'cl'): 0.2345678,
          ('FE', 'cl'): -3.3},
         ('Test2', 'Boolean'): {('B', 'cl'): True,
          ('BE', 'cl'): True,
          ('BM', 'cl'): True,
          ('F', 'cl'): False,
          ('FE', 'cl'): True}})
         
        tbl = self.doc.add_dataframe(df, incl_indx=True, sig_figures=5)
        
        that.assert_equal(tbl.rows[0].cells[2].text.strip(), 'Test1')
        that.assert_equal(tbl.rows[1].cells[2].text.strip(), 'Energy')
        that.assert_equal(float(tbl.rows[2].cells[4].text), 123460)
        that.assert_equal(float(tbl.rows[3].cells[4].text), -12345)
        that.assert_equal(float(tbl.rows[4].cells[4].text), 1)
        that.assert_equal(float(tbl.rows[5].cells[4].text), 0.23457)
        that.assert_equal(float(tbl.rows[6].cells[4].text), -3.3)
        that.assert_equal(tbl.rows[5].cells[3].text.strip(), 'False')
        that.assert_equal(tbl.rows[6].cells[3].text.strip(), 'True')
         