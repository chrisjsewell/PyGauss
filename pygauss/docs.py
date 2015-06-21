# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:52:53 2015

@author: chris sewell
"""
from math import log10, floor
from io import BytesIO

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_ORIENT
from docx.shared import Cm

from pandas import Index, MultiIndex

class MSDocument(object):
    """a class to output a Microsoft Word Document
    
    NB: docx.api.Document can't be directly inherited as it is a function which 
    returns various classes dependent on the *docx* parameter 
    """
    def __init__(self, docx=None):
        """a class to output a Microsoft Word Document

        inherited api details for :py:class:`docx.document.Document` can be 
        found at; https://python-docx.readthedocs.org/en/latest/api/document.html       
        
        Parameters
        ----------
        docx : str or file-like object
            can be either a path to a .docx file (a string) or a file-like object. 
            If docx is missing or None, the built-in default document “template” 
            is loaded.
        
        """
        self._docx = Document(docx=docx)
    
    def __getattr__(self, name):
        """ required to get docx.Document methods """
        return getattr(self._docx, name)
    
    def __dir__(self):
        """ required to have docx.Document methods in ipython tab completion"""
        dirlist = self.__class__.__dict__.keys() + self._docx.__class__.__dict__.keys()
        return sorted(dirlist)           
    
    def add_markdown(self, text='', style='Body Text'):
        """adds a paragraph to the document, allowing for 
        markdown style **bold** and *italic* text
        """
        
        p = self._docx.add_paragraph(style=style)
        if not text:
            return p
            
        bold_split = text.split('**')
        for i, text in enumerate(bold_split):
            if i % 2 == 0:
                italic_split = text.split('*')
                for j, other in enumerate(italic_split):
                    if j % 2 == 0:
                        p.add_run(other)
                    else:
                        p.add_run(other).italic = True
            else:
                p.add_run(text).bold = True
        
        return p
    
    def add_list(self, text_list=[], numbered=False):
        """adds a list """
        if numbered:
            style='List Number'
        else:
            style='List Bullet'
            
        return [self._docx.add_paragraph(tx, style=style) for tx in text_list]
    
    def add_mpl(self, fig, dpi=None, width=None, height=None, pad_inches=0.2):
        """add matplotlib figure to the document, width/height in cm 
        Amount of padding around the figure """
        stream = BytesIO()
        fig.savefig(stream, format='png', dpi=dpi,
                    bbox_inches='tight', pad_inches=pad_inches,
                    transparent=True)
        
        width = Cm(width) if width else None
        height = Cm(height) if height else None
        
        pic = self._docx.add_picture(stream, width=width, height=height)
        self._docx.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
        
        return pic
           
    def add_dataframe(self, df, incl_indx=True, autofit=True, sig_figures=5,
                      style='Medium List 1 Accent 1'):
        """add dataframe as a table
        
        """
        df = df.fillna('-')
        rows, cols = df.shape 

        if type(df.columns) == Index:
            hrows = 1
        elif type(df.columns) == MultiIndex:
            hrows = len(df.columns.levels)
        else:
            raise TypeError('df does not have standard or multi column index')
        
        if incl_indx:
            if type(df.index) == Index:
                icols = 1
            elif type(df.index) == MultiIndex:
                icols = len(df.index.levels)
            else:
                raise TypeError('df does not have standard or multi row index')
        else:
            icols = 0
            
        table = self._docx.add_table(rows=rows+hrows, cols=cols+icols, 
                                     style=style)
        table.alignment = WD_TABLE_ALIGNMENT.CENTER
        table.autofit = autofit
         
        if type(df.columns) == MultiIndex:
            
            merge_cells=[None, None]
            prev_val = None

            for col, vals in enumerate(df.columns.tolist()):
                for hrw, val in enumerate(vals):
                    cell = table.rows[hrw].cells[col+icols]
                    
                    if not hrw==0:
                        p=cell.add_paragraph(str(val))
                        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                    elif not merge_cells[0]:
                        merge_cells[0] = cell
                        prev_val = val
                    elif not merge_cells[1] and not val==prev_val:
                        p=merge_cells[0].add_paragraph(str(prev_val))
                        p.alignment = WD_ALIGN_PARAGRAPH.CENTER                        
                        merge_cells[0] = cell
                        prev_val = val
                    elif not merge_cells[1]:
                        merge_cells[1] = cell
                        prev_val = val
                    elif val == prev_val:
                        merge_cells[1] = cell
                        prev_val = val
                    else:
                        m = merge_cells[0].merge(merge_cells[1])
                        m.add_paragraph(str(prev_val))
                        m.alignment = WD_ALIGN_PARAGRAPH.CENTER
                        merge_cells=[cell, None]
                        prev_val = val
            
            if merge_cells[0] != None and merge_cells[1] != None:
                m = merge_cells[0].merge(merge_cells[1])
                m.add_paragraph(str(prev_val))
                m.alignment = WD_ALIGN_PARAGRAPH.CENTER
            else:
                p=merge_cells[0].add_paragraph(str(prev_val))
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER                        
                                    
        else:
            for col, head in enumerate(df.keys()):
                p=table.rows[hrows-1].cells[col+icols].add_paragraph(str(head))
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER

        if incl_indx and type(df.index) == MultiIndex:
           for col, name in enumerate(df.index.names):
               p=table.rows[hrows-1].cells[col].add_paragraph(str(name))
               p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        
        def rnd(val):
            try:
                if val >= 0.:
                    return round(val, -int(floor(log10(val))) + (sig_figures - 1))
                else:
                    return -round(-val, -int(floor(log10(-val))) + (sig_figures - 1))
            except Exception:
                return val
                
        for row, id_series in enumerate(df.iterrows()):

            if incl_indx and type(df.index) == Index:
                p=table.rows[row+hrows].cells[hrows-1].add_paragraph(str(df.index[row]))
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER                
                
            if incl_indx and type(df.index) == MultiIndex:
                for col, val in enumerate(df.index.tolist()[row]):
                    p=table.rows[row+hrows].cells[col].add_paragraph(str(val))
                    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
                    
            for col, item in enumerate(id_series[1].iteritems()):
                p=table.rows[row+hrows].cells[col+icols].add_paragraph(str(rnd(item[1])))
                p.alignment = WD_ALIGN_PARAGRAPH.CENTER

        return table
    

        
        