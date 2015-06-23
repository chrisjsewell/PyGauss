# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:52:53 2015

@author: chris sewell
"""
from math import log10, floor
from io import BytesIO
import re

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.section import WD_ORIENT
from docx.shared import Cm

from pandas import DataFrame, Index, MultiIndex

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
        """ required to get :py:class:`docx.document.Document` methods """
        return getattr(self._docx, name)
    
    def __dir__(self):
        """ required to have :py:class:`docx.document.Document` methods in 
        :py:mod:`IPython` tab completion"""
        dirlist = self.__class__.__dict__.keys() + self._docx.__class__.__dict__.keys()
        return sorted(dirlist)           
    

    _MARKUPS = {
            'italic':('*','*'),
            'bold':('**', '**'),
            'subscript':('_{', '}'),
            'superscript':('^{', '}'),
            'strike':('~~','~~'),
            'math': ('$', '$')
            }

    def _get_markup(self, para):
        """get markup """

        df = DataFrame(self._MARKUPS, index=['Enter', 'Exit']).T
        df['In']=False
        
        sects=[]
        place=0
        while place > -1:
            place = -1
            markup = None
            estr = None
            for mark, enter in df[df.In==False].Enter.iterkv():
                find = para.find(enter)
                if find > -1 and (find<=place or place==-1):
                    if find == place and len(enter) < len(estr):
                        continue
                    place = find
                    markup = mark
                    estr = enter
            for mark, exit in df[df.In==True].Exit.iterkv():
                find = para.find(exit)
                if find > -1 and (find<=place or place==-1):
                    if find == place and len(exit) < len(estr):
                        continue
                    place = find
                    markup = mark
                    estr = exit
        
            if place > -1:
                sects.append([para[:place], df[df.In==True].index.tolist()])
                df.set_value(markup, 'In', not df.get_value(markup, 'In'))
                para = para[place+len(estr):]

        if df.In.any():
            raise ValueError(
                'the markup does not exit from;\n{}'.format(df[df.In==True]))
            
        sects.append([para, []])
                         
        return sects

    def add_markdown(self, text='', style='Body Text', para=None):
        r"""adds a paragraph to the document, allowing for
        paragraph/font styling akin to a stripped down version of
        markdown text::
        
            - bullet list
            1. numbered list           
            **bold** 
            *italic* 
            _{subscript} 
            ^{superscript} 
            ~~strikethrough~~ 
            $mathML$
        
        Parameters
        ----------
        text : str
            the text to add
        style : str
            the style to apply (overriden if a header or list)
        para : docx.text.paragraph.Paragraph
            a pre-existing paragraph to add the text to
        
        Returns
        -------
        para : docx.text.paragraph.Paragraph
            a paragraph added to the document
            
        """
        list_pattern = re.compile('^[-+]\s')
        number_pattern = re.compile('^\d+[.]\s')
        head_pattern = re.compile('^#+\s')
        level=0
        
        if re.match(list_pattern, text):
            style = 'List Bullet'
            text = text[len(re.findall(list_pattern, text)[0]):]
        elif re.match(number_pattern, text):
            style = 'List Number'
            text = text[len(re.findall(number_pattern, text)[0]):]
        elif re.match(head_pattern, text):
            level = len(re.findall(head_pattern, text)[0]) - 1
            text = text[level+1:]
 
        if not para:
            if level:
                para = self._docx.add_heading(level=level)
            else:
                para = self._docx.add_paragraph(style=style)
        if not text:
            return para
        
        sects = self._get_markup(text)
        for txt, markups in sects:
            run = para.add_run(txt)
            font = run.font
            for markup in markups:
                setattr(font, markup, True)

        return para


    def _split_special_paras(self, text):
        """split text into paras if a header or list,
        denominated by; # heading, - bullet or 1. numbered
        """
        patterns = ['[-+]', '\d+[.]', '#+']

        for pattern in patterns:
            if re.match(re.compile('^{}\s'.format(pattern)), text):
                starts = re.findall(re.compile('\n\s*{}\s'.format(pattern)), '\n'+text)
                texts = re.split(re.compile('\n\s*{}\s'.format(pattern)), '\n'+text)
                return [s[1:]+t for s, t in zip(starts, texts[1:])]
        
        return [text]
        

    def add_docstring(self, docstring, style='Body Text',
                      markdown=True):
        """adds a doctring to the document
            
        this function will split text into paragraphs 
        (denominated by a separating blank line)
        remove new-line characters and add to document, allowing for 
        markdown style text designated in 
        :py:func:`pygauss.docs.MSDocument.add_markdown`
        
        Parameters
        ----------
        text : str
            the text to add
        style : str
            the style to apply for each paragraph
        markdown : bool
            whether to apply markdown to the text
        
        Returns
        -------
        paras : docx.text.paragraph.Paragraph
            a list of paragraphs added to the document        

        """
        docx_paras = []
        para_pattern = re.compile('\n[\s]*\n')
        paras = re.split(para_pattern, docstring)

        # remove initial linespace if present
        if paras[0][:1] == '\n':
            paras[0] = paras[0][1:]

        for para in paras:
            if markdown:
                para = para.strip()
                for p in self._split_special_paras(para):
                    p = p.replace('\n', ' ').strip()
                    docx_paras.append(self.add_markdown(p, style=style))
            else:
                para = para.replace('\n', ' ').strip()
                docx_paras.append(self._docx.add_paragraph(para, style=style))
    
        return docx_paras
        
    def add_list(self, text_list=[], numbered=False):
        """adds a list """
        if numbered:
            style='List Number'
        else:
            style='List Bullet'
            
        return [self._docx.add_paragraph(tx, style=style) for tx in text_list]
    
    def add_mpl(self, fig, dpi=None, width=None, height=None, pad_inches=0.2):
        """add matplotlib figure to the document 
        
        Parameters
        ----------
        fig : matplotlib.figure.Figure
            a matplotlib figure
        dpi : int
            Dots per inch
        width : float
            width of image in document
        height : float
            width of image in document
        pad_inches : float
            amount of padding around the figure

        Returns
        -------
        pic : docx.shape.InlineShape
            an inline picture added to the document        

        """
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
        """add dataframe as a table to the document
        
        Parameters
        ----------
        df : pandas.DataFrame
            a pandas dataframe
        incl_indx : bool
            include dataframes index in table
        autofit : bool
            autfit table in document 
        sig_figures : int
            number of significant figures for numbers in table
        style : str
            MS Word table style

        Returns
        -------
        pic : docx.table.Table
            a table added to the document        
        
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
    

        
        
