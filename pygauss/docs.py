# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:52:53 2015

@author: chris sewell
"""
import docx
#from docx import Document
#from docx.shared import Cm

class MSDoc(docx.Document):
    """a class to output a Microsoft Word Document
    """
    def __init__(self, docx=None):
        """a class to output a Microsoft Word Document

        inherited api details for :py:class:`docx.api.Document` can be found at;
        https://python-docx.readthedocs.org/en/latest/api/document.html       
        
        Parameters
        ----------
        docx : str or file-like object
            can be either a path to a .docx file (a string) or a file-like object. 
            If docx is missing or None, the built-in default document “template” 
            is loaded.
        
        """
        super(MSDoc, self).__init__(docx=docx)
    
    

        
        