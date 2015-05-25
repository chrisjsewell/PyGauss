# -*- coding: utf-8 -*-
"""
Created on Mon May 18 21:01:25 2015

@author: chris sewell
"""
import os, glob
import socket
import errno
import re
import getpass
from io import BytesIO

import paramiko

class Folder(object):
    """ an object intended to act as an entry point to a folder path 
    
    it will act identical whether the folder is local or on a server
    
    """
    def __init__(self, path, 
                 server=None, username=None, passwrd=None):
        """an object intended to act as an entry point to a folder path
        
        Parameters
        ----------
        path : str
            the path to the folder (absolute or relative)
        server : str
            the server name
        username : str
            the username to connect to the server
        passwrd : str
            server password, if not present it will be asked for during initialisation
        
        """
        assert type(path) is str
        self._path = path
        
        if not server:
            self._local = True
        else: 
            self._local = False
            self._server = server
            self._username = username
            #TODO encrypt?
            self._passwrd = passwrd
            
        # set folder and test it exists
        if self._local:
            if not os.path.exists(self._path):
                raise IOError(
                    'the folder path does not exist: {}'.format(self._path))
        else:
            ssh_failed = False
            try:
                ssh = self._connect_ssh(server, username, passwrd)
            except:
                ssh_failed = True
            
            if ssh_failed:
                if not type(self._username) is str:                   
                    self._passwrd = getpass.getuser()
                if not type(self._passwrd) is str:                   
                    self._passwrd = getpass.getpass('Please enter server password: ')

                ssh = self._connect_ssh(server, username, self._passwrd)

            sftp = ssh.open_sftp()
            try:
                sftp.stat(path)
            except IOError, e:
                ssh.close()
                if e.errno == errno.ENOENT:
                    raise IOError("{0} does not exist on server: {1}".format(path, 
                                                                          server))
                else:
                    IOError('error trying to validate folder \n {0}'.format(e))
            
            ssh.close()      

            self._ssh = None 
            self._sftp = None                 

    def get_path(self):
        return self._path
    
    def islocal(self):
        return self._local
    
    def active(self):
        if self._local:
            return True
        elif self._ssh:
            return True
        else:
            return False
        
    def _connect_ssh(self, ssh_server, ssh_username, ssh_passwrd):
        """ connect and verify ssh connection """
            
        ssh = paramiko.SSHClient() 
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        try:
            ssh.connect(ssh_server, username=ssh_username, password=ssh_passwrd)
        except socket.error, e:
            raise IOError(
            'could not connect to the ssh server: \n {0} \n {1}'.format(ssh_server, e))
        except paramiko.ssh_exception.AuthenticationException, e:
            raise IOError(
            'username or password authentication error \n {0}'.format(e))
        except Exception, e:
            raise IOError('error connecting to server: \n {0}'.format(e))

        return ssh
        
    def __enter__(self):
        """ use with statement to open ssh connection once """
        if self._local:
            return self
        ssh = self._connect_ssh(self._server, 
                  self._username, self._passwrd)
        self._ssh = ssh
        self._sftp = ssh.open_sftp()
        return self
        
    def __exit__(self, type, value, traceback):
        """ use with statement to open ssh connection once """
        if self._local:
            return
        try:
            self._ssh.close()
        except:
            pass
        self._ssh = None
        self._sftp = None

    def list_files(self, pattern=None, one_file=False):
        """ list files in folder 

        pattern : str
            a pattern the file must match that can include * wildcards
            
        """
        if self._local:

            if not pattern:
                pattern = '*'
            filepaths = glob.glob(os.path.join(self._path, pattern))
            files = [os.path.basename(f) for f in filepaths]

        else:
            
            if not self._ssh:
                ssh = self._connect_ssh(self._server, 
                          self._username, self._passwrd)
                sftp = ssh.open_sftp()
                files = sftp.listdir(self._path)
                ssh.close()

            else:
                files = self._sftp.listdir(self._path)

            if pattern:
                pattern = "".join(
                [ c if c.isalnum() or c=='*' else "["+c+"]" for c in pattern]
                ).replace('*', '.*')
                files = filter(lambda x: re.match(pattern,x), files)
        
        if not one_file:
            return files
        
        if not files:
            raise IOError(
                'no files of format {0} in path: {1}'.format(pattern, self._path))
        if len(files)>1:
            raise IOError(
            'multiple files found conforming to format {0} in path: {1}'.format(
            pattern, self._path))
            
        if self._local:
            return os.path.basename(files[0])
        else:
            return files[0]
    
    def read_file(self, file_name):
        """ """
        mode='rb'
        
        file_name = self.list_files(file_name, one_file=True)
        
        if self._local:
            return open(os.path.join(self._path, file_name), mode)

        #assume it is a unix server (so '/' is path seperator)
        #otherwise if you use os.path.join on a windows os it will not find
        if not self._path[-1] == '/':
            file_path = self._path + '/' + file_name
        else:
            file_path = self._path + file_name
        
        if not self._sftp:  
            raise IOError('must have an open ssh connection (use `with` statement)')
              
        return self._sftp.file(file_path, mode)

    def write_file(self, file_name, overwrite=False):
        """ """
        mode = 'w'
        
        if not overwrite:
            f = None
            try:
                f = self.list_files(file_name, one_file=True)                
            except:
                pass
            if f:
                raise IOError('the file {0} already exists'.format(file_name))
        
        if self._local:
            return open(os.path.join(self._path, file_name), mode)
        
        #assume it is a unix server (so '/' is path seperator)
        #otherwise if you use os.path.join on a windows os it will not find
        if not self._path[-1] == '/':
            file_path = self._path + '/' + file_name
        else:
            file_path = self._path + file_name

        if not self._sftp:  
            raise IOError('must have an open ssh connection (use `with` statement)')
        
        return self._sftp.file(file_path, mode)

    #TODO write save_mplfig for non-local      
    def save_mplfig(self, fig, fig_name, dpi=256, format='png'):
        """a function for outputing a matplotlib figure to a file
        
        fig : Matplotlib.figure.Figure
            a Matplotlib figure
        fig_name : str
            the desired name of the file
        
        """
        try:
            fig.get_figwidth()
        except AttributeError:
            raise ValueError('the fig is not a Matplotlib figure')
        
        if not os.path.splitext(fig_name)[1]:
            fig_name += os.path.extsep + 'png'
        
        if self._local:
            full_path = os.path.join(self._path, fig_name)
            fig.savefig(full_path, dpi=dpi, 
                        bbox_inches='tight')
            return os.path.abspath(full_path)

        else:
            raise NotImplementedError
    
    #TODO write save_ipyimg for non-local
    def save_ipyimg(self, img, img_name):
        """a function for outputing an IPython Image to a file
        
        img : IPython.display.Image
            an IPyton image
        img_name : str
            the desired name of the file
        
        """
        try:
            data = img.data
        except AttributeError:
            raise ValueError('the img is not an IPython Image')
            
        #_PNG = b'\x89PNG\r\n\x1a\n'
        _JPEG = b'\xff\xd8'
        ext = 'png'
        if data[:2] == _JPEG:
            ext = 'jpg'
        
        if self._local:
            
            full_path = os.path.join(self._path, img_name)+ os.path.extsep + ext
            with open(full_path, "wb") as f:
                f.write(data)
        
            return os.path.abspath(full_path)

        else:
            raise NotImplementedError
        
    #TODO write save_pilimg
    def save_pilimg(self, img, img_name):
        raise NotImplementedError

class NoOutputFolder(Folder):
    """ a folder object which will not output any data """
    def __init__(self, *args, **kwargs):
        super(NoOutputFolder, self).__init__(*args, **kwargs)
        
    def write_file(self, *arg, **kwargs):
        return BytesIO()
    def save_ipyimg(self, *arg, **kwargs):
        return ''
    def save_mplfig(self, *arg, **kwargs):
        return ''
    def save_pilimg(self, *arg, **kwargs):
        return ''
    
    
        
            

