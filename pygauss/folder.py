# -*- coding: utf-8 -*-
"""
Created on Mon May 18 21:01:25 2015

@author: chris sewell
"""
import os, glob
import socket
import errno
import re

import paramiko

from IPython.core.display import clear_output

class Folder:
    """ an object intended to act as an entry point to a folder path 
    
    it will act identical whether the folder is local or on a server
    
    """
    def __init__(self, path, 
                 server=None, username=None, passwrd=None):
        """
        
        path : str
            the path to the folder (absolute or relative)
        server : str
            the server name
        username : str
            the username to connect to the server
        passwrd : str
            the password to connect to the server
        
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
                if not type(self._passwrd) is str:
                    self._passwrd = raw_input('Please enter server password:')
                    try:
                        clear_output()    
                    except:
                        pass
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
    
    def open_file(self, file_name, mode='rb'):
        """ """
        file_name = self.list_files(file_name, one_file=True)
        
        if self._local:
            return open(os.path.join(self._path, file_name), mode)

        #assume it is a linux server (so '/' is path seperator)
        #otherwise if you use os.path.join on a windows os it will not find
        if not self._path[-1] == '/':
            file_path = self._path + '/' + file_name
        else:
            file_path = self._path + file_name
        
        if not self._ssh:  
            raise IOError('must have an open ssh connection (use with statement)')
              
        return self._sftp.file(file_path, mode)
            
            
            

