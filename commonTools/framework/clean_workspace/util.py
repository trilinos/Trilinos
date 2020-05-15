#!/usr/bin/env python
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
# pylint: disable=line-too-long
# pylint: disable=wrong-import-position
# pylint: disable=too-many-lines
"""
Useful utility routines

Generic utility routines for use by all.

$Id: util.py,v 1.148 2009/02/09 22:47:53 mhamilt Exp $
"""
from __future__ import print_function
import sys
sys.dont_write_bytecode = True

import os
import types
import subprocess

COMMON_SSH_AUTH_OPTIONS = ["-o", "BatchMode=yes",
                           "-o", "StrictHostKeyChecking=no",
                           "-o", "LogLevel=ERROR"]
OPENSSH_SPECIFIC_AUTH_OPTIONS = ["-o", "GSSAPIDelegateCredentials=yes"]


def _find_files(file_name, dir_list=None):
    """
    A generator to find and return all all instances of
    file_name in the passed dir_list.
    """
    # import pyfind
    if not dir_list:
        dir_list = os.getcwd()
    if not isinstance(dir_list, list):
        dir_list = [dir_list]
    if file_name:
        # Find all the instances of find_name in each subdirectory.
        for subdir in dir_list:
            for root, dirs, files in os.walk(subdir):
                if file_name in files:
                    yield os.path.join(root, file_name)


def _find_binary(binary_name, dir_list=None,
                 version_option='-v',
                 fail_on_error=False,
                 newest=True,
                 return_all=False):
    """
    Find a binary in the passed dir_list that returns a 0 exit
    status when passed the version_option.
    To return the newest matching binary, pass
        newest=True
    To return the first matching binary w/o testing it, pass
        version_option=None
    To return a list of all matching binaries (newest first if
    newest=True,) pass
        return_all=True

    If no matching binary is found, None will be returned.

    This function differs from util.find_file_in_hierarchy in that
    it searches down from the passed directory(ies) instead of up.

    This function differs from util.which in that it uses a passed
    list of directories, rather than using PATH, and it searches
    the entire hierarchy (rooted at each passed directory) for the
    binary rather than just the passed directories.

    This function assumes that the binary may appear more than once
    in subdirectories of each directory in dir_list, so pyfind will
    be used to search for all instances of binary_name in each directory.
    """
    if not binary_name:
        if fail_on_error:
            msg = "No executable name passed."
            sys.exit(97)
        else:
            return None
    if not dir_list:
        dir_list = os.getcwd()
    if not isinstance(dir_list, list):
        dir_list = [dir_list]
    # Find all the instances of binary_name in each subdirectory.
    # Check them, and return the first one which is executable and
    # returns 0 when executed with version_option.
    gen = _find_files(binary_name, dir_list)
    if newest:
        # If the caller wants the newest or all matching binaries, get the
        # list of matching files  from the generator and sort it by
        # modification time.
        file_list = [fname for fname in gen]
        file_list.sort(key=lambda n: os.stat(n).st_mtime, reverse=True)
        # Replace the generator with the sorted list to pass to the loop.
        gen = file_list
    file_list = []
    for fname in gen:
        if os.access(fname, os.X_OK):
            # Test the version only if an option was passed
            if version_option:
                subp = subprocess.Popen(fname + ' %s' % version_option,
                                        shell=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
                (subp_stdout, subp_stderr) = subp.communicate()
                estat = subp.wait()
                if estat != 0:
                    # If this one failed, try the next.
                    continue
            if not return_all:
                # The binary we found passed or we didn't test the version.
                # Either way terminate the generator and return the binary.
                if isinstance(gen, types.GeneratorType) and \
                   hasattr(types.GeneratorType, 'close'):
                    gen.close()
                return fname
            else:
                file_list.append(fname)

    # If all was requested and we found files, return them.
    if file_list and return_all:
        return file_list

    if fail_on_error:
        gen = _find_files(binary_name, dir_list)
        msg = None
        for fname in gen:
            if not msg:
                msg = "Failed to find a %s executable for this system." \
                      " Checked these binaries:" % binary_name
            msg += "\n  %s" % fname
        if not msg:
            msg = "No %s executable found in %s" % (binary_name, dir_list)
        sys.exit(97)

    return None


def find_first_binary(
        binary_name, dir_list=None, version_option='-v', fail_on_error=False):
    """ Get the first binary found; use when dir_list is a list
        rather than a single directory. """
    return _find_binary(binary_name, dir_list,
                        version_option=version_option,
                        fail_on_error=fail_on_error,
                        newest=False)


def find_file_in_list(fname, path_list, tmp_dir=None):
    """ Search for a file in a list of paths. """

    # Test for an absolute path to the file.
    if fname[0] == os.sep and os.path.exists(fname):
        return fname

    # Test for paths relative to those in path_list
    for a_dir in path_list:
        if tmp_dir:
            a_dir = os.path.join(a_dir, tmp_dir)
        tmp_file = os.path.join(a_dir, fname)
        if os.path.exists(tmp_file):
            return tmp_file

    return None

# These directories (other than '/') are automount points on the
# ESPHC lan and various cluster machines. Searching for non-existent
# files in them will cause errors to be logged in the system message file.
# At the request of sysadmin, this will stop searching at these points.
SIERRA_STOP_DIRS = ['/', '/home', '/u0',
                    '/ascldap/users', '/projects', '/sierra']

PROJECT_FILE_NAME = 'SNTools.project'
PROJECT_SEARCH_STOP_FILE_NAME = '.project_search_stop'
# For testing purposes the first element must the the project file;
# all others may be directories or files, but only files can contain a '/'
# (Ie, this should only stop on top-level directories, and files one or
# more levels down.)
PROJECT_STOP_FILES = [PROJECT_FILE_NAME,
                      'Nbtools',
                      'Sierra/Sierra.xml',
                      'sntools/engine']


def which(fname, path=None, find_file_func=None):
    """
    This function will find executable files in PATH.

    Given fname and an optional os.pathsep separated string of paths,
    return the first location of an executable fname within path as
    an absolute pathname.

    This is equivalent to the bash 'type' or csh 'which' command.
    The path argument defaults to the user's PATH.
    """

    # If fname is None or empty string, return None
    if not fname:
        return None

    if not find_file_func:
        def is_executable(fname):
            return os.access(fname, os.X_OK)

        find_file_func = is_executable

    # Do not put this in the function definition above.
    # It will not get executed at the correct time, and path will be wrong.
    if path is None:
        path = os.environ["PATH"]

    if isinstance(path, str):
        path = path.split(os.pathsep)

    if os.path.isabs(fname):
        if find_file_func(fname):
            return os.path.normpath(fname)
        else:
            return None

    if not path:
        path = []
    for idir in path:
        fname_path = os.path.join(idir, fname)
        if find_file_func(fname_path):
            return os.path.normpath(fname_path)
    return None
