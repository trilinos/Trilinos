#!/usr/bin/env python
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
# pylint: disable=wrong-import-position
# pylint: disable=relative-import
"""Fuctions for managing Environment Modules in Sierra python scripts."""
from __future__ import print_function
import sys
sys.dont_write_bytecode = True
import os
import subprocess

from util import which, find_first_binary, find_file_in_list


class Module(object):
    """Class to hold module configuration and provide module commands"""
    def __init__(self, name=None):
        self.paths = self._setup_paths()
        self.command_name = os.environ.get('LMOD_CMD', 'modulecmd')
        self.command = self._setup_command()
        self.init_file = self._module_setup(name)
        if name is None and self.init_file:
            execfile(self.init_file)

    @staticmethod
    def _setup_paths():
        """Set paths to modules"""
        paths = ["/usr/local/modules/default",
                 "/usr/local/Modules/default",
                 "/usr/netpub/modules/default",
                 "/opt/modules/default",
                 "/projects/modules/default",
                 "/usr/share/modules",
                 "/usr/share/Modules"]

        if 'MODULESHOME' in os.environ:
            paths.append(os.path.join(os.environ['MODULESHOME']))
        return paths

    def _setup_command(self):
        """Set command for module"""
        command = which(self.command_name)
        if command:
            self.paths.append(os.path.dirname(command))
        self.command_name = os.path.basename(self.command_name)

        return command

    def _module_setup(self, name=None):
        """Set up modules for the Sierra nightly processes."""
        if not name:
            name = 'python'
        # Look for the modules init files for python.
        # Look for a modulecmd script that we can execute on this machine,
        # and thatis paired with an init/python script.
        modules_init_file = None
        if os.environ.get('LMOD_CMD') and name == 'python':
            self.command = os.environ.get('LMOD_CMD')
        else:
            self.command = os.environ.get('LMOD_CMD',
                                          find_first_binary("modulecmd",
                                                            self.paths, '-V'))

            if self.command:
                # Look for the executable module init file relative to the
                # command binary.
                fname = os.path.join('init', name)
                grandparent_dir = os.path.dirname(
                    os.path.dirname(self.command))
                modules_init_file = which(fname, grandparent_dir)
                if not modules_init_file:
                    modules_init_file = which(fname + '.py', grandparent_dir)
                if not modules_init_file:
                    # If we didn't find it relative to the command binary, look
                    # for it in the standard paths. Again, this is a hack for
                    # Modules/3.1.6
                    modules_init_file = find_file_in_list(fname, self.paths)
                    if not modules_init_file:
                        modules_init_file = find_file_in_list(
                            fname + '.py', self.paths)
                    if not modules_init_file:
                        raise RuntimeError('Unable to find modules '
                                           'init/python file in {}.'.
                                           format(' '.join(self.paths)))
                # We need to set the MODULESHOME environment variable if
                # it is not yet set.
                if modules_init_file and 'MODULESHOME' not in os.environ:
                    os.environ['MODULESHOME'] = os.path.dirname(
                        os.path.dirname(modules_init_file))
        return modules_init_file

    def _echo_module(self, command, *arguments):
        """A function to execute module commands from a python script,
           and return the exit code, stdout, and stderr to the caller."""
        # If the modulecmd executable wasn't found when this script was
        # initialized, or if MODULESHOME is no longer in the environment, look
        # for it in PATH. This is needed because in 3.1.6 a 'module purge'
        # deletes the MODULESHOME variable, so we can't always use it to find
        # modulecmd.
        try:
            cmd = os.environ.get('LMOD_CMD')
            if not cmd:
                cmd = '%s/bin/%s' % (
                    os.environ['MODULESHOME'], self.command_name)
        except Exception:  # pylint: disable=broad-except
            cmd = self.command
        if not cmd or not os.path.exists(cmd):
            cmd = self.command
            if not which(cmd):
                print("Unable to load modules; no {} found.".format(
                    self.command_name), file=sys.stderr)
                return '', '', ''
        cmdline = '%s python %s %s' % (cmd, command, ' '.join(arguments))

        # Get all the output at once.
        subp = subprocess.Popen(cmdline,
                                shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, )
        (stdout, stderr) = subp.communicate()
        errcode = subp.wait()

        return errcode, stdout, stderr

    def module(self, command, *arguments):
        """A function to execute module commands from a python script."""
        # Normally this function returns nothing as its value. However, if the
        # special string 'no-stderr' is found in the arguments list this
        # function will return the lines collected from stderr.
        # See the _echo_module function for more information how this
        # argument is used.

        # Look for the special argument 'no-stderr'. If we find it, flag and
        # remove it. Anything written to stderr by command will not be
        # written to sys.stderr if this flag is passed. It will be
        # returned as the function value, however.
        no_stderr = 'no-stderr' in arguments
        if no_stderr:
            arguments = list(arguments)
            arguments.remove('no-stderr')

        (_, output, stderr) = self._echo_module(command, *arguments)

        # The modules return a string of python commands to be executed on
        # stdout.  Execute those commands now.
        if output:
            exec (output)  # pylint: disable=exec-used

        # Check stderr for anything that looks like an error.
        if ":ERROR:" in stderr:
            raise RuntimeError(stderr)
        elif stderr:
            # Do not print stderr if no-stderr was passed
            if not no_stderr:
                print(stderr, file=sys.stderr)
        return stderr

    def module_list(self):
        """A wrapper around 'module list' that sorts the output."""
        # This function differs from 'module list' in that it will sort the
        # module names reported to ensure a consistent order, then write
        # the resulting sorted list to stderr.
        mnames = self.module('--terse', 'list', 'no-stderr')
        mlist = []
        if mnames:
            mlist = mnames.splitlines()

        # Sort the module names list (from 1 to the end) to maintain a
        # consistent order. We skip line 1 because that is the header line.
        if len(mlist) > 1:
            names = mlist[1:]
            del mlist[1:]
            names.sort()
            mlist.extend(names)
        # Now re-write the sorted names to stderr.
        sys.stderr.writelines('\n'.join(mlist) + '\n')


def module(command, *arguments):
    """Legacy support of module call"""
    return Module().module(command, *arguments)


def module_list():
    """Legacy support of module_list call"""
    Module().module_list()


def get_init_file():
    """Called when script is executed to print name of init file"""
    # If this module is invoked directly, print the modules_init_file found
    # for the passed name.
    # Note that $1 must be a Modules shell init file name.
    if len(sys.argv) != 2:
        raise RuntimeError('A single init file name parameter is required'
                           ' by Modules.py when invoked directly')
    init_file = Module(sys.argv[1]).init_file
    if not init_file:
        raise RuntimeError('Unable to determine init file for {}'
                           ' in Modules.py'.format(sys.argv[1]))
    print(init_file)


if __name__ == '__main__':  # pragma: no cover
    get_init_file()
