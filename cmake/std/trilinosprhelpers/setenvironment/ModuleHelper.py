#!/usr/bin/env python
"""
Module helper for *nix systems

This module contains a helper for using the modules subsystem on
*nix systems.  It will attempt to load the env_modules_python
module and if that does not exist then we generate our own call
to the `module()` function.

Todo:
    Check to find what package provides env_modules_python and document.
"""
from __future__ import print_function

import os
import subprocess
import sys

if "MODULESHOME" in os.environ.keys():                                             # pragma: no cover
    sys.path.insert(1, os.path.join(os.environ['MODULESHOME'], 'init'))            # pragma: no cover
else:                                                                              # pragma: no cover
    print("WARNING: The environment variable 'MODULESHOME' was not found.")        # pragma: no cover
    print("         ModuleHelper may not be able to locate modules.")              # pragma: no cover

try:

    from env_modules_python import module

except ImportError:
    # If importing module from env_modules_python fails, we roll our own
    # version of that function.

    def module(*args):
        """
        Function that enables operations on environment modules in
        the system.

        Args:


        Raises:
            FileNotFoundError: This is thrown if `modulecmd` is not found.

        Todo:
            Update documentation for this function to list args and how its called.
        """
        try:
            import distutils.spawn
            modulecmd = distutils.spawn.find_executable("modulecmd")
            if modulecmd is None:
                raise FileNotFoundError("Unable to find modulecmd")          # pragma: no cover
        except:
            modulecmd = "/usr/bin/modulecmd"

        if type(args[0]) == type([]):
            args = args[0]
        else:
            args = list(args)

        cmd = [modulecmd, "python"] + args

        # Execute the command by loading the instructions from
        # the modulefile and run them.
        proc = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
        (output,stderr) = proc.communicate()
        errcode = proc.returncode

        if errcode:
            print("Failed to execute the module command: {}.".format(" ".join(args)))
            print("- Returned {} exit status.".format(errcode))
        else:
            exec(output)

        # Convert the bytes into UTF-8 strings
        output = output.decode('utf-8')
        stderr = stderr.decode('utf-8')

        # Check for success... module tends to return 0 regardless of outcome
        # so we'll have a look at the stderr for indications that there was a
        # problem.
        if "ERROR:" in stderr:
            errcode = 1
        elif "_mlstatus = False" in stderr:
            errcode = 1

        if errcode:
            print("module output> {}".format(output))
            print("module stderr> {}".format(stderr))

        # Uncomment this if we want to throw an error rather than exit with nonzero code
        #if errcode != 0:
        #    raise OSError("Failed to execute module command: {}.".format(" ".join(args)))

        return errcode


