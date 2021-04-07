#!/usr/bin/env python
"""
Module helper for *nix systems

This module contains a helper for using the modules subsystem on
*nix systems.  It will attempt to load the env_modules_python
module and if that does not exist then we generate our own call
to the `module()` function.
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

    # Try to import the LMOD version of the module() function.
    # See: https://github.com/TACC/Lmod/blob/master/init/env_modules_python.py.in
    import env_modules_python
    print("NOTICE> [ModuleHelper.py] Using the lmod based `env_modules_python` module handler.")


    def module(command, *arguments):
        """
        `module` wrapper for the LMOD based modules function.

        Args:
            command (str): The `module` command that we're executing. i.e., `load`, `unload`, `swap`, etc.
            *arguments   : Variable length argument list.

        Returns:
            int: status indicating success or failure. 0 = success, nonzero for failure

            For now, because the LMOD `module()` command returns nothing (`NoneType`)
            and provides no way to determine success or failure, we will only return 0.

        Raises:
            Any exception thrown by env_modules_python.module() will get caught and passed along.

        """
        output = 0

        try:

            status = env_modules_python.module(command, *arguments)

        except BaseException as error:
            print("")
            print("An ERROR occurred during execution of module command")
            print("")
            raise error

        return output



except ImportError:
    print("NOTICE> [ModuleHelper.py] Using the modulecmd based environment modules handler.")


    def module(command, *arguments):
        """
        Function that enables operations on environment modules in
        the system.

        Args:
            command (str): The `module` command that we're executing. i.e., `load`, `unload`, `swap`, etc.
            *arguments   : Variable length argument list.

        Returns:
            int: status indicating success or failure.  0 = success, nonzero for failure.

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

        numArgs = len(arguments)

        cmd = [ modulecmd, "python", command ]

        if (numArgs == 1):
            cmd += arguments[0].split()
        else:
            cmd += list(arguments)

        # Execute the `modules` command (i.e., $ module <op> <module name(s)>)
        # but we don't actually set the environment. If successful, the STDOUT will
        # contain a sequence of Python commands that we can later execute to set up
        # the proper environment for the module operation
        proc = subprocess.Popen( cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE )
        (output,stderr) = proc.communicate()
        errcode = proc.returncode

        # Convert the bytes into UTF-8 strings
        output = output.decode()
        stderr = stderr.decode()

        stderr_ok = True
        if "ERROR:" in stderr:
            print("")
            print("An error occurred in modulecmd:")
            print("")
            stderr_ok = False
            errcode = 1

        if stderr_ok and errcode != 0:
            print("")
            print("Failed to execute the module command: {}".format(" ".join(cmd[2:])))
            print("- Returned {} exit status.".format(errcode))

        if stderr_ok and errcode == 0:
            try:
                # This is where we _actually_ execute the module command body.
                exec(output)

                # Check for _mlstatus = True/False (set by some versions of modulecmd)
                if "_mlstatus" in locals():
                    if locals()["_mlstatus"] == False:
                        print("")
                        print("modulecmd set _mlstatus == False, command failed")
                        print("")
                        errcode = 1

            except BaseException as error:
                print("")
                print("An ERROR occurred during execution of module commands")
                print("")
                raise error

        if errcode != 0:
            print("")
            print("[module output start]\n{}\n[module output end]".format(output))
            print("[module stderr start]\n{}\n[module stderr end]".format(stderr))
            print("")
            sys.stdout.flush()

        # Uncomment this if we want to throw an error rather than exit with nonzero code
        #if errcode != 0:
        #    raise OSError("Failed to execute module command: {}.".format(" ".join(args)))

        if errcode is None:
            raise TypeError("ERROR: the errorcode can not be `None`")

        return errcode

