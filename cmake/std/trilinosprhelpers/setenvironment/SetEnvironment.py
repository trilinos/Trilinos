#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module handles setting up an environment using modules
and environment variables, etc. by processing a configuration
file (.ini) using the Python module configparser

See: https://docs.python.org/3/library/configparser.html

Todo:
     Add an option to throw an error if any of the actions fail during an
     apply() step.

     Expand the envvar language to allow a couple new things:
         1. envvar-prepend VARNAME: NEW_VALUE
         2. envvar-append  VARNAME: NEW_VALUE
         3. envvar-set-if-missing VARNAME: NEW_VALUE
    As with all the envvar stuff, not being able to have multiple 'keys' with
    the same value could be problematic since we can't do something like
        envvar-append PATH: foo
        envvar-append PATH: bar    # would fail b/c key is already in the dict
        envvar-append PATH: baz    # would fail because ^^^
    One way to fix this might be to allow a third thing in a tuple to uniquify
    the key... so we might do somethign like:
        envvar-append PATH 1: foo
        envvar-append PATH 2: bar
        envvar-append PATH 3: baz
    If we ignore the 3rd tuple from a space-separated key then something like
    this could work.  We might apply this to other commands too.
"""
from __future__ import print_function

import configparser
import os
import pprint
import re
import sys

#from .ModuleHelper import module
from . import ModuleHelper


class SetEnvironment(object):
    """
    Configure environment based on a configuration from a .ini file.

    Args:
        filename (str): The filename of the .ini file to load.
        profile  (str): The profile or <section> in the .ini file to
            process for an action list.

    Attributes:
        config  ()    : The data from the .ini file as loaded by configparser.
        profile (str) : The profile section from the .ini file that is loaded.
        actions (dict): The actions that would be processed when apply() is called.

    """

    def __init__(self, filename, profile):
        self._file        = filename
        self._profile     = profile
        self._actions     = None
        self._config      = None


    @property
    def config(self):
        """
        The config property returns the configuration that is specified
        by the configuration file. If the file has not been loaded yet,
        then we will load it.

        Returns:
            ConfigParser object containing the contents of the configuration
            that is loaded from the .ini file. This does not contain the
            specific profile data from the file.
        """
        if self._config is None:
            self._config = configparser.ConfigParser()

            # Prevent ConfigParser from lowercasing keys.
            self._config.optionxform = str

            try:
                with open(self._file) as ifp:
                    self._config.read_file(ifp)

            except IOError as err:
                msg = "+" + "="*78 + "+\n" + \
                      "|   ERROR: Unable to load configuration file\n" + \
                      "|   - Requested file: {}\n".format(self._file) + \
                      "|   - CWD: {}\n".format(os.getcwd()) + \
                      "+" + "="*78 + "+\n"
                sys.exit(msg)

        return self._config


    @property
    def profile(self):
        """
        Returns the profile.

        Returns:
            String: containing the profile section that was specified in the c'tor.
        """
        return self._profile


    @property
    def actions(self):
        """
        The actions property encodes the module and environment actions that we
        are queueing up to do.

        The value to assign to the 'actions' property. It should conform to this:

            { 'setenv': []
                'unsetenv': []
                'module-op': []
                'module-list': {}
            }

        Returns:
            dict: containing the list of actions that were specified in the .ini
                file and will be executed.
        """
        if self._actions is None:
            self._load_configuration()
        return self._actions


    @actions.setter
    def actions(self, value):
        if not isinstance(value, dict):
            raise TypeError("Invalid type provided.")
        if "setenv" not in value.keys():
            raise KeyError("Configuration is missing the 'setenv' key.")
        if "unsetenv" not in value.keys():
            raise KeyError("Configuration is missing the 'unsetenv' key.")
        if "module-op" not in value.keys():
            raise KeyError("Configuration is missing the 'module-op' key.")
        if "module-list" not in value.keys():
            raise KeyError("Configuration is missing the 'module-list' key.")

        self._actions = value


    def apply(self, throw_on_error=False):
        """
        Apply the settings configured in the `actions` property.

        Args:
            throw_on_error (bool): Throw an error if an instruction fails.
                                   Default: False

        Returns:
            An integer corresponding to success or failure.
            0 = success, anything else indicates failure.

        Raises:
            IndexError: Invalid number of parameters for 'module' command
        """
        status = 0
        print("Module Operations")
        print("+---------------+")
        for m in self.actions["module-op"]:
            op     = m[0]
            params = m[1:]
            cmd_status = 0
            print("[module] module {} {}".format(op, " ".join(params)), end="")
            if   1 == len(params) and op in ["purge"]:
                cmd_status = max(ModuleHelper.module(op), status)
            elif 1 == len(params):
                cmd_status = max(ModuleHelper.module(op, params[0]), status)
            elif 2 == len(params):
                cmd_status = max(ModuleHelper.module(op, params[0], params[1]), status)
            else:
                print(" ERROR")
                msg =  "Invalid number of parameters for `module` command.\n"
                msg += "- Received {} parameters to a `{}` command.".format(len(params), op)
                raise IndexError(msg)

            status = max(cmd_status, status)

            if 0 == cmd_status:
                print(" OK")
            else:
                print(" FAILED")

        print("")
        print("Environment Vars")
        print("+--------------+")

        for envvar in self.actions["setenv"]:
            k = envvar['key']
            v = envvar['value']
            print("[envvar] export {}={}".format(k,v))
            v = self._expand_envvars_in_string(v)
            os.environ[k] = v
            print("[envvar] export {}={} (actual)".format(k,v))
        for k in self.actions["unsetenv"]:
            del os.environ[k]
            print("[envvar] unset {}".format(k))

        if status != 0 and throw_on_error:
            sys.stdout.flush()
            raise Exception("FATAL ERROR in SetEnvironment.apply()")

        return status


    def pretty_print_envvars(self, envvar_filter=None, filtered_keys_only=False):
        """
        Print out a filtered list of environment variables. This routine provides
        a simplified view of the environment variables on a system since sometimes
        just printing out the contents of all envvars on a system can result in very
        cluttered console logs that are difficult to read through.

        `envvar_filter` contains a list of strings that is matched against envvar
        key names, if the envvar key contains one of the strings in envvar_filter
        then we print both the envvar and its value. Otherwise, we print just the
        envvar name.

        If we set filtered_keys_only to True then we ONLY print out envvars that
        matched envvar_filter.

        Arguments:
            envvar_filter (list): a list of keys to print out the value.
                all envvar values are printed if omitted.
                Default: None
            filtered_keys_only (bool): If true, we only display envvars that contain
                a match to an entry in `envvar_filter`.
                If false, we display the keys of all envvars and key+value of envvars
                that matched a filter in the `envvar_filter` list.
                Default: False

        Returns:
            int 0

        Todo:
            This could probably be moved outside of this class since it's not
            really directly relevant to the environment setting stuff, but it
            might also be useful to have here for checking our work.
        """
        if envvar_filter is not None:
            assert isinstance(envvar_filter, list)

        print("+" + "="*68 + "+")
        print("|   P R I N T   E N V I R O N M E N T   V A R S")
        print("+" + "="*68 + "+")
        print("--- ")
        for k,v in os.environ.items():
            matched_key = False
            if envvar_filter is not None:
                for f in envvar_filter:
                    if f in k:
                        matched_key = True
                        break
            else:
                filtered_keys_only = False

            if filtered_keys_only == False or matched_key:
                print("--- {}".format(k), end="")
                if envvar_filter is None:
                    print(" = {}".format(v), end="")
                elif matched_key:
                    print(" = {}".format(v), end="")
                print("")
        print("")
        return 0


    def pretty_print(self):
        """
        Pretty print the list of actions.

        Raises:
            IndexError: Invalid number of parameters for 'module' command
        """
        print("")
        print("Module Operations")
        print("+---------------+")
        for m in self.actions["module-op"]:
            op     = m[0]
            params = m[1:]
            print("[module] module {} {}".format(op, " ".join(params)) )

            if len(params) not in [1,2] and op not in ["purge"]:
                msg =  "Invalid number of parameters for `module` command.\n"
                msg += "- Received {} parameters to a `{}` command.".format(len(params), op)
                raise IndexError(msg)

        print("")
        print("Environment Vars")
        print("+--------------+")
        for envvar in self.actions["setenv"]:
            # Note: we can't print the expansion here because the source
            #       envvar might not exist yet.
            k = envvar['key']
            v = envvar['value']
            print("[envvar] export {}={}".format(k,v))
        for k in self.actions["unsetenv"]:
            print("[envvar] unset {}".format(k))

        return 0


    def _expand_envvars_in_string(self, string_in):
        """
        Take an input string that may contain environment variables in the style
        of BASH shell environment vars (i.e., "${foobar}") and replace them with
        the actual environment variables.

        Returns:
            A string that contains the contents of any `${ENVVAR}` entries expanded
            inline into the string.

        Raises:
            KeyError: Required environment variable does not exist.
        """
        regexp = re.compile(r"(\$\{(\S*)\})")
        string_out = string_in
        for m in re.finditer(regexp, string_out):
            #v = m.group(1)  # The full ENVVAR sequence: ${VARNAME}
            s = m.group(2)  # Just the ENVVAR itself: VARNAME
            if(s in os.environ.keys()):
                string_out = re.sub(regexp, os.environ[s], string_in)
            else:
                msg = "Required environment variable `{}` does not exist.".format(s)
                raise KeyError(msg)
        return string_out


    def _load_configuration(self):
        """
        Loads the environment data from the provided `profile` and saves the
        information to the `actions` property.

        This is a recursive operation (recurses on the `use <profile>:` entries
        in the configuration file). We use a DFS recursion strategy with back-edge
        checking to prevent cyclic recursion.

        Returns:
            A dict containing the contents of self.actions.
        """
        self.actions = { "setenv": [],         # Envvars that we'll set
                         "unsetenv": [],       # Envvars to explicitly unset (after setting)
                         "module-op": [],      # module operations (use, load, unload, etc.)
                         "module-list": {}     # Keys = modules already set to 'load'
                        }

        self.actions = self._load_configuration_r(self.profile, actions=self.actions)

        return self.actions


    def _load_configuration_r(self, config_section, actions=None, processed_secs=None):
        """
        Recursive handler for _load_configuration.
        Recursion is called when the operation is a 'use' operation in the key.

        Arguments:
            config_section (str): The section (profile) name that is getting processed by
                this recursive call. This matches the [section name] entries from the .ini
                file.
            actions (dict): The actions structure that is constructed by this routine.
                This is the ultimate output of the routine.
            processed_secs (dict): A dictionary where the keys indicate sections that
                have already been processed. It is used to prevent cycles in the recursion.

        Raises:
            KeyError: If `config_section` section is not contained in `config`
            Exception: If a section isn't loaded for a reason other than a KeyError.
        """
        # Load the section and throw if it's not found.
        sec = None
        try:
            sec = self.config[config_section]
        except KeyError as err:
            message = "ERROR: No section named '{}' was found in the configuration file.".format(config_section)
            #raise Exception(message) from err  # (This is ok in Python3 but not Python2)
            raise KeyError(message)

        # Verify that we got a section, if it's None then raise an exception.
        if sec is None:                                                                             # pragma: no cover
            raise Exception("ERROR: Unable to load section '{}'".format(config_section) +           # pragma: no cover
                            " from configuration for unknown reason.")                              # pragma: no cover

        # If processed_secs is the default parameter then we seed it with {}
        # Note: We cannot set the default processed_secs to {} or we'll have
        #       errors.
        if processed_secs is None:
            processed_secs = {}

        # Exit if we've already processed this section (recursion breaker)
        if config_section in processed_secs.keys():
            print("WARNING: Cyclic dependencies encountered when loading {} ... {}".format(self.profile, config_section))
            return actions

        processed_secs[config_section] = True

        for op2,value in sec.items():
            #print("> k=`{}` v=`{}`".format(k,v))
            value  = str(value.strip('"'))
            oplist = op2.split()
            op1    = str(oplist[0])

            if(2==len(oplist)):
                op2 = str(oplist[1])

            if "use" == op1:
                new_profile = op2
                self._load_configuration_r(new_profile, actions, processed_secs)

            elif "module-purge" == op1:
                actions["module-op"].append(["purge", ''])

            elif "module-use" == op1:
                actions["module-op"].append(["use", value])

            elif "module-load" == op1:
                if 0 == len(value) or value is None:
                    actions["module-op"].append(["load", op2])
                else:
                    actions["module-list"][op2] = True
                    actions["module-op"].append(["load", "{}/{}".format(op2,value)])

            elif "module-remove" == op1:
                if op2 in actions["module-list"].keys():
                    del actions["module-list"][op2]

                regexp = re.compile(r"^{}/.*".format(op2))
                new_modules = list(filter(lambda x : not regexp.search(x[1]), actions["module-op"]))
                actions["module-op"] = new_modules

            elif "module-unload" == op1:
                actions["module-op"].append(["unload", op2])
                # If the module is in the list of keys we have, we remove it from
                # module-list since the net effect of a load/unload would be to not
                # have the module.
                # Note: This does NOT remove it from the list of operations... just from
                # module-list entry, which is used for bookkeeping.
                if(op2 in actions['module-list'].keys()):
                    del actions['module-list'][op2]

            elif "module-swap" == op1:
                actions["module-op"].append(["swap", op2, value])

            elif "setenv" == op1:
                actions["setenv"].append( {'key': op2.upper(), 'value': value} )

            elif "unsetenv" == op1:
                actions["unsetenv"].append(op2.upper())

        return actions


