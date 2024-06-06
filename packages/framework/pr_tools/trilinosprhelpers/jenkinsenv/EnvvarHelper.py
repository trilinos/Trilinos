#!/usr/bin/env python3
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
"""
"""
import os



class EnvvarHelper(object):
    """
    Envvar helper class.  This can be used as a base class for specific environment
    variables for set environments (i.e., Jenkins sets some standard envvars for its
    jobs).

    This class implements functions that are useful fetching or creating envvars.
    """

    def __init__(self):
        """
        Constructor.
        """
        pass


    def get_envvar_str(self, envvar_name, error_if_missing=False):
        """
        Get the value of an environment variable if it exists and return
        it as a string. If the envvar does not exist, return None.

        Args:
            envvar_name (str): The environment variable name.
            error_if_missing (bool): If True then throw a KeyError if the envvar does not exist.
                                     If False, we return None if the envvar is missing. Default: False

        Returns:
            The value of the envvar as a string if it exists.

            If error_if_missing is false, we return `None` if the envvar is missing
            or we throw a KeyError if it's missing.

        Throws:
            KeyError: If error_is_missing is True and the envvar does not exist.
        """
        assert isinstance(envvar_name, str)
        assert isinstance(error_if_missing, bool)

        output = None
        if envvar_name in os.environ:
            output = str( os.environ[envvar_name] )
        elif error_if_missing:
            raise KeyError("ERROR: Missing required envvar '{}'".format(envvar_name))
        return output


    def get_or_create_if_missing(self, envvar_name, default_value=""):
        """
        Attempt to retrieve an envvar but if it does not exist then set it with
        a default value.

        Args:
            envvar_name (str)  : The name of the envvar we're searching for.
            default_value (str): The default value to set if the envvar doesn't exist.

        Returns:
            str value of the envvar that we either set or read.
        """
        output = str(default_value)
        try:
            output = self.get_envvar_str(envvar_name, error_if_missing=True)
        except KeyError:
            os.environ[envvar_name] = str(default_value)
        return output
