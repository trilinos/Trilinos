#!/usr/bin/env python
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-


import re
import subprocess



class GitUtility(object):
    """
    """
    def __init__(self):
        self._version_str = None
        self._version = None


    @property
    def version_str(self):
        if self._version_str is None:
            self._version_str = subprocess.check_output(['git', '--version'])
            self._version_str = self._version_str.decode('utf-8').strip() 
        return self._version_str


    @property
    def version(self):
        if self._version is None:
            matches = re.findall(r"\d+", self.version_str)
            self._version = { "major": int(matches[0]),
                              "minor": int(matches[1]),
                              "patch": int(matches[2]) 
                            }
        return self._version


    def check_minimum_version(self, req_major, req_minor=0):
        """
        Verify the version of git is greater than req_major.req_minor
        """
        if not isinstance(req_major, int):
            raise TypeError("Required parameter 'req_major' must be an integer.")
        if not isinstance(req_minor, int):
            raise TypeError("Required parameter 'req_minor' must be an integer.")

        det_major = self.version["major"]
        det_minor = self.version["minor"]

        str_err  = "Required Git version must be greater or equal to {}.{}\n".format(req_major, req_minor)
        str_err += "but we detected version {}.{}".format(det_major, det_minor)

        if det_major < req_major:
            raise SystemExit(str_err)
        elif det_major == req_major and det_minor < req_minor:
            raise SystemExit(str_err)
        else:
            print("Git minimum version requirement: OK")

        return 0


    def pretty_print(self):
        print("Git Version Detected: {}".format( self.version_str ))



