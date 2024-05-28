#!/usr/bin/env python3
# -*- mode: python; py-indent-offset: 4; py-continuation-offset: 4 -*-
"""
System information helpers are contained in this file.

This class attemts to load the psutil module first but if
it does not exist then we try and work around it by falling
back to pulling information from the /proc/meminfo file.

Note:
    `/proc/meminfo` will exist on *nix systems but OSX or Windows
    systems will not have this file. `psutil` is the best way
    but sometimes that package isn't available on build systems.
"""

import multiprocessing

try:
    # Preferred: we can get our information from psutil.
    import psutil                             # pragma: no cover
    have_psutil_vm = True                     # pragma: no cover
except:                                       # pragma: no cover
    # Fallback and load /proc/meminfo (only works on *nix systems)
    have_psutil_vm = False                    # pragma: no cover

    # Check that /proc/meminfo actually exists. If not, fail early.
    import os                                 # pragma: no cover
    if not os.path.exists("/proc/meminfo"):   # pragma: no cover
        msg_error = "Import psutil.virtual_memory failed. /proc/meminfo not found.\n" + \
                    "Testing cannot proceed because we can't determine system information.\n" + \
                    "Try running '$ pip install psutil'"
        raise IOError(msg_error)




class SysInfo(object):
    """
    This class handles getting system information for Trilinos Pull Request
    jobs. The main thing we need to determine about the system we're running
    a PR test on is how many cores we can use in our parallel builds.

    This is generally due to the memory high-water mark required for some of
    the modules. For example, say MueLu requires 5 GB / core to effectively
    compile. We should select the # of cores we need in a parallel make to
    ensure we don't violate the high-water mark.

    One constraint we have on this class is that there are two ways to get the
    memory avilable on a system.
        1. The preferred way is to use the ``psutil`` module in Python.
        2. If ``psutil`` isn't available, we attempt to use a fall-back
           by reading in the ``/proc/meminfo`` file, which should exist
           on *nix systems.
        3. If ``psutil`` and ``/proc/meminfo`` aren't found then we will
           fail during ``import`` time by raising ``IOError``.

    Attributes:
        have_psutil: Returns true if the psutil module is available
            on the system.
        meminfo: Returns a dictionary containing information about the
            memory available on the system.
    """
    def __init__(self):
        self._meminfo = None
        self.have_psutil = have_psutil_vm


    @property
    def have_psutil(self):
        """
        If the psutil module is available, this will be True.

        Returns:
            bool: True if 'import psutil' is successful, False otherwise.
        """
        return self._have_psutil


    @have_psutil.setter
    def have_psutil(self, value):
        """
        Setter for the 'have_psutil' property.

        Args:
            value (bool): A boolean value that indicates whether or not the
                          algorithm will use psutil or not. Note: even if
                          psutil is available, setting this explicitly to False
                          prior to the first call to 'meminfo' will force the
                          alternative method (/proc/meminfo).

        Returns:
            bool: The new value of 'have_psutil'

        Raises:
            TypeError: Raises `TypeError` if the value isn't a boolean.
        """
        if not isinstance(value, bool):
            raise TypeError("unable to convert value to a bool")

        self._have_psutil = value
        return self._have_psutil


    @property
    def meminfo(self):
        """
        Returns a dict containing information about the memory available on the system.

        On first run we generate the value and cache it for later use.

        Returns:
            dict: { 'mem_kb': <system memory in kb>, 'mem_gb': <system memory in gb> }
        """
        if self._meminfo is None:
            self._get_meminfo()
        return self._meminfo


    def compute_num_usable_cores(self, req_mem_gb_per_core=3.0, max_cores_allowed=32):
        """
        Compute the number of usable cores given memory and size
        constraints for the problem.

        Steps:
            1) compute max number of CPUs by cpu constraints: min(sys cores, max_cores_allowed)
            2) compute max cores by mem constraints: sys_mem_GB / req_mem_gb_per_core
            3) return min of (1) and (2)

        Args:
            req_mem_gb_per_core (float): The minimum memory (GB) per core required. Default=3.0
            max_cores_allowed   (int)  : The maximum number of cores to use. Default=32

        Returns:
            int: Returns the minimum of (1) and (2)
        """
        if max_cores_allowed < 1:
            max_cores_allowed = 1

        output = 1

        n_cpu = multiprocessing.cpu_count()
        mem_G = self.meminfo["mem_gb"]

        max_cpu_by_param_constraint = min(n_cpu, max_cores_allowed)
        max_cpu_by_mem_constraint   = int(mem_G / req_mem_gb_per_core)

        output = min(max_cpu_by_mem_constraint, max_cpu_by_param_constraint)

        return output


    def _get_psutil_vm_total(self):
        '''
        Testing has issues when psutil doesn't exist, putting the call
        in this function makes it easier to just mock out this call
        rather than trying to mock psutil itself.
        '''
        return psutil.virtual_memory().total        # pragma: no cover


    def _get_meminfo(self):
        """
        Get Memory Information

        Attempts to load psutil to extract the information. If that fails then we will
        attempt to read this information from /proc/meminfo. Note: backup method will only
        work on *nix systems.

        Returns: dictionary = ["mem_gb": <Gigabytes of RAM>, "mem_kb": <kilobytes of RAM>]
        """
        mem_kb = None

        if self.have_psutil:
            #mem_total = psutil.virtual_memory().total
            mem_total = self._get_psutil_vm_total()
            mem_kb = mem_total / 1024.0
        else:
            with open('/proc/meminfo') as f_ptr:
                meminfo = dict((i.split()[0].rstrip(':'), int(i.split()[1])) for i in f_ptr.readlines())
            mem_kb = meminfo['MemTotal']

        self._meminfo = {}
        self._meminfo["mem_kb"] = float(mem_kb)
        self._meminfo["mem_gb"] = mem_kb / (1024*1024.0)

        return self._meminfo
