# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

macro(tribits_set_cmake_policy_if_exists  policyName  oldOrNew)
  if (POLICY ${policyName})
    cmake_policy(SET ${policyName} ${oldOrNew})
  endif()
endmacro()

# Define policies for CMake
tribits_set_cmake_policy_if_exists(CMP0003 NEW) # Don't split up full lib paths to linker args
tribits_set_cmake_policy_if_exists(CMP0007 NEW) # Don't ignore empty list items
tribits_set_cmake_policy_if_exists(CMP0053 NEW) # Make var references much faster
tribits_set_cmake_policy_if_exists(CMP0054 NEW) # Avoid quoted strings lookup variables
tribits_set_cmake_policy_if_exists(CMP0057 NEW) # Support if ( ... IN_LIST ... )
tribits_set_cmake_policy_if_exists(CMP0082 NEW) # Install rules follow order install() called in subdirs
tribits_set_cmake_policy_if_exists(CMP0144 NEW) # find_package() use <PACKAGENAME>_ROOT env var
