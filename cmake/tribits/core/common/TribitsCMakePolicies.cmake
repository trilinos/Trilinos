# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

# Define policies for CMake
# It is assumed that the project has already called CMAKE_MINIMUM_REQUIRED.
cmake_policy(SET CMP0003 NEW) # Don't split up full lib paths to linker args
cmake_policy(SET CMP0007 NEW) # Don't ignore empty list items
cmake_policy(SET CMP0053 NEW) # Make var references much faster
cmake_policy(SET CMP0054 NEW) # Avoid quoted strings lookup variables
cmake_policy(SET CMP0057 NEW) # Support if ( ... IN_LIST ... )
cmake_policy(SET CMP0082 NEW) # Install rules follow order install() called in subdirs
