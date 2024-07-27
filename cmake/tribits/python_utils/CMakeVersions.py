# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#                      Copyright (c) 2013-2016 NTESS
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

# This file designates the officially supported versions of CMake for Trilinos.
# It is updated periodically as new builds of CMake are available.
#
# The variables defined at the bottom of this file are as follows:
#
# cmake_version_min:
#   Minimum required -- updated manually, synched with cmake_minimum_required
#   Should match Trilinos/CMakeLists.txt cmake_minimum_required
#
# cmake_version_release:
#   Latest official release -- updated manually
#   May be the same as cmake_version_min.
#
# cmake_version_rc:
#   Latest available release candidate for vdir "v3.10" -- detected automatically
#   May be the same as the latest official release if there are no release
#   candidates currently available.
#
# cmake_version_dev:
#   Latest available build for vdir "vCVS" -- detected automatically
#

cmake_version_min = "3.10.0" # manual_update

cmake_version_release = "3.10.0" # manual_update

cmake_version_rc = "3.12.0" # auto_update v2.8

cmake_version_dev = "3.12.0" # auto_update vCVS
