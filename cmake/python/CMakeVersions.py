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
#   Latest available release candidate for vdir "v2.8" -- detected automatically
#   May be the same as the latest official release if there are no release
#   candidates currently available.
#
# cmake_version_dev:
#   Latest available build for vdir "vCVS" -- detected automatically
#

cmake_version_min = "2.8.0" # manual_update

cmake_version_release = "2.8.4" # manual_update

cmake_version_rc = "2.8.5" # auto_update v2.8

cmake_version_dev = "2.8.4.20110707-g0eecf" # auto_update vCVS
