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
# The value of cmake_version_min below is set to the first available build of
# CMake v2.7 with three pre-built binaries: for Linux, Mac and Windows.
# 20090924 was the very first day that Kitware produced all three binaries
# in the vCVS directory... prior to that, only Windows binaries were available
# nightly. This is the only available date for pre-built binaries on all three
# platforms for v2.7. On 9/25, the CMake-2-8 branch was created, and Kitware
# started producing v2.9 binaries in vCVS.
#

cmake_version_min = "2.7.20090924" # manual_update

cmake_version_release = "2.8.1" # manual_update

cmake_version_rc = "2.8.1" # auto_update v2.8

cmake_version_dev = "2.9.20100317" # auto_update vCVS
