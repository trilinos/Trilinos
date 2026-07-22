# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

# This file designates the officially supported versions of CMake for Trilinos.
# It is updated periodically as new builds of CMake are available.

from CMakeVersions import *


# All CMake binaries are downloaded from some path under this base url:
#
cmake_baseurl = "http://www.cmake.org/files"


# These are the directories under the base url where we look for new builds.
# Each entry in cmake_vdirs will be scanned for the most-recent available build
# when download-cmake.py is run (unless detection is skipped)
#
# You could also add v2.4 and v2.6 to this list. Installers for the big three
# platforms have been available since v2.4. Left out here because Trilinos
# requires at least a CMake v2.7...
#
cmake_vdirs = (\
  "v2.8",\
  "vCVS",\
  )


# Minimum required -- updated manually, synched with cmake_minimum_required
#
  #
  # version_min should match Trilinos/CMakeLists.txt cmake_minimum_required
  #
  # This one is the very first day of Kitware producing all three binaries on a
  # nightly basis in the vCVS directory... prior to that, only Windows
  # binaries were available nightly. This is the only available date for
  # pre-built binaries on all three platforms for v2.7. On 9/25, the CMake-2-8
  # branch was created, and we started producing v2.9 binaries in vCVS.

cmake_min_binaries = (\
  ('linux2', cmake_baseurl + "/vCVS/cmake-" + cmake_version_min + "-Linux-i386.tar.gz"),\
  ('darwin', cmake_baseurl + "/vCVS/cmake-" + cmake_version_min + "-Darwin-universal.tar.gz"),\
  ('win32', cmake_baseurl + "/vCVS/cmake-" + cmake_version_min + "-win32-x86.zip"),\
  )


# Latest official release -- updated manually
#

cmake_release_binaries = (\
  ('linux2', cmake_baseurl + "/v2.8/cmake-" + cmake_version_release + "-Linux-i386.tar.gz"),\
  ('darwin', cmake_baseurl + "/v2.8/cmake-" + cmake_version_release + "-Darwin-universal.tar.gz"),\
  ('win32', cmake_baseurl + "/v2.8/cmake-" + cmake_version_release + "-win32-x86.zip"),\
  )


# Latest release candidate for vdir "v2.8" -- detected automatically
#

cmake_rc_binaries = (\
  ('linux2', cmake_baseurl + "/v2.8/cmake-" + cmake_version_rc + "-Linux-i386.tar.gz"),\
  ('darwin', cmake_baseurl + "/v2.8/cmake-" + cmake_version_rc + "-Darwin-universal.tar.gz"),\
  ('win32', cmake_baseurl + "/v2.8/cmake-" + cmake_version_rc + "-win32-x86.zip"),\
  )


# Latest verified build for vdir "vCVS" -- detected automatically
#

cmake_dev_binaries = (\
  ('linux2', cmake_baseurl + "/vCVS/cmake-" + cmake_version_dev + "-Linux-i386.tar.gz"),\
  ('darwin', cmake_baseurl + "/vCVS/cmake-" + cmake_version_dev + "-Darwin-universal.tar.gz"),\
  ('win32', cmake_baseurl + "/vCVS/cmake-" + cmake_version_dev + "-win32-x86.zip"),\
  )
