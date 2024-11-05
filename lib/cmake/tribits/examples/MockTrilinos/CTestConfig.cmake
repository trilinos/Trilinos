# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

include(SetDefaultAndFromEnv)

# Should never be used to submit anywhere!  This is just for the unit test
# mock project!

# Must match what is in CDash project 'Trilinos'
set(CTEST_NIGHTLY_START_TIME "00:00:00 UTC")

# Set actual CTest/CDash settings

if (NOT DEFINED CTEST_DROP_METHOD)
  set_default_and_from_env(CTEST_DROP_METHOD "http")
endif()

if (CTEST_DROP_METHOD STREQUAL "http")
  set_default_and_from_env(CTEST_DROP_SITE "dummy.com")
  set_default_and_from_env(CTEST_PROJECT_NAME "MockProjectName")
  set_default_and_from_env(CTEST_DROP_LOCATION "/cdash/submit.php?project=MockProjectName")
  set_default_and_from_env(CTEST_TRIGGER_SITE "")
  set_default_and_from_env(CTEST_DROP_SITE_CDASH TRUE)
endif()
