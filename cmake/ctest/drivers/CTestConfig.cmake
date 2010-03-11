#=============================================================================
# CMake - Cross Platform Makefile Generator
# Copyright 2000-2009 Kitware, Inc., Insight Software Consortium
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
set(CTEST_PROJECT_NAME "Trilinos_Driver")
set(CTEST_NIGHTLY_START_TIME "21:00:00 EDT")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "trilinos.sandia.gov")
set(CTEST_DROP_LOCATION "/cdash/submit.php?project=TrilinosDriver")
set(CTEST_DROP_SITE_CDASH TRUE)
set(CTEST_CDASH_VERSION "1.5")
set(CTEST_CDASH_QUERY_VERSION TRUE)
