# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

# Checks for the isnan() and isinf() functions needed by Trilinos.
# These functions are not supported reliably across all platforms.
# Even if they are supported, they sometimes don't have a prototype
# defined in a header, making it useless in c++.

# We check for a predefined version and use that.  If not, then we
# fall back on checks that may or may not work depending on the
# platforms's compliance with IEEE standards.

include(CheckCXXSourceRuns)
include(CheckCXXSourceCompiles)

#############################################################
# isnan
#############################################################

# Some machines have isnan() in the global namespace and some put it
# in the std:: namespace.  We will check for both.

set(SOURCE_GLOBAL_ISNAN
  "
#include <cmath>
int main()
{
  double x = 1.0;
  isnan(x);
  return 0;
}
  "
  )

check_cxx_source_compiles("${SOURCE_GLOBAL_ISNAN}"
  FINITE_VALUE_HAVE_GLOBAL_ISNAN)

set(SOURCE_STD_ISNAN
  "
#include <cmath>
int main()
{
  double x = 1.0;
  std::isnan(x);
  return 0;
}
  "
  )

check_cxx_source_compiles("${SOURCE_STD_ISNAN}"
  FINITE_VALUE_HAVE_STD_ISNAN)

if (CMAKE_VERBOSE_MAKEFILE)
  if (NOT FINITE_VALUE_HAVE_GLOBAL_ISNAN AND
      NOT FINITE_VALUE_HAVE_STD_ISNAN )
    message("****************************************************")
    message("** NOTE: Your compiler doesn't support isnan() or")
    message("** std::isnan()")
    message("** We will supply a default checker but it is ")
    message("** *NOT* guaranteed to work on your platform")
    message("** unless your machine is IEEE 748/754 compliant.")
    message("****************************************************")
  endif()
endif()

#############################################################
# isinf
#############################################################

set(SOURCE_GLOBAL_ISINF
  "
#include <cmath>
int main()
{
  double x = 1.0;
  isinf(x);
  return 0;
}
  "
  )

check_cxx_source_compiles("${SOURCE_GLOBAL_ISINF}"
  FINITE_VALUE_HAVE_GLOBAL_ISINF)

set(SOURCE_STD_ISINF
  "
#include <cmath>
int main()
{
  double x = 1.0;
  std::isinf(x);
  return 0;
}
  "
  )

check_cxx_source_compiles("${SOURCE_STD_ISINF}"
  FINITE_VALUE_HAVE_STD_ISINF)

if (CMAKE_VERBOSE_MAKEFILE)
  if (NOT FINITE_VALUE_HAVE_GLOBAL_ISINF AND
      NOT FINITE_VALUE_HAVE_STD_ISINF )
    message("****************************************************")
    message("** NOTE: Your compiler doesn't support isinf() or")
    message("** std::isinf()")
    message("** We will supply a default checker but it is ")
    message("** *NOT* guaranteed to work on your platform")
    message("** unless your machine is IEEE 748/754 compliant.")
    message("****************************************************")
  endif()
endif()

#############################################################
#############################################################
