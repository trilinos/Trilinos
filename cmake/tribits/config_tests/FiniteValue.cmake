# @HEADER
# ************************************************************************
#
#            TriBITS: Tribial Build, Integrate, and Test System
#                    Copyright 2013 Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ************************************************************************
# @HEADER

# Checks for the isnan() and isinf() functions needed by Trilinos.
# These functions are not supported reliably across all platforms.
# Even if they are supported, they sometimes don't have a prototype
# defined in a header, making it useless in c++. 

# We check for a predefined version and use that.  If not, then we
# fall back on checks that may or may not work depending on the
# platforms's compliance with IEEE standards.

INCLUDE(CheckCXXSourceRuns)
INCLUDE(CheckCXXSourceCompiles)

#############################################################
# isnan 
#############################################################

# Some machines have isnan() in the global namespace and some put it
# in the std:: namespace.  We will check for both.

SET(SOURCE_GLOBAL_ISNAN
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

CHECK_CXX_SOURCE_COMPILES("${SOURCE_GLOBAL_ISNAN}" 
  FINITE_VALUE_HAVE_GLOBAL_ISNAN)

SET(SOURCE_STD_ISNAN
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

CHECK_CXX_SOURCE_COMPILES("${SOURCE_STD_ISNAN}" 
  FINITE_VALUE_HAVE_STD_ISNAN)

IF (CMAKE_VERBOSE_MAKEFILE)
  IF (NOT FINITE_VALUE_HAVE_GLOBAL_ISNAN AND 
      NOT FINITE_VALUE_HAVE_STD_ISNAN )
    message("****************************************************")
    message("** Warning: Your compiler doesn't support isnan() or")
    message("** std::isnan()")
    message("** We will supply a default checker but it is ")
    message("** *NOT* guaranteed to work on your platform")
    message("** unless your machine is IEEE 748/754 compliant.")
    message("****************************************************")
  ENDIF()
ENDIF()

#############################################################
# isinf
#############################################################

SET(SOURCE_GLOBAL_ISINF
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

CHECK_CXX_SOURCE_COMPILES("${SOURCE_GLOBAL_ISINF}" 
  FINITE_VALUE_HAVE_GLOBAL_ISINF)

SET(SOURCE_STD_ISINF
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

CHECK_CXX_SOURCE_COMPILES("${SOURCE_STD_ISINF}" 
  FINITE_VALUE_HAVE_STD_ISINF)

IF (CMAKE_VERBOSE_MAKEFILE)
  IF (NOT FINITE_VALUE_HAVE_GLOBAL_ISINF AND 
      NOT FINITE_VALUE_HAVE_STD_ISINF )
    message("****************************************************")
    message("** Warning: Your compiler doesn't support isinf() or")
    message("** std::isinf()")
    message("** We will supply a default checker but it is ")
    message("** *NOT* guaranteed to work on your platform")
    message("** unless your machine is IEEE 748/754 compliant.")
    message("****************************************************")
  ENDIF()
ENDIF()

#############################################################
#############################################################
