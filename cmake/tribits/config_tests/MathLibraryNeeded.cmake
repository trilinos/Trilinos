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


INCLUDE(CheckCSourceCompiles)

IF (NOT MATH_LIBRARY_IS_SUPPLIED AND NOT MATH_LIBRARY_IS_SET)

  SET(CMAKE_REQUIRED_LIBRARIES ${${PROJECT_NAME}_EXTRA_LINK_FLAGS})
  
  CHECK_C_SOURCE_COMPILES(
    "
#include <math.h>
int main()
{
  double val1 = sqrt(1.0);
  double val2 = log10(2.0);
  double val3 = log(2.0);
  return 0;
} 
    "
    MATH_LIBRARY_IS_SUPPLIED
    )
  
  SET(CMAKE_REQUIRED_LIBRARIES)
  
  IF (NOT MATH_LIBRARY_IS_SUPPLIED)

    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE(STATUS "Searching for -lm ...")
    ENDIF()
  
    SET(MATH_LIBRARY NOTFOUND)
    FIND_LIBRARY(MATH_LIBRARY m)

    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE(STATUS "MATH_LIBRARY = ${MATH_LIBRARY}")
    ENDIF()

    IF (MATH_LIBRARY)
      IF (NOT MATH_LIBRARY_IS_SET)
        MESSAGE(STATUS "Appending math library ${MATH_LIBRARY} to link line ...")
        SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS ${${PROJECT_NAME}_EXTRA_LINK_FLAGS} ${MATH_LIBRARY}
          CACHE STRING ""  FORCE)
        GLOBAL_SET(MATH_LIBRARY_IS_SET ON)
        # NOTE: Only do this once and not over and over or you will relink
        # everything after each configure!
      ENDIF()
    ELSE()
      MESSAGE(SEND_ERROR
        "Error, the math library for C programs could not be found!"
        )
    ENDIF()
  
  ENDIF()

ENDIF()
