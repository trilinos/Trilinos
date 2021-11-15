# @HEADER
# ************************************************************************
#
#            TriBITS: Tribal Build, Integrate, and Test System
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


include(CheckCSourceCompiles)

if (NOT MATH_LIBRARY_IS_SUPPLIED AND NOT MATH_LIBRARY_IS_SET)

  set(CMAKE_REQUIRED_LIBRARIES ${${PROJECT_NAME}_EXTRA_LINK_FLAGS})

  check_c_source_compiles(
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

  set(CMAKE_REQUIRED_LIBRARIES)

  if (NOT MATH_LIBRARY_IS_SUPPLIED)

    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message(STATUS "Searching for -lm ...")
    endif()

    set(MATH_LIBRARY NOTFOUND)
    find_library(MATH_LIBRARY m)

    if (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      message(STATUS "MATH_LIBRARY = ${MATH_LIBRARY}")
    endif()

    if (MATH_LIBRARY)
      if (NOT MATH_LIBRARY_IS_SET)
        message(STATUS "Appending math library ${MATH_LIBRARY} to link line ...")
        set(${PROJECT_NAME}_EXTRA_LINK_FLAGS ${${PROJECT_NAME}_EXTRA_LINK_FLAGS} ${MATH_LIBRARY}
          CACHE STRING ""  FORCE)
        global_set(MATH_LIBRARY_IS_SET ON)
        # NOTE: Only do this once and not over and over or you will relink
        # everything after each configure!
      endif()
    else()
      message(SEND_ERROR
        "Error, the math library for C programs could not be found!"
        )
    endif()

  endif()

endif()
