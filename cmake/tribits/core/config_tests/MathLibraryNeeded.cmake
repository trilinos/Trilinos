# @HEADER
# *****************************************************************************
#            TriBITS: Tribal Build, Integrate, and Test System
#
# Copyright 2013-2016 NTESS and the TriBITS contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
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
