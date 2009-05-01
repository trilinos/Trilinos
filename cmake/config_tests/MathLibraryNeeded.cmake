
INCLUDE(CheckCSourceCompiles)

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

  FIND_LIBRARY(MATH_LIBRARY m)

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
