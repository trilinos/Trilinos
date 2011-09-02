ADD_DEFINITIONS(-DBuild64 -DADDC_)

IF (${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcray-pointer -fdefault-real-8 -fdefault-integer-8")
ELSE()
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -i8")
ENDIF()

