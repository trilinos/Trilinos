INCLUDE(TrilinosCreateClientTemplateHeaders)

# NTS: The funny thing is that I need to include this in order to
# get the link/include directories correct.  Go figure.
SET(HEADERS_EXEC "")
SET(SOURCES_EXEC "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

IF(TPL_ENABLE_MATLAB)
  TRIBITS_INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIRS})

  # Must glob the binary dir last to get all of the auto-generated headers
  TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})

  # Force the code into R2017b compatibility mode
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMX_COMPAT_64 -fPIC -DMATLAB_MEXCMD_RELEASE=R2017b")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMX_COMPAT_64 -fPIC -DMATLAB_MEXCMD_RELEASE=R2017b")
  # Debugging information
  MESSAGE("MEX_COMPILER    = " ${MEX_COMPILER})
  MESSAGE("MEX_MEXEXT      = " ${MEX_MEXEXT})
  MESSAGE("MEX_EXTENSION   = " ${MEX_EXTENSION})
  MESSAGE("MATLAB_MEX_DIR   = " ${MATLAB_MEX_DIR})

  APPEND_SET(HEADERS_EXEC muemex.h)
  APPEND_SET(SOURCES_EXEC muemex.cpp)

  # Manually drop in options from the mex script (R2009b) on a linux platform.
  # g++ -O -pthread -shared -Wl,--version-script,/usr/local/matlab/7.9/extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined -o  "mlmex.mexa64"   "mlmex-mlmex.o"  -lm -Wl,-rpath-link,/usr/local/matlab/7.9/bin/glnxa64 -L/usr/local/matlab/7.9/bin/glnxa64 -lmx -lmex -lmat -lm
  IF(NOT APPLE)
    SET(MLMEX_OPTS1 "-pthread;-shared;-Wl,--version-script,${MATLAB_ROOT}/extern/lib/${MATLAB_ARCH}/mexFunction.map;-Wl,--no-undefined")
    SET(MLMEX_OPTS2 "-Wl,-rpath-link,${MATLAB_ROOT}/bin/${MATLAB_ARCH}")
  ELSE()
    SET(MLMEX_OPTS1 "-shared")
    SET(MLMEX_OPTS2 "")
  ENDIF()

  # Use TARGET_LINK_LIBRARIES and the C++ compiler to link the mlmex.cpp file to the rest of Trilinos & the mex libs.
  # This code is extremely fragile and probably won't work on any system but GNU/Linux with gcc.
  # This is because Cmake will not allow me to *just call the mex linker*, and so I have to do this the hard way.
  #

  LINK_DIRECTORIES(${MATLAB_LIBRARY_DIRS})
  TRIBITS_ADD_EXECUTABLE(
    muemex.${MEX_EXTENSION}
    SOURCES ${SOURCES_EXEC}
    NOEXEPREFIX
    NOEXESUFFIX
  )
  SET_TARGET_PROPERTIES(muemex.${MEX_EXTENSION} PROPERTIES SUFFIX "") # remove .exe extension
  TARGET_LINK_LIBRARIES(
    muemex.${MEX_EXTENSION}
    PUBLIC ${MLMEX_OPTS1}
           ${LINK_LIBS}
           ${MLMEX_OPTS2}
           ${${PROJECT_NAME}_EXTRA_LINK_FLAGS}
  )

  # Copy these files over from src tree
  TRIBITS_COPY_FILES_TO_BINARY_DIR(Matlab_Output_cp SOURCE_FILES muelu.m setup.m)

  CONFIGURE_FILE(matlab ${CMAKE_CURRENT_BINARY_DIR}/matlab @ONLY)

ENDIF()
