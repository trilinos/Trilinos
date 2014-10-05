# NTS: The funny thing is that I need to include this in order to 
# get the link/include directories correct.  Go figure.
INCLUDE(TribitsAddExecutable)

if(TPL_ENABLE_MATLAB)
  INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIRS})

  # Debugging information  
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("MEX_COMPILER    = " ${MEX_COMPILER})
    MESSAGE("MEX_MEXEXT      = " ${MEX_MEXEXT})
    MESSAGE("MEX_EXTENSION   = " ${MEX_EXTENSION})
  ENDIF()
  
  # First, grab the package's own libraries
  SET(LINK_LIBS)
  APPEND_SET(LINK_LIBS ${${PACKAGE_NAME}_LIBRARIES})

  # Third, add test dependent package libraries
  TRIBITS_GATHER_ENABLED_ITEMS(${PACKAGE_NAME} TEST PACKAGES ALL_DEP_PACKAGES)
  TRIBITS_SORT_AND_APPEND_INCLUDE_AND_LINK_DIRS_AND_LIBS("${${PROJECT_NAME}_REVERSE_PACKAGES}"
    "${ALL_DEP_PACKAGES}" "" LINK_LIBS)
  
  # Fourth, add dependent test TPL libraries
  TRIBITS_GATHER_ENABLED_ITEMS(${PACKAGE_NAME} TEST TPLS ALL_TPLS)
  TRIBITS_SORT_AND_APPEND_INCLUDE_AND_LINK_DIRS_AND_LIBS("${${PROJECT_NAME}_REVERSE_TPLS}" "${ALL_TPLS}"
    TPL_ LINK_LIBS)

  # Fifth, add matlab-specific libs
  SET(LINK_LIBS ${LINK_LIBS} "mx" "mex" "mat")
 
  # Last, add last_lib to get extra link options on the link line
  IF (${PROJECT_NAME}_EXTRA_LINK_FLAGS)
    APPEND_SET(LINK_LIBS ${last_lib})
  ENDIF()
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(LINK_LIBS)
  ENDIF()

  # Manually drop in options from the mex script (R2009b) on a linux platform.
  # g++ -O -pthread -shared -Wl,--version-script,/usr/local/matlab/7.9/extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined -o  "mlmex.mexa64"   "mlmex-mlmex.o"  -lm -Wl,-rpath-link,/usr/local/matlab/7.9/bin/glnxa64 -L/usr/local/matlab/7.9/bin/glnxa64 -lmx -lmex -lmat -lm
  SET(MLMEX_OPTS1 "-pthread;-shared;-Wl,--version-script,${MATLAB_ROOT}/extern/lib/${MATLAB_ARCH}/mexFunction.map;-Wl,--no-undefined")
  
  SET(MLMEX_OPTS2 "-Wl,-rpath-link,${MATLAB_ROOT}/bin/${MATLAB_ARCH}")
    
  # Use TARGET_LINK_LIBRARIES and the C++ compiler to link the mlmex.cpp file to the rest of Trilinos & the mex libs.
  # This code is extremely fragile and probably won't work on any system but GNU/Linux with gcc.
  # This is because Cmake will not allow me to *just call the mex linker*, and so I have to do this the hard way.
  #
  LINK_DIRECTORIES(${MATLAB_LIBRARY_DIRS})
  ADD_EXECUTABLE(mlmex.${MEX_EXTENSION} mlmex.cpp)
  SET_TARGET_PROPERTIES(mlmex.${MEX_EXTENSION} PROPERTIES SUFFIX "") # remove .exe extension
  TARGET_LINK_LIBRARIES(mlmex.${MEX_EXTENSION} ${MLMEX_OPTS1} ${LINK_LIBS} ${MLMEX_OPTS2})

  # Copy over the ml.m file from src
  CONFIGURE_FILE("ml.m" "ml.m" COPYONLY)  
ENDIF()
  
