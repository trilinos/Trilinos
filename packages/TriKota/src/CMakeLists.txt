
INCLUDE(PackageLibraryMacros)

#
IF (NOT EXISTS ${PACKAGE_SOURCE_DIR}/Dakota/configure)
  MESSAGE(FATAL_ERROR "\nTriKota Fatal Error: Dakota needs to 
 be untarred into the Trikota Directory before running cmake!")
ENDIF()

#
# Configure Dakota libraries
#   Skip configuring if dakota_config.h already exists,
#   unless -DTriKota_ENABLE_DakotaReconfigure:BOOL=ON
#
SET(${PACKAGE_NAME}_ENABLE_DakotaReconfigure OFF CACHE BOOL
    "Whether to force reconfiguration of Dakota")

FILE(MAKE_DIRECTORY ${PACKAGE_BINARY_DIR}/Dakota)

IF(NOT ${PACKAGE_NAME}_ENABLE_DakotaReconfigure 
  AND EXISTS ${PACKAGE_BINARY_DIR}/Dakota/dakota_config.h)

EXECUTE_PROCESS(
  COMMAND echo "   TriKota: Dakota already configured, skipping reconfigure."
  WORKING_DIRECTORY ${PACKAGE_BINARY_DIR}/Dakota
)

ELSE()

EXECUTE_PROCESS(
  COMMAND echo "-->TriKota: Launching Dakota configure script -- takes several minutes."
  WORKING_DIRECTORY ${PACKAGE_BINARY_DIR}/Dakota
)

EXECUTE_PROCESS(
  WORKING_DIRECTORY ${PACKAGE_BINARY_DIR}/Dakota
  COMMAND ${PACKAGE_SOURCE_DIR}/Dakota/configure --without-graphics
    CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} MPICXX=${CMAKE_CXX_COMPILER}
    FC=${CMAKE_Fortran_COMPILER} F77=${CMAKE_Fortran_COMPILER}
    --prefix=${PACKAGE_BINARY_DIR}/Dakota/install
    --with-teuchos-include=${PACKAGE_SOURCE_DIR}/../teuchos/src
    --with-teuchos-lib=${PACKAGE_BINARY_DIR}/../teuchos/src
    --with-boost=${PACKAGE_SOURCE_DIR}/Dakota/packages/boost
    --without-dot
    OUTPUT_FILE ${PACKAGE_BINARY_DIR}/Dakota/configure.out
  )

EXECUTE_PROCESS(
  COMMAND echo "<--TriKota: Finished Dakota configure script."
  WORKING_DIRECTORY ${PACKAGE_BINARY_DIR}/Dakota
)

ENDIF()

# Customize the build process, by adding a custom make
# target, which makes Dakota (in source tree) and then
# uses Dakota's PREFIX to install it in the build tree.
#   Reverting to "make" from "make -j 4" which might
#   have been causing problems by being launched within
#   another make -j 4 at the Trilinos level

ADD_CUSTOM_TARGET( TriKota_Dakota_libs ALL
   COMMAND make -j4
   WORKING_DIRECTORY ${PACKAGE_BINARY_DIR}/Dakota
  )

ADD_CUSTOM_TARGET( TriKota_Dakota_install ALL
   COMMAND make install
   WORKING_DIRECTORY ${PACKAGE_BINARY_DIR}/Dakota
  )

## Make must finish before make install, for ampl
ADD_DEPENDENCIES(TriKota_Dakota_install  TriKota_Dakota_libs)

#
# Build the TriKota library
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../Dakota/src
  ${CMAKE_CURRENT_BINARY_DIR}/../Dakota
  )

APPEND_SET(HEADERS
    TriKota_ConfigDefs.hpp
    TriKota_DirectApplicInterface.hpp
    TriKota_Driver.hpp
    TriKota_Version.hpp
  ) 

APPEND_SET(SOURCES
    TriKota_DirectApplicInterface.cpp
    TriKota_Driver.cpp
    TriKota_Version.cpp
  )

SET(DAKOTA_LIBS
    dakota pecos fftw3 lhs surfpack 
    conmin ddace fsudace jega nlpql npsol opt
    psuade newmat ncsuopt quadrature coliny colin
    pebbl utilib 3po nappspack appspack conveyor
    shared cdd cport amplsolver lhs
  )

LINK_DIRECTORIES(${PACKAGE_BINARY_DIR}/Dakota/install/lib)

# Compiling code that includes Dakota headers requires
# -DHAVE_CONFIG_H.
ADD_DEFINITIONS(-DHAVE_CONFIG_H)

PACKAGE_ADD_LIBRARY(
  trikota
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  IMPORTEDLIBS ${DAKOTA_LIBS}
  )
ADD_DEPENDENCIES(trikota TriKota_Dakota_install)

# Custom "make install" command to copy contents of Dakota lib&include
INSTALL(DIRECTORY ${PACKAGE_BINARY_DIR}/Dakota/install/lib/ DESTINATION lib)
INSTALL(DIRECTORY ${PACKAGE_BINARY_DIR}/Dakota/install/include/ DESTINATION include)
