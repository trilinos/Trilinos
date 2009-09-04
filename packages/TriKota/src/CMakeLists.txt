
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

IF(NOT ${PACKAGE_NAME}_ENABLE_DakotaReconfigure 
  AND EXISTS ${PACKAGE_SOURCE_DIR}/Dakota/dakota_config.h)

EXECUTE_PROCESS(
  COMMAND echo "   TriKota: Dakota already configured, skipping reconfigure."
  WORKING_DIRECTORY ${PACKAGE_SOURCE_DIR}/Dakota
)

ELSE()

EXECUTE_PROCESS(
  COMMAND echo "-->TriKota: Launching Dakota configure script -- takes several minutes."
  WORKING_DIRECTORY ${PACKAGE_SOURCE_DIR}/Dakota
)

EXECUTE_PROCESS(
  WORKING_DIRECTORY ${PACKAGE_SOURCE_DIR}/Dakota
  COMMAND ./configure --without-graphics
    CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} MPICXX=${CMAKE_CXX_COMPILER}
    FC=${CMAKE_Fortran_COMPILER} F77=${CMAKE_Fortran_COMPILER}
    --prefix=${PACKAGE_BINARY_DIR}/Dakota
    --with-teuchos-include=${PACKAGE_SOURCE_DIR}/../teuchos/src
    --with-teuchos-lib=${PACKAGE_BINARY_DIR}/../teuchos/src
    OUTPUT_FILE ${PACKAGE_SOURCE_DIR}/Dakota/configure.out
  )

EXECUTE_PROCESS(
  COMMAND echo "<--TriKota: Finished Dakota configure script."
  WORKING_DIRECTORY ${PACKAGE_SOURCE_DIR}/Dakota
)

ENDIF()

# Customize the build process, by adding a custom make
# target, which makes Dakota (in source tree) and then
# uses Dakota's PREFIX to install it in the build tree.

ADD_CUSTOM_TARGET( TriKota_Dakota_libs ALL
   COMMAND make -j 4
   WORKING_DIRECTORY ${PACKAGE_SOURCE_DIR}/Dakota
  )

ADD_CUSTOM_TARGET( TriKota_Dakota_install ALL
   COMMAND make install
   WORKING_DIRECTORY ${PACKAGE_SOURCE_DIR}/Dakota
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
  ${CMAKE_CURRENT_SOURCE_DIR}/../Dakota
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
    dakota pecos fftw3 lhs evidence surfpack 
    conmin ddace fsudace jega nlpql npsol opt
    psuade newmat ncsuopt quadrature coliny colin
    pebbl utilib 3po nappspack appspack conveyor
    shared cdd cport amplsolver lhs
  )

LINK_DIRECTORIES(${PACKAGE_BINARY_DIR}/Dakota/lib)

# Compiling code that includes Dakota headers requires
# -DHAVE_CONFIG_H.
ADD_DEFINITIONS(-DHAVE_CONFIG_H)

PACKAGE_ADD_LIBRARY(
  trikota
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  IMPORTEDLIBS ${DAKOTA_LIBS}
  )

# Custom "make install" command to copy contents of Dakota lib
FILE(GLOB_RECURSE dakota_lib_files ${PACKAGE_BINARY_DIR}/Dakota/lib/*)
INSTALL(FILES ${dakota_lib_files} DESTINATION lib)

# Custom "make install" command for dakota includes,
# preserving directory structure
SET(DAKINC  ${PACKAGE_BINARY_DIR}/Dakota/include)

# Q: need *hpp.inl  as well?
FILE(GLOB dakota_header_files ${DAKINC}/*.h  ${DAKINC}/*.H  ${DAKINC}/*.hpp)
INSTALL(FILES ${dakota_header_files} DESTINATION include)

FILE(GLOB tpo_header_files ${DAKINC}/3po/*.h)
INSTALL(FILES ${tpo_header_files} DESTINATION include/3po)

# Q: need ampl/*.hd as well?
FILE(GLOB ampl_header_files ${DAKINC}/ampl/*.h)
INSTALL(FILES ${ampl_header_files} DESTINATION include/ampl)

# Q: need boost subdirectories as well?
FILE(GLOB boost_header_files ${DAKINC}/boost/*.hpp)
INSTALL(FILES ${boost_header_files} DESTINATION include/boost)

FILE(GLOB colin_header_files ${DAKINC}/colin/*.h)
INSTALL(FILES ${colin_header_files} DESTINATION include/colin)

FILE(GLOB coliny_header_files ${DAKINC}/coliny/*.h)
INSTALL(FILES ${coliny_header_files} DESTINATION include/coliny)

FILE(GLOB coliny_header_files ${DAKINC}/coliny/*.h)
INSTALL(FILES ${coliny_header_files} DESTINATION include/coliny)

FILE(GLOB pebbl_header_files ${DAKINC}/pebbl/*.h)
INSTALL(FILES ${pebbl_header_files} DESTINATION include/pebbl)

FILE(GLOB pecos_header_files ${DAKINC}/pecos/*.h ${DAKINC}/pecos/*.hpp)
INSTALL(FILES ${pecos_header_files} DESTINATION include/pecos)

FILE(GLOB utilib_header_files ${DAKINC}/utilib/*.h)
INSTALL(FILES ${utilib_header_files} DESTINATION include/utilib)
