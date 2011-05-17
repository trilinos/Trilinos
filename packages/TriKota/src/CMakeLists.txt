
INCLUDE(PackageLibraryMacros)

IF(NOT TriKota_ENABLE_DakotaCMake)

#
IF (NOT EXISTS ${PACKAGE_SOURCE_DIR}/Dakota/configure)
  MESSAGE(FATAL_ERROR "\nTriKota Fatal Error: Dakota needs to 
 be untarred into the Trikota Directory before running cmake!")
ENDIF()

# Dakota requires on Boost. Check if the user set a
# path first, otherwise point into Dakota's boost.
IF (NOT EXISTS ${Boost_INCLUDE_DIRS}/boost/type.hpp)

#  Explain what is going on with Boost paths: (too verbose)
##    MESSAGE("   TriKota: Boost headers not found in Boost_INCLUDE_DIRS
##    path <${Boost_INCLUDE_DIRS}>. Resetting path to Boost version distributed
##    within Dakota: ${PACKAGE_SOURCE_DIR}/Dakota/packages/boost")

  IF (EXISTS ${PACKAGE_SOURCE_DIR}/Dakota/packages/boost/boost/type.hpp)
    SET(Boost_INCLUDE_DIRS ${PACKAGE_SOURCE_DIR}/Dakota/packages/boost)
    SET(Boost_LIBRARY_DIRS ${PACKAGE_BINARY_DIR}/Dakota/install/lib)
    SET(TriKota_Dakota_Boost ON)
  ELSE()
    MESSAGE(FATAL_ERROR "\nTriKota Fatal Error: Dakota requires
     Boost, and the headers were not found in ${Boost_INCLUDE_DIRS}
     nor inside of Dakota at ${PACKAGE_SOURCE_DIR}/Dakota/packages/boost. ")
  ENDIF()
ELSE()
  IF (NOT EXISTS ${Boost_LIBRARY_DIRS}/libboost_signals.a)
      MESSAGE(FATAL_ERROR "\nTriKota Fatal Error: Dakota (Colin) now requires
       Boost Signals library: looking for   
       ${Boost_LIBRARY_DIRS}/libboost_signals.a ")
  ENDIF()
  SET(TriKota_Dakota_Boost OFF)
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

SET(TriKota_Rebuild_Dakota OFF)

EXECUTE_PROCESS(
  COMMAND echo "   TriKota: Dakota already configured, skipping reconfigure."
  WORKING_DIRECTORY ${PACKAGE_BINARY_DIR}/Dakota
)

ELSE()

SET(TriKota_Rebuild_Dakota ON)

SET(DAKOTA_CONFIGURE_COMMAND ${PACKAGE_SOURCE_DIR}/Dakota/configure
    CC=${CMAKE_C_COMPILER} 
    CXX=${CMAKE_CXX_COMPILER} 
    FC=${CMAKE_Fortran_COMPILER} 
    F77=${CMAKE_Fortran_COMPILER}
    --prefix=${PACKAGE_BINARY_DIR}/Dakota/install
    --with-teuchos-include=${PACKAGE_SOURCE_DIR}/../teuchos/src
    --with-teuchos-lib=${PACKAGE_BINARY_DIR}/../teuchos/src
    --disable-tests
    --with-plugin 
    --without-graphics
    --without-dot 
    --without-nlpql 
    --without-npsol)
IF(TPL_ENABLE_MPI)
  SET(DAKOTA_CONFIGURE_COMMAND ${DAKOTA_CONFIGURE_COMMAND}
    MPICXX=${CMAKE_CXX_COMPILER})
ELSE()
  SET(DAKOTA_CONFIGURE_COMMAND ${DAKOTA_CONFIGURE_COMMAND}
    --disable-mpi)
ENDIF()
IF(BUILD_SHARED_LIBS)
  SET(DAKOTA_CONFIGURE_COMMAND ${DAKOTA_CONFIGURE_COMMAND}
    --enable-shared
    --disable-static)
ELSE()
  SET(DAKOTA_CONFIGURE_COMMAND ${DAKOTA_CONFIGURE_COMMAND}
    --disable-shared
    --enable-static)
ENDIF()
IF(NOT TriKota_Dakota_Boost)
  SET(DAKOTA_CONFIGURE_COMMAND ${DAKOTA_CONFIGURE_COMMAND}
    --with-boost=${Boost_INCLUDE_DIRS}
    --with-boost-libdir=${Boost_LIBRARY_DIRS})
ENDIF()

EXECUTE_PROCESS(
  COMMAND echo "-->TriKota: Launching Dakota configure script -- takes several minutes."
  WORKING_DIRECTORY ${PACKAGE_BINARY_DIR}/Dakota
)

EXECUTE_PROCESS(
  WORKING_DIRECTORY ${PACKAGE_BINARY_DIR}/Dakota
  COMMAND ${DAKOTA_CONFIGURE_COMMAND}
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

IF(TriKota_Rebuild_Dakota)

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

ENDIF()

SET(Dakota_INCLUDE_DIRS 
  ${CMAKE_CURRENT_BINARY_DIR}/../Dakota
  ${CMAKE_CURRENT_BINARY_DIR}/../Dakota/install/include)

SET(DAKOTA_LIBS
    dakota pecos lhs surfpack conmin
    ddace fsudace jega opt
    psuade newmat ncsuopt sparsegrid 
    colin interfaces momh scolib pebbl tinyxml
    utilib 3po hopspack 
    cport amplsolver nidr dl 
    boost_signals
  )

LINK_DIRECTORIES(${PACKAGE_BINARY_DIR}/Dakota/install/lib
                 ${Boost_LIBRARY_DIRS})


ENDIF() # IF(NOT TriKota_ENABLE_DakotaCMake)

#
# Build the TriKota library
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(HEADERS "")
SET(SOURCES "")

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${Dakota_INCLUDE_DIRS}
  )

APPEND_SET(HEADERS
    TriKota_ConfigDefs.hpp
    TriKota_DirectApplicInterface.hpp
    TriKota_ThyraDirectApplicInterface.hpp
    TriKota_Driver.hpp
    TriKota_Version.hpp
  ) 

APPEND_SET(SOURCES
    TriKota_DirectApplicInterface.cpp
    TriKota_ThyraDirectApplicInterface.cpp
    TriKota_Driver.cpp
    TriKota_Version.cpp
  )


# Create a cache variable to store definitions needed when including Dakota
# headers
IF(TriKota_ENABLE_DakotaCMake)
  GET_PROPERTY(DAKOTA_DEFINITIONS_TMP DIRECTORY ../Dakota/src PROPERTY
  COMPILE_DEFINITIONS)
  # CMake apparently removes all of the -D's
  FOREACH(DAKOTA_DEF IN LISTS DAKOTA_DEFINITIONS_TMP)
    LIST(APPEND DAKOTA_DEFINITIONS "-D${DAKOTA_DEF}")
  ENDFOREACH()
ELSE()
  SET(DAKOTA_DEFINITIONS -DHAVE_CONFIG_H)
ENDIF()
SET(TRIKOTA_DAKOTA_DEFINITIONS ${DAKOTA_DEFINITIONS} CACHE INTERNAL 
"Variable used to communicate the ever-changing list of compiler definitions 
required when including Dakota headers")
ADD_DEFINITIONS(${TRIKOTA_DAKOTA_DEFINITIONS})

IF(TriKota_ENABLE_DakotaCMake)

# The DEPLIBS part is to tell CMake that trikota depends on the specified
# Dakota libraries.  The IMPORTEDLIBS part is to get the PACKAGE_ADD_LIBRARY
# macro to put the Dakota libraries TriKota_LIBRARIES when writing
# TriKotaConfig.cmake and Makefile.export.TriKota.  There should be a better
# way of doing this...
PACKAGE_ADD_LIBRARY(
  trikota
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS ${Dakota_EXTRA_LIBRARIES}
  IMPORTEDLIBS ${Dakota_LIBRARIES}
  )

ELSE()


  LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}
${PACKAGE_BINARY_DIR}/Dakota/install/lib/)

PACKAGE_ADD_LIBRARY(
  trikota
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS ${DAKOTA_LIBS}
  )

IF(TriKota_Rebuild_Dakota)

ADD_DEPENDENCIES(trikota TriKota_Dakota_install)

ENDIF()

# Custom "make install" command to copy contents of Dakota lib&include
INSTALL(DIRECTORY ${PACKAGE_BINARY_DIR}/Dakota/install/lib/ DESTINATION lib)
INSTALL(DIRECTORY ${PACKAGE_BINARY_DIR}/Dakota/install/include/ DESTINATION include)
# Temporary fix, to be removed by boostsignals TPL
INSTALL(FILES ${Boost_LIBRARY_DIRS}/libboost_signals.a  DESTINATION lib)

ENDIF() # IF(TriKota_ENABLE_DakotaCMake)
