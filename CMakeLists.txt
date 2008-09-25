# $Header$


#
# A) Basic top-level Trilinos stuff
#

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(Trilinos)
SET(Trilinos_VERSION "9.0d")

FIND_PACKAGE(Perl)

SET(TRILINOS_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(TRILINOS_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
SET(TRILINOS_TEST_CATEGORY "FRAMEWORK")

SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")


#
# B) Define the list of all Trilinos packages
#

# Setup paqckages to include ...
SET( Trilinos_PACKAGES
 Teuchos
 Epetra
 #Anasazi
 RBGen
)


#
# C) Set up user options and global variables that will be used throughout
#

INCLUDE(TrilinosGlobalHelpers)

OPTION(VERBOSE_CONFIGURE "Make the configure process verbose." OFF)

SET( Trilinos_ENABLE_ALL_PACKAGES "DEFAULT" CACHE STRING
  "Enable all packages (set to ON, OFF, or leave DEFAULT)." )
SET( Trilinos_ENABLE_TESTS "DEFAULT" CACHE STRING
  "Enable tests in all packages  (set to ON, OFF, or leave DEFAULT)." )
SET( ${PROJECT_NAME}_ENABLE_EXAMPLES "DEFAULT" CACHE STRING
  "Enable examples in all packages  (set to ON, OFF, or leave DEFAULT)." )

FOREACH(PACKAGE ${Trilinos_PACKAGES})
  TRILINOS_INSERT_STANDARD_PACAKGE_VARIABLES(${PACKAGE})
ENDFOREACH()

OPTION(BUILD_SHARED_LIBS "Build shared libraries." OFF)
OPTION(TRILINOS_ENABLE_MPI "Build packages using MPI for parallel codes." OFF)

SET(TRILINOS_INSTALL_INCLUDE_DIR "include/trilinos-${Trilinos_VERSION}")
SET(TRILINOS_INSTALL_LIB_INCLUDE_DIR "lib/trilinos-${Trilinos_VERSION}/include")

SET(Trilinos_INCLUDE_DIRS "" CACHE INTERNAL "Trilinos include directories")
SET(Trilinos_LIBRARY_DIRS "" CACHE INTERNAL "Trilinos library directories")
SET(Trilinos_LIBRARIES "" CACHE INTERNAL "Trilinos libraries")

MARK_AS_ADVANCED(BUILD_SHARED_LIBS)
MARK_AS_ADVANCED(BUILD_TESTING)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
MARK_AS_ADVANCED(CMAKE_BUILD_TYPE)
MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)
MARK_AS_ADVANCED(DART_TESTING_TIMEOUT)
MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH)
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
MARK_AS_ADVANCED(CMAKE_OSX_ARCHITECTURES)
MARK_AS_ADVANCED(CMAKE_OSX_SYSROOT)
#MARK_AS_ADVANCED(TRILINOS_BLAS_LIBRARY)
#MARK_AS_ADVANCED(TRILINOS_LAPACK_LIBRARY)

# 2008/09/25: rabartl: Above: I have not made the options for the BLAS and
# LAPACK libraries advanced since you almost always have to specify these and
# we don't want to hide that from users.

FOREACH(PACKAGE ${Trilinos_PACKAGES})
  TRILINOS_POSTPROCESS_STANDARD_PACAKGE_VARIABLES(${PACKAGE})
ENDFOREACH()


#
# D) Probe the environment on this computer
#

INCLUDE(TrilinosProbeEnv)

INCLUDE(TrilinosFortranMangling)


#
# E) Set up for testing with CTest and Trilinos test harness
#

INCLUDE(CTest)

ADD_CUSTOM_TARGET(
  runtests-serial
   ${PERL_EXECUTABLE} ${TRILINOS_HOME_DIR}/commonTools/test/utilities/runtests
  --trilinos-dir=${TRILINOS_HOME_DIR}
  --comm=serial
  --build-dir=${TRILINOS_BUILD_DIR}
  --category=${TRILINOS_TEST_CATEGORY}
  --output-dir=${TRILINOS_BUILD_DIR}/runtests-results
  )

ADD_CUSTOM_TARGET(
  runtests-mpi
   ${PERL_EXECUTABLE} ${TRILINOS_HOME_DIR}/commonTools/test/utilities/runtests
  --trilinos-dir=${TRILINOS_HOME_DIR}
  --comm=mpi
  --mpi-go="${TRILINOS_MPI_GO}"
  --max-proc=${MPIEXEC_MAX_NUMPROCS}
  --build-dir=${TRILINOS_BUILD_DIR}
  --category=${TRILINOS_TEST_CATEGORY}
  --output-dir=${TRILINOS_BUILD_DIR}/runtests-results
  )


#
# F) Setup for packaging and distribution
#

SET(CPACK_PACKAGE_DESCRIPTION "Trilinos provides algorithms and technologies for the solution of large-scale, complex multi-physics engineering and scientific problems.")
SET(CPACK_PACKAGE_FILE_NAME "trilinos-setup-${Trilinos_VERSION}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Trilinos ${Trilinos_VERSION}")
SET(CPACK_PACKAGE_REGISTRY_KEY "Trilinos ${Trilinos_VERSION}")
SET(CPACK_PACKAGE_NAME "trilinos")
SET(CPACK_PACKAGE_VENDOR "Sandia National Laboratories")
SET(CPACK_PACKAGE_VERSION "${Trilinos_VERSION}")
SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2")
SET(CPACK_SOURCE_FILE_NAME "trilinos-source-${Trilinos_VERSION}")

IF(WIN32)
  SET(CPACK_GENERATOR "NSIS")
ENDIF()

INCLUDE(CPack)


#
# G) Configure the individual pacakges
# 

FOREACH(PACKAGE ${Trilinos_PACKAGES})
  IF(${Trilinos_ENABLE_${PACKAGE}} STREQUAL "ON")
    STRING(TOLOWER ${PACKAGE} PACKAGE_DIR)
    # 2008/06/09: rabartl: Above, you can not assume that all package names are lower
    #   case.  This will not work for PyTrilinos, ForTrilinos, ThreadPool, etc.  You should
    #   just name the packages with their directory names!
    MESSAGE(STATUS "Trilinos_ENABLE_${PACKAGE}=${Trilinos_ENABLE_${PACKAGE}}")
    ADD_SUBDIRECTORY(packages/${PACKAGE_DIR})
    INCLUDE(${CMAKE_CURRENT_BINARY_DIR}/packages/${PACKAGE_DIR}/${PACKAGE}Config.cmake)
    LIST(APPEND Trilinos_INCLUDE_DIRS ${${PACKAGE}_INCLUDE_DIRS})
    LIST(APPEND Trilinos_LIBRARY_DIRS ${${PACKAGE}_LIBRARY_DIRS})
    LIST(APPEND Trilinos_LIBRARIES ${${PACKAGE}_LIBRARIES})
  ENDIF()
ENDFOREACH()

LIST(REMOVE_DUPLICATES Trilinos_INCLUDE_DIRS)
LIST(REMOVE_DUPLICATES Trilinos_LIBRARY_DIRS)
LIST(REMOVE_DUPLICATES Trilinos_LIBRARIES)


#
# H) Install commands
#

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TrilinosConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfig.cmake )


