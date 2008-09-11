# $Header$

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(Trilinos)
SET(${PROJECT_NAME}_VERSION "9.0d")

FIND_PACKAGE(Perl)

SET(TRILINOS_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(TRILINOS_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
SET(TRILINOS_TEST_CATEGORY "FRAMEWORK")

#2008/07/10: ejguill: Updated to use the same command 
#that the perl test harness uses to get the name of the host
EXECUTE_PROCESS(
  COMMAND uname -n
  OUTPUT_VARIABLE TRILINOS_HOSTNAME
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

SET(${PROJECT_NAME}_INCLUDE_DIRS "" CACHE INTERNAL "${PROJECT_NAME} include directories")
SET(${PROJECT_NAME}_LIBRARY_DIRS "" CACHE INTERNAL "${PROJECT_NAME} library directories")
SET(${PROJECT_NAME}_LIBRARIES "" CACHE INTERNAL "${PROJECT_NAME} libraries")

SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Setup paqckages to include ...
SET( ${PROJECT_NAME}_PACKAGES
 Teuchos
 Epetra
 #Anasazi
 RBGen
)

# Setup user options ...
FOREACH(PACKAGE ${${PROJECT_NAME}_PACKAGES})
  OPTION(${PROJECT_NAME}_ENABLE_${PACKAGE} "Build the ${PACKAGE} package." ON)
ENDFOREACH(PACKAGE ${PROJECT_NAME}_PACKAGES)
OPTION(BUILD_SHARED_LIBS "Build shared libraries." OFF)
OPTION(TRILINOS_ENABLE_MPI "Build packages using MPI for parallel codes." OFF)
#OPTION(TRILINOS_INSTALL_BLAS "Include BLAS library in installations / packages." OFF)
#OPTION(TRILINOS_INSTALL_LAPACK "Include LAPACK library in installations / packages." OFF)

# Setup global requirements ...
FIND_LIBRARY(TRILINOS_BLAS_LIBRARY NAMES blas blas_win32 DOC "Path to the BLAS implementation")
FIND_LIBRARY(TRILINOS_LAPACK_LIBRARY NAMES lapack lapack_win32 DOC "Path to the LAPACK implementation")

IF(CMAKE_SIZEOF_VOID_P GREATER 4)
  ADD_DEFINITIONS(-DEPETRA_ADDRESS64BIT)
ENDIF()

# Set global options that will be inserted into config.h files ...
INCLUDE(TrilinosProbeEnv)

# MPI related configuration
IF(TRILINOS_ENABLE_MPI)
  SET(HAVE_MPI TRUE)
  # 2008/06/09: rabartl: Above, HAVE_MPI gets set in ConfigureMPI.cmake!
  # You should not have to set it up again.
  INCLUDE(ConfigureMPI)
ENDIF(TRILINOS_ENABLE_MPI)

# Setup Fortran name mangling
INCLUDE(TrilinosFortranMangling)

# Setup platform-specific compiler options ...
IF(WIN32)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_DEPRECATE)
ENDIF()

# Setup testing ...
INCLUDE(CTest)

# Handle some common installation stuff ...

# This is the path where packages should install their header files.
SET(TRILINOS_INSTALL_INCLUDE_DIR "include/trilinos-${${PROJECT_NAME}_VERSION}")
# This is the path where packages should install their header configuration files ...
SET(TRILINOS_INSTALL_LIB_INCLUDE_DIR "lib/trilinos-${${PROJECT_NAME}_VERSION}/include")

# Disable these two because they're problematic.  The issue is that the procedure
# to install these libraries correctly changes radically based on platform and
# whether the libraries are static or shared.
# 
# While cmake knows how to install the target it builds, it doesn't know anything
# about these external libraries - they're just files to cmake!

#IF(TRILINOS_INSTALL_BLAS)
#  INSTALL(FILES ${TRILINOS_BLAS_LIBRARY} DESTINATION lib)
#ENDIF(TRILINOS_INSTALL_BLAS)

#IF(TRILINOS_INSTALL_LAPACK)
#  INSTALL(FILES ${TRILINOS_LAPACK_LIBRARY} DESTINATION lib)
#ENDIF(TRILINOS_INSTALL_LAPACK)


FOREACH(PACKAGE ${${PROJECT_NAME}_PACKAGES})
  SET(${PACKAGE}_INCLUDE_DIRS "" CACHE INTERNAL "${PACKAGE} include directories")
  SET(${PACKAGE}_LIBRARY_DIRS "" CACHE INTERNAL "${PACKAGE} library directories")
  SET(${PACKAGE}_LIBRARIES "" CACHE INTERNAL "${PACKAGE} libraries")
ENDFOREACH()

# Setup individual packages ...
ADD_SUBDIRECTORY(packages)

# Setup packaging and distribution ...
SET(CPACK_PACKAGE_DESCRIPTION "Trilinos provides algorithms and technologies for the solution of large-scale, complex multi-physics engineering and scientific problems.")
SET(CPACK_PACKAGE_FILE_NAME "trilinos-setup-${${PROJECT_NAME}_VERSION}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Trilinos ${${PROJECT_NAME}_VERSION}")
SET(CPACK_PACKAGE_REGISTRY_KEY "Trilinos ${${PROJECT_NAME}_VERSION}")
SET(CPACK_PACKAGE_NAME "trilinos")
SET(CPACK_PACKAGE_VENDOR "Sandia National Laboratories")
SET(CPACK_PACKAGE_VERSION "${${PROJECT_NAME}_VERSION}")
SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2")
SET(CPACK_SOURCE_FILE_NAME "trilinos-source-${${PROJECT_NAME}_VERSION}")

IF(WIN32)
  SET(CPACK_GENERATOR "NSIS")
ENDIF()

INCLUDE(CPack)

# Hide advanced options ...
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
MARK_AS_ADVANCED(TRILINOS_BLAS_LIBRARY)
MARK_AS_ADVANCED(TRILINOS_LAPACK_LIBRARY)

# 2008/06/09: rabartl: Above: I would not make the options for the BLAS and
# LAPACK libraries advanced since you almost always have to specify these and
# we don't want to hide that from users.

FOREACH(PACKAGE ${${PROJECT_NAME}_PACKAGES})
  IF(${${PROJECT_NAME}_ENABLE_${PACKAGE}} STREQUAL "ON")
    STRING(TOLOWER ${PACKAGE} PACKAGE_DIR)
    # 2008/06/09: rabartl: Above, you can not assume that all package names are lower
    #   case.  This will not work for PyTrilinos, ForTrilinos, ThreadPool, etc.  You should
    #   just name the packages with their directory names!
    MESSAGE(STATUS "${PROJECT_NAME}_ENABLE_${PACKAGE}=${${PROJECT_NAME}_ENABLE_${PACKAGE}}")
    ADD_SUBDIRECTORY(packages/${PACKAGE_DIR})
    INCLUDE(${CMAKE_CURRENT_BINARY_DIR}/packages/${PACKAGE_DIR}/${PACKAGE}Config.cmake)
    LIST(APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${${PACKAGE}_INCLUDE_DIRS})
    LIST(APPEND ${PROJECT_NAME}_LIBRARY_DIRS ${${PACKAGE}_LIBRARY_DIRS})
    LIST(APPEND ${PROJECT_NAME}_LIBRARIES ${${PACKAGE}_LIBRARIES})
  ENDIF()
ENDFOREACH()

LIST(APPEND ${PROJECT_NAME}_LIBRARIES ${TRILINOS_BLAS_LIBRARY} ${TRILINOS_LAPACK_LIBRARY})

LIST(REMOVE_DUPLICATES ${PROJECT_NAME}_INCLUDE_DIRS)
LIST(REMOVE_DUPLICATES ${PROJECT_NAME}_LIBRARY_DIRS)
LIST(REMOVE_DUPLICATES ${PROJECT_NAME}_LIBRARIES)

##############################################################################
# Configure <PACKAGE>Config.cmake for the build tree.
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
               ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake)
##############################################################################

# 2008/06/09: rabartl: General Comment: A lot of the logic that is duplicated
# in the autotools system is aggregated here at the base level.  What this
# means is that packages are not really independent entities but instead must
# be delt with at the top Trilinos level.  This is a shift in the logical
# relationship between packages and the Trilinos framework but I suspect that
# this is a more accurate description of what is really going on.  However, it
# is very dangerous to have this logic at the top level and then letting
# individual package developers modify it or not letting them change it.  A
# top-heavy build model is what SIERRA has and it does not provide enough
# freedom at the package level (i.e. product level for SIERRA).

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

ADD_CUSTOM_TARGET(
  tests
  DEPENDS
)
