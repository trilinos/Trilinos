

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


#
# A) Basic top-level Trilinos stuff
#

MESSAGE("")
MESSAGE("Configuring Trilinos build directory")
MESSAGE("")

PROJECT(Trilinos)

SET(Trilinos_VERSION "9.0d")
SET(TRILINOS_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(TRILINOS_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
SET(TRILINOS_TEST_CATEGORY "FRAMEWORK")

SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

INCLUDE(TrilinosGlobalHelpers)


#
# B) Define the list of all Trilinos packages
#

OPTION(Trilinos_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING "" OFF)
MARK_AS_ADVANCED(Trilinos_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING)

IF (Trilinos_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING)

  SET( Trilinos_PACKAGES
   Teuchos
   RTOp
   Epetra
   Triutils
   EpetraExt
   Thyra
   RBGen
  )

ELSE()

  SET( Trilinos_PACKAGES
   Teuchos
   RTOp
   Epetra
   #Anasazi
   RBGen
  )

ENDIF()

SET( Trilinos_PACKAGE_DIRS "")
FOREACH( PACKAGE ${Trilinos_PACKAGES})
  STRING(TOLOWER ${PACKAGE} PACKAGE_DIR)
  LIST(APPEND Trilinos_PACKAGE_DIRS ${PACKAGE_DIR})
ENDFOREACH()
#MESSAGE("Trilinos_PACKAGE_DIRS = ${Trilinos_PACKAGE_DIRS}")

LIST(LENGTH Trilinos_PACKAGES Trilinos_NUM_PACKAGES)
#MESSAGE("Trilinos_NUM_PACKAGES = ${Trilinos_NUM_PACKAGES}")


#
# C) Read in the package dependancies
#

MESSAGE("")
MESSAGE("Reading package dependencies ...")
MESSAGE("")

FOREACH(PACKAGE ${Trilinos_PACKAGES})

  IF (VERBOSE_CONFIGURE)
    MESSAGE("Reading package dependencies: ${PACKAGE}")
  ENDIF()

  SET(${PACKAGE}_FORWARD_LIB_REQUIRED_DEP_PACKAGES "")
  SET(${PACKAGE}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES "")
  SET(${PACKAGE}_FORWARD_TEST_REQUIRED_DEP_PACKAGES "")
  SET(${PACKAGE}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES "")

  STRING(TOLOWER ${PACKAGE} PACKAGE_DIR)

  INCLUDE(packages/${PACKAGE_DIR}/cmake/Dependencies.cmake)

  SET(${PACKAGE}_LIB_REQUIRED_DEP_PACKAGES ${LIB_REQUIRED_DEP_PACKAGES})
  SET(${PACKAGE}_LIB_OPTIONAL_DEP_PACKAGES ${LIB_OPTIONAL_DEP_PACKAGES})
  SET(${PACKAGE}_TEST_REQUIRED_DEP_PACKAGES ${TEST_REQUIRED_DEP_PACKAGES})
  SET(${PACKAGE}_TEST_OPTIONAL_DEP_PACKAGES ${TEST_OPTIONAL_DEP_PACKAGES})

  IF (VERBOSE_CONFIGURE)
    MESSAGE("${PACKAGE}_LIB_REQUIRED_DEP_PACKAGES = ${${PACKAGE}_LIB_REQUIRED_DEP_PACKAGES}")
    MESSAGE("${PACKAGE}_LIB_OPTIONAL_DEP_PACKAGES = ${${PACKAGE}_LIB_OPTIONAL_DEP_PACKAGES}")
    MESSAGE("${PACKAGE}_TEST_REQUIRED_DEP_PACKAGES = ${${PACKAGE}_TEST_REQUIRED_DEP_PACKAGES}")
    MESSAGE("${PACKAGE}_TEST_OPTIONAL_DEP_PACKAGES = ${${PACKAGE}_TEST_OPTIONAL_DEP_PACKAGES}")
    MESSAGE("")
  ENDIF()

  TRILINOS_APPEND_FORWARD_DEP_PACKAGES(${PACKAGE} LIB_REQUIRED_DEP_PACKAGES)
  TRILINOS_APPEND_FORWARD_DEP_PACKAGES(${PACKAGE} LIB_OPTIONAL_DEP_PACKAGES)
  TRILINOS_APPEND_FORWARD_DEP_PACKAGES(${PACKAGE} TEST_REQUIRED_DEP_PACKAGES)
  TRILINOS_APPEND_FORWARD_DEP_PACKAGES(${PACKAGE} TEST_OPTIONAL_DEP_PACKAGES)

ENDFOREACH()

IF (VERBOSE_CONFIGURE)

  MESSAGE("")
  MESSAGE("Printing forward package dependenies ...")
  MESSAGE("")

  FOREACH(PACKAGE ${Trilinos_PACKAGES})

    MESSAGE("Forward dependencies: ${PACKAGE}")

    MESSAGE("${PACKAGE}_FORWARD_LIB_REQUIRED_DEP_PACKAGES = ${${PACKAGE}_FORWARD_LIB_REQUIRED_DEP_PACKAGES}")
    MESSAGE("${PACKAGE}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES = ${${PACKAGE}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES}")
    MESSAGE("${PACKAGE}_FORWARD_TEST_REQUIRED_DEP_PACKAGES = ${${PACKAGE}_FORWARD_TEST_REQUIRED_DEP_PACKAGES}")
    MESSAGE("${PACKAGE}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES = ${${PACKAGE}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES}")
    MESSAGE("")

  ENDFOREACH()

ENDIF()


#
# D) Set up user options and global variables that will be used throughout
#

OPTION(VERBOSE_CONFIGURE "Make the configure process verbose." OFF)

SET( Trilinos_ENABLE_ALL_PACKAGES "" CACHE BOOL
  "Enable all packages (set to ON, OFF, or leave empty)." )
SET( Trilinos_ENABLE_TESTS "" CACHE BOOL
  "Enable tests in all packages  (set to ON, OFF, or leave empty)." )
SET( ${PROJECT_NAME}_ENABLE_EXAMPLES "" CACHE BOOL
  "Enable examples in all packages  (set to ON, OFF, or leave empty)." )

FOREACH(PACKAGE ${Trilinos_PACKAGES})
  TRILINOS_INSERT_STANDARD_PACKAGE_OPTIONS(${PACKAGE})
  TRILINOS_ADD_OPTIONAL_PACKAGE_ENABLES(${PACKAGE})
ENDFOREACH()

IF(WIN32 AND NOT CYGWIN)
  SET(Trilinos_ENABLE_FORTRAN_DEFAULT OFF)
ELSE()
  SET(Trilinos_ENABLE_FORTRAN_DEFAULT ON)
ENDIF()

OPTION(Trilinos_ENABLE_FORTRAN
  "Build Fortran library for this package."
  ${Trilinos_ENABLE_FORTRAN_DEFAULT}
  )

OPTION(BUILD_SHARED_LIBS "Build shared libraries." OFF)
OPTION(TRILINOS_ENABLE_MPI "Build packages using MPI for parallel codes." OFF)

SET(TRILINOS_INSTALL_INCLUDE_DIR "include/trilinos-${Trilinos_VERSION}")
SET(TRILINOS_INSTALL_LIB_INCLUDE_DIR "lib/trilinos-${Trilinos_VERSION}/include")

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


#
# E) Apply logic to enable Trilinos pacakges and tests
#

# Disable all forward packages and features based on package disables
# explicitly set by the user

MESSAGE("")
MESSAGE("Disabling forward pacakges that have a required dependancy on explicitly disabled packages ...")
MESSAGE("")
FOREACH(PACKAGE ${Trilinos_PACKAGES})
  IF (VERBOSE_CONFIGURE)
     MESSAGE(STATUS "Status: Trilinos_ENABLE_${PACKAGE} = ${Trilinos_ENABLE_${PACKAGE}}")
  ENDIF()
  IF (NOT ${Trilinos_ENABLE_${PACKAGE}})
    #MESSAGE(STATUS "Disabling packages depending on ${PACKAGE} ...")
    TRILINOS_DISABLE_FORWARD_REQUIRED_DEP_PACKAGES(${PACKAGE})
  ENDIF()
ENDFOREACH()

IF(Trilinos_ENABLE_ALL_PACKAGES)
  MESSAGE("")
  MESSAGE("Enabling all packages that are not currently disabled ...")
  MESSAGE("")
  FOREACH(PACKAGE ${Trilinos_PACKAGES})
    TRILINOS_APPLY_ALL_PACKAGE_ENABLES(${PACKAGE})
  ENDFOREACH()
ENDIF()

# Enable all optional forward packages if asked
IF(Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES)
  # ToDo: Implement!
  SET(Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES ON)
ENDIF()

# Enable all optional packages in the current set of enabled Trilinos packages
IF(Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES)
  # ToDo: Implement!
ENDIF()

IF (Trilinos_ENABLE_TESTS OR Trilinos_EANBLE_EXAMPLES)
  MESSAGE("")
  MESSAGE("Enabling all tests and examples that can be enabled ...")
  MESSAGE("")
  FOREACH(PACKAGE ${Trilinos_PACKAGES})
    TRILINOS_APPLY_TEST_EXAMPLE_EANBLES(${PACKAGE})
  ENDFOREACH()
ENDIF()
# NOTE: Above, we enable tests and examples here, before the remaining required
# packages so that we don't enable tests that don't need to be enabled based
# on the use of the option Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES.

# Enable all required dependent packages for the current set of enabled
# packages.

# ToDo: Implement!

MESSAGE("")
MESSAGE("Enabling all optional package dependencies that can be if both sets of packages are enabled ...")
MESSAGE("")
FOREACH(PACKAGE ${Trilinos_PACKAGES})
  TRILINOS_POSTPROCESS_OPTIONAL_PACKAGE_ENABLES(${PACKAGE})
ENDFOREACH()

IF (Trilinos_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING)
  MESSAGE("")
  MESSAGE(SEND_ERROR "Shortcircuiting after dependency tracking ...")
  RETURN()
ENDIF()


#
# F) Probe the environment on this computer
#

MESSAGE("")
MESSAGE("Probing the environment ...")
MESSAGE("")

INCLUDE(TrilinosProbeEnv)

IF(Trilinos_ENABLE_FORTRAN)
  ENABLE_LANGUAGE(Fortran)
ENDIF()

INCLUDE(TrilinosFortranMangling)


#
# G) Set up for testing with CTest and Trilinos test harness
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
# H) Setup for packaging and distribution
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
# I) Configure individual pacakges
# 

MESSAGE("")
MESSAGE("Configuring individual Trilinos packages ...")
MESSAGE("")

SET(Trilinos_INCLUDE_DIRS "" CACHE INTERNAL "")
SET(Trilinos_LIBRARY_DIRS "" CACHE INTERNAL "")
SET(Trilinos_LIBRARIES "" CACHE INTERNAL "")

FOREACH(PACKAGE ${Trilinos_PACKAGES})
  IF(Trilinos_ENABLE_${PACKAGE})
    STRING(TOLOWER ${PACKAGE} PACKAGE_DIR)
    # 2008/06/09: rabartl: Above, you can not assume that all package names are lower
    #   case.  This will not work for PyTrilinos, ForTrilinos, ThreadPool, etc.  You should
    #   just name the packages with their directory names!
    SET(PACKAGE_NAME_GLOBAL ${PACKAGE}) # For checking
    ADD_SUBDIRECTORY(packages/${PACKAGE_DIR})
    LIST(APPEND Trilinos_INCLUDE_DIRS ${${PACKAGE}_INCLUDE_DIRS})
    LIST(APPEND Trilinos_LIBRARY_DIRS ${${PACKAGE}_LIBRARY_DIRS})
    LIST(APPEND Trilinos_LIBRARIES ${${PACKAGE}_LIBRARIES})
  ENDIF()
ENDFOREACH()

LIST(REMOVE_DUPLICATES Trilinos_INCLUDE_DIRS)
LIST(REMOVE_DUPLICATES Trilinos_LIBRARY_DIRS)
LIST(REMOVE_DUPLICATES Trilinos_LIBRARIES)


#
# J) Install-related commands
#

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TrilinosConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfig.cmake )
