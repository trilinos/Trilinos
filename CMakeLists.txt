

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


#
# A) Basic top-level Trilinos stuff
#

MESSAGE("")
MESSAGE("Configuring Trilinos build directory")
MESSAGE("")

PROJECT(Trilinos C)

SET(Trilinos_VERSION "9.0d")
SET(TRILINOS_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(TRILINOS_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
SET(TRILINOS_TEST_CATEGORY "FRAMEWORK")

SET(CMAKE_MODULE_PATH "${TRILINOS_HOME_DIR}/cmake")

SET(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

INCLUDE(TrilinosGlobalHelpers)
INCLUDE(Advanced_Set)
INCLUDE(Advanced_Option)

OPTION(Trilinos_VERBOSE_CONFIGURE "Make the Trilinos configure process verbose." OFF)


#
# B) Define the list of all Trilinos packages
#

# Define the Trilinos package names and directories
SET( Trilinos_PACKAGES_AND_DIRS
 Teuchos               teuchos
 RTOp                  rtop
 Epetra                epetra
 Zoltan                zoltan
 Tpetra                tpetra
 Triutils              triutils
 EpetraExt             epetraext
 Isorropia             isorropia
 Thyra                 thyra
 AztecOO               aztecoo
 Galeri                galeri
 Ifpack                ifpack
 Komplex               komplex
 #Anasazi               anasazi
 RBGen                 rbgen
)

# Hack to allow me to set the exact list of packages so that I can have stable
# unit tests as packages are added and dependencies change
ADVANCED_OPTION(Trilinos_PACKAGES_UNITTEST_OVERRIDE "" OFF)
IF (Trilinos_PACKAGES_UNITTEST_OVERRIDE)
  SET( Trilinos_PACKAGES_AND_DIRS
    Teuchos          teuchos
    Epetra           epetra
    Zoltan           zoltan
    RTOp             rtop
    Triutils         triutils
    EpetraExt        epetraext
    Isorropia        isorropia
    Thyra            thyra
    AztecOO          aztecoo
    Galeri           galeri
    Amesos           amesos
    Ifpack           ifpack
    ML               ml
    Belos            belos
    RBGen            rbgen
    Stratimikos      stratimikos
    )
ENDIF()

# Temporary hack to allow me to state more package dependencies that we have
# full CMake support for.  This will go away once all of Trilinos is
# CMakeified.
ADVANCED_OPTION(Trilinos_PACKAGES_HACK_OVERRIDE "" OFF)
IF (Trilinos_PACKAGES_HACK_OVERRIDE)
  SET( Trilinos_PACKAGES_AND_DIRS
    Teuchos               teuchos
    RTOp                  rtop
    Kokkos                kokkos
    Epetra                epetra
    Stokhos               stokhos
    Sacado                sacado
    Zoltan                zoltan
    Shards                shards
    Intrepid              intrepid
    Tpetra                tpetra
    Triutils              triutils
    EpetraExt             epetraext
    Thyra                 thyra
    Isorropia             isorropia
    Pliris                pliris
    Claps                 claps
    AztecOO               aztecoo
    Galeri                galeri
    Amesos                amesos
    Ifpack                ifpack
    Komplex               komplex
    ML                    ml
    Belos                 belos
    Stratimikos           stratimikos
    Meros                 meros
    FEI                   fei
    RBGen                 rbgen
    Anasazi               anasazi
    ThreadPool            ThreadPool
    Phalanx               phalanx
    Pamgen                pamgen
    Phdmesh               phdmesh
    NOX                   nox
    Moertel               moertel
    TrilinosCouplings     trilinoscouplings
    Rythmos               rythmos
    MOOCHO                moocho
    Aristos               aristos
    Sundance              Sundance
    TriKota               TriKota
    CTrilinos             CTrilinos
    ForTrilinos           ForTrilinos
    PyTrilinos            PyTrilinos
    WebTrilinos           WebTrilinos
    Didasko               didasko
    NewPackage            new_package
    )
ENDIF()

# Process the packages lists and dependencies

TRILINOS_PROCESS_PACKAGES_AND_DIRS_LISTS()

FOREACH(PACKAGE_IDX RANGE ${Trilinos_LAST_PACKAGE_IDX})
  LIST(GET Trilinos_PACKAGES ${PACKAGE_IDX} PACKAGE)
  LIST(GET Trilinos_PACKAGE_DIRS ${PACKAGE_IDX} PACKAGE_DIR)
  TRILINOS_READ_PACKAGE_DEPENDENCIES(${PACKAGE} ${PACKAGE_DIR})
ENDFOREACH()

ADVANCED_OPTION(Trilinos_DUMP_PACKAGE_DEPENDENCIES
  "Dump the package dependency information." OFF)

IF (Trilinos_VERBOSE_CONFIGURE OR Trilinos_DUMP_PACKAGE_DEPENDENCIES)
  MESSAGE("")
  MESSAGE("Printing package dependencies ...")
  MESSAGE("")
  PRINT_VAR(Trilinos_PACKAGES)
  MESSAGE("")
  FOREACH(PACKAGE ${Trilinos_PACKAGES})
    TRILINOS_PRINT_PACKAGE_DEPENDENCIES(${PACKAGE})
    MESSAGE("")
  ENDFOREACH()
ENDIF()


#
# C) Define the list of offical Trilinos TPLs
#

SET(Trilinos_TPLs
   BLAS         ON
   LAPACK       ON
   MPI          OFF
   Boost        OFF
   ParMETIS     OFF
   Scotch       OFF
   PaToH        OFF
   )

# ToDo: Set up user cache varaibles


#
# D) Set up user options and global variables that will be used throughout
#

MESSAGE("")
MESSAGE("Setting up major user options ...")
MESSAGE("")

SET( Trilinos_ENABLE_ALL_PACKAGES OFF CACHE BOOL
  "Enable all packages (set to ON, OFF, or leave empty)." )

SET(Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES OFF CACHE BOOL
  "Recursively enable all packages that have required or optional dependencies for set of enabled packages." )

SET(Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES OFF CACHE BOOL
  "Recursively enable all optional packages for set of enabled packages." )

SET( Trilinos_ENABLE_TESTS "" CACHE STRING
  "Enable tests in all packages  (set to ON, OFF, or leave empty)." )

SET(Trilinos_ENABLE_EXAMPLES "" CACHE STRING
  "Enable examples in all packages  (set to ON, OFF, or leave empty).  If left empty, then this will be set to ON if Trilinos_ENABLE_TESTS=ON" )

SET(Trilinos_ENABLE_PERFORMANCE_TESTS "" CACHE BOOL
  "Enable performance tests in all packages that support them" )

IF (Trilinos_ENABLE_TESTS AND Trilinos_ENABLE_EXAMPLES STREQUAL "")
  MESSAGE(STATUS "Setting Trilinos_ENABLE_EXAMPLES=ON because Trilinos_ENABLE_TESTS=ON")
  SET(Trilinos_ENABLE_EXAMPLES ON)
ENDIF()

FOREACH(PACKAGE ${Trilinos_PACKAGES})
  TRILINOS_INSERT_STANDARD_PACKAGE_OPTIONS(${PACKAGE})
  TRILINOS_ADD_OPTIONAL_PACKAGE_ENABLES(${PACKAGE})
ENDFOREACH()

ADVANCED_OPTION(Trilinos_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING
  "Shortcircut after dependency handling is complete"
  OFF)

IF (Trilinos_PACKAGES_UNITTEST_OVERRIDE OR Trilinos_PACKAGES_HACK_OVERRIDE)
  SET(Trilinos_DEPS_XML_OUTPUT_FILE_DEFAULT
    "${TRILINOS_BUILD_DIR}/cmake/TrilinosPackageDependencies.xml" )
ELSE()
  SET(Trilinos_DEPS_XML_OUTPUT_FILE_DEFAULT
    "${TRILINOS_HOME_DIR}/cmake/python/data/TrilinosPackageDependencies.xml" )
ENDIF()
ADVANCED_SET(Trilinos_DEPS_XML_OUTPUT_FILE
  ${Trilinos_DEPS_XML_OUTPUT_FILE_DEFAULT}
  CACHE STRING
  "XML file that Trilinos dependenices used by tools are written to (if not empty)." )

IF(WIN32 AND NOT CYGWIN)
  SET(Trilinos_DEPS_HTML_OUTPUT_FILE_DEFAULT "")
ELSE()
  SET(Trilinos_DEPS_HTML_OUTPUT_FILE_DEFAULT
    "${TRILINOS_BUILD_DIR}/TrilinosPackageDependenciesTable.html" )
ENDIF()
ADVANCED_SET(Trilinos_DEPS_HTML_OUTPUT_FILE
  ${Trilinos_DEPS_HTML_OUTPUT_FILE_DEFAULT}
  CACHE STRING
  "HTML Trilinos dependenices file that will be written to (if not empty)." )

ADVANCED_OPTION(Trilinos_ENABLE_CXX
  "Enable the C++ compiler and related code"
  ON )

IF(WIN32 AND NOT CYGWIN)
  SET(Trilinos_ENABLE_FORTRAN_DEFAULT OFF)
ELSE()
  SET(Trilinos_ENABLE_FORTRAN_DEFAULT ON)
ENDIF()

OPTION(Trilinos_ENABLE_FORTRAN
  "Enable the Fortran compiler and related code"
  ${Trilinos_ENABLE_FORTRAN_DEFAULT} )

ADVANCED_OPTION(BUILD_SHARED_LIBS "Build shared libraries." OFF)

OPTION(Trilinos_ENABLE_MPI "Build packages using MPI for parallel codes." OFF)

SET(TRILINOS_INSTALL_INCLUDE_DIR "include/trilinos-${Trilinos_VERSION}")
SET(TRILINOS_INSTALL_LIB_INCLUDE_DIR "lib/trilinos-${Trilinos_VERSION}/include")

MARK_AS_ADVANCED(BUILD_TESTING)
MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
MARK_AS_ADVANCED(CMAKE_BUILD_TYPE)
MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)
MARK_AS_ADVANCED(DART_TESTING_TIMEOUT)
MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH)
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
MARK_AS_ADVANCED(CMAKE_OSX_ARCHITECTURES)
MARK_AS_ADVANCED(CMAKE_OSX_SYSROOT)
MARK_AS_ADVANCED(Trilinos_BLAS_LIBRARY)
MARK_AS_ADVANCED(Trilinos_LAPACK_LIBRARY)


#
# E) Apply logic to enable Trilinos packages and tests
#

TRILINOS_PRINT_ENABLED_PACKAGE_LIST("\nExplicitly enabled packages on input")

TRILINOS_ADJUST_PACKAGE_ENABLES()

TRILINOS_PRINT_ENABLED_PACKAGE_LIST("\nFinal set of enabled packages")

IF (Trilinos_DEPS_XML_OUTPUT_FILE)
  MESSAGE("" )
  MESSAGE("Dumping the XML dependencies file ${Trilinos_DEPS_XML_OUTPUT_FILE} ..." )
  MESSAGE("")
  TRILINOS_DUMP_DEPS_XML_FILE()
ENDIF()

IF (Trilinos_DEPS_HTML_OUTPUT_FILE)
  MESSAGE("" )
  MESSAGE("Dumping the HTML dependencies webpage file ${Trilinos_DEPS_HTML_OUTPUT_FILE} ..." )
  MESSAGE("" )
  EXECUTE_PROCESS(
    COMMAND ${TRILINOS_HOME_DIR}/cmake/python/dump-package-dep-table.py
      --input-xml-deps-file=${Trilinos_DEPS_XML_OUTPUT_FILE}
      --output-html-deps-file=${Trilinos_DEPS_HTML_OUTPUT_FILE} )
ENDIF()

IF (Trilinos_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING)
  MESSAGE("")
  MESSAGE(SEND_ERROR "Shortcircuiting after dependency tracking ...")
  RETURN()
ENDIF()



#
# F) Set up the environment on this computer
#

MESSAGE("")
MESSAGE("Probing the environment ...")
MESSAGE("")

INCLUDE(TrilinosProbeEnv)


#
# G) Set up for testing with CTest and Trilinos test harness
#

INCLUDE(CTest)

IF (WIN32 AND NOT CYGWIN)
  SET(Trilinos_ENABLE_NATIVE_TEST_HARNESS_DEFAULT OFF)
ELSE()
  SET(Trilinos_ENABLE_NATIVE_TEST_HARNESS_DEFAULT ON)
ENDIF()

ADVANCED_OPTION(Trilinos_ENABLE_NATIVE_TEST_HARNESS
  "Enable the native Trilinos perl-based test harness."
  ${Trilinos_ENABLE_NATIVE_TEST_HARNESS_DEFAULT} )

IF (Trilinos_ENABLE_NATIVE_TEST_HARNESS)

  ADD_CUSTOM_TARGET(
    runtests-serial
     ${PERL_EXECUTABLE} ${TRILINOS_HOME_DIR}/commonTools/test/utilities/runtests
    --trilinos-dir=${TRILINOS_HOME_DIR}
    --comm=serial
    --build-dir=${TRILINOS_BUILD_DIR}
    --category=${TRILINOS_TEST_CATEGORY}
    --output-dir=${TRILINOS_BUILD_DIR}/runtests-results
    )

  IF (Trilinos_ENABLE_MPI)
  
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

  ENDIF()

ENDIF()

IF (WIN32)
  SET(Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS_DEFAULT OFF)
ELSE()
  SET(Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS_DEFAULT ON)
ENDIF()

# 2008/10/17: rabartl: Above, I can not turn these tests on by default
# with cygwin because the custom script target is not working for some
# reason.

ADVANCED_OPTION(Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS
  "Enable dependency unit tests."
  ${Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS_DEFAULT}
  )

IF (Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS)
  ADD_SUBDIRECTORY(cmake/DependencyUnitTests)
  ADD_SUBDIRECTORY(cmake/python/UnitTests)
ENDIF()


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
  SET(CPACK_NSIS_MODIFY_PATH ON)
ENDIF()

INCLUDE(CPack)


#
# I) Configure individual packages
# 

MESSAGE("")
MESSAGE("Configuring individual Trilinos packages ...")
MESSAGE("")

SET(Trilinos_INCLUDE_DIRS "" CACHE INTERNAL "")
SET(Trilinos_LIBRARY_DIRS "" CACHE INTERNAL "")
SET(Trilinos_LIBRARIES "" CACHE INTERNAL "")

FOREACH(PACKAGE_IDX RANGE ${Trilinos_LAST_PACKAGE_IDX})
  LIST(GET Trilinos_PACKAGES ${PACKAGE_IDX} PACKAGE)
  LIST(GET Trilinos_PACKAGE_DIRS ${PACKAGE_IDX} PACKAGE_DIR)
  IF(Trilinos_ENABLE_${PACKAGE})
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
