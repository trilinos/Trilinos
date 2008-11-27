

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


#
# A) Basic top-level Trilinos stuff
#

MESSAGE("")
MESSAGE("Configuring Trilinos build directory")
MESSAGE("")

PROJECT(Trilinos C) # Across all of Trilinos C++ and Fortran are actually optional!

SET(Trilinos_VERSION "9.0d")
SET(TRILINOS_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(TRILINOS_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
SET(TRILINOS_TEST_CATEGORY "FRAMEWORK")

SET(CMAKE_MODULE_PATH ${TRILINOS_HOME_DIR}/cmake ${TRILINOS_HOME_DIR}/cmake/utils)

SET(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)

INCLUDE(TrilinosGlobalHelpers)
INCLUDE(Advanced_Set)
INCLUDE(Advanced_Option)

OPTION(Trilinos_VERBOSE_CONFIGURE "Make the Trilinos configure process verbose." OFF)


#
# B) Define the list of all Trilinos packages and TPLs
#

# Define the Trilinos package names and directories
SET( Trilinos_PACKAGES_AND_DIRS
 Teuchos               teuchos
 RTOp                  rtop
 Epetra                epetra
 Zoltan                zoltan
 Triutils              triutils
 Tpetra                tpetra
 EpetraExt             epetraext
 Isorropia             isorropia
 Thyra                 thyra
 AztecOO               aztecoo
 Galeri                galeri
 Amesos                amesos
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
    Triutils              triutils
    Tpetra                tpetra
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

TRILINOS_PROCESS_PACKAGES_AND_DIRS_LISTS()

SET(Trilinos_TPLS_ENABLED
  BLAS           ""
  LAPACK         ""
  MPI            ""
  Boost          OFF
  #METIS          OFF
  #ParMETIS       OFF
  #Scotch         OFF
  #PaToH          OFF
  )
  # Above, the TPLs that are not required for "Primary Stable" code areexplicitly
  # turned off by default.  That will result in all packages that have required
  # dependencies on these to be disabled by default.  We want the default configure
  # to enable only "Primary Stable" code.  You have explicitly enable code that is
  # not "Primary Stable" code.

TRILINOS_PROCESS_TPLS_LISTS()


#
# C) Process the package and TPL dependencies
#

TRILINOS_READ_ALL_PACKAGE_DEPENDENCIES()


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

IF (Trilinos_ENABLE_TESTS AND Trilinos_ENABLE_EXAMPLES STREQUAL "")
  MESSAGE(STATUS "Setting Trilinos_ENABLE_EXAMPLES=ON because Trilinos_ENABLE_TESTS=ON")
  SET(Trilinos_ENABLE_EXAMPLES ON)
ENDIF()

SET(Trilinos_ENABLE_PERFORMANCE_TESTS "" CACHE BOOL
  "Enable performance tests in all packages that support them" )

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


#
# E) Apply logic to enable Trilinos packages and tests
#

TRILINOS_PRINT_ENABLED_PACKAGE_LIST(
  "\nExplicitly enabled packages on input" ON)
TRILINOS_PRINT_ENABLED_PACKAGE_LIST(
  "\nExplicitly disabled packages on input" OFF)

TRILINOS_ADJUST_PACKAGE_ENABLES()

TRILINOS_PRINT_ENABLED_PACKAGE_LIST(
  "\nFinal set of enabled packages" ON)

TRILINOS_PRINT_ENABLED_TPL_LIST(
  "\nFinal set of enabled TPLs" ON)

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
# G) Go get the information for all enabled TPLS
#

MESSAGE("")
MESSAGE("Getting information for all enabled TPLs ...")
MESSAGE("")

TRILINOS_PROCESS_ENABLED_TPLS()


#
# H) Set up for testing with CTest and Trilinos test harness
#

MESSAGE("")
MESSAGE("Setting up testing support ...")
MESSAGE("")

TRILINOS_SETUP_TESTING_SUPPORT()


#
# I) Setup for packaging and distribution
#

TRILINOS_DEFINE_PACKAGING()


#
# J) Configure individual packages
# 

MESSAGE("")
MESSAGE("Configuring individual enabled Trilinos packages ...")
MESSAGE("")

TRILINOS_CONFIGURE_ENABLED_PACKAGES()


#
# K) Install-related commands
#

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TrilinosConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfig.cmake )
