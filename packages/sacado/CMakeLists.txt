INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)
INCLUDE(CombinedOption)

#
# A) Define the package
#

PACKAGE(Sacado)

#
# B) Set up package-specific options
#
ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_COMPLEX
  HAVE_SACADO_COMPLEX
  "Enable Sacado features for complex data type."
  ON
  )

OPTION(${PACKAGE_NAME}_ENABLE_PARAMLIB
  "Enable parameter library support"
  ON )

OPTION(${PACKAGE_NAME}_ENABLE_NESTEDFAD
  "Enable parameter library support"
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  SACADO_DEBUG
  "Enable debug code in sacado"
  ${Trilinos_ENABLE_DEBUG} )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_ADOLC
  HAVE_ADOLC
  "Enable ADOLC support for optional sacado tests"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_ADIC
  HAVE_ADIC
  "Enable ADIC support for optional sacado tests"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_RAD_NO_USING_STDCC
  RAD_NO_USING_STDCC
  "Enable if the compiler is confused by std::sin, etc., within namespace"
  OFF )

OPTION(${PACKAGE_NAME}_ENABLE_PERFORMANCE_TESTS
  "Override enable ${PACKAGE_NAME} performance tests."
  ${Trilinos_ENABLE_PERFORMANCE_TESTS}
  )

COMBINED_OPTION( Sacado_ENABLE_FEApp
  DEP_OPTIONS_NAMES Sacado_ENABLE_Epetra Sacado_ENABLE_EpetraExt
  DOCSTR "Enables compiling FEApp example finite element application"
    " that demonstrates incorporating Sacado into a applcaitons."
    "  Defaults to ON if Epetra, and EpetraExt are enabled," 
    " otherwise defaults to OFF.  Only takes affect when"
    " Sacado_ENABLE_EXAMPLES=ON."
  )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

# We are treating FEApp as a library, not an example
# This way packages that wish to use FEApp do not need to enable
# Sacado examples
ADD_SUBDIRECTORY(example/FEApp)

# This adds directories for non-performance tests as well
PACKAGE_ADD_PERFORMANCE_TEST_DIRECTORIES(
  test/performance)

# Must be after the performance tests so dependency on test/utils is picked
# up.
PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()


# ToDo: Add tradop tests use shell scripts in autotools to generate
# files.  Need to find a simpler way to do this for the CMake build
# system.
