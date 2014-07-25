INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(CombinedOption)

#
# A) Define the package
#

TRIBITS_PACKAGE(Sacado)

#
# B) Set up package-specific options
#
TRIBITS_SET_ST_FOR_DEV_MODE(ST_FOR_DEV_PT_FOR_RELEASE)

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_COMPLEX
  HAVE_SACADO_COMPLEX
  "Enable Sacado features for complex data type."
  ${ST_FOR_DEV_PT_FOR_RELEASE}
  )

OPTION(${PACKAGE_NAME}_ENABLE_PARAMLIB
  "Enable parameter library support"
  ON )

OPTION(${PACKAGE_NAME}_ENABLE_NESTEDFAD
  "Enable parameter library support"
  ON )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  SACADO_DEBUG
  "Enable debug code in sacado"
  ${${PROJECT_NAME}_ENABLE_DEBUG} )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_UNINIT
  HAVE_SACADO_UNINIT
  "Enable use of uninit_f2c for debugging Rad"
  OFF )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_VIEW_SPEC
  HAVE_SACADO_VIEW_SPEC
  "Enable Sacado specializations for Kokkos::View<...>"
  ON )

SET(Sacado_ENABLE_ADOLC_default OFF)
IF(TPL_ENABLE_ADOLC)
  SET(Sacado_ENABLE_ADOLC_default ON)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(Sacado_ENABLE_ADOLC
  HAVE_ADOLC
  "Enable ADOLC support for optional sacado tests"
  ${Sacado_ENABLE_ADOLC_default} )

SET(Sacado_ENABLE_ADIC_default OFF)
IF(TPL_ENABLE_ADIC)
  SET(Sacado_ENABLE_ADIC_default ON)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(Sacado_ENABLE_ADIC
  HAVE_ADIC
  "Enable ADIC support for optional sacado tests"
  ${Sacado_ENABLE_ADIC_default} )

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_RAD_NO_USING_STDCC
  RAD_NO_USING_STDCC
  "Enable if the compiler is confused by std::sin, etc., within namespace"
  OFF )

# Checks for C99 TR1 math functions such as acosh, asinh, atanh ...
INCLUDE(CheckCXXSourceCompiles)
SET(SOURCE_C99_TR1_CMATH
  "
#include <cmath>
int main()
{
  double x = 1.0;
  double y = std::acosh(x);
  return 0;
} 
  "
  )
CHECK_CXX_SOURCE_COMPILES("${SOURCE_C99_TR1_CMATH}" 
  HAS_C99_TR1_CMATH)

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

# Must be after the performance tests so dependency on test/utils is picked
# up.
TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_EXCLUDE_FILES(
  example/FEApp/experimental
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()


# ToDo: Add tradop tests use shell scripts in autotools to generate
# files.  Need to find a simpler way to do this for the CMake build
# system.
