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
PACKAGE_ARCH_SET_SS_FOR_DEV_MODE(SS_FOR_DEV_PS_FOR_RELEASE)

ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_COMPLEX
  HAVE_SACADO_COMPLEX
  "Enable Sacado features for complex data type."
  ${SS_FOR_DEV_PS_FOR_RELEASE}
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

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_UNINIT
  HAVE_SACADO_UNINIT
  "Enable use of uninit_f2c for debugging Rad"
  OFF )

SET(Sacado_ENABLE_ADOLC_default OFF)
IF(TPL_ENABLE_ADOLC)
  SET(Sacado_ENABLE_ADOLC_default ON)
ENDIF()
ADD_OPTION_AND_DEFINE(Sacado_ENABLE_ADOLC
  HAVE_ADOLC
  "Enable ADOLC support for optional sacado tests"
  ${Sacado_ENABLE_ADOLC_default} )

SET(Sacado_ENABLE_ADIC_default OFF)
IF(TPL_ENABLE_ADIC)
  SET(Sacado_ENABLE_ADIC_default ON)
ENDIF()
ADD_OPTION_AND_DEFINE(Sacado_ENABLE_ADIC
  HAVE_ADIC
  "Enable ADIC support for optional sacado tests"
  ${Sacado_ENABLE_ADIC_default} )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_RAD_NO_USING_STDCC
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
PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES()

PACKAGE_ARCH_EXCLUDE_FILES(
  example/FEApp/experimental
  )

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()


# ToDo: Add tradop tests use shell scripts in autotools to generate
# files.  Need to find a simpler way to do this for the CMake build
# system.
