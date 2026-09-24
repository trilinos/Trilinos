INCLUDE(CombinedOption)

#
# A) Define the package
#

TRIBITS_PACKAGE(Sacado)

#
# B) Set up package-specific options
#
TRIBITS_SET_ST_FOR_DEV_MODE(ST_FOR_DEV_PT_FOR_RELEASE)

SET(Sacado_ENABLE_DEPRECATED_CODE_DEFAULT ON)
TRIBITS_ADD_OPTION_AND_DEFINE(
  Sacado_ENABLE_DEPRECATED_CODE
  SACADO_ENABLE_DEPRECATED_CODE
  "Whether Sacado enables deprecated code (that is, anything marked with the Sacado_DEPRECATED macro) at compile time.  Default is ON (deprecated code enabled)."
  ${Sacado_ENABLE_DEPRECATED_CODE_DEFAULT}
)

IF (Trilinos_ENABLE_COMPLEX_DOUBLE)
  SET(${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT  TRUE)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT  FALSE)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_COMPLEX
  HAVE_SACADO_COMPLEX
  "Enable Sacado features for complex data type."
  ${${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT}
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_SFAD_INIT_DEFAULT_CONSTRUCTOR
  SACADO_SFAD_INIT_DEFAULT_CONSTRUCTOR
  "Force SFad (in the new design) to initialize value and derivative components in the default constructor (adds additional runtime cost, but protects against uninitialized use)."
  OFF
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_HIERARCHICAL
  SACADO_GPU_HIERARCHICAL
  "Enable hierarchical parallelism, in which a team's vector dimension is used to partition the derivative array (CUDA, HIP and SYCL)."
  FALSE
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_HIERARCHICAL_DFAD
  SACADO_GPU_HIERARCHICAL_DFAD
  "Enable hierarchical parallelism using an approach specific to DFad (CUDA, HIP and SYCL).  A DFad whose derivative array belongs to a Kokkos::View is supported everywhere; building a DFad value inside a kernel additionally needs a device heap, which means CUDA with UVM and is unavailable on SYCL."
  FALSE
  )

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_MEMORY_POOL
  SACADO_KOKKOS_USE_MEMORY_POOL
  "Enable use of memory pool for DFad in Kokkos kernels (currently only useful for hierarchical parallelism)."
  FALSE
  )

IF(${PACKAGE_NAME}_ENABLE_HIERARCHICAL AND
    ${PACKAGE_NAME}_ENABLE_HIERARCHICAL_DFAD)
message(FATAL_ERROR "Cannot simultaneously enable ${PACKAGE_NAME}_ENABLE_HIERARCHICAL and ${PACKAGE_NAME}_ENABLE_HIERARCHICAL_DFAD")
ENDIF()

# ToDo: Above, split ${PACKAGE_NAME}_ENABLE_COMPLEX into
# ${PACKAGE_NAME}_ENABLE_COMPLEX_FLOAT and
# ${PACKAGE_NAME}_ENABLE_COMPLEX_DOUBLE.  See Trilinos GitHub #362.

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

# Narrowest vector width a SYCL TeamPolicy launch can deliver.
#
# Kokkos may hand a team a narrower vector width than it asked for, and the
# clamp that depends on the built kernel cannot be seen from the host.  Sacado
# therefore sizes partitioned Fad storage for the narrowest width the device
# could produce rather than the width that was requested.  Lowering this value
# costs registers; raising it above what the device can deliver is incorrect.
#
# This changes a type, so it must be the same for every translation unit in a
# build -- hence a configured value rather than a macro a source file can set.
#
# 8 is the smallest SIMD width Intel GPUs compile to and is safe anywhere.
# Ponte Vecchio reports sub-group sizes of 16 and 32, so 16 is safe there and
# halves the derivative storage per thread.
IF (Kokkos_ARCH_INTEL_PVC)
  SET(Sacado_SYCL_MIN_VECTOR_LENGTH_DEFAULT 16)
  SET(Sacado_SYCL_MIN_VECTOR_LENGTH_WHY "from Kokkos_ARCH_INTEL_PVC")
ELSE()
  SET(Sacado_SYCL_MIN_VECTOR_LENGTH_DEFAULT 8)
  SET(Sacado_SYCL_MIN_VECTOR_LENGTH_WHY "default")
ENDIF()

SET(Sacado_SYCL_MIN_VECTOR_LENGTH "${Sacado_SYCL_MIN_VECTOR_LENGTH_DEFAULT}"
  CACHE STRING
  "Narrowest vector width a SYCL team launch can deliver.  Sacado sizes partitioned Fad storage for this width, so lowering it costs registers and raising it above what the device can deliver is incorrect.")

IF (NOT "${Sacado_SYCL_MIN_VECTOR_LENGTH}" STREQUAL
    "${Sacado_SYCL_MIN_VECTOR_LENGTH_DEFAULT}")
  SET(Sacado_SYCL_MIN_VECTOR_LENGTH_WHY "user-specified")
ENDIF()

# The name #cmakedefine tests in Sacado_config.h.in is the macro's own, so
# mirror the cache variable onto it.
SET(SACADO_SYCL_MIN_VECTOR_LENGTH ${Sacado_SYCL_MIN_VECTOR_LENGTH})

IF (Kokkos_ENABLE_SYCL)
  MESSAGE(STATUS
    "Sacado: partitioned Fad sized for a minimum SYCL vector length of "
    "${Sacado_SYCL_MIN_VECTOR_LENGTH} (${Sacado_SYCL_MIN_VECTOR_LENGTH_WHY})")
ENDIF()

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
