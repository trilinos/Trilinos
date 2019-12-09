TRIBITS_SUBPACKAGE(TSQR)

#TRIBITS_ADD_DEBUG_OPTION()

# Whether to enable TSQR's complex arithmetic (std::complex<T>)
# support.
#
# Enabled by default (unless disabled explicitly at the command line)
# if Teuchos is built with complex arithmetic support.
TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosTSQR_ENABLE_Complex
  HAVE_KOKKOSTSQR_COMPLEX
  "Enable complex arithmetic (std::complex<T>) support for TSQR.  This is currently ON if Teuchos_ENABLE_COMPLEX is ON.  The default behavior may change as we migrate TSQR to depend on new Kokkos.  New Kokkos does not currently support complex arithmetic, but this will change."
  "${Teuchos_ENABLE_COMPLEX}" 
  )

# KokkosTSQR_config.h gets created in the src/ subdirectory.
ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()
