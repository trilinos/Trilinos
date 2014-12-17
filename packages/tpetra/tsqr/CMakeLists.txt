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

# Whether to build TbbTsqr and related classes.
#
# Enabled by default (unless disabled explicitly at the command line)
# if Trilinos is built with the TBB (Intel's Threading Building
# Blocks) TPL (third-party library) enabled.
TRIBITS_ADD_OPTION_AND_DEFINE(
  KokkosTSQR_ENABLE_TBB
  HAVE_KOKKOSTSQR_TBB
  "Enable Intel Threading Building Blocks (TBB) intranode parallelization of TSQR.  This option is enabled by default if you are building Trilinos with TBB enabled as a 'third-party library' (TPL), so you should not have to enable this option manually.  TSQR will work without this, but enabling it gives another parallelization option for TSQR."
  "${TPL_ENABLE_TBB}" 
  )

# KokkosTSQR_config.h gets created in the src/ subdirectory.
ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()
