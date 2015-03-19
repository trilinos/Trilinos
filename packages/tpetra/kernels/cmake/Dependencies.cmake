# FIXME (mfh 18 Dec 2014) MKL and CUSPARSE are optional TPLs, in the
# sense that one need not include the MKL resp. CUSPARSE header files
# if those TPLs aren't enabled.  However, we should still list them as
# optional TPLs.

TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_REQUIRED_PACKAGES KokkosCore KokkosContainers
  LIB_OPTIONAL_TPLS quadmath
  TEST_REQUIRED_PACKAGES Gtest
)
