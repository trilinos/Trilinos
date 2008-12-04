#
# Define the list of TPLs
#
# Here, the TPLs must be listed in the order of increasing
# dependencies (if such dependencies exist).  
#
# NOTE: The TPLs that are not required for "Primary Stable" code are
# explicitly turned off by default.  That will result in all packages
# that have required dependencies on these to be disabled by default.
# We want the default configure to enable only "Primary Stable" code.
# You have to explicitly enable code that is not "Primary Stable"
# code.
#

SET(Trilinos_TPLS_ENABLED
  BLAS           ""
  LAPACK         ""
  MPI            ""
  Boost          OFF
  METIS          OFF
  ParMETIS       OFF
  Scotch         OFF
  PaToH          OFF
  CppUnit        OFF
  ADOLC          OFF
  )
