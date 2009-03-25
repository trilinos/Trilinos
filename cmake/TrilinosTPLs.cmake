#
# Define the list of TPLs
#
# Here, the TPLs must be listed in the order of increasing
# dependencies (if such dependencies exist).  
#
# NOTE: The TPLs that are not required for "Primary Stable" code are
# explicitly turned 'OFF' by default.  That will result in all packages
# that have required dependencies on these to be disabled by default.
# We want the default configure to enable only "Primary Stable" code.
# You have to explicitly enable code that is not "Primary Stable"
# code and therefore, you have to explicitly enabled these extra
# TPLs.
#

SET(Trilinos_TPLS_ENABLED
  MPI            ""
  BLAS           ""
  LAPACK         ""
  Boost          OFF
  Scotch         OFF
  METIS          OFF
  ParMETIS       OFF
  PaToH          OFF
  CppUnit        OFF
  ADOLC          OFF
  TVMET          OFF
  MF             OFF
  ExodusII       OFF
  ZoltanTpl      OFF
  y12m           OFF
  SuperLUDist    OFF
  )

# NOTES:
#
# (*) ParMETIS must be listed after Scotch because the
#     ParMETIS include directories must come before the
#     Scotch include directories.
