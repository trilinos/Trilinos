# TPL for linking to matlab libraries (no mex or engine)
#
# Unlike the MATLAB TPL, this TPL seems to work with MPI enabled.


TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( MATLABLib
  REQUIRED_HEADERS mat.h matrix.h
  REQUIRED_LIBS_NAMES mx mat
  )
