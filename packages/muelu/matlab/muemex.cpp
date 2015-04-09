#include <stdio.h>
#include <string.h>

#include "MueLu_ConfigDefs.hpp"


#ifdef HAVE_ML_MATLAB


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
  printf("Hello world.\n");
}


#else
#error "Do not have MATLAB"
#endif

