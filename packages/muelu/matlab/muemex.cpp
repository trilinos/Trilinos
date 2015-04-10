#include <stdio.h>
#include <string.h>


#include "MueLu_config.hpp"

#ifdef HAVE_MUELU_MATLAB
#include "muemex.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){
  printf("Hello world.\n");
}


#else
#error "Do not have MATLAB"
#endif

