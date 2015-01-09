#include <Kokkos_Core.hpp>
//#include "util.hpp"

//#include "crs_matrix_base.hpp"
//#include "graph_helper_scotch.hpp"

// matlab headers
#include <mex.h>
#include <matrix.h>

using namespace std;

typedef double  value_type;
typedef mwIndex ordinal_type;
typedef mwIndex size_type;

//typedef Kokkos::Serial space_type;

//typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBase;
//typedef Example::CrsMatrixView<CrsMatrixBase> CrsMatrixView;

//typedef Example::Uplo Uplo;
//typedef Example::Algo Algo;

#undef  __FUNCT__ 
#define __FUNCT__ "scotch"
void mexFunction(int nlhs,
                 mxArray *plhs [],
                 int nrhs,
                 const mxArray *prhs []) {
  const mxArray *matrix = prhs[0];

  if (nrhs < 1) {
    mexErrMsgTxt (" Incorrect number of arguments\n") ;
    return;
  }

  const size_type *ap = mxGetJc(matrix);
  const ordinal_type *aj = mxGetIr(matrix);
  const value_type *ax = mxGetPr(matrix);

  const ordinal_type nrow = mxGetM(matrix);
  const ordinal_type ncol = mxGetN(matrix);

  size_type nnz  = mxGetNumberOfElements(matrix);

  Kokkos::initialize();
  string name = typeid(Kokkos::DefaultExecutionSpace).name();
  mexPrintf("\n Kokkos is initialized with a default spaceL: %s\n", name.c_str());

  Kokkos::finalize();
}
