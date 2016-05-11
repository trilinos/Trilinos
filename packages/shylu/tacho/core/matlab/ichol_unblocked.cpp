#include <Kokkos_Core.hpp>
#include "util.hpp"

#include "crs_matrix_base.hpp"
#include "crs_matrix_view.hpp"
#include "crs_row_view.hpp"

#include "ichol.hpp"

// matlab headers
#include <mex.h>
#include <matrix.h>

using namespace std;

typedef double  value_type;
typedef mwSignedIndex ordinal_type;
typedef mwIndex size_type;

typedef Kokkos::Serial space_type;

typedef Example::CrsMatrixBase<value_type,ordinal_type,size_type,space_type> CrsMatrixBase;
typedef Example::CrsMatrixView<CrsMatrixBase> CrsMatrixView;

typedef Example::Uplo Uplo;
typedef Example::Algo Algo;

#undef  __FUNCT__ 
#define __FUNCT__ "ichol_unblocked"
void mexFunction(int nlhs,
                 mxArray *plhs [],
                 int nrhs,
                 const mxArray *prhs []) {
  if (nrhs < 1) {
    mexErrMsgTxt (" Incorrect number of arguments\n") ;
    return;
  }

  Kokkos::initialize();
  string name = typeid(Kokkos::DefaultExecutionSpace).name();
  mexPrintf("\n Kokkos is initialized with a default spaceL: %s\n", name.c_str());

  // input sparse matrix 
  // - assume that the matrix stores both symmetric parts
  // - matrix is square and non-singular spd
  // - matlab stores sparse matrix in CSC not CSR
  // - since ichol only deal with symmetric matrix, it does not matter
  const mxArray *Amat = prhs[0];

  const ordinal_type m = mxGetM(Amat);
  const ordinal_type n = mxGetN(Amat);

  const size_type nnz = mxGetNumberOfElements(Amat);

  CrsMatrixBase::size_type_array    ap((size_type*)(mxGetJc(Amat)), m+1);
  CrsMatrixBase::ordinal_type_array aj((ordinal_type*)(mxGetIr(Amat)), nnz);
  CrsMatrixBase::value_type_array   ax((value_type*)(mxGetPr(Amat)), nnz);

  CrsMatrixBase A("CrsMatrixBase::Matlab::A",
                  m, n, nnz,
                  ap, aj, ax);

  mxArray *Umat = mxCreateSparse(m, n, nnz, mxREAL);

  CrsMatrixBase::size_type_array    up((size_type*)(mxGetJc(Umat)), m+1);
  CrsMatrixBase::ordinal_type_array uj((ordinal_type*)(mxGetIr(Umat)), nnz);
  CrsMatrixBase::value_type_array   ux((value_type*)(mxGetPr(Umat)), nnz);


  // output sparse matrix
  CrsMatrixBase U("CrsMatrixBase::Matlab::U",
                  m, n, nnz,
                  up, uj, ux);

  U.copy(Uplo::Upper, A);

  int r_val = Example::IChol<Uplo::Upper,Algo::RightUnblockedOpt1>::invoke(CrsMatrixView(U));
  if (r_val != 0)  {
    stringstream ss;
    ss << " Error in " << __FILE__ << ", " << __LINE__ << ", " __FUNCT__ << endl;
    ss << " IChol encounters negative diagonal: " << r_val << endl;
    mexErrMsgTxt (ss.str().c_str()) ;
  }

  // return the output matrix
  plhs[0] = Umat;

  Kokkos::finalize();
}
