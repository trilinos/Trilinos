#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>  

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
#define __FUNCT__ "ichol_blocked"
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

  Kokkos::Impl::Timer time;

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

  mxArray *Lmat = mxCreateSparse(m, n, nnz, mxREAL);

  CrsMatrixBase::size_type_array    lp((size_type*)(mxGetJc(Lmat)), m+1);
  CrsMatrixBase::ordinal_type_array lj((ordinal_type*)(mxGetIr(Lmat)), nnz);
  CrsMatrixBase::value_type_array   lx((value_type*)(mxGetPr(Lmat)), nnz);


  // output sparse matrix
  CrsMatrixBase L("CrsMatrixBase::Matlab::L",
                  m, n, nnz,
                  lp, lj, lx);

  L.copy(Uplo::Lower, A);

  Example::IChol<Uplo::Lower,Algo::LeftBlocked>::blocksize = 64;
  int r_val = Example::IChol<Uplo::Lower,Algo::LeftBlocked>::invoke(CrsMatrixView(L));
  if (r_val != 0)  {
    stringstream ss;
    ss << " Error in " << __FILE__ << ", " << __LINE__ << ", " __FUNCT__ << endl;
    ss << " IChol encounters negative diagonal: " << r_val << endl;
    mexErrMsgTxt (ss.str().c_str()) ;
  }

  // return the output matrix
  plhs[0] = Lmat;
  plhs[1] = mxCreateDoubleScalar(time.seconds()); 

  Kokkos::finalize();
}
