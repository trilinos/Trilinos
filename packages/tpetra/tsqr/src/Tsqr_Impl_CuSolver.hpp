#ifndef TSQR_IMPL_CUSOLVER_HPP
#define TSQR_IMPL_CUSOLVER_HPP

#include "TpetraTSQR_config.h"
#if defined(HAVE_TPETRATSQR_CUBLAS) && defined(HAVE_TPETRATSQR_CUSOLVER)
#include "Tsqr_Impl_CuBlasHandle.hpp"
#include "Tsqr_Impl_CuSolverHandle.hpp"
#if defined(HAVE_TPETRATSQR_COMPLEX)
#  include <complex>
#endif // HAVE_TPETRATSQR_COMPLEX

namespace TSQR {
namespace Impl {

template<class Scalar>
class CuSolver {
public:
  CuSolver (CuSolverHandle handle, int* const info);

  int
  geqrfBufferSize (const int nrows,
                   const int ncols,
                   Scalar A_raw[],
                   const int lda);

  void
  geqrf (const int nrows,
         const int ncols,
         Scalar A[],
         const int lda,
         Scalar tau[],
         Scalar work[],
         const int lwork);

  int
  unmqrBufferSize (const char side,
                   const char trans,
                   const int nrows,
                   const int ncols_C,
                   const int ncols_Q,
                   const Scalar Q[],
                   const int ldq,
                   const Scalar tau[],
                   const Scalar C[],
                   const int ldc);

  void
  unmqr (const char side,
         const char trans,
         const int nrows,
         const int ncols_C,
         const int ncols_Q,
         const Scalar Q[],
         const int ldq,
         const Scalar tau[],
         Scalar C[],
         const int ldc,
         Scalar work[],
         const int lwork);

private:
  CuSolverHandle handle_;
  int* info_; // DEVICE MEMORY
};

extern template class CuSolver<double>;
extern template class CuSolver<float>;
#if defined(HAVE_TPETRATSQR_COMPLEX)
extern template class CuSolver<std::complex<double>>;
extern template class CuSolver<std::complex<float>>;
#endif // defined(HAVE_TPETRATSQR_COMPLEX)

} // namespace Impl
} // namespace TSQR

#endif // HAVE_TPETRATSQR_CUBLAS && HAVE_TPETRATSQR_CUSOLVER
#endif // TSQR_IMPL_CUSOLVER_HPP
