// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tsqr_Impl_CuSolver.hpp"
#if defined(HAVE_TPETRATSQR_CUBLAS) && defined(HAVE_TPETRATSQR_CUSOLVER)
#include "Tsqr_Impl_CuSolverHandle.hpp"
#include "Tsqr_Impl_CuTypes.hpp"
#include "Teuchos_Assert.hpp"

namespace TSQR {
namespace Impl {

template<class T>
class RawCuSolver {};

template<>
class RawCuSolver<double> {
public:
  using impl_scalar_type = double;

  static cusolverStatus_t
  compute_QR_lwork (cusolverDnHandle_t handle,
                    int m,
                    int n,
                    impl_scalar_type* A,
                    int lda,
                    int *lwork)
  {
    return cusolverDnDgeqrf_bufferSize (handle, m, n, A, lda, lwork);
  }

  static cusolverStatus_t
  compute_QR (cusolverDnHandle_t handle,
              int m,
              int n,
              impl_scalar_type* A,
              int lda,
              impl_scalar_type* tau,
              impl_scalar_type* work,
              int lwork,
              int* info)
  {
    return cusolverDnDgeqrf (handle, m, n, A, lda, tau,
                             work, lwork, info);
  }

  static cusolverStatus_t
  apply_Q_factor_lwork (cusolverDnHandle_t handle,
                        cublasSideMode_t side,
                        cublasOperation_t trans,
                        int m,
                        int n,
                        int k,
                        const impl_scalar_type* A,
                        int lda,
                        const impl_scalar_type* tau,
                        const impl_scalar_type* C,
                        int ldc,
                        int *lwork)
  {
    return cusolverDnDormqr_bufferSize (handle, side, trans,
                                        m, n, k, A, lda, tau,
                                        C, ldc, lwork);
  }

  static cusolverStatus_t
  apply_Q_factor (cusolverDnHandle_t handle,
                  cublasSideMode_t side,
                  cublasOperation_t trans,
                  int m,
                  int n,
                  int k,
                  const impl_scalar_type* A,
                  int lda,
                  const impl_scalar_type* tau,
                  impl_scalar_type* C,
                  int ldc,
                  impl_scalar_type* work,
                  int lwork,
                  int* devInfo)
  {
    return cusolverDnDormqr (handle, side, trans, m, n, k,
                             A, lda, tau, C, ldc,
                             work, lwork, devInfo);
  }

  static cusolverStatus_t
  compute_explicit_Q_lwork (cusolverDnHandle_t handle,
                            int m,
                            int n,
                            int k,
                            const impl_scalar_type *A,
                            int lda,
                            const impl_scalar_type *tau,
                            int *lwork)
  {
    return cusolverDnDorgqr_bufferSize(handle, m, n, k, A, lda,
                                       tau, lwork);
  }

  static cusolverStatus_t
  compute_explicit_Q (cusolverDnHandle_t handle,
                      int m,
                      int n,
                      int k,
                      impl_scalar_type *A,
                      int lda,
                      const impl_scalar_type *tau,
                      impl_scalar_type *work,
                      int lwork,
                      int *devInfo)
  {
    return cusolverDnDorgqr(handle, m, n, k, A, lda, tau,
                            work, lwork, devInfo);
  }
};

template<>
class RawCuSolver<float> {
public:
  using impl_scalar_type = float;

  static cusolverStatus_t
  compute_QR_lwork (cusolverDnHandle_t handle,
                    int m,
                    int n,
                    impl_scalar_type* A,
                    int lda,
                    int *lwork)
  {
    return cusolverDnSgeqrf_bufferSize (handle, m, n, A, lda, lwork);
  }

  static cusolverStatus_t
  compute_QR (cusolverDnHandle_t handle,
              int m,
              int n,
              impl_scalar_type* A,
              int lda,
              impl_scalar_type* tau,
              impl_scalar_type* work,
              int lwork,
              int* info)
  {
    return cusolverDnSgeqrf (handle, m, n, A, lda, tau,
                             work, lwork, info);
  }

  static cusolverStatus_t
  apply_Q_factor_lwork (cusolverDnHandle_t handle,
                        cublasSideMode_t side,
                        cublasOperation_t trans,
                        int m,
                        int n,
                        int k,
                        const impl_scalar_type* A,
                        int lda,
                        const impl_scalar_type* tau,
                        const impl_scalar_type* C,
                        int ldc,
                        int *lwork)
  {
    return cusolverDnSormqr_bufferSize (handle, side, trans,
                                        m, n, k, A, lda, tau,
                                        C, ldc, lwork);
  }

  static cusolverStatus_t
  apply_Q_factor (cusolverDnHandle_t handle,
                  cublasSideMode_t side,
                  cublasOperation_t trans,
                  int m,
                  int n,
                  int k,
                  const impl_scalar_type* A,
                  int lda,
                  const impl_scalar_type* tau,
                  impl_scalar_type* C,
                  int ldc,
                  impl_scalar_type* work,
                  int lwork,
                  int* devInfo)
  {
    return cusolverDnSormqr (handle, side, trans, m, n, k,
                             A, lda, tau, C, ldc,
                             work, lwork, devInfo);
  }

  static cusolverStatus_t
  compute_explicit_Q_lwork (cusolverDnHandle_t handle,
                            int m,
                            int n,
                            int k,
                            const impl_scalar_type *A,
                            int lda,
                            const impl_scalar_type *tau,
                            int *lwork)
  {
    return cusolverDnSorgqr_bufferSize(handle, m, n, k, A, lda,
                                       tau, lwork);
  }

  static cusolverStatus_t
  compute_explicit_Q (cusolverDnHandle_t handle,
                      int m,
                      int n,
                      int k,
                      impl_scalar_type *A,
                      int lda,
                      const impl_scalar_type *tau,
                      impl_scalar_type *work,
                      int lwork,
                      int *devInfo)
  {
    return cusolverDnSorgqr(handle, m, n, k, A, lda, tau,
                            work, lwork, devInfo);
  }
};

#if defined(HAVE_TPETRATSQR_COMPLEX)
template<>
class RawCuSolver<CudaValue<std::complex<double>>::type> {
public:
  using impl_scalar_type = CudaValue<std::complex<double>>::type;

  static cusolverStatus_t
  compute_QR_lwork (cusolverDnHandle_t handle,
                    int m,
                    int n,
                    impl_scalar_type* A,
                    int lda,
                    int *lwork)
  {
    return cusolverDnZgeqrf_bufferSize (handle, m, n, A, lda, lwork);
  }

  static cusolverStatus_t
  compute_QR (cusolverDnHandle_t handle,
              int m,
              int n,
              impl_scalar_type* A,
              int lda,
              impl_scalar_type* tau,
              impl_scalar_type* work,
              int lwork,
              int* info)
  {
    return cusolverDnZgeqrf (handle, m, n, A, lda, tau,
                             work, lwork, info);
  }

  static cusolverStatus_t
  apply_Q_factor_lwork (cusolverDnHandle_t handle,
                        cublasSideMode_t side,
                        cublasOperation_t trans,
                        int m,
                        int n,
                        int k,
                        const impl_scalar_type* A,
                        int lda,
                        const impl_scalar_type* tau,
                        const impl_scalar_type* C,
                        int ldc,
                        int *lwork)
  {
    return cusolverDnZunmqr_bufferSize (handle, side, trans,
                                        m, n, k, A, lda, tau,
                                        C, ldc, lwork);
  }

  static cusolverStatus_t
  apply_Q_factor (cusolverDnHandle_t handle,
                  cublasSideMode_t side,
                  cublasOperation_t trans,
                  int m,
                  int n,
                  int k,
                  const impl_scalar_type* A,
                  int lda,
                  const impl_scalar_type* tau,
                  impl_scalar_type* C,
                  int ldc,
                  impl_scalar_type* work,
                  int lwork,
                  int* devInfo)
  {
    return cusolverDnZunmqr (handle, side, trans, m, n, k,
                             A, lda, tau, C, ldc,
                             work, lwork, devInfo);
  }

  static cusolverStatus_t
  compute_explicit_Q_lwork (cusolverDnHandle_t handle,
                            int m,
                            int n,
                            int k,
                            const impl_scalar_type *A,
                            int lda,
                            const impl_scalar_type *tau,
                            int *lwork)
  {
    return cusolverDnZungqr_bufferSize(handle, m, n, k, A, lda,
                                       tau, lwork);
  }

  static cusolverStatus_t
  compute_explicit_Q (cusolverDnHandle_t handle,
                      int m,
                      int n,
                      int k,
                      impl_scalar_type *A,
                      int lda,
                      const impl_scalar_type *tau,
                      impl_scalar_type *work,
                      int lwork,
                      int *devInfo)
  {
    return cusolverDnZungqr(handle, m, n, k, A, lda, tau,
                            work, lwork, devInfo);
  }
};

template<>
class RawCuSolver<CudaValue<std::complex<float>>::type> {
public:
  using impl_scalar_type = CudaValue<std::complex<float>>::type;

  static cusolverStatus_t
  compute_QR_lwork (cusolverDnHandle_t handle,
                    int m,
                    int n,
                    impl_scalar_type* A,
                    int lda,
                    int *lwork)
  {
    return cusolverDnCgeqrf_bufferSize (handle, m, n, A, lda, lwork);
  }

  static cusolverStatus_t
  compute_QR (cusolverDnHandle_t handle,
              int m,
              int n,
              impl_scalar_type* A,
              int lda,
              impl_scalar_type* tau,
              impl_scalar_type* work,
              int lwork,
              int* info)
  {
    return cusolverDnCgeqrf (handle, m, n, A, lda, tau,
                             work, lwork, info);
  }

  static cusolverStatus_t
  apply_Q_factor_lwork (cusolverDnHandle_t handle,
                        cublasSideMode_t side,
                        cublasOperation_t trans,
                        int m,
                        int n,
                        int k,
                        const impl_scalar_type* A,
                        int lda,
                        const impl_scalar_type* tau,
                        const impl_scalar_type* C,
                        int ldc,
                        int *lwork)
  {
    return cusolverDnCunmqr_bufferSize (handle, side, trans,
                                        m, n, k, A, lda, tau,
                                        C, ldc, lwork);
  }

  static cusolverStatus_t
  apply_Q_factor (cusolverDnHandle_t handle,
                  cublasSideMode_t side,
                  cublasOperation_t trans,
                  int m,
                  int n,
                  int k,
                  const impl_scalar_type* A,
                  int lda,
                  const impl_scalar_type* tau,
                  impl_scalar_type* C,
                  int ldc,
                  impl_scalar_type* work,
                  int lwork,
                  int* devInfo)
  {
    return cusolverDnCunmqr (handle, side, trans, m, n, k,
                             A, lda, tau, C, ldc,
                             work, lwork, devInfo);
  }

  static cusolverStatus_t
  compute_explicit_Q_lwork (cusolverDnHandle_t handle,
                            int m,
                            int n,
                            int k,
                            const impl_scalar_type *A,
                            int lda,
                            const impl_scalar_type *tau,
                            int *lwork)
  {
    return cusolverDnCungqr_bufferSize(handle, m, n, k, A, lda,
                                       tau, lwork);
  }

  static cusolverStatus_t
  compute_explicit_Q (cusolverDnHandle_t handle,
                      int m,
                      int n,
                      int k,
                      impl_scalar_type *A,
                      int lda,
                      const impl_scalar_type *tau,
                      impl_scalar_type *work,
                      int lwork,
                      int *devInfo)
  {
    return cusolverDnCungqr(handle, m, n, k, A, lda, tau,
                            work, lwork, devInfo);
  }
};
#endif // defined(HAVE_TPETRATSQR_COMPLEX)

template<class Scalar>
CuSolver<Scalar>::CuSolver (int* const info) :
  handle_ (getCuSolverHandleSingleton ()), info_ (info)
{}

template<class Scalar>
CuSolver<Scalar>::
CuSolver (const std::shared_ptr<CuSolverHandle>& handle,
          int* const info) :
  handle_ (handle), info_ (info)
{}

template<class Scalar>
int
CuSolver<Scalar>::
compute_QR_lwork (const int nrows,
                  const int ncols,
                  Scalar A[],
                  const int lda) const
{
  cusolverDnHandle_t rawHandle = handle_->getHandle ();
  int lwork = 0;

  using IST = typename CudaValue<Scalar>::type;
  IST* A_raw = reinterpret_cast<IST*> (A);

  using impl_type = RawCuSolver<IST>;
  const auto status =
    impl_type::compute_QR_lwork (rawHandle, nrows, ncols,
                                 A_raw, lda, &lwork);
  TEUCHOS_ASSERT( status == CUSOLVER_STATUS_SUCCESS );
  return lwork;
}

template<class Scalar>
void
CuSolver<Scalar>::
compute_QR (const int nrows,
            const int ncols,
            Scalar A[],
            const int lda,
            Scalar tau[],
            Scalar work[],
            const int lwork) const
{
  cusolverDnHandle_t rawHandle = handle_->getHandle ();

  using IST = typename CudaValue<Scalar>::type;
  IST* A_raw = reinterpret_cast<IST*> (A);
  IST* tau_raw = reinterpret_cast<IST*> (tau);
  IST* work_raw = reinterpret_cast<IST*> (work);

  using impl_type = RawCuSolver<IST>;
  const auto status =
    impl_type::compute_QR (rawHandle, nrows, ncols, A_raw, lda,
                           tau_raw, work_raw, lwork, info_);
  TEUCHOS_ASSERT( status == CUSOLVER_STATUS_SUCCESS );
}

template<class Scalar>
int
CuSolver<Scalar>::
apply_Q_factor_lwork (const char side,
                      const char trans,
                      const int nrows,
                      const int ncols_C,
                      const int ncols_Q,
                      const Scalar Q[],
                      const int ldq,
                      const Scalar tau[],
                      Scalar C[],
                      const int ldc) const
{
  cusolverDnHandle_t rawHandle = handle_->getHandle ();
  const cublasSideMode_t cuSide = cuBlasSide (side);
  const cublasOperation_t cuTrans = cuBlasTrans (trans);
  int lwork = 0;

  using IST = typename CudaValue<Scalar>::type;
  const IST* Q_raw = reinterpret_cast<const IST*> (Q);
  const IST* tau_raw = reinterpret_cast<const IST*> (tau);
  const IST* C_raw = reinterpret_cast<const IST*> (C);

  using impl_type = RawCuSolver<IST>;
  const auto status =
    impl_type::apply_Q_factor_lwork (rawHandle, cuSide, cuTrans,
                                     nrows, ncols_C, ncols_Q,
                                     Q_raw, ldq, tau_raw,
                                     C_raw, ldc, &lwork);
  TEUCHOS_ASSERT( status == CUSOLVER_STATUS_SUCCESS );
  return lwork;
}

template<class Scalar>
void
CuSolver<Scalar>::
apply_Q_factor (const char side,
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
                const int lwork) const
{
  cusolverDnHandle_t rawHandle = handle_->getHandle ();
  const cublasSideMode_t cuSide = cuBlasSide (side);
  const cublasOperation_t cuTrans = cuBlasTrans (trans);

  using IST = typename CudaValue<Scalar>::type;
  const IST* Q_raw = reinterpret_cast<const IST*> (Q);
  const IST* tau_raw = reinterpret_cast<const IST*> (tau);
  IST* C_raw = reinterpret_cast<IST*> (C);
  IST* work_raw = reinterpret_cast<IST*> (work);

  using impl_type = RawCuSolver<IST>;
  const auto status =
    impl_type::apply_Q_factor (rawHandle, cuSide, cuTrans,
                               nrows, ncols_C, ncols_Q,
                               Q_raw, ldq, tau_raw, C_raw, ldc,
                               work_raw, lwork, info_);
  TEUCHOS_ASSERT( status == CUSOLVER_STATUS_SUCCESS );
}

template<class Scalar>
int
CuSolver<Scalar>::
compute_explicit_Q_lwork(const int m, const int n, const int k,
                         Scalar A[], const int lda,
                         const Scalar tau[]) const
{
  cusolverDnHandle_t rawHandle = handle_->getHandle ();
  int lwork = 0;

  using IST = typename CudaValue<Scalar>::type;
  const IST* A_raw = reinterpret_cast<const IST*> (A);
  const IST* tau_raw = reinterpret_cast<const IST*> (tau);

  using impl_type = RawCuSolver<IST>;
  const auto status =
    impl_type::compute_explicit_Q_lwork (rawHandle, m, n, k,
                                         A_raw, lda, tau_raw, &lwork);
  TEUCHOS_ASSERT( status == CUSOLVER_STATUS_SUCCESS );
  return lwork;
}

template<class Scalar>
void
CuSolver<Scalar>::
compute_explicit_Q(const int m, const int n, const int k,
                   Scalar A[], const int lda,
                   const Scalar tau[],
                   Scalar work[], const int lwork) const
{
  cusolverDnHandle_t rawHandle = handle_->getHandle ();
  using IST = typename CudaValue<Scalar>::type;
  IST* A_raw = reinterpret_cast<IST*> (A);
  const IST* tau_raw = reinterpret_cast<const IST*> (tau);
  IST* work_raw = reinterpret_cast<IST*> (work);

  using impl_type = RawCuSolver<IST>;
  const auto status =
    impl_type::compute_explicit_Q (rawHandle, m, n, k, A_raw, lda,
                                   tau_raw, work_raw, lwork, info_);
  TEUCHOS_ASSERT( status == CUSOLVER_STATUS_SUCCESS );
}

template class CuSolver<double>;
template class CuSolver<float>;
#if defined(HAVE_TPETRATSQR_COMPLEX)
template class CuSolver<std::complex<double>>;
template class CuSolver<std::complex<float>>;
#endif // defined(HAVE_TPETRATSQR_COMPLEX)

} // namespace Impl
} // namespace TSQR

#endif // HAVE_TPETRATSQR_CUBLAS && HAVE_TPETRATSQR_CUSOLVER
