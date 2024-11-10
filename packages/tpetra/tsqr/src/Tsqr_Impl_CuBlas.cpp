// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tsqr_Impl_CuBlas.hpp"
#if defined(HAVE_TPETRATSQR_CUBLAS)
#include "Tsqr_Impl_CuBlasHandle.hpp"
#include "Tsqr_Impl_CuTypes.hpp"
#include "Teuchos_Assert.hpp"

namespace TSQR {
namespace Impl {

template<class T>
class RawCuBlas {};

template<>
class RawCuBlas<double> {
public:
  using impl_scalar_type = double;

  static cublasStatus_t
  gemm (cublasHandle_t handle,
        cublasOperation_t transa,
        cublasOperation_t transb,
        const int m, const int n, const int k,
        const impl_scalar_type* alpha,
        const impl_scalar_type* A, const int lda,
        const impl_scalar_type* B, const int ldb,
        const impl_scalar_type* beta,
        impl_scalar_type* C, const int ldc)
  {
    return cublasDgemm (handle, transa, transb, m, n, k,
                        alpha, A, lda, B, ldb, beta, C, ldc);
  }
};

template<>
class RawCuBlas<float> {
public:
  using impl_scalar_type = float;

  static cublasStatus_t
  gemm (cublasHandle_t handle,
        cublasOperation_t transa,
        cublasOperation_t transb,
        const int m, const int n, const int k,
        const impl_scalar_type* alpha,
        const impl_scalar_type* A, const int lda,
        const impl_scalar_type* B, const int ldb,
        const impl_scalar_type* beta,
        impl_scalar_type* C, const int ldc)
  {
    return cublasSgemm (handle, transa, transb, m, n, k,
                        alpha, A, lda, B, ldb, beta, C, ldc);
  }
};

#if defined(HAVE_TPETRATSQR_COMPLEX)
template<>
class RawCuBlas<CudaValue<std::complex<double>>::type> {
public:
  using impl_scalar_type = CudaValue<std::complex<double>>::type;

  static cublasStatus_t
  gemm (cublasHandle_t handle,
        cublasOperation_t transa,
        cublasOperation_t transb,
        const int m, const int n, const int k,
        const impl_scalar_type* alpha,
        const impl_scalar_type* A, const int lda,
        const impl_scalar_type* B, const int ldb,
        const impl_scalar_type* beta,
        impl_scalar_type* C, const int ldc)
  {
    return cublasZgemm (handle, transa, transb, m, n, k,
                        alpha, A, lda, B, ldb, beta, C, ldc);
  }
};

template<>
class RawCuBlas<CudaValue<std::complex<float>>::type> {
public:
  using impl_scalar_type = CudaValue<std::complex<float>>::type;

  static cublasStatus_t
  gemm (cublasHandle_t handle,
        cublasOperation_t transa,
        cublasOperation_t transb,
        const int m, const int n, const int k,
        const impl_scalar_type* alpha,
        const impl_scalar_type* A, const int lda,
        const impl_scalar_type* B, const int ldb,
        const impl_scalar_type* beta,
        impl_scalar_type* C, const int ldc)
  {
    return cublasCgemm (handle, transa, transb, m, n, k,
                        alpha, A, lda, B, ldb, beta, C, ldc);
  }
};
#endif // defined(HAVE_TPETRATSQR_COMPLEX)

template<class Scalar>
CuBlas<Scalar>::CuBlas () :
  handle_ (getCuBlasHandleSingleton ()) {}

template<class Scalar>
CuBlas<Scalar>::CuBlas (const std::shared_ptr<CuBlasHandle>& handle) :
  handle_ (handle) {}

template<class Scalar>
void
CuBlas<Scalar>::
gemm (const char transa,
      const char transb,
      const int m, const int n, const int k,
      const Scalar alpha,
      const Scalar* A, const int lda,
      const Scalar* B, const int ldb,
      const Scalar beta,
      Scalar* C, const int ldc)
{
  cublasHandle_t rawHandle = handle_->getHandle ();
  const cublasOperation_t cuTransa = cuBlasTrans (transa);
  const cublasOperation_t cuTransb = cuBlasTrans (transb);

  using IST = typename CudaValue<Scalar>::type;
  const IST alpha_raw = CudaValue<Scalar>::makeValue (alpha);
  const IST* A_raw = reinterpret_cast<const IST*> (A);
  const IST* B_raw = reinterpret_cast<const IST*> (B);
  const IST beta_raw = CudaValue<Scalar>::makeValue (beta);
  IST* C_raw = reinterpret_cast<IST*> (C);

  using impl_type = RawCuBlas<IST>;
  // https://docs.nvidia.com/cuda/cublas/index.html#cublas-lt-t-gt-gemm
  // says that alpha and beta may be host or device pointers.
  const auto status =
    impl_type::gemm (rawHandle, cuTransa, cuTransb,
                     m, n, k,
                     &alpha_raw, A_raw, lda,
                     B_raw, ldb,
                     &beta_raw, C_raw, ldc);
  TEUCHOS_ASSERT( status == CUBLAS_STATUS_SUCCESS );
}

template class CuBlas<double>;
template class CuBlas<float>;
#if defined(HAVE_TPETRATSQR_COMPLEX)
template class CuBlas<std::complex<double>>;
template class CuBlas<std::complex<float>>;
#endif // defined(HAVE_TPETRATSQR_COMPLEX)

} // namespace Impl
} // namespace TSQR

#endif // HAVE_TPETRATSQR_CUBLAS
