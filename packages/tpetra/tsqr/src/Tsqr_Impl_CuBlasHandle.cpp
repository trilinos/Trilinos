// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tsqr_Impl_CuBlasHandle.hpp"

#ifdef HAVE_TPETRATSQR_CUBLAS
#include "Kokkos_Core.hpp"
#include "Teuchos_Assert.hpp"

namespace TSQR {
namespace Impl {

cublasHandle_t cuBlasRawHandle_ = nullptr;

CuBlasHandle::CuBlasHandle (cublasHandle_t handle) :
  handle_ (handle)
{}

cublasHandle_t
CuBlasHandle::getHandle () const {
  return handle_;
}

std::shared_ptr<CuBlasHandle> getCuBlasHandleSingleton ()
{
  static std::shared_ptr<CuBlasHandle> singleton_;
  if (singleton_.get () == nullptr) {
    auto finalizer = [] () {
      if (cuBlasRawHandle_ != nullptr) {
        (void) cublasDestroy (cuBlasRawHandle_);
        cuBlasRawHandle_ = nullptr;
      }
    };
    Kokkos::push_finalize_hook (finalizer);
    auto status = cublasCreate (&cuBlasRawHandle_);
    TEUCHOS_ASSERT( status == CUBLAS_STATUS_SUCCESS );

    singleton_ = std::shared_ptr<CuBlasHandle>
      (new CuBlasHandle (cuBlasRawHandle_));
  }
  TEUCHOS_ASSERT( cuBlasRawHandle_ != nullptr );
  TEUCHOS_ASSERT( singleton_.get () != nullptr );
  return singleton_;
}

} // namespace Impl
} // namespace TSQR
#endif // HAVE_TPETRATSQR_CUBLAS
