// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tsqr_Impl_CuSolverHandle.hpp"

#ifdef HAVE_TPETRATSQR_CUSOLVER
#include "Kokkos_Core.hpp"
#include "Teuchos_Assert.hpp"

namespace TSQR {
namespace Impl {

cusolverDnHandle_t cuSolverRawHandle_ = nullptr;

CuSolverHandle::CuSolverHandle (cusolverDnHandle_t handle) :
  handle_ (handle)
{}

cusolverDnHandle_t
CuSolverHandle::getHandle () const {
  return handle_;
}

std::shared_ptr<CuSolverHandle> getCuSolverHandleSingleton ()
{
  static std::shared_ptr<CuSolverHandle> singleton_;
  if (singleton_.get () == nullptr) {
    auto finalizer = [] () {
      if (cuSolverRawHandle_ != nullptr) {
        (void) cusolverDnDestroy (cuSolverRawHandle_);
        cuSolverRawHandle_ = nullptr;
      }
    };
    Kokkos::push_finalize_hook (finalizer);
    auto status = cusolverDnCreate (&cuSolverRawHandle_);
    TEUCHOS_ASSERT( status == CUSOLVER_STATUS_SUCCESS );

    singleton_ = std::shared_ptr<CuSolverHandle>
      (new CuSolverHandle (cuSolverRawHandle_));
  }
  TEUCHOS_ASSERT( cuSolverRawHandle_ != nullptr );
  TEUCHOS_ASSERT( singleton_.get () != nullptr );
  return singleton_;
}

} // namespace Impl
} // namespace TSQR
#endif // HAVE_TPETRATSQR_CUSOLVER
