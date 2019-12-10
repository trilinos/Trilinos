#include "Tsqr_Impl_CuSolverHandle.hpp"

#ifdef HAVE_TPETRATSQR_CUSOLVER
#include "Kokkos_Core.hpp"
#include "Teuchos_Assert.hpp"
#include <cusolverDn.h>

namespace TSQR {
namespace Impl {

cusolverDnHandle_t cuSolverRawHandle_ = nullptr;

CuSolverHandle::CuSolverHandle (void* handle) :
  handle_ (handle)
{}

CuSolverHandle CuSolverHandle::getSingleton ()
{
  static int called_before = 0;
  if (called_before == 0) {
    auto finalizer = [] () {
      if (cuSolverRawHandle_ != nullptr) {
        (void) cusolverDnDestroy (cuSolverRawHandle_);
        cuSolverRawHandle_ = nullptr;
      }
    };
    Kokkos::push_finalize_hook (finalizer);
    auto status = cusolverDnCreate (&cuSolverRawHandle_);
    TEUCHOS_ASSERT( status == CUSOLVER_STATUS_SUCCESS );
    called_before = 1;
  }
  TEUCHOS_ASSERT( cuSolverRawHandle_ != nullptr );
  return CuSolverHandle (cuSolverRawHandle_);
}

} // namespace Impl
} // namespace TSQR
#endif // HAVE_TPETRATSQR_CUSOLVER
