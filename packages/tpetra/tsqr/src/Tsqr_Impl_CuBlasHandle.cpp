#include "Tsqr_Impl_CuBlasHandle.hpp"

#ifdef HAVE_TPETRATSQR_CUBLAS
#include "Kokkos_Core.hpp"
#include "Teuchos_Assert.hpp"
#include <cublas_v2.h>

namespace TSQR {
namespace Impl {

cublasHandle_t cuBlasRawHandle_ = nullptr;

CuBlasHandle::CuBlasHandle (void* handle) :
  handle_ (handle)
{}

CuBlasHandle CuBlasHandle::getSingleton ()
{
  static int called_before = 0;
  if (called_before == 0) {
    auto finalizer = [] () {
      if (cuBlasRawHandle_ != nullptr) {
        (void) cublasDestroy (cuBlasRawHandle_);
        cuBlasRawHandle_ = nullptr;
      }
    };
    Kokkos::push_finalize_hook (finalizer);
    auto status = cublasCreate (&cuBlasRawHandle_);
    TEUCHOS_ASSERT( status == CUBLAS_STATUS_SUCCESS );
    called_before = 1;
  }
  TEUCHOS_ASSERT( cuBlasRawHandle_ != nullptr );
  return CuBlasHandle (cuBlasRawHandle_);
}

} // namespace Impl
} // namespace TSQR
#endif // HAVE_TPETRATSQR_CUBLAS
