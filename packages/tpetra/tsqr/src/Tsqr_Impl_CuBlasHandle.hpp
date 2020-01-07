#ifndef TSQR_IMPL_CUBLASHANDLE_HPP
#define TSQR_IMPL_CUBLASHANDLE_HPP

#include "TpetraTSQR_config.h"
#ifdef HAVE_TPETRATSQR_CUBLAS

namespace TSQR {
namespace Impl {

class CuBlasHandle {
private:
  // This is actually a cublasHandle_t, which is a pointer type.
  void* handle_ {nullptr};

  CuBlasHandle (void* handle);

public:
  static CuBlasHandle getSingleton ();

  // This is not really encapsulation, because the "handle" type is
  // just a pointer.  However, it lets us define cuBlas wrapper
  // functions without needing to make them friends of CuBlasHandle.
  void* getHandle () const {
    return handle_;
  }
};

} // namespace Impl
} // namespace TSQR

#endif // HAVE_TPETRATSQR_CUBLAS

#endif // TSQR_IMPL_CUBLASHANDLE_HPP
