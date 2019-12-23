#ifndef TSQR_IMPL_CUSOLVERHANDLE_HPP
#define TSQR_IMPL_CUSOLVERHANDLE_HPP

#include "TpetraTSQR_config.h"
#ifdef HAVE_TPETRATSQR_CUSOLVER

namespace TSQR {
namespace Impl {

class CuSolverHandle {
private:
  // This is actually a cusolverDnHandle_t, which is a pointer type.
  void* handle_ {nullptr};

  CuSolverHandle (void* handle);

public:
  static CuSolverHandle getSingleton ();

  // This is not really encapsulation, because the "handle" type is
  // just a pointer.  However, it lets us define cuSolver wrapper
  // functions without needing to make them friends of CuSolverHandle.
  void* getHandle () const {
    return handle_;
  }
};

} // namespace Impl
} // namespace TSQR

#endif // HAVE_TPETRATSQR_CUSOLVER

#endif // TSQR_IMPL_CUSOLVERHANDLE_HPP
