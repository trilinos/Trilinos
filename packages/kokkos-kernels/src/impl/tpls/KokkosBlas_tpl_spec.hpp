#ifndef KOKKOSBLAS_TPL_SPEC_HPP_
#define KOKKOSBLAS_TPL_SPEC_HPP_

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUBLAS
#include "cuda_runtime.h"
#include "cublas_v2.h"

namespace KokkosBlas {
namespace Impl {

struct CudaBlasSingleton {
  cublasHandle_t handle;

  CudaBlasSingleton();

  static CudaBlasSingleton & singleton();
};

}
}
#endif

#endif