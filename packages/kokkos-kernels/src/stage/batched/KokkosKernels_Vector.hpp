#ifndef __KOKKOSKERNELS_VECTOR_HPP__
#define __KOKKOSKERNELS_VECTOR_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosKernels_Util.hpp"

namespace KokkosKernels {
  template<typename T>
  class Vector;
}

#include "KokkosKernels_Vector_SIMD.hpp"
#include "KokkosKernels_Vector_AVX256D.hpp"
#include "KokkosKernels_Vector_AVX512D.hpp"

#endif
