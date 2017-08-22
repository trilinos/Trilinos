#ifndef __KOKKOSBATCHED_VECTOR_HPP__
#define __KOKKOSBATCHED_VECTOR_HPP__

/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosBatched_Util.hpp"

namespace KokkosBatched {
  namespace Experimental {
    template<typename T>
    class Vector;
  }
}

#include "KokkosBatched_Vector_SIMD.hpp"
#include "KokkosBatched_Vector_AVX256D.hpp"
#include "KokkosBatched_Vector_AVX512D.hpp"

#endif
