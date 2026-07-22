// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_CUDA_WARP_SHUFFLE_HPP
#define STOKHOS_CUDA_WARP_SHUFFLE_HPP

#include "Kokkos_Core.hpp"

#ifdef __CUDA_ARCH__
#  if (__CUDA_ARCH__ >= 300)
#    define HAVE_CUDA_SHUFFLE 1
#  else
#    define HAVE_CUDA_SHUFFLE 0
#  endif
#else
#  define HAVE_CUDA_SHUFFLE 0
#endif

namespace Stokhos {

template<typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar shfl_down(const Scalar &val, const int& delta, const int& width){
  return Kokkos::shfl_down(val, delta, width);
}

template<typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar shfl_up(const Scalar &val, const int& delta, const int& width){
  return Kokkos::shfl_up(val, delta, width);
}

template<typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar shfl_down(const Scalar &val, const int& delta, const int& width,
		 const int& mask){
#ifdef __CUDA_ARCH__
  return __shfl_down_sync(mask, val, delta, width);
#else
  (void) delta;
  (void) width;
  (void) mask;
  return val;
#endif
}

template<typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar shfl_up(const Scalar &val, const int& delta, const int& width,
	       const int& mask){
#ifdef __CUDA_ARCH__
  return __shfl_up_sync(mask, val, delta, width);
#else
  (void) delta;
  (void) width;
  (void) mask;
  return val;
#endif
}

KOKKOS_INLINE_FUNCTION
void sync_warp(const int& mask) {
#ifdef __CUDA_ARCH__
  __syncwarp(mask);
#endif
}

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_WARP_SHUFFLE_HPP */
