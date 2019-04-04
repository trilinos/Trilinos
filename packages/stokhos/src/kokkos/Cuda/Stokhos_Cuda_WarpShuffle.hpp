// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
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
  return KOKKOS_IMPL_CUDA_SHFL_DOWN_MASK(mask, val, delta, width);
}

template<typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar shfl_up(const Scalar &val, const int& delta, const int& width,
	       const int& mask){
  return KOKKOS_IMPL_CUDA_SHFL_UP_MASK(mask, val, delta, width);
}

KOKKOS_INLINE_FUNCTION
void sync_warp(const int& mask) {
  KOKKOS_IMPL_CUDA_SYNCWARP_MASK(mask);
}

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_WARP_SHUFFLE_HPP */
