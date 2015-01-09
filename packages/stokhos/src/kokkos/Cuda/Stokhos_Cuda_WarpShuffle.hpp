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
  return val;
}

template<typename Scalar>
KOKKOS_INLINE_FUNCTION
Scalar shfl_up(const Scalar &val, const int& delta, const int& width){
  return val;
}


#if HAVE_CUDA_SHUFFLE

KOKKOS_INLINE_FUNCTION
unsigned int shfl_down(
  const unsigned int &val, const int& delta, const int& width) {
  unsigned int tmp1 = val;
  int tmp = *reinterpret_cast<int*>(&tmp1);
  tmp = __shfl_down(tmp,delta,width);
  return *reinterpret_cast<unsigned int*>(&tmp);
}

KOKKOS_INLINE_FUNCTION
int shfl_down(const int &val, const int& delta, const int& width) {
  return __shfl_down(val,delta,width);
}

KOKKOS_INLINE_FUNCTION
float shfl_down(const float &val, const int& delta, const int& width) {
  return __shfl_down(val,delta,width);
}

KOKKOS_INLINE_FUNCTION
double shfl_down(const double &val, const int& delta, const int& width) {
  int lo = __double2loint(val);
  int hi = __double2hiint(val);
  lo = __shfl_down(lo,delta,width);
  hi = __shfl_down(hi,delta,width);
  return __hiloint2double(hi,lo);
}

KOKKOS_INLINE_FUNCTION
unsigned int shfl_up(
  const unsigned int &val, const int& delta, const int& width) {
  unsigned int tmp1 = val;
  int tmp = *reinterpret_cast<int*>(&tmp1);
  tmp = __shfl_up(tmp,delta,width);
  return *reinterpret_cast<unsigned int*>(&tmp);
}

KOKKOS_INLINE_FUNCTION
int shfl_up(const int &val, const int& delta, const int& width) {
  return __shfl_up(val,delta,width);
}

KOKKOS_INLINE_FUNCTION
float shfl_up(const float &val, const int& delta, const int& width) {
  return __shfl_up(val,delta,width);
}

KOKKOS_INLINE_FUNCTION
double shfl_up(const double &val, const int& delta, const int& width) {
  int lo = __double2loint(val);
  int hi = __double2hiint(val);
  lo = __shfl_up(lo,delta,width);
  hi = __shfl_up(hi,delta,width);
  return __hiloint2double(hi,lo);
}

#endif // #if HAVE_CUDA_SHUFFLE

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_WARP_SHUFFLE_HPP */
