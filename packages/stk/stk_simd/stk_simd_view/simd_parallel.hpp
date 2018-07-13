// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
 // 
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
 // 
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef STK_SIMD_PARALLEL_H
#define STK_SIMD_PARALLEL_H

#include <stk_simd_view/simd_view.hpp>
#include <stk_simd_view/has_typedef.hpp>
#include <stk_simd_view/simd_index.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <typeinfo>

namespace stk {
namespace simd {

template <typename Func>
inline constexpr bool is_gpu() {
  typedef typename
    Kokkos::Impl::FunctorPolicyExecutionSpace<Func, void>::execution_space execution_space;
#ifdef KOKKOS_ENABLE_CUDA
  return std::is_same<execution_space, Kokkos::Cuda>::value;
#else
  return false;
#endif
}


template <typename T, typename Func>
STK_INLINE
int get_simd_loop_size(int N) {
  return is_gpu<Func>() ? N : simd_pad<T>(N) / SimdSizeTraits<T>::simd_width;
}


#ifdef KOKKOS_ENABLE_CUDA

template <typename T=double, typename Func, typename PolicyTag>
inline void
parallel_for(std::string forName, const Kokkos::RangePolicy<PolicyTag>& range, const Func& func) {
  assert(range.begin() == 0);
  const int simdLoopSize = get_simd_loop_size<T, Func>(range.end());
  Kokkos::RangePolicy<PolicyTag> simdRange(0, simdLoopSize);
  Kokkos::parallel_for(forName, simdRange, KOKKOS_LAMBDA(const PolicyTag& tag, const int i) {
    func( tag, simd::DeviceIndex(i) );
  });
}

template <typename T=double, typename Func>
inline void
parallel_for(std::string forName, const Kokkos::RangePolicy<void>& range, const Func& func) {
  assert(range.begin() == 0); // only supporting ranges starting at 0
  const int simdLoopSize = get_simd_loop_size<T, Func>(range.end());
  Kokkos::RangePolicy<void> simdRange(0, simdLoopSize);
  Kokkos::parallel_for(forName, simdRange, KOKKOS_LAMBDA(const int i) {
    func( simd::DeviceIndex(i) );
  });
}

template <typename T=double, typename Func>
inline void
parallel_for(std::string forName, int N, const Func& func) {
  parallel_for<T>(forName, Kokkos::RangePolicy<void>(0,N), func);
}

#else

template <typename T, typename Func, typename PolicyTag>
inline void
parallel_for(std::string forName, const Kokkos::RangePolicy<PolicyTag>& range, const Func& func) {
  assert(range.begin() == 0);
  const int simdLoopSize = get_simd_loop_size<T, Func>(range.end());
  Kokkos::RangePolicy<PolicyTag> simdRange(0, simdLoopSize);
  Kokkos::parallel_for(forName, simdRange, KOKKOS_LAMBDA(const PolicyTag& tag, const int i) {
    func( tag, simd::DeviceIndex(i) );
  });
}

template <typename T, typename Func>
inline void
parallel_for(std::string forName, const Kokkos::RangePolicy<void>& range, const Func& func) {
  assert(range.begin() == 0); // only supporting ranges starting at 0
  const int simdLoopSize = get_simd_loop_size<T, Func>(range.end());
  Kokkos::RangePolicy<void> simdRange(0, simdLoopSize);
  Kokkos::parallel_for(forName, simdRange, KOKKOS_LAMBDA(const int i) {
    func( simd::DeviceIndex(i) );
  });
}

template <typename T, typename Func>
inline void
parallel_for(std::string forName, int N, const Func& func) {
  parallel_for<T>(forName, Kokkos::RangePolicy<void>(0,N), func);
}

// specialize/overoad for the case where no type is specified: default to double

template <typename Func, typename PolicyTag>
inline void
parallel_for(std::string forName, const Kokkos::RangePolicy<PolicyTag>& range, const Func& func) {
  parallel_for<double>(forName, range, func);
}

template <typename Func>
inline void
parallel_for(std::string forName, int N, const Func& func) {
  parallel_for<double>(forName, Kokkos::RangePolicy<void>(0,N), func);
}

#endif


// ThreadTeamRange
template <typename Func>
KOKKOS_INLINE_FUNCTION void
parallel_for(const Kokkos::TeamPolicy<>::member_type& thread, int N, const Func& func) {
  const int simdLoopSize = get_simd_loop_size<double, Func>(N);
  Kokkos::parallel_for(Kokkos::TeamThreadRange(thread, simdLoopSize), [&](const int i) {
    func( simd::DeviceIndex(i) );
  });
}


template <typename T, typename Func>
STK_INLINE void
for_each(int N, const Func& func) {
  const int simdLoopSize = get_simd_loop_size<T, Func>(N);
  for (int i=0; i < simdLoopSize; ++i) {
    func( simd::Index(i) );
  }
}

template <typename Func>
STK_INLINE void
for_each(int N, const Func& func) {
  const int simdLoopSize = get_simd_loop_size<double, Func>(N);
  for (int i=0; i < simdLoopSize; ++i) {
    func( DeviceIndex(i) );
  }
}


template <typename Func, typename RealType, typename PolicyTag>
struct SimdReduct {
  SimdReduct(const Func& func_, const int& simdLoopSizeWithoutRemainder_, const RealType& remainderMask_)
    : func(func_), simdLoopSizeWithoutRemainder(simdLoopSizeWithoutRemainder_), remainderMask(remainderMask_) {
  }

  typedef RealType value_type;
  
  static void init(PolicyTag, value_type& initialVal) {
    initialVal = 0.0;
  }

  STK_INLINE void operator() (PolicyTag tag, const int n, value_type& reductee) const {
    if (n != simdLoopSizeWithoutRemainder) {
      func(tag, simd::DeviceIndex(n), reductee);
    } else {
      value_type tmp(0.0);
      func(tag, simd::DeviceIndex(n), tmp);
      reductee += tmp*remainderMask;
    }
  }
  
 private:
  const Func func;
  const int simdLoopSizeWithoutRemainder;
  const value_type remainderMask;
};


template <typename Func, typename RealType>
struct SimdReduct<Func, RealType, void> {
  SimdReduct(const Func& func_, const int& simdLoopSizeWithoutRemainder_, const RealType& remainderMask_)
    : func(func_), simdLoopSizeWithoutRemainder(simdLoopSizeWithoutRemainder_), remainderMask(remainderMask_) {}

  typedef RealType value_type;
  
  static void init(value_type& initialVal) {
    initialVal = 0.0;
  }

  STK_INLINE void operator() (const int n, value_type& reductee) const {
    if (n != simdLoopSizeWithoutRemainder) {
      func(simd::DeviceIndex(n), reductee);
    } else {
      value_type tmp(0.0);
      func(simd::DeviceIndex(n), tmp);
      reductee += tmp*remainderMask;
    }
  }
  
 private:
  const Func func;
  const int simdLoopSizeWithoutRemainder;
  const value_type remainderMask;
};

template <typename Func, typename PolicyTag, typename ScalarType>
typename std::enable_if<!is_gpu<Func>()>::type
parallel_reduce_sum(std::string reduceName, Kokkos::RangePolicy<PolicyTag> policy, const Func& func, ScalarType& reductee) {
  typedef typename stk::Traits<ScalarType>::simd_type SimdType;
  constexpr int simdLength = stk::Traits<SimdType>::length;
  SimdType reducteeSimd(reductee);
  const int N = policy.end();
  const int simdLoopSize = get_simd_loop_size<ScalarType, Func>(N);
  const int simdLoopSizeWithoutRemainder = simdLoopSize*simdLength==N ? simdLoopSize : simdLoopSize-1;
  ScalarType ones[simdLength];
  for (int i=0; i<simdLength; ++i) ones[i] = 1.0;
  SimdType remainderMask = stk::simd::load_part(ones, N%simdLength);
  SimdReduct<Func,SimdType,PolicyTag> reductor(func, simdLoopSizeWithoutRemainder, remainderMask);
  Kokkos::RangePolicy<PolicyTag> policySimd(0, simdLoopSize);
  Kokkos::parallel_reduce(reduceName, policySimd, reductor, reducteeSimd);
  reductee = stk::simd::reduce_sum(reducteeSimd);
}

template <typename Func, typename ScalarType>
typename std::enable_if<!is_gpu<Func>()>::type
parallel_reduce_sum(std::string reduceName, const int N, const Func& func, ScalarType& reductee) {
  parallel_reduce_sum(reduceName, Kokkos::RangePolicy<void>(0,N), func, reductee);
}

template <typename Func, typename PolicyTag, typename ScalarType>
typename std::enable_if<is_gpu<Func>()>::type
parallel_reduce_sum(std::string reduceName, Kokkos::RangePolicy<PolicyTag> policy, const Func& func, ScalarType& reductee) {
  Kokkos::parallel_reduce(reduceName, policy, KOKKOS_LAMBDA(PolicyTag tag, const int i, ScalarType& v) {
    func(tag, i, v);
  }, 
  reductee);
}

template <typename Func, typename ScalarType>
typename std::enable_if<is_gpu<Func>()>::type
parallel_reduce_sum(std::string reduceName, const int N, const Func& func, ScalarType& reductee) {
  Kokkos::parallel_reduce(reduceName, Kokkos::RangePolicy<void>(0,N), KOKKOS_LAMBDA(const int i, ScalarType& v) {
    func(i, v);
  }, 
  reductee);
}


template <typename Func> inline
double reduce_sum_each(int N, const Func& func) {
  typedef typename DeviceTraits<double>::simd_type SimdDouble;
  constexpr int length = stk::Traits<SimdDouble>::length;
  SimdDouble sum = 0.0;

  const int remainder = N%length;
  const int simdLoopSize = get_simd_loop_size<double, Func>(N);
  if (remainder==0) {
    for (int i=0; i < simdLoopSize; ++i) { sum += func( DeviceIndex(i) ); }
  } else {
    for (int i=0; i < simdLoopSize-1; ++i) { sum += func( DeviceIndex(i) ); }
    SimdDouble mask=0;
    for (int i=0; i < remainder; ++i) { stk::simd::set_data(mask, i, 1.0); }
    sum += mask*func( DeviceIndex(simdLoopSize-1) );
  }
  return stk::simd::reduce_sum(sum);
}

}} // stk::simd

#endif
