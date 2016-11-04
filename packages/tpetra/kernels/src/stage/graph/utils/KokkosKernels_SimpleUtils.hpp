/*
//@HEADER
// ************************************************************************
//
//          KokkosKernels: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef _KOKKOSKERNELS_SIMPLEUTILS_HPP
#define _KOKKOSKERNELS_SIMPLEUTILS_HPP
#include "Kokkos_Core.hpp"
#include "Kokkos_Atomic.hpp"
#include "impl/Kokkos_Timer.hpp"

#define KOKKOSKERNELS_MACRO_MIN(x,y) ((x) < (y) ? (x) : (y))
#define KOKKOSKERNELS_MACRO_ABS(x) ((x) > (0) ? (x): (-x))

namespace KokkosKernels{

namespace Experimental{

namespace Util{

template <typename view_t>
struct ExclusiveParallelPrefixSum{
  typedef typename view_t::value_type idx;
  view_t array_sum;
  ExclusiveParallelPrefixSum(view_t arr_): array_sum(arr_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, size_t& update, const bool final) const {

    idx val = array_sum(ii);
    if (final) {
      array_sum(ii) = idx (update);
    }
    update += val;
  }
};

template <typename array_type>
struct InclusiveParallelPrefixSum{
  typedef typename array_type::value_type idx;
  array_type array_sum;
  InclusiveParallelPrefixSum(array_type arr_): array_sum(arr_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t ii, size_t& update, const bool final) const {
    update += array_sum(ii);
    if (final) {
      array_sum(ii) = idx (update);
    }
  }
};


/***
 * \brief Function performs the exclusive parallel prefix sum. That is each entry holds the sum
 * until itself.
 * \param num_elements: size of the array
 * \param arr: the array for which the prefix sum will be performed.
 */
template <typename view_t, typename MyExecSpace>
inline void kk_exclusive_parallel_prefix_sum(typename view_t::value_type num_elements, view_t arr){
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_scan( my_exec_space(0, num_elements), ExclusiveParallelPrefixSum<view_t>(arr));
}




/***
 * \brief Function performs the inclusive parallel prefix sum. That is each entry holds the sum
 * until itself including itself.
 * \param num_elements: size of the array
 * \param arr: the array for which the prefix sum will be performed.
 */
template <typename forward_array_type, typename MyExecSpace>
void kk_inclusive_parallel_prefix_sum(typename forward_array_type::value_type num_elements, forward_array_type arr){
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_scan( my_exec_space(0, num_elements), InclusiveParallelPrefixSum<forward_array_type>(arr));
}


template<typename view_type1, typename view_type2, typename eps_type>
struct IsIdenticalFunctor{
  view_type1 view1;
  view_type2 view2;
  eps_type eps;


  IsIdenticalFunctor(view_type1 view1_, view_type2 view2_, eps_type eps_):
    view1(view1_), view2(view2_), eps(eps_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i, size_t &is_equal) const {
    auto val_diff = view1(i) - view2(i);

    if (KOKKOSKERNELS_MACRO_ABS (val_diff) > eps) {
      is_equal+=1;
    }
  }
};

template <typename view_type1, typename view_type2, typename eps_type, typename MyExecSpace>
bool kk_is_identical_view(view_type1 view1, view_type2 view2, eps_type eps){

  if (view1.dimension_0() != view2.dimension_0()){
    return false;
  }

  size_t num_elements = view1.dimension_0();

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  size_t issame = 0;
  Kokkos::parallel_reduce( my_exec_space(0,num_elements),
      IsIdenticalFunctor<view_type1, view_type2, eps_type>(view1, view2, eps), issame);
  MyExecSpace::fence();
  if (issame > 0){
    return false;
  }
  else {
    return true;
  }
}

}
}
}
#endif
