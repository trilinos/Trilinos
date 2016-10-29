#ifndef _KOKKOSKERNELS_SIMPLEUTILS_HPP
#define _KOKKOSKERNELS_SIMPLEUTILS_HPP
#include "Kokkos_Core.hpp"
#include "Kokkos_Atomic.hpp"
#include "impl/Kokkos_Timer.hpp"

#define KOKKOSKERNELS_MACRO_MIN(x,y) ((x) < (y) ? (x) : (y))

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
}
}
}
#endif
