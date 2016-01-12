#include "KokkosKernels_SortKeyValue_impl.hpp"

#ifndef _KOKKOSKERNELSSORTING_HPP
#define _KOKKOSKERNELSSORTING_HPP
namespace KokkosKernels{

namespace Experimental{

namespace KokkosKernelsSorting{


template <typename view1,  typename view2, typename MyExecSpace>
void sort_key_value_views (view1 key, view2 value){
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  typedef typename Impl::DefaultBinOp2D<view1, view2> CompType;
  Kokkos::SortImpl::min_max<typename view1::non_const_value_type> val;
  Kokkos::parallel_reduce(my_exec_space(0,key.dimension_0()),Kokkos::SortImpl::min_max_functor<view1>(key),val);
  Impl::BinSort2D<view1,view2, CompType>
      bin_sort(key,value, CompType(key.dimension_0()/2,val.min,val.max),true);
  bin_sort.create_permute_vector();
  bin_sort.sort(key);
  bin_sort.sort(value);
}
}
}
}

#endif
