#include <Kokkos_Sort.hpp>

#ifndef _KOKKOSKERNELSSORTINGIMP_HPP
#define _KOKKOSKERNELSSORTINGIMP_HPP


namespace KokkosKernels{

namespace Experimental{

namespace KokkosKernelsSorting{


namespace Impl{

template<class KeyViewType, class ValueViewType>
struct DefaultBinOp2D {
  const int max_bins_;
  const double mul_;
  typename KeyViewType::const_value_type range_k;

  typename KeyViewType::const_value_type min_k;

  //Construct BinOp with number of bins, minimum value and maxuimum value
  DefaultBinOp2D(int max_bins__,
      typename KeyViewType::const_value_type min_k_,
      typename KeyViewType::const_value_type max_k_)
  :max_bins_(max_bins__+1),mul_(1.0*max_bins__/(max_k_-min_k_)),range_k(max_k_-min_k_),min_k(min_k_) {}

  //Determine bin index from key value
  template<class KViewType, class VViewType>
  KOKKOS_INLINE_FUNCTION
  int bin(KViewType& keys, VViewType& values,  const int& i) const {
    return int(mul_*(keys(i)-min_k));
  }

  //Return maximum bin index + 1
  KOKKOS_INLINE_FUNCTION
  int max_bins() const {
    return max_bins_;
  }

  //Compare to keys within a bin if true new_val will be put before old_val
  template<class KViewType, class VViewType, typename iType1, typename iType2>
  KOKKOS_INLINE_FUNCTION
  bool operator()(KViewType& keys, VViewType& values, iType1& i1, iType2& i2) const {
    return keys(i1)<keys(i2) || (keys(i1) == keys(i2) && values(i1) < values(i2)) ;
  }
};



template<class KeyViewType, class SecondKeyViewType, class BinSortOp, class ExecutionSpace = typename KeyViewType::execution_space,
         class SizeType = typename KeyViewType::memory_space::size_type>
class BinSort2D {


public:
  template<class ValuesViewType, class PermuteViewType, class CopyOp>
  struct bin_sort_sort_functor {
    typedef ExecutionSpace execution_space;
    typedef typename ValuesViewType::non_const_type values_view_type;
    typedef typename ValuesViewType::const_type const_values_view_type;
    Kokkos::View<typename values_view_type::const_data_type,typename values_view_type::array_layout,
                 typename values_view_type::memory_space,Kokkos::MemoryTraits<Kokkos::RandomAccess> > values;
    values_view_type sorted_values;
    typename PermuteViewType::const_type sort_order;
    bin_sort_sort_functor(const_values_view_type values_, values_view_type  sorted_values_, PermuteViewType sort_order_):
       values(values_),sorted_values(sorted_values_),sort_order(sort_order_) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const int& i)  const {
      //printf("Sort: %i %i\n",i,sort_order(i));
      CopyOp::copy(sorted_values,i,values,sort_order(i));
    }
  };

  typedef ExecutionSpace execution_space;
  typedef BinSortOp bin_op_type;

  struct bin_count_tag {};
  struct bin_offset_tag {};
  struct bin_binning_tag {};
  struct bin_sort_bins_tag {};

public:
  typedef SizeType size_type;
  typedef size_type value_type;

  typedef Kokkos::View<size_type*, execution_space> offset_type;
  typedef Kokkos::View<const int*, execution_space> bin_count_type;


  typedef Kokkos::View<typename KeyViewType::const_data_type,
                       typename KeyViewType::array_layout,
                       typename KeyViewType::memory_space> const_key_view_type;
  typedef Kokkos::View<typename KeyViewType::const_data_type,
                       typename KeyViewType::array_layout,
                       typename KeyViewType::memory_space,
                       Kokkos::MemoryTraits<Kokkos::RandomAccess> > const_rnd_key_view_type;


  typedef Kokkos::View<typename SecondKeyViewType::const_data_type,
                       typename SecondKeyViewType::array_layout,
                       typename SecondKeyViewType::memory_space> const_2ndkey_view_type;
  typedef Kokkos::View<typename SecondKeyViewType::const_data_type,
                       typename SecondKeyViewType::array_layout,
                       typename SecondKeyViewType::memory_space,
                       Kokkos::MemoryTraits<Kokkos::RandomAccess> > const_rnd_2ndkey_view_type;



  typedef typename KeyViewType::non_const_value_type non_const_key_scalar;
  typedef typename KeyViewType::const_value_type     const_key_scalar;

private:
  const_key_view_type keys;
  const_2ndkey_view_type vals;

  const_rnd_key_view_type keys_rnd;
  const_rnd_2ndkey_view_type vals_rnd;
public:
  BinSortOp bin_op;

  offset_type bin_offsets;

  Kokkos::View<int*, ExecutionSpace, Kokkos::MemoryTraits<Kokkos::Atomic> > bin_count_atomic;
  bin_count_type bin_count_const;

  offset_type sort_order;

  bool sort_within_bins;

public:

  // Constructor: takes the keys, the binning_operator and optionally whether to sort within bins (default false)
  BinSort2D(const_key_view_type keys_, const_2ndkey_view_type vals_, BinSortOp bin_op_,
          bool sort_within_bins_ = false)
     :keys(keys_), vals(vals_), keys_rnd(keys_), vals_rnd(vals_), bin_op(bin_op_) {

    bin_count_atomic = Kokkos::View<int*, ExecutionSpace >("Kokkos::SortImpl::BinSortFunctor::bin_count",bin_op.max_bins());
    bin_count_const =  bin_count_atomic;
    bin_offsets =      offset_type("Kokkos::SortImpl::BinSortFunctor::bin_offsets",bin_op.max_bins());
    sort_order =       offset_type("PermutationVector",keys.dimension_0());
    sort_within_bins = sort_within_bins_;
  }

  // Create the permutation vector, the bin_offset array and the bin_count array. Can be called again if keys changed
  void create_permute_vector() {
    Kokkos::parallel_for (Kokkos::RangePolicy<ExecutionSpace,bin_count_tag>    (0,keys.dimension_0()),*this);
    Kokkos::parallel_scan(Kokkos::RangePolicy<ExecutionSpace,bin_offset_tag>   (0,bin_op.max_bins()) ,*this);

    Kokkos::deep_copy(bin_count_atomic,0);
    Kokkos::parallel_for (Kokkos::RangePolicy<ExecutionSpace,bin_binning_tag>  (0,keys.dimension_0()),*this);

    if(sort_within_bins)
      Kokkos::parallel_for (Kokkos::RangePolicy<ExecutionSpace,bin_sort_bins_tag>(0,bin_op.max_bins()) ,*this);
  }

  // Sort a view with respect ot the first dimension using the permutation array
  template<class ValuesViewType>
  void sort(ValuesViewType values) {
    ValuesViewType sorted_values = ValuesViewType("Copy",
           values.dimension_0(),
           values.dimension_1(),
           values.dimension_2(),
           values.dimension_3(),
           values.dimension_4(),
           values.dimension_5(),
           values.dimension_6(),
           values.dimension_7());

    parallel_for(values.dimension_0(),
        bin_sort_sort_functor<ValuesViewType, offset_type,
                              Kokkos::SortImpl::CopyOp<ValuesViewType> >(values,sorted_values,sort_order));

    deep_copy(values,sorted_values);
  }

  // Get the permutation vector
  KOKKOS_INLINE_FUNCTION
  offset_type get_permute_vector() const { return sort_order;}

  // Get the start offsets for each bin
  KOKKOS_INLINE_FUNCTION
  offset_type get_bin_offsets() const { return bin_offsets;}

  // Get the count for each bin
  KOKKOS_INLINE_FUNCTION
  bin_count_type get_bin_count() const {return bin_count_const;}

public:
  KOKKOS_INLINE_FUNCTION
  void operator() (const bin_count_tag& tag, const int& i) const {
    bin_count_atomic(bin_op.bin(keys, vals,i))++;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const bin_offset_tag& tag, const int& i, value_type& offset, const bool& final)  const {
    if(final) {
      bin_offsets(i) = offset;
    }
    offset+=bin_count_const(i);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const bin_binning_tag& tag, const int& i)  const {
    const int bin = bin_op.bin(keys,vals, i);
    const int count = bin_count_atomic(bin)++;

    sort_order(bin_offsets(bin) + count) = i;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const bin_sort_bins_tag& tag, const int&i )  const {
    bool sorted = false;
    int upper_bound = bin_offsets(i)+bin_count_const(i);
    while(!sorted) {
      sorted = true;
      int old_idx = sort_order(bin_offsets(i));
      int new_idx;
      for(int k=bin_offsets(i)+1; k<upper_bound; k++) {
        new_idx = sort_order(k);

        if(!bin_op(keys_rnd,vals_rnd, old_idx,new_idx)) {
          sort_order(k-1) = new_idx;
          sort_order(k) = old_idx;
          sorted = false;
        } else {
          old_idx = new_idx;
        }
      }
      upper_bound--;
    }
  }
};


}
}
}
}
#endif
