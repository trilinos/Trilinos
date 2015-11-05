
#include <Kokkos_Core.hpp>
#include <Kokkos_Atomic.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <Kokkos_MemoryTraits.hpp>

#ifndef _KOKKOSKERNELSUTILS_HPP
#define _KOKKOSKERNELSUTILS_HPP


namespace Experimental{
namespace KokkosKernels{
namespace Util{

template <typename forward_map_type, typename reverse_map_type>
struct Reverse_Map_Init{
  typedef typename forward_map_type::value_type forward_type;
  typedef typename reverse_map_type::value_type reverse_type;
  forward_map_type forward_map;
  reverse_map_type reverse_map_xadj;
  Reverse_Map_Init(
      forward_map_type forward_map_,
      reverse_map_type reverse_xadj_):
        forward_map(forward_map_), reverse_map_xadj(reverse_xadj_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const reverse_type &ii) const {
    forward_type fm = forward_map[ii];
    Kokkos::atomic_fetch_add( &(reverse_map_xadj(fm)), 1);
  }

  /*
  KOKKOS_INLINE_FUNCTION
  void operator()(const forward_type ii, size_t& update, const bool final) const {
    update += reverse_map_xadj(ii);
    if (final) {
      reverse_map_xadj(ii) = reverse_type (update);
    }
  }
  */
};





template <typename forward_map_type, typename reverse_map_type>
struct Fill_Reverse_Map{
  typedef typename forward_map_type::value_type forward_type;
  typedef typename reverse_map_type::value_type reverse_type;
  forward_map_type forward_map;
  reverse_map_type reverse_map_xadj;
  reverse_map_type reverse_map_adj;


  Fill_Reverse_Map(
      forward_map_type forward_map_,
      reverse_map_type reverse_map_xadj_,
      reverse_map_type reverse_map_adj):
        forward_map(forward_map_), reverse_map_xadj(reverse_map_xadj_), reverse_map_adj(reverse_map_adj){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const reverse_type &ii) const {
    forward_type c = forward_map[ii];
    const reverse_type future_index = Kokkos::atomic_fetch_add( &(reverse_map_xadj(c - 1)), 1);
    reverse_map_adj(future_index) = ii;
  }
};


template <typename array_type>
struct InclusiveParallelPrefixSum{
  typedef typename array_type::value_type idx;
  array_type array_sum;
  InclusiveParallelPrefixSum(array_type arr_): array_sum(arr_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const idx ii, size_t& update, const bool final) const {
    update += array_sum(ii);
    if (final) {
      array_sum(ii) = idx (update);
    }
  }
};

template <typename forward_array_type, typename MyExecSpace>
void inclusive_parallel_prefix_sum(typename forward_array_type::value_type num_elements, forward_array_type arr){
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_scan( my_exec_space(0, num_elements), InclusiveParallelPrefixSum<forward_array_type>(arr));
}

/**
 * \brief Utility function to obtain a reverse map given a map.
 * Input is a map with the number of elements within the map.
 * forward_map[c] = i, where c is a forward elements and forward_map has a size of num_forward_elements.
 * i is the value that c is mapped in the forward map, and the range of that is num_reverse_elements.
 * Output is the reverse_map_xadj and reverse_map_adj such that,
 * all c, forward_map[c] = i, will appear in  reverse_map_adj[ reverse_map_xadj[i]: reverse_map_xadj[i+1])
 * \param: num_forward_elements: the number of elements in the forward map, the size of the forward map.
 * \param: num_reverse_elements: the number of elements that forward map is mapped to. It is the value of max i.
 * \param: forward_map: input forward_map, where forward_map[c] = i.
 * \param: reverse_map_xadj: reverse map xadj, that is it will hold the beginning and
 * end indices on reverse_map_adj such that all values mapped to i will be [ reverse_map_xadj[i]: reverse_map_xadj[i+1])
 * its size will be num_reverse_elements + 1.
 * \param: reverse_map_adj: reverse map adj, holds the values of reverse maps. Its size will be num_forward_elements.
 *
 */
template <typename forward_array_type, typename reverse_array_type, typename MyExecSpace>
void create_reverse_map(
    const typename reverse_array_type::value_type &num_forward_elements, //num_vertices
    const typename forward_array_type::value_type &num_reverse_elements, //num_colors

    const forward_array_type &forward_map, //vertex to colors
    reverse_array_type &reverse_map_xadj, // colors to vertex xadj
    reverse_array_type &reverse_map_adj){ //colros to vertex adj

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;

  reverse_map_xadj = reverse_array_type("Reverse Map Xadj", num_reverse_elements + 1);
  reverse_map_adj = reverse_array_type(Kokkos::ViewAllocateWithoutInitializing("REVERSE_ADJ"), num_forward_elements);
  reverse_array_type tmp_color_xadj (Kokkos::ViewAllocateWithoutInitializing("TMP_REVERSE_XADJ"), num_reverse_elements + 1);

  Reverse_Map_Init<forward_array_type, reverse_array_type> rmi(forward_map, reverse_map_xadj);
  Kokkos::parallel_for (my_exec_space (0, num_forward_elements) , rmi);
  inclusive_parallel_prefix_sum<reverse_array_type, MyExecSpace>(num_reverse_elements + 1, reverse_map_xadj);
  //Kokkos::parallel_scan (my_exec_space (0, num_reverse_elements + 1) , rmi);
  Kokkos::deep_copy (tmp_color_xadj, reverse_map_xadj);
  Fill_Reverse_Map<forward_array_type, reverse_array_type> frm (forward_map, tmp_color_xadj, reverse_map_adj);
  Kokkos::parallel_for (my_exec_space (0, num_forward_elements) , frm);
}


template <typename value_array_type, typename idx_array_type>
struct PermuteVector{
  typedef typename idx_array_type::value_type idx;
  value_array_type old_vector;
  value_array_type new_vector;
  idx_array_type old_to_new_mapping;
  PermuteVector(
      value_array_type old_vector_,
      value_array_type new_vector_,
      idx_array_type old_to_new_mapping_):
        old_vector(old_vector_), new_vector(new_vector_),old_to_new_mapping(old_to_new_mapping_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const idx &ii) const {

//    /std::cout << "ii:" << ii << " old_to_new_mapping[ii]:" << old_to_new_mapping[ii] << std::endl;
    new_vector[old_to_new_mapping[ii]] = old_vector[ii];
  }
};

template <typename value_array_type, typename idx_array_type, typename MyExecSpace>
void permute_vector(
    typename idx_array_type::value_type num_elements,
    idx_array_type &old_to_new_index_map,
    value_array_type &old_vector,
    value_array_type &new_vector
    ){
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;

  Kokkos::parallel_for( my_exec_space(0,num_elements),
      PermuteVector<value_array_type, idx_array_type>(old_vector, new_vector, old_to_new_index_map));

}

template <typename value_array_type>
struct ZeroVector{
  value_array_type myview;
  ZeroVector(value_array_type myview_): myview(myview_){}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i ) const {
    myview(i) = 0;
  }
};

template <typename value_array_type, typename MyExecSpace>
void zero_vector(
    typename value_array_type::value_type num_elements,
    value_array_type &vector
    ){
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_for( my_exec_space(0,num_elements), ZeroVector<value_array_type>(vector));

}

}
}
}
#endif
