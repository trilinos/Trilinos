
#include <Kokkos_Core.hpp>
#include <Kokkos_Atomic.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <KokkosKernels_SortKeyValue.hpp>
#include <iostream>

#include <Kokkos_UnorderedMap.hpp>

#ifndef _KOKKOSKERNELSUTILS_HPP
#define _KOKKOSKERNELSUTILS_HPP

namespace KokkosKernels{

namespace Experimental{

namespace Util{

template <typename idx, typename ExecutionSpace>
void get_suggested_vector_team_size(
    int max_allowed_team_size,
    int &suggested_vector_size_,
    int &suggested_team_size_,
    idx nr, idx nnz){

#if defined( KOKKOS_HAVE_SERIAL )
  if (Kokkos::Impl::is_same< Kokkos::Serial , ExecutionSpace >::value){
    suggested_vector_size_ =  1;
    suggested_team_size_ = max_allowed_team_size;
    return;
  }
#endif

#if defined( KOKKOS_HAVE_PTHREAD )
  if (Kokkos::Impl::is_same< Kokkos::Threads , ExecutionSpace >::value){
    suggested_vector_size_ =  1;
    suggested_team_size_ =  max_allowed_team_size;
    return;
  }
#endif

#if defined( KOKKOS_HAVE_OPENMP )
  if (Kokkos::Impl::is_same< Kokkos::OpenMP, ExecutionSpace >::value){
    suggested_vector_size_ =  1;
    suggested_team_size_ = max_allowed_team_size;
  }
#endif

#if defined( KOKKOS_HAVE_CUDA )
  if (Kokkos::Impl::is_same<Kokkos::Cuda, ExecutionSpace >::value){

    suggested_vector_size_ = nnz / double (nr) + 0.5;

    if (suggested_vector_size_ <= 3){
      suggested_vector_size_ = 2;
    }
    else if (suggested_vector_size_ <= 6){
      suggested_vector_size_ = 4;
    }
    else if (suggested_vector_size_ <= 12){
      suggested_vector_size_ = 8;
    }
    else if (suggested_vector_size_ <= 24){
      suggested_vector_size_ = 16;
    }
    else {
      suggested_vector_size_ = 32;
    }
    if (max_allowed_team_size < 32){
      std::cerr << "max_allowed_team_size:" << max_allowed_team_size << std::endl;
    }

    suggested_team_size_ = max_allowed_team_size / suggested_vector_size_;
  }
#endif

#if defined( KOKKOS_HAVE_QTHREAD)
  if (Kokkos::Impl::is_same< Kokkos::Qthread, ExecutionSpace >::value){
    suggested_vector_size_ = 1;
    suggested_team_size_ = max_allowed_team_size;
  }
#endif

}


template <typename idx_array_type,
          typename idx_edge_array_type,
          typename idx_out_edge_array_type,
          typename team_member>
struct FillSymmetricEdges{
  typedef typename idx_array_type::value_type idx;
  idx num_rows;
  idx nnz;
  idx_array_type xadj;
  idx_edge_array_type adj;

  idx_out_edge_array_type srcs;
  idx_out_edge_array_type dsts;

  FillSymmetricEdges(
    typename idx_array_type::value_type num_rows_,
    idx_array_type xadj_,
    idx_edge_array_type adj_,

    idx_out_edge_array_type srcs_,
    idx_out_edge_array_type dsts_
    ):num_rows(num_rows_),nnz(adj_.dimension_0()), xadj(xadj_), adj(adj_), srcs(srcs_), dsts(dsts_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member & teamMember) const {
    idx ii = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
    if (ii >= num_rows) return;
    idx row_begin = xadj[ii];
    idx row_end = xadj[ii + 1];

    Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(teamMember, row_end - row_begin),
        [&] (idx i) {
      idx adjind = i + row_begin;
      idx colIndex = adj[adjind];
      if (colIndex < num_rows){
        srcs[adjind] = ii + 1;
        dsts[adjind] = colIndex + 1;
        if (colIndex != ii){
          srcs[adjind + nnz] = colIndex + 1;
          dsts[adjind + nnz] = ii + 1;
        }
      }

    });

  }
};


template <typename idx_array_type>
void print_1Dview(idx_array_type view, bool print_all = false){

  typedef typename idx_array_type::HostMirror host_type;
  typedef typename idx_array_type::size_type idx;
  host_type host_view = Kokkos::create_mirror_view (view);
  Kokkos::deep_copy (host_view , view);
  idx nr = host_view.dimension_0();
  if (!print_all){

    idx n = nr/10 > 10 ? 10:nr/10;
    for (idx i = 0; i < n; ++i){
      std::cout << host_view(i) << " ";
    }
    std::cout << "... ... ... ";

    for (idx i = nr-n; i < nr; ++i){
      std::cout << host_view(i) << " ";
    }
    std::cout << std::endl;
  }
  else {
    for (idx i = 0; i < nr; ++i){
      std::cout << host_view(i) << " ";
    }
    std::cout << std::endl;
  }
}


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

template <typename array_type>
struct ExclusiveParallelPrefixSum{
  typedef typename array_type::value_type idx;
  array_type array_sum;
  ExclusiveParallelPrefixSum(array_type arr_): array_sum(arr_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const idx ii, size_t& update, const bool final) const {

    idx val = array_sum(ii);
    if (final) {
      array_sum(ii) = idx (update);
    }
    update += val;
  }
};

template <typename forward_array_type, typename MyExecSpace>
void exclusive_parallel_prefix_sum(typename forward_array_type::value_type num_elements, forward_array_type arr){
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_scan( my_exec_space(0, num_elements), ExclusiveParallelPrefixSum<forward_array_type>(arr));
}

template <typename array_type>
struct PropogataMaxValstoZeros{
  typedef typename array_type::value_type idx;
  array_type array_sum;
  PropogataMaxValstoZeros(array_type arr_): array_sum(arr_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const idx ii, idx& update, const bool final) const {

    idx value = array_sum(ii);
    if (value != 0) {
      update = value;
    }
    else if (final ){
      array_sum(ii) = idx (update);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile       idx & update
           , volatile const idx & input ) const {
    if (input > update) update = input;
  }


};

template <typename forward_array_type, typename MyExecSpace>
void remove_zeros_in_xadj_vector(typename forward_array_type::value_type num_elements, forward_array_type arr){
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_scan( my_exec_space(0, num_elements), PropogataMaxValstoZeros<forward_array_type>(arr));
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
  MyExecSpace::fence();
  Kokkos::deep_copy (tmp_color_xadj, reverse_map_xadj);
  MyExecSpace::fence();
  Fill_Reverse_Map<forward_array_type, reverse_array_type> frm (forward_map, tmp_color_xadj, reverse_map_adj);
  Kokkos::parallel_for (my_exec_space (0, num_forward_elements) , frm);
}


template <typename value_array_type, typename out_value_array_type, typename idx_array_type>
struct PermuteVector{
  typedef typename idx_array_type::value_type idx;
  value_array_type old_vector;
  out_value_array_type new_vector;
  idx_array_type old_to_new_mapping;
  idx mapping_size;
  PermuteVector(
      value_array_type old_vector_,
      out_value_array_type new_vector_,
      idx_array_type old_to_new_mapping_):
        old_vector(old_vector_), new_vector(new_vector_),old_to_new_mapping(old_to_new_mapping_), mapping_size(old_to_new_mapping_.dimension_0()){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const idx &ii) const {

    idx mapping = ii;
    if (ii < mapping_size) mapping = old_to_new_mapping[ii];
    new_vector[mapping] = old_vector[ii];
  }
};

template <typename value_array_type, typename out_value_array_type, typename idx_array_type, typename MyExecSpace>
void permute_vector(
    typename idx_array_type::value_type num_elements,
    idx_array_type &old_to_new_index_map,
    value_array_type &old_vector,
    out_value_array_type &new_vector
    ){
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;

  Kokkos::parallel_for( my_exec_space(0,num_elements),
      PermuteVector<value_array_type, out_value_array_type, idx_array_type>(old_vector, new_vector, old_to_new_index_map));

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


template <typename v1, typename v2, typename v3>
struct MarkDuplicateSortedKeyValuePairs{
  v1 keys;
  v2 vals;
  v3 prefix_sum;
  typename v1::size_type overall_size;
  MarkDuplicateSortedKeyValuePairs(v1 keys_,v2 vals_, v3 prefix_sum_, typename v1::size_type overall_size_):
    keys(keys_), vals(vals_), prefix_sum(prefix_sum_), overall_size(overall_size_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i, typename v3::value_type &num_result) const {
    typename v1::value_type my_key = keys(i);
    typename v2::value_type my_val = vals(i);

    if ((my_key != 0 && my_val != 0) && ((i + 1 >= overall_size) || (my_key != keys(i + 1) || my_val != vals(i + 1)))){
      prefix_sum(i) = 1;
      num_result += 1;
    }
  }


};

template <typename v1, typename v2, typename v3, typename v4, typename v5>
struct FillSymmetricCSR{
  v1 keys;
  v2 vals;
  v3 prefix_sum;
  typename v3::size_type array_size;
  v4 out_xadj;
  v5 out_adj;
  FillSymmetricCSR(
      v1 keys_,v2 vals_, v3 prefix_sum_, typename v3::size_type array_size_,
      v4 out_xadj_, v5 out_adj_):
        keys(keys_), vals(vals_), prefix_sum(prefix_sum_), array_size(array_size_),
        out_xadj(out_xadj_), out_adj(out_adj_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int &i) const {
    typename v3::value_type my_pos = prefix_sum(i);

    if (i + 1 >= array_size){
      typename v2::value_type my_val = vals(i);
      typename v1::value_type my_key = keys(i);
      out_adj(my_pos) = my_val - 1;
      out_xadj(my_key) = my_pos + 1;
    }
    else {
      typename v3::value_type next_pos = prefix_sum(i + 1);
      if (my_pos != next_pos){

        typename v2::value_type my_val = vals(i);
        typename v1::value_type my_key = keys(i);
        typename v1::value_type next_key = keys(i + 1);
        out_adj(my_pos) = my_val - 1;
        if (my_key != next_key){
          out_xadj(my_key) = my_pos + 1;

        }

      }
    }
  }


};

/*
template <typename idx_array_type, typename idx_edge_array_type, typename hashmap, typename team_member>
struct InsertEdgesToMap{
  idx_array_type xadj;
  idx_array_type adj;
  hashmap my_map;
  InsertEdgesToMap(idx_array_type xadj_, idx_array_type adj_, hashmap my_map_):
    xadj(xadj_), adj(adj_), my_map(my_map_) {}

  KOKKOS_INLINE_FUNCTION
    void operator()(const team_member & teamMember, bool &failure) const {
      idx ii = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (ii >= num_rows) return;
      idx row_begin = xadj[ii];
      idx row_end = xadj[ii + 1];

      bool success = true;
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, row_end - row_begin),
          [&] (idx i) {
        idx adjind = i + row_begin;
        idx colIndex = adj[adjind];
        if (colIndex != ii){
          success = success && my_map.insert(ii, colIndex).success();
        }
      });

      failure = failure || !success;
    }


  KOKKOS_INLINE_FUNCTION
  void init (bool & dst) const {
    dst = false;
  }

  KOKKOS_INLINE_FUNCTION
  void join (volatile value_type& dst, const volatile value_type& src) const {
    dst = dst || src;
  }
};

template <typename idx_array_type, typename idx_edge_array_type, typename hashmap, typename team_member>
struct InsertEdgesToMap{
  idx_array_type xadj;
  idx_array_type adj;
  hashmap my_map;
  InsertEdgesToMap(idx_array_type xadj_, idx_array_type adj_, hashmap my_map_):
    xadj(xadj_), adj(adj_), my_map(my_map_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member & teamMember, bool &failure) const {
    idx ii = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
    if (ii >= num_rows) return;
    idx row_begin = xadj[ii];
    idx row_end = xadj[ii + 1];

    bool success = true;
    Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(teamMember, row_end - row_begin),
        [&] (idx i) {
      idx adjind = i + row_begin;
      idx colIndex = adj[adjind];
      if (colIndex != ii){
        success = success && my_map.insert(ii, colIndex).success();
      }
    });

    failure = failure || !success;
  }


  KOKKOS_INLINE_FUNCTION
  void init (bool & dst) const {
    dst = false;
  }

  KOKKOS_INLINE_FUNCTION
  void join (volatile value_type& dst, const volatile value_type& src) const {
    dst = dst || src;
  }
};


template <typename idx_array_type, typename idx_edge_array_type, typename hashmap, typename team_member>
struct CheckReverseEdgesInMap{
  idx_array_type xadj;
  idx_array_type adj;
  hashmap my_map;
  InsertEdgesToMap(idx_array_type xadj_, idx_array_type adj_, hashmap my_map_):
    xadj(xadj_), adj(adj_), my_map(my_map_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member & teamMember, bool &failure) const {
    idx ii = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
    if (ii >= num_rows) return;
    idx row_begin = xadj[ii];
    idx row_end = xadj[ii + 1];

    bool success = true;
    Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(teamMember, row_end - row_begin),
        [&] (idx i) {
      idx adjind = i + row_begin;
      idx colIndex = adj[adjind];
      if (colIndex != ii){
        success = success && my_map.insert(ii, colIndex).success();
        int index = cmap.find(colIndex, ii);
      }
    });

    failure = failure || !success;
  }


  KOKKOS_INLINE_FUNCTION
  void init (bool & dst) const {
    dst = false;
  }

  KOKKOS_INLINE_FUNCTION
  void join (volatile value_type& dst, const volatile value_type& src) const {
    dst = dst || src;
  }
};

template <typename idx_array_type, typename idx_edge_array_type, typename MyExecSpace>
bool is_graph_symmetric_symbolic(
    typename idx_array_type::value_type num_rows,
    typename idx_array_type::value_type num_cols,
    typename idx_array_type::value_type nnz,
    idx_array_type xadj,
    idx_edge_array_type adj,
    ){

  if (num_rows != num_cols){
    return false;
  }

  typedef typename idx_array_type::value_type idx;
  typedef Kokkos::UnorderedMap< idx, idx, ExecSpace > edge_map;
  edge_map EdgeMap(nnz);



  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy ;
  typedef typename team_policy::member_type team_member ;

  InsertEdgesToMap<idx_array_type, idx_edge_array_type,edge_map, team_member> ietm(xadj, adj, EdgeMap);

  int teamSizeMax = 0;
  int vector_size = 0;
  int max_allowed_team_size = team_policy::team_size_max(ietm);

  get_suggested_vector_team_size<idx, MyExecSpace>(
      max_allowed_team_size,
      vector_size,
      teamSizeMax,
      num_rows, nnz);


  bool insert_failure = false;
  Kokkos::parallel_for(
            team_policy(num_rows / teamSizeMax + 1 , teamSizeMax, vector_size), ietm, insert_failure);
  MyExecSpace::fence();
  if (insert_failure){
    return false;
  }



}
*/
template <typename idx_array_type, typename idx_edge_array_type,
          typename idx_out_array_type, typename idx_out_edge_array_type, typename MyExecSpace>
void symmetrize_graph_symbolic(
    typename idx_array_type::value_type num_rows_to_symmetrize,
    idx_array_type xadj,
    idx_edge_array_type adj,
    idx_out_array_type &sym_xadj,
    idx_out_edge_array_type &sym_adj
    ){

  typedef typename idx_array_type::value_type idx;


  idx nnz = adj.dimension_0();

  idx_out_edge_array_type tmp_srcs("tmpsrc", nnz * 2);
  idx_out_edge_array_type tmp_dsts("tmpdst",nnz * 2);

  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy ;
  typedef typename team_policy::member_type team_member ;

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;

  FillSymmetricEdges <idx_array_type,idx_edge_array_type,idx_out_edge_array_type, team_member> fse(
      num_rows_to_symmetrize,
      xadj,
      adj,
      tmp_srcs,
      tmp_dsts
      );



  int teamSizeMax = 0;
  int vector_size = 0;
  int max_allowed_team_size = team_policy::team_size_max(fse);

  get_suggested_vector_team_size<idx, MyExecSpace>(
      max_allowed_team_size,
      vector_size,
      teamSizeMax,
      xadj.dimension_0() - 1, nnz);
  //std::cout << "max_allowed_team_size:" << max_allowed_team_size << " vs:" << vector_size << " tsm:" << teamSizeMax<< std::endl;

  Kokkos::parallel_for(
            team_policy(num_rows_to_symmetrize / teamSizeMax + 1 , teamSizeMax, vector_size),
            fse);
  MyExecSpace::fence();

  KokkosKernelsSorting::sort_key_value_views <idx_out_edge_array_type, idx_out_edge_array_type, MyExecSpace>(tmp_srcs, tmp_dsts);

  MyExecSpace::fence();

  idx_out_edge_array_type pps("PPS", nnz * 2);

  typename idx_out_edge_array_type::non_const_value_type num_symmetric_edges = 0;
  if (nnz > 0)
  Kokkos::parallel_reduce(
            my_exec_space(0, nnz * 2),
            MarkDuplicateSortedKeyValuePairs<idx_out_edge_array_type, idx_out_edge_array_type, idx_out_edge_array_type>(
                tmp_srcs, tmp_dsts, pps, nnz * 2), num_symmetric_edges);

  Kokkos::fence();
  if (nnz > 0)
  exclusive_parallel_prefix_sum<idx_out_edge_array_type, MyExecSpace>(nnz * 2, pps);

  MyExecSpace::fence();
  sym_xadj = idx_out_array_type("sym_xadj", num_rows_to_symmetrize + 1);
  sym_adj = idx_out_edge_array_type("sym_adj", num_symmetric_edges);

  MyExecSpace::fence();
  Kokkos::parallel_for(
        my_exec_space(0, nnz * 2),
        FillSymmetricCSR<idx_out_edge_array_type, idx_out_edge_array_type, idx_out_edge_array_type, idx_out_array_type, idx_out_edge_array_type>
  (tmp_srcs, tmp_dsts, pps, nnz * 2, sym_xadj, sym_adj));

  MyExecSpace::fence();
  remove_zeros_in_xadj_vector<idx_out_array_type, MyExecSpace>(num_rows_to_symmetrize + 1, sym_xadj);
  MyExecSpace::fence();
}

template <typename from_vector, typename to_vector>
struct CopyVector{
  from_vector from;
  to_vector to;

  CopyVector(from_vector &from_, to_vector to_): from(from_), to(to_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i) const {
    to[i] = from[i];
  }
};
template <typename from_vector, typename to_vector, typename MyExecSpace>
void copy_vector(
                size_t num_elements,
                from_vector from, to_vector to){

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  Kokkos::parallel_for( my_exec_space(0,num_elements), CopyVector<from_vector, to_vector>(from, to));

}

}
}
}
#endif
