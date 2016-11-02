#ifndef _KOKKOSKERNELS_SPARSEUTILS_HPP
#define _KOKKOSKERNELS_SPARSEUTILS_HPP
#include "Kokkos_Core.hpp"
#include "Kokkos_Atomic.hpp"
#include "impl/Kokkos_Timer.hpp"
#include "KokkosKernels_SimpleUtils.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"
namespace KokkosKernels{

namespace Experimental{

namespace Util{
template <typename in_row_view_t,
          typename in_nnz_view_t,
          typename in_scalar_view_t,
          typename out_row_view_t,
          typename out_nnz_view_t,
          typename out_scalar_view_t,
          typename tempwork_row_view_t,
          typename MyExecSpace>
struct TransposeMatrix{

  struct CountTag{};
  struct FillTag{};

  typedef struct CountTag CountTag;
  typedef struct FillTag FillTag;

  typedef Kokkos::TeamPolicy<CountTag, MyExecSpace> team_count_policy_t ;
  typedef Kokkos::TeamPolicy<FillTag, MyExecSpace> team_fill_policy_t ;

  typedef Kokkos::TeamPolicy<CountTag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> > dynamic_team_count_policy_t ;
  typedef Kokkos::TeamPolicy<FillTag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> > dynamic_team_fill_policy_t ;


  typedef typename team_count_policy_t::member_type team_count_member_t ;
  typedef typename team_fill_policy_t::member_type team_fill_member_t ;

  typedef typename in_nnz_view_t::non_const_value_type nnz_lno_t;
  typedef typename in_row_view_t::non_const_value_type size_type;


  typename in_nnz_view_t::non_const_value_type num_rows;
  typename in_nnz_view_t::non_const_value_type num_cols;
  in_row_view_t xadj;
  in_nnz_view_t adj;
  in_scalar_view_t vals;
  out_row_view_t t_xadj; //allocated
  out_nnz_view_t t_adj;  //allocated
  out_nnz_view_t t_vals;  //allocated
  tempwork_row_view_t tmp_txadj;
  bool transpose_values;
  nnz_lno_t team_work_size;

  TransposeMatrix(
      nnz_lno_t num_rows_,
      nnz_lno_t num_cols_,
      in_row_view_t xadj_,
      in_nnz_view_t adj_,
      in_scalar_view_t vals_,
      out_row_view_t t_xadj_,
      out_nnz_view_t t_adj_,
      out_nnz_view_t t_vals_,
      tempwork_row_view_t tmp_txadj_,
      bool transpose_values_,
      nnz_lno_t team_row_work_size_):
        num_rows(num_rows_), num_cols(num_cols_),
        xadj(xadj_), adj(adj_), vals(vals_),
        t_xadj(t_xadj_),  t_adj(t_adj_), t_vals(t_vals_),
        tmp_txadj(tmp_txadj_), transpose_values(transpose_values_), team_work_size(team_row_work_size_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const CountTag&, const team_count_member_t & teamMember) const {

    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, num_rows);
    //TODO we dont need to go over rows
    //just go over nonzeroes.
    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,team_row_begin,team_row_end), [&] (const nnz_lno_t& row_index) {
      const size_type col_begin = xadj[row_index];
      const size_type col_end = xadj[row_index + 1];
      const nnz_lno_t left_work = col_end - col_begin;
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, left_work),
          [&] (nnz_lno_t i) {
        const size_type adjind = i + col_begin;
        const nnz_lno_t colIndex = adj[adjind];
        Kokkos::atomic_fetch_add(&(t_xadj(colIndex)),1);
      });
    });
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const FillTag&, const team_fill_member_t & teamMember) const {
    const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, num_rows);


    Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember,team_row_begin,team_row_end), [&] (const nnz_lno_t& row_index) {
    const nnz_lno_t teamsize = teamMember.team_size();
    //for (nnz_lno_t row_index = team_row_begin + teamMember.team_rank(); row_index < team_row_end; row_index += teamsize){
      const size_type col_begin = xadj[row_index];
      const size_type col_end = xadj[row_index + 1];
      const nnz_lno_t left_work = col_end - col_begin;
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, left_work),
          [&] (nnz_lno_t i) {
        const size_type adjind = i + col_begin;
        const nnz_lno_t colIndex = adj[adjind];
        const size_type pos = Kokkos::atomic_fetch_add(&(tmp_txadj(colIndex)),1);

        t_adj(pos) = row_index;
        if (transpose_values){
          t_vals(pos) = vals[adjind];
        }

      });
    //}
    });
  }
};

/**
 * \brief function returns transpose of the given graph.
 * \param num_rows: num rows in input graph
 * \param num_cols: num cols in input graph
 * \param xadj: row pointers of the input graph
 * \param adj: column indices of the input graph
 * \param t_xadj: output, the row indices of the output graph. MUST BE INITIALIZED WITH ZEROES.
 * \param t_adj: output, column indices. No need for initializations.
 * \param vector_size: suggested vector size, optional. if -1, kernel will decide.
 * \param suggested_team_size: suggested team size, optional. if -1, kernel will decide.
 * \param team_work_chunk_size: suggested work size of a team, optional. if -1, kernel will decide.
 * \param use_dynamic_scheduling: whether to use dynamic scheduling. Default is true.
 */
template <typename in_row_view_t,
          typename in_nnz_view_t,
          typename out_row_view_t,
          typename out_nnz_view_t,
          typename tempwork_row_view_t,
          typename MyExecSpace>
inline void kk_transpose_graph(
    typename in_nnz_view_t::non_const_value_type num_rows,
    typename in_nnz_view_t::non_const_value_type num_cols,
    in_row_view_t xadj,
    in_nnz_view_t adj,
    out_row_view_t t_xadj, //pre-allocated -- initialized with 0
    out_nnz_view_t t_adj,  //pre-allocated -- no need for initialize
    int vector_size = -1,
    int suggested_team_size = -1,
    typename in_nnz_view_t::non_const_value_type team_work_chunk_size = -1,
    bool use_dynamic_scheduling = true
    ){

  //allocate some memory for work for row pointers
  tempwork_row_view_t tmp_row_view(Kokkos::ViewAllocateWithoutInitializing("tmp_row_view"), num_cols + 1);

  in_nnz_view_t tmp1;
  out_nnz_view_t tmp2;

  //create the functor for tranpose.
  typedef TransposeMatrix <
      in_row_view_t, in_nnz_view_t, in_nnz_view_t,
      out_row_view_t, out_nnz_view_t, out_nnz_view_t,
      tempwork_row_view_t, MyExecSpace>  TransposeFunctor_t;

  TransposeFunctor_t tm ( num_rows, num_cols, xadj, adj, tmp1,
                          t_xadj, t_adj, tmp2,
                          tmp_row_view,
                          false,
                          team_work_chunk_size);

  typedef typename TransposeFunctor_t::team_count_policy_t count_tp_t;
  typedef typename TransposeFunctor_t::team_fill_policy_t fill_tp_t;
  typedef typename TransposeFunctor_t::dynamic_team_count_policy_t d_count_tp_t;
  typedef typename TransposeFunctor_t::dynamic_team_fill_policy_t d_fill_tp_t;

  typename in_row_view_t::non_const_value_type nnz = adj.dimension_0();

  //set the vector size, if not suggested.
  if (vector_size == -1)
    vector_size = kk_get_suggested_vector_size(num_rows, nnz, kk_get_exec_space_type<MyExecSpace>());

  //set the team size, if not suggested.
  if (suggested_team_size == -1)
    suggested_team_size = kk_get_suggested_team_size(vector_size, kk_get_exec_space_type<MyExecSpace>());

  //set the chunk size, if not suggested.
  if (team_work_chunk_size == -1)
    team_work_chunk_size = suggested_team_size;



  if (use_dynamic_scheduling){
    Kokkos::parallel_for(  d_count_tp_t(num_rows  / team_work_chunk_size + 1 , suggested_team_size, vector_size), tm);
  }
  else {
    Kokkos::parallel_for(  count_tp_t(num_rows  / team_work_chunk_size + 1 , suggested_team_size, vector_size), tm);
  }
  MyExecSpace::fence();

  kk_exclusive_parallel_prefix_sum<out_row_view_t, MyExecSpace>(num_cols+1, t_xadj);
  MyExecSpace::fence();

  Kokkos::deep_copy(tmp_row_view, t_xadj);
  MyExecSpace::fence();


  if (use_dynamic_scheduling){
    Kokkos::parallel_for(  fill_tp_t(num_rows  / team_work_chunk_size + 1 , suggested_team_size, vector_size), tm);
  }
  else {
    Kokkos::parallel_for(  d_fill_tp_t(num_rows  / team_work_chunk_size + 1 , suggested_team_size, vector_size), tm);
  }
  MyExecSpace::fence();
}

template <typename forward_map_type, typename reverse_map_type>
struct Fill_Reverse_Scale_Functor{

  struct CountTag{};
  struct FillTag{};

  typedef struct CountTag CountTag;
  typedef struct FillTag FillTag;


  typedef typename forward_map_type::value_type forward_type;
  typedef typename reverse_map_type::value_type reverse_type;
  forward_map_type forward_map;
  reverse_map_type reverse_map_xadj;
  reverse_map_type reverse_map_adj;

  const reverse_type multiply_shift_for_scale;
  const reverse_type division_shift_for_bucket;


  Fill_Reverse_Scale_Functor(
      forward_map_type forward_map_,
      reverse_map_type reverse_map_xadj_,
      reverse_map_type reverse_map_adj_,
      reverse_type multiply_shift_for_scale_,
      reverse_type division_shift_for_bucket_):
        forward_map(forward_map_), reverse_map_xadj(reverse_map_xadj_), reverse_map_adj(reverse_map_adj_),
        multiply_shift_for_scale(multiply_shift_for_scale_),
        division_shift_for_bucket(division_shift_for_bucket_){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const CountTag&, const size_t &ii) const {
    forward_type fm = forward_map[ii];
    fm = fm << multiply_shift_for_scale;
    fm += ii >> division_shift_for_bucket;
    Kokkos::atomic_fetch_add( &(reverse_map_xadj(fm)), 1);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const FillTag&, const size_t &ii) const {
    forward_type fm = forward_map[ii];

    fm = fm << multiply_shift_for_scale;
    fm += ii >> division_shift_for_bucket;
    const reverse_type future_index = Kokkos::atomic_fetch_add( &(reverse_map_xadj(fm )), 1);
    reverse_map_adj(future_index) = ii;
  }
};


template <typename from_view_t, typename to_view_t>
struct StridedCopy1{
  const from_view_t from;
  to_view_t to;
  const size_t stride;
  StridedCopy1(
      const from_view_t from_,
      to_view_t to_,
      size_t stride_):from(from_), to (to_), stride(stride_){}


  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &ii) const {
    to[ii] = from[(ii) * stride];
  }
};

template <typename forward_map_type, typename reverse_map_type>
struct Reverse_Map_Functor{

  struct CountTag{};
  struct FillTag{};

  typedef struct CountTag CountTag;
  typedef struct FillTag FillTag;


  typedef typename forward_map_type::value_type forward_type;
  typedef typename reverse_map_type::value_type reverse_type;
  forward_map_type forward_map;
  reverse_map_type reverse_map_xadj;
  reverse_map_type reverse_map_adj;


  Reverse_Map_Functor(
      forward_map_type forward_map_,
      reverse_map_type reverse_map_xadj_,
      reverse_map_type reverse_map_adj):
        forward_map(forward_map_), reverse_map_xadj(reverse_map_xadj_), reverse_map_adj(reverse_map_adj){}

  KOKKOS_INLINE_FUNCTION
  void operator()(const CountTag&, const size_t &ii) const {
    forward_type fm = forward_map[ii];
    Kokkos::atomic_fetch_add( &(reverse_map_xadj(fm)), 1);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const FillTag&, const size_t &ii) const {
    forward_type c = forward_map[ii];
    const reverse_type future_index = Kokkos::atomic_fetch_add( &(reverse_map_xadj(c)), 1);
    reverse_map_adj(future_index) = ii;
  }
};


/**
 * \brief Utility function to obtain a reverse map given a map.
 * Input is a map with the number of elements within the map.
 * forward_map[c] = i, where c is a forward element and forward_map has a size of num_forward_elements.
 * i is the value that c is mapped in the forward map, and the range of that is num_reverse_elements.
 * Output is the reverse_map_xadj and reverse_map_adj such that,
 * all c, forward_map[c] = i, will appear in  reverse_map_adj[ reverse_map_xadj[i]: reverse_map_xadj[i+1])
 * \param: num_forward_elements: the number of elements in the forward map, the size of the forward map.
 * \param: num_reverse_elements: the number of elements that forward map is mapped to. It is the value of max i.
 * \param: forward_map: input forward_map, where forward_map[c] = i.
 * \param: reverse_map_xadj: reverse map xadj, that is it will hold the beginning and
 * end indices on reverse_map_adj such that all values mapped to i will be [ reverse_map_xadj[i]: reverse_map_xadj[i+1])
 * its size will be num_reverse_elements + 1. NO NEED TO INITIALIZE.
 * \param: reverse_map_adj: reverse map adj, holds the values of reverse maps. Its size is num_forward_elements.
 *
 */
template <typename forward_array_type, typename reverse_array_type, typename MyExecSpace>
void kk_create_reverse_map(
    const typename reverse_array_type::value_type &num_forward_elements, //num_vertices
    const typename forward_array_type::value_type &num_reverse_elements, //num_colors

    const forward_array_type &forward_map, //vertex to colors
    const reverse_array_type &reverse_map_xadj, // colors to vertex xadj
    const reverse_array_type &reverse_map_adj){ //colros to vertex adj

  typedef typename reverse_array_type::value_type lno_t;
  typedef typename forward_array_type::value_type reverse_lno_t;

  const lno_t  MINIMUM_TO_ATOMIC = 128;

  //typedef Kokkos::TeamPolicy<CountTag, MyExecSpace> team_count_policy_t ;
  //typedef Kokkos::TeamPolicy<FillTag, MyExecSpace> team_fill_policy_t ;

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;

  //IF There are very few reverse elements, atomics are likely to create contention.
  if (num_reverse_elements < MINIMUM_TO_ATOMIC){
    const lno_t scale_size = 1024;
    const lno_t multiply_shift_for_scale = 10;

    //there will be 1024 buckets
    const lno_t division_shift_for_bucket =
          lno_t (ceil(log(double (num_forward_elements) / scale_size)/log(2)));

    //coloring indices are base-1. we end up using not using element 1.
    const reverse_lno_t tmp_reverse_size =
        (num_reverse_elements + 1) << multiply_shift_for_scale;

    typename reverse_array_type::non_const_type
        tmp_color_xadj ("TMP_REVERSE_XADJ", tmp_reverse_size + 1);

    typedef Fill_Reverse_Scale_Functor<forward_array_type, reverse_array_type> frsf;
    typedef typename frsf::CountTag cnt_tag;
    typedef typename frsf::FillTag fill_tag;
    typedef Kokkos::RangePolicy<cnt_tag, MyExecSpace> my_cnt_exec_space;
    typedef Kokkos::RangePolicy<fill_tag, MyExecSpace> my_fill_exec_space;

    frsf frm (forward_map, tmp_color_xadj, reverse_map_adj,
            multiply_shift_for_scale, division_shift_for_bucket);

    Kokkos::parallel_for (my_cnt_exec_space (0, num_forward_elements) , frm);
    MyExecSpace::fence();


    //kk_inclusive_parallel_prefix_sum<reverse_array_type, MyExecSpace>(tmp_reverse_size + 1, tmp_color_xadj);
    kk_exclusive_parallel_prefix_sum<reverse_array_type, MyExecSpace>
      (tmp_reverse_size + 1, tmp_color_xadj);
    MyExecSpace::fence();

    Kokkos::parallel_for (
        my_exec_space (0, num_reverse_elements + 1) ,
        StridedCopy1<reverse_array_type, reverse_array_type>
          (tmp_color_xadj, reverse_map_xadj, scale_size));
    MyExecSpace::fence();
    Kokkos::parallel_for (my_fill_exec_space (0, num_forward_elements) , frm);
    MyExecSpace::fence();
  }
  else
  //atomic implementation.
  {
    reverse_array_type tmp_color_xadj ("TMP_REVERSE_XADJ", num_reverse_elements + 1);

    typedef Reverse_Map_Functor<forward_array_type, reverse_array_type> rmp_functor_type;
    typedef typename rmp_functor_type::CountTag cnt_tag;
    typedef typename rmp_functor_type::FillTag fill_tag;
    typedef Kokkos::RangePolicy<cnt_tag, MyExecSpace> my_cnt_exec_space;
    typedef Kokkos::RangePolicy<fill_tag, MyExecSpace> my_fill_exec_space;

    rmp_functor_type frm (forward_map, tmp_color_xadj, reverse_map_adj);

    Kokkos::parallel_for (my_cnt_exec_space (0, num_forward_elements) , frm);
    MyExecSpace::fence();

    //kk_inclusive_parallel_prefix_sum<reverse_array_type, MyExecSpace>(num_reverse_elements + 1, reverse_map_xadj);
    kk_exclusive_parallel_prefix_sum<reverse_array_type, MyExecSpace>
      (num_reverse_elements + 1, tmp_color_xadj);
    MyExecSpace::fence();

    Kokkos::deep_copy (reverse_map_xadj, tmp_color_xadj);
    MyExecSpace::fence();

    Kokkos::parallel_for (my_fill_exec_space (0, num_forward_elements) , frm);
    MyExecSpace::fence();
  }
}

}
}
}
#endif
