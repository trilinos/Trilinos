
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_Core.hpp>
#include "KokkosKernels_Utils.hpp"
#ifndef _GRAPHCOLORHANDLE_HPP
#define _GRAPHCOLORHANDLE_HPP
//#define VERBOSE
namespace KokkosKernels{

namespace Experimental{


namespace Graph{

enum ColoringAlgorithm{COLORING_DEFAULT, COLORING_SERIAL, COLORING_VB, COLORING_VBBIT, COLORING_VBCS, COLORING_EB};

enum ConflictList{COLORING_NOCONFLICT, COLORING_ATOMIC, COLORING_PPS};

enum ColoringType {Distance1, Distance2};

template <class lno_row_view_t_, class nonconst_color_view_t_, class lno_nnz_view_t_,
      class ExecutionSpace, class TemporaryMemorySpace, class PersistentMemorySpace>
class GraphColoringHandle{

public:
  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;


  typedef lno_row_view_t_ in_lno_row_view_t;
  //typedef row_index_view_type_ idx_array_type;

  //typedef nonzero_index_view_type_ idx_edge_array_type;
  typedef lno_nnz_view_t_ in_lno_nnz_view_t;

  //typedef typename row_index_view_type::value_type idx;
  typedef typename in_lno_row_view_t::non_const_value_type row_lno_t;

  //typedef typename row_index_view_type::array_layout idx_array_layout;
  typedef typename in_lno_row_view_t::array_layout row_lno_view_array_layout;

  //typedef typename row_index_view_type::device_type idx_device_type;
  typedef typename in_lno_row_view_t::device_type row_lno_view_device_t;

  //typedef typename row_index_view_type::memory_traits idx_memory_traits;
  typedef typename in_lno_row_view_t::memory_traits row_lno_view_memory_traits;

  //typedef typename row_index_view_type::HostMirror host_view_type;
  typedef typename in_lno_row_view_t::HostMirror row_lno_host_view_t; //Host view type
  //typedef typename idx_memory_traits::MemorySpace MyMemorySpace;

  //typedef typename nonzero_index_view_type::non_const_value_type idx_edge;
  typedef typename in_lno_nnz_view_t::non_const_value_type nnz_lno_t;

  //typedef typename nonzero_index_view_type::array_layout idx_edge_array_layout;
  typedef typename in_lno_nnz_view_t::array_layout nnz_lno_view_array_layout;

  //typedef typename nonzero_index_view_type::device_type idx_edge_device_type;
  typedef typename in_lno_nnz_view_t::device_type nnz_lno_view_device_t;

  //typedef typename nonzero_index_view_type::memory_traits idx_edge_memory_traits;
  typedef typename in_lno_nnz_view_t::memory_traits nnz_lno_view_memory_traits;

  //typedef typename nonzero_index_view_type::HostMirror host_edge_view_type; //Host view type
  typedef typename in_lno_nnz_view_t::HostMirror nnz_lno_host_view_t; //Host view type



  typedef nonconst_color_view_t_ color_view_t;

  typedef typename color_view_t::non_const_value_type color_t;
  typedef typename color_view_t::array_layout color_view_array_layout;
  typedef typename color_view_t::device_type color_view_device_t;
  typedef typename color_view_t::memory_traits color_view_memory_traits;
  typedef typename color_view_t::HostMirror color_host_view_t; //Host view type



  typedef typename in_lno_row_view_t::const_data_type const_row_lno_t;
  typedef typename in_lno_row_view_t::non_const_data_type non_const_row_lno_t;
  //typedef typename in_row_index_view_type::memory_space row_view_memory_space;
  //typedef typename Kokkos::View<const_row_lno_t, row_view_array_layout,
  //    row_view_device_t, row_view_memory_traits> const_row_index_view_type;

  typedef typename in_lno_row_view_t::const_type const_lno_row_view_t;

  //typedef typename Kokkos::View<non_const_row_lno_t, row_view_array_layout,
  //    row_view_device_t, row_view_memory_traits> non_const_row_index_view_type;
  typedef typename in_lno_row_view_t::non_const_type non_const_lno_row_view_t;



  typedef typename in_lno_nnz_view_t::const_data_type const_nnz_lno_t;
  typedef typename in_lno_nnz_view_t::non_const_data_type non_const_nnz_lno_t;
  //typedef typename in_nonzero_index_view_type::memory_space nonzero_index_view_memory_space;
  //typedef typename Kokkos::View<const_nnz_lno_t, nnz_lno_view_array_layout,
  //    nnz_lno_view_device_t, nnz_lno_view_memory_traits> const_nonzero_index_view_type;

  typedef typename in_lno_nnz_view_t::const_type const_lno_nnz_view_t;

  //typedef typename Kokkos::View<non_const_nnz_lno_t, nnz_lno_view_array_layout,
  //    nnz_lno_view_device_t, nnz_lno_view_memory_traits> non_const_nonzero_index_view_type;
  typedef typename in_lno_nnz_view_t::non_const_type non_const_lno_nnz_view_t;




  //typedef typename Kokkos::View<row_index_type *, HandleTempMemorySpace> idx_temp_work_array_type;
  typedef typename Kokkos::View<row_lno_t *, HandleTempMemorySpace> row_lno_temp_work_view_t;
  //typedef typename row_index_persistent_work_view_type idx_persistent_work_array_type;
  typedef typename Kokkos::View<row_lno_t *, PersistentMemorySpace> row_lno_persistent_work_view_t;

  //typedef typename row_index_persistent_work_view_type::HostMirror host_idx_persistent_view_type; //Host view type
  typedef typename row_lno_persistent_work_view_t::HostMirror row_lno_persistent_work_host_view_t; //Host view type

  typedef typename Kokkos::View<nnz_lno_t *, HandleTempMemorySpace> nnz_lno_temp_work_view_t;
  typedef typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace> nnz_lno_persistent_work_view_t;
  typedef typename nnz_lno_persistent_work_view_t::HostMirror nnz_lno_persistent_work_host_view_t; //Host view type
  typedef Kokkos::TeamPolicy<HandleExecSpace> team_policy_t ;

  typedef typename team_policy_t::member_type team_member_t ;
private:

  ColoringType GraphColoringType;
  //Parameters
  ColoringAlgorithm coloring_algorithm_type; //VB, VBBIT or EB.
  ConflictList conflict_list_type;  // whether to use a conflict list or not, and
                                    // if using it wheter to create it with atomic or parallel prefix sum.

  double min_reduction_for_conflictlist;
                      //if used pps is selected to create conflict list, what min percantage should be the vertex list
                      //be reduced, to create the new vertexlist. If it is reduced less than this percantage, use the
                      //previous array.

  int min_elements_for_conflictlist;
                            //minimum number of elements to create a new conflict list.
                            //if current conflict list size is smaller than this number,
                            //than we do not need to create a new conflict list.

  bool serial_conflict_resolution;//perform parallel greedy coloring once, then resolve conflict serially.
  bool tictoc; //print time at every step

  bool vb_edge_filtering;  //whether to do edge filtering or not in vertex based algorithms. Swaps on the ad error.

  int vb_chunk_size;  //the (minimum) size of the consecutive works that a thread will be assigned to.
  int max_number_of_iterations; //maximum allowed number of phases

  int eb_num_initial_colors; //the number of colors to assign at the beginning of the edge-based algorithm

  //STATISTICS
  double overall_coloring_time; //the overall time that it took to color the graph. In the case of the iterative calls.
  double coloring_time; //the time that it took to color the graph
  int num_phases; //


  row_lno_t size_of_edge_list;
  row_lno_persistent_work_view_t lower_triangle_src;
  row_lno_persistent_work_view_t lower_triangle_dst;

  color_view_t vertex_colors;
  bool is_coloring_called_before;
  row_lno_t num_colors;



  public:


  /**
   * \brief Default constructor.
   */
  GraphColoringHandle():
    GraphColoringType(Distance1),
    coloring_algorithm_type(COLORING_DEFAULT),
    conflict_list_type(COLORING_ATOMIC),
    min_reduction_for_conflictlist(0.35),
    min_elements_for_conflictlist(1000 /*5000*/),
    serial_conflict_resolution(false),
    tictoc(false),
    vb_edge_filtering(false),
    vb_chunk_size(8),
    max_number_of_iterations(200), eb_num_initial_colors(1),
    overall_coloring_time(0),
    coloring_time(0),
    num_phases(0), size_of_edge_list(0), lower_triangle_src(), lower_triangle_dst(),
    vertex_colors(), is_coloring_called_before(false), num_colors(0){
    this->choose_default_algorithm();
    this->set_defaults(this->coloring_algorithm_type);
  }

  /** \brief Sets the graph coloring type. Whether it is distance-1 or distance-2 coloring.
   *  \param col_type: Coloring Type: KokkosKernels::Experimental::Graph::ColoringType which can be
   *        either KokkosKernels::Experimental::Graph::Distance1 or KokkosKernels::Experimental::Graph::Distance2
   */
  void set_coloring_type(const ColoringType &col_type){
    this->GraphColoringType = col_type;
  }

  /** \brief Gets the graph coloring type. Whether it is distance-1 or distance-2 coloring.
   *  returns Coloring Type: KokkosKernels::Experimental::Graph::ColoringType which can be
   *        either KokkosKernels::Experimental::Graph::Distance1 or KokkosKernels::Experimental::Graph::Distance2
   */
  ColoringType get_coloring_type(){
    return this->GraphColoringType;
  }

  /** \brief Changes the graph coloring algorithm.
   *  \param col_algo: Coloring algorithm: one of COLORING_VB, COLORING_VBBIT, COLORING_VBCS, COLORING_EB
   *  \param set_default_parameters: whether or not to reset the default parameters for the given algorithm.
   */
  void set_algorithm(const ColoringAlgorithm &col_algo, bool set_default_parameters = true){
    if (col_algo == COLORING_DEFAULT){
      this->choose_default_algorithm();
    }
    else {
      this->coloring_algorithm_type = col_algo;
    }
    if (set_default_parameters){
      this->set_defaults(this->coloring_algorithm_type);
    }
  }

  /** \brief Chooses best algorithm based on the execution space. COLORING_EB if cuda, COLORING_VB otherwise.
   */
  void choose_default_algorithm(){

#if defined( KOKKOS_HAVE_SERIAL )
    if (Kokkos::Impl::is_same< Kokkos::Serial , ExecutionSpace >::value){
      this->coloring_algorithm_type = COLORING_SERIAL;
#ifdef VERBOSE
      std::cout << "Serial Execution Space, Default Algorithm: COLORING_VB" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_PTHREAD )
    if (Kokkos::Impl::is_same< Kokkos::Threads , ExecutionSpace >::value){
      this->coloring_algorithm_type = COLORING_VB;
#ifdef VERBOSE
      std::cout << "PTHREAD Execution Space, Default Algorithm: COLORING_VB" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_OPENMP )
    if (Kokkos::Impl::is_same< Kokkos::OpenMP, ExecutionSpace >::value){
      this->coloring_algorithm_type = COLORING_VB;
#ifdef VERBOSE
      std::cout << "OpenMP Execution Space, Default Algorithm: COLORING_VB" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_CUDA )
    if (Kokkos::Impl::is_same<Kokkos::Cuda, ExecutionSpace >::value){
      this->coloring_algorithm_type = COLORING_EB;
#ifdef VERBOSE
      std::cout << "Qthread Execution Space, Default Algorithm: COLORING_VB" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_QTHREAD)
    if (Kokkos::Impl::is_same< Kokkos::Qthread, ExecutionSpace >::value){
      this->coloring_algorithm_type = COLORING_VB;
#ifdef VERBOSE
      std::cout << "Qthread Execution Space, Default Algorithm: COLORING_VB" << std::endl;
#endif
    }
#endif
  }

  template<typename v1, typename v2, typename v3>
  struct CountLowerTriangle{
    row_lno_t nv;
    v1 xadj;
    v2 adj;
    v3 lower_xadj_counts;

    CountLowerTriangle(
        row_lno_t nv_,
        v1 xadj_,
        v2 adj_,
        v3 lower_xadj_counts_
        ): nv(nv_),
            xadj(xadj_), adj(adj_),
            lower_xadj_counts(lower_xadj_counts_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const row_lno_t &i, row_lno_t &new_num_edge) const {
      row_lno_t xadj_begin = xadj(i);
      row_lno_t xadj_end = xadj(i + 1);

      row_lno_t new_edge_count = 0;
      for (row_lno_t j = xadj_begin; j < xadj_end; ++j){
        row_lno_t n = adj(j);
        if (i < n && n < nv){
          new_edge_count += 1;
        }
      }
      lower_xadj_counts(i + 1) = new_edge_count;
      new_num_edge += new_edge_count;
    }
  };

  template<typename v1, typename v2, typename v3>
  struct CountLowerTriangleTeam{

    row_lno_t nv;
    v1 xadj;
    v2 adj;
    v3 lower_xadj_counts;

    CountLowerTriangleTeam(
        row_lno_t nv_,
        v1 xadj_,
        v2 adj_,
        v3 lower_xadj_counts_
        ): nv(nv_),
            xadj(xadj_), adj(adj_),
            lower_xadj_counts(lower_xadj_counts_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember/*, row_lno_t &new_num_edge*/) const {


      row_lno_t ii = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (ii >= nv) {
        return;
      }

      row_lno_t xadj_begin = xadj(ii);
      row_lno_t xadj_end = xadj(ii + 1);

      row_lno_t new_edge_count = 0;


      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(teamMember, xadj_end - xadj_begin),
          [&] (row_lno_t i, row_lno_t &numEdges) {

        row_lno_t adjind = i + xadj_begin;
        row_lno_t n = adj[adjind];
        if (ii < n && n < nv){
          numEdges += 1;
        }
      }, new_edge_count);

      Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
        lower_xadj_counts(ii + 1) = new_edge_count;
      });
    }

  };

  template<typename v1, typename v2, typename v3, typename v4>
  struct FillLowerTriangleTeam{
    row_lno_t nv;
    v1 xadj;
    v2 adj;
    v3 lower_xadj_counts;
    v4 lower_srcs;
    v4 lower_dsts;

    FillLowerTriangleTeam(
        row_lno_t nv_,
        v1 xadj_,
        v2 adj_,
        v3 lower_xadj_counts_,
        v4 lower_srcs_,
        v4 lower_dsts_
        ):  nv(nv_),
            xadj(xadj_), adj(adj_),
            lower_xadj_counts(lower_xadj_counts_),
            lower_srcs(lower_srcs_), lower_dsts(lower_dsts_) {}


    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember) const {


      row_lno_t ii = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (ii >= nv) {
        return;
      }

      row_lno_t xadj_begin = xadj(ii);
      row_lno_t xadj_end = xadj(ii + 1);

      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, xadj_end - xadj_begin),
          [&] (row_lno_t i) {

        row_lno_t adjind = i + xadj_begin;
        row_lno_t n = adj[adjind];
        if (ii < n && n < nv){
          row_lno_t position =
              Kokkos::atomic_fetch_add( &(lower_xadj_counts(ii)), 1);
          lower_srcs(position) = ii;
          lower_dsts(position) = n;
        }
      });
    }
  };


  template<typename v1, typename v2, typename v3, typename v4>
  struct FillLowerTriangle{
    row_lno_t nv;
    v1 xadj;
    v2 adj;
    v3 lower_xadj_counts;
    v4 lower_srcs;
    v4 lower_dsts;

    FillLowerTriangle(
        row_lno_t nv_,
        v1 xadj_,
        v2 adj_,
        v3 lower_xadj_counts_,
        v4 lower_srcs_,
        v4 lower_dsts_
        ):  nv(nv_),
            xadj(xadj_), adj(adj_),
            lower_xadj_counts(lower_xadj_counts_),
            lower_srcs(lower_srcs_), lower_dsts(lower_dsts_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const row_lno_t &i) const{

      row_lno_t xadj_begin = xadj[i];
      row_lno_t xadj_end = xadj[i + 1];

      for (row_lno_t j = xadj_begin; j < xadj_end; ++j){
        row_lno_t n = adj(j);
        if (i < n && n < nv){
          row_lno_t position = lower_xadj_counts(i)++;
          lower_srcs(position) = i;
          lower_dsts(position) = n;
        }
      }
    }
  };
  template <typename row_index_view_type, typename nonzero_view_type>
  void symmetrize_and_calculate_lower_diagonal_edge_list(
      row_lno_t nv,
      row_index_view_type xadj, nonzero_view_type adj){

    KokkosKernels::Experimental::Util::symmetrize_and_get_lower_diagonal_edge_list
    <row_index_view_type, nonzero_view_type, row_lno_persistent_work_view_t, HandleExecSpace>
      (
        nv,
        xadj,
        adj,
        lower_triangle_src,
        lower_triangle_dst);

    size_of_edge_list = lower_triangle_src.dimension_0();

  }

  template <typename row_index_view_type, typename nonzero_view_type>
  void get_lower_diagonal_edge_list(
      row_lno_t nv, row_lno_t ne,
      row_index_view_type xadj, nonzero_view_type adj,
      row_lno_t  &num_out_edges,
      row_lno_persistent_work_view_t &src,
      row_lno_persistent_work_view_t &dst){

    if (size_of_edge_list > 0){
      num_out_edges = size_of_edge_list;
      //src = Kokkos::View<idx *, HandlePersistentMemorySpace> (this->lower_triangle_src);
      //dst = Kokkos::View<idx *, HandlePersistentMemorySpace> (this->lower_triangle_dst);
      src = (this->lower_triangle_src);
      dst = (this->lower_triangle_dst);
    }
    else {

      row_lno_temp_work_view_t lower_count("LowerXADJ", nv + 1);
      row_lno_t new_num_edge = 0;
      typedef Kokkos::RangePolicy<HandleExecSpace> my_exec_space;


      if (
#if defined( KOKKOS_HAVE_CUDA )
          Kokkos::Impl::is_same<Kokkos::Cuda, ExecutionSpace >::value ||
#endif
          0){


        int teamSizeMax = 0;
        int vector_size = 0;

        CountLowerTriangleTeam<row_index_view_type, nonzero_view_type, row_lno_temp_work_view_t> clt (nv, xadj, adj, lower_count);
        int max_allowed_team_size = team_policy_t::team_size_max(clt);
        KokkosKernels::Experimental::Util::get_suggested_vector_team_size<row_lno_t, HandleExecSpace>(
            max_allowed_team_size,
            vector_size,
            teamSizeMax,
            nv, ne);

        //std::cout << "teamSizeMax:" << teamSizeMax << " vector_size:" << vector_size << std::endl;
        //Kokkos::parallel_reduce(



        Kokkos::parallel_for(
            team_policy_t(nv / teamSizeMax + 1 , teamSizeMax, vector_size),
            clt//, new_num_edge
        );

        KokkosKernels::Experimental::Util::inclusive_parallel_prefix_sum<row_lno_temp_work_view_t, HandleExecSpace>
        (nv+1, lower_count);
        //Kokkos::parallel_scan (my_exec_space(0, nv + 1), PPS<row_lno_temp_work_view_t>(lower_count));
        HandleExecSpace::fence();
        auto lower_total_count = Kokkos::subview(lower_count, nv);
        auto hlower = Kokkos::create_mirror_view (lower_total_count);
        Kokkos::deep_copy (hlower, lower_total_count);

        new_num_edge = hlower();
        row_lno_persistent_work_view_t half_src = row_lno_persistent_work_view_t(Kokkos::ViewAllocateWithoutInitializing("HALF SRC"),new_num_edge);
        row_lno_persistent_work_view_t half_dst = row_lno_persistent_work_view_t(Kokkos::ViewAllocateWithoutInitializing("HALF DST"),new_num_edge);
        Kokkos::parallel_for(
            team_policy_t(nv / teamSizeMax + 1 , teamSizeMax, vector_size),
            FillLowerTriangleTeam
            <row_index_view_type, nonzero_view_type,
           row_lno_temp_work_view_t,row_lno_persistent_work_view_t> (nv, xadj, adj, lower_count, half_src, half_dst));

        src = lower_triangle_src = half_src;
        dst = lower_triangle_dst = half_dst;
        num_out_edges = size_of_edge_list = new_num_edge;
      }
      else {
        if (nv > 0) {
          Kokkos::parallel_reduce(my_exec_space(0,nv),
              CountLowerTriangle<row_index_view_type, nonzero_view_type, row_lno_temp_work_view_t> (nv, xadj, adj, lower_count), new_num_edge);
        }

        //Kokkos::parallel_scan (my_exec_space(0, nv + 1), PPS<row_lno_temp_work_view_t>(lower_count));

        KokkosKernels::Experimental::Util::inclusive_parallel_prefix_sum<row_lno_temp_work_view_t, HandleExecSpace>
        (nv+1, lower_count);
        row_lno_persistent_work_view_t half_src = row_lno_persistent_work_view_t(Kokkos::ViewAllocateWithoutInitializing("HALF SRC"),new_num_edge);
        row_lno_persistent_work_view_t half_dst = row_lno_persistent_work_view_t(Kokkos::ViewAllocateWithoutInitializing("HALF DST"),new_num_edge);

        Kokkos::parallel_for(my_exec_space(0,nv), FillLowerTriangle
            <row_index_view_type, nonzero_view_type,
            row_lno_temp_work_view_t,row_lno_persistent_work_view_t> (nv, xadj, adj, lower_count, half_src, half_dst));

        src = lower_triangle_src = half_src;
        dst = lower_triangle_dst = half_dst;
        num_out_edges = size_of_edge_list = new_num_edge;
      }

    }
  }

  struct ReduceMaxFunctor{
    color_view_t colors;
    ReduceMaxFunctor(color_view_t cat):colors(cat){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const row_lno_t &i, color_t & color_max) const {
      if (color_max < colors(i) ) color_max = colors(i);
    }

    KOKKOS_INLINE_FUNCTION
    void join (volatile color_t& dst , const volatile color_t& src) const { // max -plus semiring equivalent of "plus"
      if (dst < src) {
        dst = src;
      }
    }

    KOKKOS_INLINE_FUNCTION
    void init (color_t& dst) const {
      dst = 0;
    }
  };


  row_lno_t get_num_colors(){
    if (num_colors == 0){
      typedef typename Kokkos::RangePolicy<ExecutionSpace> my_exec_space;
      Kokkos::parallel_reduce(my_exec_space(0, vertex_colors.dimension_0()),
          ReduceMaxFunctor(vertex_colors) ,num_colors);
    }
    return num_colors;
  }

  /** \brief Sets Default Parameter settings for the given algorithm.
   */
  void set_defaults (const ColoringAlgorithm &col_algo){
    switch (col_algo){
    case COLORING_VB:
    case COLORING_VBBIT:
    case COLORING_VBCS:
    case COLORING_SERIAL:
      this->conflict_list_type = COLORING_ATOMIC;
      this->min_reduction_for_conflictlist = 0.35;
      this->min_elements_for_conflictlist = 1000;
      this->serial_conflict_resolution = false;
      this->tictoc = false;
      this->vb_edge_filtering = false;
      this->vb_chunk_size = 8;
      this->max_number_of_iterations = 200;
      this->eb_num_initial_colors = 1;
      break;
    case COLORING_EB:
      this->conflict_list_type = COLORING_PPS;
      this->min_reduction_for_conflictlist = 0.35;
      this->min_elements_for_conflictlist = 5000;
      this->serial_conflict_resolution = false;
      this->tictoc = false;
      this->vb_edge_filtering = false;
      this->vb_chunk_size = 8;
      this->max_number_of_iterations = 20000;
      this->eb_num_initial_colors = 1;
      break;
    default:
      std::cerr << "Unknown Coloring Algorithm" << std::endl;
      break;
    }
  }

  virtual ~GraphColoringHandle(){};

  //getters
  ColoringAlgorithm get_coloring_algo_type() const {return this->coloring_algorithm_type;}
  ConflictList get_conflict_list_type() const {return this->conflict_list_type;}
  double get_min_reduction_for_conflictlist() const{return this->min_reduction_for_conflictlist;}
  int get_min_elements_for_conflictlist() const{ return this->min_elements_for_conflictlist;}
  bool get_serial_conflict_resolution() const{return this->serial_conflict_resolution;}
  bool get_tictoc() const{return this->tictoc;}
  bool get_vb_edge_filtering() const{return this-vb_edge_filtering;}
  int get_vb_chunk_size() const{return this->vb_chunk_size;}
  int get_max_number_of_iterations() const{return this->max_number_of_iterations;}
  int get_eb_num_initial_colors() const{return this->eb_num_initial_colors;}

  double get_overall_coloring_time() const { return this->overall_coloring_time;}
  double get_coloring_time() const { return this->coloring_time;}
  int get_num_phases() const { return this->num_phases;}
  color_view_t get_vertex_colors() const {return this->vertex_colors;}
  bool is_coloring_called() const {return this->is_coloring_called_before;}
  //setters
  void set_coloring_algo_type(const ColoringAlgorithm &col_algo){this->coloring_algorithm_type = col_algo;}
  void set_conflict_list_type(const ConflictList &cl){this->conflict_list_type = cl;}
  void set_min_reduction_for_conflictlist(const double &min_reduction){this->min_reduction_for_conflictlist = min_reduction;}
  void set_min_elements_for_conflictlist(const int &min_elements){ this->min_elements_for_conflictlist = min_elements;}
  void set_serial_conflict_resolution(const bool &use_serial_conflist_resolution){this->serial_conflict_resolution = use_serial_conflist_resolution;}
  void set_tictoc(const bool use_tictoc){this->tictoc = use_tictoc;}
  void set_vb_edge_filtering(const bool  &use_vb_edge_filtering){this-vb_edge_filtering = use_vb_edge_filtering;}
  void set_vb_chunk_size(const int &chunksize){this->vb_chunk_size = chunksize;}
  void set_max_number_of_iterations(const int &max_phases){this->max_number_of_iterations = max_phases;}
  void set_eb_num_initial_colors(const int &num_initial_colors){this->eb_num_initial_colors = num_initial_colors;}
  void add_to_overall_coloring_time(const double &coloring_time_){this->overall_coloring_time += coloring_time_;}
  void set_coloring_time(const double &coloring_time_){this->coloring_time = coloring_time_;}
  void set_num_phases(const double &num_phases_){this->num_phases = num_phases_;}
  void set_vertex_colors( const color_view_t vertex_colors_){
    this->vertex_colors = vertex_colors_;
    this->is_coloring_called_before = true;
    this->num_colors = 0;
  }


};
}
}
}

#endif // _GRAPHCOLORHANDLE_HPP
