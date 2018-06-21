/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_Utils.hpp>

#ifndef _GRAPHCOLORHANDLE_HPP
#define _GRAPHCOLORHANDLE_HPP

//#define VERBOSE
namespace KokkosGraph {

enum ColoringAlgorithm { COLORING_DEFAULT,
                         COLORING_SERIAL,
                         COLORING_VB,
                         COLORING_VBBIT,
                         COLORING_VBCS,
                         COLORING_EB,
                         COLORING_SERIAL2,
                         COLORING_SPGEMM,
                         COLORING_D2_MATRIX_SQUARED,          // Distance-2 Graph Coloring (Brian's Code)
                         COLORING_D2                          // Distance-2 Graph Coloring (WCMCLEN)
                       };

enum ConflictList{COLORING_NOCONFLICT, COLORING_ATOMIC, COLORING_PPS};

enum ColoringType {Distance1, Distance2};

template <class size_type_, class color_t_, class lno_t_, 
         //class lno_row_view_t_, class nonconst_color_view_t_, class lno_nnz_view_t_,
          class ExecutionSpace, class TemporaryMemorySpace, class PersistentMemorySpace>
class GraphColoringHandle
{

public:
  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;


  typedef typename std::remove_const<size_type_>::type  size_type;
  typedef const size_type const_size_type;

  typedef typename std::remove_const<lno_t_>::type  nnz_lno_t;
  typedef const nnz_lno_t const_nnz_lno_t;

  typedef typename std::remove_const<color_t_>::type  color_t;
  typedef const color_t const_color_t;



  typedef typename Kokkos::View<color_t *, HandlePersistentMemorySpace> color_view_t;

  typedef typename color_view_t::array_layout color_view_array_layout;
  typedef typename color_view_t::device_type color_view_device_t;
  typedef typename color_view_t::memory_traits color_view_memory_traits;
  typedef typename color_view_t::HostMirror color_host_view_t; //Host view type


  typedef typename Kokkos::View<size_type *, HandleTempMemorySpace> size_type_temp_work_view_t;
  typedef typename Kokkos::View<size_type *, HandlePersistentMemorySpace> size_type_persistent_work_view_t;

  typedef typename size_type_persistent_work_view_t::HostMirror size_type_persistent_work_host_view_t; //Host view type

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


  size_type size_of_edge_list;
  nnz_lno_persistent_work_view_t lower_triangle_src;
  nnz_lno_persistent_work_view_t lower_triangle_dst;

  color_view_t vertex_colors;
  bool is_coloring_called_before;
  nnz_lno_t num_colors;



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
  void choose_default_algorithm()
  {
#if defined( KOKKOS_ENABLE_SERIAL )
    if (Kokkos::Impl::is_same< Kokkos::Serial , ExecutionSpace >::value){
      this->coloring_algorithm_type = COLORING_SERIAL;
#ifdef VERBOSE
      std::cout << "Serial Execution Space, Default Algorithm: COLORING_VB" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_ENABLE_THREADS )
    if (Kokkos::Impl::is_same< Kokkos::Threads , ExecutionSpace >::value){
      this->coloring_algorithm_type = COLORING_VB;
#ifdef VERBOSE
      std::cout << "PTHREAD Execution Space, Default Algorithm: COLORING_VB" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_ENABLE_OPENMP )
    if (Kokkos::Impl::is_same< Kokkos::OpenMP, ExecutionSpace >::value){
      this->coloring_algorithm_type = COLORING_VB;
#ifdef VERBOSE
      std::cout << "OpenMP Execution Space, Default Algorithm: COLORING_VB" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_ENABLE_CUDA )
    if (Kokkos::Impl::is_same<Kokkos::Cuda, ExecutionSpace >::value){
      this->coloring_algorithm_type = COLORING_EB;
#ifdef VERBOSE
      std::cout << "Cuda Execution Space, Default Algorithm: COLORING_VB" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_ENABLE_QTHREAD)
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
    nnz_lno_t nv;
    v1 xadj;
    v2 adj;
    v3 lower_xadj_counts;

    CountLowerTriangle(
        nnz_lno_t nv_,
        v1 xadj_,
        v2 adj_,
        v3 lower_xadj_counts_
        ): nv(nv_),
            xadj(xadj_), adj(adj_),
            lower_xadj_counts(lower_xadj_counts_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t &i, size_type &new_num_edge) const {
      size_type xadj_begin = xadj(i);
      size_type xadj_end = xadj(i + 1);

      size_type new_edge_count = 0;
      for (size_type j = xadj_begin; j < xadj_end; ++j){
        nnz_lno_t n = adj(j);
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

    nnz_lno_t nv;
    v1 xadj;
    v2 adj;
    v3 lower_xadj_counts;

    CountLowerTriangleTeam(
        nnz_lno_t nv_,
        v1 xadj_,
        v2 adj_,
        v3 lower_xadj_counts_
        ): nv(nv_),
            xadj(xadj_), adj(adj_),
            lower_xadj_counts(lower_xadj_counts_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t & teamMember/*, row_lno_t &new_num_edge*/) const {


      nnz_lno_t ii = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (ii >= nv) {
        return;
      }

      size_type xadj_begin = xadj(ii);
      size_type xadj_end = xadj(ii + 1);

      size_type new_edge_count = 0;


      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(teamMember, xadj_end - xadj_begin),
          [&] (size_type i, size_type &numEdges) {

        size_type adjind = i + xadj_begin;
        nnz_lno_t n = adj[adjind];
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
    nnz_lno_t nv;
    v1 xadj;
    v2 adj;
    v3 lower_xadj_counts;
    v4 lower_srcs;
    v4 lower_dsts;

    FillLowerTriangleTeam(
        nnz_lno_t nv_,
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
      typedef typename std::remove_reference< decltype( lower_xadj_counts(0) ) >::type atomic_incr_type;

      nnz_lno_t ii = teamMember.league_rank()  * teamMember.team_size()+ teamMember.team_rank();
      if (ii >= nv) {
        return;
      }

      size_type xadj_begin = xadj(ii);
      size_type xadj_end = xadj(ii + 1);

      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(teamMember, xadj_end - xadj_begin),
          [&] (size_type i) {

        size_type adjind = i + xadj_begin;
        nnz_lno_t n = adj[adjind];
        if (ii < n && n < nv){
          size_type position =
              Kokkos::atomic_fetch_add( &(lower_xadj_counts(ii)), atomic_incr_type(1));
          lower_srcs(position) = ii;
          lower_dsts(position) = n;
        }
      });
    }
  };



  template<typename v1, typename v2, typename v3, typename v4>
  struct FillLowerTriangle{
    nnz_lno_t nv;
    v1 xadj;
    v2 adj;
    v3 lower_xadj_counts;
    v4 lower_srcs;
    v4 lower_dsts;

    FillLowerTriangle(
        nnz_lno_t nv_,
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
    void operator()(const nnz_lno_t &i) const{

      size_type xadj_begin = xadj[i];
      size_type xadj_end = xadj[i + 1];

      for (size_type j = xadj_begin; j < xadj_end; ++j){
        nnz_lno_t n = adj(j);
        if (i < n && n < nv){
          size_type position = lower_xadj_counts(i)++;
          lower_srcs(position) = i;
          lower_dsts(position) = n;
        }
      }
    }
  };



  template <typename row_index_view_type, typename nonzero_view_type>
  void symmetrize_and_calculate_lower_diagonal_edge_list(
      nnz_lno_t nv,
      row_index_view_type xadj, nonzero_view_type adj){

    KokkosKernels::Impl::symmetrize_and_get_lower_diagonal_edge_list
    <row_index_view_type, nonzero_view_type, nnz_lno_persistent_work_view_t, HandleExecSpace>
      (
        nv,
        xadj,
        adj,
        lower_triangle_src,
        lower_triangle_dst);

    size_of_edge_list = lower_triangle_src.extent(0);

  }



  template <typename row_index_view_type, typename nonzero_view_type>
  void get_lower_diagonal_edge_list(
      nnz_lno_t nv, size_type ne,
      row_index_view_type xadj, nonzero_view_type adj,
      size_type  &num_out_edges,
      nnz_lno_persistent_work_view_t &src,
      nnz_lno_persistent_work_view_t &dst){

    if (size_of_edge_list > 0){
      num_out_edges = size_of_edge_list;
      //src = Kokkos::View<idx *, HandlePersistentMemorySpace> (this->lower_triangle_src);
      //dst = Kokkos::View<idx *, HandlePersistentMemorySpace> (this->lower_triangle_dst);
      src = (this->lower_triangle_src);
      dst = (this->lower_triangle_dst);
    }
    else {

      size_type_temp_work_view_t lower_count("LowerXADJ", nv + 1);
      size_type new_num_edge = 0;
      typedef Kokkos::RangePolicy<HandleExecSpace> my_exec_space;

      if (
#if defined( KOKKOS_ENABLE_CUDA )
          Kokkos::Impl::is_same<Kokkos::Cuda, ExecutionSpace >::value ||
#endif
          0){


        int teamSizeMax = 0;
        int vector_size = 0;

        CountLowerTriangleTeam<row_index_view_type, nonzero_view_type, size_type_temp_work_view_t> clt (nv, xadj, adj, lower_count);
        int max_allowed_team_size = team_policy_t::team_size_max(clt);
        KokkosKernels::Impl::get_suggested_vector_team_size<size_type, HandleExecSpace>(
            max_allowed_team_size,
            vector_size,
            teamSizeMax,
            nv, ne);

        Kokkos::parallel_for("KokkosGraph::CountLowerTriangleTeam",
            team_policy_t(nv / teamSizeMax + 1 , teamSizeMax, vector_size),
            clt//, new_num_edge
        );

        KokkosKernels::Impl::inclusive_parallel_prefix_sum<size_type_temp_work_view_t, HandleExecSpace>
        (nv+1, lower_count);
        //Kokkos::parallel_scan (my_exec_space(0, nv + 1), PPS<row_lno_temp_work_view_t>(lower_count));
        HandleExecSpace::fence();
        auto lower_total_count = Kokkos::subview(lower_count, nv);
        auto hlower = Kokkos::create_mirror_view (lower_total_count);
        Kokkos::deep_copy (hlower, lower_total_count);

        new_num_edge = hlower();
        nnz_lno_persistent_work_view_t half_src (Kokkos::ViewAllocateWithoutInitializing("HALF SRC"),new_num_edge);
        nnz_lno_persistent_work_view_t half_dst (Kokkos::ViewAllocateWithoutInitializing("HALF DST"),new_num_edge);
        Kokkos::parallel_for("KokkosGraph::FillLowerTriangleTeam",
            team_policy_t(nv / teamSizeMax + 1 , teamSizeMax, vector_size),
            FillLowerTriangleTeam
            <row_index_view_type, nonzero_view_type,
           size_type_temp_work_view_t,nnz_lno_persistent_work_view_t> (nv, xadj, adj, lower_count, half_src, half_dst));

        src = lower_triangle_src = half_src;
        dst = lower_triangle_dst = half_dst;
        num_out_edges = size_of_edge_list = new_num_edge;
      }
      else {
        if (nv > 0) {
          Kokkos::parallel_reduce("KokkosGraph::CountLowerTriangleTeam",my_exec_space(0,nv),
              CountLowerTriangle<row_index_view_type, nonzero_view_type, size_type_temp_work_view_t> (nv, xadj, adj, lower_count), new_num_edge);
        }

        //Kokkos::parallel_scan (my_exec_space(0, nv + 1), PPS<row_lno_temp_work_view_t>(lower_count));

        KokkosKernels::Impl::inclusive_parallel_prefix_sum<size_type_temp_work_view_t, HandleExecSpace>
        (nv+1, lower_count);
        nnz_lno_persistent_work_view_t half_src (Kokkos::ViewAllocateWithoutInitializing("HALF SRC"),new_num_edge);
        nnz_lno_persistent_work_view_t half_dst (Kokkos::ViewAllocateWithoutInitializing("HALF DST"),new_num_edge);

        Kokkos::parallel_for("KokkosGraph::FillLowerTriangleTeam",my_exec_space(0,nv), FillLowerTriangle
            <row_index_view_type, nonzero_view_type,
            size_type_temp_work_view_t,nnz_lno_persistent_work_view_t> (nv, xadj, adj, lower_count, half_src, half_dst));

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
    void operator()(const nnz_lno_t &i, color_t & color_max) const {
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



  nnz_lno_t get_num_colors(){
    if (num_colors == 0){
      typedef typename Kokkos::RangePolicy<ExecutionSpace> my_exec_space;
      Kokkos::parallel_reduce("KokkosKernels::FindMax", my_exec_space(0, vertex_colors.extent(0)),
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
    case COLORING_SERIAL2:
    case COLORING_SPGEMM:
    case COLORING_D2_MATRIX_SQUARED:
    case COLORING_D2:
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
      throw std::runtime_error ("Unknown Coloring Algorithm\n");
      //break;
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
  bool get_vb_edge_filtering() const{return this->vb_edge_filtering;}
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
  void set_vb_edge_filtering(const bool  &use_vb_edge_filtering){this->vb_edge_filtering = use_vb_edge_filtering;}
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

}   // namespace KokkosGraph

#endif // _GRAPHCOLORHANDLE_HPP
