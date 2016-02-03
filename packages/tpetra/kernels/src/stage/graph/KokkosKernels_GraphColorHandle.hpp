
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

template <class row_index_view_type_, class nonconst_color_array_type_, class nonzero_index_view_type_,
      class ExecutionSpace, class TemporaryMemorySpace, class PersistentMemorySpace>
class GraphColoringHandle{

public:
  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;


  typedef row_index_view_type_ in_row_index_view_type;
  //typedef row_index_view_type_ idx_array_type;

  //typedef nonzero_index_view_type_ idx_edge_array_type;
  typedef nonzero_index_view_type_ in_nonzero_index_view_type;

  //typedef typename row_index_view_type::value_type idx;
  typedef typename in_row_index_view_type::non_const_value_type row_index_type;

  //typedef typename row_index_view_type::array_layout idx_array_layout;
  typedef typename in_row_index_view_type::array_layout row_view_array_layout;

  //typedef typename row_index_view_type::device_type idx_device_type;
  typedef typename in_row_index_view_type::device_type row_view_device_type;

  //typedef typename row_index_view_type::memory_traits idx_memory_traits;
  typedef typename in_row_index_view_type::memory_traits row_view_memory_traits;

  //typedef typename row_index_view_type::HostMirror host_view_type;
  typedef typename in_row_index_view_type::HostMirror row_host_view_type; //Host view type
  //typedef typename idx_memory_traits::MemorySpace MyMemorySpace;

  //typedef typename nonzero_index_view_type::non_const_value_type idx_edge;
  typedef typename in_nonzero_index_view_type::non_const_value_type nonzero_index_type;

  //typedef typename nonzero_index_view_type::array_layout idx_edge_array_layout;
  typedef typename in_nonzero_index_view_type::array_layout nonzero_index_view_array_layout;

  //typedef typename nonzero_index_view_type::device_type idx_edge_device_type;
  typedef typename in_nonzero_index_view_type::device_type nonzero_index_view_device_type;

  //typedef typename nonzero_index_view_type::memory_traits idx_edge_memory_traits;
  typedef typename in_nonzero_index_view_type::memory_traits nonzero_index_view_memory_traits;

  //typedef typename nonzero_index_view_type::HostMirror host_edge_view_type; //Host view type
  typedef typename in_nonzero_index_view_type::HostMirror nonzero_index_host_view_type; //Host view type



  typedef nonconst_color_array_type_ color_view_type;

  typedef typename color_view_type::non_const_value_type color_type;
  typedef typename color_view_type::array_layout color_view_array_layout;
  typedef typename color_view_type::device_type color_view_device_type;
  typedef typename color_view_type::memory_traits color_view_memory_traits;
  typedef typename color_view_type::HostMirror color_host_view_type; //Host view type



  typedef typename in_row_index_view_type::const_data_type const_row_data_type;
  typedef typename in_row_index_view_type::non_const_data_type non_const_row_data_type;
  typedef typename in_row_index_view_type::memory_space row_view_memory_space;
  typedef typename Kokkos::View<const_row_data_type, row_view_array_layout,
      row_view_memory_space, row_view_memory_traits> const_row_index_view_type;
  typedef typename Kokkos::View<non_const_row_data_type, row_view_array_layout,
      row_view_memory_space, row_view_memory_traits> non_const_row_index_view_type;



  typedef typename in_nonzero_index_view_type::const_data_type const_nonzero_index_data_type;
  typedef typename in_nonzero_index_view_type::non_const_data_type non_const_nonzero_index_data_type;
  typedef typename in_nonzero_index_view_type::memory_space nonzero_index_view_memory_space;
  typedef typename Kokkos::View<const_nonzero_index_data_type, nonzero_index_view_array_layout,
      nonzero_index_view_memory_space, nonzero_index_view_memory_traits> const_nonzero_index_view_type;
  typedef typename Kokkos::View<non_const_nonzero_index_data_type, nonzero_index_view_array_layout,
      nonzero_index_view_memory_space, nonzero_index_view_memory_traits> non_const_nonzero_index_view_type;




  //typedef typename Kokkos::View<row_index_type *, HandleTempMemorySpace> idx_temp_work_array_type;
  typedef typename Kokkos::View<row_index_type *, HandleTempMemorySpace> row_index_temp_work_view_type;
  //typedef typename row_index_persistent_work_view_type idx_persistent_work_array_type;
  typedef typename Kokkos::View<row_index_type *, PersistentMemorySpace> row_index_persistent_work_view_type;

  //typedef typename row_index_persistent_work_view_type::HostMirror host_idx_persistent_view_type; //Host view type
  typedef typename row_index_persistent_work_view_type::HostMirror row_index_persistent_host_view_type; //Host view type

private:
  //Parameters
  ColoringAlgorithm coloring_type; //VB, VBBIT or EB.
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


  row_index_type size_of_edge_list;
  row_index_persistent_work_view_type lower_triangle_src;
  row_index_persistent_work_view_type lower_triangle_dst;

  color_view_type vertex_colors;
  bool is_coloring_called_before;
  row_index_type num_colors;



  public:


  /**
   * \brief Default constructor.
   */
  GraphColoringHandle():
    coloring_type(COLORING_DEFAULT),
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
    this->set_defaults(this->coloring_type);
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
      this->coloring_type = col_algo;
    }
    if (set_default_parameters){
      this->set_defaults(this->coloring_type);
    }
  }

  /** \brief Chooses best algorithm based on the execution space. COLORING_EB if cuda, COLORING_VB otherwise.
   */
  void choose_default_algorithm(){

#if defined( KOKKOS_HAVE_SERIAL )
    if (Kokkos::Impl::is_same< Kokkos::Serial , ExecutionSpace >::value){
      this->coloring_type = COLORING_SERIAL;
#ifdef VERBOSE
      std::cout << "Serial Execution Space, Default Algorithm: COLORING_VB" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_PTHREAD )
    if (Kokkos::Impl::is_same< Kokkos::Threads , ExecutionSpace >::value){
      this->coloring_type = COLORING_VB;
#ifdef VERBOSE
      std::cout << "PTHREAD Execution Space, Default Algorithm: COLORING_VB" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_OPENMP )
    if (Kokkos::Impl::is_same< Kokkos::OpenMP, ExecutionSpace >::value){
      this->coloring_type = COLORING_VB;
#ifdef VERBOSE
      std::cout << "OpenMP Execution Space, Default Algorithm: COLORING_VB" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_CUDA )
    if (Kokkos::Impl::is_same<Kokkos::Cuda, ExecutionSpace >::value){
      this->coloring_type = COLORING_EB;
#ifdef VERBOSE
      std::cout << "Qthread Execution Space, Default Algorithm: COLORING_VB" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_QTHREAD)
    if (Kokkos::Impl::is_same< Kokkos::Qthread, ExecutionSpace >::value){
      this->coloring_type = COLORING_VB;
#ifdef VERBOSE
      std::cout << "Qthread Execution Space, Default Algorithm: COLORING_VB" << std::endl;
#endif
    }
#endif
  }

  template<typename v1, typename v2, typename v3>
  struct CountLowerTriangle{
    row_index_type nv;
    v1 xadj;
    v2 adj;
    v3 lower_xadj_counts;

    CountLowerTriangle(
        row_index_type nv_,
        v1 xadj_,
        v2 adj_,
        v3 lower_xadj_counts_
        ): nv(nv_),
            xadj(xadj_), adj(adj_),
            lower_xadj_counts(lower_xadj_counts_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const row_index_type &i, row_index_type &new_num_edge) const {
      row_index_type xadj_begin = xadj(i);
      row_index_type xadj_end = xadj(i + 1);

      row_index_type new_edge_count = 0;
      for (row_index_type j = xadj_begin; j < xadj_end; ++j){
        row_index_type n = adj(j);
        if (i < n && n < nv){
          new_edge_count += 1;
        }
      }
      lower_xadj_counts(i + 1) = new_edge_count;
      new_num_edge += new_edge_count;
    }
  };

  template <typename view_type>
  struct PPS{
    view_type lower_xadj_counts;
    PPS(view_type lower_xadj_counts_):
        lower_xadj_counts(lower_xadj_counts_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const row_index_type &ii, size_t& update, const bool final) const{
      update += lower_xadj_counts(ii);
      if (final) {
        lower_xadj_counts(ii)  = update;
      }
    }
  };

  template<typename v1, typename v2, typename v3, typename v4>
  struct FillLowerTriangle{
    row_index_type nv;
    v1 xadj;
    v2 adj;
    v3 lower_xadj_counts;
    v4 lower_srcs;
    v4 lower_dsts;

    FillLowerTriangle(
        row_index_type nv_,
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
    void operator()(const row_index_type &i) const{

      row_index_type xadj_begin = xadj[i];
      row_index_type xadj_end = xadj[i + 1];

      for (row_index_type j = xadj_begin; j < xadj_end; ++j){
        row_index_type n = adj(j);
        if (i < n && n < nv){
          row_index_type position = lower_xadj_counts(i)++;
          lower_srcs(position) = i;
          lower_dsts(position) = n;
        }
      }
    }
  };

  template <typename row_index_view_type, typename nonzero_view_type>
  void symmetrize_and_calculate_lower_diagonal_edge_list(
      row_index_type nv,
      row_index_view_type xadj, nonzero_view_type adj){

    KokkosKernels::Experimental::Util::symmetrize_and_get_lower_diagonal_edge_list
    <row_index_view_type, nonzero_view_type, row_index_persistent_work_view_type, HandleExecSpace>
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
      row_index_type nv, row_index_type ne,
      row_index_view_type xadj, nonzero_view_type adj,
      row_index_type  &num_out_edges,
      row_index_persistent_work_view_type &src,
      row_index_persistent_work_view_type &dst){

    if (size_of_edge_list > 0){
      num_out_edges = size_of_edge_list;
      //src = Kokkos::View<idx *, HandlePersistentMemorySpace> (this->lower_triangle_src);
      //dst = Kokkos::View<idx *, HandlePersistentMemorySpace> (this->lower_triangle_dst);
      src = (this->lower_triangle_src);
      dst = (this->lower_triangle_dst);
    }
    else {

      row_index_temp_work_view_type lower_count("LowerXADJ", nv + 1);
      row_index_type new_num_edge = 0;
      typedef Kokkos::RangePolicy<HandleExecSpace> my_exec_space;
      if (nv > 0) {
        Kokkos::parallel_reduce(my_exec_space(0,nv),
            CountLowerTriangle<row_index_view_type, nonzero_view_type, row_index_temp_work_view_type> (nv, xadj, adj, lower_count), new_num_edge);
      }

      //std::cout << "nv:" << nv << " ne:" << ne << " new_num_edge:" << new_num_edge << std::endl;

      row_index_persistent_work_view_type half_src = row_index_persistent_work_view_type("HALF SRC",new_num_edge);
      row_index_persistent_work_view_type half_dst = row_index_persistent_work_view_type("HALF DST",new_num_edge);
      Kokkos::parallel_scan (my_exec_space(0, nv + 1), PPS<row_index_temp_work_view_type>(lower_count));
      Kokkos::parallel_for(my_exec_space(0,nv), FillLowerTriangle
          <row_index_view_type, nonzero_view_type,
          row_index_temp_work_view_type,row_index_persistent_work_view_type> (nv, xadj, adj, lower_count, half_src, half_dst));
      src = lower_triangle_src = half_src;
      dst = lower_triangle_dst = half_dst;
      num_out_edges = size_of_edge_list = new_num_edge;
    }
  }

  struct ReduceMaxFunctor{
    color_view_type colors;
    ReduceMaxFunctor(color_view_type cat):colors(cat){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const row_index_type &i, color_type & color_max) const {
      if (color_max < colors(i) ) color_max = colors(i);
    }

    KOKKOS_INLINE_FUNCTION
    void join (volatile color_type& dst , const volatile color_type& src) const { // max -plus semiring equivalent of "plus"
      if (dst < src) {
        dst = src;
      }
    }

    KOKKOS_INLINE_FUNCTION
    void init (color_type& dst) const {
      dst = 0;
    }
  };


  row_index_type get_num_colors(){
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
  ColoringAlgorithm get_coloring_type() const {return this->coloring_type;}
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
  color_view_type get_vertex_colors() const {return this->vertex_colors;}
  bool is_coloring_called() const {return this->is_coloring_called_before;}
  //setters
  void set_coloring_type(const ColoringAlgorithm &col_algo){this->coloring_type = col_algo;}
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
  void set_vertex_colors( const color_view_type vertex_colors_){
    this->vertex_colors = vertex_colors_;
    this->is_coloring_called_before = true;
    this->num_colors = 0;
  }


};
}
}
}

#endif // _GRAPHCOLORHANDLE_HPP
