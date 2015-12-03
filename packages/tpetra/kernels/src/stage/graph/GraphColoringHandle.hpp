
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_Core.hpp>

#ifndef _GRAPHCOLORHANDLE_HPP
#define _GRAPHCOLORHANDLE_HPP
//#define VERBOSE
namespace KokkosKernels{

namespace Experimental{


namespace Graph{

enum ColoringAlgorithm{COLORING_DEFAULT, COLORING_SERIAL, COLORING_VB, COLORING_VBBIT, COLORING_VBCS, COLORING_EB};

enum ConflictList{COLORING_NOCONFLICT, COLORING_ATOMIC, COLORING_PPS};

template <class idx_array_type_, class color_array_type_, class idx_edge_array_type_,
      class ExecutionSpace, class TemporaryMemorySpace, class PersistentMemorySpace>
class GraphColoringHandle{

public:
  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;

  typedef idx_array_type_ idx_array_type;
  typedef idx_edge_array_type_ idx_edge_array_type;
  typedef color_array_type_ color_array_type;


  typedef typename idx_array_type::value_type idx;
  typedef typename idx_array_type::array_layout idx_array_layout;
  typedef typename idx_array_type::device_type idx_device_type;
  typedef typename idx_array_type::memory_traits idx_memory_traits;
  typedef typename idx_array_type::HostMirror host_view_type; //Host view type
  //typedef typename idx_memory_traits::MemorySpace MyMemorySpace;

  typedef typename idx_edge_array_type::value_type idx_edge;
  typedef typename idx_edge_array_type::array_layout idx_edge_array_layout;
  typedef typename idx_edge_array_type::device_type idx_edge_device_type;
  typedef typename idx_edge_array_type::memory_traits idx_edge_memory_traits;
  typedef typename idx_edge_array_type::HostMirror host_edge_view_type; //Host view type
  //typedef typename idx_edge_memory_traits::MemorySpace MyEdgeMemorySpace;

  typedef typename color_array_type::value_type color_type;
  typedef typename color_array_type::array_layout color_type_array_layout;
  typedef typename color_array_type::device_type color_type_device_type;
  typedef typename color_array_type::memory_traits color_type_memory_traits;
  typedef typename color_array_type::HostMirror host_color_view_type; //Host view type


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


  idx size_of_edge_list;
  Kokkos::View<idx *, HandlePersistentMemorySpace> lower_triangle_src;
  Kokkos::View<idx *, HandlePersistentMemorySpace> lower_triangle_dst;

  color_array_type vertex_colors;
  bool is_coloring_called_before;
  idx num_colors;



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
    v1 xadj;
    v2 adj;
    v3 lower_xadj_counts;

    CountLowerTriangle(
        v1 xadj_,
        v2 adj_,
        v3 lower_xadj_counts_
        ): xadj(xadj_), adj(adj_), lower_xadj_counts(lower_xadj_counts_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &i, idx &new_num_edge) const {
      idx xadj_begin = xadj(i);
      idx xadj_end = xadj(i + 1);

      idx new_edge_count = 0;
      for (idx j = xadj_begin; j < xadj_end; ++j){
        idx n = adj(j);
        if (i < n){
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
    void operator()(const idx &ii, size_t& update, const bool final) const{
      update += lower_xadj_counts(ii);
      if (final) {
        lower_xadj_counts(ii)  = update;
      }
    }
  };

  template<typename v1, typename v2, typename v3, typename v4>
  struct FillLowerTriangle{
    v1 xadj;
    v2 adj;
    v3 lower_xadj_counts;
    v4 lower_srcs;
    v4 lower_dsts;

    FillLowerTriangle(
        v1 xadj_,
        v2 adj_,
        v3 lower_xadj_counts_,
        v4 lower_srcs_,
        v4 lower_dsts_
        ): xadj(xadj_), adj(adj_), lower_xadj_counts(lower_xadj_counts_),
            lower_srcs(lower_srcs_), lower_dsts(lower_dsts_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &i) const{

      idx xadj_begin = xadj[i];
      idx xadj_end = xadj[i + 1];

      for (idx j = xadj_begin; j < xadj_end; ++j){
        idx n = adj(j);
        if (i < n){
          idx position = lower_xadj_counts(i)++;
          lower_srcs(position) = i;
          lower_dsts(position) = n;
        }
      }
    }
  };


  void get_lower_diagonal_edge_list(
      idx nv, idx ne,
      idx_array_type xadj, idx_edge_array_type adj,
      idx  &num_out_edges,
      Kokkos::View<idx *, HandlePersistentMemorySpace> &src,
      Kokkos::View<idx *, HandlePersistentMemorySpace> &dst){

    if (size_of_edge_list > 0){
      num_out_edges = size_of_edge_list;
      //src = Kokkos::View<idx *, HandlePersistentMemorySpace> (this->lower_triangle_src);
      //dst = Kokkos::View<idx *, HandlePersistentMemorySpace> (this->lower_triangle_dst);
      src = (this->lower_triangle_src);
      dst = (this->lower_triangle_dst);
    }
    else {
      typedef typename Kokkos::View<idx *, HandleTempMemorySpace> idx_temp_view;
      typedef typename Kokkos::RangePolicy<ExecutionSpace> my_exec_space;


      idx_temp_view lower_count("LowerXADJ", nv + 1);
      idx new_num_edge = 0;
      Kokkos::parallel_reduce(my_exec_space(0,nv), CountLowerTriangle<idx_array_type, idx_edge_array_type, idx_temp_view>(xadj, adj, lower_count), new_num_edge);

      Kokkos::View<idx*, HandlePersistentMemorySpace> half_src =
          Kokkos::View<idx *, HandlePersistentMemorySpace>("HALF SRC",new_num_edge);
      Kokkos::View<idx *, HandlePersistentMemorySpace> half_dst =
          Kokkos::View<idx *, HandlePersistentMemorySpace>("HALF DST",new_num_edge);
      Kokkos::parallel_scan (my_exec_space(0, nv + 1), PPS<idx_temp_view>(lower_count));
      Kokkos::parallel_for(my_exec_space(0,nv),
          FillLowerTriangle<idx_array_type, idx_edge_array_type, idx_temp_view, Kokkos::View<idx*, HandlePersistentMemorySpace> >
          (xadj, adj, lower_count, half_src, half_dst));
      src = lower_triangle_src = half_src;
      dst = lower_triangle_dst = half_dst;
      num_out_edges = size_of_edge_list = new_num_edge;
    }
  }

  struct ReduceMaxFunctor{
    color_array_type colors;
    ReduceMaxFunctor(color_array_type cat):colors(cat){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &i, color_type & color_max) const {
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


  idx get_num_colors(){
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
  color_array_type get_vertex_colors() const {return this->vertex_colors;}
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
  void set_vertex_colors( const color_array_type vertex_colors_){
    this->vertex_colors = vertex_colors_;
    this->is_coloring_called_before = true;
    this->num_colors = 0;
  }


};
}
}
}

#endif // _GRAPHCOLORHANDLE_HPP
