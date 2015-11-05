#include <Kokkos_Core.hpp>
#include <Kokkos_Atomic.hpp>
#include <impl/Kokkos_Timer.hpp>
#ifdef SKIPFORNOW
#include <Kokkos_Sort.hpp>
#endif
#include <Kokkos_MemoryTraits.hpp>
#include <GraphColoringHandle.hpp>

#ifndef _KOKKOSCOLORINGIMP_HPP
#define _KOKKOSCOLORINGIMP_HPP


#define EBCOLORING_HIGHER_QUALITY        //suggested

namespace Experimental{
namespace KokkosKernels{
namespace Graph{


namespace Impl{

#define VB_COLORING_FORBIDDEN_SIZE 64
#define VBBIT_COLORING_FORBIDDEN_SIZE 64
/*! \brief Base class for graph coloring purposes.
 *  Each color represents the set of the vertices that are independent,
 *  e.g. no vertex having same color shares an edge.
 *  General aim is to find the minimum number of colors, minimum number of independent sets.
 */
template <typename HandleType>
class GraphColor {
public:

  typedef typename HandleType::idx_array_type idx_array_type;
  typedef typename HandleType::idx_edge_array_type idx_edge_array_type;
  typedef typename HandleType::color_array_type color_array_type;

  typedef typename HandleType::idx idx;
  typedef typename HandleType::idx_array_layout idx_array_layout;
  typedef typename HandleType::idx_device_type idx_device_type;
  typedef typename HandleType::idx_memory_traits idx_memory_traits;
  typedef typename HandleType::host_view_type host_view_type; //Host view type
  //typedef typename idx_memory_traits::MemorySpace MyMemorySpace;

  typedef typename HandleType::idx_edge idx_edge;
  typedef typename HandleType::idx_edge_array_layout idx_edge_array_layout;
  typedef typename HandleType::idx_edge_device_type idx_edge_device_type;
  typedef typename HandleType::idx_edge_memory_traits idx_edge_memory_traits;
  typedef typename HandleType::host_edge_view_type host_edge_view_type; //Host view type
  //typedef typename idx_edge_memory_traits::MemorySpace MyEdgeMemorySpace;

  typedef typename HandleType::color_type color_type;
  typedef typename HandleType::color_type_array_layout color_type_array_layout;
  typedef typename HandleType::color_type_device_type color_type_device_type;
  typedef typename HandleType::color_type_memory_traits color_type_memory_traits;
  typedef typename HandleType::host_color_view_type host_color_view_type; //Host view type

  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;


protected:
  idx nv; //# vertices
  idx ne; //# edges
  idx_array_type xadj; //rowmap
  idx_edge_array_type adj; // entries
  idx_edge_array_type kok_src, kok_dst; //Edge list storage of the graph
  HandleType *cp;

public:
  /**
   * \brief GraphColor constructor.
   * \param nv_: number of vertices in the graph
   * \param ne_: number of edges in the graph
   * \param row_map: the xadj array of the graph. Its size is nv_ +1
   * \param entries: adjacency array of the graph. Its size is ne_
   * \param coloring_handle: GraphColoringHandle object that holds the specification about the graph coloring,
   *    including parameters.
   */
  GraphColor(
      idx nv_,
      idx ne_,
      idx_array_type row_map,
      idx_edge_array_type entries,
      HandleType *coloring_handle):
        nv (nv_), ne(ne_),xadj(row_map), adj (entries),
        kok_src(), kok_dst(), cp(coloring_handle){}

  /** \brief GraphColor destructor.
   */
  virtual ~GraphColor(){}


  /** \brief Function to color the vertices of the graphs. This is the base class,
   * therefore, it only performs sequential coloring on the host device, ignoring the execution space.
   * \param colors is the output array corresponding the color of each vertex.Size is this->nv.
   *   Attn: Color array must be nonnegative numbers. If there is no initial colors,
   *   it should be all initialized with zeros. Any positive value in the given array, will make the
   *   algorithm to assume that the color is fixed for the corresponding vertex.
   * \param num_phases: The number of iterations (phases) that algorithm takes to converge.
   */
  virtual void color_graph(
      color_array_type &d_colors,
      int &num_phases){

    num_phases = 1;

    host_color_view_type colors = Kokkos::create_mirror_view (d_colors);
    host_view_type h_xadj = Kokkos::create_mirror_view (this->xadj);
    host_edge_view_type h_adj = Kokkos::create_mirror_view (this->adj);
    Kokkos::deep_copy (h_xadj, this->xadj);
    Kokkos::deep_copy (h_adj, this->adj);

    MyExecSpace::fence();



    //create a ban color array to keep track of
    //which colors have been taking by the neighbor vertices.
    idx *banned_colors = new idx[this->nv];
    for (idx i = 0; i < this->nv; ++i) banned_colors[i] = 0;

    color_type max_color = 0;
    //traverse vertices greedily
    for (idx i = 0; i < this->nv; ++i){
      idx nbegin = h_xadj(i);
      idx nend = h_xadj(i + 1);
      //std::cout << "nb:" << nbegin << " ne:" << nend << std::endl;
      //check the colors of neighbors
      for (idx j = nbegin; j < nend; ++j){
        idx n = h_adj(j);
        //set the banned_color of the color of the neighbor vertex to my vertex index.
        //the entries in the banned_color array that has my vertex index will be the set of prohibeted colors.
        banned_colors[colors(n)] = i;
      }
      //check the prohibeted colors, and pick the first available one.
      for (color_type j = 1; j <= max_color; ++j) {
        if(banned_colors[j] != i){
          colors(i) = j;
          break;
        }
      }
      //if no color is available, pick a new color.
      if (colors(i) == 0) colors(i) = ++max_color;
    }
    delete [] banned_colors;

    Kokkos::deep_copy (d_colors, colors); // Copy from host to device.
  }

#ifdef SKIPFORNOW

  struct FILL_COLOR_XADJ{
    color_array_type _d_colors;
    idx_array_type _color_xadj;


    FILL_COLOR_XADJ(
        color_array_type d_colors,
        idx_array_type color_xadj):
          _d_colors(d_colors), _color_xadj(color_xadj){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &ii) const {
      color_type c = _d_colors[ii];
      const idx future_index = Kokkos::atomic_fetch_add( &(_color_xadj(c)), 1);

    }
  };
  struct XADJ_PPS{
    idx_array_type _color_xadj;


    XADJ_PPS(
        idx_array_type color_xadj):
          _color_xadj(color_xadj){}


    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii, size_t& update, const bool final) const {
      update += _color_xadj(ii);
      if (final) {
        _color_xadj(ii) = idx (update);
      }
    }
  };

  struct FILL_COLOR_ADJ{
    color_array_type _d_colors;
    idx_array_type _color_xadj;
    idx_array_type _color_adj;


    FILL_COLOR_ADJ(
        color_array_type d_colors,
        idx_array_type color_xadj,
        idx_array_type color_adj):
          _d_colors(d_colors), _color_xadj(color_xadj), _color_adj(color_adj){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &ii) const {
      color_type c = _d_colors[ii];
      if (c <= 0) std::cout << "ii:" << ii << " c:" << c << std::endl;
      const idx future_index = Kokkos::atomic_fetch_add( &(_color_xadj(c - 1)), 1);
      _color_adj(future_index) = ii;
    }
  };
  void create_color_adjacency(
      idx numColors,
      color_array_type d_colors,
      idx_array_type &color_xadj,
      idx_array_type &color_adj){
    typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;

    color_xadj = idx_array_type("COLOR_XADJ", numColors + 1);

    idx_array_type tmp_color_xadj (Kokkos::ViewAllocateWithoutInitializing("TMP_COLOR_XADJ"), numColors + 1);
    color_adj = idx_array_type(Kokkos::ViewAllocateWithoutInitializing("COLOR_ADJ"), this->nv);

    Kokkos::parallel_for (my_exec_space (0, this->nv) , FILL_COLOR_XADJ(d_colors, color_xadj));

    Kokkos::parallel_scan (my_exec_space (0, numColors + 1) , XADJ_PPS(color_xadj));

    Kokkos::deep_copy (tmp_color_xadj, color_xadj);

    Kokkos::parallel_for (my_exec_space (0, this->nv) , FILL_COLOR_ADJ(d_colors, tmp_color_xadj, color_adj));
  }


  //void setParameters(){};



  struct CountLowerTriangle{

    Kokkos::View<idx *, MyMemorySpace> xadj;

    Kokkos::View<idx *, MyEdgeMemorySpace> adj;
    Kokkos::View<idx *, MyMemorySpace> lower_xadj_counts;

    CountLowerTriangle(
        Kokkos::View<idx *, MyMemorySpace> xadj_,
        Kokkos::View<idx *, MyEdgeMemorySpace> adj_,
        Kokkos::View<idx *, MyMemorySpace> lower_xadj_counts_
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


  struct FillLowerTriangle{

    Kokkos::View<idx *, MyMemorySpace> xadj;
    Kokkos::View<idx *, MyEdgeMemorySpace> adj;
    Kokkos::View<idx *, MyMemorySpace> lower_xadj_counts;
    Kokkos::View<idx *, MyEdgeMemorySpace> lower_srcs;
    Kokkos::View<idx *, MyEdgeMemorySpace> lower_dsts;

    FillLowerTriangle(
        Kokkos::View<idx *, MyMemorySpace> xadj_,
        Kokkos::View<idx *, MyEdgeMemorySpace> adj_,
        Kokkos::View<idx *, MyMemorySpace> lower_xadj_counts_,
        Kokkos::View<idx *, MyEdgeMemorySpace> lower_srcs_,
        Kokkos::View<idx *, MyEdgeMemorySpace> lower_dsts_
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


  struct PPS{
    Kokkos::View<idx *, MyMemorySpace> lower_xadj_counts;
    PPS(Kokkos::View<idx *, MyMemorySpace> lower_xadj_counts_):
      lower_xadj_counts(lower_xadj_counts_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &ii, size_t& update, const bool final) const{
      update += lower_xadj_counts(ii);
      if (final) {
        lower_xadj_counts(ii)  = update;
      }
    }
  };

  void CreateHalfAdjacencyList(
      idx nv,
      idx ne,
      Kokkos::View<idx *, MyMemorySpace> xadj,
      Kokkos::View<idx *, MyEdgeMemorySpace> adj,
      Kokkos::View<idx *, MyEdgeMemorySpace> &half_src,
      Kokkos::View<idx *, MyEdgeMemorySpace> &half_dst
  ){
    typedef Kokkos::View<idx *, MyEdgeMemorySpace> idx_edge_array_type;
    typedef Kokkos::View<idx *, MyMemorySpace> idx_array_type;
    typedef typename idx_array_type::HostMirror host_view_type; //Host view type
    typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;

    idx_array_type lower_count("LowerXADJ", nv + 1);

    idx new_num_edge = 0;
    Kokkos::parallel_reduce(my_exec_space(0,nv), CountLowerTriangle(xadj, adj, lower_count), new_num_edge);
    std::cout << "NE:" << this->ne << " new_num_edge:" << new_num_edge * 2 << std::endl;
    this->ne = new_num_edge * 2;
    half_src = idx_array_type("HALF SRC", ne / 2);
    half_dst = idx_array_type("HALF dst", ne / 2);
    Kokkos::parallel_scan (my_exec_space(0, nv + 1), PPS(lower_count));
    Kokkos::parallel_for(my_exec_space(0,nv), FillLowerTriangle(xadj, adj, lower_count, half_src, half_dst));
  }
#endif
};


/*! \brief Class for the vertex based graph coloring algorithms.
 *  They work better on CPUs and Xeon Phis, but edge-based ones are better on GPUs.
 *  Includes 3 algorithms:
 *  VB: Speculative parallel vertex based algorithm using a forbidden array of size 64 per thread.
 *  Best on Xeon Phi among the vertex based algorithms.
 *  VBBIT: Speculative parallel vertex based using a long integer for forbidden colors per vertex.
 *  Best on GPUs among the vertex based algorithms.
 *  VBCS: Speculative parallel vertex based using color set implementation.
 */
template <typename HandleType>
class GraphColor_VB:public GraphColor <HandleType>{
public:


  typedef typename HandleType::idx_array_type idx_array_type;
  typedef typename HandleType::idx_edge_array_type idx_edge_array_type;
  typedef typename HandleType::color_array_type color_array_type;

  typedef long long int ban_type;

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

  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;

  typedef typename Kokkos::View<idx, idx_device_type> idx_type;
  typedef typename Kokkos::View<idx, Kokkos::MemoryUnmanaged> um_array_type;

  typedef typename idx_type::HostMirror idx_host_type; //Host view type
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;

  typedef typename Kokkos::View<idx *, idx_array_layout, MyTempMemorySpace> idx_temp_work_array_type;
  typedef typename Kokkos::View<idx *, idx_array_layout, MyPersistentMemorySpace> idx_persistent_work_array_type;
protected:



  //typedef Kokkos::View<idx, /*Kokkos::Serial::array_layout,*//* Kokkos::Serial,*/ Kokkos::MemoryUnmanaged> um_array_type;


  bool _serialConflictResolution; //if true use serial conflict resolution
  bool _ticToc; //if true print info in each step
  char _conflictlist; //0 for no conflictlist, 1 for atomic, 2 for pps

  double _pps_ratio; //the minimum number of reduction on the size of the conflictlist to create a new conflictlist
  idx _min_vertex_cut_off; //minimum number of vertices to reduce the conflictlist further.
  bool _edge_filtering; //if true, edge-filtering is applied by swaps on adjacency array.
  int _chunkSize; //the size of the minimum work unit assigned to threads. Changes the convergence on GPUs
  char _use_color_set; //the VB algorithm type.
                        // 0 for VB:
                        // 1: for VBCS
                        // 2: for VBBIT

  int _max_num_iterations;

public:
  /**
   * \brief GraphColor_VB constructor .
   * \param nv_: number of vertices in the graph
   * \param ne_: number of edges in the graph
   * \param row_map: the xadj array of the graph. Its size is nv_ +1
   * \param entries: adjacency array of the graph. Its size is ne_
   * \param coloring_handle: GraphColoringHandle object that holds the specification about the graph coloring,
   *    including parameters.
   */
  GraphColor_VB(idx nv_, idx ne_, idx_array_type row_map, idx_edge_array_type entries, HandleType *coloring_handle):
    GraphColor<HandleType>(nv_, ne_, row_map, entries, coloring_handle),
    _serialConflictResolution(coloring_handle->get_serial_conflict_resolution()),
    _ticToc(coloring_handle->get_tictoc()),
    _conflictlist(),
    _pps_ratio(coloring_handle->get_min_reduction_for_conflictlist()),
    _min_vertex_cut_off(coloring_handle->get_min_elements_for_conflictlist()),
    _edge_filtering(coloring_handle->get_vb_edge_filtering()),
    _chunkSize(coloring_handle->get_vb_chunk_size()),
    _use_color_set(),
    _max_num_iterations(coloring_handle->get_max_number_of_iterations())
    {
      switch (coloring_handle->get_coloring_type()){
      case Experimental::KokkosKernels::Graph::COLORING_VB:
        this->_use_color_set = 0;
        break;
      case Experimental::KokkosKernels::Graph::COLORING_VBBIT:
        this->_use_color_set = 2;
        break;
      case Experimental::KokkosKernels::Graph::COLORING_VBCS:
        this->_use_color_set = 1;
        break;
      default: //cannnot get in here.
        this->_use_color_set = 0;
        break;

      }

      switch (coloring_handle->get_conflict_list_type()){
      case Experimental::KokkosKernels::Graph::COLORING_NOCONFLICT:
        this->_conflictlist = 0;
        break;
      case Experimental::KokkosKernels::Graph::COLORING_ATOMIC:
        this->_conflictlist = 1;
        break;
      case Experimental::KokkosKernels::Graph::COLORING_PPS:
        this->_conflictlist = 2;
        break;
      }
    }

  /** \brief GraphColor_VB destructor.
    */
  virtual ~GraphColor_VB(){}

  /** \brief Function to color the vertices of the graphs. Performs a vertex-based coloring.
   * \param colors is the output array corresponding the color of each vertex.Size is this->nv.
   *   Attn: Color array must be nonnegative numbers. If there is no initial colors,
   *   it should be all initialized with zeros. Any positive value in the given array, will make the
   *   algorithm to assume that the color is fixed for the corresponding vertex.
   * \param num_phases: The number of iterations (phases) that algorithm takes to converge.
   */
  virtual void color_graph(color_array_type &colors,int &num_loops){
    if (this->_ticToc) {
      std::cout
          << "\tVB params:" << std::endl
          << "\tuseConflictList:" << int (this->_conflictlist) << std::endl
          << "\talgorithm:" << (int)this->_use_color_set << std::endl
          << "\tserialConflictResolution:"  << (int) this->_serialConflictResolution << std::endl
          << "\tticToc:" << (int) this->_ticToc << std::endl
          << "\tuse_color_set:" << (int) this->_use_color_set << std::endl
          << "\tpps_ratio:" << this->_pps_ratio << std::endl
          << "\tmin_vertex_cut_off:" << this->_min_vertex_cut_off << std::endl
          << "\tedge_filtering:" << (int)this->_edge_filtering << std::endl
          << "\tmax_num_iterations:" << this->_max_num_iterations << std::endl
          << "\tchunkSize:" << this->_chunkSize << std::endl;
    }

    //if the edge filtering is selected, then we do swaps on the adj array.
    //to not to touch the given one, we copy the adj array.
    idx_edge_array_type adj_copy = this->adj;

    //if we use edge-filtering, we perform swaps.
    //We need to copy the adjacency array so that we dont harm the original one.
    if (this->_edge_filtering){
      adj_copy = idx_edge_array_type("adj copy", this->ne);
      Kokkos::deep_copy(adj_copy, this->adj);
    }

    //if color set algorithm is used, we need one more array to represent the range.
    idx_temp_work_array_type vertex_color_set;
    if (this->_use_color_set == 1){
      vertex_color_set = idx_temp_work_array_type("colorset", this->nv);
    }

    //the conflictlist
    idx_temp_work_array_type current_vertexList =
        idx_temp_work_array_type(Kokkos::ViewAllocateWithoutInitializing("vertexList"), this->nv);

    //init vertexList sequentially.
    Kokkos::parallel_for(my_exec_space(0, this->nv), functorInitList (current_vertexList));


    // the next iteration's conflict list
    idx_temp_work_array_type next_iteration_recolorList;
    // the size of the current conflictlist

    //the size of the next iteration's conflictlist
    idx_type next_iteration_recolorListLength;
    //if parallel prefix sum is selected instead of atomic operations,
    //we need one more work array to do the prefix sum.
    idx_temp_work_array_type pps_work_view;

    // if a conflictlist is used
    if (this->_conflictlist > 0){
      // Vertices to recolor. Will swap with vertexList.
      next_iteration_recolorList = idx_temp_work_array_type(Kokkos::ViewAllocateWithoutInitializing("recolorList"), this->nv);
      next_iteration_recolorListLength = idx_type("recolorListLength");
      if (this->_conflictlist == 2) {
        pps_work_view = idx_temp_work_array_type("pps_view", this->nv);
      }
    }

    idx numUncolored = this->nv;
    idx current_vertexListLength = this->nv;


    double t, total=0.0;
    Kokkos::Impl::Timer timer;


    int iter=0;
    for (; (iter < this->_max_num_iterations) && (numUncolored>0); iter++){

      // First color greedy speculatively,
      //some conflicts expected
      this->colorGreedy(
          this->xadj,
          adj_copy,
          colors,
          vertex_color_set,
          current_vertexList,
          current_vertexListLength
      );

      MyExecSpace::fence();

      if (this->_ticToc){
        t = timer.seconds();
        total += t;
        std::cout << "\tTime speculative greedy phase " << iter << " : " << t << std::endl;
        timer.reset();
      }

      bool swap_work_arrays = true;

      numUncolored = this->findConflicts(
          swap_work_arrays,
          this->xadj, adj_copy,
          colors, vertex_color_set,
          current_vertexList, current_vertexListLength,
          next_iteration_recolorList, next_iteration_recolorListLength,
          pps_work_view
      );

      MyExecSpace::fence();

      if (_ticToc){
        t = timer.seconds();
        total += t;
        std::cout << "\tTime conflict detection " << iter << " : " << t << std::endl;
        timer.reset();
      }

      if (this->_serialConflictResolution) break; // Break after first iteration.
      if (this->_conflictlist && swap_work_arrays && (iter + 1< this->_max_num_iterations)){
        // Swap recolorList and vertexList
        idx_temp_work_array_type temp = current_vertexList;
        current_vertexList = next_iteration_recolorList;
        next_iteration_recolorList = temp;
        current_vertexListLength = numUncolored;
        next_iteration_recolorListLength = idx_type("recolorListLength");
      }

    }

    //if VBCS algorithm is used, the colors are converted back to original form.
    if (this->_use_color_set == 1){
      Kokkos::parallel_for(my_exec_space(0, this->nv), set_final_colors (colors, vertex_color_set));
    }
    if (numUncolored > 0){

      // Resolve conflicts by recoloring in serial
      this->resolveConflicts(
          this->nv,
          this->xadj, adj_copy,
          colors,
          current_vertexList, current_vertexListLength
      );
      MyExecSpace::fence();
      if (_ticToc){
        t = timer.seconds();
        total += t;
        std::cout << "\tTime serial conflict resolution: " << t << std::endl;
      }
    }
    num_loops = iter;
  }





private:
  /** \brief Performs speculative coloring based on the given colorings.
   *  \param xadj_: row map of the graph
   *  \param adj_: entries, columns of the graph
   *  \param vertex_colors_: colors corresponding to each vertex
   *  \param vertex_color_set: if VBCS is used, color set of each vertex
   *  \param current_vertexList_: current conflictlist
   *  \param current_vertexListLength_: size of current conflictlist
   */
  void  colorGreedy(
      idx_array_type xadj_,
      idx_edge_array_type adj_,
      color_array_type vertex_colors_,
      idx_array_type vertex_color_set,
      idx_temp_work_array_type current_vertexList_,
      idx current_vertexListLength_) {

    idx chunkSize_ = this->_chunkSize; // Process chunkSize vertices in one chunk

    if (current_vertexListLength_ < 100*chunkSize_) chunkSize_ = 1;

    //if the algorithm VBBIT
    if (this->_use_color_set == 2) {

      //If edge filtering is applied
      if (this->_edge_filtering){
        functorGreedyColor_IMPLOG_EF gc(
            xadj_, adj_,
            vertex_colors_, current_vertexList_,
            current_vertexListLength_, chunkSize_);
        Kokkos::parallel_for(my_exec_space(0, current_vertexListLength_/chunkSize_+1), gc);
      }
      else {
        functorGreedyColor_IMPLOG gc(
            xadj_, adj_,
            vertex_colors_,  current_vertexList_,
            current_vertexListLength_, chunkSize_);
        Kokkos::parallel_for(my_exec_space(0, current_vertexListLength_/chunkSize_+1), gc);
      }
    }
    // VBCS algorithm
    else if (this->_use_color_set == 1){
      if (this->_edge_filtering){
        functorGreedyColor_IMP_EF gc(
            xadj_, adj_,
            vertex_colors_, vertex_color_set, current_vertexList_,
            current_vertexListLength_, chunkSize_);
        Kokkos::parallel_for(my_exec_space(0, current_vertexListLength_/chunkSize_+1), gc);
      }
      else {
        functorGreedyColor_IMP gc(
            xadj_, adj_,
            vertex_colors_, vertex_color_set, current_vertexList_,
            current_vertexListLength_, chunkSize_);
        Kokkos::parallel_for(my_exec_space(0, current_vertexListLength_/chunkSize_+1), gc);
      }
    }
    //VB algorithm
    else if (this->_use_color_set == 0)
    {
      if (this->_edge_filtering){
        functorGreedyColor_EF  gc(
            xadj_, adj_,
            vertex_colors_,
            current_vertexList_, current_vertexListLength_, chunkSize_);
        Kokkos::parallel_for(my_exec_space(0, current_vertexListLength_/chunkSize_+1), gc);
      }
      else{
        functorGreedyColor  gc(
            xadj_, adj_,
            vertex_colors_,
            current_vertexList_, current_vertexListLength_, chunkSize_);
        Kokkos::parallel_for(my_exec_space(0, current_vertexListLength_/chunkSize_+1), gc);
      }
    }
  }

  /** \brief Performs conflict resolution
   *  \param swap_work_arrays: An output parameter whether indicating whether the work arrays should be swapped or not.
   *  \param xadj_: row map of the graph
   *  \param adj_: entries, columns of the graph
   *  \param vertex_colors_: colors corresponding to each vertex
   *  \param vertex_color_set: if VBCS is used, color set of each vertex
   *  \param current_vertexList_: current conflictlist
   *  \param current_vertexListLength_: size of current conflictlist
   *  \param next_iteration_recolorList_: current conflictlist
   *  \param next_iteration_recolorListLength_: size of current conflictlist
   *  \param pps_work_view: size of current conflictlist
   */
  idx findConflicts(
      bool &swap_work_arrays,
      idx_array_type xadj_,
      idx_edge_array_type adj_,
      color_array_type vertex_colors_, idx_array_type vertex_color_set_,
      idx_array_type current_vertexList_, idx current_vertexListLength_,
      idx_array_type next_iteration_recolorList_, idx_type next_iteration_recolorListLength_,
      idx_array_type pps_work_view) {

    swap_work_arrays = true;
    idx numUncolored = 0;
    if (this->_conflictlist == 0){
      if (this->_use_color_set == 0 || this->_use_color_set == 2){
        functorFindConflicts_No_Conflist conf( xadj_, adj_, vertex_colors_);
        Kokkos::parallel_reduce(my_exec_space(0, current_vertexListLength_), conf, numUncolored);
      }
      else {
        functorFindConflicts_No_Conflist_IMP conf(xadj_, adj_,vertex_colors_, vertex_color_set_);
        Kokkos::parallel_reduce(my_exec_space(0, current_vertexListLength_), conf, numUncolored);
      }
    }
    else if (this->_conflictlist == 2){ //IF PPS
      if (this->_use_color_set == 0 || this->_use_color_set == 2){
        // Check for conflicts. Compute numUncolored == numConflicts.
        functorFindConflicts_PPS conf(xadj_, adj_,vertex_colors_,current_vertexList_,next_iteration_recolorList_);
        Kokkos::parallel_reduce(my_exec_space(0, current_vertexListLength_), conf, numUncolored);
      }
      else {
        functorFindConflicts_PPS_IMP conf(
            xadj_, adj_,vertex_colors_, vertex_color_set_,
            current_vertexList_,next_iteration_recolorList_);
        Kokkos::parallel_reduce(my_exec_space(0, current_vertexListLength_), conf, numUncolored);
      }


      if( numUncolored && (current_vertexListLength_ >= this->_min_vertex_cut_off) &&
          (double (numUncolored) / current_vertexListLength_  <  (1.0 - this->_pps_ratio))){
        if (this->_ticToc){
          std::cout << "\tcreating work array with pps current_vertexListLength_:" <<
              current_vertexListLength_ << " params->min_vertex_cut_off:" << this->_min_vertex_cut_off << std::endl;
        }
        um_array_type h_numUncolored(&numUncolored);
        Kokkos::deep_copy (next_iteration_recolorListLength_, h_numUncolored);

        MyExecSpace::fence();

        Kokkos::parallel_scan (my_exec_space(0, current_vertexListLength_),
            parallel_prefix_sum(current_vertexList_, next_iteration_recolorList_, pps_work_view));

        MyExecSpace::fence();
        Kokkos::parallel_for (my_exec_space(0, current_vertexListLength_),
            create_new_work_array(current_vertexList_, next_iteration_recolorList_, pps_work_view));
      }
      else {
        swap_work_arrays = false;
      }
    }
    else { //IF ATOMIC
      if (this->_use_color_set == 0 || this->_use_color_set == 2){
        // Check for conflicts. Compute numUncolored == numConflicts.
        functorFindConflicts_Atomic conf(
            xadj_, adj_,vertex_colors_,current_vertexList_,
            next_iteration_recolorList_, next_iteration_recolorListLength_);
        Kokkos::parallel_reduce(my_exec_space(0, current_vertexListLength_), conf, numUncolored);
      }
      else {
        functorFindConflicts_Atomic_IMP conf(
            xadj_, adj_,vertex_colors_, vertex_color_set_,
            current_vertexList_,next_iteration_recolorList_, next_iteration_recolorListLength_);
        Kokkos::parallel_reduce(my_exec_space(0, current_vertexListLength_), conf, numUncolored);
      }
    }
    if (this->_ticToc){
      std::cout << "\tnumUncolored:" << numUncolored << std::endl;
    }
    return numUncolored;
  }

  /** \brief Performs serial conflict resolution
   *  \param nv: Number of vertices
   *  \param xadj_: row map of the graph
   *  \param adj_: entries, columns of the graph
   *  \param vertex_colors_: colors corresponding to each vertex
   *  \param current_vertexList_: current conflictlist
   *  \param current_vertexListLength_: size of current conflictlist
   */
  void  resolveConflicts(
      idx nv,
      idx_array_type xadj_,
      idx_edge_array_type adj_,
      color_array_type vertex_colors_,
      idx_array_type current_vertexList_,
      idx current_vertexListLength_) {

    color_type *forbidden = new color_type[nv];
    idx i=0;
    idx end = nv;
    host_view_type h_recolor_list;

    if (this->_conflictlist){
      end = current_vertexListLength_;
      h_recolor_list = Kokkos::create_mirror_view (current_vertexList_);
      Kokkos::deep_copy (h_recolor_list, current_vertexList_);
    }
    host_color_view_type h_colors = Kokkos::create_mirror_view (vertex_colors_);
    host_view_type h_idx = Kokkos::create_mirror_view (xadj_);
    host_edge_view_type h_adj = Kokkos::create_mirror_view (adj_);


    Kokkos::deep_copy (h_colors, vertex_colors_);
    Kokkos::deep_copy (h_idx, xadj_);
    Kokkos::deep_copy (h_adj, adj_);

    for (idx k=0; k <end; k++){
      if (this->_conflictlist){
        i = h_recolor_list(k);
      }
      else {
        // Check for uncolored vertices
        i = k;
      }
      if (h_colors(i) > 0) continue;
      for (idx j=h_idx(i); j<h_idx(i+1); j++){
        if (h_adj(j) == i) continue; // Skip self-loops
        forbidden[h_colors(h_adj(j))] = i;
      }
      // color vertex i with smallest available color
      int c=1;
      while ((forbidden[c]==i)) c++;
      h_colors(i) = c;
    }
    Kokkos::deep_copy (vertex_colors_, h_colors);
    delete [] forbidden;
  }

public:
  //Speculative Coloring Functors

  /**
   * Functor for VBBIT algorithms speculative coloring with edge filtering.
   */
  struct functorGreedyColor_IMPLOG_EF {

    idx_array_type _idx; //rowmap
    idx_edge_array_type _adj; //entries
    color_array_type _colors; // vertex colors
    idx_array_type _vertexList; // conflictlist
    idx _vertexListLength;
    idx _chunkSize;

    functorGreedyColor_IMPLOG_EF(
        idx_array_type xadj,
        idx_edge_array_type adj,
        color_array_type colors,
        idx_array_type vertexList,
        idx vertexListLength,
        idx chunkSize) :_idx(xadj), _adj(adj), _colors(colors),
      _vertexList(vertexList), _vertexListLength(vertexListLength), _chunkSize(chunkSize){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii) const {
      idx i = 0;

      //outer loop is on chunks, a thread is assigned as many vertex as the chunksize
      for (idx ichunk=0; ichunk<_chunkSize; ichunk++){
        if (ii*_chunkSize +ichunk < _vertexListLength){
          i = _vertexList(ii*_chunkSize +ichunk);
        }
        else{
          continue;
        }

        if (_colors(i) > 0) continue; // Already colored this vertex

        idx my_xadj_end = _idx(i+1);
        idx xadjbegin = _idx(i);

        // Do multiple passes if array is too small.
        idx degree = my_xadj_end-xadjbegin; // My degree
        idx offset = 0;

        //we parse the neigborlist multiple times,
        //each time we look for a certain range of colors.
        for (; (offset <= degree + VBBIT_COLORING_FORBIDDEN_SIZE); offset += VBBIT_COLORING_FORBIDDEN_SIZE){

          // Forbidden colors
          //we use a single (long) int for forbidden colors
          ban_type forbidden = 0;

          // Check nbors, fill forbidden array.
          for (idx j=xadjbegin; j<my_xadj_end; ++j){
            idx n  = _adj(j);
            if (n == i) continue; // Skip self-loops
            color_type c = _colors(n);

            color_type color_offset = c-offset;
            //if color is within the current range, or if its color is in a previously traversed range
            if (c && color_offset <= VBBIT_COLORING_FORBIDDEN_SIZE){
              //apply edge filtering, place it to front of the adjacency list,
              //so that we wont see that anymore.
              if (j > xadjbegin){
                _adj(j) = _adj(xadjbegin);
                _adj(xadjbegin) = n;
              }
              ++xadjbegin;

              //if it is in the current range, then add the color to banned colors
              if (c > offset){
                //convert color to bit representation.
                ban_type ban_color_bit = 1;
                ban_color_bit = ban_color_bit << (color_offset - 1);
                //add it to forbidden colors
                forbidden = forbidden | ban_color_bit;
                //if there are no available colors in this range,
                //early exit, no need to traverse the rest.
                if (~forbidden == 0) {
                  break;
                }
              }
            }
          }

          forbidden = ~(forbidden);
          //check if an available color exits.
          if (forbidden){
            //if there is an available color, choose the first color,
            //using 2s complement.
            ban_type my_new_color = forbidden & (-forbidden);
            color_type val = 1;
            //convert it back to decimal color.
            while ((my_new_color & 1) == 0) {
              ++val;
              my_new_color = my_new_color >> 1;
            }
            _colors(i) = val + offset;
            break;
          }
        }
      }
    }
  };

  /**
   * Functor for VBBIT algorithms speculative coloring without edge filtering.
   */
  struct functorGreedyColor_IMPLOG {
    idx_array_type _idx;
    idx_edge_array_type _adj;
    color_array_type _colors;
    idx_array_type _vertexList;
    idx _vertexListLength;
    idx _chunkSize;

    functorGreedyColor_IMPLOG(
        idx_array_type xadj,
        idx_edge_array_type adj,
        color_array_type colors,
        idx_array_type vertexList,
        idx vertexListLength,
        idx chunkSize) :
          _idx(xadj), _adj(adj), _colors(colors),
          _vertexList(vertexList), _vertexListLength(vertexListLength), _chunkSize(chunkSize){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii) const {
      idx i = 0;
      for (idx ichunk=0; ichunk<_chunkSize; ichunk++){
        if (ii*_chunkSize +ichunk < _vertexListLength)
          i = _vertexList(ii*_chunkSize +ichunk);
        else
          continue;

        if (_colors(i) > 0) continue; // Already colored this vertex

        idx my_xadj_end = _idx(i+1);
        idx xadjbegin = _idx(i);
        // Do multiple passes if array is too small.
        idx degree = my_xadj_end - xadjbegin; // My degree
        idx offset = 0;

        for (; (offset <= degree + VBBIT_COLORING_FORBIDDEN_SIZE); offset += VBBIT_COLORING_FORBIDDEN_SIZE){

          ban_type forbidden = 0; // Forbidden colors

          // Check nbors, fill forbidden array.
          for (idx j=xadjbegin; j<my_xadj_end; ++j){
            idx n  = _adj(j);
            if (n == i) continue; // Skip self-loops
            color_type c = _colors(n);
            color_type color_offset = c-offset;
            //if color is in the current range
            //convert it to binary and add it to forbidden
            if (color_offset <= VBBIT_COLORING_FORBIDDEN_SIZE && c > offset){

              ban_type ban_color_bit = 1;
              ban_color_bit = ban_color_bit << (color_offset - 1);

              forbidden = forbidden | ban_color_bit;
              if (~forbidden == 0) {
                break;
              }
            }
          }
          forbidden = (~forbidden);

          if (forbidden){
            ban_type my_new_color = forbidden & (-forbidden);

            color_type val = 1;

            while ((my_new_color & 1) == 0) {
              ++val;
              my_new_color = my_new_color >> 1;
            }
            _colors(i) = val + offset;
            break;
          }
        }
      }
    }

  };

  /**
   * Functor for VBCS algorithms speculative coloring with edge filtering.
   */
  struct functorGreedyColor_IMP_EF {
    idx_array_type _xadj;
    idx_edge_array_type _adj;
    color_array_type _colors;
    idx_array_type _color_set ;
    idx_array_type _vertexList;
    idx _vertexListLength;
    idx _chunkSize;

    functorGreedyColor_IMP_EF(
        idx_array_type xadj,
        idx_edge_array_type adj,
        color_array_type colors, idx_array_type color_set,
        idx_array_type vertexList,
        idx vertexListLength,
        idx chunkSize):
          _xadj(xadj), _adj(adj),
          _colors(colors), _color_set(color_set),
          _vertexList(vertexList), _vertexListLength(vertexListLength),
          _chunkSize(chunkSize){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &ii) const {
      idx i = 0;
      for (idx ichunk=0; ichunk<_chunkSize; ichunk++){
        if (ii*_chunkSize +ichunk < _vertexListLength)
          i = _vertexList(ii*_chunkSize +ichunk);
        else
          continue;

        if (_colors(i) > 0) continue; // Already colored this vertex
        idx xadj_end = _xadj(i+1);
        idx xadj_begin = _xadj(i);

        //my color set starts from zero, but if we are leaving vertices
        //that cannot be colored in this iteration, we retrieve it from their previous color_sets.
        idx my_color_set = 0;
        while (1){
          color_type ban_colors = 0;

          for (idx j = xadj_begin; j < xadj_end && ~ban_colors; ++j){
            idx n = _adj(j);
            if (n == i) continue; // Skip self-loops

            idx neighbor_color_set = _color_set(n);
            //only if the neigbor has the same color set
            if (neighbor_color_set <= my_color_set ){
              color_type ncolor = _colors(n);
              if (ncolor){
                if(j > xadj_begin){
                  _adj(j) = _adj(xadj_begin);
                  _adj(xadj_begin) = n;
                }
                ++xadj_begin;
                if (neighbor_color_set == my_color_set){
                  ban_colors = ban_colors | ncolor;
                }
              }
            }
          }

          ban_colors = ~(ban_colors);
          if (ban_colors){
            color_type my_color = ban_colors & (-ban_colors);
            _color_set(i) = my_color_set;
            _colors(i) = my_color;
            break;
          }
          else{
            my_color_set += 1;
          }
        }
      }
    }
  };

  /**
   * Functor for VBCS algorithms speculative coloring without edge filtering.
   */
  struct functorGreedyColor_IMP {
    idx_array_type _xadj;
    idx_edge_array_type _adj;
    color_array_type _colors;
    idx_array_type _color_set ;
    idx_array_type _vertexList;
    idx _vertexListLength;
    idx _chunkSize;


    functorGreedyColor_IMP(
        const idx_array_type xadj,
        const idx_edge_array_type adj,
        color_array_type colors, idx_array_type color_set,
        idx_array_type vertexList,
        idx vertexListLength,
        idx chunkSize
    ) :
      _xadj(xadj), _adj(adj),
      _colors(colors), _color_set(color_set),
      _vertexList(vertexList), _vertexListLength(vertexListLength),
      _chunkSize(chunkSize){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &ii) const {
      idx i = 0;
      for (idx ichunk=0; ichunk<_chunkSize; ichunk++){
        if (ii*_chunkSize +ichunk < _vertexListLength)
          i = _vertexList(ii*_chunkSize +ichunk);
        else
          continue;

        if (_colors(i) > 0) continue; // Already colored this vertex
        idx xadj_end = _xadj(i+1);
        idx xadj_begin = _xadj(i);

        //my color set starts from zero, but if we are leaving vertices
        //that cannot be colored in this iteration, we retrieve it from their previous color_sets.
        idx my_color_set =  0;
        //idx degree = xadj_end - xadj_begin;
        for (; ;){
          color_type ban_colors = 0;
          for (idx j = xadj_begin; j < xadj_end ;++j){
            idx n = _adj(j);
            if (n == i) continue; // Skip self-loops
            if ( my_color_set == _color_set(n)){
              ban_colors = ban_colors | _colors(n);
              if (~ban_colors == 0){
                break;
              }
            }
          }

          ban_colors = ~(ban_colors);

          if (ban_colors){
            color_type my_color = ban_colors & (-ban_colors);
            _color_set(i) = my_color_set;
            _colors(i) = my_color;
            break;
          }
          else{
            my_color_set += 1;
          }
        }
      }
    }
  };

  /**
   * Functor for VB algorithm speculative coloring with edge filtering.
   */
  struct functorGreedyColor_EF {
    idx_array_type _idx;
    idx_edge_array_type _adj;
    color_array_type _colors;
    idx_array_type _vertexList;
    idx _vertexListLength;
    idx _chunkSize;

    functorGreedyColor_EF(
        const idx_array_type xadj,
        const idx_edge_array_type adj,
        color_array_type colors,
        idx_array_type vertexList,
        idx vertexListLength,
        idx chunkSize
    ) :
      _idx(xadj), _adj(adj), _colors(colors),
      _vertexList(vertexList), _vertexListLength(vertexListLength), _chunkSize(chunkSize){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii) const {
      // Color vertex i with smallest available color.
      //
      // Each thread colors a chunk of vertices to prevent all
      // vertices getting the same color.
      //
      // This version uses a bool array of size FORBIDDEN_SIZE.
      // TODO: With chunks, the forbidden array should be char/int
      //       and reused for all vertices in the chunk.
      //
      idx i = 0;
      for (idx ichunk=0; ichunk<_chunkSize; ichunk++){
        if (ii*_chunkSize +ichunk < _vertexListLength)
          i = _vertexList(ii*_chunkSize +ichunk);
        else
          continue;

        if (_colors(i) > 0) continue; // Already colored this vertex

        bool foundColor = false; // Have we found a valid color?

        // Use forbidden array to find available color.
        // This array should be small enough to fit in fast memory (use Kokkos memoryspace?)
        bool forbidden[VB_COLORING_FORBIDDEN_SIZE]; // Forbidden colors

        // Do multiple passes if array is too small.
        idx degree = _idx(i+1)-_idx(i); // My degree
        idx my_xadj_end = _idx(i+1);
        idx offset = 0;
        idx xadjbegin = _idx(i);

        for (; (offset <= degree + VB_COLORING_FORBIDDEN_SIZE) && (!foundColor); offset += VB_COLORING_FORBIDDEN_SIZE){
          // initialize
          for (int j=0; j< VB_COLORING_FORBIDDEN_SIZE; j++){
            forbidden[j] = false;
          }
          if (offset == 0) forbidden[0] = true; // by convention, start at 1

          // Check nbors, fill forbidden array.
          for (idx j=xadjbegin; j<my_xadj_end; ++j){
            idx n  = _adj(j);
            if (n == i) {
              continue; // Skip self-loops
            }
            color_type c= _colors(n);
            // Removed option to leave potentially conflicted vertices uncolored.
            //if (c== -1){ // Nbor is being colored at same time
            //  _colors[i] = 0; // Neutral color, skip and recolor later
            //  foundColor = true;
            //  return;
            //}
            if ((c>= offset) && (c-offset < VB_COLORING_FORBIDDEN_SIZE)){
              forbidden[c-offset] = true;
            }
            if (c && c-offset < VB_COLORING_FORBIDDEN_SIZE){
              if (j > xadjbegin){
                _adj(j) = _adj(xadjbegin);
                _adj(xadjbegin) = n;
              }
              ++xadjbegin;
            }

          }

          // color vertex i with smallest available color (FirstFit)
          // TODO: Add options for other color choices (Random, LeastUsed)
          for (int c=0; c< VB_COLORING_FORBIDDEN_SIZE; c++){
            if (!forbidden[c]){
              _colors(i) = offset+c;
              //_colors[i] += (i&1); // RandX strategy to reduce conflicts
              foundColor = true;
              break;
            }
          }
        }
      }
    }

  };


  /**
   * Functor for VB algorithm speculative coloring without edge filtering.
   */
  struct functorGreedyColor {
    idx_array_type _idx;
    idx_edge_array_type _adj;
    color_array_type _colors;
    idx_array_type _vertexList;
    idx _vertexListLength;
    idx _chunkSize;

    functorGreedyColor(
        const idx_array_type xadj,
        const idx_edge_array_type adj,
        color_array_type colors,
        idx_array_type vertexList,
        idx vertexListLength,
        idx chunkSize
    ) :
      _idx(xadj), _adj(adj), _colors(colors),
      _vertexList(vertexList), _vertexListLength(vertexListLength), _chunkSize(chunkSize){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii) const {
      // Color vertex i with smallest available color.
      //
      // Each thread colors a chunk of vertices to prevent all
      // vertices getting the same color.
      //
      // This version uses a bool array of size FORBIDDEN_SIZE.
      // TODO: With chunks, the forbidden array should be char/int
      //       and reused for all vertices in the chunk.
      //
      idx i = 0;
      for (idx ichunk=0; ichunk<_chunkSize; ichunk++){
        if (ii*_chunkSize +ichunk < _vertexListLength)
          i = _vertexList(ii*_chunkSize +ichunk);
        else
          continue;

        if (_colors(i) > 0) continue; // Already colored this vertex

        bool foundColor = false; // Have we found a valid color?

        // Use forbidden array to find available color.
        // This array should be small enough to fit in fast memory (use Kokkos memoryspace?)
            bool forbidden[VB_COLORING_FORBIDDEN_SIZE]; // Forbidden colors

        // Do multiple passes if array is too small.
        idx degree = _idx(i+1)-_idx(i); // My degree
        idx offset = 0;
        for (; (offset <= degree + VB_COLORING_FORBIDDEN_SIZE) && (!foundColor); offset += VB_COLORING_FORBIDDEN_SIZE){
          // initialize
          for (int j=0; j< VB_COLORING_FORBIDDEN_SIZE; j++){
            forbidden[j] = false;
          }
          if (offset == 0) forbidden[0] = true; // by convention, start at 1

          // Check nbors, fill forbidden array.
          for (idx j=_idx(i); j<_idx(i+1); j++){
            if (_adj(j) == i) continue; // Skip self-loops
            color_type c= _colors(_adj(j));
            // Removed option to leave potentially conflicted vertices uncolored.
            //if (c== -1){ // Nbor is being colored at same time
            //  _colors[i] = 0; // Neutral color, skip and recolor later
            //  foundColor = true;
            //  return;
            //}
            if ((c>= offset) && (c-offset < VB_COLORING_FORBIDDEN_SIZE))
              forbidden[c-offset] = true;
          }

          // color vertex i with smallest available color (FirstFit)
          // TODO: Add options for other color choices (Random, LeastUsed)
          for (int c=0; c< VB_COLORING_FORBIDDEN_SIZE; c++){
            if (!forbidden[c]){
              _colors(i) = offset+c;
              //_colors[i] += (i&1); // RandX strategy to reduce conflicts
              foundColor = true;
              break;
            }
          }
        }
      }
    }

  };


  //Conflict find and worklist creation functors.

  /**
   * Finds conflicts without creating a new worklist
   */
  struct functorFindConflicts_No_Conflist {
    idx_array_type _idx;
    idx_edge_array_type _adj;
    color_array_type _colors;

    functorFindConflicts_No_Conflist(
        const idx_array_type xadj,
        const idx_edge_array_type adj,
        color_array_type colors) :
          _idx(xadj), _adj(adj),_colors(colors)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii, idx &numConflicts) const {

      color_type my_color = _colors(ii);
      idx xadjend = _idx(ii+1);
      idx j=_idx(ii);
#ifdef DEGREECOMP
      idx myDegree = xadjend - j;
#endif
      for (; j< xadjend; j++){
        idx neighbor = _adj(j);

        if (
#ifndef DEGREECOMP
            ii < neighbor &&
#endif
            _colors(neighbor) == my_color
#ifdef DEGREECOMP
            && (myDegree < _idx(neighbor + 1) - _idx(neighbor) ||
                (myDegree == _idx(neighbor + 1) - _idx(neighbor) && ii < neighbor))
#endif
        ) {
          //std::cout << "me:" << ii << " n:" << neighbor << " color:" << my_color << std::endl;
          _colors(ii) = 0; // Uncolor vertex i
          numConflicts += 1;
          break; // Once i is uncolored and marked conflict
        }
      }
    }
  };


  /**
   * Finds conflicts by marking the work vertices to be used later for creation of new worklist with PPS
   */
  struct functorFindConflicts_PPS {
    idx_array_type _idx;
    idx_edge_array_type _adj;
    color_array_type _colors;
    idx_array_type _vertexList;
    idx_array_type _recolorList;



    functorFindConflicts_PPS(
        const idx_array_type xadj,
        const idx_edge_array_type adj,
        color_array_type colors,
        idx_array_type vertexList,
        idx_array_type recolorList) :
          _idx(xadj), _adj(adj), _colors(colors),
          _vertexList(vertexList),
          _recolorList(recolorList){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii, idx &numConflicts) const {
      idx i = _vertexList(ii);
      color_type my_color = _colors(i);
      _recolorList(i) = 0;
      // check vertex i conflicts

      idx xadjend = _idx(i+1);
      idx j=_idx(i);
#ifdef DEGREECOMP
      idx myDegree = xadjend - j;
#endif
      for (; j<xadjend; j++){
        idx neighbor = _adj(j);
        if (
#ifndef DEGREECOMP
            i < neighbor &&
#endif
            _colors(neighbor) == my_color
#ifdef DEGREECOMP
            && (myDegree < _idx(neighbor + 1) - _idx(neighbor)
                                                     || (myDegree == _idx(neighbor + 1) - _idx(neighbor) && i < neighbor))
#endif
        ) {
          _colors(i) = 0; // Uncolor vertex i
          _recolorList(i) = 1;
          numConflicts += 1;
          break; // Once i is uncolored and marked conflict
        }
      }
    }
  };


  /**
   * Finds conflicts and creates new worklist using atomic operations.
   */
  struct functorFindConflicts_Atomic {
    idx_array_type _idx;
    idx_edge_array_type _adj;
    color_array_type _colors;
    idx_array_type _vertexList;
    idx_array_type _recolorList;
    idx_type _recolorListLength;


    functorFindConflicts_Atomic(
        const idx_array_type xadj,
        const idx_edge_array_type adj,
        color_array_type colors,
        idx_array_type vertexList,
        idx_array_type recolorList,
        idx_type recolorListLength
    ) :
      _idx(xadj), _adj(adj), _colors(colors),
      _vertexList(vertexList),
      _recolorList(recolorList),
      _recolorListLength(recolorListLength){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii, idx &numConflicts) const {

      idx i = _vertexList(ii);
      color_type my_color = _colors(i);

      idx xadjend = _idx(i+1);
      idx j=_idx(i);
#ifdef DEGREECOMP
      idx myDegree = xadjend - j;
#endif

      for (; j < xadjend; j++){
        idx neighbor = _adj(j);
        if (
#ifndef DEGREECOMP
            i < neighbor &&
#endif
            _colors(neighbor) == my_color
#ifdef DEGREECOMP
            && (myDegree < _idx(neighbor + 1) - _idx(neighbor) ||
                (myDegree == _idx(neighbor + 1) - _idx(neighbor) && i < neighbor))
#endif
        ) {
          _colors(i) = 0; // Uncolor vertex i
          // Atomically add vertex i to recolorList
          const idx k = Kokkos::atomic_fetch_add( &_recolorListLength(), 1);
          _recolorList(k) = i;
          numConflicts += 1;
          break; // Once i is uncolored and marked conflict
        }
      }
    }
  };


  /**
   * VBCS:  Finds conflicts without creating a new worklist
   */
  struct functorFindConflicts_No_Conflist_IMP {

    idx_array_type _xadj;
    idx_edge_array_type _adj;
    color_array_type _colors;
    idx_array_type _color_sets;


    functorFindConflicts_No_Conflist_IMP(
        const idx_array_type xadj,
        const idx_edge_array_type adj,
        color_array_type colors,
        idx_array_type color_sets
    ) :
      _xadj(xadj), _adj(adj), _colors(colors), _color_sets(color_sets){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii, idx &numConflicts) const {
      color_type my_color = _colors(ii);
      if (my_color == 0){
        // this should only happen when one_color_set_per_iteration is set to true.
        numConflicts++;
      }
      else {
        idx my_color_set = _color_sets(ii);
        idx my_xadj_end = _xadj(ii+1);
        // check vertex i conflicts

        idx j=_xadj(ii);
#ifdef DEGREECOMP
        idx myDegree = my_xadj_end - j;
#endif

        for (; j<my_xadj_end; j++){
          idx neighbor = _adj(j);
          if (
#ifndef DEGREECOMP
              ii < neighbor &&
#endif
              _colors(neighbor) == my_color && my_color_set == _color_sets(neighbor)
#ifdef DEGREECOMP
              && (myDegree < _xadj(neighbor + 1) - _xadj(neighbor)||
                  (myDegree == _xadj(neighbor + 1) - _xadj(neighbor) && ii < neighbor))
#endif
          ) {
            _colors(ii) = 0; // Uncolor vertex i
            _color_sets(ii) = 0;
            numConflicts++;
            break; // Once i is uncolored and marked conflict
          }
        }

      }
    }
  };


  /**
   * VBCS: Finds conflicts by marking the work vertices to be used later for creation of new worklist with PPS
   */
  struct functorFindConflicts_PPS_IMP {
    idx_array_type _xadj;
    idx_edge_array_type _adj;
    color_array_type _colors;
    idx_array_type _color_sets;
    idx_array_type _vertexList;
    idx_array_type _recolorList;

    functorFindConflicts_PPS_IMP(
        const idx_array_type xadj,
        const idx_edge_array_type adj,
        color_array_type colors,
        idx_array_type color_sets,
        idx_array_type vertexList,
        idx_array_type recolorList
    ) :
      _xadj(xadj), _adj(adj), _colors(colors), _color_sets(color_sets),
      _vertexList(vertexList),
      _recolorList(recolorList){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii, idx &numConflicts) const {
      idx i = _vertexList(ii);
      _recolorList(i) = 0;
      color_type my_color = _colors(i);
      if (my_color == 0){
        _recolorList(i) = 1;
        numConflicts++;
      }
      else {
        idx my_color_set = _color_sets(i);
        idx my_xadj_end = _xadj(i+1);
        // check vertex i conflicts

        idx j=_xadj(i);
#ifdef DEGREECOMP
        idx myDegree = my_xadj_end - j;
#endif
        for (; j<my_xadj_end; j++){
          idx neighbor = _adj(j);
          if (
#ifndef DEGREECOMP
              i < neighbor &&
#endif
              _colors(neighbor) == my_color && my_color_set == _color_sets(neighbor)
#ifdef DEGREECOMP
              && (myDegree < _xadj(neighbor + 1) - _xadj(neighbor)||
                  (myDegree == _xadj(neighbor + 1) - _xadj(neighbor) && i < neighbor))
#endif
          ) {
            _colors(i) = 0; // Uncolor vertex i
            _color_sets(i) = 0;
            _recolorList(i) = 1;
            numConflicts++;
            break; // Once i is uncolored and marked conflict
          }
        }
      }
    }
  };


  /**
   * VBCS:Finds conflicts and creates new worklist using atomic operations.
   */
  struct functorFindConflicts_Atomic_IMP {

    idx_array_type _xadj;
    idx_edge_array_type _adj;
    color_array_type _colors;
    idx_array_type _color_sets;
    idx_array_type _vertexList;
    idx_array_type _recolorList;
    idx_type _recolorListLength;


    functorFindConflicts_Atomic_IMP(
        const idx_array_type xadj,
        const idx_edge_array_type adj,
        color_array_type colors,
        idx_array_type color_sets,
        idx_array_type vertexList,
        idx_array_type recolorList,
        idx_type recolorListLength
    ) :
      _xadj(xadj), _adj(adj), _colors(colors), _color_sets(color_sets),
      _vertexList(vertexList),
      _recolorList(recolorList),
      _recolorListLength(recolorListLength){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii, idx &numConflicts) const {
      idx i = _vertexList(ii);
      color_type my_color = _colors(i);
      if (my_color == 0){
        // this should only happen when one_color_set_per_iteration is set to true.
        const idx k = Kokkos::atomic_fetch_add( &_recolorListLength(), 1);
        _recolorList(k) = i;
        numConflicts++;
      }
      else {
        idx my_color_set = _color_sets(i);
        idx my_xadj_end = _xadj(i+1);
        // check vertex i conflicts

        idx j=_xadj(i);
#ifdef DEGREECOMP
        idx myDegree = my_xadj_end - j;
#endif
        for (; j< my_xadj_end; j++){
          idx neighbor = _adj(j);
          if (
#ifndef DEGREECOMP
              i < neighbor &&
#endif
              _colors(neighbor) == my_color && my_color_set == _color_sets(neighbor)
#ifdef DEGREECOMP
              && (myDegree < _xadj(neighbor + 1) - _xadj(neighbor)||
                  (myDegree == _xadj(neighbor + 1) - _xadj(neighbor) && i < neighbor))
#endif
          ) {
            _colors(i) = 0; // Uncolor vertex i
            _color_sets(i) = 0;
            // Atomically add vertex i to recolorList
            const idx k = Kokkos::atomic_fetch_add( &_recolorListLength(), 1);
            _recolorList(k) = i;
            numConflicts++;
            break; // Once i is uncolored and marked conflict
          }
        }
      }
    }
  };

  //Helper Functors
  /**
   * Functor to init a list sequentialy, that is list[i] = i
   */
  struct functorInitList{
    idx_array_type _vertexList;
    functorInitList (idx_array_type vertexList) : _vertexList(vertexList) { }
    KOKKOS_INLINE_FUNCTION
    void operator()(const idx i) const {
      // Natural order
      _vertexList(i) = i;
    }
  };

  /**
   * Functor for parallel prefix sum
   */
  struct parallel_prefix_sum{
    idx_array_type _vertexList;
    idx_array_type _recolorList;
    idx_array_type _pps_view;

    parallel_prefix_sum(
        idx_array_type vertexList,
        idx_array_type recolorList,
        idx_array_type pps_view):
          _vertexList(vertexList),_recolorList(recolorList),_pps_view(pps_view){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii, size_t& update, const bool final) const {
      idx w = _vertexList(ii);
      update += _recolorList(w);
      if (final) {
        _pps_view(w) = (update);
      }
    }
  };

  /**
   * Functor for creating new worklist using pps
   */
  struct create_new_work_array{
    idx_array_type _vertexList;
    idx_array_type _recolorList;
    idx_array_type _pps_view;

    create_new_work_array(
        idx_array_type vertexList,
        idx_array_type recolorList,
        idx_array_type pps_view):
          _vertexList(vertexList),_recolorList(recolorList),_pps_view(pps_view){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii) const {
      idx w = _vertexList(ii);
      idx left_work = 0;
      if (ii > 0){
        left_work = _pps_view(_vertexList(ii - 1));
      }
      idx pps_current = _pps_view(w);
      if(pps_current != left_work){
        idx future_index = pps_current;
        _recolorList(future_index - 1) = w;
      }
    }
  };


  /**
   * Converting VBCS colors to final colors.
   */
  struct set_final_colors{
    color_array_type kokcol;
    idx_array_type kokcolset; //the colors that are represented with bits, and the colors set that the color is in.
    color_type color_size;

    /** \brief functor constructor.
     * \param kokcol_  the colors of the vertices. Represented with bits.
     * \param kokcolset_  the color set of the vertices. kokcolors_ and color_set_ together
     *      is used to represent the colors e.g. color_set_(v) * (numbits_in_idx-1) + set_bit_position_in_kokcolors_(v)
     */
    set_final_colors(color_array_type kokcol_, idx_array_type kokcolset_
    ): kokcol(kokcol_),kokcolset(kokcolset_), color_size ( sizeof(color_type) * 8){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &ii) const {

      color_type val = kokcol(ii);
      if (val){
        //find the position in the bit.
        idx i = 1;
        while ((val & 1) == 0) {
          ++i;
          val = val >> 1;
        }

        //idx i = log2(val) + 1;
        //set the final color.
        kokcol(ii) = i + kokcolset(ii) * color_size;
      }
    }
  };


};


/*! \brief Class for modular parallel graph coloring using Kokkos.
 *  Performs a edge_base coloring, with the hope of better load balance
 *  as well as better memory accesses on GPUs.
 */
template <typename HandleType>
class GraphColor_EB:public GraphColor <HandleType>{
public:

  typedef typename HandleType::idx_array_type idx_array_type;
  typedef typename HandleType::idx_edge_array_type idx_edge_array_type;
  typedef typename HandleType::color_array_type color_array_type;

  typedef long long int ban_type;
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




  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;

  typedef typename Kokkos::View<idx, idx_device_type> idx_type;
  typedef typename Kokkos::View<idx, Kokkos::MemoryUnmanaged> um_array_type;

  typedef typename idx_type::HostMirror idx_host_type; //Host view type
  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;

  typedef typename Kokkos::View<idx *, idx_array_layout, MyTempMemorySpace> idx_temp_work_array_type;

  typedef typename Kokkos::View<idx *, MyPersistentMemorySpace> idx_persistent_work_array_type;


  typedef typename Kokkos::View<color_type *, color_type_array_layout, MyTempMemorySpace> color_temp_work_array_type;




  typedef Kokkos::View<char *, MyTempMemorySpace> char_array_type;
  typedef typename char_array_type::HostMirror host_char_view_type; //Host view type


public:

  /**
   * \brief GraphColor_EB constructor.
   * \param nv_ number of vertices in the graph
   * \param ne_ number of edges in the graph
   * \param xadj_ the xadj array of the graph. Its size is nv_ +1
   * \param adj_ adjacency array of the graph. Its size is ne_
   */
  GraphColor_EB(idx nv_, idx ne_, idx_array_type row_map, idx_edge_array_type entries, HandleType *coloring_handle):
    GraphColor<HandleType>(nv_, ne_, row_map, entries, coloring_handle)
    {}

  /**
   * \brief Class Destructor.
   */
  virtual ~GraphColor_EB(){}



  /** \brief function to color the vertices of the graphs. Performs an edge based graph coloring.
   *  the algorithm uses kokkos, so it is modular.
   * \param colors is the output array corresponding the color of each vertex.Size is this->nv.
   * \param num_loops is the output for the number of phases that the algorithm took to converge.
   */
  virtual void color_graph(color_array_type &kok_colors, int &num_loops ){

    //get EB parameters
    color_type numInitialColors = this->cp->get_eb_num_initial_colors();
    double pps_cutoff = this->cp->get_min_reduction_for_conflictlist();
    idx ps_min = this->cp->get_min_elements_for_conflictlist();
    bool use_pps = (this->cp->get_conflict_list_type() == Experimental::KokkosKernels::Graph::COLORING_PPS);



    bool tictoc = this->cp->get_tictoc();
    Kokkos::Impl::Timer *timer;

    if (tictoc){
      timer = new Kokkos::Impl::Timer();
      std::cout << "\tRewriting EB params. num_initial_colors:" << numInitialColors
          << " prefix_sum_shrink_min:"  << ps_min
          << " ps_cutoff:" << pps_cutoff
          << std::endl;
    }

    idx numEdges = 0;
    idx_persistent_work_array_type kok_src, kok_dst;
    this->cp->get_lower_diagonal_edge_list
        (this->nv, this->ne, this->xadj, this->adj, numEdges, kok_src, kok_dst);
    idx num_work_edges = numEdges;

    //allocate memory for vertex ban colors, and tentative bans
    color_temp_work_array_type color_ban (Kokkos::ViewAllocateWithoutInitializing("color_ban"), this->nv);
    color_temp_work_array_type tentative_color_ban(Kokkos::ViewAllocateWithoutInitializing("tentative_color_ban"), this->nv);//views are initialized with zero
    //allocate memory for vertex color set shifts.
    idx_temp_work_array_type color_set ("color_set", this->nv); //initialized with zero.
    //initialize colors, color bans
    Kokkos::parallel_for (my_exec_space (0, this->nv) , init_colors (kok_colors, color_ban, numInitialColors));

    //worklist
    idx_temp_work_array_type edge_conflict_indices
    (Kokkos::ViewAllocateWithoutInitializing("edge_conflict_indices"), num_work_edges);
    //next iterations conflict list
    idx_temp_work_array_type new_edge_conflict_indices
    (Kokkos::ViewAllocateWithoutInitializing("new_edge_conflict_indices"), num_work_edges);

    idx_temp_work_array_type
    pps(Kokkos::ViewAllocateWithoutInitializing("prefix_sum"), num_work_edges);

    char_array_type edge_conflict_marker
    (Kokkos::ViewAllocateWithoutInitializing("edge_conflict_marker"), num_work_edges);


    //initialize the worklist sequentiall, and markers as 1.
    Kokkos::parallel_for (
        my_exec_space (0, num_work_edges),
        init_work_arrays(edge_conflict_indices, edge_conflict_marker)
    );
    MyExecSpace::fence();

    double inittime = 0;
    if (tictoc){
      inittime = timer->seconds();
      timer->reset();
    }
    double mc_time = 0, cnt_time = 0, ban_time = 0, expand_ban_time = 0, color_time = 0, pps_time = 0;

    idx i = 0;
    while(1){
      ++i;
      //conflict detection mark conflicts as color 0.
      //update their bans
      Kokkos::parallel_for(
          my_exec_space(0,num_work_edges),
          halfedge_mark_conflicts (
              kok_src, kok_dst,
              kok_colors, color_set,
              color_ban, tentative_color_ban
              ,edge_conflict_indices
          )
      );
      MyExecSpace::fence();
      if (tictoc){
        mc_time += timer->seconds();
        timer->reset();
      }

      idx num_conflict_reduction = 0;
      //count conflicts, and mark the edges that does not need to be processed.
      Kokkos::parallel_reduce(
          my_exec_space(0, num_work_edges),
          halfedge_conflict_count(
              kok_src, kok_dst,
              kok_colors, color_set,
              edge_conflict_indices, edge_conflict_marker
          ),num_conflict_reduction);

      MyExecSpace::fence();

      if (tictoc){
        cnt_time += timer->seconds();
        timer->reset();
      }


      if (num_work_edges - num_conflict_reduction == 0) break;

      //if the reduction is good enough w.r.t. parameters, create new worklist.
      if (num_work_edges > ps_min && num_conflict_reduction / double (num_work_edges) > pps_cutoff)
      {
        //use_pps = false;
        if (use_pps){
          //calculate new positions of the edges in new worklist
          Kokkos::parallel_scan (
              my_exec_space(0, num_work_edges),
              parallel_prefix_sum(edge_conflict_indices, edge_conflict_marker, pps)
          );
          MyExecSpace::fence();

          //write the edge indices to new worklist.
          Kokkos::parallel_for (
              my_exec_space(0, num_work_edges),
              create_new_work_array(edge_conflict_indices, edge_conflict_marker, pps, new_edge_conflict_indices));
          MyExecSpace::fence();
        }
        else {
          //create new worklist
          idx_type new_index = idx_type("recolorListLength");;
          Kokkos::parallel_for (
              my_exec_space(0, num_work_edges),
              atomic_create_new_work_array(new_index, edge_conflict_indices, edge_conflict_marker, new_edge_conflict_indices));
          MyExecSpace::fence();
        }

        //swap old and new worklist
        idx_temp_work_array_type tmp = new_edge_conflict_indices;
        new_edge_conflict_indices =  edge_conflict_indices;
        edge_conflict_indices = tmp;
        num_work_edges -= num_conflict_reduction;
        num_conflict_reduction = 0;

        if (tictoc){
          pps_time += timer->seconds();
          timer->reset();
        }
      }

      //create ban colors using the colored neighbors
      Kokkos::parallel_for (
          my_exec_space(0,num_work_edges),
          halfedge_ban_colors(
              kok_src, kok_dst,
              kok_colors, color_set,
              color_ban,edge_conflict_indices, edge_conflict_marker
          )
      );

      MyExecSpace::fence();

      if (tictoc){
        ban_time += timer->seconds();
        timer->reset();
      }


      //create tentative ban using the uncolored neighbors.
      Kokkos::parallel_for (
          my_exec_space(0,num_work_edges),
          halfedge_expand_ban_for_unmatched_neighbors(
              kok_src, kok_dst,
              kok_colors, color_set,
              color_ban,
              tentative_color_ban,edge_conflict_indices)
      );

      if (tictoc){
        expand_ban_time += timer->seconds();
        timer->reset();
      }



      //chose a color based on the ban arrays.
      //if all colors in current set are taken, increase the color set, try again in the next iteration.
      Kokkos::parallel_for(my_exec_space(0,this->nv), choose_colors(kok_colors, color_set, color_ban, tentative_color_ban));
      if (tictoc){
        color_time += timer->seconds();
        timer->reset();
      }
    }
    if (tictoc){
      std::cout <<
          "\tinit_time:" << inittime <<
          " mc:" << mc_time <<
          " cnt_time:" << cnt_time <<
          " ban_time:" << ban_time <<
          " expand ban time:" << expand_ban_time <<
          " pps time:" << pps_time <<
          " color time:" << color_time << std::endl<< std::endl;
    }

    //set the final colors.
    Kokkos::parallel_for(my_exec_space(0,this->nv), set_final_colors (kok_colors, color_set));

    num_loops = i;

    if (tictoc){
      delete timer;
    }
  }


  /*! \brief Functor to initialize the colors of the vertices randomly,
   *  with the hope that it will reduce the conflict in parallel execution.
   *  It also initializes the color bans.
   */
  struct init_colors{

    color_array_type kokcolors;
    color_temp_work_array_type color_ban; //colors
    color_type hash; //the number of colors to be assigned initially.

    //the value to initialize the color_ban_. We avoid using the first bit representing the sign.
    //Therefore if idx is int, it can represent 32-1 colors. Use color_set to represent more.
    color_type color_ban_init_val;


    init_colors (color_array_type colors,color_temp_work_array_type color_ban_,color_type hash_):
      kokcolors(colors), color_ban(color_ban_), hash(hash_){
      color_type tmp = 1;
      color_ban_init_val = tmp <<( sizeof(color_type) * 8 -1);
    };

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &ii) const {
      //set colors based on their indices.
      color_type tmp1 = 1;
      kokcolors(ii) = tmp1 << (ii % hash);
      color_ban(ii) = color_ban_init_val;
    }
  };


  /*! \brief Functor to initialize the worklist
     */
  struct init_work_arrays{
    idx_edge_array_type _edge_conflict_indices;
    char_array_type _edge_conflict_marker;

    init_work_arrays (
        idx_edge_array_type edge_conflict_indices,
        char_array_type edge_conflict_marker):
          _edge_conflict_indices(edge_conflict_indices),
          _edge_conflict_marker(edge_conflict_marker){};

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &ii) const {
      _edge_conflict_indices(ii)= ii; //every edge needs to be worked on initially.
      _edge_conflict_marker(ii) = 1; //every edge is a conflict initially.
    }
  };

  /**\brief Functor to mark the conflicts.
   * It goes to all the edges, and checks if two ends of an edge have the same color.
   * If they do, then it marks the one with the larger index as conflict.
   */
  struct halfedge_mark_conflicts {
    //edge list, source and destinations of the edge list.
    idx_persistent_work_array_type srcs, dsts;
    color_array_type kokcolors;

    idx_temp_work_array_type color_set; //the colors that are represented with bits, and the colors set that the color is in.
    color_temp_work_array_type color_ban, tentative_color_ban; //color ban for each vertex represented with bit, as well as tentative color ban.
    idx_temp_work_array_type edge_conf_indices;


    //      idx color_ban_init_val;  //the value to initialize the color_ban_. We avoid using the first bit representing the sign.
    //                                //Therefore if idx is int, it can represent 32-1 colors. Use color_set to represent more.

    /** \brief functor constructor.
     * \param srcs_ sources of the edgelist
     * \param dsts_ destinations of the edgelist
     * \param kokcolors_  the colors of the vertices. Represented with bits.
     * \param color_set_  the color set of the vertices. kokcolors_ and color_set_ together
     *      is used to represent the colors e.g. color_set_(v) * (numbits_in_idx-1) + set_bit_position_in_kokcolors_(v)
     * \param color_ban_ the bit representation of the neighbor colors that are in the same color_set.
     *                   color_ban_ only includes the colors of the neighbors that have been colored correctly.
     * \param tentative_color_ban_ :bit reprensentation of the neighbor vertex colors, that have been colored in the current iteration.
     *                              it is tentative, because coloring might have conflicts.
     * \param edge_conf_indices_ : The worklist for the edges.
     */
    halfedge_mark_conflicts (
        idx_persistent_work_array_type srcs_,
        idx_persistent_work_array_type dsts_,
        color_array_type kokcolors_,
        idx_temp_work_array_type color_set_,
        color_temp_work_array_type color_ban_,
        color_temp_work_array_type tentative_color_ban_,
        idx_temp_work_array_type edge_conf_indices_):
      srcs(srcs_), dsts (dsts_),
      kokcolors(kokcolors_), color_set(color_set_),
      color_ban(color_ban_), tentative_color_ban(tentative_color_ban_),
      edge_conf_indices(edge_conf_indices_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &ii) const {
      idx work_index = edge_conf_indices(ii);
      //traverse edges,
      idx src_id = srcs(work_index);
      idx dst_id = dsts(work_index);


      color_type source_color = kokcolors(src_id);
      color_type dst_color = kokcolors(dst_id);

      //if the source and destionation have the same color, e.g. same color and same color_set.
      //then we have a conflict.
      char is_conflicted = (source_color != 0 && (source_color == dst_color) && (color_set(src_id) == color_set(dst_id)));
      if (is_conflicted){
        //this functor works both sides, although there is a reverse edge that will be encountered.
        //we want to mark the conflicts as soon as possible, so that those conflicted vertices neighbors wont have unnecessary conflicts.
        //idx conflict_ver = (src_id < dst_id) ? src_id : dst_id;
        //
        //TODO: dst_id seems to reduce the num colors, without increaisng runtime
        kokcolors(dst_id) = 0;
        tentative_color_ban(dst_id) = 0;
      }
    }
  };

  /**\brief Functor to count the number of conflicts
   * Also, as a side effect, it marks edge_conflict_marker and
   * remove those that are not needed to be looked further
   */
  struct halfedge_conflict_count{
    idx_persistent_work_array_type _kok_src;
    idx_persistent_work_array_type _kok_dst;
    color_array_type _kok_colors;
    idx_temp_work_array_type _color_set; //the colors that are represented with bits, and the colors set that the color is in.
    idx_temp_work_array_type _edge_conflict_indices;
    char_array_type _edge_conflict_marker;

    halfedge_conflict_count(
        idx_persistent_work_array_type kok_src,
        idx_persistent_work_array_type kok_dst,
        color_array_type kok_colors,
        idx_temp_work_array_type color_set,
        idx_temp_work_array_type edge_conflict_indices,
        char_array_type edge_conflict_marker):
          _kok_src( kok_src),
          _kok_dst( kok_dst),
          _kok_colors( kok_colors),
          _color_set( color_set),
          _edge_conflict_indices( edge_conflict_indices),
          _edge_conflict_marker( edge_conflict_marker){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &ii, idx &sum) const {

      idx w = _edge_conflict_indices(ii);

      if (_edge_conflict_marker(w) == 0){
        sum += 1;
      }
      else{
        idx d = this->_kok_dst(w);
        idx s = this->_kok_src(w);

        color_type dc = _kok_colors(d);
        color_type sc = _kok_colors(s);

        if ( (dc && sc) || //if both colored
            (sc && (_color_set(d) > _color_set(s))) || //if source is colored, and destination color set is larger than source
            (dc && (_color_set(s) > _color_set(d))) //or if destionation is colored, and the source color set is larger
        ){
          //then no need to look at this edge anymore.
          _edge_conflict_marker(w) = 0;
          sum += 1;
        }
      }
    }
  };


  /**
   * \brief Functor to perform parallel prefix sum for edges so that the position
   * on the next conflictlist is calculated.
   */
  struct parallel_prefix_sum{
    idx_temp_work_array_type _edge_conflict_indices;
    char_array_type _edge_conflict_marker;
    idx_temp_work_array_type _pps_view;

    parallel_prefix_sum(
        idx_temp_work_array_type edge_conflict_indices,
        char_array_type edge_conflict_marker,
        idx_temp_work_array_type pps_view):
          _edge_conflict_indices(edge_conflict_indices),
          _edge_conflict_marker(edge_conflict_marker),
          _pps_view(pps_view){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii, size_t& update, const bool final) const {
      idx w = _edge_conflict_indices(ii);
      if (final) {
        _pps_view(w) = idx (update);
      }
      update += _edge_conflict_marker(w);
    }
  };

  /**
   * \brief Functor to create the new work array.
   */
  struct create_new_work_array{
    idx_temp_work_array_type _edge_conflict_indices;
    char_array_type _edge_conflict_marker;
    idx_temp_work_array_type _pps_view;
    idx_temp_work_array_type _new_edge_conflict_indices;

    create_new_work_array(
        idx_temp_work_array_type edge_conflict_indices,
        char_array_type edge_conflict_marker,
        idx_temp_work_array_type pps_view,
        idx_temp_work_array_type new_edge_conflict_indices):
          _edge_conflict_indices(edge_conflict_indices),
          _edge_conflict_marker(edge_conflict_marker),
          _pps_view(pps_view),
          _new_edge_conflict_indices(new_edge_conflict_indices){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii) const {
      idx w = _edge_conflict_indices(ii);
      if(_edge_conflict_marker(w)){
        idx future_index = _pps_view(w);
        _new_edge_conflict_indices(future_index) = w;
      }
    }
  };

  /**
   * \brief Functor to create the new work array with atomic operations.
   */
  struct atomic_create_new_work_array{
    idx_type _new_index;
    idx_temp_work_array_type _edge_conflict_indices;
    char_array_type _edge_conflict_marker;
    idx_temp_work_array_type _new_edge_conflict_indices;

    atomic_create_new_work_array(
        idx_type new_index,
        idx_temp_work_array_type edge_conflict_indices,
        char_array_type edge_conflict_marker,
        idx_temp_work_array_type new_edge_conflict_indices):
          _new_index(new_index),
          _edge_conflict_indices(edge_conflict_indices),
          _edge_conflict_marker(edge_conflict_marker),
          _new_edge_conflict_indices(new_edge_conflict_indices){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx ii) const {
      idx w = _edge_conflict_indices(ii);
      if(_edge_conflict_marker(w)){
        const idx future_index = Kokkos::atomic_fetch_add( &_new_index(), 1);
        _new_edge_conflict_indices(future_index) = w;
      }
    }

  };


  /** \brief Functor for creating the ban colors for uncolored vertices.
   * It only creates ban_color, based on the certain information.
   * That is it works on the edges where one end is colored (for good)
   * and the other part is not colored.
   */
  struct halfedge_ban_colors {
    idx_persistent_work_array_type srcs, dsts;  //edge list, source and destinations of the edge list.
    color_array_type kokcolors;
    idx_temp_work_array_type color_set; //the colors that are represented with bits, and the colors set that the color is in.
    color_temp_work_array_type color_ban; //color ban for each vertex represented with bit
    idx_temp_work_array_type conflict_indices;
    char_array_type edge_conflict_marker;

    /** \brief Functor constructor.
     * \param srcs_ sources of the edgelist
     * \param dsts_ destinations of the edgelist
     * \param kokcolors_  the colors of the vertices. Represented with bits.
     * \param color_set_  the color set of the vertices. kokcolors_ and color_set_ together
     *      is used to represent the colors e.g. color_set_(v) * (numbits_in_idx-1) + set_bit_position_in_kokcolors_(v)
     * \param color_ban_ the bit representation of the neighbor colors that are in the same color_set.
     *                   color_ban_ only includes the colors of the neighbors that have been colored correctly.
     */
    halfedge_ban_colors (
        idx_persistent_work_array_type srcs_, idx_persistent_work_array_type dsts_,
        color_array_type kokcolors_, idx_temp_work_array_type  color_set_,
        color_temp_work_array_type color_ban_,
        idx_temp_work_array_type conflict_indices_, char_array_type edge_conflict_marker_):
          srcs(srcs_), dsts (dsts_),
          kokcolors(kokcolors_), color_set(color_set_),
          color_ban(color_ban_),
          conflict_indices(conflict_indices_), edge_conflict_marker(edge_conflict_marker_){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &ii) const {
      idx work_index = conflict_indices(ii);
      idx dst_id = dsts(work_index);
      color_type dst_col = kokcolors(dst_id);
      idx src_id = srcs(work_index);
      color_type src_col = kokcolors(src_id);

      //check destionation color.
      //continue only if it is not colored
      if ((!dst_col  && src_col) || (!src_col  && dst_col)){
        //check src color, send its color to ban colors only if it is colored.
        idx dest_col_set = color_set (dst_id);
        idx src_col_set = color_set (src_id);
        //check if they are in the same color set.
        //if they are not, we do not ban the color, as it represents a different color.
        if (src_col_set == dest_col_set){
          //atomic or, as no threads owns 'dst' (neither src)
          idx uncolored_vertex = dst_col? src_id: dst_id;
          Kokkos::atomic_fetch_or<color_type>(&(color_ban(uncolored_vertex)), src_col | dst_col);
          edge_conflict_marker(work_index) = 0;
        }
      }
    }
  };


  /**
   * \brief Functor to tentatively color vertices. It propogates the color information
   * to other end.
   */
  struct halfedge_expand_ban_for_unmatched_neighbors{
    idx_persistent_work_array_type srcs, dsts; //edge list, source and destinations of the edge list.
    color_array_type kokcolors;
    idx_temp_work_array_type color_set; //the colors that are represented with bits, and the colors set that the color is in.
    color_temp_work_array_type color_ban, tentative_color_ban; //color ban for each vertex represented with bit, as well as tentative color ban.
    color_type first_digit;
    idx_temp_work_array_type conflict_indices;


    /** \brief functor constructor.
     * \param srcs_ sources of the edgelist
     * \param dsts_ destinations of the edgelist
     * \param kokcolors_  the colors of the vertices. Represented with bits.
     * \param color_set_  the color set of the vertices. kokcolors_ and color_set_ together
     *      is used to represent the colors e.g. color_set_(v) * (numbits_in_idx-1) + set_bit_position_in_kokcolors_(v)
     * \param color_ban_ the bit representation of the neighbor colors that are in the same color_set.
     *                   color_ban_ only includes the colors of the neighbors that have been colored correctly.
     * \param tentative_color_ban_ :bit reprensentation of the neighbor vertex colors, that have been colored in the current iteration.
     *                              it is tentative, because coloring might have conflicts.
     */
    halfedge_expand_ban_for_unmatched_neighbors (
        idx_persistent_work_array_type srcs_, idx_persistent_work_array_type dsts_,
        color_array_type kokcolors_, idx_temp_work_array_type  color_set_ ,
        color_temp_work_array_type color_ban_, color_temp_work_array_type tentative_color_ban_,
        idx_temp_work_array_type conflict_indices_):
          srcs(srcs_), dsts (dsts_),
          kokcolors(kokcolors_), color_set(color_set_),
          color_ban(color_ban_), tentative_color_ban(tentative_color_ban_),
          conflict_indices(conflict_indices_){
      color_type tmp = 1;
      first_digit = tmp <<( sizeof(color_type) * 8 -1);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &ii) const {
      idx work_index = conflict_indices(ii);
      idx dst_id = dsts(work_index);
      color_type dst_col = kokcolors(dst_id);

      //if the destionation is colored already, we have nothing to do.
      //otherwise, if destionation is uncolored, or if its color < 0 (it has been tentatively colored)
      //then we need to check the source.
      if (dst_col == 0 || (dst_col & first_digit) ){
        idx src_id = srcs(work_index);
        color_type src_col = kokcolors(src_id);
        //if source is colored, again we have nothing to do.
        //if it is tentatively colored or uncolored, then we have work to do.
        if (src_col == 0 || (src_col & first_digit)){
          //check their colors sets, if they are on different color sets,
          //we dont need to care about the prohibted colors on each other -- at least in this iteration.
          idx dest_col_set = color_set (dst_id);
          idx src_col_set = color_set (src_id);

          if (src_col_set == dest_col_set){
            if ((dst_col & first_digit)  && (src_col & first_digit)){
              //if both ends are tentatively colored, we can check for the conflict here,
              //and fix it as below. However doing so increased the number of colors,
              //so, it has been turned of for now.

              if (src_col == dst_col && dest_col_set == src_col_set){
                //idx smaller_index = (src_id > dst_id) ? src_id : dst_id;
                idx smaller_index = dst_id; //TODO which one is better? this seems to be not much changing
                //idx smaller_index = src_id;
                //then both have been colored tentavitely. propoagate the color of src to dst.
                Kokkos::atomic_fetch_or<color_type>(&(tentative_color_ban(smaller_index)), -src_col);
                idx banned_colors = ~(color_ban(smaller_index) | tentative_color_ban(smaller_index));
                idx larger_col = banned_colors & (-banned_colors);
                kokcolors(smaller_index) = -(larger_col);
              }
            }
            else if(src_col != 0) {
              //if src is tentavily colored, and dst is not colored,
              //then we send the color information to dst's tentative_ban.

              //Kokkos::atomic_fetch_or<color_type>(&(color_ban(dst_id)), -src_col);
              Kokkos::atomic_fetch_or<color_type>(&(tentative_color_ban(dst_id)), -src_col);
            }
            else if (dst_col != 0){
              //if it is dst tentatively colors, but src is not colored,
              //then we send the dst color info to src's tentative_ban

              //Kokkos::atomic_fetch_or<color_type>(&(color_ban(src_id)), -dst_col);
              Kokkos::atomic_fetch_or<color_type>(&(tentative_color_ban(src_id)), -dst_col);
            }
            else {
              //idx smaller_index = src_id < dst_id > 0 ? src_id: dst_id;
              //idx larger_index = src_id < dst_id > 0 ? dst_id : src_id;
#ifndef TOOHIGHQUALITY
              idx smaller_index = src_id;
              idx larger_index = dst_id;
#endif
#ifdef TOOHIGHQUALITY
              idx smaller_index = dst_id;
              idx larger_index = src_id;
#endif

              //idx smaller_col =  src_id < dst_id > 0 ? src_col: dst_col;
              //if both ends are uncolored, tentatively color the the source if its index is smaller than dst.
              //make an 'bitwise or' of color_ban and tentative_color_ban to get the all prohibited colors.
              //we need to find the right most zero here. it is easier to find right most 1, so we do a not of the result color ban.
              color_type banned_colors = ~(color_ban(smaller_index) | tentative_color_ban(smaller_index));
              //the 'bitwise and' of banned_colors with two's complement result in only the rightmost 1 to be set, which is our color.
              src_col = banned_colors & (-banned_colors);
              //set it to minus of the color, as it is tentative coloring.
              kokcolors(smaller_index) = -(src_col);
              //send the color information to dst's tentative color ban.
              Kokkos::atomic_fetch_or<color_type>(&(tentative_color_ban(larger_index)), src_col);
              //Kokkos::atomic_fetch_or<color_type>(&(color_ban(dst_id)), src_col);
            }
          }
        }
      }
    }
  };





  /** \brief Functor responsible for choosing a color for each uncolored vertex,
   * given the color_ban and tentative_color_ban
   */
  struct choose_colors {
    color_array_type kokcolors;
    idx_temp_work_array_type color_set; //the colors that are represented with bits, and the colors set that the color is in.
    color_temp_work_array_type color_ban, tentative_color_ban;  //color ban for each vertex represented with bit, as well as tentative color ban.
    color_type color_ban_init_val;  //the value to initialize the color_ban_. We avoid using the first bit representing the sign.
    //Therefore if idx is int, it can represent 32-1 colors. Use color_set to represent more

    /** \brief functor constructor.
     * \param kokcolors_  the colors of the vertices. Represented with bits.
     * \param color_set_  the color set of the vertices. kokcolors_ and color_set_ together
     *      is used to represent the colors e.g. color_set_(v) * (numbits_in_idx-1) + set_bit_position_in_kokcolors_(v)
     * \param color_ban_ the bit representation of the neighbor colors that are in the same color_set.
     *                   color_ban_ only includes the colors of the neighbors that have been colored correctly.
     * \param tentative_color_ban_ :bit reprensentation of the neighbor vertex colors, that have been colored in the current iteration.
     *                              it is tentative, because coloring might have conflicts.
     */
    choose_colors ( color_array_type kokcolors_, idx_array_type  color_set_,
        color_array_type color_ban_,  color_array_type tentative_color_ban_):
          kokcolors(kokcolors_), color_set(color_set_),
          color_ban(color_ban_), tentative_color_ban(tentative_color_ban_){
      //color_ban should always have 1 at the first bit, so that that color is not allowed.
      color_type tmp = 1;
      color_ban_init_val = tmp <<( sizeof(color_type) * 8 -1);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &ii) const {
      //if the vertex is uncolored, we will choose a new color for the vertex.
      if (kokcolors[ii] == 0){
        color_type certain_info = color_ban(ii);
        //get the banned_color_set by taking 'bitwise or' or color ban and tentative_color_ban
        color_type banned_colors = ~(certain_info | tentative_color_ban(ii));
        //my color is the first non set bit in the banned_colors. We perform a not operation,
        //and make a 'bitwise and' with its 2's complement to find the first zero bit.
        color_type my_color = banned_colors & (-banned_colors);
        if (my_color == 0){
#ifdef EBCOLORING_HIGHER_QUALITY
          //if my color is zero, that is all the available colors in this set has been taken by the neighbors
          //then I might need to change my color set. But we need to be sure about this, so we check if
          //color_ban is full as well, since tentative_color_ban might be too pessimist.
          banned_colors = ~(certain_info);
          my_color = banned_colors & (-banned_colors);
          if (my_color == 0){
#endif
            //if all colors are taken by the neighbors (certainly excluding tentative colors), then
            //I need to change my color set.
            //if there are still available colors w.r.t. color_ban, then try one more time without increasing the color_set.
            color_set(ii) += 1; // increase color set.
            color_ban(ii) = color_ban_init_val; //set the color ban to its initial value.
#ifdef EBCOLORING_HIGHER_QUALITY
          }
#endif
          //in each case we cannot color this vertex. set the tentative_color_ban to 0
          //try to color it at the next iteration.
          tentative_color_ban(ii) = 0;
          //color_ban(ii) = color_ban_init_val; //set the color ban to its initial value.
        }
        else {
          kokcolors(ii) = my_color;
        }
      }
      else if (kokcolors(ii) & color_ban_init_val){
        kokcolors(ii) = -kokcolors(ii);
      }
    }
  };

  /** Functor responsible for setting the final color for each vertex,
   *  The color of the vertex is found ast color_set * (sizeof(color_type) * 8 -1) + log2(color)
   */
  struct set_final_colors{
    color_array_type kokcol;
    idx_temp_work_array_type kokcolset; //the colors that are represented with bits, and the colors set that the color is in.
    color_type color_size;

    /** \brief functor constructor.
     * \param kokcol_  the colors of the vertices. Represented with bits.
     * \param kokcolset_  the color set of the vertices. kokcolors_ and color_set_ together
     *      is used to represent the colors e.g. color_set_(v) * (numbits_in_idx-1) + set_bit_position_in_kokcolors_(v)
     */
    set_final_colors(color_array_type kokcol_, idx_temp_work_array_type  kokcolset_): kokcol(kokcol_),kokcolset(kokcolset_), color_size ( sizeof(color_type) * 8 -1){}

    KOKKOS_INLINE_FUNCTION
    void operator()(const idx &ii) const {
      idx i = 0;
      color_type val = kokcol(ii);
      //if check below is necessary.
      // this happens when a vertices all neighbors are colored,
      //so the information from all neighbors are taken, no edge to be processed by this vertex.
      //the algorithm works on the number of edges, if the edges are all consumed, the loop
      //might terminate with an early exit without coloring this vertex.
      //this happens when all neighbors consumes all the colors in the current vertex set,
      //and the vertex left to be colored in the next iteration.
      //but the vertex couldnt be colored, because there is no more edge left to be worked on.

      if(val == 0) val = 1;

      //find the position in the bit.
      while (val != 0) {
        ++i;
        val = val >> 1;
      }
      //set the final color.
      kokcol(ii) = i + kokcolset(ii) * color_size;
    }
  };
};

}
}
}
}

#endif //_KOKKOSCOLORINGIMP_HPP
