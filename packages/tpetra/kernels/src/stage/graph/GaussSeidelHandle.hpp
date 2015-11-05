
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_Core.hpp>

#ifndef _GAUSSSEIDELHANDLE_HPP
#define _GAUSSSEIDELHANDLE_HPP
//#define VERBOSE
namespace Experimental{

namespace KokkosKernels{
namespace Graph{

enum GSAlgorithm{GS_DEFAULT, GS_PERMUTED, GS_TEAM};

template <class idx_array_type_, class idx_edge_array_type_, class value_array_type_, class ExecutionSpace, class TemporaryMemorySpace, class PersistentMemorySpace>
class GaussSeidelHandle{


public:
  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;

  typedef idx_array_type_ idx_array_type;
  typedef idx_edge_array_type_ idx_edge_array_type;
  typedef value_array_type_ value_array_type;

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

  typedef typename value_array_type::value_type value_type;
  typedef typename value_array_type::array_layout value_type_array_layout;
  typedef typename value_array_type::device_type value_type_device_type;
  typedef typename value_array_type::memory_traits value_type_memory_traits;
  typedef typename value_array_type::HostMirror host_value_view_type; //Host view type

  typedef typename Kokkos::View<idx *, HandleTempMemorySpace> idx_temp_work_array_type;
  typedef typename Kokkos::View<idx *, HandlePersistentMemorySpace> idx_persistent_work_array_type;
  typedef typename idx_persistent_work_array_type::HostMirror host_idx_persistent_view_type; //Host view type

  typedef typename Kokkos::View<value_type *, HandleTempMemorySpace> value_temp_work_array_type;
  typedef typename Kokkos::View<value_type *, HandlePersistentMemorySpace> value_persistent_work_array_type;

private:
  bool owner_of_coloring;
  GSAlgorithm algorithm_type;

  host_idx_persistent_view_type color_set_xadj;
  idx_persistent_work_array_type color_sets;
  idx numColors;

  idx_persistent_work_array_type permuted_xadj;
  idx_persistent_work_array_type permuted_adj;
  value_persistent_work_array_type permuted_adj_vals;
  idx_persistent_work_array_type old_to_new_map;

  bool called_symbolic;


  value_persistent_work_array_type permuted_y_vector;
  value_persistent_work_array_type permuted_x_vector;

  int suggested_vector_size;
  int suggested_team_size;
  public:

  /**
   * \brief Default constructor.
   */
  GaussSeidelHandle(GSAlgorithm gs = GS_DEFAULT):
    owner_of_coloring(false),
    algorithm_type(gs),
    color_set_xadj(), color_sets(), numColors(0),
    permuted_xadj(),  permuted_adj(), permuted_adj_vals(), old_to_new_map(),
    called_symbolic(false), permuted_y_vector(), permuted_x_vector(),
    suggested_vector_size(0), suggested_team_size(0)
    {
    if (gs == GS_DEFAULT){
      this->choose_default_algorithm();
    }


  }

    /** \brief Chooses best algorithm based on the execution space. COLORING_EB if cuda, COLORING_VB otherwise.
   */
  void choose_default_algorithm(){
#if defined( KOKKOS_HAVE_SERIAL )
    if (Kokkos::Impl::is_same< Kokkos::Serial , ExecutionSpace >::value){
      this->algorithm_type = GS_PERMUTED;
#ifdef VERBOSE
      std::cout << "Serial Execution Space, Default Algorithm: GS_PERMUTED" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_PTHREAD )
    if (Kokkos::Impl::is_same< Kokkos::Threads , ExecutionSpace >::value){
      this->algorithm_type = GS_PERMUTED;
#ifdef VERBOSE
      std::cout << "PTHREAD Execution Space, Default Algorithm: GS_PERMUTED" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_OPENMP )
    if (Kokkos::Impl::is_same< Kokkos::OpenMP, ExecutionSpace >::value){
      this->algorithm_type = GS_PERMUTED;
#ifdef VERBOSE
      std::cout << "OpenMP Execution Space, Default Algorithm: GS_PERMUTED" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_CUDA )
    if (Kokkos::Impl::is_same<Kokkos::Cuda, ExecutionSpace >::value){
      this->algorithm_type = GS_TEAM;
#ifdef VERBOSE
      std::cout << "Qthread Execution Space, Default Algorithm: GS_TEAM" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_QTHREAD)
    if (Kokkos::Impl::is_same< Kokkos::Qthread, ExecutionSpace >::value){
      this->algorithm_type = GS_PERMUTED;
#ifdef VERBOSE
      std::cout << "Qthread Execution Space, Default Algorithm: GS_PERMUTED" << std::endl;
#endif
    }
#endif
  }

  virtual ~GaussSeidelHandle(){};

  //getters
  GSAlgorithm get_algorithm_type() const {return this->algorithm_type;}
  bool is_owner_of_coloring() const {return this->owner_of_coloring;}

  host_idx_persistent_view_type get_color_xadj() {
    return this->color_set_xadj;
  }
  idx_persistent_work_array_type get_color_adj() {
    return this->color_sets;
  }
  idx get_num_colors() {
    return this->numColors;
  }

  idx_persistent_work_array_type get_new_xadj() {
    return this->permuted_xadj;
  }
  idx_persistent_work_array_type get_new_adj() {
    return this->permuted_adj;
  }
  value_persistent_work_array_type get_new_adj_val() {
    return this->permuted_adj_vals;
  }
  idx_persistent_work_array_type get_old_to_new_map() {
    return this->old_to_new_map;
  }

  bool is_symbolic_called(){return this->called_symbolic;}

  //setters
  void set_algorithm_type(const GSAlgorithm &sgs_algo){this->algorithm_type = sgs_algo;}
  void set_owner_of_coloring(bool owner = true){this->owner_of_coloring = owner;}

  void set_call_symbolic(bool call = true){this->called_symbolic = call;}

  void set_color_set_xadj(const host_idx_persistent_view_type &color_set_xadj_) {
    this->color_set_xadj = color_set_xadj_;
  }
  void set_color_set_adj(const idx_persistent_work_array_type &color_sets_) {
    this->color_sets = color_sets_;
  }
  void set_num_colors(const idx &numColors_) {
    this->numColors = numColors_;
  }

  void set_new_xadj(const idx_persistent_work_array_type &xadj_) {
    this->permuted_xadj = xadj_;
  }
  void set_new_adj(const idx_persistent_work_array_type &adj_) {
    this->permuted_adj = adj_;
  }
  void set_new_adj_val(const value_persistent_work_array_type &adj_vals_) {
    this->permuted_adj_vals = adj_vals_;
  }
  void set_old_to_new_map(const idx_persistent_work_array_type &old_to_new_map_) {
    this->old_to_new_map = old_to_new_map_;
  }

  void allocate_x_y_vectors(){
    if(permuted_y_vector.dimension_0() != permuted_xadj.dimension_0() - 1){
      permuted_y_vector = value_persistent_work_array_type("PERMUTED Y VECTOR", permuted_xadj.dimension_0() - 1);
    }
    if(permuted_x_vector.dimension_0() != permuted_xadj.dimension_0() - 1){
      permuted_x_vector = value_persistent_work_array_type("PERMUTED X VECTOR", permuted_xadj.dimension_0() - 1);
    }
  }

  value_persistent_work_array_type get_permuted_y_vector (){return this->permuted_y_vector;}
  value_persistent_work_array_type get_permuted_x_vector (){return this->permuted_x_vector;}
  void vector_team_size(int max_allowed_team_size, int &suggested_vector_size_, int &suggested_team_size_){
    //suggested_team_size_ =  this->suggested_team_size = 1;
    //suggested_vector_size_=this->suggested_vector_size = 1;
    //return;
    if (this->suggested_team_size && this->suggested_vector_size) {
      suggested_vector_size_ = this->suggested_vector_size;
      suggested_team_size_ = this->suggested_team_size;
      return;
    }

#if defined( KOKKOS_HAVE_SERIAL )
    if (Kokkos::Impl::is_same< Kokkos::Serial , ExecutionSpace >::value){
      suggested_vector_size_ = this->suggested_vector_size = 1;
      suggested_team_size_ = this->suggested_team_size = max_allowed_team_size;
      return;
    }
#endif

#if defined( KOKKOS_HAVE_PTHREAD )
    if (Kokkos::Impl::is_same< Kokkos::Threads , ExecutionSpace >::value){
      suggested_vector_size_ = this->suggested_vector_size = 1;
      suggested_team_size_ = this->suggested_team_size = max_allowed_team_size;
      return;
    }
#endif

#if defined( KOKKOS_HAVE_OPENMP )
    if (Kokkos::Impl::is_same< Kokkos::OpenMP, ExecutionSpace >::value){
      suggested_vector_size_ = this->suggested_vector_size = 1;
      suggested_team_size_ = this->suggested_team_size = max_allowed_team_size;
    }
#endif

#if defined( KOKKOS_HAVE_CUDA )
    if (Kokkos::Impl::is_same<Kokkos::Cuda, ExecutionSpace >::value){

      this->suggested_vector_size = this->entries.dimension_0() / double (this->row_map.dimension_0()) + 0.5;

      if (this->suggested_vector_size <= 3){
        this->suggested_vector_size = 2;
      }
      else if (this->suggested_vector_size <= 6){
        this->suggested_vector_size = 4;
      }
      else if (this->suggested_vector_size <= 12){
        this->suggested_vector_size = 8;
      }
      else if (this->suggested_vector_size <= 24){
        this->suggested_vector_size = 16;
      }
      else {
        this->suggested_vector_size = 32;
      }
      if (max_allowed_team_size < 32){
        std::cerr << "max_allowed_team_size:" << max_allowed_team_size << std::endl;
      }
      suggested_vector_size_ = this->suggested_vector_size;
      this->suggested_team_size= suggested_team_size_ = max_allowed_team_size / this->suggested_vector_size;
    }
#endif

#if defined( KOKKOS_HAVE_QTHREAD)
    if (Kokkos::Impl::is_same< Kokkos::Qthread, ExecutionSpace >::value){
      suggested_vector_size_ = this->suggested_vector_size = 1;
      suggested_team_size_ = this->suggested_team_size = max_allowed_team_size;
    }
#endif

  }

};
}
}
}

#endif
