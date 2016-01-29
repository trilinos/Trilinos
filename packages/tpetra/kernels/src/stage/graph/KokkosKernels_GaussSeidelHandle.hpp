
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_Core.hpp>

#ifndef _GAUSSSEIDELHANDLE_HPP
#define _GAUSSSEIDELHANDLE_HPP
//#define VERBOSE

namespace KokkosKernels{

namespace Experimental{

namespace Graph{

enum GSAlgorithm{GS_DEFAULT, GS_PERMUTED, GS_TEAM};

template <class row_index_view_type_,
          class nonzero_index_view_type_,
          class nonzero_value_view_type_,
          class ExecutionSpace,
          class TemporaryMemorySpace,
          class PersistentMemorySpace>
class GaussSeidelHandle{
public:
  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;

  typedef row_index_view_type_ in_row_index_view_type;
  typedef nonzero_index_view_type_ in_nonzero_index_view_type;
  typedef nonzero_value_view_type_ in_nonzero_value_view_type;

  typedef typename in_row_index_view_type::non_const_value_type row_index_type;
  typedef typename in_row_index_view_type::array_layout row_view_array_layout;
  typedef typename in_row_index_view_type::device_type row_view_device_type;
  typedef typename in_row_index_view_type::memory_traits row_view_memory_traits;
  typedef typename in_row_index_view_type::HostMirror row_host_view_type; //Host view type
  //typedef typename idx_memory_traits::MemorySpace MyMemorySpace;

  typedef typename in_nonzero_index_view_type::non_const_value_type nonzero_index_type;
  typedef typename in_nonzero_index_view_type::array_layout nonzero_index_view_array_layout;
  typedef typename in_nonzero_index_view_type::device_type nonzero_index_view_device_type;
  typedef typename in_nonzero_index_view_type::memory_traits nonzero_index_view_memory_traits;
  typedef typename in_nonzero_index_view_type::HostMirror nonzero_index_host_view_type; //Host view type
  //typedef typename idx_edge_memory_traits::MemorySpace MyEdgeMemorySpace;

  typedef typename in_nonzero_value_view_type::non_const_value_type nonzero_value_type;
  typedef typename in_nonzero_value_view_type::array_layout nonzero_value_view_array_layout;
  typedef typename in_nonzero_value_view_type::device_type nonzero_value_view_device_type;
  typedef typename in_nonzero_value_view_type::memory_traits nonzero_value_view_memory_traits;
  typedef typename in_nonzero_value_view_type::HostMirror nonzero_value_host_view_type; //Host view type



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



  typedef typename in_nonzero_value_view_type::const_data_type const_nonzero_value_data_type;
  typedef typename in_nonzero_value_view_type::non_const_data_type non_const_nonzero_value_data_type;
  typedef typename in_nonzero_value_view_type::memory_space nonzero_value_view_memory_space;
  typedef typename Kokkos::View<const_nonzero_value_data_type, nonzero_value_view_array_layout,
        nonzero_value_view_memory_space, nonzero_value_view_memory_traits> const_nonzero_value_view_type;
    typedef typename Kokkos::View<non_const_nonzero_value_data_type, nonzero_value_view_array_layout,
        nonzero_value_view_memory_space, nonzero_value_view_memory_traits> non_const_nonzero_value_view_type;



  //typedef typename Kokkos::View<row_index_type *, HandleTempMemorySpace> idx_temp_work_array_type;
  typedef typename Kokkos::View<row_index_type *, HandleTempMemorySpace> row_index_temp_work_view_type;
  //typedef typename Kokkos::View<row_index_type *, HandlePersistentMemorySpace> idx_persistent_work_array_type;
  typedef typename Kokkos::View<row_index_type *, HandlePersistentMemorySpace> row_index_persistent_work_view_type;

  //typedef typename row_index_persistent_work_view_type::HostMirror host_idx_persistent_view_type; //Host view type
  typedef typename row_index_persistent_work_view_type::HostMirror row_index_persistent_host_view_type; //Host view type

  //typedef typename Kokkos::View<nonzero_value_type *, HandleTempMemorySpace> value_temp_work_array_type;
  typedef typename Kokkos::View<nonzero_value_type *, HandleTempMemorySpace> nonzero_value_temp_work_view_type;

  //typedef typename Kokkos::View<nonzero_value_type *, HandlePersistentMemorySpace> value_persistent_work_array_type;
  typedef typename Kokkos::View<nonzero_value_type *, HandlePersistentMemorySpace> nonzero_value_persistent_work_view_type;

private:
  bool owner_of_coloring;
  GSAlgorithm algorithm_type;

  row_index_persistent_host_view_type color_set_xadj;
  row_index_persistent_work_view_type color_sets;
  row_index_type numColors;

  row_index_persistent_work_view_type permuted_xadj;
  row_index_persistent_work_view_type permuted_adj;
  nonzero_value_persistent_work_view_type permuted_adj_vals;
  row_index_persistent_work_view_type old_to_new_map;

  bool called_symbolic;
  bool called_numeric;


  nonzero_value_persistent_work_view_type permuted_y_vector;
  nonzero_value_persistent_work_view_type permuted_x_vector;

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
    called_symbolic(false), called_numeric(false), permuted_y_vector(), permuted_x_vector(),
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

  row_index_persistent_host_view_type get_color_xadj() {
    return this->color_set_xadj;
  }
  row_index_persistent_work_view_type get_color_adj() {
    return this->color_sets;
  }
  row_index_type get_num_colors() {
    return this->numColors;
  }

  row_index_persistent_work_view_type get_new_xadj() {
    return this->permuted_xadj;
  }
  row_index_persistent_work_view_type get_new_adj() {
    return this->permuted_adj;
  }
  nonzero_value_persistent_work_view_type get_new_adj_val() {
    return this->permuted_adj_vals;
  }
  row_index_persistent_work_view_type get_old_to_new_map() {
    return this->old_to_new_map;
  }

  bool is_symbolic_called(){return this->called_symbolic;}
  bool is_numeric_called(){return this->called_numeric;}

  //setters
  void set_algorithm_type(const GSAlgorithm &sgs_algo){this->algorithm_type = sgs_algo;}
  void set_owner_of_coloring(bool owner = true){this->owner_of_coloring = owner;}

  void set_call_symbolic(bool call = true){this->called_symbolic = call;}
  void set_call_numeric(bool call = true){this->called_numeric = call;}

  void set_color_set_xadj(const row_index_persistent_host_view_type &color_set_xadj_) {
    this->color_set_xadj = color_set_xadj_;
  }
  void set_color_set_adj(const row_index_persistent_work_view_type &color_sets_) {
    this->color_sets = color_sets_;
  }
  void set_num_colors(const row_index_type &numColors_) {
    this->numColors = numColors_;
  }

  void set_new_xadj(const row_index_persistent_work_view_type &xadj_) {
    this->permuted_xadj = xadj_;
  }
  void set_new_adj(const row_index_persistent_work_view_type &adj_) {
    this->permuted_adj = adj_;
  }
  void set_new_adj_val(const nonzero_value_persistent_work_view_type &adj_vals_) {
    this->permuted_adj_vals = adj_vals_;
  }
  void set_old_to_new_map(const row_index_persistent_work_view_type &old_to_new_map_) {
    this->old_to_new_map = old_to_new_map_;
  }

  void allocate_x_y_vectors(row_index_type num_rows, row_index_type num_cols){
    if(permuted_y_vector.dimension_0() != num_rows){
      permuted_y_vector = nonzero_value_persistent_work_view_type("PERMUTED Y VECTOR", num_rows);
    }
    if(permuted_x_vector.dimension_0() != num_cols){
      permuted_x_vector = nonzero_value_persistent_work_view_type("PERMUTED X VECTOR", num_cols);
    }
  }

  nonzero_value_persistent_work_view_type get_permuted_y_vector (){return this->permuted_y_vector;}
  nonzero_value_persistent_work_view_type get_permuted_x_vector (){return this->permuted_x_vector;}
  void vector_team_size(
      int max_allowed_team_size,
      int &suggested_vector_size_,
      int &suggested_team_size_,
      row_index_type nr, row_index_type nnz){
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

      this->suggested_vector_size = nnz / double (nr) + 0.5;

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
