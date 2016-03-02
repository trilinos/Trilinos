
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_Core.hpp>
#include <KokkosKernels_Utils.hpp>
#ifndef _GAUSSSEIDELHANDLE_HPP
#define _GAUSSSEIDELHANDLE_HPP
//#define VERBOSE

namespace KokkosKernels{

namespace Experimental{

namespace Graph{

enum GSAlgorithm{GS_DEFAULT, GS_PERMUTED, GS_TEAM};

template <class lno_row_view_t_,
          class lno_nnz_view_t_,
          class scalar_nnz_view_t_,
          class ExecutionSpace,
          class TemporaryMemorySpace,
          class PersistentMemorySpace>
class GaussSeidelHandle{
public:
  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;

  typedef lno_row_view_t_ in_lno_row_view_t;
  typedef lno_nnz_view_t_ in_lno_nnz_view_t;
  typedef scalar_nnz_view_t_ in_scalar_nnz_view_t;

  typedef typename in_lno_row_view_t::non_const_value_type row_lno_t;
  typedef typename in_lno_row_view_t::array_layout row_lno_view_array_layout;
  typedef typename in_lno_row_view_t::device_type row_lno_view_device_t;
  typedef typename in_lno_row_view_t::memory_traits row_lno_view_memory_traits;
  typedef typename in_lno_row_view_t::HostMirror row_lno_host_view_t; //Host view type
  //typedef typename idx_memory_traits::MemorySpace MyMemorySpace;

  typedef typename in_lno_nnz_view_t::non_const_value_type nnz_lno_t;
  typedef typename in_lno_nnz_view_t::array_layout nnz_lno_view_array_layout;
  typedef typename in_lno_nnz_view_t::device_type nnz_lno_view_device_t;
  typedef typename in_lno_nnz_view_t::memory_traits nnz_lno_view_memory_traits;
  typedef typename in_lno_nnz_view_t::HostMirror nnz_lno_host_view_t; //Host view type
  //typedef typename idx_edge_memory_traits::MemorySpace MyEdgeMemorySpace;

  typedef typename in_scalar_nnz_view_t::non_const_value_type nnz_scalar_t;
  typedef typename in_scalar_nnz_view_t::array_layout nnz_scalar_view_array_layout;
  typedef typename in_scalar_nnz_view_t::device_type nnz_scalar_view_device_t;
  typedef typename in_scalar_nnz_view_t::memory_traits nnz_scalar_view_memory_traits;
  typedef typename in_scalar_nnz_view_t::HostMirror nnz_scalar_view_t; //Host view type



  typedef typename in_lno_row_view_t::const_data_type const_row_lno_t;
  typedef typename in_lno_row_view_t::non_const_data_type non_const_row_lno_t;
  typedef typename in_lno_row_view_t::const_type const_lno_row_view_t;
  typedef typename in_lno_row_view_t::non_const_type non_const_lno_row_view_t;




  typedef typename in_lno_nnz_view_t::const_data_type const_nnz_lno_t;
  typedef typename in_lno_nnz_view_t::non_const_data_type non_const_nnz_lno_t;
  typedef typename in_lno_nnz_view_t::const_type const_lno_nnz_view_t;
  typedef typename in_lno_nnz_view_t::non_const_type non_const_lno_nnz_view_t;



  typedef typename in_scalar_nnz_view_t::const_data_type const_nnz_scalar_t;
  typedef typename in_scalar_nnz_view_t::non_const_data_type non_const_nnz_scalar_t;
  typedef typename in_scalar_nnz_view_t::const_type const_scalar_nnz_view_t;
  typedef typename in_scalar_nnz_view_t::non_const_type non_const_scalar_nnz_view_t;


  //typedef typename Kokkos::View<row_index_type *, HandleTempMemorySpace> idx_temp_work_array_type;
  typedef typename Kokkos::View<row_lno_t *, HandleTempMemorySpace> row_lno_temp_work_view_t;
  //typedef typename Kokkos::View<row_index_type *, HandlePersistentMemorySpace> idx_persistent_work_array_type;
  typedef typename Kokkos::View<row_lno_t *, HandlePersistentMemorySpace> row_lno_persistent_work_view_t;

  //typedef typename row_index_persistent_work_view_type::HostMirror host_idx_persistent_view_type; //Host view type
  typedef typename row_lno_persistent_work_view_t::HostMirror row_lno_persistent_work_host_view_t; //Host view type

  //typedef typename Kokkos::View<nonzero_value_type *, HandleTempMemorySpace> value_temp_work_array_type;
  typedef typename Kokkos::View<nnz_scalar_t *, HandleTempMemorySpace> scalar_temp_work_view_t;

  //typedef typename Kokkos::View<nonzero_value_type *, HandlePersistentMemorySpace> value_persistent_work_array_type;
  typedef typename Kokkos::View<nnz_scalar_t *, HandlePersistentMemorySpace> scalar_persistent_work_view_t;


  typedef typename Kokkos::View<nnz_lno_t *, HandleTempMemorySpace> nnz_lno_temp_work_view_t;
  typedef typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace> nnz_lno_persistent_work_view_t;
  typedef typename nnz_lno_persistent_work_view_t::HostMirror nnz_lno_persistent_work_host_view_t; //Host view type


private:
  bool owner_of_coloring;
  GSAlgorithm algorithm_type;

  row_lno_persistent_work_host_view_t color_set_xadj;
  row_lno_persistent_work_view_t color_sets;
  row_lno_t numColors;

  row_lno_persistent_work_view_t permuted_xadj;
  row_lno_persistent_work_view_t permuted_adj;
  scalar_persistent_work_view_t permuted_adj_vals;
  row_lno_persistent_work_view_t old_to_new_map;

  bool called_symbolic;
  bool called_numeric;


  scalar_persistent_work_view_t permuted_y_vector;
  scalar_persistent_work_view_t permuted_x_vector;

  int suggested_vector_size;
  int suggested_team_size;

  scalar_persistent_work_view_t permuted_diagonals;
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
    suggested_vector_size(0), suggested_team_size(0), permuted_diagonals()
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

  row_lno_persistent_work_host_view_t get_color_xadj() {
    return this->color_set_xadj;
  }
  row_lno_persistent_work_view_t get_color_adj() {
    return this->color_sets;
  }
  row_lno_t get_num_colors() {
    return this->numColors;
  }

  row_lno_persistent_work_view_t get_new_xadj() {
    return this->permuted_xadj;
  }
  row_lno_persistent_work_view_t get_new_adj() {
    return this->permuted_adj;
  }
  scalar_persistent_work_view_t get_new_adj_val() {
    return this->permuted_adj_vals;
  }
  row_lno_persistent_work_view_t get_old_to_new_map() {
    return this->old_to_new_map;
  }

  bool is_symbolic_called(){return this->called_symbolic;}
  bool is_numeric_called(){return this->called_numeric;}

  //setters
  void set_algorithm_type(const GSAlgorithm &sgs_algo){this->algorithm_type = sgs_algo;}
  void set_owner_of_coloring(bool owner = true){this->owner_of_coloring = owner;}

  void set_call_symbolic(bool call = true){this->called_symbolic = call;}
  void set_call_numeric(bool call = true){this->called_numeric = call;}

  void set_color_set_xadj(const row_lno_persistent_work_host_view_t &color_set_xadj_) {
    this->color_set_xadj = color_set_xadj_;
  }
  void set_color_set_adj(const row_lno_persistent_work_view_t &color_sets_) {
    this->color_sets = color_sets_;
  }
  void set_num_colors(const row_lno_t &numColors_) {
    this->numColors = numColors_;
  }

  void set_new_xadj(const row_lno_persistent_work_view_t &xadj_) {
    this->permuted_xadj = xadj_;
  }
  void set_new_adj(const row_lno_persistent_work_view_t &adj_) {
    this->permuted_adj = adj_;
  }
  void set_new_adj_val(const scalar_persistent_work_view_t &adj_vals_) {
    this->permuted_adj_vals = adj_vals_;
  }
  void set_old_to_new_map(const row_lno_persistent_work_view_t &old_to_new_map_) {
    this->old_to_new_map = old_to_new_map_;
  }
  void set_permuted_diagonals (const scalar_persistent_work_view_t permuted_diagonals_){
    this->permuted_diagonals = permuted_diagonals_;
  }

  scalar_persistent_work_view_t get_permuted_diagonals (){
    return this->permuted_diagonals;
  }

  void allocate_x_y_vectors(row_lno_t num_rows, row_lno_t num_cols){
    if(permuted_y_vector.dimension_0() != size_t(num_rows)){
      permuted_y_vector = scalar_persistent_work_view_t("PERMUTED Y VECTOR", num_rows);
    }
    if(permuted_x_vector.dimension_0() != size_t(num_cols)){
      permuted_x_vector = scalar_persistent_work_view_t("PERMUTED X VECTOR", num_cols);
    }
  }

  scalar_persistent_work_view_t get_permuted_y_vector (){return this->permuted_y_vector;}
  scalar_persistent_work_view_t get_permuted_x_vector (){return this->permuted_x_vector;}
  void vector_team_size(
      int max_allowed_team_size,
      int &suggested_vector_size_,
      int &suggested_team_size_,
      row_lno_t nr, row_lno_t nnz){
    //suggested_team_size_ =  this->suggested_team_size = 1;
    //suggested_vector_size_=this->suggested_vector_size = 1;
    //return;
    if (this->suggested_team_size && this->suggested_vector_size) {
      suggested_vector_size_ = this->suggested_vector_size;
      suggested_team_size_ = this->suggested_team_size;
      return;
    }
    else {
      KokkosKernels::Experimental::Util::get_suggested_vector_team_size<row_lno_t, ExecutionSpace>(
          max_allowed_team_size, suggested_vector_size_, suggested_team_size_, nr, nnz);
      this->suggested_team_size = suggested_vector_size_;
      this->suggested_vector_size = suggested_vector_size_;

    }
  }

};
}
}
}

#endif
