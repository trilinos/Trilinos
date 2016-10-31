
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_Core.hpp>
#include <iostream>

//#define KERNELS_HAVE_CUSPARSE

#ifdef KERNELS_HAVE_CUSPARSE
#include "cusparse.h"
#endif
#ifndef _SPGEMMHANDLE_HPP
#define _SPGEMMHANDLE_HPP
//#define VERBOSE

namespace KokkosKernels{

namespace Experimental{

namespace Graph{

enum SPGEMMAlgorithm{SPGEMM_DEFAULT, SPGEMM_CUSPARSE, SPGEMM_SERIAL, SPGEMM_CUSP, SPGEMM_MKL, SPGEMM_VIENNA,
                     SPGEMM_KK1, SPGEMM_KK2, SPGEMM_KK3, SPGEMM_KK4,
                     SPGEMM_KK_SPEED, SPGEMM_KK_MEMORY, SPGEMM_KK_COLOR, SPGEMM_KK_MULTICOLOR, SPGEMM_KK_MULTICOLOR2, SPGEMM_KK_MEMSPEED};

template <class lno_row_view_t_,
          class lno_nnz_view_t_,
          class scalar_nnz_view_t_,
          class ExecutionSpace,
          class TemporaryMemorySpace,
          class PersistentMemorySpace>
class SPGEMMHandle{
public:
  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;

  typedef lno_row_view_t_ in_lno_row_view_t;
  typedef lno_nnz_view_t_ in_lno_nnz_view_t;
  typedef scalar_nnz_view_t_ in_scalar_nnz_view_t;

  typedef typename in_lno_row_view_t::non_const_value_type size_type;
  /*
  typedef typename in_lno_row_view_t::array_layout row_lno_view_array_layout;
  typedef typename in_lno_row_view_t::device_type row_lno_view_device_t;
  typedef typename in_lno_row_view_t::memory_traits row_lno_view_memory_traits;
  typedef typename in_lno_row_view_t::HostMirror row_lno_host_view_t; //Host view type
  */
  //typedef typename idx_memory_traits::MemorySpace MyMemorySpace;

  typedef typename in_lno_nnz_view_t::non_const_value_type nnz_lno_t;
  /*
  typedef typename in_lno_nnz_view_t::array_layout nnz_lno_view_array_layout;
  typedef typename in_lno_nnz_view_t::device_type nnz_lno_view_device_t;
  typedef typename in_lno_nnz_view_t::memory_traits nnz_lno_view_memory_traits;
  typedef typename in_lno_nnz_view_t::HostMirror nnz_lno_host_view_t; //Host view type
  */
  //typedef typename idx_edge_memory_traits::MemorySpace MyEdgeMemorySpace;

  typedef typename in_scalar_nnz_view_t::non_const_value_type nnz_scalar_t;
  /*
  typedef typename in_scalar_nnz_view_t::array_layout nnz_scalar_view_array_layout;
  typedef typename in_scalar_nnz_view_t::device_type nnz_scalar_view_device_t;
  typedef typename in_scalar_nnz_view_t::memory_traits nnz_scalar_view_memory_traits;
  typedef typename in_scalar_nnz_view_t::HostMirror nnz_scalar_view_t; //Host view type
  */


  /*
  typedef typename in_lno_row_view_t::const_data_type const_row_lno_t;
  typedef typename in_lno_row_view_t::non_const_data_type non_const_row_lno_t;
  */

  typedef typename in_lno_row_view_t::const_type const_lno_row_view_t;
  typedef typename in_lno_row_view_t::non_const_type non_const_lno_row_view_t;

  /*
  typedef typename in_lno_nnz_view_t::const_data_type const_nnz_lno_t;
  typedef typename in_lno_nnz_view_t::non_const_data_type non_const_nnz_lno_t;
  */
  typedef typename in_lno_nnz_view_t::const_type const_lno_nnz_view_t;
  typedef typename in_lno_nnz_view_t::non_const_type non_const_lno_nnz_view_t;

  /*
  typedef typename in_scalar_nnz_view_t::const_data_type const_nnz_scalar_t;
  typedef typename in_scalar_nnz_view_t::non_const_data_type non_const_nnz_scalar_t;
  */
  typedef typename in_scalar_nnz_view_t::const_type const_scalar_nnz_view_t;
  typedef typename in_scalar_nnz_view_t::non_const_type non_const_scalar_nnz_view_t;


  typedef typename Kokkos::View<size_type *, HandleTempMemorySpace> row_lno_temp_work_view_t;
  typedef typename Kokkos::View<size_type *, HandlePersistentMemorySpace> row_lno_persistent_work_view_t;
  typedef typename row_lno_persistent_work_view_t::HostMirror row_lno_persistent_work_host_view_t; //Host view type

  typedef typename Kokkos::View<nnz_scalar_t *, HandleTempMemorySpace> scalar_temp_work_view_t;
  typedef typename Kokkos::View<nnz_scalar_t *, HandlePersistentMemorySpace> scalar_persistent_work_view_t;


  typedef typename Kokkos::View<nnz_lno_t *, HandleTempMemorySpace> nnz_lno_temp_work_view_t;
  typedef typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace> nnz_lno_persistent_work_view_t;
  typedef typename nnz_lno_persistent_work_view_t::HostMirror nnz_lno_persistent_work_host_view_t; //Host view type


#ifdef KERNELS_HAVE_CUSPARSE
  struct cuSparseHandleType{
    cusparseHandle_t handle;
    cusparseOperation_t transA;
    cusparseOperation_t transB;
    cusparseMatDescr_t a_descr;
    cusparseMatDescr_t b_descr;
    cusparseMatDescr_t c_descr;
    cuSparseHandleType(bool transposeA, bool transposeB){
      cusparseStatus_t status;
      status= cusparseCreate(&handle);
      if (status != CUSPARSE_STATUS_SUCCESS) {
        std::cerr << ("cusparseCreate ERROR") << std::endl;
        return;
      }
      cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST);

      if (transposeA){
        transA = CUSPARSE_OPERATION_TRANSPOSE;
      }
      else {
        transA  = CUSPARSE_OPERATION_NON_TRANSPOSE;
      }
      if (transposeB){
        transB = CUSPARSE_OPERATION_TRANSPOSE;
      }
      else {
        transB  = CUSPARSE_OPERATION_NON_TRANSPOSE;
      }


      status = cusparseCreateMatDescr(&a_descr);
      if (status != CUSPARSE_STATUS_SUCCESS) {
        std::cerr << "cusparseCreateMatDescr a_descr ERROR" << std::endl;
        return;
      }
      cusparseSetMatType(a_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
      cusparseSetMatIndexBase(a_descr,CUSPARSE_INDEX_BASE_ZERO);

      status = cusparseCreateMatDescr(&b_descr);
      if (status != CUSPARSE_STATUS_SUCCESS) {
        std::cerr << ("cusparseCreateMatDescr b_descr ERROR") << std::endl;
        return;
      }
      cusparseSetMatType(b_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
      cusparseSetMatIndexBase(b_descr,CUSPARSE_INDEX_BASE_ZERO);

      status = cusparseCreateMatDescr(&c_descr);
      if (status != CUSPARSE_STATUS_SUCCESS) {
        std::cerr << ("cusparseCreateMatDescr  c_descr ERROR") << std::endl;
        return;
      }
      cusparseSetMatType(c_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
      cusparseSetMatIndexBase(c_descr,CUSPARSE_INDEX_BASE_ZERO);
    }
    ~cuSparseHandleType(){
      cusparseDestroyMatDescr(a_descr);
      cusparseDestroyMatDescr(b_descr);
      cusparseDestroyMatDescr(c_descr);
      cusparseDestroy(handle);
    }
  };

  typedef cuSparseHandleType SPGEMMcuSparseHandleType;
#endif
private:
  SPGEMMAlgorithm algorithm_type;
  size_type result_nnz_size;

  bool called_symbolic;
  bool called_numeric;

  int suggested_vector_size;
  int suggested_team_size;
  nnz_lno_t max_nnz_inresult;
  nnz_lno_t max_nnz_compressed_result;

  row_lno_temp_work_view_t compressed_b_rowmap;// compressed_b_set_begins, compressed_b_set_nexts;
  nnz_lno_temp_work_view_t compressed_b_set_indices, compressed_b_sets;
  row_lno_temp_work_view_t compressed_c_rowmap;


  row_lno_temp_work_view_t tranpose_a_xadj, tranpose_b_xadj, tranpose_c_xadj;
  nnz_lno_temp_work_view_t tranpose_a_adj, tranpose_b_adj, tranpose_c_adj;

  bool transpose_a,transpose_b, transpose_c_symbolic;


  nnz_lno_t num_colors;
  nnz_lno_persistent_work_host_view_t color_xadj;
  nnz_lno_persistent_work_view_t color_adj, vertex_colors;
  nnz_lno_t num_multi_colors, num_used_colors;

  double multi_color_scale;


#ifdef KERNELS_HAVE_CUSPARSE
  SPGEMMcuSparseHandleType *cuSPARSEHandle;
#endif
  public:
  void set_color_xadj(
      nnz_lno_t num_colors_,
      nnz_lno_persistent_work_host_view_t color_xadj_,
      nnz_lno_persistent_work_view_t color_adj_,
      nnz_lno_persistent_work_view_t vertex_colors_,
      nnz_lno_t num_multi_colors_, nnz_lno_t num_used_colors_){

    num_colors = num_colors_;
    color_xadj = color_xadj_;
    color_adj = color_adj_;
    vertex_colors = vertex_colors_;

    num_multi_colors = num_multi_colors_;
    num_used_colors = num_used_colors_;
  }

  /**
   * \brief sets the result nnz size.
   * \param result_nnz_size: size of the output matrix.
   */
  void set_c_nnz(size_type result_nnz_size_){
    this->result_nnz_size = result_nnz_size_;
  }
  /**
   * \brief returns the result nnz size.
   */
  void get_c_nnz(){
    return this->result_nnz_size;
  }

  void set_multi_color_scale(double multi_color_scale_){
    this->multi_color_scale = multi_color_scale_;
  }

  double get_multi_color_scale(){
    return this->multi_color_scale;
  }

  void get_color_xadj(
      nnz_lno_t &num_colors_,
      nnz_lno_persistent_work_host_view_t &color_xadj_,
      nnz_lno_persistent_work_view_t &color_adj_,
      nnz_lno_persistent_work_view_t &vertex_colors_,
      nnz_lno_t &num_multi_colors_, nnz_lno_t &num_used_colors_){
    num_colors_ = num_colors;
    color_xadj_ = color_xadj;
    color_adj_ = color_adj ;
    num_multi_colors_ = num_multi_colors;
    num_used_colors_ = num_used_colors ;
    vertex_colors_ = vertex_colors;
  }

  void set_compressed_c(
      row_lno_temp_work_view_t compressed_c_rowmap_){
    compressed_c_rowmap = compressed_c_rowmap_;
  }

  void get_compressed_c(
      row_lno_temp_work_view_t &compressed_c_rowmap_){
    compressed_c_rowmap_ = compressed_c_rowmap;
  }

  //TODO: store transpose here.
  void get_c_transpose_symbolic(){}





  void get_compressed_b(
      row_lno_temp_work_view_t &compressed_b_rowmap_,
      nnz_lno_temp_work_view_t &compressed_b_set_indices_,
      nnz_lno_temp_work_view_t &compressed_b_sets_,
      row_lno_temp_work_view_t &compressed_b_set_begins_,
      row_lno_temp_work_view_t &compressed_b_set_nexts_){
    compressed_b_rowmap_ = compressed_b_rowmap;
    compressed_b_set_indices_ = compressed_b_set_indices;
    compressed_b_sets_ = compressed_b_sets;
  }

  /**
   * \brief Default constructor.
   */
  SPGEMMHandle(SPGEMMAlgorithm gs = SPGEMM_DEFAULT):
    algorithm_type(gs), result_nnz_size(0),
    called_symbolic(false), called_numeric(false),
    suggested_vector_size(0), suggested_team_size(0), max_nnz_inresult(0),
    tranpose_a_xadj(), tranpose_b_xadj(), tranpose_c_xadj(),
    tranpose_a_adj(), tranpose_b_adj(), tranpose_c_adj(),
    transpose_a(false),transpose_b(false), transpose_c_symbolic(false),
    num_colors(0),
    color_xadj(), color_adj(), vertex_colors(), num_multi_colors(0),num_used_colors(0), multi_color_scale(1)
#ifdef KERNELS_HAVE_CUSPARSE
  ,cuSPARSEHandle(NULL)
#endif
  {
    if (gs == SPGEMM_DEFAULT){
      this->choose_default_algorithm();
    }
  }


  virtual ~SPGEMMHandle(){

#ifdef KERNELS_HAVE_CUSgPARSE
    this->destroy_cuSPARSE_Handle();
#endif
  };

#ifdef KERNELS_HAVE_CUSPARSE
  void create_cuSPARSE_Handle(bool transA, bool transB){
    this->destroy_cuSPARSE_Handle();
    this->cuSPARSEHandle = new cuSparseHandleType(transA, transB);
  }
  void destroy_cuSPARSE_Handle(){
    if (this->cuSPARSEHandle != NULL){
      delete this->cuSPARSEHandle;
      this->cuSPARSEHandle = NULL;
    }
  }

  SPGEMMcuSparseHandleType *get_cuSparseHandle(){
    return this->cuSPARSEHandle;
  }
#endif
    /** \brief Chooses best algorithm based on the execution space. COLORING_EB if cuda, COLORING_VB otherwise.
   */
  void choose_default_algorithm(){
#if defined( KOKKOS_HAVE_SERIAL )
    if (Kokkos::Impl::is_same< Kokkos::Serial , ExecutionSpace >::value){
      this->algorithm_type = SPGEMM_SERIAL;
#ifdef VERBOSE
      std::cout << "Serial Execution Space, Default Algorithm: SPGEMM_SERIAL" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_PTHREAD )
    if (Kokkos::Impl::is_same< Kokkos::Threads , ExecutionSpace >::value){
      this->algorithm_type = SPGEMM_SERIAL;
#ifdef VERBOSE
      std::cout << "PTHREAD Execution Space, Default Algorithm: SPGEMM_SERIAL" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_OPENMP )
    if (Kokkos::Impl::is_same< Kokkos::OpenMP, ExecutionSpace >::value){
      this->algorithm_type = SPGEMM_SERIAL;
#ifdef VERBOSE
      std::cout << "OpenMP Execution Space, Default Algorithm: SPGEMM_SERIAL" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_CUDA )
    if (Kokkos::Impl::is_same<Kokkos::Cuda, ExecutionSpace >::value){
      this->algorithm_type = SPGEMM_CUSPARSE;
#ifdef VERBOSE
      std::cout << "Cuda Execution Space, Default Algorithm: SPGEMM_CUSPARSE" << std::endl;
#endif
    }
#endif

#if defined( KOKKOS_HAVE_QTHREAD)
    if (Kokkos::Impl::is_same< Kokkos::Qthread, ExecutionSpace >::value){
      this->algorithm_type = SPGEMM_SERIAL;
#ifdef VERBOSE
      std::cout << "Qthread Execution Space, Default Algorithm: SPGEMM_SERIAL" << std::endl;
#endif
    }
#endif
  }






  //getters
  SPGEMMAlgorithm get_algorithm_type() const {return this->algorithm_type;}

  bool is_symbolic_called(){return this->called_symbolic;}
  bool is_numeric_called(){return this->called_numeric;}


  nnz_lno_t get_max_result_nnz() const{
    return this->max_nnz_inresult ;
  }

  nnz_lno_t get_max_compresed_result_nnz() const{
    return this->max_nnz_compressed_result ;
  }


  //setters
  void set_algorithm_type(const SPGEMMAlgorithm &sgs_algo){this->algorithm_type = sgs_algo;}
  void set_call_symbolic(bool call = true){this->called_symbolic = call;}
  void set_call_numeric(bool call = true){this->called_numeric = call;}

  void set_max_result_nnz(nnz_lno_t num_result_nnz_){
    this->max_nnz_inresult = num_result_nnz_;
  }

  void set_max_compresed_result_nnz(nnz_lno_t num_result_nnz_){
    this->max_nnz_compressed_result = num_result_nnz_;
  }

  void vector_team_size(
      int max_allowed_team_size,
      int &suggested_vector_size_,
      int &suggested_team_size_,
      size_type nr, size_type nnz){
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
