
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

enum SPGEMMAlgorithm{SPGEMM_DEFAULT, SPGEMM_CUSPARSE, SPGEMM_SERIAL, SPGEMM_CUSP, SPGEMM_MKL};

template <class row_index_view_type_,
          class nonzero_index_view_type_,
          class nonzero_value_view_type_,
          class ExecutionSpace,
          class TemporaryMemorySpace,
          class PersistentMemorySpace>
class SPGEMMHandle{
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


  typedef typename Kokkos::View<row_index_type *, HandleTempMemorySpace> row_index_temp_work_view_type;
  typedef typename Kokkos::View<row_index_type *, HandlePersistentMemorySpace> row_index_persistent_work_view_type;
  typedef typename row_index_persistent_work_view_type::HostMirror row_index_persistent_host_view_type; //Host view type

  typedef typename Kokkos::View<nonzero_value_type *, HandleTempMemorySpace> nonzero_value_temp_work_view_type;
  typedef typename Kokkos::View<nonzero_value_type *, HandlePersistentMemorySpace> nonzero_value_persistent_work_view_type;

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

  bool called_symbolic;
  bool called_numeric;

  int suggested_vector_size;
  int suggested_team_size;

#ifdef KERNELS_HAVE_CUSPARSE
  SPGEMMcuSparseHandleType *cuSPARSEHandle;
#endif
  public:
  /**
   * \brief Default constructor.
   */
  SPGEMMHandle(SPGEMMAlgorithm gs = SPGEMM_DEFAULT):
    algorithm_type(gs),
    called_symbolic(false), called_numeric(false),
    suggested_vector_size(0), suggested_team_size(0)
#ifdef KERNELS_HAVE_CUSPARSE
  ,cuSPARSEHandle(NULL)
#endif
  {
    if (gs == SPGEMM_DEFAULT){
      this->choose_default_algorithm();
    }
  }


  virtual ~SPGEMMHandle(){

#ifdef KERNELS_HAVE_CUSPARSE
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

  //setters
  void set_algorithm_type(const SPGEMMAlgorithm &sgs_algo){this->algorithm_type = sgs_algo;}
  void set_call_symbolic(bool call = true){this->called_symbolic = call;}
  void set_call_numeric(bool call = true){this->called_numeric = call;}

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
