
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

enum SPGEMMAlgorithm{SPGEMM_DEFAULT, SPGEMM_CUSPARSE, SPGEMM_SERIAL, SPGEMM_CUSP};

template <class idx_array_type_,
          class idx_edge_array_type_,
          class value_array_type_,
          class ExecutionSpace,
          class TemporaryMemorySpace,
          class PersistentMemorySpace>
class SPGEMMHandle{
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
      idx nr, idx nnz){
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
