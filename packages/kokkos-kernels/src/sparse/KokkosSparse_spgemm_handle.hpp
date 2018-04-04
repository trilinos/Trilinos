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
#include <iostream>
#include <string>

//#define KOKKOSKERNELS_ENABLE_TPL_CUSPARSE

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "cusparse.h"
#endif
#ifndef _SPGEMMHANDLE_HPP
#define _SPGEMMHANDLE_HPP
//#define VERBOSE

namespace KokkosSparse{

//TODO:SPGEMM_KK_MEMORY2 option is for testing in openmp.
//it wont work on cuda, not bind to a test program.
//hidden parameter for StringToSPGEMMAlgorithm for now.
enum SPGEMMAlgorithm{
		/*DEFAULT*/SPGEMM_KK, SPGEMM_KK_DENSE, SPGEMM_KK_MEMORY, SPGEMM_KK_LP, //KKVARIANTS
		SPGEMM_CUSPARSE,  SPGEMM_CUSP, SPGEMM_MKL, SPGEMM_MKL2PHASE, SPGEMM_VIENNA, //TPLS

		//TRIANGLE COUNTING SPECIALIZED
		SPGEMM_KK_TRIANGLE_AI, //SPGEMM_KK_TRIANGLE_DEFAULT, SPGEMM_KK_TRIANGLE_MEM, SPGEMM_KK_TRIANGLE_DENSE,
		SPGEMM_KK_TRIANGLE_IA_UNION, //SPGEMM_KK_TRIANGLE_DEFAULT_IA_UNION, SPGEMM_KK_TRIANGLE_MEM_IA_UNION, SPGEMM_KK_TRIANGLE_DENSE_IA_UNION,
		SPGEMM_KK_TRIANGLE_IA,//SPGEMM_KK_TRIANGLE_IA_DEFAULT, SPGEMM_KK_TRIANGLE_IA_MEM, SPGEMM_KK_TRIANGLE_IA_DENSE,
		SPGEMM_KK_TRIANGLE_LL,
		SPGEMM_KK_TRIANGLE_LU,

		//below research code.
		SPGEMM_KK_MULTIMEM, SPGEMM_KK_OUTERMULTIMEM,
		SPGEMM_DEFAULT, SPGEMM_DEBUG, SPGEMM_SERIAL,
		SPGEMM_KK_CUCKOO, //USE CUCKOO HASHING
		SPGEMM_KK_TRACKED_CUCKOO, //USE SCALED CUCKOO HASHING
		SPGEMM_KK_TRACKED_CUCKOO_F, //USE SCALED FANCY CUCKOO HASHING
		SPGEMM_KK_SPEED,  // DENSE ACCUMULATOR SAME AS SPEED
		SPGEMM_KK_MEMORY_SORTED,
		SPGEMM_KK_MEMORY_TEAM,
		SPGEMM_KK_MEMORY_BIGTEAM,
		SPGEMM_KK_MEMORY_SPREADTEAM,
		SPGEMM_KK_MEMORY_BIGSPREADTEAM,
		SPGEMM_KK_MEMORY2,
		SPGEMM_KK_COLOR,
		SPGEMM_KK_MULTICOLOR,
		SPGEMM_KK_MULTICOLOR2,
		SPGEMM_KK_MEMSPEED};

enum SPGEMMAccumulator{
  SPGEMM_ACC_DEFAULT, SPGEMM_ACC_DENSE, SPGEMM_ACC_SPARSE,
};
template <class size_type_, class lno_t_, class scalar_t_,
          class ExecutionSpace,
          class TemporaryMemorySpace,
          class PersistentMemorySpace>
class SPGEMMHandle{
public:
  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;


  typedef typename std::remove_const<size_type_>::type  size_type;
  typedef const size_type const_size_type;

  typedef typename std::remove_const<lno_t_>::type  nnz_lno_t;
  typedef const nnz_lno_t const_nnz_lno_t;

  typedef typename std::remove_const<scalar_t_>::type  nnz_scalar_t;
  typedef const nnz_scalar_t const_nnz_scalar_t;


  typedef typename Kokkos::View<size_type *, HandleTempMemorySpace> row_lno_temp_work_view_t;
  typedef typename Kokkos::View<size_type *, HandlePersistentMemorySpace> row_lno_persistent_work_view_t;
  typedef typename row_lno_persistent_work_view_t::HostMirror row_lno_persistent_work_host_view_t; //Host view type

  typedef typename Kokkos::View<nnz_scalar_t *, HandleTempMemorySpace> scalar_temp_work_view_t;
  typedef typename Kokkos::View<nnz_scalar_t *, HandlePersistentMemorySpace> scalar_persistent_work_view_t;


  typedef typename Kokkos::View<nnz_lno_t *, HandleTempMemorySpace> nnz_lno_temp_work_view_t;
  typedef typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace> nnz_lno_persistent_work_view_t;
  typedef typename nnz_lno_persistent_work_view_t::HostMirror nnz_lno_persistent_work_host_view_t; //Host view type




#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
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
        throw std::runtime_error ("cusparseCreate ERROR\n");
        //return;
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
        throw std::runtime_error ("cusparseCreateMatDescr a_descr ERROR\n");
        //return;
      }
      cusparseSetMatType(a_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
      cusparseSetMatIndexBase(a_descr,CUSPARSE_INDEX_BASE_ZERO);

      status = cusparseCreateMatDescr(&b_descr);
      if (status != CUSPARSE_STATUS_SUCCESS) {
        throw std::runtime_error ("cusparseCreateMatDescr b_descr ERROR\n");
        //return;
      }
      cusparseSetMatType(b_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
      cusparseSetMatIndexBase(b_descr,CUSPARSE_INDEX_BASE_ZERO);

      status = cusparseCreateMatDescr(&c_descr);
      if (status != CUSPARSE_STATUS_SUCCESS) {
        throw std::runtime_error ("cusparseCreateMatDescr  c_descr ERROR\n");
        //return;
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
  SPGEMMAccumulator accumulator_type;
  size_type result_nnz_size;

  bool called_symbolic;
  bool called_numeric;

  int suggested_vector_size;
  int suggested_team_size;
  nnz_lno_t max_nnz_inresult;
  nnz_lno_t max_nnz_compressed_result;

  size_type compressed_b_size;
  row_lno_temp_work_view_t compressed_b_rowmap;// compressed_b_set_begins, compressed_b_set_nexts;
  nnz_lno_temp_work_view_t compressed_b_set_indices, compressed_b_sets;

  row_lno_temp_work_view_t compressed_c_rowmap;

  nnz_lno_temp_work_view_t c_column_indices;

  row_lno_temp_work_view_t tranpose_a_xadj, tranpose_b_xadj, tranpose_c_xadj;
  nnz_lno_temp_work_view_t tranpose_a_adj, tranpose_b_adj, tranpose_c_adj;

  bool transpose_a,transpose_b, transpose_c_symbolic;


  nnz_lno_t num_colors;
  nnz_lno_persistent_work_host_view_t color_xadj;
  nnz_lno_persistent_work_view_t color_adj, vertex_colors;
  nnz_lno_t num_multi_colors, num_used_colors;
  nnz_lno_persistent_work_view_t min_result_row_for_each_row;


  bool create_lower_triangular;
  int sort_lower_triangular; //0 - do not sort // 1 - sort // 2 - Algorithm decides (default)
  int sort_option ;
  nnz_lno_persistent_work_view_t lower_triangular_permutation;

  row_lno_persistent_work_view_t lower_triangular_matrix_rowmap;
  nnz_lno_persistent_work_view_t lower_triangular_matrix_entries;


  row_lno_persistent_work_view_t incidence_matrix_row_map;
  nnz_lno_persistent_work_view_t incidence_matrix_entries;
  bool compress_second_matrix;


  double multi_color_scale;
  int mkl_sort_option;
  bool calculate_read_write_cost;

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  SPGEMMcuSparseHandleType *cuSPARSEHandle;
#endif
  public:

  std::string coloring_input_file;
  std::string coloring_output_file;

  int min_hash_size_scale;
  double compression_cut_off;
  double first_level_hash_cut_off;
  size_t original_max_row_flops, original_overall_flops;
  row_lno_persistent_work_view_t row_flops;

  size_t compressed_max_row_flops, compressed_overall_flops;

  void set_first_level_hash_cut_off(double first_level_hash_cut_off_){
    this->first_level_hash_cut_off = first_level_hash_cut_off_;
  }

  double get_first_level_hash_cut_off(){
    return this->first_level_hash_cut_off;
  }

  void set_compression_cut_off(double compression_cut_off_){
    this->compression_cut_off = compression_cut_off_;
  }

  double get_compression_cut_off(){
    return this->compression_cut_off;
  }
  void set_min_hash_size_scale(int scale){
    min_hash_size_scale = scale;
  }
  int get_min_hash_size_scale(){
    return min_hash_size_scale;
  }
  void set_read_write_cost_calc(bool read_write_cost_cal){
    this->calculate_read_write_cost = read_write_cost_cal;
  }
  int get_read_write_cost_calc(){
    return this->calculate_read_write_cost;
  }

  typename Kokkos::View<int *, HandlePersistentMemorySpace> persistent_c_xadj, persistent_a_xadj, persistent_b_xadj, persistent_a_adj, persistent_b_adj;
  size_t MaxColDenseAcc;
  bool mkl_keep_output;
  bool mkl_convert_to_1base;
  bool is_compression_single_step;

  void set_mkl_sort_option(int mkl_sort_option_){
    this->mkl_sort_option = mkl_sort_option_;
  }
  int get_mkl_sort_option(){
    return this->mkl_sort_option;
  }
  void set_c_column_indices(nnz_lno_temp_work_view_t c_col_indices_){
    this->c_column_indices = c_col_indices_;
  }

  nnz_lno_temp_work_view_t get_c_column_indices(){
    return this->c_column_indices;
  }

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
  size_type get_c_nnz(){
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

  
  void set_sort_lower_triangular(int option){
    this->sort_lower_triangular = option;
  }
  int get_sort_lower_triangular(){
    return this->sort_lower_triangular;
  }

  void set_sort_option(int option){
    this->sort_option = option;
  }
  int get_sort_option(){
    return this->sort_option;
  }

  void set_create_lower_triangular(bool option){
    this->create_lower_triangular = option;
  }
  bool get_create_lower_triangular(){
    return this->create_lower_triangular;
  }

  void set_lower_triangular_permutation(nnz_lno_persistent_work_view_t ltp_){
    this->lower_triangular_permutation = ltp_;
  }

  nnz_lno_persistent_work_view_t get_lower_triangular_permutation(){
    return this->lower_triangular_permutation;
  }

  void set_lower_triangular_matrix(
      row_lno_persistent_work_view_t lower_triangular_matrix_rowmap_,
      nnz_lno_persistent_work_view_t lower_triangular_matrix_entries_){
    this->lower_triangular_matrix_rowmap = lower_triangular_matrix_rowmap_;
    this->lower_triangular_matrix_entries = lower_triangular_matrix_entries_;
  }
  void get_lower_triangular_matrix(
      row_lno_persistent_work_view_t &lower_triangular_matrix_rowmap_,
      nnz_lno_persistent_work_view_t &lower_triangular_matrix_entries_){
    lower_triangular_matrix_rowmap_ = this->lower_triangular_matrix_rowmap;
    lower_triangular_matrix_entries_ = this->lower_triangular_matrix_entries;
  }


  void set_compressed_b(
      size_type b_nnz_size,
      row_lno_temp_work_view_t compressed_b_rowmap_,
      nnz_lno_temp_work_view_t compressed_b_set_indices_,
      nnz_lno_temp_work_view_t compressed_b_sets_){
    compressed_b_size = b_nnz_size;
    compressed_b_rowmap = compressed_b_rowmap_;
    compressed_b_set_indices = compressed_b_set_indices_;
    compressed_b_sets = compressed_b_sets_;
  }


  void get_compressed_b(
      size_type &b_nnz_size,
      row_lno_temp_work_view_t &compressed_b_rowmap_,
      nnz_lno_temp_work_view_t &compressed_b_set_indices_,
      nnz_lno_temp_work_view_t &compressed_b_sets_){
    b_nnz_size = compressed_b_size;
    compressed_b_rowmap_ = compressed_b_rowmap;
    compressed_b_set_indices_ = compressed_b_set_indices;
    compressed_b_sets_ = compressed_b_sets;
  }

  /**
   * \brief Default constructor.
   */
  SPGEMMHandle(SPGEMMAlgorithm gs = SPGEMM_DEFAULT):
    algorithm_type(gs), accumulator_type(SPGEMM_ACC_DEFAULT), result_nnz_size(0),
    called_symbolic(false), called_numeric(false),
    suggested_vector_size(0), suggested_team_size(0), max_nnz_inresult(0),
    c_column_indices(),
    tranpose_a_xadj(), tranpose_b_xadj(), tranpose_c_xadj(),
    tranpose_a_adj(), tranpose_b_adj(), tranpose_c_adj(),
    transpose_a(false),transpose_b(false), transpose_c_symbolic(false),
    num_colors(0),
    color_xadj(), color_adj(), vertex_colors(), num_multi_colors(0),num_used_colors(0),
    min_result_row_for_each_row(),

    create_lower_triangular(false),
    sort_lower_triangular(2),
    sort_option (-1),
    lower_triangular_permutation(),
    lower_triangular_matrix_rowmap(),
    lower_triangular_matrix_entries(),
    incidence_matrix_row_map(),
    incidence_matrix_entries(),compress_second_matrix(true),

    multi_color_scale(1), mkl_sort_option(7), calculate_read_write_cost(false),
	coloring_input_file(""),
	coloring_output_file(""), min_hash_size_scale(1), compression_cut_off(0.85), first_level_hash_cut_off(0.50),
	original_max_row_flops(std::numeric_limits<size_t>::max()), original_overall_flops(std::numeric_limits<size_t>::max()),
    persistent_a_xadj(), persistent_b_xadj(), persistent_a_adj(), persistent_b_adj(), MaxColDenseAcc(250001),
    mkl_keep_output(true),
    mkl_convert_to_1base(true), is_compression_single_step(false)
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
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

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
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

#if defined( KOKKOS_ENABLE_CUDA )
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


  void set_compression(bool compress_second_matrix_){
    this->compress_second_matrix = compress_second_matrix_;
  }

  bool get_compression (){
    return this->compress_second_matrix;
  }

  SPGEMMAccumulator get_accumulator_type() const {return this->accumulator_type;}
  void set_accumulator_type(const SPGEMMAccumulator &acc_type){this->accumulator_type = acc_type;}


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

#if defined( KOKKOS_ENABLE_CUDA )
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

  void set_compression_steps(bool isCompressionSingleStep){
    this->is_compression_single_step = isCompressionSingleStep;
  }

  void set_min_col_of_row(nnz_lno_persistent_work_view_t min_result_row_for_each_row_){
    this->min_result_row_for_each_row = min_result_row_for_each_row_;
  }

  nnz_lno_persistent_work_view_t get_min_col_of_row(){
    return this->min_result_row_for_each_row;
  }

  bool get_compression_step(){
    return is_compression_single_step;
  }
};


  inline SPGEMMAlgorithm StringToSPGEMMAlgorithm(std::string & name) {
    if(name=="SPGEMM_DEFAULT")             return SPGEMM_KK;
    else if(name=="SPGEMM_KK")       	   return SPGEMM_KK;
    else if(name=="SPGEMM_KK_MEMORY")      return SPGEMM_KK_MEMORY;
    else if(name=="SPGEMM_KK_DENSE")       return SPGEMM_KK_DENSE;
    else if(name=="SPGEMM_KK_LP")  		   return SPGEMM_KK_LP;
    else if(name=="SPGEMM_KK_MEMSPEED")    return SPGEMM_KK;

    else if(name=="SPGEMM_DEBUG")          return SPGEMM_SERIAL;
    else if(name=="SPGEMM_SERIAL")         return SPGEMM_SERIAL;
    else if(name=="SPGEMM_CUSPARSE")       return SPGEMM_CUSPARSE;
    else if(name=="SPGEMM_CUSP")           return SPGEMM_CUSP;
    else if(name=="SPGEMM_MKL")            return SPGEMM_MKL;
    else if(name=="SPGEMM_VIENNA")         return SPGEMM_VIENNA;
    else
      throw std::runtime_error("Invalid SPGEMMAlgorithm name");
  }



}

#endif
