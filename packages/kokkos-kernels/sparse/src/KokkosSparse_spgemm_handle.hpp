//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#ifndef _SPGEMMHANDLE_HPP
#define _SPGEMMHANDLE_HPP

#include <KokkosKernels_config.h>
#include <KokkosKernels_Controls.hpp>
#include <KokkosSparse_Utils.hpp>
#include <Kokkos_Core.hpp>
#include <iostream>
#include <string>
// #define VERBOSE

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
#include "KokkosSparse_Utils_rocsparse.hpp"
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "KokkosSparse_Utils_cusparse.hpp"
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
#include "KokkosSparse_Utils_mkl.hpp"
#endif

namespace KokkosSparse {

// TODO:SPGEMM_KK_MEMORY2 option is for testing in openmp.
// it wont work on cuda, not bind to a test program.
// hidden parameter for StringToSPGEMMAlgorithm for now.
enum SPGEMMAlgorithm {
  /*DEFAULT*/ SPGEMM_KK,
  SPGEMM_KK_DENSE,
  SPGEMM_KK_MEMORY,
  SPGEMM_KK_LP,  // KKVARIANTS
  SPGEMM_CUSPARSE [[deprecated("cuSPARSE is now used automatically in all "
                               "supported SpGEMM calls, if enabled.")]],
  SPGEMM_MKL [[deprecated("MKL is now used automatically in all supported "
                          "SpGEMM calls, if enabled.")]],
  SPGEMM_MKL2PHASE [[deprecated("MKL is now used automatically in all "
                                "supported SpGEMM calls, if enabled.")]],
  SPGEMM_ROCSPARSE [[deprecated("rocSPARSE is now used automatically in all "
                                "supported SpGEMM calls, if enabled.")]],

  // TRIANGLE COUNTING SPECIALIZED
  SPGEMM_KK_TRIANGLE_AI,        // SPGEMM_KK_TRIANGLE_DEFAULT, SPGEMM_KK_TRIANGLE_MEM,
                                // SPGEMM_KK_TRIANGLE_DENSE,
  SPGEMM_KK_TRIANGLE_IA_UNION,  // SPGEMM_KK_TRIANGLE_DEFAULT_IA_UNION,
                                // SPGEMM_KK_TRIANGLE_MEM_IA_UNION,
                                // SPGEMM_KK_TRIANGLE_DENSE_IA_UNION,
  SPGEMM_KK_TRIANGLE_IA,        // SPGEMM_KK_TRIANGLE_IA_DEFAULT,
                                // SPGEMM_KK_TRIANGLE_IA_MEM,
                                // SPGEMM_KK_TRIANGLE_IA_DENSE,
  SPGEMM_KK_TRIANGLE_LL,
  SPGEMM_KK_TRIANGLE_LU,

  // below research code.
  SPGEMM_KK_MULTIMEM,
  SPGEMM_KK_OUTERMULTIMEM,
  SPGEMM_DEFAULT,
  SPGEMM_DEBUG,
  SPGEMM_SERIAL,
  SPGEMM_KK_CUCKOO,            // USE CUCKOO HASHING
  SPGEMM_KK_TRACKED_CUCKOO,    // USE SCALED CUCKOO HASHING
  SPGEMM_KK_TRACKED_CUCKOO_F,  // USE SCALED FANCY CUCKOO HASHING
  SPGEMM_KK_SPEED,             // DENSE ACCUMULATOR SAME AS SPEED
  SPGEMM_KK_MEMORY_SORTED,
  SPGEMM_KK_MEMORY_TEAM,
  SPGEMM_KK_MEMORY_BIGTEAM,
  SPGEMM_KK_MEMORY_SPREADTEAM,
  SPGEMM_KK_MEMORY_BIGSPREADTEAM,
  SPGEMM_KK_MEMORY2,
  SPGEMM_KK_MEMSPEED
};

enum SPGEMMAccumulator {
  SPGEMM_ACC_DEFAULT,
  SPGEMM_ACC_DENSE,
  SPGEMM_ACC_SPARSE,
};
template <class size_type_, class lno_t_, class scalar_t_, class ExecutionSpace, class TemporaryMemorySpace,
          class PersistentMemorySpace>
class SPGEMMHandle {
 public:
  typedef ExecutionSpace HandleExecSpace;
  typedef TemporaryMemorySpace HandleTempMemorySpace;
  typedef PersistentMemorySpace HandlePersistentMemorySpace;

  typedef typename std::remove_const<size_type_>::type size_type;
  typedef const size_type const_size_type;

  typedef typename std::remove_const<lno_t_>::type nnz_lno_t;
  typedef const nnz_lno_t const_nnz_lno_t;

  typedef typename std::remove_const<scalar_t_>::type nnz_scalar_t;
  typedef const nnz_scalar_t const_nnz_scalar_t;

  typedef typename Kokkos::View<size_type *, HandleTempMemorySpace> row_lno_temp_work_view_t;
  typedef typename Kokkos::View<size_type *, HandlePersistentMemorySpace> row_lno_persistent_work_view_t;
  typedef typename row_lno_persistent_work_view_t::HostMirror row_lno_persistent_work_host_view_t;  // Host view type

  typedef typename Kokkos::View<nnz_scalar_t *, HandleTempMemorySpace> scalar_temp_work_view_t;
  typedef typename Kokkos::View<nnz_scalar_t *, HandlePersistentMemorySpace> scalar_persistent_work_view_t;

  typedef typename Kokkos::View<nnz_lno_t *, HandleTempMemorySpace> nnz_lno_temp_work_view_t;
  typedef typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace> nnz_lno_persistent_work_view_t;
  typedef typename nnz_lno_persistent_work_view_t::HostMirror nnz_lno_persistent_work_host_view_t;  // Host view type

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
  struct rocSparseSpgemmHandleType {
    KokkosKernels::Experimental::Controls kkControls;  // give a singleton rocsparse handle
    rocsparse_handle rocsparseHandle;
    rocsparse_operation opA, opB;
    rocsparse_mat_descr descr_A, descr_B, descr_C, descr_D;
    rocsparse_mat_info info_C;
    size_t bufferSize;
    void *buffer;

    rocSparseSpgemmHandleType(bool transposeA, bool transposeB) {
      opA = transposeA ? rocsparse_operation_transpose : rocsparse_operation_none;
      opB = transposeB ? rocsparse_operation_transpose : rocsparse_operation_none;

      bufferSize = 0;
      buffer     = nullptr;
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_mat_descr(&descr_A));
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_mat_descr(&descr_B));
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_mat_descr(&descr_C));
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_mat_descr(&descr_D));
      KOKKOS_ROCSPARSE_SAFE_CALL_IMPL(rocsparse_create_mat_info(&info_C));
      rocsparseHandle = kkControls.getRocsparseHandle();
    }

    ~rocSparseSpgemmHandleType() {
      rocsparse_destroy_mat_info(info_C);
      rocsparse_destroy_mat_descr(descr_A);
      rocsparse_destroy_mat_descr(descr_B);
      rocsparse_destroy_mat_descr(descr_C);
      rocsparse_destroy_mat_descr(descr_D);
    }
  };

#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#if (CUDA_VERSION >= 11000)
  struct cuSparseSpgemmHandleType {
    KokkosKernels::Experimental::Controls kkControls;
    cusparseHandle_t cusparseHandle;
    cusparseOperation_t opA, opB;
    cusparseSpMatDescr_t descr_A, descr_B, descr_C;
    cusparseSpGEMMDescr_t spgemmDescr;

    cudaDataType scalarType;
    cusparseSpGEMMAlg_t alg;

    // buffer1~2 are not needed in the numeric phase,
    // so we don't have them as member variables
    size_t bufferSize3, bufferSize4, bufferSize5;
    void *buffer3, *buffer4, *buffer5;

    cuSparseSpgemmHandleType(bool transposeA, bool transposeB) {
      opA        = transposeA ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE;
      opB        = transposeB ? CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE;
      scalarType = Impl::cuda_data_type_from<nnz_scalar_t>();

      alg         = CUSPARSE_SPGEMM_DEFAULT;
      bufferSize3 = bufferSize4 = bufferSize5 = 0;
      buffer3 = buffer4 = buffer5 = nullptr;

      cusparseHandle = kkControls.getCusparseHandle();
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_createDescr(&spgemmDescr));
    }

    ~cuSparseSpgemmHandleType() {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(descr_A));
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(descr_B));
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(descr_C));
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpGEMM_destroyDescr(spgemmDescr));
      KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(buffer3));
      KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(buffer4));
      KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(buffer5));
    }
  };
#else
  struct cuSparseSpgemmHandleType {
    cusparseHandle_t cusparseHandle;
    // Descriptor for any general matrix with index base 0
    cusparseMatDescr_t generalDescr;

    cuSparseSpgemmHandleType(bool /* transposeA */, bool /* transposeB */) {
      KokkosKernels::Experimental::Controls kkControls;
      // Get singleton cusparse handle from default controls
      cusparseHandle = kkControls.getCusparseHandle();

      KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreateMatDescr(&generalDescr));
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatType(generalDescr, CUSPARSE_MATRIX_TYPE_GENERAL));
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(generalDescr, CUSPARSE_INDEX_BASE_ZERO));
    }
    ~cuSparseSpgemmHandleType() { KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroyMatDescr(generalDescr)); }
  };
#endif
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
  struct mklSpgemmHandleType {
    // Allow mkl_sparse_sp2m (in SPARSE_STAGE_NNZ_COUNT mode) to construct C.
    // Then this assumes ownership of it and will destroy it later.
    mklSpgemmHandleType(const sparse_matrix_t &C_) : C(C_) {}

    ~mklSpgemmHandleType() { KOKKOSKERNELS_MKL_SAFE_CALL(mkl_sparse_destroy(C)); }

    sparse_matrix_t C;
  };
#endif

 private:
  SPGEMMAlgorithm algorithm_type;
  SPGEMMAccumulator accumulator_type;
  size_type result_nnz_size;

  bool called_symbolic;
  bool computed_rowptrs;
  bool computed_rowflops;
  bool computed_entries;
  bool called_numeric;

  int suggested_vector_size;
  int suggested_team_size;
  nnz_lno_t max_nnz_inresult;  // C max nonzeros per row.
  bool computed_max_nnz_inresult;
  nnz_lno_t max_nnz_compressed_result;

  size_type compressed_b_size;
  row_lno_temp_work_view_t compressed_b_rowmap;  // compressed_b_set_begins, compressed_b_set_nexts;
  nnz_lno_temp_work_view_t compressed_b_set_indices, compressed_b_sets;

  row_lno_temp_work_view_t compressed_c_rowmap;

  nnz_lno_temp_work_view_t c_column_indices;

  row_lno_temp_work_view_t tranpose_a_xadj, tranpose_b_xadj, tranpose_c_xadj;
  nnz_lno_temp_work_view_t tranpose_a_adj, tranpose_b_adj, tranpose_c_adj;

  bool transpose_a, transpose_b, transpose_c_symbolic;

  nnz_lno_t num_colors;
  nnz_lno_persistent_work_host_view_t color_xadj;
  nnz_lno_persistent_work_view_t color_adj, vertex_colors;
  nnz_lno_t num_multi_colors, num_used_colors;
  nnz_lno_persistent_work_view_t min_result_row_for_each_row;

  bool create_lower_triangular;
  int sort_lower_triangular;  // 0 - do not sort // 1 - sort // 2 - Algorithm
                              // decides (default)
  nnz_lno_persistent_work_view_t lower_triangular_permutation;

  row_lno_persistent_work_view_t lower_triangular_matrix_rowmap;
  nnz_lno_persistent_work_view_t lower_triangular_matrix_entries;

  row_lno_persistent_work_view_t incidence_matrix_row_map;
  nnz_lno_persistent_work_view_t incidence_matrix_entries;
  bool compress_second_matrix;

  double multi_color_scale;
  int mkl_sort_option;
  bool calculate_read_write_cost;

 public:
  std::string coloring_input_file;
  std::string coloring_output_file;

  int min_hash_size_scale;
  double compression_cut_off;
  double first_level_hash_cut_off;
  size_t original_max_row_flops, original_overall_flops;
  row_lno_persistent_work_view_t row_flops;

  size_t compressed_max_row_flops, compressed_overall_flops;

  void set_first_level_hash_cut_off(double first_level_hash_cut_off_) {
    this->first_level_hash_cut_off = first_level_hash_cut_off_;
  }

  double get_first_level_hash_cut_off() { return this->first_level_hash_cut_off; }

  void set_compression_cut_off(double compression_cut_off_) { this->compression_cut_off = compression_cut_off_; }

  double get_compression_cut_off() { return this->compression_cut_off; }
  void set_min_hash_size_scale(int scale) { min_hash_size_scale = scale; }
  int get_min_hash_size_scale() { return min_hash_size_scale; }
  void set_read_write_cost_calc(bool read_write_cost_cal) { this->calculate_read_write_cost = read_write_cost_cal; }
  int get_read_write_cost_calc() { return this->calculate_read_write_cost; }

  typename Kokkos::View<int *, HandlePersistentMemorySpace> persistent_c_xadj, persistent_a_xadj, persistent_b_xadj,
      persistent_a_adj, persistent_b_adj;
  size_t MaxColDenseAcc;
  bool mkl_keep_output;
  bool mkl_convert_to_1base;
  bool is_compression_single_step;

  void set_mkl_sort_option(int mkl_sort_option_) { this->mkl_sort_option = mkl_sort_option_; }
  int get_mkl_sort_option() { return this->mkl_sort_option; }

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
 private:
  rocSparseSpgemmHandleType *rocsparse_spgemm_handle;

 public:
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
 private:
  cuSparseSpgemmHandleType *cusparse_spgemm_handle;

 public:
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
 private:
  mklSpgemmHandleType *mkl_spgemm_handle;

 public:
#endif

  void set_c_column_indices(nnz_lno_temp_work_view_t c_col_indices_) { this->c_column_indices = c_col_indices_; }

  nnz_lno_temp_work_view_t get_c_column_indices() { return this->c_column_indices; }

  void set_color_xadj(nnz_lno_t num_colors_, nnz_lno_persistent_work_host_view_t color_xadj_,
                      nnz_lno_persistent_work_view_t color_adj_, nnz_lno_persistent_work_view_t vertex_colors_,
                      nnz_lno_t num_multi_colors_, nnz_lno_t num_used_colors_) {
    num_colors    = num_colors_;
    color_xadj    = color_xadj_;
    color_adj     = color_adj_;
    vertex_colors = vertex_colors_;

    num_multi_colors = num_multi_colors_;
    num_used_colors  = num_used_colors_;
  }

  /// \brief sets the result nnz size.
  /// \param result_nnz_size_ size of the output matrix.
  void set_c_nnz(size_type result_nnz_size_) { this->result_nnz_size = result_nnz_size_; }
  /**
   * \brief returns the result nnz size.
   */
  size_type get_c_nnz() { return this->result_nnz_size; }

  void set_multi_color_scale(double multi_color_scale_) { this->multi_color_scale = multi_color_scale_; }

  double get_multi_color_scale() { return this->multi_color_scale; }

  void get_color_xadj(nnz_lno_t &num_colors_, nnz_lno_persistent_work_host_view_t &color_xadj_,
                      nnz_lno_persistent_work_view_t &color_adj_, nnz_lno_persistent_work_view_t &vertex_colors_,
                      nnz_lno_t &num_multi_colors_, nnz_lno_t &num_used_colors_) {
    num_colors_       = num_colors;
    color_xadj_       = color_xadj;
    color_adj_        = color_adj;
    num_multi_colors_ = num_multi_colors;
    num_used_colors_  = num_used_colors;
    vertex_colors_    = vertex_colors;
  }

  void set_compressed_c(row_lno_temp_work_view_t compressed_c_rowmap_) { compressed_c_rowmap = compressed_c_rowmap_; }

  void get_compressed_c(row_lno_temp_work_view_t &compressed_c_rowmap_) { compressed_c_rowmap_ = compressed_c_rowmap; }

  void set_sort_lower_triangular(int option) { this->sort_lower_triangular = option; }
  int get_sort_lower_triangular() { return this->sort_lower_triangular; }

  void set_create_lower_triangular(bool option) { this->create_lower_triangular = option; }
  bool get_create_lower_triangular() { return this->create_lower_triangular; }

  void set_lower_triangular_permutation(nnz_lno_persistent_work_view_t ltp_) {
    this->lower_triangular_permutation = ltp_;
  }

  nnz_lno_persistent_work_view_t get_lower_triangular_permutation() { return this->lower_triangular_permutation; }

  void set_lower_triangular_matrix(row_lno_persistent_work_view_t lower_triangular_matrix_rowmap_,
                                   nnz_lno_persistent_work_view_t lower_triangular_matrix_entries_) {
    this->lower_triangular_matrix_rowmap  = lower_triangular_matrix_rowmap_;
    this->lower_triangular_matrix_entries = lower_triangular_matrix_entries_;
  }
  void get_lower_triangular_matrix(row_lno_persistent_work_view_t &lower_triangular_matrix_rowmap_,
                                   nnz_lno_persistent_work_view_t &lower_triangular_matrix_entries_) {
    lower_triangular_matrix_rowmap_  = this->lower_triangular_matrix_rowmap;
    lower_triangular_matrix_entries_ = this->lower_triangular_matrix_entries;
  }

  void set_compressed_b(size_type b_nnz_size, row_lno_temp_work_view_t compressed_b_rowmap_,
                        nnz_lno_temp_work_view_t compressed_b_set_indices_,
                        nnz_lno_temp_work_view_t compressed_b_sets_) {
    compressed_b_size        = b_nnz_size;
    compressed_b_rowmap      = compressed_b_rowmap_;
    compressed_b_set_indices = compressed_b_set_indices_;
    compressed_b_sets        = compressed_b_sets_;
  }

  void get_compressed_b(size_type &b_nnz_size, row_lno_temp_work_view_t &compressed_b_rowmap_,
                        nnz_lno_temp_work_view_t &compressed_b_set_indices_,
                        nnz_lno_temp_work_view_t &compressed_b_sets_) {
    b_nnz_size                = compressed_b_size;
    compressed_b_rowmap_      = compressed_b_rowmap;
    compressed_b_set_indices_ = compressed_b_set_indices;
    compressed_b_sets_        = compressed_b_sets;
  }

  /**
   * \brief Default constructor.
   */
  SPGEMMHandle(SPGEMMAlgorithm gs = SPGEMM_DEFAULT)
      : algorithm_type(gs),
        accumulator_type(SPGEMM_ACC_DEFAULT),
        result_nnz_size(0),
        called_symbolic(false),
        computed_rowptrs(false),
        computed_rowflops(false),
        computed_entries(false),
        called_numeric(false),
        suggested_vector_size(0),
        suggested_team_size(0),
        max_nnz_inresult(0),
        computed_max_nnz_inresult(false),
        c_column_indices(),
        tranpose_a_xadj(),
        tranpose_b_xadj(),
        tranpose_c_xadj(),
        tranpose_a_adj(),
        tranpose_b_adj(),
        tranpose_c_adj(),
        transpose_a(false),
        transpose_b(false),
        transpose_c_symbolic(false),
        num_colors(0),
        color_xadj(),
        color_adj(),
        vertex_colors(),
        num_multi_colors(0),
        num_used_colors(0),
        min_result_row_for_each_row(),

        create_lower_triangular(false),
        sort_lower_triangular(2),
        lower_triangular_permutation(),
        lower_triangular_matrix_rowmap(),
        lower_triangular_matrix_entries(),
        incidence_matrix_row_map(),
        incidence_matrix_entries(),
        compress_second_matrix(true),

        multi_color_scale(1),
        mkl_sort_option(7),
        calculate_read_write_cost(false),
        coloring_input_file(""),
        coloring_output_file(""),
        min_hash_size_scale(1),
        compression_cut_off(0.85),
        first_level_hash_cut_off(0.50),
        original_max_row_flops(std::numeric_limits<size_t>::max()),
        original_overall_flops(std::numeric_limits<size_t>::max()),
        persistent_a_xadj(),
        persistent_b_xadj(),
        persistent_a_adj(),
        persistent_b_adj(),
        MaxColDenseAcc(250001),
        mkl_keep_output(true),
        mkl_convert_to_1base(true),
        is_compression_single_step(false)
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
        ,
        rocsparse_spgemm_handle(nullptr)
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
        ,
        cusparse_spgemm_handle(nullptr)
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
        ,
        mkl_spgemm_handle(nullptr)
#endif
  {
    if (gs == SPGEMM_DEFAULT) {
      this->choose_default_algorithm();
    }
  }

  virtual ~SPGEMMHandle() {
#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
    this->destroy_rocsparse_spgemm_handle();
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    this->destroy_cusparse_spgemm_handle();
#endif
#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
    this->destroy_mkl_spgemm_handle();
#endif
  };

#ifdef KOKKOSKERNELS_ENABLE_TPL_ROCSPARSE
  void create_rocsparse_spgemm_handle(bool transA, bool transB) {
    this->destroy_rocsparse_spgemm_handle();
    this->rocsparse_spgemm_handle = new rocSparseSpgemmHandleType(transA, transB);
  }

  rocSparseSpgemmHandleType *get_rocsparse_spgemm_handle() { return this->rocsparse_spgemm_handle; }

  void destroy_rocsparse_spgemm_handle() {
    if (this->rocsparse_spgemm_handle != nullptr) {
      delete this->rocsparse_spgemm_handle;
      this->rocsparse_spgemm_handle = nullptr;
    }
  }

#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  void create_cusparse_spgemm_handle(bool transA, bool transB) {
    this->destroy_cusparse_spgemm_handle();
    this->cusparse_spgemm_handle = new cuSparseSpgemmHandleType(transA, transB);
  }
  void destroy_cusparse_spgemm_handle() {
    if (this->cusparse_spgemm_handle != nullptr) {
      delete this->cusparse_spgemm_handle;
      this->cusparse_spgemm_handle = nullptr;
    }
  }

  cuSparseSpgemmHandleType *get_cusparse_spgemm_handle() { return this->cusparse_spgemm_handle; }
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_MKL
  void create_mkl_spgemm_handle(sparse_matrix_t C) {
    this->destroy_mkl_spgemm_handle();
    this->mkl_spgemm_handle = new mklSpgemmHandleType(C);
  }
  void destroy_mkl_spgemm_handle() {
    if (this->mkl_spgemm_handle != nullptr) {
      delete this->mkl_spgemm_handle;
      this->mkl_spgemm_handle = nullptr;
    }
  }

  mklSpgemmHandleType *get_mkl_spgemm_handle() { return this->mkl_spgemm_handle; }
#endif

  void choose_default_algorithm() {
#if defined(KOKKOS_ENABLE_SERIAL)
    if (std::is_same<Kokkos::Serial, ExecutionSpace>::value) {
      this->algorithm_type = SPGEMM_SERIAL;
#ifdef VERBOSE
      std::cout << "Serial Execution Space, Default Algorithm: SPGEMM_SERIAL" << std::endl;
#endif
    }
#endif

#if defined(KOKKOS_ENABLE_THREADS)
    if (std::is_same<Kokkos::Threads, ExecutionSpace>::value) {
      this->algorithm_type = SPGEMM_SERIAL;
#ifdef VERBOSE
      std::cout << "THREADS Execution Space, Default Algorithm: SPGEMM_SERIAL" << std::endl;
#endif
    }
#endif

#if defined(KOKKOS_ENABLE_OPENMP)
    if (std::is_same<Kokkos::OpenMP, ExecutionSpace>::value) {
      this->algorithm_type = SPGEMM_SERIAL;
#ifdef VERBOSE
      std::cout << "OpenMP Execution Space, Default Algorithm: SPGEMM_SERIAL" << std::endl;
#endif
    }
#endif

#if defined(KOKKOS_ENABLE_CUDA)
    if (std::is_same<Kokkos::Cuda, ExecutionSpace>::value) {
      this->algorithm_type = SPGEMM_KK;
#ifdef VERBOSE
      std::cout << "Cuda Execution Space, Default Algorithm: SPGEMM_KK" << std::endl;
#endif
    }
#endif

#if defined(KOKKOS_ENABLE_HIP)
    if (std::is_same<Kokkos::HIP, ExecutionSpace>::value) {
      this->algorithm_type = SPGEMM_KK;
#ifdef VERBOSE
      std::cout << "HIP Execution Space, Default Algorithm: SPGEMM_KK" << std::endl;
#endif
    }
#endif

#if defined(KOKKOS_ENABLE_SYCL)
    if (std::is_same<Kokkos::Experimental::SYCL, ExecutionSpace>::value) {
      this->algorithm_type = SPGEMM_KK;
#ifdef VERBOSE
      std::cout << "SYCL Execution Space, Default Algorithm: SPGEMM_KK" << std::endl;
#endif
    }
#endif
  }

  void set_compression(bool compress_second_matrix_) { this->compress_second_matrix = compress_second_matrix_; }

  bool get_compression() { return this->compress_second_matrix; }

  SPGEMMAccumulator get_accumulator_type() const { return this->accumulator_type; }
  void set_accumulator_type(const SPGEMMAccumulator &acc_type) { this->accumulator_type = acc_type; }

  // getters
  SPGEMMAlgorithm get_algorithm_type() const { return this->algorithm_type; }

  bool is_symbolic_called() { return this->called_symbolic; }
  bool are_rowptrs_computed() { return this->computed_rowptrs; }
  bool are_rowflops_computed() { return this->computed_rowflops; }
  bool are_entries_computed() { return this->computed_entries; }
  bool is_numeric_called() { return this->called_numeric; }

  template <typename c_row_view_t>
  nnz_lno_t get_max_result_nnz(const c_row_view_t &row_mapC) {
    if (!this->computed_max_nnz_inresult) {
      this->max_nnz_inresult = KokkosSparse::Impl::graph_max_degree<HandleExecSpace, size_type, c_row_view_t>(row_mapC);
      this->computed_max_nnz_inresult = true;
    }
    return this->max_nnz_inresult;
  }

  nnz_lno_t get_max_compresed_result_nnz() const { return this->max_nnz_compressed_result; }

  // setters
  void set_algorithm_type(const SPGEMMAlgorithm &sgs_algo) { this->algorithm_type = sgs_algo; }
  void set_call_symbolic(bool call = true) { this->called_symbolic = call; }
  void set_computed_rowptrs() { this->computed_rowptrs = true; }
  void set_computed_rowflops() { this->computed_rowflops = true; }
  void set_computed_entries() { this->computed_entries = true; }
  void set_call_numeric(bool call = true) { this->called_numeric = call; }

  void set_max_result_nnz(nnz_lno_t nz) {
    this->max_nnz_inresult          = nz;
    this->computed_max_nnz_inresult = true;
  }

  void set_max_compresed_result_nnz(nnz_lno_t num_result_nnz_) { this->max_nnz_compressed_result = num_result_nnz_; }

  void vector_team_size(int max_allowed_team_size, int &suggested_vector_size_, int &suggested_team_size_, size_type nr,
                        size_type nnz) {
    // suggested_team_size_ =  this->suggested_team_size = 1;
    // suggested_vector_size_=this->suggested_vector_size = 1;
    // return;
    if (this->suggested_team_size && this->suggested_vector_size) {
      // already set in the handle
      suggested_vector_size_ = this->suggested_vector_size;
      suggested_team_size_   = this->suggested_team_size;
      return;
    }

    // otherwise, recompute team_size/vector_size based on heuristic and save
    // them in the handle
    suggested_vector_size_ = KokkosKernels::Impl::kk_get_suggested_vector_size(
        nr, nnz, KokkosKernels::Impl::kk_get_exec_space_type<ExecutionSpace>());
    if (KokkosKernels::Impl::kk_is_gpu_exec_space<ExecutionSpace>())
      suggested_team_size_ = max_allowed_team_size / suggested_vector_size_;
    else
      suggested_team_size = max_allowed_team_size;
    this->suggested_vector_size = suggested_vector_size_;
    this->suggested_team_size   = suggested_vector_size_;
  }

  void set_compression_steps(bool isCompressionSingleStep) {
    this->is_compression_single_step = isCompressionSingleStep;
  }

  void set_min_col_of_row(nnz_lno_persistent_work_view_t min_result_row_for_each_row_) {
    this->min_result_row_for_each_row = min_result_row_for_each_row_;
  }

  nnz_lno_persistent_work_view_t get_min_col_of_row() { return this->min_result_row_for_each_row; }

  bool get_compression_step() { return is_compression_single_step; }

 private:
  // An SpGEMM handle can be reused for multiple products C = A*B, but only if
  // the sparsity patterns of A and B do not change. Enforce this (in debug
  // builds only) by recording hashes of the graphs, and then checking they
  // match in later calls.
  bool computedInputHashes = false;
  uint32_t a_graph_hash    = 0U;
  uint32_t b_graph_hash    = 0U;

 public:
  template <typename a_rowptrs_t, typename a_entries_t, typename b_rowptrs_t, typename b_entries_t>
  bool checkMatrixIdentitiesSymbolic(const a_rowptrs_t &a_rowptrsIn, const a_entries_t &a_entriesIn,
                                     const b_rowptrs_t &b_rowptrsIn, const b_entries_t &b_entriesIn) {
#ifndef NDEBUG
    // If this is the first symbolic call, assign the handle's CRS pointers to
    // check against later
    if (!computedInputHashes) {
      a_graph_hash        = KokkosKernels::Impl::hashView(a_rowptrsIn) ^ KokkosKernels::Impl::hashView(a_entriesIn);
      b_graph_hash        = KokkosKernels::Impl::hashView(b_rowptrsIn) ^ KokkosKernels::Impl::hashView(b_entriesIn);
      computedInputHashes = true;
    } else {
      if (a_graph_hash != (KokkosKernels::Impl::hashView(a_rowptrsIn) ^ KokkosKernels::Impl::hashView(a_entriesIn)))
        return false;
      if (b_graph_hash != (KokkosKernels::Impl::hashView(b_rowptrsIn) ^ KokkosKernels::Impl::hashView(b_entriesIn)))
        return false;
    }
#else
    (void)a_rowptrsIn;
    (void)a_entriesIn;
    (void)b_rowptrsIn;
    (void)b_entriesIn;
#endif
    return true;
  }

  template <typename a_rowptrs_t, typename a_entries_t, typename b_rowptrs_t, typename b_entries_t>
  bool checkMatrixIdentitiesNumeric(const a_rowptrs_t &a_rowptrsIn, const a_entries_t &a_entriesIn,
                                    const b_rowptrs_t &b_rowptrsIn, const b_entries_t &b_entriesIn) {
#ifndef NDEBUG
    if (a_graph_hash != (KokkosKernels::Impl::hashView(a_rowptrsIn) ^ KokkosKernels::Impl::hashView(a_entriesIn)))
      return false;
    if (b_graph_hash != (KokkosKernels::Impl::hashView(b_rowptrsIn) ^ KokkosKernels::Impl::hashView(b_entriesIn)))
      return false;
#else
    (void)a_rowptrsIn;
    (void)a_entriesIn;
    (void)b_rowptrsIn;
    (void)b_entriesIn;
#endif
    return true;
  }
};

inline SPGEMMAlgorithm StringToSPGEMMAlgorithm(std::string &name) {
  if (name == "SPGEMM_DEFAULT")
    return SPGEMM_KK;
  else if (name == "SPGEMM_KK")
    return SPGEMM_KK;
  else if (name == "SPGEMM_KK_MEMORY")
    return SPGEMM_KK_MEMORY;
  else if (name == "SPGEMM_KK_DENSE")
    return SPGEMM_KK_DENSE;
  else if (name == "SPGEMM_KK_LP")
    return SPGEMM_KK_LP;
  else if (name == "SPGEMM_KK_MEMSPEED")
    return SPGEMM_KK;
  else if (name == "SPGEMM_DEBUG")
    return SPGEMM_SERIAL;
  else if (name == "SPGEMM_SERIAL")
    return SPGEMM_SERIAL;
  else if (name == "SPGEMM_CUSPARSE")
    throw std::runtime_error(
        "Enum value SPGEMM_CUSPARSE is deprecated. cuSPARSE is automatically "
        "used in all supported SpGEMM calls.");
  else if (name == "SPGEMM_MKL")
    throw std::runtime_error(
        "Enum value SPGEMM_MKL is deprecated. MKL is automatically used in all "
        "supported SpGEMM calls.");
  else if (name == "SPGEMM_MKL2PHASE")
    throw std::runtime_error(
        "Enum value SPGEMM_MKL2PHASE is deprecated. MKL is automatically used "
        "in all supported SpGEMM calls.");
  else if (name == "SPGEMM_ROCSPARSE")
    throw std::runtime_error(
        "Enum value SPGEMM_ROCSPARSE is deprecated. rocSPARSE is automatically "
        "used in all supported SpGEMM calls.");
  else
    throw std::runtime_error("Invalid SPGEMMAlgorithm name");
}

}  // namespace KokkosSparse

#endif
