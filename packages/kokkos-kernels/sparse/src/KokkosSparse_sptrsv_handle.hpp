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

#include <Kokkos_Core.hpp>
#include <iostream>
#include <string>

#ifndef KOKKOSSPARSE_SPTRSVHANDLE_HPP
#define KOKKOSSPARSE_SPTRSVHANDLE_HPP

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#include "KokkosSparse_Utils_cusparse.hpp"
#endif

#if defined(KOKKOS_ENABLE_CUDA) && 10000 < CUDA_VERSION && defined(KOKKOSKERNELS_ENABLE_EXP_CUDAGRAPH)
#define KOKKOSKERNELS_SPTRSV_CUDAGRAPHSUPPORT
#endif

// Enable supernodal sptrsv
#ifdef KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV
#include <KokkosSparse_CrsMatrix.hpp>
#endif

namespace KokkosSparse {
namespace Experimental {

// TODO TP2 algorithm had issues with some offset-ordinal combo to be addressed
// when compiled in Trilinos...
enum class SPTRSVAlgorithm {
  SEQLVLSCHD_RP,
  SEQLVLSCHD_TP1 /*, SEQLVLSCHED_TP2*/,
  SEQLVLSCHD_TP1CHAIN,
  SPTRSV_CUSPARSE,
  SUPERNODAL_NAIVE,
  SUPERNODAL_ETREE,
  SUPERNODAL_DAG,
  SUPERNODAL_SPMV,
  SUPERNODAL_SPMV_DAG
};

template <class size_type_, class lno_t_, class scalar_t_, class ExecutionSpace, class TemporaryMemorySpace,
          class PersistentMemorySpace>
class SPTRSVHandle {
 public:
  using HandleExecSpace             = ExecutionSpace;
  using HandleTempMemorySpace       = TemporaryMemorySpace;
  using HandlePersistentMemorySpace = PersistentMemorySpace;

  using execution_space = ExecutionSpace;
  using memory_space    = HandlePersistentMemorySpace;

  using TeamPolicy  = Kokkos::TeamPolicy<execution_space>;
  using RangePolicy = Kokkos::RangePolicy<execution_space>;

  using size_type       = typename std::remove_const<size_type_>::type;
  using const_size_type = const size_type;

  using nnz_lno_t       = typename std::remove_const<lno_t_>::type;
  using const_nnz_lno_t = const nnz_lno_t;

  using scalar_t           = typename std::remove_const<scalar_t_>::type;
  using const_nnz_scalar_t = const scalar_t;

  // Row_map type (managed memory)
  using nnz_row_view_temp_t = typename Kokkos::View<size_type *, HandleTempMemorySpace>;
  using nnz_row_view_t      = typename Kokkos::View<size_type *, HandlePersistentMemorySpace>;
  using host_nnz_row_view_t = typename nnz_row_view_t::HostMirror;
  using int_row_view_t      = typename Kokkos::View<int *, HandlePersistentMemorySpace>;
  using int64_row_view_t    = typename Kokkos::View<int64_t *, HandlePersistentMemorySpace>;
  // typedef typename row_lno_persistent_work_view_t::HostMirror
  // row_lno_persistent_work_host_view_t; //Host view type
  using nnz_row_unmanaged_view_t =
      typename Kokkos::View<const size_type *, HandlePersistentMemorySpace,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;  // for rank1 subviews

  // values type (managed memory)
  using nnz_scalar_view_temp_t = typename Kokkos::View<scalar_t *, HandleTempMemorySpace>;
  using nnz_scalar_view_t      = typename Kokkos::View<scalar_t *, HandlePersistentMemorySpace>;
  using host_nnz_scalar_view_t = typename nnz_scalar_view_t::HostMirror;
  using nnz_scalar_unmanaged_view_t =
      typename Kokkos::View<const scalar_t *, HandlePersistentMemorySpace,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;  // for rank1 subviews

  // entries type (managed memory)
  using nnz_lno_view_temp_t      = typename Kokkos::View<nnz_lno_t *, HandleTempMemorySpace>;
  using nnz_lno_view_t           = typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace>;
  using hostspace_nnz_lno_view_t = typename Kokkos::View<nnz_lno_t *, Kokkos::HostSpace>;
  using host_nnz_lno_view_t      = typename nnz_lno_view_t::HostMirror;
  using nnz_lno_unmanaged_view_t =
      typename Kokkos::View<const nnz_lno_t *, HandlePersistentMemorySpace,
                            Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::RandomAccess>>;  // for rank1 subviews
  // typedef typename nnz_lno_persistent_work_view_t::HostMirror
  // nnz_lno_persistent_work_host_view_t; //Host view type

  using signed_integral_t = typename std::make_signed<typename nnz_row_view_t::non_const_value_type>::type;
  using signed_nnz_lno_view_t =
      Kokkos::View<signed_integral_t *, typename nnz_row_view_t::array_layout, typename nnz_row_view_t::device_type,
                   typename nnz_row_view_t::memory_traits>;

  using host_signed_nnz_lno_view_t = typename signed_nnz_lno_view_t::HostMirror;

  using mtx_scalar_view_t = typename Kokkos::View<scalar_t **, HandlePersistentMemorySpace>;

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
#if (CUDA_VERSION >= 11030)
  struct cuSparseHandleType {
    cusparseHandle_t handle;
    cusparseOperation_t transpose;
    cusparseSpMatDescr_t matDescr;
    cusparseDnVecDescr_t vecBDescr, vecBDescr_dummy;
    cusparseDnVecDescr_t vecXDescr, vecXDescr_dummy;
    cusparseSpSVDescr_t spsvDescr;
    void *pBuffer{nullptr};

    cuSparseHandleType(bool transpose_, bool /*is_lower*/) {
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseCreate(&handle));

      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST));

      if (transpose_) {
        transpose = CUSPARSE_OPERATION_TRANSPOSE;
      } else {
        transpose = CUSPARSE_OPERATION_NON_TRANSPOSE;
      }

      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpSV_createDescr(&spsvDescr));
    }

    ~cuSparseHandleType() {
      if (pBuffer != nullptr) {
        KOKKOS_IMPL_CUDA_SAFE_CALL(cudaFree(pBuffer));
        pBuffer = nullptr;
      }
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroySpMat(matDescr));
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseSpSV_destroyDescr(spsvDescr));
      KOKKOS_CUSPARSE_SAFE_CALL(cusparseDestroy(handle));
    }
  };
#else  // CUDA_VERSION < 11030
  struct cuSparseHandleType {
    cusparseHandle_t handle;
    cusparseOperation_t transpose;
    csrsv2Info_t info{0};
    cusparseMatDescr_t descr;
    cusparseSolvePolicy_t policy;
    void *pBuffer{nullptr};

    cuSparseHandleType(bool transpose_, bool is_lower) {
      cusparseStatus_t status;
      status = cusparseCreate(&handle);
      if (status != CUSPARSE_STATUS_SUCCESS) {
        throw std::runtime_error("cusparseCreate ERROR\n");
      }
      cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST);

      if (transpose_) {
        transpose = CUSPARSE_OPERATION_TRANSPOSE;
      } else {
        transpose = CUSPARSE_OPERATION_NON_TRANSPOSE;
      }

      status = cusparseCreateMatDescr(&descr);
      if (status != CUSPARSE_STATUS_SUCCESS) {
        throw std::runtime_error("cusparseCreateMatDescr descr ERROR\n");
      }
      cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
      cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

      if (is_lower)
        cusparseSetMatFillMode(descr, CUSPARSE_FILL_MODE_LOWER);
      else
        cusparseSetMatFillMode(descr, CUSPARSE_FILL_MODE_UPPER);

      cusparseSetMatDiagType(descr, CUSPARSE_DIAG_TYPE_NON_UNIT);

      policy = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
    }

    ~cuSparseHandleType() {
      if (pBuffer != nullptr) {
        cudaFree(pBuffer);
        pBuffer = nullptr;
      }
      cusparseDestroyMatDescr(descr);
      cusparseDestroy(handle);
    }
  };
#endif

  using SPTRSVcuSparseHandleType = cuSparseHandleType;
#endif

#ifdef KOKKOSKERNELS_SPTRSV_CUDAGRAPHSUPPORT
  bool cudagraphCreated;  // Move this later
  struct cudaGraphWrapperType {
    cudaGraph_t cudagraph;
    cudaGraphExec_t cudagraphinstance;
    cudaStream_t stream;

    // cudaGraphWrapperType() { }
    //~cudaGraphWrapperType() { }
  };

  using SPTRSVcudaGraphWrapperType = cudaGraphWrapperType;

  void create_SPTRSVcudaGraphWrapperType() {
    destroy_SPTRSVcudaGraphWrapperType();
    sptrsvCudaGraph = new SPTRSVcudaGraphWrapperType;
    cudaStreamCreate(&sptrsvCudaGraph->stream);
  }

  void destroy_SPTRSVcudaGraphWrapperType() {
    if (sptrsvCudaGraph != nullptr) {
      // cudaGraphExecDestroy(sptrsvCudaGraph->cudagraphinstance);
      // cudaGraphDestroy(sptrsvCudaGraph->cudagraph);
      cudaStreamDestroy(sptrsvCudaGraph->stream);
      delete sptrsvCudaGraph;
      sptrsvCudaGraph = nullptr;
    }
  }

  SPTRSVcudaGraphWrapperType *get_sptrsvCudaGraph() { return sptrsvCudaGraph; }
#endif

#ifdef KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV
  using supercols_memory_space = TemporaryMemorySpace;

  using supercols_host_execution_space = Kokkos::DefaultHostExecutionSpace;
  using supercols_host_memory_space    = typename supercols_host_execution_space::memory_space;

  using integer_view_t      = Kokkos::View<int *, supercols_memory_space>;
  using integer_view_host_t = Kokkos::View<int *, supercols_host_memory_space>;

  using workspace_t = typename Kokkos::View<scalar_t *, Kokkos::Device<execution_space, supercols_memory_space>>;

  //
  using host_crsmat_t = KokkosSparse::CrsMatrix<scalar_t, nnz_lno_t, supercols_host_execution_space, void, size_type>;
  using crsmat_t = KokkosSparse::CrsMatrix<scalar_t, nnz_lno_t, Kokkos::Device<execution_space, PersistentMemorySpace>,
                                           void, size_type>;

  //
  using host_graph_t = typename host_crsmat_t::StaticCrsGraphType;
  using graph_t      = typename crsmat_t::StaticCrsGraphType;

  //
  using crsmat_list_t = typename std::vector<crsmat_t>;
#endif

 private:
#ifdef KOKKOSKERNELS_SPTRSV_CUDAGRAPHSUPPORT
  SPTRSVcudaGraphWrapperType *sptrsvCudaGraph;
#endif

  size_type nrows;

  bool lower_tri;

  SPTRSVAlgorithm algm;

  // Symbolic: Level scheduling data
  signed_nnz_lno_view_t level_list;
  nnz_lno_view_t nodes_per_level;
  hostspace_nnz_lno_view_t hnodes_per_level;  // NEW
  nnz_lno_view_t nodes_grouped_by_level;
  hostspace_nnz_lno_view_t hnodes_grouped_by_level;  // NEW
  size_type nlevel;
  size_type block_size;  // block_size > 0 implies BSR

  int team_size;
  int vector_size;

  bool stored_diagonal;
  nnz_lno_view_t diagonal_offsets;
  nnz_scalar_view_t diagonal_values;  // inserted by rowid

  host_nnz_lno_view_t hdiagonal_offsets;
  host_nnz_scalar_view_t hdiagonal_values;  // inserted by rowid

  // Symbolic: Single-block chain data
  host_signed_nnz_lno_view_t h_chain_ptr;
  size_type num_chain_entries;
  signed_integral_t chain_threshold;

  bool symbolic_complete;
  bool numeric_complete;
  bool require_symbolic_lvlsched_phase;
  bool require_symbolic_chain_phase;

  void set_if_algm_require_symb_lvlsched() {
    if (algm == SPTRSVAlgorithm::SEQLVLSCHD_RP ||
        algm == SPTRSVAlgorithm::SEQLVLSCHD_TP1
        /*|| algm == SPTRSVAlgorithm::SEQLVLSCHED_TP2*/
        || algm == SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN
#ifdef KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV
        || algm == SPTRSVAlgorithm::SUPERNODAL_NAIVE || algm == SPTRSVAlgorithm::SUPERNODAL_ETREE ||
        algm == SPTRSVAlgorithm::SUPERNODAL_DAG || algm == SPTRSVAlgorithm::SUPERNODAL_SPMV ||
        algm == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG
#endif
    ) {
      require_symbolic_lvlsched_phase = true;
    } else {
      require_symbolic_lvlsched_phase = false;
    }
  }

  void set_if_algm_require_symb_chain() {
    if (algm == KokkosSparse::Experimental::SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN) {
      require_symbolic_chain_phase = true;
    } else {
      require_symbolic_chain_phase = false;
    }
  }

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  SPTRSVcuSparseHandleType *cuSPARSEHandle;
  int_row_view_t tmp_int_rowmap;
  int64_row_view_t tmp_int64_rowmap;
#endif

#ifdef KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV
  // specify if unit diagonal
  bool unit_diag;

  // stored either in CSR or CSC
  bool col_major;

  // number of supernodal columns
  signed_integral_t nsuper;

  // etree, parent's id for each supernodal column
  integer_view_host_t etree_host;

  // dag
  host_graph_t dag_host;

  // map from supernode to column id, i.e., superdols[s] = the first column id
  // of s-th supernode
  integer_view_host_t supercols_host;  // on the host
  integer_view_t supercols;            // on the default host/device

  // workspace size
  signed_integral_t lwork;
  workspace_t work;
  // offset to workspace for each supernodal column
  integer_view_host_t work_offset_host;
  integer_view_t work_offset;

  //
  bool merge_supernodes;
  bool invert_diagonal;
  bool invert_offdiagonal;
  int *etree;

  // type of kernels used at each level
  bool trmm_on_device;
  int sup_size_unblocked;
  int sup_size_blocked;
  integer_view_host_t diag_kernel_type_host;
  integer_view_t diag_kernel_type;
  integer_view_host_t kernel_type_host;
  integer_view_t kernel_type;

  // permutation
  bool perm_avail;
  integer_view_t perm;

  // graphs
  host_graph_t original_graph_host;  // graph on host before merge (only if merged)
  host_graph_t graph_host;           // mirror of graph on host
  graph_t graph;

  // crsmat
  crsmat_t crsmat;

  // for supernodal spmv
  bool spmv_trans;
  crsmat_list_t sub_crsmats;
  crsmat_list_t diag_blocks;

  int num_streams;
#if defined(KOKKOS_ENABLE_CUDA)
  cudaStream_t *cuda_streams;
#endif

  // verbose
  bool verbose;
#endif

 public:
  SPTRSVHandle(SPTRSVAlgorithm choice, const size_type nrows_, bool lower_tri_, const size_type block_size_ = 0,
               bool symbolic_complete_ = false, bool numeric_complete_ = false)
      :
#ifdef KOKKOSKERNELS_SPTRSV_CUDAGRAPHSUPPORT
        cudagraphCreated(false),
        sptrsvCudaGraph(nullptr),
#endif
        nrows(nrows_),
        lower_tri(lower_tri_),
        algm(choice),
        level_list(),
        nodes_per_level(),
        hnodes_per_level(),
        nodes_grouped_by_level(),
        hnodes_grouped_by_level(),
        nlevel(0),
        block_size(block_size_),
        team_size(-1),
        vector_size(-1),
        stored_diagonal(false),
        diagonal_offsets(),
        diagonal_values(),  // inserted by rowid
        hdiagonal_offsets(),
        hdiagonal_values(),
        h_chain_ptr(),
        num_chain_entries(0),
        chain_threshold(-1),
        symbolic_complete(symbolic_complete_),
        numeric_complete(numeric_complete_),
        require_symbolic_lvlsched_phase(false),
        require_symbolic_chain_phase(false)
#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
        ,
        cuSPARSEHandle(nullptr),
        tmp_int_rowmap(),
        tmp_int64_rowmap()
#endif
#ifdef KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV
        ,
        unit_diag(false),
        merge_supernodes(false),
        invert_diagonal(true),
        invert_offdiagonal(false),
        etree(nullptr),
        trmm_on_device(true),
        sup_size_unblocked(100),
        sup_size_blocked(200),
        perm_avail(false),
        spmv_trans(true),
        verbose(false)
#endif
  {
    this->set_if_algm_require_symb_lvlsched();
    this->set_if_algm_require_symb_chain();

    // Check a few prerequisites before allowing users
    // to run with the cusparse implementation of sptrsv.
    if (algm == SPTRSVAlgorithm::SPTRSV_CUSPARSE) {
#if !defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
      throw(
          std::runtime_error("sptrsv handle: SPTRSV_CUSPARSE requested but "
                             "cuSPARSE TPL not enabled."));
#else
      if (!std::is_same_v<HandleExecSpace, Kokkos::Cuda>) {
        throw(
            std::runtime_error("sptrsv handle: SPTRSV_CUSPARSE requested but "
                               "HandleExecSpace is not Kokkos::CUDA."));
      }
#endif
    }

#ifdef KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV
    if (lower_tri) {
      // lower-triangular is stored in CSC
      col_major = true;
    } else {
      // upper-triangular is stored in CSR
      col_major = false;
    }
#endif
  }

#ifdef KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV
  // set nsuper and supercols (# of supernodes, and map from supernode to column
  // id
  template <class input_int_type>
  void set_supernodes(signed_integral_t nsuper_, input_int_type *supercols_, int *etree_) {
    // set etree (just wrap etree in a view)
    this->etree_host = integer_view_host_t(etree_, nsuper_);
    // set supernodes (make a copy, from input_int_type to int)
    integer_view_host_t supercols_view = integer_view_host_t("supercols", 1 + nsuper_);
    for (signed_integral_t i = 0; i <= nsuper_; i++) {
      supercols_view(i) = supercols_[i];
    }

    set_supernodes(nsuper_, supercols_view, etree_);
  }

  void set_supernodes(signed_integral_t nsuper_, integer_view_host_t supercols_, int *etree_) {
    // set etree
    this->etree_host = integer_view_host_t(etree_, nsuper_);
    // set supernodes
    set_supernodes(nsuper_, supercols_);
  }

  void set_supernodes(signed_integral_t nsuper_, integer_view_host_t supercols_view) {
    this->nsuper = nsuper_;

    // supercols
    integer_view_host_t supercols_subview(supercols_view.data(), 1 + nsuper_);
    this->supercols_host = integer_view_host_t("supercols_host", 1 + nsuper_);
    Kokkos::deep_copy(this->supercols_host, supercols_subview);

    this->supercols = integer_view_t("supercols", 1 + nsuper_);
    Kokkos::deep_copy(this->supercols, this->supercols_host);

    // workspace offset
    this->work_offset_host = integer_view_host_t("workoffset_host", nsuper_);
    this->work_offset      = integer_view_t("workoffset", nsuper_);

    // kernel type
    this->diag_kernel_type_host = integer_view_host_t("diag_kernel_type_host", nsuper_);
    this->diag_kernel_type      = integer_view_t("diag_kernel_type", nsuper_);
    this->kernel_type_host      = integer_view_host_t("kernel_type_host", nsuper_);
    this->kernel_type           = integer_view_t("kernel_type", nsuper_);

    // number of streams
    this->num_streams = 0;
  }

  // set lower/upper triangular
  void set_lower_tri(bool lower_tri_) { lower_tri = lower_tri_; }

  // set supernodal dag
  void set_supernodal_dag(host_graph_t dag_) { this->dag_host = dag_; }

  // return number of supernodes
  signed_integral_t get_num_supernodes() { return this->nsuper; }

  // return map to supernode to column id
  const int *get_supercols() { return this->supercols.data(); }

  const int *get_supercols_host() { return this->supercols_host.data(); }

  // return parents info in etree of supernodes
  const int *get_etree_parents() { return this->etree_host.data(); }

  // return parents info in etree of supernodes
  host_graph_t get_supernodal_dag() { return this->dag_host; }

  // workspace size
  void set_workspace_size(signed_integral_t lwork_) {
    this->lwork = lwork_;
    this->work  = workspace_t("work", lwork);
  }
  signed_integral_t get_workspace_size() { return this->lwork; }

  // workspace
  KOKKOS_INLINE_FUNCTION
  workspace_t get_workspace() const { return this->work; }

  // workspace
  KOKKOS_INLINE_FUNCTION
  integer_view_t get_work_offset() const { return this->work_offset; }

  integer_view_host_t get_work_offset_host() const { return this->work_offset_host; }

  // specify whether too run KokkosKernels::trmm on device or not
  void set_trmm_on_device(bool flag) { this->trmm_on_device = flag; }

  bool get_trmm_on_device() { return trmm_on_device; }

  // supernode size tolerance to pick right kernel type
  int get_supernode_size_unblocked() { return this->sup_size_unblocked; }

  int get_supernode_size_blocked() { return this->sup_size_blocked; }

  void set_supernode_size_unblocked(int size_unblocked) { this->sup_size_unblocked = size_unblocked; }

  void set_supernode_size_blocked(int size_blocked) { this->sup_size_blocked = size_blocked; }

  // specify to merge supernodes
  void set_merge_supernodes(bool flag) { this->merge_supernodes = flag; }

  bool get_merge_supernodes() { return this->merge_supernodes; }

  // specify etree
  void set_etree(int *etree_) {
    // NOTE: make a copy?
    this->etree = etree_;
  }

  int *get_etree() { return this->etree; }

  // specify to invertt diagonal
  void set_invert_diagonal(bool flag) { this->invert_diagonal = flag; }

  bool get_invert_diagonal() { return this->invert_diagonal; }

  // specify to apply the inverse of diagonal to the offdiagonal blocks
  void set_invert_offdiagonal(bool flag) { this->invert_offdiagonal = flag; }

  bool get_invert_offdiagonal() { return this->invert_offdiagonal; }

  // kernel type
  integer_view_host_t get_kernel_type_host() { return this->kernel_type_host; }

  integer_view_host_t get_diag_kernel_type_host() { return this->diag_kernel_type_host; }

  KOKKOS_INLINE_FUNCTION
  integer_view_t get_kernel_type() { return this->kernel_type; }

  KOKKOS_INLINE_FUNCTION
  integer_view_t get_diag_kernel_type() { return this->diag_kernel_type; }

  // permutation vector
  void set_perm(int *perm_) {
    this->perm     = integer_view_t("PermView", nrows);
    auto perm_host = Kokkos::create_mirror_view(this->perm);

    // copy perm to device
    for (int i = 0; i < nrows; i++) {
      perm_host[i] = perm_[i];
    }
    Kokkos::deep_copy(this->perm, perm_host);
    this->perm_avail = true;
  }

  bool has_perm() { return this->perm_avail; }

  // graph on host (before merge)
  void set_original_graph_host(host_graph_t graph_host_) { this->original_graph_host = graph_host_; }

  host_graph_t get_original_graph_host() { return this->original_graph_host; }

  // graph on host
  void set_graph_host(host_graph_t graph_host_) { this->graph_host = graph_host_; }

  host_graph_t get_graph_host() { return this->graph_host; }

  // graph on device
  void set_graph(graph_t graph_) { this->graph = graph_; }

  graph_t get_graph() { return this->graph; }

  // set if unit diagonal
  void set_unit_diagonal(bool unit_diag_) { this->unit_diag = unit_diag_; }
  bool is_unit_diagonal() { return this->unit_diag; }

  // set CSR or CSC format
  void set_column_major(bool col_major_) { this->col_major = col_major_; }

  bool is_column_major() { return this->col_major; }

  // crsmat
  void set_crsmat(crsmat_t crsmat_) { this->crsmat = crsmat_; }

  crsmat_t get_crsmat() { return this->crsmat; }

  // submatrices
  void set_submatrices(crsmat_list_t subcrsmats) { this->sub_crsmats = subcrsmats; }

  crsmat_t get_submatrix(int i) { return this->sub_crsmats[i]; }

  // diagonal subblocks
  void set_diagblocks(crsmat_list_t subcrsmats) { this->diag_blocks = subcrsmats; }

  crsmat_t get_diagblock(int i) { return this->diag_blocks[i]; }

  // spmv option
  void set_transpose_spmv(bool spmv_trans_) { this->spmv_trans = spmv_trans_; }

  bool transpose_spmv() { return this->spmv_trans; }

  // verbose
  void set_verbose(bool verbose_) { this->verbose = verbose_; }

#if defined(KOKKOS_ENABLE_CUDA)
  // streams
  void setNumStreams(int num_streams_) {
    this->num_streams = num_streams_;
    if (num_streams_ > 0) {
      this->cuda_streams = (cudaStream_t *)malloc(num_streams_ * sizeof(cudaStream_t));
      for (int i = 0; i < num_streams_; i++) {
        cudaStreamCreate(&(this->cuda_streams[i]));
      }
    }
  }

  cudaStream_t *getStream(int id) { return &(this->cuda_streams[id]); }
#endif
#endif

  // Requires nrows_ input
  // Allocates all views
  void new_init_handle(const size_type nrows_) {
    // set_nrows(nrows_);
    nrows = nrows_;
    // Assumed that level scheduling occurs during symbolic phase for all
    // algorithms, for now

    // TODO: Set sizes differently/smaller, resize during symbolic to save space
    if (this->require_symbolic_lvlsched_phase == true) {
      set_num_levels(0);
      level_list = signed_nnz_lno_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "level_list"), nrows_);
      Kokkos::deep_copy(level_list, signed_integral_t(-1));
      // The host side views need to be initialized, but the device-side views
      // don't. Symbolic computes on the host (and requires these are 0
      // initialized), and then copies to device.
      hnodes_per_level        = hostspace_nnz_lno_view_t("host nodes_per_level", nrows_);
      hnodes_grouped_by_level = hostspace_nnz_lno_view_t("host nodes_grouped_by_level", nrows_);
      nodes_per_level = nnz_lno_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "nodes_per_level"), nrows_);
      nodes_grouped_by_level =
          nnz_lno_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "nodes_grouped_by_level"), nrows_);

#if 0
      std::cout << "  newinit_handle: level schedule allocs" << std::endl;
      std::cout << "  ll.extent = " << level_list.extent(0) << std::endl;
      std::cout << "  npl.extent = " << nodes_per_level.extent(0) << std::endl;
      std::cout << "  hnpl.extent = " << hnodes_per_level.extent(0) << std::endl;
      std::cout << "  ngbl.extent = " << nodes_grouped_by_level.extent(0) << std::endl;
      std::cout << "  hngbl.extent = " << hnodes_grouped_by_level.extent(0) << std::endl;
#endif
    }

    if (stored_diagonal) {
      diagonal_offsets  = nnz_lno_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "diagonal_offsets"), nrows_);
      diagonal_values   = nnz_scalar_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "diagonal_values"),
                                            nrows_);  // inserted by rowid
      hdiagonal_offsets = Kokkos::create_mirror_view(diagonal_offsets);
      hdiagonal_values  = Kokkos::create_mirror_view(diagonal_values);
    }

    if (this->require_symbolic_chain_phase == true) {
      if (this->chain_threshold == -1) {
        // Need default if chain_threshold not set
        // 0 means every level, regardless of number of nodes, is launched
        // within a kernel
        if (team_size == -1) {
          this->chain_threshold = 0;
          h_chain_ptr           = host_signed_nnz_lno_view_t("h_chain_ptr", this->nrows);
        } else {
          std::cout << "  Warning: chain_threshold was not set - will default "
                       "to team_size = "
                    << this->team_size << "  chain_threshold = " << this->chain_threshold << std::endl;
          this->chain_threshold = this->team_size;
          h_chain_ptr           = host_signed_nnz_lno_view_t("h_chain_ptr", this->nrows);
        }
      } else {
        if (this->team_size >= this->chain_threshold) {
          h_chain_ptr = host_signed_nnz_lno_view_t("h_chain_ptr", this->nrows);
        } else if (this->team_size == -1 && chain_threshold > 0) {
          std::cout << "  Warning: team_size was not set; chain_threshold = " << this->chain_threshold << std::endl;
          std::cout << "  Automatically setting team_size to chain_threshold - "
                       "if this exceeds the hardware limitations relaunch with "
                       "reduced chain_threshold or set a valid team_size"
                    << std::endl;
          this->team_size = this->chain_threshold;
          h_chain_ptr     = host_signed_nnz_lno_view_t("h_chain_ptr", this->nrows);
        } else {
          std::cout << "  EXPERIMENTAL: team_size less than chain size. team_size = " << this->team_size
                    << "  chain_threshold = " << this->chain_threshold << std::endl;
          h_chain_ptr = host_signed_nnz_lno_view_t("h_chain_ptr", this->nrows);
        }
      }
    } else {
      h_chain_ptr           = host_signed_nnz_lno_view_t();
      this->chain_threshold = -1;
    }

#ifdef KOKKOSKERNELS_SPTRSV_CUDAGRAPHSUPPORT
    create_SPTRSVcudaGraphWrapperType();
#endif

    set_num_chain_entries(0);
    set_symbolic_incomplete();
  }

  virtual ~SPTRSVHandle() {
#ifdef KOKKOSKERNELS_SPTRSV_CUDAGRAPHSUPPORT
    destroy_SPTRSVcudaGraphWrapperType();
#endif

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
    this->destroy_cuSPARSE_Handle();
#endif
  };

#ifdef KOKKOSKERNELS_ENABLE_TPL_CUSPARSE
  void create_cuSPARSE_Handle(bool transpose, bool is_lower) {
    this->destroy_cuSPARSE_Handle();
    this->cuSPARSEHandle = new cuSparseHandleType(transpose, is_lower);
  }
  void destroy_cuSPARSE_Handle() {
    if (this->cuSPARSEHandle != nullptr) {
      delete this->cuSPARSEHandle;
      this->cuSPARSEHandle = nullptr;
    }
  }

  SPTRSVcuSparseHandleType *get_cuSparseHandle() { return this->cuSPARSEHandle; }

  void allocate_tmp_int_rowmap(size_type N) {
    tmp_int_rowmap = int_row_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "tmp_int_rowmap"), N);
  }
  template <typename RowViewType>
  int_row_view_t get_int_rowmap_view_copy(const RowViewType &rowmap) {
    Kokkos::deep_copy(tmp_int_rowmap, rowmap);
    return tmp_int_rowmap;
  }
  template <typename RowViewType>
  int *get_int_rowmap_ptr_copy(const RowViewType &rowmap) {
    Kokkos::deep_copy(tmp_int_rowmap, rowmap);
    Kokkos::fence();
    return tmp_int_rowmap.data();
  }
  int_row_view_t get_int_rowmap_view() { return tmp_int_rowmap; }
  int *get_int_rowmap_ptr() { return tmp_int_rowmap.data(); }

  void allocate_tmp_int64_rowmap(size_type N) {
    tmp_int64_rowmap = int64_row_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "tmp_int64_rowmap"), N);
  }
  template <typename RowViewType>
  int64_t *get_int64_rowmap_ptr_copy(const RowViewType &rowmap) {
    Kokkos::deep_copy(tmp_int64_rowmap, rowmap);
    Kokkos::fence();
    return tmp_int64_rowmap.data();
  }
  int64_t *get_int64_rowmap_ptr() { return tmp_int64_rowmap.data(); }
#endif

  bool algm_requires_symb_lvlsched() const { return require_symbolic_lvlsched_phase; }

  bool algm_requires_symb_chain() const { return require_symbolic_chain_phase; }

  // Can change the algorithm to a "Compatible algorithms" - for ease in some
  // testing cases
  void set_algorithm(SPTRSVAlgorithm choice) {
    if (algm != choice) {
      algm = choice;
    }
  }

  KOKKOS_INLINE_FUNCTION
  SPTRSVAlgorithm get_algorithm() { return algm; }

  KOKKOS_INLINE_FUNCTION
  signed_nnz_lno_view_t get_level_list() const { return level_list; }

  inline host_signed_nnz_lno_view_t get_host_level_list() const {
    auto hlevel_list = Kokkos::create_mirror_view(this->level_list);
    Kokkos::deep_copy(hlevel_list, this->level_list);
    return hlevel_list;
  }

  void set_stored_diagonal(const bool stored_diagonal_) { stored_diagonal = stored_diagonal_; }

  KOKKOS_INLINE_FUNCTION
  nnz_lno_view_t get_diagonal_offsets() const { return diagonal_offsets; }

  KOKKOS_INLINE_FUNCTION
  nnz_scalar_view_t get_diagonal_values() const { return diagonal_values; }

  KOKKOS_INLINE_FUNCTION
  host_nnz_lno_view_t get_host_diagonal_offsets() const { return hdiagonal_offsets; }

  KOKKOS_INLINE_FUNCTION
  host_nnz_scalar_view_t get_host_diagonal_values() const { return hdiagonal_values; }

  inline host_signed_nnz_lno_view_t get_host_chain_ptr() const { return h_chain_ptr; }

  KOKKOS_INLINE_FUNCTION
  nnz_lno_view_t get_nodes_per_level() const { return nodes_per_level; }

  inline hostspace_nnz_lno_view_t get_host_nodes_per_level() const { return hnodes_per_level; }

  KOKKOS_INLINE_FUNCTION
  nnz_lno_view_t get_nodes_grouped_by_level() const { return nodes_grouped_by_level; }

  inline hostspace_nnz_lno_view_t get_host_nodes_grouped_by_level() const { return hnodes_grouped_by_level; }

  KOKKOS_INLINE_FUNCTION
  size_type get_nrows() const { return nrows; }
  void set_nrows(const size_type nrows_) { this->nrows = nrows_; }

  void reset_chain_threshold(const signed_integral_t threshold) {
    if (threshold != this->chain_threshold || h_chain_ptr.span() == 0) {
      this->chain_threshold = threshold;
      if (this->team_size >= this->chain_threshold) {
        //  h_chain_ptr = host_signed_nnz_lno_view_t("h_chain_ptr",
        //  this->nrows);
      } else if (this->team_size == -1 && chain_threshold > 0) {
        // std::cout << "  Warning: team_size was not set  team_size = " <<
        // this->team_size << "  chain_threshold = " << this->chain_threshold <<
        // std::endl; std::cout << "  Automatically setting team_size to
        // chain_threshold - if this exceeds the hardware limitation a runtime
        // error will occur during kernel launch - reduce chain_threshold in
        // that case" << std::endl;
        this->team_size = this->chain_threshold;
        //  h_chain_ptr = host_signed_nnz_lno_view_t("h_chain_ptr",
        //  this->nrows);
      } else {
        std::cout << "  EXPERIMENTAL: team_size < chain_size: team_size = " << this->team_size
                  << "  chain_threshold = " << this->chain_threshold << std::endl;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  signed_integral_t get_chain_threshold() const { return this->chain_threshold; }

  bool is_lower_tri() const { return lower_tri; }
  bool is_upper_tri() const { return !lower_tri; }

  bool is_symbolic_complete() const { return symbolic_complete; }
  bool is_numeric_complete() const { return numeric_complete; }

  bool is_stored_diagonal() const { return stored_diagonal; }

  KOKKOS_INLINE_FUNCTION
  size_type get_num_levels() const { return nlevel; }

  void set_num_levels(size_type nlevels_) { this->nlevel = nlevels_; }

  KOKKOS_INLINE_FUNCTION
  size_type get_block_size() const { return block_size; }

  KOKKOS_INLINE_FUNCTION
  void set_block_size(const size_type block_size_) { this->block_size = block_size_; }

  bool is_block_enabled() const { return block_size > 0; }
  void set_symbolic_complete() { this->symbolic_complete = true; }
  void set_symbolic_incomplete() { this->symbolic_complete = false; }

  void set_numeric_complete() { this->numeric_complete = true; }
  void set_numeric_incomplete() { this->numeric_complete = false; }

  KOKKOS_INLINE_FUNCTION
  int get_team_size() const { return this->team_size; }
  // Called by user at setup - should only set a value, no alloc
  void set_team_size(const int ts) { this->team_size = ts; }

  KOKKOS_INLINE_FUNCTION
  int get_vector_size() const { return this->vector_size; }
  // Called by user at setup - should only set a value, no alloc
  void set_vector_size(const int vs) { this->vector_size = vs; }

  KOKKOS_INLINE_FUNCTION
  int get_num_chain_entries() const { return this->num_chain_entries; }
  void set_num_chain_entries(const int nce) { this->num_chain_entries = nce; }

  void print_algorithm() {
    if (algm == SPTRSVAlgorithm::SEQLVLSCHD_RP) std::cout << "SEQLVLSCHD_RP" << std::endl;
    ;

    if (algm == SPTRSVAlgorithm::SEQLVLSCHD_TP1) std::cout << "SEQLVLSCHD_TP1" << std::endl;
    ;
    /*
        if ( algm == SPTRSVAlgorithm::SEQLVLSCHED_TP2 ) {
          std::cout << "SEQLVLSCHED_TP2" << std::endl;;
          std::cout << "WARNING: With CUDA this is currently only reliable with
       int-int ordinal-offset pair" << std::endl;
        }
    */
    if (algm == SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN) std::cout << "SEQLVLSCHD_TP1CHAIN" << std::endl;
    ;

    if (algm == SPTRSVAlgorithm::SPTRSV_CUSPARSE) std::cout << "SPTRSV_CUSPARSE" << std::endl;
    ;

    if (algm == SPTRSVAlgorithm::SUPERNODAL_NAIVE) std::cout << "SUPERNODAL_NAIVE" << std::endl;

    if (algm == SPTRSVAlgorithm::SUPERNODAL_ETREE) std::cout << "SUPERNODAL_ETREE" << std::endl;

    if (algm == SPTRSVAlgorithm::SUPERNODAL_DAG) std::cout << "SUPERNODAL_DAG" << std::endl;

    if (algm == SPTRSVAlgorithm::SUPERNODAL_SPMV) std::cout << "SUPERNODAL_SPMV" << std::endl;

    if (algm == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG) std::cout << "SUPERNODAL_SPMV_DAG" << std::endl;
  }

  std::string return_algorithm_string() {
    std::string ret_string;

    if (algm == SPTRSVAlgorithm::SEQLVLSCHD_RP) ret_string = "SEQLVLSCHD_RP";

    if (algm == SPTRSVAlgorithm::SEQLVLSCHD_TP1) ret_string = "SEQLVLSCHD_TP1";
    /*
        if ( algm == SPTRSVAlgorithm::SEQLVLSCHED_TP2 )
          ret_string = "SEQLVLSCHED_TP2";
    */
    if (algm == SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN) ret_string = "SEQLVLSCHD_TP1CHAIN";

    if (algm == SPTRSVAlgorithm::SPTRSV_CUSPARSE) ret_string = "SPTRSV_CUSPARSE";

    return ret_string;
  }

  inline SPTRSVAlgorithm StringToSPTRSVAlgorithm(std::string &name) {
    if (name == "SPTRSV_DEFAULT")
      return SPTRSVAlgorithm::SEQLVLSCHD_RP;
    else if (name == "SPTRSV_RANGEPOLICY")
      return SPTRSVAlgorithm::SEQLVLSCHD_RP;
    else if (name == "SPTRSV_TEAMPOLICY1")
      return SPTRSVAlgorithm::SEQLVLSCHD_TP1;
    /*else if(name=="SPTRSV_TEAMPOLICY2")       return
     * SPTRSVAlgorithm::SEQLVLSCHED_TP2;*/
    else if (name == "SPTRSV_TEAMPOLICY1CHAIN")
      return SPTRSVAlgorithm::SEQLVLSCHD_TP1CHAIN;
    else if (name == "SPTRSV_CUSPARSE")
      return SPTRSVAlgorithm::SPTRSV_CUSPARSE;
    else
      throw std::runtime_error("Invalid SPTRSVAlgorithm name");
  }
};

}  // namespace Experimental
}  // namespace KokkosSparse

#endif
