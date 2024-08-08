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
#include <KokkosKernels_ExecSpaceUtils.hpp>
#include "KokkosGraph_Distance1ColorHandle.hpp"
#include "KokkosGraph_Distance2ColorHandle.hpp"
#include "KokkosSparse_gauss_seidel_handle.hpp"
#include "KokkosSparse_spgemm_handle.hpp"
#include "KokkosSparse_spadd_handle.hpp"
#include "KokkosSparse_sptrsv_handle.hpp"
#include "KokkosSparse_spiluk_handle.hpp"
#include "KokkosSparse_par_ilut_handle.hpp"
#include "KokkosSparse_gmres_handle.hpp"
#include "KokkosKernels_default_types.hpp"

#ifndef _KOKKOSKERNELHANDLE_HPP
#define _KOKKOSKERNELHANDLE_HPP

namespace KokkosKernels {

namespace Experimental {

template <class size_type_, class lno_t_, class scalar_t_, class ExecutionSpace, class TemporaryMemorySpace,
          class PersistentMemorySpace>
class KokkosKernelsHandle {
  static_assert(std::is_signed<lno_t_>::value,
                "KokkosKernelsHandle requires that lno_t_ (ordinal) is a "
                "signed integer type.");

 public:
  typedef typename ExecutionSpace::execution_space HandleExecSpace;
  typedef typename TemporaryMemorySpace::memory_space HandleTempMemorySpace;
  typedef typename PersistentMemorySpace::memory_space HandlePersistentMemorySpace;
  // typedef
  // Kokkos::Device<ExecutionSpace::execution_space,ExecutionSpace::memory_space>
  // HandleExecSpace; typedef
  // Kokkos::Device<TemporaryMemorySpace::execution_space,TemporaryMemorySpace::memory_space>
  // HandleTempMemorySpace; typedef
  // Kokkos::Device<PersistentMemorySpace::execution_space,PersistentMemorySpace::memory_space>
  // HandlePersistentMemorySpace;

  typedef typename std::remove_const<size_type_>::type size_type;
  typedef const size_type const_size_type;

  typedef typename std::remove_const<lno_t_>::type nnz_lno_t;
  typedef const nnz_lno_t const_nnz_lno_t;

  typedef typename std::remove_const<scalar_t_>::type nnz_scalar_t;
  typedef const nnz_scalar_t const_nnz_scalar_t;

  template <typename right_size_type_, typename right_lno_t_, typename right_scalar_t_, typename right_ExecutionSpace,
            typename right_TemporaryMemorySpace, typename right_PersistentMemorySpace>
  // KokkosKernelsHandle<const_size_type,const_nnz_lno_t, const_nnz_scalar_t,
  // HandleExecSpace, HandleTempMemorySpace, HandlePersistentMemorySpace>
  // &operator=
  KokkosKernelsHandle(KokkosKernelsHandle<right_size_type_, right_lno_t_, right_scalar_t_, right_ExecutionSpace,
                                          right_TemporaryMemorySpace, right_PersistentMemorySpace> &right_side_handle) {
    static_assert(std::is_same<size_type_, const_size_type>::value,
                  "Kernel handle left hand side should have const size type in "
                  "assignment");
    static_assert(std::is_same<lno_t_, const_nnz_lno_t>::value,
                  "Kernel handle left hand side should have const lno type in "
                  "assignment");
    static_assert(std::is_same<scalar_t_, const_nnz_scalar_t>::value,
                  "Kernel handle left hand side should have const scalar type "
                  "in assignment");

    static_assert(std::is_same<ExecutionSpace, HandleExecSpace>::value,
                  "Kernel handle left hand side should have execution space in "
                  "assignment");
    static_assert(std::is_same<TemporaryMemorySpace, HandleTempMemorySpace>::value,
                  "Kernel handle left hand side should have temp memory space in "
                  "assignment");
    static_assert(std::is_same<PersistentMemorySpace, HandlePersistentMemorySpace>::value,
                  "Kernel handle left hand side should have persistent memory space in "
                  "assignment");

    typedef typename std::remove_const<right_size_type_>::type nonconst_right_size_type;
    typedef const nonconst_right_size_type const_right_size_type;

    typedef typename std::remove_const<right_lno_t_>::type nonconst_right_lno_t;
    typedef const nonconst_right_lno_t const_right_lno_t;

    typedef typename std::remove_const<right_scalar_t_>::type nonconst_right_scalar_t;
    typedef const nonconst_right_scalar_t const_right_scalar_t;

    static_assert(std::is_same<size_type_, const_right_size_type>::value,
                  "Kernel handle left and right sides should have same size "
                  "type in assignment");
    static_assert(std::is_same<lno_t_, const_right_lno_t>::value,
                  "Kernel handle left and right sides should have same lno "
                  "type in assignment");
    static_assert(std::is_same<scalar_t_, const_right_scalar_t>::value,
                  "Kernel handle left and right sides should have same scalar "
                  "type in assignment");

    static_assert(
        std::is_same<typename ExecutionSpace::execution_space, typename right_ExecutionSpace::execution_space>::value,
        "Kernel handle left and right sides should have same execution_space "
        "in assignment");
    /*
    static_assert (std::is_same<typename TemporaryMemorySpace::execution_space,
    typename right_TemporaryMemorySpace::execution_space>::value, "Kernel handle
    left and right sides should have same TemporaryMemorySpace in assignment");
    static_assert (std::is_same<typename PersistentMemorySpace::execution_space,
    typename right_PersistentMemorySpace::execution_space>::value, "Kernel
    handle left and right sides should have same PersistentMemorySpace in
    assignment");
      */
    static_assert(
        std::is_same<typename ExecutionSpace::memory_space, typename right_ExecutionSpace::memory_space>::value,
        "Kernel handle left and right sides should have same ExecutionSpace in "
        "assignment");
    static_assert(std::is_same<typename TemporaryMemorySpace::memory_space,
                               typename right_TemporaryMemorySpace::memory_space>::value,
                  "Kernel handle left and right sides should have same "
                  "TemporaryMemorySpace in assignment");
    static_assert(std::is_same<typename PersistentMemorySpace::memory_space,
                               typename right_PersistentMemorySpace::memory_space>::value,
                  "Kernel handle left and right sides should have same "
                  "PersistentMemorySpace in assignment");

    this->gcHandle    = right_side_handle.get_graph_coloring_handle();
    this->gcHandle_d2 = right_side_handle.get_distance2_graph_coloring_handle();

    this->gsHandle = right_side_handle.get_gs_handle();
    // ---------------------------------------- //
    // Handles for Classical GS (inner SpTRSV)
    this->gs_sptrsvLHandle = right_side_handle.get_gs_sptrsvL_handle();
    this->gs_sptrsvUHandle = right_side_handle.get_gs_sptrsvU_handle();

    this->spgemmHandle = right_side_handle.get_spgemm_handle();
    this->spaddHandle  = right_side_handle.get_spadd_handle();

    this->sptrsvHandle   = right_side_handle.get_sptrsv_handle();
    this->spilukHandle   = right_side_handle.get_spiluk_handle();
    this->par_ilutHandle = right_side_handle.get_par_ilut_handle();
    this->gmresHandle    = right_side_handle.get_gmres_handle();

    this->team_work_size      = right_side_handle.get_set_team_work_size();
    this->shared_memory_size  = right_side_handle.get_shmem_size();
    this->suggested_team_size = right_side_handle.get_set_suggested_team_size();

    this->my_exec_space          = right_side_handle.get_handle_exec_space();
    this->use_dynamic_scheduling = right_side_handle.is_dynamic_scheduling();
    this->KKVERBOSE              = right_side_handle.get_verbose();
    this->vector_size            = right_side_handle.get_set_suggested_vector_size();

    is_owner_of_the_gc_handle = false;
    // ---------------------------------------- //
    // Handles for Classical GS (inner SpTRSV)
    is_owner_of_the_gs_sptrsvL_handle = false;
    is_owner_of_the_gs_sptrsvU_handle = false;
    // ---------------------------------------- //
    is_owner_of_the_d2_gc_handle    = false;
    is_owner_of_the_gs_handle       = false;
    is_owner_of_the_spgemm_handle   = false;
    is_owner_of_the_spadd_handle    = false;
    is_owner_of_the_sptrsv_handle   = false;
    is_owner_of_the_spiluk_handle   = false;
    is_owner_of_the_par_ilut_handle = false;
    is_owner_of_the_gmres_handle    = false;
    // return *this;
  }

  typedef typename KokkosGraph::GraphColoringHandle<const_size_type, const_nnz_lno_t, const_nnz_lno_t, HandleExecSpace,
                                                    HandleTempMemorySpace, HandlePersistentMemorySpace>
      GraphColoringHandleType;

  typedef typename KokkosGraph::GraphColorDistance2Handle<const_size_type, const_nnz_lno_t, const_nnz_lno_t,
                                                          HandleExecSpace, HandleTempMemorySpace,
                                                          HandlePersistentMemorySpace>
      GraphColorDistance2HandleType;

  typedef typename KokkosSparse::GaussSeidelHandle<const_size_type, const_nnz_lno_t, const_nnz_scalar_t,
                                                   HandleExecSpace, HandleTempMemorySpace, HandlePersistentMemorySpace>
      GaussSeidelHandleType;

  // These are subclasses of GaussSeidelHandleType.
  typedef
      typename KokkosSparse::PointGaussSeidelHandle<const_size_type, const_nnz_lno_t, const_nnz_scalar_t,
                                                    HandleExecSpace, HandleTempMemorySpace, HandlePersistentMemorySpace>
          PointGaussSeidelHandleType;
  typedef typename KokkosSparse::ClusterGaussSeidelHandle<const_size_type, const_nnz_lno_t, const_nnz_scalar_t,
                                                          HandleExecSpace, HandleTempMemorySpace,
                                                          HandlePersistentMemorySpace>
      ClusterGaussSeidelHandleType;
  // ---------------------------------------- //
  // These are for Two-stage Gauss-Seidel
  typedef typename KokkosSparse::TwoStageGaussSeidelHandle<const_size_type, const_nnz_lno_t, const_nnz_scalar_t,
                                                           HandleExecSpace, HandleTempMemorySpace,
                                                           HandlePersistentMemorySpace>
      TwoStageGaussSeidelHandleType;
  typedef KokkosKernelsHandle<const_size_type, const_nnz_lno_t, const_nnz_scalar_t, HandleExecSpace,
                              HandleTempMemorySpace, HandleTempMemorySpace>
      TwoStageGaussSeidelSPTRSVHandleType;
  // ---------------------------------------- //

  typedef typename KokkosSparse::SPGEMMHandle<const_size_type, const_nnz_lno_t, const_nnz_scalar_t, HandleExecSpace,
                                              HandleTempMemorySpace, HandlePersistentMemorySpace>
      SPGEMMHandleType;

  typedef typename Kokkos::View<nnz_scalar_t *, HandleTempMemorySpace> in_scalar_nnz_view_t;

  typedef typename Kokkos::View<size_type *, HandleTempMemorySpace> row_lno_temp_work_view_t;
  typedef typename Kokkos::View<size_type *, HandleTempMemorySpace> size_type_temp_work_view_t;
  typedef typename Kokkos::View<size_type *, HandlePersistentMemorySpace> row_lno_persistent_work_view_t;
  typedef typename Kokkos::View<size_type *, HandlePersistentMemorySpace> size_type_persistent_work_view_t;
  typedef typename row_lno_persistent_work_view_t::HostMirror row_lno_persistent_work_host_view_t;  // Host view type
  typedef
      typename size_type_persistent_work_view_t::HostMirror size_type_persistent_work_host_view_t;  // Host view type
  typedef typename Kokkos::View<nnz_scalar_t *, HandleTempMemorySpace> scalar_temp_work_view_t;
  typedef typename Kokkos::View<nnz_scalar_t *, HandlePersistentMemorySpace> scalar_persistent_work_view_t;
  typedef typename Kokkos::View<nnz_scalar_t **, default_layout, HandlePersistentMemorySpace>
      scalar_persistent_work_view2d_t;
  typedef typename Kokkos::View<nnz_lno_t *, HandleTempMemorySpace> nnz_lno_temp_work_view_t;
  typedef typename Kokkos::View<nnz_lno_t *, HandlePersistentMemorySpace> nnz_lno_persistent_work_view_t;
  typedef typename nnz_lno_persistent_work_view_t::HostMirror nnz_lno_persistent_work_host_view_t;  // Host view type
  typedef typename Kokkos::View<bool *, HandlePersistentMemorySpace> bool_persistent_view_t;
  typedef typename Kokkos::View<bool *, HandleTempMemorySpace> bool_temp_view_t;

  typedef typename KokkosSparse::SPADDHandle<row_lno_temp_work_view_t, nnz_lno_temp_work_view_t,
                                             scalar_temp_work_view_t, HandleExecSpace, HandleTempMemorySpace>
      SPADDHandleType;

  typedef typename KokkosSparse::Experimental::SPTRSVHandle<const_size_type, const_nnz_lno_t, const_nnz_scalar_t,
                                                            HandleExecSpace, HandleTempMemorySpace,
                                                            HandlePersistentMemorySpace>
      SPTRSVHandleType;
#ifdef KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV
  using integer_view_host_t = typename SPTRSVHandleType::integer_view_host_t;
#endif

  typedef typename KokkosSparse::Experimental::SPILUKHandle<const_size_type, const_nnz_lno_t, const_nnz_scalar_t,
                                                            HandleExecSpace, HandleTempMemorySpace,
                                                            HandlePersistentMemorySpace>
      SPILUKHandleType;

  typedef typename KokkosSparse::Experimental::PAR_ILUTHandle<const_size_type, const_nnz_lno_t, const_nnz_scalar_t,
                                                              HandleExecSpace, HandleTempMemorySpace,
                                                              HandlePersistentMemorySpace>
      PAR_ILUTHandleType;

  typedef typename KokkosSparse::Experimental::GMRESHandle<const_size_type, const_nnz_lno_t, const_nnz_scalar_t,
                                                           HandleExecSpace, HandleTempMemorySpace,
                                                           HandlePersistentMemorySpace>
      GMRESHandleType;

 private:
  GraphColoringHandleType *gcHandle;
  GraphColorDistance2HandleType *gcHandle_d2;

  GaussSeidelHandleType *gsHandle;
  // ---------------------------------------- //
  // Handles for Classical GS (inner SpTRSV)
  // NOTE: move these handles inside GS handle
  TwoStageGaussSeidelSPTRSVHandleType *gs_sptrsvLHandle;
  TwoStageGaussSeidelSPTRSVHandleType *gs_sptrsvUHandle;
  // ---------------------------------------- //
  SPGEMMHandleType *spgemmHandle;
  SPADDHandleType *spaddHandle;
  SPTRSVHandleType *sptrsvHandle;
  SPILUKHandleType *spilukHandle;
  PAR_ILUTHandleType *par_ilutHandle;
  GMRESHandleType *gmresHandle;

  int team_work_size;
  size_t shared_memory_size;
  int suggested_team_size;

  KokkosKernels::Impl::ExecSpaceType my_exec_space;
  bool use_dynamic_scheduling;
  bool KKVERBOSE;
  int vector_size;

  bool is_owner_of_the_gc_handle;
  bool is_owner_of_the_d2_gc_handle;
  bool is_owner_of_the_gs_handle;
  // ---------------------------------------- //
  // Handles for Classical GS (inner SpTRSV)
  bool is_owner_of_the_gs_sptrsvL_handle;
  bool is_owner_of_the_gs_sptrsvU_handle;
  // ---------------------------------------- //
  bool is_owner_of_the_spgemm_handle;
  bool is_owner_of_the_spadd_handle;
  bool is_owner_of_the_sptrsv_handle;
  bool is_owner_of_the_spiluk_handle;
  bool is_owner_of_the_par_ilut_handle;
  bool is_owner_of_the_gmres_handle;

 public:
  KokkosKernelsHandle()
      : gcHandle(NULL),
        gcHandle_d2(NULL),
        gsHandle(NULL)
        // Handles for Classical GS (inner SpTRSV)
        ,
        gs_sptrsvLHandle(NULL),
        gs_sptrsvUHandle(NULL)
        // ---------------------------------------- //
        ,
        spgemmHandle(NULL),
        spaddHandle(NULL),
        sptrsvHandle(NULL),
        spilukHandle(NULL),
        par_ilutHandle(NULL),
        gmresHandle(NULL),
        team_work_size(-1),
        shared_memory_size(16128),
        suggested_team_size(-1),
        my_exec_space(KokkosKernels::Impl::kk_get_exec_space_type<HandleExecSpace>()),
        use_dynamic_scheduling(true),
        KKVERBOSE(false),
        vector_size(-1),
        is_owner_of_the_gc_handle(true),
        is_owner_of_the_d2_gc_handle(true),
        is_owner_of_the_gs_handle(true)
        // Handles for Classical GS (inner SpTRSV)
        ,
        is_owner_of_the_gs_sptrsvL_handle(true),
        is_owner_of_the_gs_sptrsvU_handle(true)
        // ---------------------------------------- //
        ,
        is_owner_of_the_spgemm_handle(true),
        is_owner_of_the_spadd_handle(true),
        is_owner_of_the_sptrsv_handle(true),
        is_owner_of_the_spiluk_handle(true),
        is_owner_of_the_par_ilut_handle(true),
        is_owner_of_the_gmres_handle(true) {}

  ~KokkosKernelsHandle() {
    this->destroy_gs_handle();
    // ---------------------------------------- //
    // Handles for Classical GS (inner SpTRSV)
    this->destroy_gs_sptrsvL_handle();
    this->destroy_gs_sptrsvU_handle();
    // ---------------------------------------- //
    this->destroy_graph_coloring_handle();
    this->destroy_distance2_graph_coloring_handle();
    this->destroy_spgemm_handle();
    this->destroy_spadd_handle();
    this->destroy_sptrsv_handle();
    this->destroy_spiluk_handle();
    this->destroy_par_ilut_handle();
    this->destroy_gmres_handle();
  }

  void set_verbose(bool verbose_) { this->KKVERBOSE = verbose_; }
  bool get_verbose() { return this->KKVERBOSE; }
  /**
   * \brief Sets the chunk size to be handled by each team.
   * Big chunks are good, they reduce the overhead that occurs at the every
   * iteration. However, too big chunks causes load-imbalances. If this is not
   * set, the algorithm will use teamsize as chunk sizes. On GPUs, this is
   * usually okay. On CPUs, we get better performance for bigger chunksizes.
   * This depends on how big is your work. On SPGEMM a single work is the
   * multiplication for a row. For laplace we had 50-500 FLOPs per row for AP
   * and RxAP. Best was to use 256 chunks. ~12800 - 128K flops per team was
   * fine. For brick we had 250 - 2500 FLOPs per  row for AP and RxAP. Best was
   * to use 256 chunks. 64K to 640K flops per team was fine. But on EMPIRE 1200
   * - 20800 FLOPs per row for AP and RxAP. Best was to use 16 chunks. 19K to
   * 330K was flops per team was fine. This bases on the load balancing issues
   * introduced. A general way is to have 100K flops per team for SPGEMM. \param
   * team_work_size_: input, the size of the chunks.
   */
  void set_team_work_size(const int team_work_size_) { this->team_work_size = team_work_size_; }
  int get_set_team_work_size() { return this->team_work_size; }
  /**
   * \brief Returns the enum type for the execution space.
   */
  KokkosKernels::Impl::ExecSpaceType get_handle_exec_space() { return this->my_exec_space; }

  /// \brief Returns the suggested team work size. If set with
  /// set_team_work_size, it will return the set value. Otherwise it will return
  /// the teamsize.
  /// \param team_size input, team size used by the kernel.
  /// \param concurrency filler for concurrency
  /// \param overall_work_size filler for overall_work_size
  int get_team_work_size(const int team_size, [[maybe_unused]] const int concurrency,
                         [[maybe_unused]] const nnz_lno_t overall_work_size) {
    if (this->team_work_size != -1) {
      return this->team_work_size;
    } else {
      if (my_exec_space == KokkosKernels::Impl::Exec_CUDA || my_exec_space == KokkosKernels::Impl::Exec_HIP) {
        return team_size;
      } else {
        return 16;
      }
    }
  }

  /**
   * \brief Sets whether to use dynamic scheduling or
   * not for those kernels where load-imbalances might occur.
   * @param is_dynamic: true or false -> dynamic or static scheduling..
   */
  void set_dynamic_scheduling(const bool is_dynamic) { this->use_dynamic_scheduling = is_dynamic; }

  /**
   * \brief Returns true or false, use dynamic scheduling or not.
   */
  bool is_dynamic_scheduling() { return this->use_dynamic_scheduling; }

  /// \brief sets the shared memory size to be used by the kernels using shared
  /// memory on GPUs.
  /// \param shared_memory_size_ input, shared memory size to be used by the
  /// kernel.
  void set_shmem_size(const size_t shared_memory_size_) { this->shared_memory_size = shared_memory_size_; }

  /**
   * \brief Returns the shared memory size suggested by the handle.
   */
  size_t get_shmem_size() { return shared_memory_size; }

  /**
   * \brief Returns the suggested vector size based on the execution space.
   * Basically, it calculates the average degree, and rounds it up to the
   * closes power of 2 on GPUs. on CPUs, it returns 1.
   * \param nr: number of rows, of vertices.
   * \param nnz: number of nonzeroes, or edges.
   */
  int get_suggested_vector_size(const size_t nr, const size_t nnz) {
    if (vector_size == -1) {
      return KokkosKernels::Impl::kk_get_suggested_vector_size(nr, nnz, my_exec_space);
    } else {
      return vector_size;
    }
  }

  void set_suggested_vector_size(int vector_size_) { this->vector_size = vector_size_; }

  int get_set_suggested_vector_size() { return this->vector_size; }
  /**
   * \brief Sets the team size to be used by the kernels. On GPUs and CPUs
   * usually the defaults are fine. But on CPUs with hyperthreads it might be
   * worth of trying different team sizes.
   * \param suggested_team_size_: team size to set.
   */
  void set_suggested_team_size(const int suggested_team_size_) { this->suggested_team_size = suggested_team_size_; }

  int get_set_suggested_team_size() { return this->suggested_team_size; }

  /// \brief Returns the team size, either set by the user or suggested by the
  /// handle.
  /// \param vector_size_ suggested vector size by the handle.
  int get_suggested_team_size(const int vector_size_) {
    if (this->suggested_team_size != -1) {
      return this->suggested_team_size;
    } else {
      return KokkosKernels::Impl::kk_get_suggested_team_size(vector_size_, my_exec_space);
    }
  }

  // SPGEM
  SPGEMMHandleType *get_spgemm_handle() { return this->spgemmHandle; }
  void create_spgemm_handle(KokkosSparse::SPGEMMAlgorithm spgemm_algo = KokkosSparse::SPGEMM_DEFAULT) {
    this->destroy_spgemm_handle();
    this->is_owner_of_the_spgemm_handle = true;
    this->spgemmHandle                  = new SPGEMMHandleType(spgemm_algo);
  }
  void destroy_spgemm_handle() {
    if (is_owner_of_the_spgemm_handle && this->spgemmHandle != NULL) {
      delete this->spgemmHandle;
      this->spgemmHandle = NULL;
    }
  }

  // Distance-1 Graph Coloring
  GraphColoringHandleType *get_graph_coloring_handle() {
    // (wcmclen): Should there be a check here to make sure we've created a GC
    // handle before
    //            handing the pointer out to something? This is disabled for now
    //            because it gets thrown in tests run by spot check. Moving
    //            forward, we should consider whether a "get the handle ptr,
    //            then allocate" vs. "only give out the handle ptr if it
    //            actually exists" model.
    // if(!this->is_owner_of_the_gc_handle)
    //{
    //  throw std::runtime_error("Graph coloring handle has not been created.");
    //}
    return this->gcHandle;
  }
  void create_graph_coloring_handle(KokkosGraph::ColoringAlgorithm coloring_type = KokkosGraph::COLORING_DEFAULT) {
    this->destroy_graph_coloring_handle();
    this->is_owner_of_the_gc_handle = true;
    this->gcHandle                  = new GraphColoringHandleType();
    this->gcHandle->set_algorithm(coloring_type, true);
    this->gcHandle->set_tictoc(KKVERBOSE);
  }
  void destroy_graph_coloring_handle() {
    if (is_owner_of_the_gc_handle && this->gcHandle != NULL) {
      delete this->gcHandle;
      this->gcHandle = NULL;
    }
  }

  // Distance-2 Graph Coloring
  GraphColorDistance2HandleType *get_distance2_graph_coloring_handle() {
    /* disabled for consistency with `get_graph_coloring_handle()`. See the
    comment there for reasons. if(!this->is_owner_of_the_d2_gc_handle)
    {
      throw std::runtime_error("D2 graph coloring handle has not been
    created.");
    }
    */
    return this->gcHandle_d2;
  }
  void create_distance2_graph_coloring_handle(
      KokkosGraph::GraphColoringAlgorithmDistance2 coloring_type = KokkosGraph::COLORING_D2_DEFAULT) {
    this->destroy_distance2_graph_coloring_handle();
    this->is_owner_of_the_d2_gc_handle = true;
    this->gcHandle_d2                  = new GraphColorDistance2HandleType();
    this->gcHandle_d2->set_algorithm(coloring_type, true);
    this->gcHandle_d2->set_tictoc(KKVERBOSE);
    this->gcHandle_d2->set_verbose(KKVERBOSE);
  }
  void destroy_distance2_graph_coloring_handle() {
    if (is_owner_of_the_d2_gc_handle && this->gcHandle_d2 != NULL) {
      delete this->gcHandle_d2;
      this->gcHandle_d2 = NULL;
    }
  }

  GaussSeidelHandleType *get_gs_handle() { return this->gsHandle; }
  PointGaussSeidelHandleType *get_point_gs_handle() {
    auto pgs = dynamic_cast<PointGaussSeidelHandleType *>(this->gsHandle);
    if (this->gsHandle && !pgs)
      throw std::runtime_error("GaussSeidelHandle exists but is not set up for point-coloring GS.");
    return pgs;
  }
  ClusterGaussSeidelHandleType *get_cluster_gs_handle() {
    auto cgs = dynamic_cast<ClusterGaussSeidelHandleType *>(this->gsHandle);
    if (this->gsHandle && !cgs)
      throw std::runtime_error(
          "GaussSeidelHandle exists but is not set up for cluster-coloring "
          "GS.");
    return cgs;
  }

  // clang-format off
  /**
   * @brief Create a gauss seidel handle object
   *
   * @param handle_exec_space The execution space instance to execute kernels on.
   * @param num_streams The number of streams to allocate memory for.
   * @param gs_algorithm Specifies which algorithm to use:
   *
   *                     KokkosSpace::GS_DEFAULT  PointGaussSeidel
   *                     KokkosSpace::GS_PERMUTED ??
   *                     KokkosSpace::GS_TEAM     ??
   *                     KokkosSpace::GS_CLUSTER  ??
   *                     KokkosSpace::GS_TWOSTAGE ??
   * @param coloring_algorithm Specifies which coloring algorithm to color the graph with:
   *
   *                           KokkosGraph::COLORING_DEFAULT  ??
   *                           KokkosGraph::COLORING_SERIAL   Serial Greedy Coloring
   *                           KokkosGraph::COLORING_VB       Vertex Based Coloring
   *                           KokkosGraph::COLORING_VBBIT    Vertex Based Coloring with bit array
   *                           KokkosGraph::COLORING_VBCS     Vertex Based Color Set
   *                           KokkosGraph::COLORING_VBD      Vertex Based Deterministic Coloring
   *                           KokkosGraph::COLORING_VBDBIT   Vertex Based Deterministic Coloring with bit array
   *                           KokkosGraph::COLORING_EB       Edge Based Coloring
   *                           KokkosGraph::COLORING_SERIAL2  Serial Distance-2 Graph Coloring (kept here for
   *                                                           backwards compatibility for SPGEMM and other use cases)
   */
  // clang-format on
  void create_gs_handle(const HandleExecSpace &handle_exec_space, int num_streams,
                        KokkosSparse::GSAlgorithm gs_algorithm            = KokkosSparse::GS_DEFAULT,
                        KokkosGraph::ColoringAlgorithm coloring_algorithm = KokkosGraph::COLORING_DEFAULT) {
    this->destroy_gs_handle();
    this->is_owner_of_the_gs_handle = true;
    // ---------------------------------------- //
    // Two-stage Gauss-Seidel
    if (gs_algorithm == KokkosSparse::GS_TWOSTAGE)
      this->gsHandle = new TwoStageGaussSeidelHandleType(handle_exec_space, num_streams);
    else
      this->gsHandle = new PointGaussSeidelHandleType(handle_exec_space, num_streams, gs_algorithm, coloring_algorithm);
  }

  // clang-format off
  /**
   * @brief Create a gauss seidel handle object
   *
   * @param gs_algorithm Specifies which algorithm to use:
   *
   *                     KokkosSpace::GS_DEFAULT  PointGaussSeidel or BlockGaussSeidel, depending on matrix type.
   *                     KokkosSpace::GS_PERMUTED Reorders rows/cols into colors to improve locality. Uses RangePolicy over rows.
   *                     KokkosSpace::GS_TEAM     Uses TeamPolicy over batches of rows with ThreadVector within rows.
   *                     KokkosSpace::GS_CLUSTER  Uses independent clusters of nodes in the graph. Within a cluster, x is updated sequentially.
   *                                              For more information, see: https://arxiv.org/pdf/2204.02934.pdf.
   *                     KokkosSpace::GS_TWOSTAGE Uses spmv to parallelize inner sweeps of x.
   *                                              For more information, see: https://arxiv.org/pdf/2104.01196.pdf.
   * @param coloring_algorithm Specifies which coloring algorithm to color the graph with:
   *
   *                           KokkosGraph::COLORING_DEFAULT  Depends on execution space:
   *                                                            COLORING_SERIAL on Kokkos::Serial;
   *                                                            COLORING_EB on GPUs;
   *                                                            COLORING_VBBIT on Kokkos::Sycl or elsewhere.
   *                           KokkosGraph::COLORING_SERIAL   Serial Greedy Coloring
   *                           KokkosGraph::COLORING_VB       Vertex Based Coloring
   *                           KokkosGraph::COLORING_VBBIT    Vertex Based Coloring with bit array
   *                           KokkosGraph::COLORING_VBCS     Vertex Based Color Set
   *                           KokkosGraph::COLORING_VBD      Vertex Based Deterministic Coloring
   *                           KokkosGraph::COLORING_VBDBIT   Vertex Based Deterministic Coloring with bit array
   *                           KokkosGraph::COLORING_EB       Edge Based Coloring
   *                           KokkosGraph::COLORING_SERIAL2  Serial Distance-2 Graph Coloring (kept here for
   *                                                           backwards compatibility for SPGEMM and other use cases)
   */
  // clang-format on
  void create_gs_handle(KokkosSparse::GSAlgorithm gs_algorithm            = KokkosSparse::GS_DEFAULT,
                        KokkosGraph::ColoringAlgorithm coloring_algorithm = KokkosGraph::COLORING_DEFAULT) {
    HandleExecSpace handle_exec_space;
    return create_gs_handle(handle_exec_space, 1, gs_algorithm, coloring_algorithm);
  }
  // ---------------------------------------- //
  // Two-stage Gauss-Seidel handle
  TwoStageGaussSeidelHandleType *get_twostage_gs_handle() {
    auto gs2 = dynamic_cast<TwoStageGaussSeidelHandleType *>(this->gsHandle);
    if (this->gsHandle && !gs2)
      throw std::runtime_error("GaussSeidelHandle exists but is not set up for two-stage GS.");
    return gs2;
  }
  // ---------------------------------------- //
  // Specify numer of outer sweeps for two-stage Gauss-Seidel
  void set_gs_set_num_outer_sweeps(int num_outer_sweeps) {
    auto gs2 = get_twostage_gs_handle();
    gs2->setNumOuterSweeps(num_outer_sweeps);
  }
  // ---------------------------------------- //
  // Specify numer of inner sweeps for two-stage Gauss-Seidel
  void set_gs_set_num_inner_sweeps(int num_inner_sweeps) {
    auto gs2 = get_twostage_gs_handle();
    gs2->setNumInnerSweeps(num_inner_sweeps);
  }
  // ---------------------------------------- //
  // Specify damping factor of inner sweeps for two-stage Gauss-Seidel
  void set_gs_set_inner_damp_factor(nnz_scalar_t damp_factor) {
    auto gs2 = get_twostage_gs_handle();
    gs2->setInnerDampFactor(damp_factor);
  }
  // ---------------------------------------- //
  // Specify to use either Two-stage or Classical (i.e., inner Jacobi-Richardson
  // or SpTrsv)
  void set_gs_twostage(bool two_stage, size_type nrows) {
    auto gs2 = get_twostage_gs_handle();
    gs2->setTwoStage(two_stage);
    if (!two_stage) {
      using namespace KokkosSparse::Experimental;
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE)
      // NOTE: we call CuSPARSE on GPU, if possible
      if (std::is_same<size_type, int>::value && std::is_same<nnz_lno_t, int>::value &&
          std::is_same<HandleExecSpace, Kokkos::Cuda>::value) {
        this->create_gs_sptrsvL_handle(SPTRSVAlgorithm::SPTRSV_CUSPARSE, nrows);
        this->create_gs_sptrsvU_handle(SPTRSVAlgorithm::SPTRSV_CUSPARSE, nrows);
      } else
#endif
      {
        this->create_gs_sptrsvL_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows);
        this->create_gs_sptrsvU_handle(SPTRSVAlgorithm::SEQLVLSCHD_TP1, nrows);
      }
    }
  }
  // ---------------------------------------- //
  // Specify to use either Compact or Classical form of recurrence
  void set_gs_twostage_compact_form(bool compact_form) {
    auto gs2 = get_twostage_gs_handle();
    gs2->setCompactForm(compact_form);
  }

  // clang-format off
  /**
   * @brief Create a gs handle object
   *
   * @param clusterAlgo Specifies which clustering algorithm to use:
   *
   *                    KokkosSparse::CLUSTER_DEFAULT           ??
   *                    KokkosSparse::CLUSTER_MIS2              ??
   *                    KokkosSparse::CLUSTER_BALLOON           ??
   *                    KokkosSparse::NUM_CLUSTERING_ALGORITHMS ??
   * @param hint_verts_per_cluster Hint how many verticies to use per cluster
   * @param coloring_algorithm Specifies which coloring algorithm to color the graph with:
   *
   *                           KokkosGraph::COLORING_DEFAULT  ??
   *                           KokkosGraph::COLORING_SERIAL   Serial Greedy Coloring
   *                           KokkosGraph::COLORING_VB       Vertex Based Coloring
   *                           KokkosGraph::COLORING_VBBIT    Vertex Based Coloring with bit array
   *                           KokkosGraph::COLORING_VBCS     Vertex Based Color Set
   *                           KokkosGraph::COLORING_VBD      Vertex Based Deterministic Coloring
   *                           KokkosGraph::COLORING_VBDBIT   Vertex Based Deterministic Coloring with bit array
   *                           KokkosGraph::COLORING_EB       Edge Based Coloring
   *                           KokkosGraph::COLORING_SERIAL2  Serial Distance-2 Graph Coloring (kept here for
   *                                                           backwards compatibility for SPGEMM and other use cases)
   */
  // clang-format on
  void create_gs_handle(KokkosSparse::ClusteringAlgorithm clusterAlgo, nnz_lno_t hint_verts_per_cluster,
                        KokkosGraph::ColoringAlgorithm coloring_algorithm = KokkosGraph::COLORING_DEFAULT) {
    this->destroy_gs_handle();
    this->is_owner_of_the_gs_handle = true;
    this->gsHandle = new ClusterGaussSeidelHandleType(clusterAlgo, hint_verts_per_cluster, coloring_algorithm);
  }
  void destroy_gs_handle() {
    if (is_owner_of_the_gs_handle && this->gsHandle != NULL) {
      delete this->gsHandle;
      this->gsHandle = NULL;
    }
  }

  // ---------------------------------------- //
  // Handles for Classical GS (inner SpTRSV)
  TwoStageGaussSeidelSPTRSVHandleType *get_gs_sptrsvL_handle() { return this->gs_sptrsvLHandle; }
  TwoStageGaussSeidelSPTRSVHandleType *get_gs_sptrsvU_handle() { return this->gs_sptrsvUHandle; }
  void create_gs_sptrsvL_handle(KokkosSparse::Experimental::SPTRSVAlgorithm algm, size_type nrows) {
    this->destroy_gs_sptrsvL_handle();
    this->is_owner_of_the_gs_sptrsvL_handle = true;
    this->gs_sptrsvLHandle                  = new TwoStageGaussSeidelSPTRSVHandleType();
    this->gs_sptrsvLHandle->create_sptrsv_handle(algm, nrows, true);
  }
  void create_gs_sptrsvU_handle(KokkosSparse::Experimental::SPTRSVAlgorithm algm, size_type nrows) {
    this->destroy_gs_sptrsvU_handle();
    this->is_owner_of_the_gs_sptrsvU_handle = true;
    this->gs_sptrsvUHandle                  = new TwoStageGaussSeidelSPTRSVHandleType();
    this->gs_sptrsvUHandle->create_sptrsv_handle(algm, nrows, false);
  }
  void destroy_gs_sptrsvL_handle() {
    if (this->is_owner_of_the_gs_sptrsvL_handle && this->gs_sptrsvLHandle != nullptr) {
      delete this->gs_sptrsvLHandle;
      this->gs_sptrsvLHandle = nullptr;
    }
  }
  void destroy_gs_sptrsvU_handle() {
    if (this->is_owner_of_the_gs_sptrsvU_handle && this->gs_sptrsvUHandle != nullptr) {
      delete this->gs_sptrsvUHandle;
      this->gs_sptrsvUHandle = nullptr;
    }
  }
  // ---------------------------------------- //

  SPADDHandleType *get_spadd_handle() { return this->spaddHandle; }
  void create_spadd_handle(bool input_sorted = false, bool input_merged = false) {
    this->destroy_spadd_handle();
    this->is_owner_of_the_spadd_handle = true;
    this->spaddHandle                  = new SPADDHandleType(input_sorted, input_merged);
  }
  void destroy_spadd_handle() {
    if (is_owner_of_the_spadd_handle && this->spaddHandle != NULL) {
      delete this->spaddHandle;
      this->spaddHandle = NULL;
    }
  }

  SPTRSVHandleType *get_sptrsv_handle() { return this->sptrsvHandle; }

  void create_sptrsv_handle(KokkosSparse::Experimental::SPTRSVAlgorithm algm, size_type nrows, bool lower_tri,
                            size_type block_size = 0) {
    this->destroy_sptrsv_handle();
    this->is_owner_of_the_sptrsv_handle = true;
    this->sptrsvHandle                  = new SPTRSVHandleType(algm, nrows, lower_tri, block_size);
    //    this->sptrsvHandle->init_handle(nrows);
    this->sptrsvHandle->set_team_size(this->team_work_size);
    this->sptrsvHandle->set_vector_size(this->vector_size);

#ifdef KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV
    // default SpMV option
    using namespace KokkosSparse::Experimental;
    if (algm == SPTRSVAlgorithm::SUPERNODAL_SPMV || algm == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG) {
      this->set_sptrsv_column_major(true);
    }
#endif
  }

#ifdef KOKKOSKERNELS_ENABLE_SUPERNODAL_SPTRSV
  void set_sptrsv_verbose(bool verbose) { this->sptrsvHandle->set_verbose(verbose); }

  void set_sptrsv_perm(int *perm) { this->sptrsvHandle->set_perm(perm); }

  void set_sptrsv_supernodes(int nsuper, integer_view_host_t supercols, int *etree) {
    this->sptrsvHandle->set_supernodes(nsuper, supercols, etree);
  }

  void set_sptrsv_diag_supernode_sizes(int unblocked, int blocked) {
    this->sptrsvHandle->set_supernode_size_unblocked(unblocked);
    this->sptrsvHandle->set_supernode_size_blocked(blocked);
  }

  void set_sptrsv_unit_diagonal(bool flag) { this->sptrsvHandle->set_unit_diagonal(flag); }

  void set_sptrsv_merge_supernodes(bool flag) { this->sptrsvHandle->set_merge_supernodes(flag); }

  void set_sptrsv_invert_diagonal(bool flag) {
    auto algo = this->sptrsvHandle->get_algorithm();
    using namespace KokkosSparse::Experimental;
    if (!flag && (algo == SPTRSVAlgorithm::SUPERNODAL_SPMV || algo == SPTRSVAlgorithm::SUPERNODAL_SPMV_DAG)) {
      std::cout << std::endl
                << " ** Supernodal SpTRSV with SpMV require diagonal inversion **" << std::endl
                << std::endl;
      this->sptrsvHandle->set_invert_diagonal(true);
      return;
    }
    this->sptrsvHandle->set_invert_diagonal(flag);
  }

  bool get_sptrsv_invert_diagonal() { return this->sptrsvHandle->get_invert_diagonal(); }

  void set_sptrsv_invert_offdiagonal(bool flag) {
    if (flag == true && !(this->is_sptrsv_column_major())) {
      std::cout << std::endl << " ** cannot invert offdiagonal in CSR **" << std::endl << std::endl;
      this->sptrsvHandle->set_invert_offdiagonal(false);
      return;
    }

    this->sptrsvHandle->set_invert_offdiagonal(flag);
  }

  bool get_sptrsv_invert_offdiagonal() { return this->sptrsvHandle->get_invert_offdiagonal(); }

  void set_sptrsv_etree(int *etree) { this->sptrsvHandle->set_etree(etree); }

  void set_sptrsv_column_major(bool col_major) {
    if (col_major == false && this->sptrsvHandle->get_invert_offdiagonal()) {
      std::cout << std::endl << " ** cannot use CSR for inverting offdiagonal **" << std::endl << std::endl;
      return;
    }
    this->sptrsvHandle->set_column_major(col_major);
  }

  bool is_sptrsv_lower_tri() { return this->sptrsvHandle->is_lower_tri(); }

  bool is_sptrsv_column_major() { return this->sptrsvHandle->is_column_major(); }

  void set_sptrsv_trmm_on_device(bool trmm_on_device) { this->sptrsvHandle->set_trmm_on_device(trmm_on_device); }
#endif
  void destroy_sptrsv_handle() {
    if (is_owner_of_the_sptrsv_handle && this->sptrsvHandle != nullptr) {
      delete this->sptrsvHandle;
      this->sptrsvHandle = nullptr;
    }
  }

  SPILUKHandleType *get_spiluk_handle() { return this->spilukHandle; }
  void create_spiluk_handle(KokkosSparse::Experimental::SPILUKAlgorithm algm, size_type nrows, size_type nnzL,
                            size_type nnzU, size_type block_size = 0) {
    this->destroy_spiluk_handle();
    this->is_owner_of_the_spiluk_handle = true;
    this->spilukHandle                  = new SPILUKHandleType(algm, nrows, nnzL, nnzU, block_size);
    this->spilukHandle->reset_handle(nrows, nnzL, nnzU, block_size);
    this->spilukHandle->set_team_size(this->team_work_size);
    this->spilukHandle->set_vector_size(this->vector_size);
  }
  void destroy_spiluk_handle() {
    if (is_owner_of_the_spiluk_handle && this->spilukHandle != nullptr) {
      delete this->spilukHandle;
      this->spilukHandle = nullptr;
    }
  }

  PAR_ILUTHandleType *get_par_ilut_handle() { return this->par_ilutHandle; }
  void create_par_ilut_handle(const size_type max_iter                                            = 20,
                              const typename PAR_ILUTHandleType::float_t residual_norm_delta_stop = 1e-2,
                              const typename PAR_ILUTHandleType::float_t fill_in_limit            = 0.75,
                              const bool async_update = false, const bool verbose = false) {
    this->destroy_par_ilut_handle();
    this->is_owner_of_the_par_ilut_handle = true;
    this->par_ilutHandle =
        new PAR_ILUTHandleType(max_iter, residual_norm_delta_stop, fill_in_limit, async_update, verbose);
    this->par_ilutHandle->set_team_size(this->team_work_size);
    this->par_ilutHandle->set_vector_size(this->vector_size);
  }
  void destroy_par_ilut_handle() {
    if (is_owner_of_the_par_ilut_handle && this->par_ilutHandle != nullptr) {
      delete this->par_ilutHandle;
      this->par_ilutHandle = nullptr;
    }
  }

  GMRESHandleType *get_gmres_handle() { return this->gmresHandle; }
  void create_gmres_handle(const size_type m = 50, const typename GMRESHandleType::float_t tol = 1e-8,
                           const size_type max_restart = 50) {
    this->destroy_gmres_handle();
    this->is_owner_of_the_gmres_handle = true;
    this->gmresHandle                  = new GMRESHandleType(m, tol, max_restart);
  }
  void destroy_gmres_handle() {
    if (is_owner_of_the_gmres_handle && this->gmresHandle != nullptr) {
      delete this->gmresHandle;
      this->gmresHandle = nullptr;
    }
  }

};  // end class KokkosKernelsHandle

}  // namespace Experimental
}  // namespace KokkosKernels

#endif  //_KOKKOSKERNELHANDLE_HPP
