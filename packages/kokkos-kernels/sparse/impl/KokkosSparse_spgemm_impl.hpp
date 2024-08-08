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

#ifndef _KOKKOSSPGEMMIMPL_HPP
#define _KOKKOSSPGEMMIMPL_HPP

#include <KokkosKernels_Utils.hpp>
#include <KokkosKernels_SimpleUtils.hpp>
#include <KokkosSparse_Utils.hpp>
#include <KokkosKernels_VectorUtils.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "KokkosKernels_HashmapAccumulator.hpp"
#include "KokkosKernels_Uniform_Initialized_MemoryPool.hpp"
#include "KokkosSparse_spgemm_handle.hpp"
#include "KokkosGraph_Distance1Color.hpp"

namespace KokkosSparse {

namespace Impl {

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
class KokkosSPGEMM {
 public:
  typedef a_row_view_t_ a_row_view_t;
  typedef a_lno_nnz_view_t_ a_in_lno_nnz_view_t;
  typedef a_scalar_nnz_view_t_ a_in_scalar_nnz_view_t;

  typedef b_lno_row_view_t_ b_in_lno_row_view_t;
  typedef b_lno_nnz_view_t_ b_in_lno_nnz_view_t;
  typedef b_scalar_nnz_view_t_ b_in_scalar_nnz_view_t;

  typedef typename a_row_view_t::non_const_value_type size_type;
  typedef typename a_row_view_t::const_value_type const_size_type;

  typedef typename a_in_lno_nnz_view_t::non_const_value_type nnz_lno_t;
  typedef typename a_in_lno_nnz_view_t::const_value_type const_nnz_lno_t;

  typedef typename a_in_scalar_nnz_view_t::non_const_value_type scalar_t;
  typedef typename a_in_scalar_nnz_view_t::const_value_type const_scalar_t;

  typedef typename a_row_view_t::const_type const_a_lno_row_view_t;
  typedef typename a_row_view_t::non_const_type non_const_a_lno_row_view_t;

  typedef typename a_in_lno_nnz_view_t::const_type const_a_lno_nnz_view_t;
  typedef typename a_in_lno_nnz_view_t::non_const_type non_const_a_lno_nnz_view_t;

  typedef typename a_in_scalar_nnz_view_t::const_type const_a_scalar_nnz_view_t;
  typedef typename a_in_scalar_nnz_view_t::non_const_type non_const_a_scalar_nnz_view_t;

  typedef typename b_in_lno_row_view_t::const_type const_b_lno_row_view_t;
  typedef typename b_in_lno_row_view_t::non_const_type non_const_b_lno_row_view_t;

  typedef typename b_in_lno_nnz_view_t::const_type const_b_lno_nnz_view_t;
  typedef typename b_in_lno_nnz_view_t::non_const_type non_const_b_lno_nnz_view_t;

  typedef typename b_in_scalar_nnz_view_t::const_type const_b_scalar_nnz_view_t;
  typedef typename b_in_scalar_nnz_view_t::non_const_type non_const_b_scalar_nnz_view_t;

  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;

  typedef typename HandleType::row_lno_temp_work_view_t row_lno_temp_work_view_t;
  typedef typename HandleType::row_lno_persistent_work_view_t row_lno_persistent_work_view_t;
  typedef
      typename HandleType::row_lno_persistent_work_host_view_t row_lno_persistent_work_host_view_t;  // Host view type

  typedef typename HandleType::nnz_lno_temp_work_view_t nnz_lno_temp_work_view_t;
  typedef typename HandleType::nnz_lno_persistent_work_view_t nnz_lno_persistent_work_view_t;
  typedef
      typename HandleType::nnz_lno_persistent_work_host_view_t nnz_lno_persistent_work_host_view_t;  // Host view type

  typedef typename HandleType::scalar_temp_work_view_t scalar_temp_work_view_t;
  typedef typename HandleType::scalar_persistent_work_view_t scalar_persistent_work_view_t;

  typedef typename HandleType::bool_persistent_view_t bool_persistent_view_t;
  typedef typename HandleType::bool_temp_view_t bool_temp_view_t;

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy_t;
  typedef typename team_policy_t::member_type team_member_t;

  struct CountTag {};
  struct GPUCountTag {};
  struct CountTag2 {};

  struct FillTag {};
  struct FillTag2 {};
  struct MultiCoreDenseAccumulatorTag {};
  struct MultiCoreDenseAccumulatorTag2 {};
  struct MultiCoreDenseAccumulatorTag3 {};
  struct NoCompressMultiCoreDenseAccumulatorTag {};
  struct NoCompressMultiCoreDenseAccumulatorTag2 {};
  struct NoCompressMultiCoreDenseAccumulatorTag3 {};
  struct MultiCoreTag {};
  struct MultiCoreTag2 {};
  struct MultiCoreTag3 {};
  struct MultiCoreTag4 {};
  struct MultiCoreTag5 {};
  struct MultiCoreTag6 {};
  struct GPUTag {};
  struct GPUTag2 {};
  struct GPUTag3 {};
  struct GPUTag4 {};
  struct GPUTag5 {};
  struct GPUTag6 {};

  struct Numeric1Tag {};
  struct Numeric2Tag {};
  struct Numeric3Tag {};

  typedef Kokkos::TeamPolicy<MultiCoreDenseAccumulatorTag, MyExecSpace> multicore_dense_team_count_policy_t;
  typedef Kokkos::TeamPolicy<MultiCoreDenseAccumulatorTag2, MyExecSpace> multicore_dense_team2_count_policy_t;
  typedef Kokkos::TeamPolicy<MultiCoreDenseAccumulatorTag3, MyExecSpace> multicore_dense_team3_count_policy_t;

  typedef Kokkos::TeamPolicy<NoCompressMultiCoreDenseAccumulatorTag, MyExecSpace>
      nc_multicore_dense_team_count_policy_t;
  typedef Kokkos::TeamPolicy<NoCompressMultiCoreDenseAccumulatorTag2, MyExecSpace>
      nc_multicore_dense_team2_count_policy_t;
  typedef Kokkos::TeamPolicy<NoCompressMultiCoreDenseAccumulatorTag3, MyExecSpace>
      nc_multicore_dense_team3_count_policy_t;

  typedef Kokkos::TeamPolicy<NoCompressMultiCoreDenseAccumulatorTag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
      nc_dynamic_multicore_dense_team_count_policy_t;
  typedef Kokkos::TeamPolicy<NoCompressMultiCoreDenseAccumulatorTag2, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
      nc_dynamic_multicore_dense_team2_count_policy_t;
  typedef Kokkos::TeamPolicy<NoCompressMultiCoreDenseAccumulatorTag3, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
      nc_dynamic_multicore_dense_team3_count_policy_t;

  typedef Kokkos::TeamPolicy<MultiCoreTag, MyExecSpace> multicore_team_policy_t;
  typedef Kokkos::TeamPolicy<MultiCoreTag2, MyExecSpace> multicore_team_policy2_t;
  typedef Kokkos::TeamPolicy<MultiCoreTag3, MyExecSpace> multicore_team_policy3_t;
  typedef Kokkos::TeamPolicy<MultiCoreTag4, MyExecSpace> multicore_team_policy4_t;
  typedef Kokkos::TeamPolicy<MultiCoreTag5, MyExecSpace> multicore_team_policy5_t;
  typedef Kokkos::TeamPolicy<MultiCoreTag6, MyExecSpace> multicore_team_policy6_t;

  typedef Kokkos::TeamPolicy<GPUTag, MyExecSpace> gpu_team_policy_t;
  typedef Kokkos::TeamPolicy<GPUTag2, MyExecSpace> gpu_team_policy2_t;
  typedef Kokkos::TeamPolicy<GPUTag3, MyExecSpace> gpu_team_policy3_t;
  typedef Kokkos::TeamPolicy<GPUTag4, MyExecSpace> gpu_team_policy4_t;
  typedef Kokkos::TeamPolicy<GPUTag5, MyExecSpace> gpu_team_policy5_t;
  typedef Kokkos::TeamPolicy<GPUTag6, MyExecSpace> gpu_team_policy6_t;

  typedef Kokkos::TeamPolicy<CountTag, MyExecSpace> team_count_policy_t;
  typedef Kokkos::TeamPolicy<CountTag2, MyExecSpace> team_count2_policy_t;

  typedef Kokkos::TeamPolicy<GPUCountTag, MyExecSpace> team_gpucount_policy_t;

  typedef Kokkos::TeamPolicy<FillTag, MyExecSpace> team_fill_policy_t;
  typedef Kokkos::TeamPolicy<FillTag2, MyExecSpace> team_fill2_policy_t;

  typedef Kokkos::TeamPolicy<Numeric1Tag, MyExecSpace> team_numeric1_policy_t;
  typedef Kokkos::TeamPolicy<Numeric2Tag, MyExecSpace> team_numeric2_policy_t;
  typedef Kokkos::TeamPolicy<Numeric3Tag, MyExecSpace> team_numeric3_policy_t;

  typedef Kokkos::TeamPolicy<MultiCoreDenseAccumulatorTag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
      dynamic_multicore_dense_team_count_policy_t;
  typedef Kokkos::TeamPolicy<MultiCoreDenseAccumulatorTag2, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
      dynamic_multicore_dense_team2_count_policy_t;
  typedef Kokkos::TeamPolicy<MultiCoreDenseAccumulatorTag3, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
      dynamic_multicore_dense_team3_count_policy_t;

  typedef Kokkos::TeamPolicy<MultiCoreTag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
      dynamic_multicore_team_policy_t;
  typedef Kokkos::TeamPolicy<MultiCoreTag2, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
      dynamic_multicore_team_policy2_t;
  typedef Kokkos::TeamPolicy<MultiCoreTag3, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
      dynamic_multicore_team_policy3_t;
  typedef Kokkos::TeamPolicy<MultiCoreTag4, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
      dynamic_multicore_team_policy4_t;
  typedef Kokkos::TeamPolicy<MultiCoreTag5, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
      dynamic_multicore_team_policy5_t;
  typedef Kokkos::TeamPolicy<MultiCoreTag6, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
      dynamic_multicore_team_policy6_t;

  typedef Kokkos::TeamPolicy<CountTag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> > dynamic_team_count_policy_t;
  typedef Kokkos::TeamPolicy<FillTag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> > dynamic_team_fill_policy_t;
  typedef Kokkos::TeamPolicy<Numeric1Tag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
      dynamic_team_numeric1_policy_t;
  typedef Kokkos::TeamPolicy<Numeric2Tag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
      dynamic_team_numeric2_policy_t;
  typedef Kokkos::TeamPolicy<Numeric3Tag, MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >
      dynamic_team_numeric3_policy_t;

  typedef Kokkos::TeamPolicy<MyExecSpace, Kokkos::Schedule<Kokkos::Dynamic> > dynamic_team_policy_t;

 protected:
  HandleType *handle;
  nnz_lno_t a_row_cnt;
  nnz_lno_t b_row_cnt;
  nnz_lno_t b_col_cnt;

  const_a_lno_row_view_t row_mapA;
  const_a_lno_nnz_view_t entriesA;
  const_a_scalar_nnz_view_t valsA;
  bool transposeA;

  const_b_lno_row_view_t row_mapB;
  const_b_lno_nnz_view_t entriesB;
  const_b_scalar_nnz_view_t valsB;
  bool transposeB;

  const size_t shmem_size;
  size_t concurrency;
  const bool use_dynamic_schedule;
  const bool KOKKOSKERNELS_VERBOSE;
  // const int KOKKOSKERNELS_VERBOSE = 1;

  const KokkosKernels::Impl::ExecSpaceType MyEnumExecSpace;
  const SPGEMMAlgorithm spgemm_algorithm;
  const SPGEMMAccumulator spgemm_accumulator;

  //////////////////////////////////////////////////////////////////////////////
  //////Function and Struct for matrix compression.
  //////Declerations are at KokkosKernels_SPGEMM_impl_compression.hpp
  //////////////////////////////////////////////////////////////////////////////

  /**
   * \brief Given a symbolic matrix (a graph), it compresses the graph using
   * bits. \param in_row_map: input row pointers. \param in_entries: input
   * column entries \param out_row_map: output row pointers of the compressed
   * matrix \param out_nnz_indices: output, column set indices of the output
   * matrix. \param out_nnz_sets: output, column sets of the output matrix.
   *
   */
  template <typename in_row_view_t, typename in_nnz_view_t, typename out_rowmap_view_t, typename out_nnz_view_t>
  bool compressMatrix(nnz_lno_t n, size_type nnz, in_row_view_t in_row_map, in_nnz_view_t in_entries,
                      out_rowmap_view_t out_row_map, out_nnz_view_t &out_nnz_indices, out_nnz_view_t &out_nnz_sets,
                      bool singleStep);

 public:
  /**
   *\brief Functor to zip the B matrix.
   */
  template <typename row_view_t, typename nnz_view_t, typename new_row_view_t, typename new_nnz_view_t,
            typename pool_memory_space>
  struct SingleStepZipMatrix;

 private:
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  ////BELOW code is for triangle count specific.
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  template <typename struct_visit_t>
  void triangle_count_ai(const int is_symbolic_or_numeric, const nnz_lno_t m, const size_type *row_mapA_,
                         const nnz_lno_t *entriesA_,

                         const size_type bnnz, const size_type *old_row_mapB, const size_type *row_mapB_,
                         const nnz_lno_t *entriesSetIndex, const nnz_lno_t *entriesSets,

                         size_type *rowmapC, nnz_lno_t *entriesC, struct_visit_t visit_applier);

 public:
  template <typename pool_memory_space, typename struct_visit_t>
  struct TriangleCount;

  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosSPGEMM_numeric_triangle(c_row_view_t rowmapC_, c_lno_nnz_view_t entriesC_, c_scalar_nnz_view_t valuesC_);

  template <typename c_row_view_t>
  void KokkosSPGEMM_symbolic_triangle(c_row_view_t rowmapC_);
  template <typename visit_struct_t>
  void KokkosSPGEMM_generic_triangle(visit_struct_t visit_apply);

  /*
  template <typename visit_struct_t>
  void KokkosSPGEMM_generic_triangle_no_compression(visit_struct_t visit_apply);

  template <typename struct_visit_t>
  void triangle_count_ai_no_compression(
          const nnz_lno_t m,
          const size_type* row_mapA_,
          const nnz_lno_t * entriesA_,

          const size_type bnnz,
          const size_type * rowmapB_begins,
          const size_type * rowmapB_ends,
          const nnz_lno_t * entriesB,
          struct_visit_t visit_applier);
  */
  void KokkosSPGEMM_symbolic_triangle_setup();

 private:
  template <typename c_row_view_t, typename c_lno_nnz_view_t>
  void KokkosSPGEMM_numeric_triangle_ai(c_row_view_t rowmapC_, c_lno_nnz_view_t entriesC_);

 public:
  //////////////////////////////////////////////////////////////////////////
  /////BELOW CODE IS TO for SPEED SPGEMM
  ////DECL IS AT _speed.hpp
  //////////////////////////////////////////////////////////////////////////
  template <typename a_row_view_t, typename a_nnz_view_t, typename a_scalar_view_t, typename b_row_view_t,
            typename b_nnz_view_t, typename b_scalar_view_t, typename c_row_view_t, typename c_nnz_view_t,
            typename c_scalar_view_t, typename mpool_type>
  struct NumericCMEM_CPU;

  template <typename a_row_view_t__, typename a_nnz_view_t__, typename a_scalar_view_t__, typename b_row_view_t__,
            typename b_nnz_view_t__, typename b_scalar_view_t__, typename c_row_view_t__, typename c_nnz_view_t__,
            typename c_scalar_view_t__, typename c_nnz_tmp_view_t>
  struct NumericCMEM;

 private:
  /**
   * \brief Numeric phase with speed method
   */
  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosSPGEMM_numeric_speed(c_row_view_t rowmapC_, c_lno_nnz_view_t entriesC_, c_scalar_nnz_view_t valuesC_,
                                  KokkosKernels::Impl::ExecSpaceType my_exec_space);

 public:
  /*
    //////////////////////////////////////////////////////////////////////////
    /////BELOW CODE IS TO for colored SPGEMM
    ////DECL IS AT _color.hpp
    //////////////////////////////////////////////////////////////////////////
    template <typename a_row_view_t__, typename a_nnz_view_t__, typename
    a_scalar_view_t__, typename b_row_view_t__, typename b_nnz_view_t__,
    typename b_scalar_view_t__, typename c_row_view_t__, typename
    c_nnz_view_t__, typename c_scalar_view_t__> struct NumericCCOLOR;
  */
 private:
  /**
   * \brief Numeric phase with speed method
   */
  /*
    template <typename c_row_view_t, typename c_lno_nnz_view_t, typename
    c_scalar_nnz_view_t> void KokkosSPGEMM_numeric_color( c_row_view_t rowmapC_,
        c_lno_nnz_view_t entriesC_,
        c_scalar_nnz_view_t valuesC_,
        SPGEMMAlgorithm spgemm_algorithm);

    template <typename c_row_view_t, typename c_nnz_view_t>
    void d2_color_c_matrix(
        c_row_view_t rowmapC,
        c_nnz_view_t entryIndicesC_,

        nnz_lno_t &original_num_colors,
        nnz_lno_persistent_work_host_view_t &h_color_xadj,
        nnz_lno_persistent_work_view_t &color_adj,
        nnz_lno_persistent_work_view_t &vertex_colors_to_store,

        nnz_lno_t &num_colors_in_one_step,
        nnz_lno_t &num_multi_color_steps,
        SPGEMMAlgorithm spgemm_algorithm);
  */
 public:
  //////////////////////////////////////////////////////////////////////////
  /////BELOW CODE IS TO for kkmem SPGEMM
  ////DECL IS AT _kkmem.hpp
  //////////////////////////////////////////////////////////////////////////
  template <typename a_row_view_t, typename a_nnz_view_t, typename a_scalar_view_t, typename b_row_view_t,
            typename b_nnz_view_t, typename b_scalar_view_t, typename c_row_view_t, typename c_nnz_view_t,
            typename c_scalar_view_t, typename pool_memory_type>
  struct PortableNumericCHASH;

 private:
  // KKMEM only difference is work memory does not use output memory for 2nd
  // level accumulator.
  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosSPGEMM_numeric_hash2(c_row_view_t rowmapC_, c_lno_nnz_view_t entriesC_, c_scalar_nnz_view_t valuesC_,
                                  KokkosKernels::Impl::ExecSpaceType my_exec_space);

  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosSPGEMM_numeric_hash(c_row_view_t rowmapC_, c_lno_nnz_view_t entriesC_, c_scalar_nnz_view_t valuesC_,
                                 KokkosKernels::Impl::ExecSpaceType my_exec_space);
#if defined(KOKKOS_ENABLE_OPENMP)
#ifdef KOKKOSKERNELS_HAVE_OUTER
 public:
  // OUTER PRODUCT CODES
  struct Triplet;

  template <typename a_col_view_t, typename a_nnz_view_t, typename a_scalar_view_t, typename b_row_view_t,
            typename b_nnz_view_t, typename b_scalar_view_t, typename flop_row_view_t>
  struct OuterProduct;

  template <typename a_row_view_t, typename b_row_view_t, typename flop_row_view_t>
  struct FlopsPerRowOuter;

 private:
  template <typename triplet_view_t>
  void sort_triplets(triplet_view_t triplets, size_t num_triplets);

  template <typename host_triplet_view_t>
  void merge_triplets_on_slow_memory(host_triplet_view_t *triplets, size_t num_blocks, size_t overall_size,
                                     host_triplet_view_t output_triplets);

  template <typename triplet_view_t, typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  size_t final_collapse_triplets_omp(triplet_view_t triplets, size_t num_triplets, c_row_view_t &rowmapC_,
                                     c_lno_nnz_view_t &entriesC_, c_scalar_nnz_view_t &valuesC_);

  template <typename triplet_view_t>
  size_t collapse_triplets(triplet_view_t triplets, size_t num_triplets);

  template <typename triplet_view_t>
  size_t collapse_triplets_omp(triplet_view_t triplets, size_t num_triplets, triplet_view_t out_triplets);

#endif
#endif

  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosSPGEMM_numeric_outer(c_row_view_t &rowmapC_, c_lno_nnz_view_t &entriesC_, c_scalar_nnz_view_t &valuesC_,
                                  KokkosKernels::Impl::ExecSpaceType my_exec_space);
  //////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

#ifdef KOKKOSKERNELS_ANALYZE_MEMORYACCESS
  //////////////////////////////////////////////////////////////////////////
  /////BELOW CODE IS TO CALCULATE MEMORY ACCESSES WITH HYPERGRAPH MODEL/////
  ////DECL IS AT _memaccess.hpp
  //////////////////////////////////////////////////////////////////////////
 public:
  // Functor to calculate how many flops is performed per row of C.
  template <typename a_row_view_t, typename a_nnz_view_t, typename b_row_view_t, typename b_nnz_view_t,
            typename c_row_view_t>
  struct FlopsPerRow;
  struct Cache;

 private:
  void create_read_write_hg(size_t &overall_flops, row_lno_temp_work_view_t &c_flop_rowmap,
                            row_lno_temp_work_view_t &c_comp_a_net_index, row_lno_temp_work_view_t &c_comp_b_net_index,
                            nnz_lno_temp_work_view_t &c_comp_row_index, nnz_lno_temp_work_view_t &c_comp_col_index);

  template <typename c_row_view_t>
  void print_read_write_cost(c_row_view_t rowmapC);

  template <typename c_row_view_t>
  void read_write_cost(nnz_lno_t num_colors, nnz_lno_t num_multi_colors, nnz_lno_t num_parallel_colors, bool isGPU,
                       int num_cores,

                       nnz_lno_t num_hyperthreads_in_core, nnz_lno_t hyper_threads_in_team,

                       int vectorlane, const int cache_line_size, const int data_size, const int cache_size,

                       nnz_lno_persistent_work_host_view_t color_xadj,
                       typename nnz_lno_persistent_work_view_t::HostMirror color_adj,
                       typename nnz_lno_persistent_work_view_t::HostMirror vertex_colors,

                       size_t overall_flops, typename row_lno_temp_work_view_t::HostMirror c_flop_rowmap,
                       typename row_lno_temp_work_view_t::HostMirror c_comp_a_net_index,
                       typename row_lno_temp_work_view_t::HostMirror c_comp_b_net_index,
                       typename nnz_lno_temp_work_view_t::HostMirror c_comp_row_index,
                       typename nnz_lno_temp_work_view_t::HostMirror c_comp_col_index, c_row_view_t rowmapC,
                       int write_type  // 0 -- KKMEM, 1-KKSPEED, 2- KKCOLOR 3-KKMULTICOLOR
                                       // 4-KKMULTICOLOR2
  );

#endif

 public:
  //////////////////////////////////////////////////////////////////////////
  /////BELOW CODE IS for public symbolic and numeric functions
  ////DECL IS AT _def.hpp
  //////////////////////////////////////////////////////////////////////////
  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosSPGEMM_numeric(c_row_view_t &rowmapC_, c_lno_nnz_view_t &entriesC_, c_scalar_nnz_view_t &valuesC_);
  // TODO: These are references only for outer product algorithm.
  // If the algorithm is removed, then remove the references.

  /**
   * \brief Symbolic phase of the SPGEMM.
   * \param rowmapC_: row pointers for the result matrix. Allocated before the
   * call with size (n+1), where n is the number of rows of first matrix.
   */
  template <typename c_row_view_t>
  void KokkosSPGEMM_symbolic(c_row_view_t rowmapC_);

  template <typename c_row_view_t, typename c_nnz_view_t>
  void write_matrix_to_plot(nnz_lno_t &num_colors, nnz_lno_persistent_work_host_view_t &h_color_xadj,
                            nnz_lno_persistent_work_view_t &color_adj, c_row_view_t &rowmapC,
                            c_nnz_view_t &entryIndicesC_);

  KokkosSPGEMM(HandleType *handle_, nnz_lno_t m_, nnz_lno_t n_, nnz_lno_t k_, const_a_lno_row_view_t row_mapA_,
               const_a_lno_nnz_view_t entriesA_, bool transposeA_, const_b_lno_row_view_t row_mapB_,
               const_b_lno_nnz_view_t entriesB_, bool transposeB_)
      : handle(handle_),
        a_row_cnt(m_),
        b_row_cnt(n_),
        b_col_cnt(k_),
        row_mapA(row_mapA_),
        entriesA(entriesA_),
        valsA(),
        transposeA(transposeA_),
        row_mapB(row_mapB_),
        entriesB(entriesB_),
        valsB(),
        transposeB(transposeB_),
        shmem_size(handle_->get_shmem_size()),
        concurrency(MyExecSpace().concurrency()),
        use_dynamic_schedule(handle_->is_dynamic_scheduling()),
        KOKKOSKERNELS_VERBOSE(handle_->get_verbose()),
        MyEnumExecSpace(this->handle->get_handle_exec_space()),
        spgemm_algorithm(this->handle->get_spgemm_handle()->get_algorithm_type()),
        spgemm_accumulator(this->handle->get_spgemm_handle()->get_accumulator_type())
  //,row_mapC(), entriesC(), valsC()
  {}

  KokkosSPGEMM(HandleType *handle_, nnz_lno_t m_, nnz_lno_t n_, nnz_lno_t k_, const_a_lno_row_view_t row_mapA_,
               const_a_lno_nnz_view_t entriesA_, const_a_scalar_nnz_view_t valsA_, bool transposeA_,
               const_b_lno_row_view_t row_mapB_, const_b_lno_nnz_view_t entriesB_, const_b_scalar_nnz_view_t valsB_,
               bool transposeB_)
      : handle(handle_),
        a_row_cnt(m_),
        b_row_cnt(n_),
        b_col_cnt(k_),
        row_mapA(row_mapA_),
        entriesA(entriesA_),
        valsA(valsA_),
        transposeA(transposeA_),
        row_mapB(row_mapB_),
        entriesB(entriesB_),
        valsB(valsB_),
        transposeB(transposeB_),
        shmem_size(handle_->get_shmem_size()),
        concurrency(MyExecSpace().concurrency()),
        use_dynamic_schedule(handle_->is_dynamic_scheduling()),
        KOKKOSKERNELS_VERBOSE(handle_->get_verbose()),
        MyEnumExecSpace(this->handle->get_handle_exec_space()),
        spgemm_algorithm(this->handle->get_spgemm_handle()->get_algorithm_type()),
        spgemm_accumulator(this->handle->get_spgemm_handle()->get_accumulator_type())
  //,row_mapB(), entriesC(), valsC()
  {}

  //////////////////////////////////////////////////////////////////////////
  /////BELOW CODE IS for symbolic phase
  ////DECL IS AT _symbolic.hpp
  //////////////////////////////////////////////////////////////////////////
 public:
  /***
   * \brief Functor to calculate the row sizes of C.
   */
  template <typename a_row_view_t, typename a_nnz_view_t, typename b_original_row_view_t,
            typename b_compressed_row_view_t, typename b_nnz_view_t,
            typename c_row_view_t,  // typename nnz_lno_temp_work_view_t,
            typename pool_memory_space>
  struct StructureC;

  template <typename a_row_view_t, typename a_nnz_view_t, typename b_original_row_view_t,
            typename b_compressed_row_view_t, typename b_nnz_view_t,
            typename c_row_view_t,  // typename nnz_lno_temp_work_view_t,
            typename pool_memory_space>
  struct StructureC_NC;

  template <typename a_row_view_t, typename a_nnz_view_t, typename b_original_row_view_t,
            typename b_compressed_row_view_t, typename b_nnz_view_t, typename c_row_view_t,
            typename nnz_lno_temp_work_view_t, typename pool_memory_space>
  struct NonzeroesC;

  /**
   * \brief Functor to calculate the max flops in a row of SPGEMM.
   *
   */
  template <typename a_row_view_t, typename a_nnz_view_t, typename b_oldrow_view_t, typename b_row_view_t>
  struct PredicMaxRowNNZ;

  struct PredicMaxRowNNZIntersection;
  struct PredicMaxRowNNZ_p;

 private:
  /**
   * \brief function return max flops for a row in the result multiplication.
   * \param m: number of rows in A
   * \param row_mapA: row pointers of A.
   * \param entriesA: column indices of A
   * \param row_pointers_begin_B: beginning of the row indices for B
   * \param row_pointers_end_B: end of the row indices for B
   */
  template <typename a_row_view_t, typename a_nnz_view_t, typename b_oldrow_view_t, typename b_row_view_t>
  size_t getMaxRoughRowNNZ(nnz_lno_t m, a_row_view_t row_mapA_, a_nnz_view_t entriesA_,

                           b_oldrow_view_t row_pointers_begin_B, b_row_view_t row_pointers_end_B,
                           size_type *flops_per_row = NULL);

  size_t getMaxRoughRowNNZ_p(const nnz_lno_t m, const size_type annz, const size_type *row_mapA_,
                             const nnz_lno_t *entriesA_,

                             const size_type *row_pointers_begin_B, const size_type *row_pointers_end_B);

  size_t getMaxRoughRowNNZIntersection_p(const nnz_lno_t m, const size_type annz, const size_type *row_mapA_,
                                         const nnz_lno_t *entriesA_,

                                         const size_type *row_pointers_begin_B, const size_type *row_pointers_end_B,
                                         nnz_lno_t *min_result_row_for_each_row);

  template <typename a_r_view_t, typename a_nnz_view_t, typename b_original_row_view_t,
            typename b_compressed_row_view_t, typename b_nnz_view_t, typename c_row_view_t>
  void symbolic_c(nnz_lno_t m, a_r_view_t row_mapA_, a_nnz_view_t entriesA_,

                  b_original_row_view_t old_row_mapB, b_compressed_row_view_t row_mapB_, b_nnz_view_t entriesSetIndex,
                  b_nnz_view_t entriesSets,

                  c_row_view_t rowmapC, nnz_lno_t maxNumRoughNonzeros);

  template <typename a_r_view_t, typename a_nnz_view_t, typename b_original_row_view_t,
            typename b_compressed_row_view_t, typename b_nnz_view_t, typename c_row_view_t>
  void symbolic_c_no_compression(nnz_lno_t m, a_r_view_t row_mapA_, a_nnz_view_t entriesA_,

                                 b_original_row_view_t b_rowmap_begin, b_compressed_row_view_t b_rowmap_end,
                                 b_nnz_view_t entriesb_, c_row_view_t rowmapC, nnz_lno_t maxNumRoughNonzeros);

  //////////////////////////////////////////////////////////////////////////
  ///// Jacobi-fused SpGEMM declarations
  //////////////////////////////////////////////////////////////////////////
 public:
  template <typename a_row_view_t, typename a_nnz_view_t, typename a_scalar_view_t, typename b_row_view_t,
            typename b_nnz_view_t, typename b_scalar_view_t, typename c_row_view_t, typename c_nnz_view_t,
            typename c_scalar_view_t, typename dinv_view_t, typename pool_memory_type>
  struct JacobiSpGEMMSparseAcc;

  template <typename a_row_view_t, typename a_nnz_view_t, typename a_scalar_view_t, typename b_row_view_t,
            typename b_nnz_view_t, typename b_scalar_view_t, typename c_row_view_t, typename c_nnz_view_t,
            typename c_scalar_view_t, typename dinv_view_t, typename mpool_type>
  struct JacobiSpGEMMDenseAcc;

  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t, typename dinv_view_t>
  void KokkosSPGEMM_jacobi_sparseacc(c_row_view_t rowmapC_, c_lno_nnz_view_t entriesC_, c_scalar_nnz_view_t valuesC_,
                                     typename c_scalar_nnz_view_t::const_value_type omega, dinv_view_t dinv,
                                     KokkosKernels::Impl::ExecSpaceType lcl_my_exec_space);

 protected:
  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t, typename dinv_view_t>
  void KokkosSPGEMM_jacobi_denseacc(c_row_view_t rowmapC_, c_lno_nnz_view_t entriesC_, c_scalar_nnz_view_t valuesC_,
                                    typename c_scalar_nnz_view_t::const_value_type omega, dinv_view_t dinv,
                                    KokkosKernels::Impl::ExecSpaceType my_exec_space);

  // Utility to compute the number of pool chunks for L2 hashmap accumulators.
  // Uses free memory query for accelerators/GPUs but assumes infinite available
  // host memory.
  //
  // chunk_bytes: bytes in each chunk
  // ideal_num_chunks: number of chunks that would give each thread/team its own
  // chunk (no contention)
  template <typename Pool>
  size_t compute_num_pool_chunks(size_t chunk_bytes, size_t ideal_num_chunks) {
    if (!KokkosKernels::Impl::kk_is_gpu_exec_space<typename Pool::execution_space>()) return ideal_num_chunks;
    size_t free_byte, total_byte;
    KokkosKernels::Impl::kk_get_free_total_memory<typename Pool::memory_space>(free_byte, total_byte);
    size_t required_size = ideal_num_chunks * chunk_bytes;
    if (KOKKOSKERNELS_VERBOSE)
      std::cout << "\tmempool required size:" << required_size << " free_byte:" << free_byte
                << " total_byte:" << total_byte << std::endl;
    size_t num_chunks = ideal_num_chunks;
    // If there is not enough memory to safely allocate ideal_num_chunks, use
    // half the free memory, rounded down
    if (required_size > free_byte / 2) {
      num_chunks = (free_byte / 2) / chunk_bytes;
    }
    // then take the largest power of 2 smaller than that
    size_t po2_num_chunks = 1;
    while (po2_num_chunks * 2 < num_chunks) {
      po2_num_chunks *= 2;
    }
    return po2_num_chunks;
  }

 public:
  void compute_row_flops() {
    auto sh                    = this->handle->get_spgemm_handle();
    nnz_lno_t maxNumRoughZeros = 0;
    size_t overall_flops       = 0;
    Kokkos::Timer timer1;
    auto new_row_mapB_begin = Kokkos::subview(row_mapB, std::make_pair(nnz_lno_t(0), b_row_cnt));
    auto new_row_mapB_end   = Kokkos::subview(row_mapB, std::make_pair(nnz_lno_t(1), b_row_cnt + 1));
    row_lno_persistent_work_view_t flops_per_row(Kokkos::view_alloc(Kokkos::WithoutInitializing, "original row flops"),
                                                 a_row_cnt);

    // get maximum row flops.
    maxNumRoughZeros = this->getMaxRoughRowNNZ(a_row_cnt, row_mapA, entriesA, new_row_mapB_begin, new_row_mapB_end,
                                               flops_per_row.data());

    // calculate overal flops.
    KokkosKernels::Impl::kk_reduce_view2<row_lno_persistent_work_view_t, MyExecSpace>(a_row_cnt, flops_per_row,
                                                                                      overall_flops);
    if (KOKKOSKERNELS_VERBOSE) {
      std::cout << "\tOriginal Max Row Flops:" << maxNumRoughZeros << std::endl;
      std::cout << "\tOriginal overall_flops Flops:" << overall_flops << std::endl;
      std::cout << "\ttOriginal Max Row Flop Calc Time:" << timer1.seconds() << std::endl;
    }
    sh->original_max_row_flops = maxNumRoughZeros;
    sh->original_overall_flops = overall_flops;
    sh->row_flops              = flops_per_row;
    sh->set_computed_rowflops();
  }
};

}  // namespace Impl
}  // namespace KokkosSparse
#include "KokkosSparse_spgemm_imp_outer.hpp"
#include "KokkosSparse_spgemm_impl_memaccess.hpp"
#include "KokkosSparse_spgemm_impl_kkmem.hpp"
#include "KokkosSparse_spgemm_impl_speed.hpp"
#include "KokkosSparse_spgemm_impl_compression.hpp"
#include "KokkosSparse_spgemm_impl_def.hpp"
#include "KokkosSparse_spgemm_impl_symbolic.hpp"
#include "KokkosSparse_spgemm_impl_triangle.hpp"
#endif
