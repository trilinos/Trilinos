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

#ifndef _KOKKOSBSPGEMMIMPL_HPP
#define _KOKKOSBSPGEMMIMPL_HPP

#include "KokkosSparse_spgemm_impl.hpp"

namespace KokkosSparse {

namespace Impl {

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
class KokkosBSPGEMM : public KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_,
                                          b_lno_row_view_t_, b_lno_nnz_view_t_, b_scalar_nnz_view_t_> {
 public:
  using Base = KokkosSparse::Impl::KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_,
                                                b_lno_row_view_t_, b_lno_nnz_view_t_, b_scalar_nnz_view_t_>;

#define USE_BASE_TYPE(type) using type = typename Base::type;

  USE_BASE_TYPE(nnz_lno_t)
  USE_BASE_TYPE(scalar_t)
  USE_BASE_TYPE(size_type)
  USE_BASE_TYPE(const_a_lno_row_view_t)
  USE_BASE_TYPE(const_a_lno_nnz_view_t)
  USE_BASE_TYPE(const_a_scalar_nnz_view_t)
  USE_BASE_TYPE(const_b_lno_row_view_t)
  USE_BASE_TYPE(const_b_lno_nnz_view_t)
  USE_BASE_TYPE(const_b_scalar_nnz_view_t)
  USE_BASE_TYPE(row_lno_persistent_work_view_t)
  USE_BASE_TYPE(nnz_lno_temp_work_view_t)
  USE_BASE_TYPE(team_member_t)

  USE_BASE_TYPE(MyExecSpace)
  USE_BASE_TYPE(MyTempMemorySpace)
  USE_BASE_TYPE(MultiCoreTag)
  USE_BASE_TYPE(MultiCoreTag4)
  USE_BASE_TYPE(GPUTag)
  USE_BASE_TYPE(GPUTag4)
  USE_BASE_TYPE(GPUTag6)
  USE_BASE_TYPE(gpu_team_policy_t)
  USE_BASE_TYPE(gpu_team_policy4_t)
  USE_BASE_TYPE(gpu_team_policy6_t)
  USE_BASE_TYPE(dynamic_multicore_team_policy_t)
  USE_BASE_TYPE(dynamic_multicore_team_policy4_t)
  USE_BASE_TYPE(multicore_team_policy_t)
  USE_BASE_TYPE(multicore_team_policy4_t)

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
  void KokkosBSPGEMM_numeric_speed(c_row_view_t rowmapC_, c_lno_nnz_view_t entriesC_, c_scalar_nnz_view_t valuesC_,
                                   KokkosKernels::Impl::ExecSpaceType my_exec_space);

 private:
  // How many extra bytes are needed to align a scalar_t after an array of
  // nnz_lno_t, in the worst case? Incurred once per hashmap, which may be per
  // team or per thread depending on algorithm
  static constexpr size_t scalarAlignPad =
      (alignof(scalar_t) > alignof(nnz_lno_t)) ? (alignof(scalar_t) - alignof(nnz_lno_t)) : 0;

  static constexpr bool exec_gpu = KokkosKernels::Impl::kk_is_gpu_exec_space<MyExecSpace>();

 private:
  nnz_lno_t block_dim;

 public:
  //////////////////////////////////////////////////////////////////////////
  /////BELOW CODE IS TO for kkmem SPGEMM
  ////DECL IS AT _kkmem.hpp
  //////////////////////////////////////////////////////////////////////////
  template <typename a_row_view_t, typename a_nnz_view_t, typename a_scalar_view_t, typename b_row_view_t,
            typename b_nnz_view_t, typename b_scalar_view_t, typename c_row_view_t, typename c_nnz_view_t,
            typename c_scalar_view_t, typename pool_memory_type>
  struct PortableNumericCHASH;

  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosBSPGEMM_numeric_hash(c_row_view_t rowmapC_, c_lno_nnz_view_t entriesC_, c_scalar_nnz_view_t valuesC_,
                                  KokkosKernels::Impl::ExecSpaceType my_exec_space);

 public:
  //////////////////////////////////////////////////////////////////////////
  /////BELOW CODE IS for public symbolic and numeric functions
  ////DECL IS AT _def.hpp
  //////////////////////////////////////////////////////////////////////////
  template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
  void KokkosBSPGEMM_numeric(c_row_view_t &rowmapC_, c_lno_nnz_view_t &entriesC_, c_scalar_nnz_view_t &valuesC_);

  KokkosBSPGEMM(HandleType *handle_, nnz_lno_t m_, nnz_lno_t n_, nnz_lno_t k_, nnz_lno_t block_dim_,
                const_a_lno_row_view_t row_mapA_, const_a_lno_nnz_view_t entriesA_, bool transposeA_,
                const_b_lno_row_view_t row_mapB_, const_b_lno_nnz_view_t entriesB_, bool transposeB_)
      : Base(handle_, m_, n_, k_, row_mapA_, entriesA_, transposeA_, row_mapB_, entriesB_, transposeB_),
        block_dim(block_dim_) {}

  KokkosBSPGEMM(HandleType *handle_, nnz_lno_t m_, nnz_lno_t n_, nnz_lno_t k_, nnz_lno_t block_dim_,
                const_a_lno_row_view_t row_mapA_, const_a_lno_nnz_view_t entriesA_, const_a_scalar_nnz_view_t valsA_,
                bool transposeA_, const_b_lno_row_view_t row_mapB_, const_b_lno_nnz_view_t entriesB_,
                const_b_scalar_nnz_view_t valsB_, bool transposeB_)
      : Base(handle_, m_, n_, k_, row_mapA_, entriesA_, valsA_, transposeA_, row_mapB_, entriesB_, valsB_, transposeB_),
        block_dim(block_dim_) {}
};

}  // namespace Impl
}  // namespace KokkosSparse
#include "KokkosSparse_bspgemm_impl_kkmem.hpp"
#include "KokkosSparse_bspgemm_impl_speed.hpp"
#include "KokkosSparse_bspgemm_impl_def.hpp"
#endif
