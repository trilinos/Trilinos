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
#ifndef _KOKKOS_TRIANGLE_HPP
#define _KOKKOS_TRIANGLE_HPP
#include "KokkosSparse_spgemm_impl.hpp"
#include "KokkosKernels_IOUtils.hpp"
#include "KokkosKernels_Handle.hpp"
namespace KokkosGraph {

namespace Experimental {
/*
template <typename KernelHandle,
typename alno_row_view_t_,
typename alno_nnz_view_t_,
typename blno_row_view_t_,
typename blno_nnz_view_t_,
typename clno_row_view_t_>
void triangle_count(
    KernelHandle *handle,
    typename KernelHandle::nnz_lno_t m,
    typename KernelHandle::nnz_lno_t n,
    typename KernelHandle::nnz_lno_t k,
    alno_row_view_t_ row_mapA,
    alno_nnz_view_t_ entriesA,
    bool transposeA,
    blno_row_view_t_ row_mapB,
    blno_nnz_view_t_ entriesB,
    bool transposeB,
    clno_row_view_t_ row_mapC){
  using namespace KokkosSparse;

  typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
  spgemmHandleType *sh = handle->get_spgemm_handle();
  switch (sh->get_algorithm_type()){
  case SPGEMM_KK_TRIANGLE_LL:
  {
    KokkosSparse::Impl::KokkosSPGEMM
    <KernelHandle,
    alno_row_view_t_, alno_nnz_view_t_, typename
KernelHandle::in_scalar_nnz_view_t, blno_row_view_t_, blno_nnz_view_t_, typename
KernelHandle::in_scalar_nnz_view_t> kspgemm (handle,m,n,k,row_mapA, entriesA,
transposeA, row_mapB, entriesB, transposeB);
    kspgemm.KokkosSPGEMM_symbolic_triangle(row_mapC);
  }
  break;

  case SPGEMM_KK_TRIANGLE_AI:
  case SPGEMM_KK_TRIANGLE_IA:
  case SPGEMM_KK_TRIANGLE_IA_UNION:
  default:
  {
    KokkosSparse::Impl::KokkosSPGEMM
    <KernelHandle,
    alno_row_view_t_, alno_nnz_view_t_, typename
KernelHandle::in_scalar_nnz_view_t, blno_row_view_t_, blno_nnz_view_t_, typename
KernelHandle::in_scalar_nnz_view_t> kspgemm (handle,m,n,k,row_mapA, entriesA,
transposeA, row_mapB, entriesB, transposeB);
    kspgemm.KokkosSPGEMM_symbolic_triangle(row_mapC);
  }
  break;


  }
  sh->set_call_symbolic();

}


template <typename KernelHandle,
typename alno_row_view_t_,
typename alno_nnz_view_t_,
typename blno_row_view_t_,
typename blno_nnz_view_t_,
typename clno_row_view_t_,
typename clno_nnz_view_t_>
void triangle_enumerate(
    KernelHandle *handle,
    typename KernelHandle::nnz_lno_t m,
    typename KernelHandle::nnz_lno_t n,
    typename KernelHandle::nnz_lno_t k,
    alno_row_view_t_ row_mapA,
    alno_nnz_view_t_ entriesA,

    bool transposeA,
    blno_row_view_t_ row_mapB,
    blno_nnz_view_t_ entriesB,
    bool transposeB,
    clno_row_view_t_ row_mapC,
    clno_nnz_view_t_ &entriesC1,
    clno_nnz_view_t_ &entriesC2 = NULL
){
  using namespace KokkosSparse;


  typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
  spgemmHandleType *sh = handle->get_spgemm_handle();
  if (!sh->is_symbolic_called()){
    triangle_count<KernelHandle,
    alno_row_view_t_, alno_nnz_view_t_,
    blno_row_view_t_, blno_nnz_view_t_,
    clno_row_view_t_>(
        handle, m, n, k,
        row_mapA, entriesA, transposeA,
        row_mapB, entriesB, transposeB,
        row_mapC
    );

    typename clno_row_view_t_::value_type c_nnz_size =
handle->get_spgemm_handle()->get_c_nnz(); if (c_nnz_size){ entriesC1 =
clno_nnz_view_t_ (Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC"),
c_nnz_size);
      //entriesC2 = clno_nnz_view_t_
(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC"), c_nnz_size);
    }
  }


  switch (sh->get_algorithm_type()){
  default:

  case SPGEMM_KK_TRIANGLE_AI:
  case SPGEMM_KK_TRIANGLE_IA:
  case SPGEMM_KK_TRIANGLE_IA_UNION:
  {
    KokkosSparse::Impl::KokkosSPGEMM
    <KernelHandle,
    alno_row_view_t_, alno_nnz_view_t_, typename
KernelHandle::in_scalar_nnz_view_t, blno_row_view_t_, blno_nnz_view_t_, typename
KernelHandle::in_scalar_nnz_view_t> kspgemm (handle,m,n,k,row_mapA, entriesA,
transposeA, row_mapB, entriesB, transposeB);
    kspgemm.KokkosSPGEMM_numeric_triangle(row_mapC, entriesC1, entriesC2);
  }
  break;
  }
}
*/

template <typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_, typename blno_row_view_t_,
          typename blno_nnz_view_t_, typename visit_struct_t>
void triangle_generic(KernelHandle *handle, typename KernelHandle::nnz_lno_t m, typename KernelHandle::nnz_lno_t n,
                      typename KernelHandle::nnz_lno_t k, alno_row_view_t_ row_mapA, alno_nnz_view_t_ entriesA,
                      bool transposeA, blno_row_view_t_ row_mapB, blno_nnz_view_t_ entriesB, bool transposeB,
                      visit_struct_t visit_struct) {
  using namespace KokkosSparse;

  typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
  spgemmHandleType *sh = handle->get_spgemm_handle();
  switch (sh->get_algorithm_type()) {
    // case SPGEMM_KK_TRIANGLE_LL:
    case SPGEMM_KK_TRIANGLE_AI:
    case SPGEMM_KK_TRIANGLE_IA:
    case SPGEMM_KK_TRIANGLE_IA_UNION:
    default: {
      KokkosSparse::Impl::KokkosSPGEMM<KernelHandle, alno_row_view_t_, alno_nnz_view_t_,
                                       typename KernelHandle::in_scalar_nnz_view_t, blno_row_view_t_, blno_nnz_view_t_,
                                       typename KernelHandle::in_scalar_nnz_view_t>
          kspgemm(handle, m, n, k, row_mapA, entriesA, transposeA, row_mapB, entriesB, transposeB);
      kspgemm.KokkosSPGEMM_generic_triangle(visit_struct);
    } break;
  }
}

template <typename KernelHandle, typename alno_row_view_t_, typename alno_nnz_view_t_, typename visit_struct_t>
void triangle_generic(KernelHandle *handle, typename KernelHandle::nnz_lno_t m, alno_row_view_t_ row_mapA,
                      alno_nnz_view_t_ entriesA, visit_struct_t visit_struct) {
  typedef typename KernelHandle::nnz_lno_t nnz_lno_t;
  typedef typename KernelHandle::size_type size_type;

  typedef typename KernelHandle::SPGEMMHandleType spgemmHandleType;
  typedef typename KernelHandle::nnz_lno_persistent_work_view_t nnz_lno_persistent_work_view_t;
  typedef typename KernelHandle::row_lno_persistent_work_view_t row_lno_persistent_work_view_t;

  typedef typename KernelHandle::HandleExecSpace ExecutionSpace;

  using namespace KokkosSparse;

  spgemmHandleType *sh = handle->get_spgemm_handle();
  Kokkos::Timer timer1;

  //////SORT BASE ON THE SIZE OF ROWS/////
  int sort_lower_triangle = sh->get_sort_lower_triangular();
  bool should_i_sort      = false;
  if (sort_lower_triangle == 1)
    should_i_sort = true;
  else if (sort_lower_triangle == 2) {
    size_type max_row_size = 0;
    KokkosKernels::Impl::kk_view_reduce_max_row_size<size_type, ExecutionSpace>(m, row_mapA.data(), row_mapA.data() + 1,
                                                                                max_row_size);

    if (max_row_size > 1000) {
      should_i_sort = true;
    }
  }

  if (should_i_sort) {
    if (sh->get_lower_triangular_permutation().data() == NULL) {
      nnz_lno_persistent_work_view_t new_indices(Kokkos::view_alloc(Kokkos::WithoutInitializing, "new_indices"), m);
      int sort_decreasing_order = 1;
      ////If true we place the largest row to top, so that largest row size will
      /// be minimized in lower triangle.
      if (sh->get_algorithm_type() == SPGEMM_KK_TRIANGLE_AI || sh->get_algorithm_type() == SPGEMM_KK_TRIANGLE_LU) {
        sort_decreasing_order = 0;
        // if false we place the largest row to bottom, so that largest column
        // is minimizedin lower triangle.
      } else if (sh->get_algorithm_type() == SPGEMM_KK_TRIANGLE_LL) {
        sort_decreasing_order = 1;
        // if 2, we do an interleaved sort.
      }
      {
        KokkosSparse::Impl::kk_sort_by_row_size<size_type, nnz_lno_t, ExecutionSpace>(
            m, row_mapA.data(), new_indices.data(), sort_decreasing_order, ExecutionSpace().concurrency());
      }
      sh->set_lower_triangular_permutation(new_indices);
    }
  }
  if (handle->get_verbose()) {
    std::cout << "Preprocess Sorting Time:" << timer1.seconds() << std::endl;
  }
  //////SORT BASE ON THE SIZE OF ROWS/////

  /////////CREATE LOWER TRIANGLE///////
  bool create_lower_triangular = sh->get_create_lower_triangular();
  row_lno_persistent_work_view_t lower_triangular_matrix_rowmap;
  nnz_lno_persistent_work_view_t lower_triangular_matrix_entries;
  timer1.reset();
  if (create_lower_triangular || sh->get_algorithm_type() == SPGEMM_KK_TRIANGLE_LL ||
      sh->get_algorithm_type() == SPGEMM_KK_TRIANGLE_LU) {
    sh->get_lower_triangular_matrix(lower_triangular_matrix_rowmap, lower_triangular_matrix_entries);
    if (lower_triangular_matrix_rowmap.data() == NULL || lower_triangular_matrix_entries.data() == NULL) {
      alno_nnz_view_t_ null_values;
      nnz_lno_persistent_work_view_t new_indices = sh->get_lower_triangular_permutation();

      KokkosSparse::Impl::kk_get_lower_triangle<alno_row_view_t_, alno_nnz_view_t_, alno_nnz_view_t_,
                                                row_lno_persistent_work_view_t, nnz_lno_persistent_work_view_t,
                                                alno_nnz_view_t_, nnz_lno_persistent_work_view_t, ExecutionSpace>(
          m, row_mapA, entriesA, null_values, lower_triangular_matrix_rowmap, lower_triangular_matrix_entries,
          null_values, new_indices, handle->is_dynamic_scheduling(),
          handle->get_team_work_size(1, ExecutionSpace().concurrency(), m));

      sh->set_lower_triangular_matrix(lower_triangular_matrix_rowmap, lower_triangular_matrix_entries);
    }
  }
  if (handle->get_verbose()) {
    std::cout << "Preprocess Create Lower Triangular Time:" << timer1.seconds() << std::endl;
  }
  timer1.reset();

  row_lno_persistent_work_view_t upper_triangular_matrix_rowmap;
  nnz_lno_persistent_work_view_t upper_triangular_matrix_entries;
  if (sh->get_algorithm_type() == SPGEMM_KK_TRIANGLE_LU) {
    sh->get_lower_triangular_matrix(lower_triangular_matrix_rowmap, lower_triangular_matrix_entries);
    alno_nnz_view_t_ null_values;
    nnz_lno_persistent_work_view_t new_indices = sh->get_lower_triangular_permutation();

    KokkosSparse::Impl::kk_get_lower_triangle<alno_row_view_t_, alno_nnz_view_t_, alno_nnz_view_t_,
                                              row_lno_persistent_work_view_t, nnz_lno_persistent_work_view_t,
                                              alno_nnz_view_t_, nnz_lno_persistent_work_view_t, ExecutionSpace>(
        m, row_mapA, entriesA, null_values, upper_triangular_matrix_rowmap, upper_triangular_matrix_entries,
        null_values, new_indices, handle->is_dynamic_scheduling(), 4, false);
  }
  if (handle->get_verbose()) {
    std::cout << "Preprocess Create Upper Triangular Time:" << timer1.seconds() << std::endl;
  }

  /////////CREATE LOWER TRIANGLE///////

  ////
  /// CREATE INCIDENCE MATRIX
  ///
  timer1.reset();
  row_lno_persistent_work_view_t incidence_transpose_rowmap;
  nnz_lno_persistent_work_view_t incidence_transpose_entries;

  row_lno_persistent_work_view_t incidence_rowmap;
  nnz_lno_persistent_work_view_t incidence_entries;
  switch (sh->get_algorithm_type()) {
    // IF it is one of below, we perform I^T x (A) or (L).
    // so create the transpose of I.
    case SPGEMM_KK_TRIANGLE_IA_UNION:
    case SPGEMM_KK_TRIANGLE_IA: {
      // these are the algorithms that requires transpose of the incidence
      // matrix.
      sh->get_lower_triangular_matrix(lower_triangular_matrix_rowmap, lower_triangular_matrix_entries);

      if (lower_triangular_matrix_rowmap.data() == NULL || lower_triangular_matrix_entries.data() == NULL) {
        std::cout << "Creating lower triangular A" << std::endl;

        alno_nnz_view_t_ null_values;
        nnz_lno_persistent_work_view_t new_indices = sh->get_lower_triangular_permutation();

        KokkosSparse::Impl::kk_get_lower_triangle<alno_row_view_t_, alno_nnz_view_t_, alno_nnz_view_t_,
                                                  row_lno_persistent_work_view_t, nnz_lno_persistent_work_view_t,
                                                  alno_nnz_view_t_, nnz_lno_persistent_work_view_t, ExecutionSpace>(
            m, row_mapA, entriesA, null_values, lower_triangular_matrix_rowmap, lower_triangular_matrix_entries,
            null_values, new_indices, handle->is_dynamic_scheduling());
      }
      KokkosSparse::Impl::kk_create_incidence_tranpose_matrix_from_lower_triangle<
          row_lno_persistent_work_view_t, nnz_lno_persistent_work_view_t, row_lno_persistent_work_view_t,
          nnz_lno_persistent_work_view_t, ExecutionSpace>(m, lower_triangular_matrix_rowmap,
                                                          lower_triangular_matrix_entries, incidence_transpose_rowmap,
                                                          incidence_transpose_entries, handle->is_dynamic_scheduling());
    } break;

    // IF it is one of below, we perform (A) or (L) x I
    // so create I.
    case SPGEMM_KK_TRIANGLE_AI: {
      // these are the algorithms that requires the incidence matrix.

      KokkosSparse::Impl::kk_create_incidence_matrix_from_original_matrix<
          alno_row_view_t_, alno_nnz_view_t_, row_lno_persistent_work_view_t, nnz_lno_persistent_work_view_t,
          nnz_lno_persistent_work_view_t, ExecutionSpace>(m, row_mapA, entriesA, incidence_rowmap, incidence_entries,
                                                          sh->get_lower_triangular_permutation(),
                                                          handle->is_dynamic_scheduling());
    } break;
    case SPGEMM_KK_TRIANGLE_LU:
    case SPGEMM_KK_TRIANGLE_LL:
    default: {
      break;
    }
  }

  if (handle->get_verbose()) {
    std::cout << "Preprocess Incidence Matrix Create Time:" << timer1.seconds() << std::endl;
  }
  ////
  /// CREATE INCIDENCE MATRIX END
  ///

  switch (sh->get_algorithm_type()) {
    default:
    case SPGEMM_KK_TRIANGLE_LL: {
      KokkosSparse::Impl::KokkosSPGEMM<KernelHandle, row_lno_persistent_work_view_t, nnz_lno_persistent_work_view_t,
                                       nnz_lno_persistent_work_view_t, row_lno_persistent_work_view_t,
                                       nnz_lno_persistent_work_view_t, nnz_lno_persistent_work_view_t>
          kspgemm(handle, m, m, m, lower_triangular_matrix_rowmap, lower_triangular_matrix_entries, false,
                  lower_triangular_matrix_rowmap, lower_triangular_matrix_entries, false);
      kspgemm.KokkosSPGEMM_generic_triangle(visit_struct);
    } break;
    case SPGEMM_KK_TRIANGLE_LU: {
      KokkosSparse::Impl::KokkosSPGEMM<KernelHandle, row_lno_persistent_work_view_t, nnz_lno_persistent_work_view_t,
                                       nnz_lno_persistent_work_view_t, row_lno_persistent_work_view_t,
                                       nnz_lno_persistent_work_view_t, nnz_lno_persistent_work_view_t>
          kspgemm(handle, m, m, m, lower_triangular_matrix_rowmap, lower_triangular_matrix_entries, false,
                  upper_triangular_matrix_rowmap, upper_triangular_matrix_entries, false);
      kspgemm.KokkosSPGEMM_generic_triangle(visit_struct);
    } break;
    case SPGEMM_KK_TRIANGLE_AI: {
      if (create_lower_triangular) {
        KokkosSparse::Impl::KokkosSPGEMM<KernelHandle, row_lno_persistent_work_view_t, nnz_lno_persistent_work_view_t,
                                         nnz_lno_persistent_work_view_t, row_lno_persistent_work_view_t,
                                         nnz_lno_persistent_work_view_t, nnz_lno_persistent_work_view_t>
            kspgemm(handle, m, m, incidence_entries.extent(0) / 2, lower_triangular_matrix_rowmap,
                    lower_triangular_matrix_entries,
                    false,  // transpose ignore.
                    incidence_rowmap, incidence_entries, false);
        kspgemm.KokkosSPGEMM_generic_triangle(visit_struct);
      } else {
        KokkosSparse::Impl::KokkosSPGEMM<KernelHandle, alno_row_view_t_, alno_nnz_view_t_,
                                         nnz_lno_persistent_work_view_t, row_lno_persistent_work_view_t,
                                         nnz_lno_persistent_work_view_t, nnz_lno_persistent_work_view_t>
            kspgemm(handle, m, m, incidence_entries.extent(0) / 2, row_mapA, entriesA,
                    false,  // transpose ignore.
                    incidence_rowmap, incidence_entries, false);
        kspgemm.KokkosSPGEMM_generic_triangle(visit_struct);
      }
    }

    break;
    case SPGEMM_KK_TRIANGLE_IA_UNION:
    case SPGEMM_KK_TRIANGLE_IA: {
      if (create_lower_triangular) {
        KokkosSparse::Impl::KokkosSPGEMM<KernelHandle, row_lno_persistent_work_view_t, nnz_lno_persistent_work_view_t,
                                         nnz_lno_persistent_work_view_t, row_lno_persistent_work_view_t,
                                         nnz_lno_persistent_work_view_t, nnz_lno_persistent_work_view_t>
            kspgemm(handle, incidence_transpose_rowmap.extent(0) - 1, m, m, incidence_transpose_rowmap,
                    incidence_transpose_entries,
                    false,  // transpose ignore.
                    lower_triangular_matrix_rowmap, lower_triangular_matrix_entries, false);
        kspgemm.KokkosSPGEMM_generic_triangle(visit_struct);
      } else {
        KokkosSparse::Impl::KokkosSPGEMM<KernelHandle, row_lno_persistent_work_view_t, nnz_lno_persistent_work_view_t,
                                         nnz_lno_persistent_work_view_t, alno_row_view_t_, alno_nnz_view_t_,
                                         nnz_lno_persistent_work_view_t>
            kspgemm(handle, incidence_transpose_rowmap.extent(0) - 1, m, m, incidence_transpose_rowmap,
                    incidence_transpose_entries,
                    false,  // transpose ignore.
                    row_mapA, entriesA, false);
        kspgemm.KokkosSPGEMM_generic_triangle(visit_struct);
      }
    } break;
  }
}

}  // namespace Experimental
}  // namespace KokkosGraph
#endif
