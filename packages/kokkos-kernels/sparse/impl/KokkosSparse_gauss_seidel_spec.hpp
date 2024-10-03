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
#ifndef KOKKOSSPARSE_IMPL_GAUSS_SEIDEL_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_GAUSS_SEIDEL_SPEC_HPP_

#include <KokkosKernels_config.h>

#include <Kokkos_Core.hpp>
// #include <Kokkos_ArithTraits.hpp>
#include "KokkosKernels_Handle.hpp"
// Include the actual functors
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include "KokkosSparse_gauss_seidel_impl.hpp"
#include "KokkosSparse_cluster_gauss_seidel_impl.hpp"
#include "KokkosSparse_twostage_gauss_seidel_impl.hpp"
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t>
struct gauss_seidel_symbolic_eti_spec_avail {
  enum : bool { value = false };
};
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t>
struct gauss_seidel_numeric_eti_spec_avail {
  enum : bool { value = false };
};
template <class KernelHandle, class a_size_view_t_, class a_lno_view_t, class a_scalar_view_t, class x_scalar_view_t,
          class y_scalar_view_t>
struct gauss_seidel_apply_eti_spec_avail {
  enum : bool { value = false };
};
}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_GAUSS_SEIDEL_SYMBOLIC_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,   \
                                                          EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                       \
  template <>                                                                                                    \
  struct gauss_seidel_symbolic_eti_spec_avail<                                                                   \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                                 \
    enum : bool { value = true };                                                                                \
  };

#define KOKKOSSPARSE_GAUSS_SEIDEL_NUMERIC_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,    \
                                                         EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                        \
  template <>                                                                                                    \
  struct gauss_seidel_numeric_eti_spec_avail<                                                                    \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                                 \
    enum : bool { value = true };                                                                                \
  };

#define KOKKOSSPARSE_GAUSS_SEIDEL_APPLY_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,      \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                          \
  template <>                                                                                                    \
  struct gauss_seidel_apply_eti_spec_avail<                                                                      \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> > > {                                                 \
    enum : bool { value = true };                                                                                \
  };

// Include the actual specialization declarations
#include <KokkosSparse_gauss_seidel_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_gauss_seidel_symbolic_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_gauss_seidel_numeric_eti_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_gauss_seidel_apply_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

template <class ExecSpaceIn, class KernelHandle, class a_size_view_t_, class a_lno_view_t,
          bool tpl_spec_avail = gauss_seidel_symbolic_tpl_spec_avail<KernelHandle, a_size_view_t_, a_lno_view_t>::value,
          bool eti_spec_avail = gauss_seidel_symbolic_eti_spec_avail<KernelHandle, a_size_view_t_, a_lno_view_t>::value>
struct GAUSS_SEIDEL_SYMBOLIC {
  static void gauss_seidel_symbolic(const ExecSpaceIn &exec_space_in, KernelHandle *handle,
                                    typename KernelHandle::const_nnz_lno_t num_rows,
                                    typename KernelHandle::const_nnz_lno_t num_cols, a_size_view_t_ row_map,
                                    a_lno_view_t entries, bool is_graph_symmetric);
};

template <class ExecSpaceIn, class KernelHandle, KokkosSparse::SparseMatrixFormat format, class a_size_view_t_,
          class a_lno_view_t, class a_scalar_view_t,
          bool tpl_spec_avail =
              gauss_seidel_numeric_tpl_spec_avail<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t>::value,
          bool eti_spec_avail =
              gauss_seidel_numeric_eti_spec_avail<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t>::value>
struct GAUSS_SEIDEL_NUMERIC {
  static void gauss_seidel_numeric(const ExecSpaceIn &exec_space_in, KernelHandle *handle,
                                   typename KernelHandle::const_nnz_lno_t num_rows,
                                   typename KernelHandle::const_nnz_lno_t num_cols, a_size_view_t_ row_map,
                                   a_lno_view_t entries, a_scalar_view_t values, bool is_graph_symmetric);

  static void gauss_seidel_numeric(const ExecSpaceIn &exec_space_in, KernelHandle *handle,
                                   typename KernelHandle::const_nnz_lno_t num_rows,
                                   typename KernelHandle::const_nnz_lno_t num_cols, a_size_view_t_ row_map,
                                   a_lno_view_t entries, a_scalar_view_t values, a_scalar_view_t given_inverse_diagonal,
                                   bool is_graph_symmetric);
};

template <class ExecSpaceIn, class KernelHandle, KokkosSparse::SparseMatrixFormat format, class a_size_view_t_,
          class a_lno_view_t, class a_scalar_view_t, class x_scalar_view_t, class y_scalar_view_t,
          bool tpl_spec_avail = gauss_seidel_apply_tpl_spec_avail<
              KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t, x_scalar_view_t, y_scalar_view_t>::value,
          bool eti_spec_avail = gauss_seidel_apply_eti_spec_avail<
              KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t, x_scalar_view_t, y_scalar_view_t>::value>
struct GAUSS_SEIDEL_APPLY {
  static void gauss_seidel_apply(const ExecSpaceIn &exec_space_in, KernelHandle *handle,
                                 typename KernelHandle::const_nnz_lno_t num_rows,
                                 typename KernelHandle::const_nnz_lno_t num_cols, a_size_view_t_ row_map,
                                 a_lno_view_t entries, a_scalar_view_t values, x_scalar_view_t x_lhs_output_vec,
                                 y_scalar_view_t y_rhs_input_vec, bool init_zero_x_vector, bool update_y_vector,
                                 typename KernelHandle::nnz_scalar_t omega, int numIter, bool apply_forward,
                                 bool apply_backward);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

template <class ExecSpaceIn, class KernelHandle, class a_size_view_t_, class a_lno_view_t_>
struct GAUSS_SEIDEL_SYMBOLIC<ExecSpaceIn, KernelHandle, a_size_view_t_, a_lno_view_t_, false,
                             KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void gauss_seidel_symbolic(const ExecSpaceIn &exec_space_in, KernelHandle *handle,
                                    typename KernelHandle::const_nnz_lno_t num_rows,
                                    typename KernelHandle::const_nnz_lno_t num_cols, a_size_view_t_ row_map,
                                    a_lno_view_t_ entries, bool is_graph_symmetric) {
    Kokkos::Profiling::pushRegion("KokkosSparse::Impl::gauss_seidel_symbolic");
    auto gsHandle = handle->get_gs_handle();
    gsHandle->set_execution_space(exec_space_in);
    if (gsHandle->get_algorithm_type() == GS_CLUSTER) {
      using SGS = typename Impl::ClusterGaussSeidel<KernelHandle, a_size_view_t_, a_lno_view_t_,
                                                    typename KernelHandle::in_scalar_nnz_view_t>;
      SGS sgs(handle, num_rows, num_cols, row_map, entries, is_graph_symmetric);
      sgs.initialize_symbolic();
    } else if (gsHandle->get_algorithm_type() == GS_TWOSTAGE) {
      using SGS = typename Impl::TwostageGaussSeidel<KernelHandle, a_size_view_t_, a_lno_view_t_,
                                                     typename KernelHandle::in_scalar_nnz_view_t>;
      SGS sgs(handle, num_rows, num_cols, row_map, entries);
      sgs.initialize_symbolic();
    } else {
      using SGS = typename Impl::PointGaussSeidel<KernelHandle, a_size_view_t_, a_lno_view_t_,
                                                  typename KernelHandle::in_scalar_nnz_view_t>;
      SGS sgs(handle, num_rows, num_cols, row_map, entries, is_graph_symmetric);
      sgs.initialize_symbolic();
    }
    Kokkos::Profiling::popRegion();
  }
};

template <class ExecSpaceIn, class KernelHandle, KokkosSparse::SparseMatrixFormat format, class a_size_view_t_,
          class a_lno_view_t, class a_scalar_view_t>
struct GAUSS_SEIDEL_NUMERIC<ExecSpaceIn, KernelHandle, format, a_size_view_t_, a_lno_view_t, a_scalar_view_t, false,
                            KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void gauss_seidel_numeric(const ExecSpaceIn &exec_space_in, KernelHandle *handle,
                                   typename KernelHandle::const_nnz_lno_t num_rows,
                                   typename KernelHandle::const_nnz_lno_t num_cols, a_size_view_t_ row_map,
                                   a_lno_view_t entries, a_scalar_view_t values, bool is_graph_symmetric) {
    Kokkos::Profiling::pushRegion("KokkosSparse::Impl::gauss_seidel_numeric");
    auto gsHandle = handle->get_gs_handle();
    gsHandle->set_execution_space(exec_space_in);
    if (gsHandle->get_algorithm_type() == GS_CLUSTER) {
      using SGS = typename Impl::ClusterGaussSeidel<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t>;
      SGS sgs(handle, num_rows, num_cols, row_map, entries, values, is_graph_symmetric);
      sgs.initialize_numeric();
    } else if (gsHandle->get_algorithm_type() == GS_TWOSTAGE) {
      using SGS = typename Impl::TwostageGaussSeidel<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t>;
      SGS sgs(handle, num_rows, num_cols, row_map, entries, values);
      sgs.initialize_numeric();
    } else {
      using SGS = typename Impl::PointGaussSeidel<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t, format>;
      SGS sgs(handle, num_rows, num_cols, row_map, entries, values, is_graph_symmetric);
      sgs.initialize_numeric();
    }
    Kokkos::Profiling::popRegion();
  }

  static void gauss_seidel_numeric(const ExecSpaceIn &exec_space_in, KernelHandle *handle,
                                   typename KernelHandle::const_nnz_lno_t num_rows,
                                   typename KernelHandle::const_nnz_lno_t num_cols, a_size_view_t_ row_map,
                                   a_lno_view_t entries, a_scalar_view_t values, a_scalar_view_t given_inverse_diagonal,
                                   bool is_graph_symmetric) {
    Kokkos::Profiling::pushRegion("KokkosSparse::Impl::gauss_seidel_numeric");
    auto gsHandle = handle->get_gs_handle();
    gsHandle->set_execution_space(exec_space_in);
    if (gsHandle->get_algorithm_type() == GS_CLUSTER) {
      using SGS = typename Impl::ClusterGaussSeidel<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t>;
      SGS sgs(handle, num_rows, num_cols, row_map, entries, values, given_inverse_diagonal, is_graph_symmetric);
      sgs.initialize_numeric();
    } else if (gsHandle->get_algorithm_type() == GS_TWOSTAGE) {
      using SGS = typename Impl::TwostageGaussSeidel<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t>;
      SGS sgs(handle, num_rows, num_cols, row_map, entries, values, given_inverse_diagonal);
      sgs.initialize_numeric();
    } else {
      using SGS = typename Impl::PointGaussSeidel<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t>;
      SGS sgs(handle, num_rows, num_cols, row_map, entries, values, given_inverse_diagonal, is_graph_symmetric);
      sgs.initialize_numeric();
    }
    Kokkos::Profiling::popRegion();
  }
};

template <class ExecSpaceIn, class KernelHandle, KokkosSparse::SparseMatrixFormat format, class a_size_view_t_,
          class a_lno_view_t, class a_scalar_view_t, class x_scalar_view_t, class y_scalar_view_t>
struct GAUSS_SEIDEL_APPLY<ExecSpaceIn, KernelHandle, format, a_size_view_t_, a_lno_view_t, a_scalar_view_t,
                          x_scalar_view_t, y_scalar_view_t, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static void gauss_seidel_apply(const ExecSpaceIn &exec_space_in, KernelHandle *handle,
                                 typename KernelHandle::const_nnz_lno_t num_rows,
                                 typename KernelHandle::const_nnz_lno_t num_cols, a_size_view_t_ row_map,
                                 a_lno_view_t entries, a_scalar_view_t values, x_scalar_view_t x_lhs_output_vec,
                                 y_scalar_view_t y_rhs_input_vec, bool init_zero_x_vector, bool update_y_vector,
                                 typename KernelHandle::nnz_scalar_t omega, int numIter, bool apply_forward,
                                 bool apply_backward) {
    Kokkos::Profiling::pushRegion("KokkosSparse::Impl::gauss_seidel_apply");
    auto gsHandle = handle->get_gs_handle();
    gsHandle->set_execution_space(exec_space_in);
    if (gsHandle->get_algorithm_type() == GS_CLUSTER) {
      using SGS = typename Impl::ClusterGaussSeidel<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t>;
      SGS sgs(handle, num_rows, num_cols, row_map, entries, values);
      sgs.apply(x_lhs_output_vec, y_rhs_input_vec, init_zero_x_vector, numIter, omega, apply_forward, apply_backward,
                update_y_vector);
    } else if (gsHandle->get_algorithm_type() == GS_TWOSTAGE) {
      using SGS = typename Impl::TwostageGaussSeidel<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t>;
      SGS sgs(handle, num_rows, num_cols, row_map, entries, values);
      sgs.apply(x_lhs_output_vec, y_rhs_input_vec, init_zero_x_vector, numIter, omega, apply_forward, apply_backward,
                update_y_vector);
    } else {
      using SGS = typename Impl::PointGaussSeidel<KernelHandle, a_size_view_t_, a_lno_view_t, a_scalar_view_t, format>;
      SGS sgs(handle, num_rows, num_cols, row_map, entries, values);
      sgs.apply(x_lhs_output_vec, y_rhs_input_vec, init_zero_x_vector, numIter, omega, apply_forward, apply_backward,
                update_y_vector);
    }
    Kokkos::Profiling::popRegion();
  }
};
#endif

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_GAUSS_SEIDEL_SYMBOLIC_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,    \
                                                         EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                        \
  extern template struct GAUSS_SEIDEL_SYMBOLIC<                                                                  \
      EXEC_SPACE_TYPE,                                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      false, true>;

#define KOKKOSSPARSE_GAUSS_SEIDEL_SYMBOLIC_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,    \
                                                         EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                        \
  template struct GAUSS_SEIDEL_SYMBOLIC<                                                                         \
      EXEC_SPACE_TYPE,                                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      false, true>;

#define KOKKOSSPARSE_GAUSS_SEIDEL_NUMERIC_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,     \
                                                        EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                         \
  extern template struct GAUSS_SEIDEL_NUMERIC<                                                                   \
      EXEC_SPACE_TYPE,                                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      KokkosSparse::SparseMatrixFormat::BSR,                                                                     \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      false, true>;                                                                                              \
  extern template struct GAUSS_SEIDEL_NUMERIC<                                                                   \
      EXEC_SPACE_TYPE,                                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      KokkosSparse::SparseMatrixFormat::CRS,                                                                     \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      false, true>;

#define KOKKOSSPARSE_GAUSS_SEIDEL_NUMERIC_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,     \
                                                        EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                         \
  template struct GAUSS_SEIDEL_NUMERIC<                                                                          \
      EXEC_SPACE_TYPE,                                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      KokkosSparse::SparseMatrixFormat::BSR,                                                                     \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      false, true>;                                                                                              \
  template struct GAUSS_SEIDEL_NUMERIC<                                                                          \
      EXEC_SPACE_TYPE,                                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      KokkosSparse::SparseMatrixFormat::CRS,                                                                     \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      false, true>;

#define KOKKOSSPARSE_GAUSS_SEIDEL_APPLY_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,       \
                                                      EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                           \
  extern template struct GAUSS_SEIDEL_APPLY<                                                                     \
      EXEC_SPACE_TYPE,                                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      KokkosSparse::SparseMatrixFormat::BSR,                                                                     \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      false, true>;                                                                                              \
  extern template struct GAUSS_SEIDEL_APPLY<                                                                     \
      EXEC_SPACE_TYPE,                                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      KokkosSparse::SparseMatrixFormat::CRS,                                                                     \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      false, true>;

#define KOKKOSSPARSE_GAUSS_SEIDEL_APPLY_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE,       \
                                                      EXEC_SPACE_TYPE, MEM_SPACE_TYPE)                           \
  template struct GAUSS_SEIDEL_APPLY<                                                                            \
      EXEC_SPACE_TYPE,                                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      KokkosSparse::SparseMatrixFormat::BSR,                                                                     \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      false, true>;                                                                                              \
  template struct GAUSS_SEIDEL_APPLY<                                                                            \
      EXEC_SPACE_TYPE,                                                                                           \
      KokkosKernels::Experimental::KokkosKernelsHandle<const OFFSET_TYPE, const ORDINAL_TYPE, const SCALAR_TYPE, \
                                                       EXEC_SPACE_TYPE, MEM_SPACE_TYPE, MEM_SPACE_TYPE>,         \
      KokkosSparse::SparseMatrixFormat::CRS,                                                                     \
      Kokkos::View<const OFFSET_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const ORDINAL_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE *, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,            \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,                 \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      Kokkos::View<const SCALAR_TYPE **, LAYOUT_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>,           \
                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >,                                                    \
      false, true>;

#include <KokkosSparse_gauss_seidel_tpl_spec_decl.hpp>

#endif  // KOKKOS_BLAS1_MV_IMPL_DOT_HPP_
