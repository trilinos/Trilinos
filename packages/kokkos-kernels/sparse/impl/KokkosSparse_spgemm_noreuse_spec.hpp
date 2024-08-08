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
#ifndef KOKKOSSPARSE_IMPL_SPGEMM_NOREUSE_SPEC_HPP_
#define KOKKOSSPARSE_IMPL_SPGEMM_NOREUSE_SPEC_HPP_

#include <KokkosKernels_config.h>

#include <Kokkos_Core.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY
#include "KokkosKernels_Handle.hpp"
#include "KokkosSparse_spgemm_symbolic.hpp"
#include "KokkosSparse_spgemm_numeric.hpp"
#endif

namespace KokkosSparse {
namespace Impl {
// Specialization struct which defines whether a specialization exists
template <class CMatrix, class AMatrix, class BMatrix>
struct spgemm_noreuse_eti_spec_avail {
  enum : bool { value = false };
};

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPGEMM_NOREUSE_ETI_SPEC_AVAIL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, EXEC_SPACE_TYPE,           \
                                                   MEM_SPACE_TYPE)                                                    \
  template <>                                                                                                         \
  struct spgemm_noreuse_eti_spec_avail<                                                                               \
      KokkosSparse::CrsMatrix<SCALAR_TYPE, ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, void,       \
                              OFFSET_TYPE>,                                                                           \
      KokkosSparse::CrsMatrix<const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,                            \
      KokkosSparse::CrsMatrix<const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>> {                          \
    enum : bool { value = true };                                                                                     \
  };

// Include the actual specialization declarations
#include <KokkosSparse_spgemm_noreuse_tpl_spec_avail.hpp>
#include <generated_specializations_hpp/KokkosSparse_spgemm_noreuse_eti_spec_avail.hpp>

namespace KokkosSparse {
namespace Impl {

// Unification layer
/// \brief Implementation of KokkosSparse::spgemm (sparse matrix - dense
///   vector multiply) for multiple vectors at a time (multivectors)
///   and possibly multiple coefficients at a time.

template <class CMatrix, class AMatrix, class BMatrix,
          bool tpl_spec_avail = spgemm_noreuse_tpl_spec_avail<CMatrix, AMatrix, BMatrix>::value,
          bool eti_spec_avail = spgemm_noreuse_eti_spec_avail<CMatrix, AMatrix, BMatrix>::value>
struct SPGEMM_NOREUSE {
  static CMatrix spgemm_noreuse(const AMatrix& A, bool transA, const BMatrix& B, bool transB);
};

#if !defined(KOKKOSKERNELS_ETI_ONLY) || KOKKOSKERNELS_IMPL_COMPILE_LIBRARY

// Unification layer
template <class CMatrix, class AMatrix, class BMatrix>
struct SPGEMM_NOREUSE<CMatrix, AMatrix, BMatrix, false, KOKKOSKERNELS_IMPL_COMPILE_LIBRARY> {
  static CMatrix spgemm_noreuse(const AMatrix& A, bool transA, const BMatrix& B, bool transB) {
    using device_t    = typename CMatrix::device_type;
    using scalar_t    = typename CMatrix::value_type;
    using ordinal_t   = typename CMatrix::ordinal_type;
    using size_type   = typename CMatrix::size_type;
    using c_rowmap_t  = typename CMatrix::row_map_type::non_const_type;
    using c_entries_t = typename CMatrix::index_type::non_const_type;
    using c_values_t  = typename CMatrix::values_type::non_const_type;
    KokkosKernels::Experimental::KokkosKernelsHandle<size_type, ordinal_t, scalar_t, typename device_t::execution_space,
                                                     typename device_t::memory_space, typename device_t::memory_space>
        kh;
    kh.create_spgemm_handle();
    // A is m*n, B is n*k, C is m*k
    ordinal_t m = A.numRows();
    ordinal_t n = B.numRows();
    ordinal_t k = B.numCols();
    c_rowmap_t row_mapC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "C rowmap"), m + 1);
    KokkosSparse::Experimental::spgemm_symbolic(&kh, m, n, k, A.graph.row_map, A.graph.entries, transA, B.graph.row_map,
                                                B.graph.entries, transB, row_mapC);
    size_type c_nnz = kh.get_spgemm_handle()->get_c_nnz();
    c_entries_t entriesC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "C entries"), c_nnz);
    c_values_t valuesC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "C values"), c_nnz);
    KokkosSparse::Experimental::spgemm_numeric(&kh, m, n, k, A.graph.row_map, A.graph.entries, A.values, transA,
                                               B.graph.row_map, B.graph.entries, B.values, transB, row_mapC, entriesC,
                                               valuesC);
    kh.destroy_spgemm_handle();
    return CMatrix("C", m, k, c_nnz, valuesC, row_mapC, entriesC);
  }
};

#endif

}  // namespace Impl
}  // namespace KokkosSparse

#define KOKKOSSPARSE_SPGEMM_NOREUSE_ETI_SPEC_DECL(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, EXEC_SPACE_TYPE,            \
                                                  MEM_SPACE_TYPE)                                                     \
  extern template struct SPGEMM_NOREUSE<                                                                              \
      KokkosSparse::CrsMatrix<SCALAR_TYPE, ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, void,       \
                              OFFSET_TYPE>,                                                                           \
      KokkosSparse::CrsMatrix<const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,                            \
      KokkosSparse::CrsMatrix<const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,                            \
      false, true>;

#define KOKKOSSPARSE_SPGEMM_NOREUSE_ETI_SPEC_INST(SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, EXEC_SPACE_TYPE,            \
                                                  MEM_SPACE_TYPE)                                                     \
  template struct SPGEMM_NOREUSE<                                                                                     \
      KokkosSparse::CrsMatrix<SCALAR_TYPE, ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, void,       \
                              OFFSET_TYPE>,                                                                           \
      KokkosSparse::CrsMatrix<const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,                            \
      KokkosSparse::CrsMatrix<const SCALAR_TYPE, const ORDINAL_TYPE, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
                              Kokkos::MemoryTraits<Kokkos::Unmanaged>, const OFFSET_TYPE>,                            \
      false, true>;

#include <KokkosSparse_spgemm_noreuse_tpl_spec_decl.hpp>

#endif  // KOKKOSSPARSE_IMPL_SPGEMM_NOREUSE_SPEC_HPP_
