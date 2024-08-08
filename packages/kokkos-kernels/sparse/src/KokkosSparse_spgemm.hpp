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
#ifndef _KOKKOS_SPGEMM_HPP
#define _KOKKOS_SPGEMM_HPP

#include "KokkosSparse_spgemm_numeric.hpp"
#include "KokkosSparse_spgemm_symbolic.hpp"
#include "KokkosSparse_spgemm_jacobi.hpp"
#include "KokkosSparse_spgemm_noreuse_spec.hpp"

namespace KokkosSparse {

///
/// @brief
///
/// @tparam KernelHandle
/// @tparam AMatrix
/// @tparam BMatrix
/// @tparam CMatrix
/// @param kh
/// @param A
/// @param Amode
/// @param B
/// @param Bmode
/// @param C
////
template <class KernelHandle, class AMatrix, class BMatrix, class CMatrix>
void spgemm_symbolic(KernelHandle& kh, const AMatrix& A, const bool Amode, const BMatrix& B, const bool Bmode,
                     CMatrix& C) {
  using row_map_type = typename CMatrix::row_map_type::non_const_type;
  using entries_type = typename CMatrix::index_type::non_const_type;
  using values_type  = typename CMatrix::values_type::non_const_type;

  row_map_type row_mapC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "non_const_lnow_row"), A.numRows() + 1);
  entries_type entriesC;
  values_type valuesC;

  KokkosSparse::Experimental::spgemm_symbolic(&kh, A.numRows(), B.numRows(), B.numCols(), A.graph.row_map,
                                              A.graph.entries, Amode, B.graph.row_map, B.graph.entries, Bmode,
                                              row_mapC);

  const size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
  if (c_nnz_size) {
    entriesC = entries_type(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC"), c_nnz_size);
    valuesC  = values_type(Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesC"), c_nnz_size);
  }

  C = CMatrix("C=AB", A.numRows(), B.numCols(), c_nnz_size, valuesC, row_mapC, entriesC);
}

///
/// @brief Symbolic phase for block SpGEMM (BSR matrices)
///
/// @tparam KernelHandle
/// @tparam AMatrixType
/// @tparam BMatrixType
/// @tparam CMatrixType
/// @param kh
/// @param A
/// @param transposeA
/// @param B
/// @param transposeB
/// @param C
///
template <class KernelHandle, class AMatrixType, class BMatrixType, class CMatrixType>
void block_spgemm_symbolic(KernelHandle& kh, const AMatrixType& A, const bool transposeA, const BMatrixType& B,
                           const bool transposeB, CMatrixType& C) {
  using row_map_type = typename CMatrixType::row_map_type::non_const_type;
  using entries_type = typename CMatrixType::index_type::non_const_type;
  using values_type  = typename CMatrixType::values_type::non_const_type;

  auto blockDim = A.blockDim();
  if (blockDim != B.blockDim()) {
    throw std::invalid_argument("Block SpGEMM must be called for matrices with the same block size");
  }

  row_map_type row_mapC(Kokkos::view_alloc(Kokkos::WithoutInitializing, "non_const_lnow_row"), A.numRows() + 1);

  KokkosSparse::Experimental::spgemm_symbolic(&kh, A.numRows(), B.numRows(), B.numCols(), A.graph.row_map,
                                              A.graph.entries, transposeA, B.graph.row_map, B.graph.entries, transposeB,
                                              row_mapC);

  entries_type entriesC;
  values_type valuesC;
  const size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
  if (c_nnz_size) {
    entriesC = entries_type(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC"), c_nnz_size);
    valuesC = values_type(Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesC"), c_nnz_size * blockDim * blockDim);
  }

  C = CMatrixType("C=AB", A.numRows(), B.numCols(), c_nnz_size, valuesC, row_mapC, entriesC, blockDim);
}

///
/// @brief
///
/// @tparam KernelHandle
/// @tparam AMatrix
/// @tparam BMatrix
/// @tparam CMatrix
/// @param kh
/// @param A
/// @param Amode
/// @param B
/// @param Bmode
/// @param C
///
template <class KernelHandle, class AMatrix, class BMatrix, class CMatrix>
void spgemm_numeric(KernelHandle& kh, const AMatrix& A, const bool Amode, const BMatrix& B, const bool Bmode,
                    CMatrix& C) {
  // using row_map_type = typename CMatrix::index_type::non_const_type;
  // using entries_type = typename CMatrix::row_map_type::non_const_type;
  // using values_type  = typename CMatrix::values_type::non_const_type;

  KokkosSparse::Experimental::spgemm_numeric(&kh, A.numRows(), B.numRows(), B.numCols(), A.graph.row_map,
                                             A.graph.entries, A.values, Amode, B.graph.row_map, B.graph.entries,
                                             B.values, Bmode, C.graph.row_map, C.graph.entries, C.values);
}

///
/// @brief
///
/// @tparam KernelHandle
/// @tparam AMatrix
/// @tparam BMatrix
/// @tparam CMatrix
/// @param kh
/// @param A
/// @param Amode
/// @param B
/// @param Bmode
/// @param C
///
template <class KernelHandle, class AMatrix, class BMatrix, class CMatrix>
void block_spgemm_numeric(KernelHandle& kh, const AMatrix& A, const bool Amode, const BMatrix& B, const bool Bmode,
                          CMatrix& C) {
  auto blockDim = A.blockDim();
  if (blockDim != B.blockDim() || blockDim != C.blockDim()) {
    throw std::invalid_argument("Block SpGEMM must be called for matrices with the same block size");
  }

  KokkosSparse::Experimental::spgemm_numeric(&kh, A.numRows(), B.numRows(), B.numCols(), A.graph.row_map,
                                             A.graph.entries, A.values, Amode, B.graph.row_map, B.graph.entries,
                                             B.values, Bmode, C.graph.row_map, C.graph.entries, C.values, blockDim);
}

///
/// @brief
///
/// @tparam CMatrix
/// @tparam AMatrix
/// @tparam BMatrix
/// @param A
/// @param Amode
/// @param B
/// @param Bmode
/// @return CMatrix
///
template <class CMatrix, class AMatrix, class BMatrix>
CMatrix spgemm(const AMatrix& A, const bool Amode, const BMatrix& B, const bool Bmode) {
  // Canonicalize the matrix types:
  //  - Make A,B have const values and entries.
  //  - Make all views in A,B unmanaged, but otherwise default memory traits
  //  - C must have managed memory since its views are allocated in this
  //  function
  using AMatrix_Internal =
      KokkosSparse::CrsMatrix<typename AMatrix::const_value_type, typename AMatrix::const_ordinal_type,
                              typename AMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                              typename AMatrix::const_size_type>;
  using BMatrix_Internal =
      KokkosSparse::CrsMatrix<typename BMatrix::const_value_type, typename BMatrix::const_ordinal_type,
                              typename BMatrix::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>,
                              typename BMatrix::const_size_type>;
  using CMatrix_Internal =
      KokkosSparse::CrsMatrix<typename CMatrix::non_const_value_type, typename CMatrix::non_const_ordinal_type,
                              typename CMatrix::device_type, void, typename CMatrix::non_const_size_type>;
  // Check now that A, B dimensions are compatible to multiply
  auto opACols = Amode ? A.numRows() : A.numCols();
  auto opBRows = Bmode ? B.numCols() : B.numRows();
  if (Amode || Bmode) throw std::invalid_argument("KokkosSparse::spgemm: transposing A and/or B is not yet supported");
  if (opACols != opBRows)
    throw std::invalid_argument(
        "KokkosSparse::spgemm: op(A) and op(B) have incompatible dimensions "
        "for multiplication");
  // Make sure C has managed memory. If its memory traits are void (default),
  // then that also means it's managed.
  if constexpr (!std::is_same<typename CMatrix::memory_traits, void>::value) {
    if (CMatrix::memory_traits::is_unmanaged)
      throw std::invalid_argument(
          "KokkosSparse::spgemm: C must not have the Unmanaged memory trait, "
          "because spgemm needs to allocate its Views");
  }
  AMatrix_Internal A_internal(A);
  BMatrix_Internal B_internal(B);
  // Intercept empty C case here so that TPL wrappers don't have to deal with it
  if (!A.numRows() || !A.numCols() || !B.numCols() || !A.nnz() || !B.nnz()) {
    auto Crows = Amode ? A.numCols() : A.numRows();
    auto Ccols = Bmode ? B.numRows() : B.numCols();
    typename CMatrix::row_map_type::non_const_type row_mapC("C rowmap", Crows + 1);
    typename CMatrix::index_type entriesC;
    typename CMatrix::values_type valuesC;
    return CMatrix("C", Crows, Ccols, 0, valuesC, row_mapC, entriesC);
  }
  return CMatrix(
      KokkosSparse::Impl::SPGEMM_NOREUSE<CMatrix_Internal, AMatrix_Internal, BMatrix_Internal>::spgemm_noreuse(
          A_internal, Amode, B_internal, Bmode));
}

}  // namespace KokkosSparse

#endif
