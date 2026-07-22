// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/// \file KokkosSparse_getDiagCopy.hpp
/// \brief Get a copy of the diagonal entries of a KokkosSparse::CrsMatrix.

#ifndef KOKKOSSPARSE_GETDIAGCOPY_HPP_
#define KOKKOSSPARSE_GETDIAGCOPY_HPP_

#include "KokkosSparse_getDiagCopyWithOffsets_impl.hpp"
#include <type_traits>

namespace KokkosSparse {

template <class DiagType, class OffsetsType, class CrsMatrixType>
void getDiagCopy(const DiagType& D, const OffsetsType& offsets, const CrsMatrixType& A) {
  static_assert(Kokkos::is_view<DiagType>::value, "The DiagType template parameter must be a Kokkos::View.");
  static_assert(static_cast<int>(DiagType::rank) == 1, "The DiagType template parameter must be a 1-D Kokkos::View.");
  static_assert(std::is_same<typename DiagType::value_type, typename DiagType::non_const_value_type>::value,
                "The DiagType template parameter must be a nonconst Kokkos::View.");
  static_assert(Kokkos::is_view<OffsetsType>::value, "The OffsetsType template parameter must be a Kokkos::View.");
  static_assert(static_cast<int>(OffsetsType::rank) == 1,
                "The OffsetsType template parameter must be a 1-D Kokkos::View.");

  typedef typename CrsMatrixType::value_type scalar_type;
  typedef typename CrsMatrixType::ordinal_type ordinal_type;
  typedef typename CrsMatrixType::device_type device_type;
  typedef typename CrsMatrixType::size_type offset_type;

  // Standardize on unmanaged Views, in order to avoid proliferation
  // of instantiations of the implementation type.
  Kokkos::View<typename DiagType::non_const_value_type*, typename DiagType::array_layout,
               typename DiagType::device_type, Kokkos::MemoryUnmanaged>
      D_internal = D;
  Kokkos::View<typename OffsetsType::const_value_type*, typename OffsetsType::array_layout,
               typename OffsetsType::device_type, Kokkos::MemoryUnmanaged>
      offsets_internal = offsets;

  typedef Impl::CrsMatrixGetDiagCopyWithOffsets<scalar_type, ordinal_type, device_type, offset_type> impl_type;
  impl_type::getDiagCopy(D_internal, offsets_internal, A);
}

}  // namespace KokkosSparse

#endif  // KOKOS_SPARSE_GETDIAGCOPY_HPP_
