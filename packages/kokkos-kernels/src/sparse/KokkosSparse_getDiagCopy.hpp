/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

/// \file Kokkos_Sparse_getDiagCopy.hpp
/// \brief Get a copy of the diagonal entries of a KokkosSparse::CrsMatrix.

#ifndef KOKKOS_SPARSE_GETDIAGCOPY_HPP_
#define KOKKOS_SPARSE_GETDIAGCOPY_HPP_

#include "KokkosSparse_getDiagCopyWithOffsets_impl.hpp"
#include <type_traits>

namespace KokkosSparse {

template<class DiagType,
         class OffsetsType,
         class CrsMatrixType>
void
getDiagCopy (const DiagType& D,
             const OffsetsType& offsets,
             const CrsMatrixType& A)
{
  static_assert (Kokkos::Impl::is_view<DiagType>::value,
                 "The DiagType template parameter must be a Kokkos::View.");
  static_assert (static_cast<int> (DiagType::rank) == 1,
                 "The DiagType template parameter must be a 1-D Kokkos::View.");
  static_assert (std::is_same<DiagType, typename DiagType::non_const_type>::value,
                 "The DiagType template parameter must be a nonconst Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<OffsetsType>::value,
                 "The OffsetsType template parameter must be a Kokkos::View.");
  static_assert (static_cast<int> (OffsetsType::rank) == 1,
                 "The OffsetsType template parameter must be a 1-D Kokkos::View.");

  typedef typename CrsMatrixType::value_type scalar_type;
  typedef typename CrsMatrixType::ordinal_type ordinal_type;
  typedef typename CrsMatrixType::device_type device_type;
  typedef typename CrsMatrixType::size_type offset_type;

  // Standardize on unmanaged Views, in order to avoid proliferation
  // of instantiations of the implementation type.
  Kokkos::View<typename DiagType::non_const_value_type*,
    typename DiagType::array_layout,
    typename DiagType::device_type,
    Kokkos::MemoryUnmanaged> D_internal = D;
  Kokkos::View<typename OffsetsType::const_value_type*,
    typename OffsetsType::array_layout,
    typename OffsetsType::device_type,
    Kokkos::MemoryUnmanaged> offsets_internal = offsets;

  typedef Impl::CrsMatrixGetDiagCopyWithOffsets<scalar_type,
    ordinal_type, device_type, offset_type> impl_type;
  impl_type::getDiagCopy (D_internal, offsets_internal, A);
}

} // namespace KokkosSparse

#endif // KOKKOS_SPARSE_GETDIAGCOPY_HPP_

