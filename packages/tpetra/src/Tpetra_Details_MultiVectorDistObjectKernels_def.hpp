/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#ifndef TPETRA_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_DEF_HPP
#define TPETRA_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_DEF_HPP

#include "Tpetra_ConfigDefs.hpp"
#if TPETRA_USE_KOKKOS_DISTOBJECT

namespace Tpetra {
namespace Details {

////////////////////////////////////////////////////////////////////////////////

template <typename Scalar, typename LocalOrdinal, typename Device>
void
PackArraySingleColumnConstantStride<Scalar,LocalOrdinal,Device>::
pack() {
  Kokkos::parallel_for( exportLIDs.size(), *this );
}

////////////////////////////////////////////////////////////////////////////////

template <typename Scalar, typename LocalOrdinal, typename Device>
void
PackArraySingleColumnOffset<Scalar,LocalOrdinal,Device>::
pack() {
  Kokkos::parallel_for( exportLIDs.size(), *this );
}

////////////////////////////////////////////////////////////////////////////////

template <typename Scalar, typename LocalOrdinal, typename Device>
void
PackArrayMultiColumnConstantStride<Scalar,LocalOrdinal,Device>::
pack() {
  Kokkos::parallel_for( exportLIDs.size(), *this );
}

////////////////////////////////////////////////////////////////////////////////

template <typename Scalar, typename LocalOrdinal, typename Device>
void
PackArrayMultiColumnVariableStride<Scalar,LocalOrdinal,Device>::
pack() {
  Kokkos::parallel_for( exportLIDs.size(), *this );
}

////////////////////////////////////////////////////////////////////////////////

template <typename Scalar, typename LocalOrdinal, typename Op, typename Device>
void
UnpackArrayMultiColumnConstantStride<Scalar,LocalOrdinal,Op,Device>::
unpack() {
  Kokkos::parallel_for( importLIDs.size(), *this );
}

////////////////////////////////////////////////////////////////////////////////

template <typename Scalar, typename LocalOrdinal, typename Op, typename Device>
void
UnpackArrayMultiColumnVariableStride<Scalar,LocalOrdinal,Op,Device>::
unpack() {
  Kokkos::parallel_for( importLIDs.size(), *this );
}

////////////////////////////////////////////////////////////////////////////////

template <typename Scalar, typename LocalOrdinal, typename Device>
void
PermuteArrayMultiColumnConstantStride<Scalar,LocalOrdinal,Device>::
permute() {
  const size_type numPermuteLIDs =
    std::min (permuteToLIDs.size (), permuteFromLIDs.size ());
  Kokkos::parallel_for( numPermuteLIDs, *this );
}

////////////////////////////////////////////////////////////////////////////////

template <typename Scalar, typename LocalOrdinal, typename Device>
void
PermuteArrayMultiColumnVariableStride<Scalar,LocalOrdinal,Device>::
permute() {
  const size_type numPermuteLIDs =
    std::min (permuteToLIDs.size (), permuteFromLIDs.size ());
  Kokkos::parallel_for( numPermuteLIDs, *this );
}

#define PACK_KERNELS_INSTANT(SC,LO,DE) \
  template class PackArraySingleColumnConstantStride< SC, LO, DE >; \
  template class PackArraySingleColumnOffset< SC, LO, DE >; \
  template class PackArrayMultiColumnConstantStride< SC, LO, DE >; \
  template class PackArrayMultiColumnVariableStride< SC, LO, DE >;

#define UNPACK_KERNELS_INSTANT(SC,LO,OP,DE) \
  template class UnpackArrayMultiColumnConstantStride< SC, LO, OP, DE >; \
  template class UnpackArrayMultiColumnVariableStride< SC, LO, OP, DE >;

#define PERM_KERNELS_INSTANT(SC,LO,DE) \
  template class PermuteArrayMultiColumnConstantStride< SC, LO, DE >; \
  template class PermuteArrayMultiColumnVariableStride< SC, LO, DE >;

#define KERNELS_INSTANT(SC,LO,DE)               \
  PACK_KERNELS_INSTANT(SC,LO,DE)                \
  PERM_KERNELS_INSTANT(SC,LO,DE)                \
  UNPACK_KERNELS_INSTANT(SC,LO,InsertOp,DE)     \
  UNPACK_KERNELS_INSTANT(SC,LO,AddOp,DE)        \
  UNPACK_KERNELS_INSTANT(SC,LO,AbsMaxOp,DE)

} // Details namespace
} // Tpetra namespace

#endif // TPETRA_USE_KOKKOS_DISTOBJECT

#endif // TPETRA_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_DEF_HPP
