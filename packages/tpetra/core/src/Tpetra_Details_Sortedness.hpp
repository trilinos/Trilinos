// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRAEXT_DETAILS_SORTEDNESS_HPP
#define TPETRAEXT_DETAILS_SORTEDNESS_HPP

namespace Tpetra::Details {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class MatrixTraits {
 public:
  static constexpr bool spgemmNeedsSortedInputs() {
    return true;
  }
};

#if defined(HAVE_TPETRA_SERIAL)
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal>
class MatrixTraits<Scalar, LocalOrdinal, GlobalOrdinal, KokkosCompat::KokkosSerialWrapperNode> {
 public:
  static constexpr bool spgemmNeedsSortedInputs() {
    return true;
  }
};
#endif

#if defined(HAVE_TPETRA_CUDA)
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal>
class MatrixTraits<Scalar, LocalOrdinal, GlobalOrdinal, KokkosCompat::KokkosCudaWrapperNode> {
 public:
  static constexpr bool spgemmNeedsSortedInputs() {
    // if Kokkos Kernels is going to use the cuSparse TPL for SpGEMM, this matrix
    // needs to be sorted
    // Try to mirror the Kokkos Kernels internal SpGEMM TPL use logic
    // inspired by https://github.com/trilinos/Trilinos/pull/11709
    // and https://github.com/kokkos/kokkos-kernels/pull/2008
#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOSKERNELS_ENABLE_TPL_CUSPARSE) && ((CUDA_VERSION < 11000) || (CUDA_VERSION >= 11040))
    return true;
#else
    return false;
#endif
  }
};
#endif

}  // namespace Tpetra::Details

#endif
