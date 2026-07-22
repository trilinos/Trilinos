// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_BLOCKHELPER_ETI_HPP
#define IFPACK2_BLOCKHELPER_ETI_HPP

#include <complex>

namespace Ifpack2::BlockTriDiContainerDetails {

///
/// impl tag to distinguish built-in types and sacado types.
///

// The BlockTriDiContainer implementation supports only the four
// basic scalar types (those for which ImplTag<Scalar>::type == ImplSimdTag).
//
// BTDC and related classes have specializations with ImplSimdTag that actually
// have the implementation, and specializations with ImplNotAvailTag where all
// the methods do nothing. This way, ETI can be performed even when
// ImplTag<Scalar>::type == ImplNotAvailTag to avoid needing special ETI logic in Stokhos.

struct ImplNotAvailTag {};
struct ImplSimdTag {};
struct ImplSacadoTag {};

template <typename T>
struct ImplTag {
  typedef ImplNotAvailTag type;
};
template <>
struct ImplTag<float> {
  typedef ImplSimdTag type;
};
template <>
struct ImplTag<double> {
  typedef ImplSimdTag type;
};
template <>
struct ImplTag<std::complex<float> > {
  typedef ImplSimdTag type;
};
template <>
struct ImplTag<std::complex<double> > {
  typedef ImplSimdTag type;
};

/// forward declaration
template <typename MatrixType>
struct ImplObject;
}  // namespace Ifpack2::BlockTriDiContainerDetails

#endif
