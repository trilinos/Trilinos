// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_LOCALCRSMATRIXOPERATOR_FWD_HPP
#define TPETRA_LOCALCRSMATRIXOPERATOR_FWD_HPP

#include "Tpetra_Details_DefaultTypes.hpp"

/// \file Tpetra_LocalOperator_fwd.hpp
/// \brief Forward declaration of Tpetra::LocalCrsMatrixOperator

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Tpetra {
template<class MultiVectorScalar = ::Tpetra::Details::DefaultTypes::scalar_type,
         class MatrixScalar = ::Tpetra::Details::DefaultTypes::scalar_type,
         class Device = ::Tpetra::Details::DefaultTypes::node_type::device_type>
class LocalCrsMatrixOperator;
} // namespace Tpetra
#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif // TPETRA_LOCALCRSMATRIXOPERATOR_FWD_HPP
