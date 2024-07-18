// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_LOCALOPERATOR_FWD_HPP
#define TPETRA_LOCALOPERATOR_FWD_HPP

#include "Tpetra_Details_DefaultTypes.hpp"

/// \file Tpetra_LocalOperator_fwd.hpp
/// \brief Forward declaration of Tpetra::LocalOperator

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Tpetra {
template<class Scalar = ::Tpetra::Details::DefaultTypes::scalar_type,
         class Device = ::Tpetra::Details::DefaultTypes::node_type::device_type>
class LocalOperator;
} // namespace Tpetra
#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif // TPETRA_LOCALOPERATOR_FWD_HPP
