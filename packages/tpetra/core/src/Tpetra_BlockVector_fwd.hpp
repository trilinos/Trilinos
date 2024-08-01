// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKVECTOR_FWD_HPP
#define TPETRA_BLOCKVECTOR_FWD_HPP

#include "Tpetra_Details_DefaultTypes.hpp"

/// \file Tpetra_BlockVector_fwd.hpp
/// \brief Forward declaration of Tpetra::BlockVector.

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Tpetra {
template<class SC = ::Tpetra::Details::DefaultTypes::scalar_type,
         class LO = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
         class GO = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
         class NT = ::Tpetra::Details::DefaultTypes::node_type>
class BlockVector;
} // namespace Tpetra
#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif // TPETRA_BLOCKVECTOR_FWD_HPP
