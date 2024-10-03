// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_PACKABLE_FWD_HPP
#define TPETRA_PACKABLE_FWD_HPP

#include "Tpetra_Details_DefaultTypes.hpp"

/// \file Tpetra_Packable_fwd.hpp
/// \brief Forward declaration of Tpetra::Packable

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Tpetra {
template<class Packet,
         class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type>
class Packable;
} // namespace Tpetra
#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif // TPETRA_PACKABLE_FWD_HPP
