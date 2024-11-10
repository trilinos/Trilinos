// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_ORDINALTRAITS_HPP
#define TPETRA_DETAILS_ORDINALTRAITS_HPP

/// \file Tpetra_Details_OrdinalTraits.hpp
/// \brief Import KokkosSparse::OrdinalTraits, a traits class for
///   "invalid" (flag) values of integer types, into the
///   Tpetra::Details namespace.

#include "Tpetra_ConfigDefs.hpp"
#include "KokkosSparse_OrdinalTraits.hpp"

namespace Tpetra {
namespace Details {

// KokkosSparse::OrdinalTraits is a traits class for "invalid" (flag)
// values of integer types.  Tpetra uses those flag values with both
// local and global indices.
using ::KokkosSparse::OrdinalTraits;

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_ORDINALTRAITS_HPP
