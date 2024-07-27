// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DISTOBJECT_FWD_HPP
#define TPETRA_DISTOBJECT_FWD_HPP

#include "Tpetra_Details_DefaultTypes.hpp"

/// \file Tpetra_DistObject_fwd.hpp
/// \brief Forward declaration of Tpetra::DistObject

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Tpetra {
template<class Packet,
         class LocalOrdinal = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
         class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
         class Node = ::Tpetra::Details::DefaultTypes::node_type>
class DistObject;
} // namespace Tpetra
#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif // TPETRA_DISTOBJECT_FWD_HPP
