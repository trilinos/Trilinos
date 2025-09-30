// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_LOCALMAP_FWD_HPP
#define TPETRA_DETAILS_LOCALMAP_FWD_HPP

#include "Tpetra_Details_DefaultTypes.hpp"

/// \file Tpetra_Details_LocalMap_fwd.hpp
/// \brief Forward declaration of Tpetra::Details::LocalMap

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Tpetra {
namespace Details {
template <class LocalOrdinal  = ::Tpetra::Details::DefaultTypes::local_ordinal_type,
          class GlobalOrdinal = ::Tpetra::Details::DefaultTypes::global_ordinal_type,
          class DeviceType    = ::Tpetra::Details::DefaultTypes::node_type::device_type>
class LocalMap;
}  // namespace Details
}  // namespace Tpetra
#endif  // DOXYGEN_SHOULD_SKIP_THIS

#endif  // TPETRA_DETAILS_LOCALMAP_FWD_HPP
