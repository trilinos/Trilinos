// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_MAP_DEF_HPP
#define XPETRA_MAP_DEF_HPP

#include "Xpetra_Map_decl.hpp"

namespace Xpetra {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Map<LocalOrdinal, GlobalOrdinal, Node>::
    ~Map() {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >
Map<LocalOrdinal, GlobalOrdinal, Node>::
    getMap() const {
  return rcpFromRef(*this);
}

}  // namespace Xpetra

#endif  // XPETRA_MAP_DEF_HPP
