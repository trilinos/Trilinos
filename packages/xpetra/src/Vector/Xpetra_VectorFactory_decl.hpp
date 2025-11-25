// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_VECTORFACTORY_DECL_HPP
#define XPETRA_VECTORFACTORY_DECL_HPP

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Vector.hpp"

#include "Xpetra_TpetraVector_decl.hpp"

#include "Xpetra_BlockedMap_decl.hpp"
#include "Xpetra_BlockedVector_decl.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

template <class Scalar /* = Vector<>::scalar_type*/,
          class LocalOrdinal /* = typename Vector<Scalar>::local_ordinal_type*/,
          class GlobalOrdinal /* = typename Vector<Scalar, LocalOrdinal>::local_ordinal_type*/,
          class Node /* = typename Vector<Scalar, LocalOrdinal, GlobalOrdinal>::node_type*/>
class VectorFactory {
#undef XPETRA_VECTORFACTORY_SHORT
#include "Xpetra_UseShortNames.hpp"

 private:
  //! Private constructor. This is a static class.
  VectorFactory() = default;

 public:
  //! Constructor specifying the number of non-zeros for all rows.
  static Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Build(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>& map, bool zeroOut = true) {
    XPETRA_MONITOR("VectorFactory::Build");

    RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>
        bmap = Teuchos::rcp_dynamic_cast<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>(map);

    if (!bmap.is_null()) {
      return rcp(new Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, zeroOut));
    }

    if (map->lib() == UseTpetra) {
      return rcp(new TpetraVector(map, zeroOut));
    }

    XPETRA_FACTORY_END;
    TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
  }

};  // class VectorFactory

#define XPETRA_VECTORFACTORY_SHORT

}  // namespace Xpetra

#define XPETRA_VECTORFACTORY_SHORT
#endif  // XPETRA_VECTORFACTORY_DECL_HPP
