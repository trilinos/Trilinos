// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_IMPORTFACTORY_HPP
#define XPETRA_IMPORTFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Import.hpp"

#include "Xpetra_TpetraImport.hpp"

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class ImportFactory {
 private:
  //! Private constructor. This is a static class.
  ImportFactory() {}

 public:
  //! Constructor specifying the number of non-zeros for all rows.
  static RCP<Import<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source,
                                                               const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target,
                                                               const Teuchos::RCP<Teuchos::ParameterList> &plist = Teuchos::null) {
    XPETRA_MONITOR("ImportFactory::Build");

    TEUCHOS_TEST_FOR_EXCEPTION(source->lib() != target->lib(), Xpetra::Exceptions::RuntimeError, "");

    if (source->lib() == UseTpetra)
      return rcp(new TpetraImport<LocalOrdinal, GlobalOrdinal, Node>(source, target, plist));

    XPETRA_FACTORY_END;
  }
};

}  // namespace Xpetra

#define XPETRA_IMPORTFACTORY_SHORT
#endif
