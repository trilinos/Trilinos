// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_EXPORTFACTORY_HPP
#define XPETRA_EXPORTFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Export.hpp"

#include "Xpetra_TpetraExport.hpp"

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class ExportFactory {
 private:
  //! Private constructor. This is a static class.
  ExportFactory() {}

 public:
  //! Constructor specifying the number of non-zeros for all rows.
  static RCP<Export<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source, const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target) {
    XPETRA_MONITOR("ExportFactory::Build");
    TEUCHOS_TEST_FOR_EXCEPTION(source->lib() != target->lib(), Xpetra::Exceptions::RuntimeError, "");

    if (source->lib() == UseTpetra)
      return rcp(new TpetraExport<LocalOrdinal, GlobalOrdinal, Node>(source, target));

    XPETRA_FACTORY_ERROR_IF_EPETRA(source->lib());
    XPETRA_FACTORY_END;
  }
};

// we need the Epetra specialization only if Epetra is enabled

// we need the Epetra specialization only if Epetra is enabled

}  // namespace Xpetra

#define XPETRA_EXPORTFACTORY_SHORT
#endif
