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

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraExport.hpp"
#endif
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraExport.hpp"
#endif

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

#ifdef HAVE_XPETRA_TPETRA
    if (source->lib() == UseTpetra)
      return rcp(new TpetraExport<LocalOrdinal, GlobalOrdinal, Node>(source, target));
#endif

    XPETRA_FACTORY_ERROR_IF_EPETRA(source->lib());
    XPETRA_FACTORY_END;
  }
};

// we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES))
template <>
class ExportFactory<int, int, EpetraNode> {
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef EpetraNode Node;

 private:
  //! Private constructor. This is a static class.
  ExportFactory() {}

 public:
  static RCP<Export<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source, const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target) {
    XPETRA_MONITOR("ExportFactory::Build");

    TEUCHOS_TEST_FOR_EXCEPTION(source->lib() != target->lib(), Xpetra::Exceptions::RuntimeError, "");

#ifdef HAVE_XPETRA_TPETRA
    if (source->lib() == UseTpetra)
      return rcp(new TpetraExport<LocalOrdinal, GlobalOrdinal, Node>(source, target));
#endif

    if (source->lib() == UseEpetra)
      return rcp(new EpetraExportT<int, Node>(source, target));

    XPETRA_FACTORY_END;
  }
};
#endif

// we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES))
template <>
class ExportFactory<int, long long, EpetraNode> {
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  typedef EpetraNode Node;

 private:
  //! Private constructor. This is a static class.
  ExportFactory() {}

 public:
  static RCP<Export<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source, const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target) {
    XPETRA_MONITOR("ExportFactory::Build");

    TEUCHOS_TEST_FOR_EXCEPTION(source->lib() != target->lib(), Xpetra::Exceptions::RuntimeError, "");

#ifdef HAVE_XPETRA_TPETRA
    if (source->lib() == UseTpetra)
      return rcp(new TpetraExport<LocalOrdinal, GlobalOrdinal, Node>(source, target));
#endif

    if (source->lib() == UseEpetra)
      return rcp(new EpetraExportT<long long, Node>(source, target));

    XPETRA_FACTORY_END;
  }
};
#endif

}  // namespace Xpetra

#define XPETRA_EXPORTFACTORY_SHORT
#endif
