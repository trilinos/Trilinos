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

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraImport.hpp"
#endif
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraImport.hpp"
#endif

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

#ifdef HAVE_XPETRA_TPETRA
    if (source->lib() == UseTpetra)
      return rcp(new TpetraImport<LocalOrdinal, GlobalOrdinal, Node>(source, target, plist));
#endif

    XPETRA_FACTORY_ERROR_IF_EPETRA(source->lib());
    XPETRA_FACTORY_END;
  }
};

// we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES))

// Specialization on LO=GO=int with serial node.
// Used for Epetra and Tpetra
// For any other node definition the general default implementation is used which allows Tpetra only
template <>
class ImportFactory<int, int, EpetraNode> {
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef EpetraNode Node;

 private:
  //! Private constructor. This is a static class.
  ImportFactory() {}

 public:
  static RCP<Import<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source,
                                                               const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target,
                                                               const Teuchos::RCP<Teuchos::ParameterList> &plist = Teuchos::null) {
    XPETRA_MONITOR("ImportFactory::Build");
    TEUCHOS_TEST_FOR_EXCEPTION(source->lib() != target->lib(), Xpetra::Exceptions::RuntimeError, "");

#ifdef HAVE_XPETRA_TPETRA
    if (source->lib() == UseTpetra)
      return rcp(new TpetraImport<LocalOrdinal, GlobalOrdinal, Node>(source, target, plist));
#endif

    if (source->lib() == UseEpetra)
      return rcp(new EpetraImportT<int, Node>(source, target));

    XPETRA_FACTORY_END;
  }
};
#endif

// we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES))
template <>
class ImportFactory<int, long long, EpetraNode> {
  typedef int LocalOrdinal;
  typedef long long GlobalOrdinal;
  typedef EpetraNode Node;

 private:
  //! Private constructor. This is a static class.
  ImportFactory() {}

 public:
  static RCP<Import<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source,
                                                               const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target,
                                                               const Teuchos::RCP<Teuchos::ParameterList> &plist = Teuchos::null) {
    XPETRA_MONITOR("ImportFactory::Build");
    TEUCHOS_TEST_FOR_EXCEPTION(source->lib() != target->lib(), Xpetra::Exceptions::RuntimeError, "");

#ifdef HAVE_XPETRA_TPETRA
    if (source->lib() == UseTpetra)
      return rcp(new TpetraImport<LocalOrdinal, GlobalOrdinal, Node>(source, target, plist));
#endif

    if (source->lib() == UseEpetra)
      return rcp(new EpetraImportT<long long, Node>(source, target));

    XPETRA_FACTORY_END;
  }
};
#endif
}  // namespace Xpetra

#define XPETRA_IMPORTFACTORY_SHORT
#endif
