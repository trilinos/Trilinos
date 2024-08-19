// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_TPETRAIMPORT_DEF_HPP
#define XPETRA_TPETRAIMPORT_DEF_HPP
#include "Xpetra_TpetraConfigDefs.hpp"

#include "Xpetra_Import.hpp"
#include "Xpetra_TpetraImport_decl.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_TpetraMap.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Distributor.hpp"

namespace Xpetra {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::TpetraImport(const Teuchos::RCP<const map_type> &source, const Teuchos::RCP<const map_type> &target)
  : import_(Teuchos::rcp(new Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(source), toTpetra(target)))) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::TpetraImport(const Teuchos::RCP<const map_type> &source, const Teuchos::RCP<const map_type> &target, const Teuchos::RCP<Teuchos::ParameterList> &plist)
  : import_(Teuchos::rcp(new Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(source), toTpetra(target), plist))) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::TpetraImport(const Import<LocalOrdinal, GlobalOrdinal, Node> &import)
  : import_(Teuchos::rcp(new Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(import)))) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::~TpetraImport() {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node> >
TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::createRemoteOnlyImport(const Teuchos::RCP<const map_type> &remoteTarget) const {
  Teuchos::RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > newImport = import_->createRemoteOnlyImport(toTpetra(remoteTarget));
  return Teuchos::rcp(new TpetraImport<LocalOrdinal, GlobalOrdinal, Node>(newImport));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getNumSameIDs() const {
  XPETRA_MONITOR("TpetraImport::getNumSameIDs");
  return import_->getNumSameIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getNumPermuteIDs() const {
  XPETRA_MONITOR("TpetraImport::getNumPermuteIDs");
  return import_->getNumPermuteIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal> TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getPermuteFromLIDs() const {
  XPETRA_MONITOR("TpetraImport::getPermuteFromLIDs");
  return import_->getPermuteFromLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal> TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getPermuteToLIDs() const {
  XPETRA_MONITOR("TpetraImport::getPermuteToLIDs");
  return import_->getPermuteToLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getNumRemoteIDs() const {
  XPETRA_MONITOR("TpetraImport::getNumRemoteIDs");
  return import_->getNumRemoteIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::setDistributorParameters(const Teuchos::RCP<Teuchos::ParameterList> params) const {
  XPETRA_MONITOR("TpetraImport::setDistributorParameters");
  import_->getDistributor().setParameterList(params);
  auto revDistor = import_->getDistributor().getReverse(false);
  if (!revDistor.is_null())
    revDistor->setParameterList(params);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal> TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getRemoteLIDs() const {
  XPETRA_MONITOR("TpetraImport::getRemoteLIDs");
  return import_->getRemoteLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getNumExportIDs() const {
  XPETRA_MONITOR("TpetraImport::getNumExportIDs");
  return import_->getNumExportIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal> TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getExportLIDs() const {
  XPETRA_MONITOR("TpetraImport::getExportLIDs");
  return import_->getExportLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const int> TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getExportPIDs() const {
  XPETRA_MONITOR("TpetraImport::getExportPIDs");
  return import_->getExportPIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getSourceMap() const {
  XPETRA_MONITOR("TpetraImport::getSourceMap");
  return toXpetra(import_->getSourceMap());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getTargetMap() const {
  XPETRA_MONITOR("TpetraImport::getTargetMap");
  return toXpetra(import_->getTargetMap());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::print(std::ostream &os) const {
  XPETRA_MONITOR("TpetraImport::print");
  import_->print(os);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::TpetraImport(const RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > &import)
  : import_(import) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > TpetraImport<LocalOrdinal, GlobalOrdinal, Node>::getTpetra_Import() const { return import_; }

}  // namespace Xpetra

#endif
