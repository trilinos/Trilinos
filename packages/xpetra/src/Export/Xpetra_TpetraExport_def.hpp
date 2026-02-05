// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_TPETRAEXPORT_DEF_HPP
#define XPETRA_TPETRAEXPORT_DEF_HPP

#include "Xpetra_TpetraExport_decl.hpp"
#include "Tpetra_Distributor.hpp"

namespace Xpetra {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraExport(const Teuchos::RCP<const map_type>& source,
                 const Teuchos::RCP<const map_type>& target)
  : export_(Teuchos::rcp(new Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(source), toTpetra(target)))) {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraExport(const Teuchos::RCP<const map_type>& source,
                 const Teuchos::RCP<const map_type>& target,
                 const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : export_(Teuchos::rcp(new Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(source), toTpetra(target), plist))) {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraExport(const Export<LocalOrdinal, GlobalOrdinal, Node>& rhs)
  : export_(Teuchos::rcp(new Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>(toTpetra(rhs)))) {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    ~TpetraExport() {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getNumSameIDs() const {
  XPETRA_MONITOR("TpetraExport::getNumSameIDs");
  return export_->getNumSameIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getNumPermuteIDs() const {
  XPETRA_MONITOR("TpetraExport::getNumPermuteIDs");
  return export_->getNumPermuteIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getPermuteFromLIDs() const {
  XPETRA_MONITOR("TpetraExport::getPermuteFromLIDs");
  return export_->getPermuteFromLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getPermuteToLIDs() const {
  XPETRA_MONITOR("TpetraExport::getPermuteToLIDs");
  return export_->getPermuteToLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getNumRemoteIDs() const {
  XPETRA_MONITOR("TpetraExport::getNumRemoteIDs");
  return export_->getNumRemoteIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getRemoteLIDs() const {
  XPETRA_MONITOR("TpetraExport::getRemoteLIDs");
  return export_->getRemoteLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getNumExportIDs() const {
  XPETRA_MONITOR("TpetraExport::getNumExportIDs");
  return export_->getNumExportIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const LocalOrdinal>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getExportLIDs() const {
  XPETRA_MONITOR("TpetraExport::getExportLIDs");
  return export_->getExportLIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayView<const int>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getExportPIDs() const {
  XPETRA_MONITOR("TpetraExport::getExportPIDs");
  return export_->getExportPIDs();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getSourceMap() const {
  XPETRA_MONITOR("TpetraExport::getSourceMap");
  return toXpetra(export_->getSourceMap());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getTargetMap() const {
  XPETRA_MONITOR("TpetraExport::getTargetMap");
  return toXpetra(export_->getTargetMap());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    setDistributorParameters(const Teuchos::RCP<Teuchos::ParameterList> params) const {
  XPETRA_MONITOR("TpetraExport::setDistributorParameters");
  export_->getDistributor().setParameterList(params);
  auto revDistor = export_->getDistributor().getReverse(false);
  if (!revDistor.is_null())
    revDistor->setParameterList(params);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    print(std::ostream& os) const {
  XPETRA_MONITOR("TpetraExport::print");
  export_->print(os);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraExport(
        const RCP<const Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>>& exp)
  : export_(exp) {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>>
TpetraExport<LocalOrdinal, GlobalOrdinal, Node>::
    getTpetra_Export() const {
  return export_;
}

}  // namespace Xpetra

#endif  // XPETRA_TPETRAEXPORT_DEF_HPP
