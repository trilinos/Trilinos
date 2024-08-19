// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_TPETRAMAP_DEF_HPP
#define XPETRA_TPETRAMAP_DEF_HPP

#include "Xpetra_TpetraMap_decl.hpp"

namespace Xpetra {

//! @name Constructor/Destructor Methods
//@{

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraMap(global_size_t numGlobalElements,
              GlobalOrdinal indexBase,
              const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
              LocalGlobal lg)
  : map_(Teuchos::rcp(new Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>(numGlobalElements,
                                                                         indexBase, comm,
                                                                         toTpetra(lg)))) {}

//! Constructor with a user-defined contiguous distribution.

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraMap(global_size_t numGlobalElements,
              size_t numLocalElements,
              GlobalOrdinal indexBase,
              const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
  : map_(Teuchos::rcp(new Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>(numGlobalElements,
                                                                         numLocalElements,
                                                                         indexBase, comm))) {}

//! Constructor with user-defined arbitrary (possibly noncontiguous) distribution.

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraMap(global_size_t numGlobalElements,
              const Teuchos::ArrayView<const GlobalOrdinal> &elementList,
              GlobalOrdinal indexBase,
              const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
  : map_(Teuchos::rcp(new Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>(numGlobalElements,
                                                                         elementList,
                                                                         indexBase,
                                                                         comm))) {}

//! Constructor with user-defined arbitrary (possibly noncontiguous) distribution passed as a Kokkos::View.
template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::
    TpetraMap(global_size_t numGlobalElements,
              const Kokkos::View<const GlobalOrdinal *, typename Node::device_type> &indexList,
              GlobalOrdinal indexBase,
              const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
  : map_(Teuchos::rcp(new Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>(numGlobalElements,
                                                                         indexList,
                                                                         indexBase,
                                                                         comm))) {}

//! Destructor.
template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::~TpetraMap() {}

//@}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getGlobalNumElements() const {
  XPETRA_MONITOR("TpetraMap::getGlobalNumElements");
  return map_->getGlobalNumElements();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getLocalNumElements() const {
  XPETRA_MONITOR("TpetraMap::getLocalNumElements");
  return map_->getLocalNumElements();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getIndexBase() const {
  XPETRA_MONITOR("TpetraMap::getIndexBase");
  return map_->getIndexBase();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getMinLocalIndex() const {
  XPETRA_MONITOR("TpetraMap::getMinLocalIndex");
  return map_->getMinLocalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getMaxLocalIndex() const {
  XPETRA_MONITOR("TpetraMap::getMaxLocalIndex");
  return map_->getMaxLocalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getMinGlobalIndex() const {
  XPETRA_MONITOR("TpetraMap::getMinGlobalIndex");
  return map_->getMinGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getMaxGlobalIndex() const {
  XPETRA_MONITOR("TpetraMap::getMaxGlobalIndex");
  return map_->getMaxGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getMinAllGlobalIndex() const {
  XPETRA_MONITOR("TpetraMap::getMinAllGlobalIndex");
  return map_->getMinAllGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getMaxAllGlobalIndex() const {
  XPETRA_MONITOR("TpetraMap::getMaxAllGlobalIndex");
  return map_->getMaxAllGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getLocalElement(GlobalOrdinal globalIndex) const {
  XPETRA_MONITOR("TpetraMap::getLocalElement");
  return map_->getLocalElement(globalIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getGlobalElement(LocalOrdinal localIndex) const {
  XPETRA_MONITOR("TpetraMap::getGlobalElement");
  return map_->getGlobalElement(localIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LookupStatus TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> &GIDList, const Teuchos::ArrayView<int> &nodeIDList, const Teuchos::ArrayView<LocalOrdinal> &LIDList) const {
  XPETRA_MONITOR("TpetraMap::getRemoteIndexList");
  return toXpetra(map_->getRemoteIndexList(GIDList, nodeIDList, LIDList));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LookupStatus TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal> &GIDList, const Teuchos::ArrayView<int> &nodeIDList) const {
  XPETRA_MONITOR("TpetraMap::getRemoteIndexList");
  return toXpetra(map_->getRemoteIndexList(GIDList, nodeIDList));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayView<const GlobalOrdinal> TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getLocalElementList() const {
  XPETRA_MONITOR("TpetraMap::getLocalElementList");
  return map_->getLocalElementList();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
typename Map<LocalOrdinal, GlobalOrdinal, Node>::global_indices_array_device_type TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getMyGlobalIndicesDevice() const {
  XPETRA_MONITOR("TpetraMap::getMyGlobalIndicesDevice");
  return map_->getMyGlobalIndicesDevice();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::isNodeLocalElement(LocalOrdinal localIndex) const {
  XPETRA_MONITOR("TpetraMap::isNodeLocalElement");
  return map_->isNodeLocalElement(localIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::isNodeGlobalElement(GlobalOrdinal globalIndex) const {
  XPETRA_MONITOR("TpetraMap::isNodeGlobalElement");
  return map_->isNodeGlobalElement(globalIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::isContiguous() const {
  XPETRA_MONITOR("TpetraMap::isContiguous");
  return map_->isContiguous();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::isDistributed() const {
  XPETRA_MONITOR("TpetraMap::isDistributed");
  return map_->isDistributed();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::isCompatible(const Map<LocalOrdinal, GlobalOrdinal, Node> &map) const {
  XPETRA_MONITOR("TpetraMap::isCompatible");
  return map_->isCompatible(toTpetra(map));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::isSameAs(const Map<LocalOrdinal, GlobalOrdinal, Node> &map) const {
  XPETRA_MONITOR("TpetraMap::isSameAs");
  return map_->isSameAs(toTpetra(map));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Teuchos::Comm<int> > TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getComm() const {
  XPETRA_MONITOR("TpetraMap::getComm");
  return map_->getComm();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
std::string TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::description() const {
  XPETRA_MONITOR("TpetraMap::description");
  return map_->description();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
  XPETRA_MONITOR("TpetraMap::describe");
  map_->describe(out, verbLevel);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::removeEmptyProcesses() const {
  return toXpetra(map_->removeEmptyProcesses());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::replaceCommWithSubset(const Teuchos::RCP<const Teuchos::Comm<int> > &newComm) const {
  return toXpetra(map_->replaceCommWithSubset(newComm));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::TpetraMap(const Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > &map)
  : map_(map) {}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
UnderlyingLib TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::lib() const { return UseTpetra; }

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getTpetra_Map() const { return map_; }

template <class LocalOrdinal, class GlobalOrdinal, class Node>
typename Map<LocalOrdinal, GlobalOrdinal, Node>::local_map_type TpetraMap<LocalOrdinal, GlobalOrdinal, Node>::getLocalMap() const {
  return map_->getLocalMap();
}

}  // namespace Xpetra

// TODO: remove?
//!  Returns true if \c map is identical to this map. Implemented in TpetraMap::isSameAs().
template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool operator==(const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> &map1, const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> &map2) {
  XPETRA_MONITOR("TpetraMap==TpetraMap");
  return map1.isSameAs(map2);
}

//! Returns true if \c map is not identical to this map. Implemented in TpetraMap::isSameAs().
template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool operator!=(const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> &map1, const Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> &map2) {
  XPETRA_MONITOR("TpetraMap!=TpetraMap");
  return !map1.isSameAs(map2);
}

#endif  // XPETRA_TPETRAMAP_DEF_HPP
