// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKEDMAP_DEF_HPP
#define TPETRA_BLOCKEDMAP_DEF_HPP

#include "Tpetra_BlockedMap_decl.hpp"

#include <stdexcept>

#include <Teuchos_as.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_TestForException.hpp>

#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"

namespace Tpetra {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMap() {
  bThyraMode_ = false;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMap(const Teuchos::RCP<const map_type>& fullmap,
               const std::vector<Teuchos::RCP<const map_type>>& maps,
               bool bThyraMode) {
  bThyraMode_ = bThyraMode;

  if (bThyraMode == false) {
    // use Xpetra-style numbering for sub-block maps
    // That is, all sub-block maps have unique GIDs which may not be contiguous and start with GIDs different than zero.

    // plausibility check
    size_t numAllElements = 0;
    for (size_t v = 0; v < maps.size(); ++v) {
      numAllElements += maps[v]->getGlobalNumElements();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(fullmap->getGlobalNumElements() != numAllElements,
                               std::logic_error,
                               "logic error. full map and sub maps have not same number of elements ("
                                   << fullmap->getGlobalNumElements() << " versus " << numAllElements
                                   << "). We cannot build BlockedMap with Xpetra-style numbering. Please make sure that you want "
                                      "Xpetra-style numbering instead of Thyra-style numbering.");

    fullmap_ = fullmap;
    maps_    = maps;
  } else {
    // use Thyra-style numbering for sub-block maps
    // That is, all sub-block maps start with zero as GID and are contiguous

    // plausibility check
    for (size_t v = 0; v < maps.size(); ++v) {
      TEUCHOS_TEST_FOR_EXCEPTION(maps[v]->getMinAllGlobalIndex() != 0,
                                 std::logic_error,
                                 "logic error. When using Thyra-style numbering all sub-block maps must start with zero as GID. Map block "
                                     << v << " starts with GID " << maps[v]->getMinAllGlobalIndex());
    }

    // store submaps in Thyra-style ordering
    thyraMaps_ = maps;

    // get offsets
    std::vector<GlobalOrdinal> gidOffsets(maps.size(), 0);
    for (size_t v = 1; v < maps.size(); ++v) {
      gidOffsets[v] = maps[v - 1]->getMaxAllGlobalIndex() + gidOffsets[v - 1] + 1;
    }

    // build submaps
    maps_.resize(maps.size());
    std::vector<GlobalOrdinal> fullMapGids;
    const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
    for (size_t v = 0; v < maps.size(); ++v) {
      size_t myNumElements = maps[v]->getLocalNumElements();
      std::vector<GlobalOrdinal> subMapGids(myNumElements, 0);
      for (LocalOrdinal l = 0; l < Teuchos::as<LocalOrdinal>(myNumElements); ++l) {
        GlobalOrdinal myGid = maps[v]->getGlobalElement(l);
        subMapGids[l]       = myGid + gidOffsets[v];
        fullMapGids.push_back(myGid + gidOffsets[v]);
      }

      Teuchos::ArrayView<const GlobalOrdinal> subMapGidsView(subMapGids.data(), subMapGids.size());

      Teuchos::RCP<const map_type> mySubMap = Teuchos::rcp(
          new map_type(INVALID, subMapGidsView, maps[v]->getIndexBase(), maps[v]->getComm()));
      maps_[v] = mySubMap;
    }

    Teuchos::ArrayView<const GlobalOrdinal> fullMapGidsView(fullMapGids.data(), fullMapGids.size());

    fullmap_ = Teuchos::rcp(
        new map_type(INVALID, fullMapGidsView, fullmap->getIndexBase(), fullmap->getComm()));

    // plausibility check
    size_t numAllElements = 0;
    for (size_t v = 0; v < maps_.size(); ++v) {
      numAllElements += maps_[v]->getGlobalNumElements();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
        fullmap_->getGlobalNumElements() != numAllElements,
        std::logic_error,
        "logic error. full map and sub maps have not same number of elements. This cannot be. Please report the bug to the Tpetra developers!");
  }

  // build importers for sub maps
  importers_.resize(maps_.size());
  for (unsigned i = 0; i < maps_.size(); ++i) {
    if (maps[i] != Teuchos::null) {
      importers_[i] = Teuchos::rcp(new import_type(fullmap_, maps_[i]));
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
      CheckConsistency() == false, std::logic_error, "logic error. full map and sub maps are inconsistently distributed over the processors.");
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMap(const std::vector<Teuchos::RCP<const map_type>>& maps,
               const std::vector<Teuchos::RCP<const map_type>>& thyramaps) {
  bThyraMode_ = true;

  // plausibility check
  TEUCHOS_TEST_FOR_EXCEPTION(thyramaps.size() != maps.size(), std::logic_error, "logic error. The number of submaps must be identical!");
  for (size_t v = 0; v < thyramaps.size(); ++v) {
    TEUCHOS_TEST_FOR_EXCEPTION(thyramaps[v]->getMinAllGlobalIndex() != 0,
                               std::logic_error,
                               "logic error. When using Thyra-style numbering all sub-block maps must start with zero as GID.");

    TEUCHOS_TEST_FOR_EXCEPTION(thyramaps[v]->getLocalNumElements() != maps[v]->getLocalNumElements(),
                               std::logic_error,
                               "logic error. The size of the submaps must be identical (same distribution, just different GIDs)");
  }

  // store user-provided maps and thyramaps
  thyraMaps_ = thyramaps;
  maps_      = maps;
  fullmap_   = this->concatenateMaps(maps);

  // plausibility check
  size_t numAllElements = 0;
  for (size_t v = 0; v < maps_.size(); ++v) {
    numAllElements += maps_[v]->getGlobalNumElements();
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
      fullmap_->getGlobalNumElements() != numAllElements,
      std::logic_error,
      "logic error. full map and sub maps have not same number of elements. This cannot be. Please report the bug to the Tpetra developers!");

  // build importers for sub maps
  importers_.resize(maps_.size());
  for (unsigned i = 0; i < maps_.size(); ++i) {
    if (maps[i] != Teuchos::null) {
      importers_[i] = Teuchos::rcp(new import_type(fullmap_, maps_[i]));
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
      CheckConsistency() == false, std::logic_error, "logic error. full map and sub maps are inconsistently distributed over the processors.");
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMap(const BlockedMap& input) {
  bThyraMode_ = input.getThyraMode();
  fullmap_    = Teuchos::null;
  maps_.resize(input.getNumMaps(), Teuchos::null);
  thyraMaps_.resize(input.getNumMaps(), Teuchos::null);
  this->assign(input);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    ~BlockedMap() {
  // make sure all RCP's are freed
  for (size_t v = 0; v < maps_.size(); ++v) {
    maps_[v] = Teuchos::null;
    if (bThyraMode_ == true)
      thyraMaps_[v] = Teuchos::null;
    importers_[v] = Teuchos::null;
  }

  fullmap_ = Teuchos::null;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalNumElements() const {
  return fullmap_->getGlobalNumElements();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::getLocalNumElements() const {
  return fullmap_->getLocalNumElements();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getIndexBase() const {
  return fullmap_->getIndexBase();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMinLocalIndex() const {
  return fullmap_->getMinLocalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMaxLocalIndex() const {
  return fullmap_->getMaxLocalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMinGlobalIndex() const {
  return fullmap_->getMinGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMaxGlobalIndex() const {
  return fullmap_->getMaxGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMinAllGlobalIndex() const {
  return fullmap_->getMinAllGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMaxAllGlobalIndex() const {
  return fullmap_->getMaxAllGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getLocalElement(GlobalOrdinal globalIndex) const {
  return fullmap_->getLocalElement(globalIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalElement(LocalOrdinal localIndex) const {
  return fullmap_->getGlobalElement(localIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LookupStatus
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal>& /* GIDList    */,
                       const Teuchos::ArrayView<int>& /* nodeIDList */,
                       const Teuchos::ArrayView<LocalOrdinal>& /* LIDList    */) const {
  throw std::runtime_error("BlockedMap::getRemoteIndexList: routine not implemented.");
  TEUCHOS_UNREACHABLE_RETURN(IDNotPresent);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LookupStatus
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal>& /* GIDList */,
                       const Teuchos::ArrayView<int>& /* nodeIDList */) const {
  throw std::runtime_error("BlockedMap::getRemoteIndexList: routine not implemented.");
  TEUCHOS_UNREACHABLE_RETURN(IDNotPresent);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayView<const GlobalOrdinal>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getLocalElementList() const {
  return fullmap_->getLocalElementList();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
typename BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::global_indices_array_device_type
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMyGlobalIndicesDevice() const {
  return fullmap_->getMyGlobalIndicesDevice();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isNodeLocalElement(LocalOrdinal localIndex) const {
  return fullmap_->isNodeLocalElement(localIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isNodeGlobalElement(GlobalOrdinal globalIndex) const {
  return fullmap_->isNodeGlobalElement(globalIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isContiguous() const {
  throw std::runtime_error("BlockedMap::isContiguous: routine not implemented.");
  TEUCHOS_UNREACHABLE_RETURN(false);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isDistributed() const {
  return fullmap_->isDistributed();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isCompatible(const map_type& map) const {
  // A plain (non-blocked) map is only compatible with a single-block BlockedMap.
  if (getNumMaps() != 1)
    return false;
  return getMap(0, bThyraMode_)->isCompatible(map);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isCompatible(const BlockedMap& rcpBMap) const {
  if (getNumMaps() != rcpBMap.getNumMaps())
    return false;

  for (size_t v = 0; v < maps_.size(); ++v) {
    bool bSame = getMap(v, false)->isCompatible(*(rcpBMap.getMap(v, false)));
    if (bSame == false)
      return false;
    if (bThyraMode_ && rcpBMap.getThyraMode()) {
      bSame = getMap(v, true)->isCompatible(*(rcpBMap.getMap(v, true)));
      if (bSame == false)
        return false;
    }
  }
  return true;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isSameAs(const map_type& map) const {
  // If this is a blocked map with > 1 blocks but "map" is a plain map they can't be the same
  if (this->getNumMaps() > 1) {
    return false;
  }

  // special case: this is a single blocked map and "map" is a plain map object
  return getMap(0, bThyraMode_)->isSameAs(map);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isSameAs(const BlockedMap& rcpBMap) const {
  if (getNumMaps() != rcpBMap.getNumMaps())
    return false;

  for (size_t v = 0; v < maps_.size(); ++v) {
    bool bSame = getMap(v, false)->isSameAs(*(rcpBMap.getMap(v, false)));
    if (bSame == false) {
      return false;
    }
    if (bThyraMode_ && rcpBMap.getThyraMode()) {
      bSame = getMap(v, true)->isSameAs(*(rcpBMap.getMap(v, true)));
      if (bSame == false) {
        return false;
      }
    }
  }
  return true;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Teuchos::Comm<int>>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getComm() const {
  return fullmap_->getComm();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>&
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
operator=(const BlockedMap& rhs) {
  assign(rhs);  // dispatch to protected virtual method
  return *this;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::getThyraMode() const {
  return bThyraMode_;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getNumMaps() const {
  return maps_.size();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::map_type>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMap(size_t i, bool bThyraMode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(i >= getNumMaps(),
                             std::runtime_error,
                             "BlockedMap::getMap: tried to access block " << i << ", but BlockedMap has only " << getNumMaps()
                                                                          << " blocks! Block indices must be between 0 and " << getNumMaps() - 1
                                                                          << ".");
  if (bThyraMode_ == true && bThyraMode == true) {
    return thyraMaps_[i];
  }

  TEUCHOS_TEST_FOR_EXCEPTION(bThyraMode_ == false && bThyraMode == true,
                             std::runtime_error,
                             "BlockedMap::getMap: cannot return sub map in Thyra-style numbering if BlockedMap object is not created using "
                             "Thyra-style numbered submaps.");
  return maps_[i];
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::import_type>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getImporter(size_t i) const {
  TEUCHOS_TEST_FOR_EXCEPTION(i >= getNumMaps(),
                             std::runtime_error,
                             "BlockedMap::getImporter: tried to access block " << i << ", but BlockedMap has only " << getNumMaps()
                                                                               << " blocks! Block indices must be between 0 and " << getNumMaps() - 1
                                                                               << ".");
  return importers_[i];
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::map_type>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getFullMap() const {
  return fullmap_;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMapIndexForGID(GlobalOrdinal gid) const {
  for (size_t i = 0; i < getNumMaps(); i++)
    if (getMap(i)->isNodeGlobalElement(gid) == true)
      return i;

  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error, "getMapIndexForGID: GID " << gid << " is not contained by a map in mapextractor.");
  return 0;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    setBlockedSubMap(size_t i, const Teuchos::RCP<const BlockedMap>& submap) {
  TEUCHOS_TEST_FOR_EXCEPTION(i >= getNumMaps(), std::runtime_error,
                             "BlockedMap::setBlockedSubMap: block " << i << " is out of range (" << getNumMaps() << " maps).");
  if (nestedMaps_.size() != getNumMaps())
    nestedMaps_.resize(getNumMaps());
  nestedMaps_[i] = submap;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getBlockedSubMap(size_t i) const {
  if (i >= nestedMaps_.size())
    return Teuchos::null;
  return nestedMaps_[i];
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isBlocked(size_t i) const {
  return !getBlockedSubMap(i).is_null();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
std::string
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    description() const {
  return std::string("BlockedMap");
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
  out << "------------- Blocked Map -----------" << std::endl;
  out << description() << std::endl;
  out << "Thyra mode: " << getThyraMode() << std::endl;
  out << "No of submaps: " << getNumMaps() << std::endl;
  Teuchos::OSTab tab(out);
  for (size_t r = 0; r < getNumMaps(); r++) {
    out << "MAP " << r << "/" << getNumMaps() - 1 << std::endl;
    getMap(r, false)->describe(out, verbLevel);
  }
  if (getThyraMode() == true) {
    for (size_t r = 0; r < getNumMaps(); r++) {
      out << "Thyra MAP " << r << "/" << getNumMaps() - 1 << std::endl;
      getMap(r, true)->describe(out, verbLevel);
    }
  }
  out << "-------------------------------------" << std::endl;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    assign(const BlockedMap& input) {
  // Tpetra::Map objects are immutable and reference counted, so a deep copy is
  // unnecessary: sharing the underlying RCPs is semantically identical to (and
  // cheaper than) rebuilding the maps as the Xpetra implementation does.
  bThyraMode_ = input.bThyraMode_;
  fullmap_    = input.fullmap_;
  maps_       = input.maps_;
  thyraMaps_  = input.thyraMaps_;
  importers_  = input.importers_;
  nestedMaps_ = input.nestedMaps_;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::map_type>
BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    concatenateMaps(const std::vector<Teuchos::RCP<const map_type>>& subMaps) {
  // merge submaps to global map
  std::vector<GlobalOrdinal> gids;
  for (size_t tt = 0; tt < subMaps.size(); ++tt) {
    Teuchos::RCP<const map_type> subMap = subMaps[tt];
    Teuchos::ArrayView<const GlobalOrdinal> subMapGids = subMap->getLocalElementList();
    gids.insert(gids.end(), subMapGids.begin(), subMapGids.end());
  }

  const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
  Teuchos::ArrayView<const GlobalOrdinal> gidsView(gids.data(), gids.size());

  Teuchos::RCP<const map_type> fullMap = Teuchos::rcp(
      new map_type(INVALID, gidsView, subMaps[0]->getIndexBase(), subMaps[0]->getComm()));

  return fullMap;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMap<LocalOrdinal, GlobalOrdinal, Node>::
    CheckConsistency() const {
  const Teuchos::RCP<const map_type> fullMap = getFullMap();

  for (size_t i = 0; i < getNumMaps(); i++) {
    const Teuchos::RCP<const map_type> map = getMap(i);

    Teuchos::ArrayView<const GlobalOrdinal> mapGids = map->getLocalElementList();
    for (typename Teuchos::ArrayView<const GlobalOrdinal>::const_iterator it = mapGids.begin(); it != mapGids.end(); it++) {
      if (fullMap->isNodeGlobalElement(*it) == false) {
        return false;  // Global ID (*it) not found locally on this proc in fullMap -> error
      }
    }
  }
  return true;
}

}  // namespace Tpetra

//
// Explicit instantiation macro
//
#define TPETRA_BLOCKEDMAP_INSTANT(LO, GO, NODE) \
  template class BlockedMap<LO, GO, NODE>;

#endif  // TPETRA_BLOCKEDMAP_DEF_HPP
