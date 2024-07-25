// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_STRIDEDMAP_DEF_HPP
#define XPETRA_STRIDEDMAP_DEF_HPP

#include "Xpetra_StridedMap.hpp"

#include <Teuchos_OrdinalTraits.hpp>

#include "Xpetra_Exceptions.hpp"
#include "Xpetra_MapFactory.hpp"

namespace Xpetra {

template <class LocalOrdinal, class GlobalOrdinal, class Node>
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    StridedMap(UnderlyingLib xlib,
               global_size_t numGlobalElements,
               GlobalOrdinal indexBase,
               std::vector<size_t>& stridingInfo,
               const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
               LocalOrdinal stridedBlockId,  // FIXME (mfh 03 Sep 2014) This breaks for unsigned LocalOrdinal
               GlobalOrdinal offset,
               LocalGlobal lg)
  : stridingInfo_(stridingInfo)
  , stridedBlockId_(stridedBlockId)
  , offset_(offset)
  , indexBase_(indexBase) {
  using MapFactory_t = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>;

  size_t blkSize = getFixedBlockSize();

  TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo.size() == 0,
                             Exceptions::RuntimeError,
                             "StridedMap::StridedMap: stridingInfo not valid: stridingInfo.size() = 0?");

  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalElements == Teuchos::OrdinalTraits<global_size_t>::invalid(),
                             std::invalid_argument,
                             "StridedMap::StridedMap: numGlobalElements is invalid");

  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalElements % blkSize != 0,
                             Exceptions::RuntimeError,
                             "StridedMap::StridedMap: stridingInfo not valid: getFixedBlockSize "
                             "is not an integer multiple of numGlobalElements.");

  if (stridedBlockId != -1) {
    TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo.size() < static_cast<size_t>(stridedBlockId),
                               Exceptions::RuntimeError,
                               "StridedTpetraMap::StridedTpetraMap: "
                               "stridedBlockId > stridingInfo.size()");
  }

  // Try to create a shortcut
  if (blkSize != 1 || offset_ != 0) {
    // check input data and reorganize map
    global_size_t numGlobalNodes = numGlobalElements / blkSize;

    // build an equally distributed node map
    RCP<Map> nodeMap            = MapFactory_t::Build(xlib, numGlobalNodes, indexBase, comm, lg);
    global_size_t numLocalNodes = nodeMap->getLocalNumElements();

    // translate local node ids to local dofs
    size_t nStridedOffset = 0;
    size_t nDofsPerNode   = blkSize;  // dofs per node for local striding block
    if (stridedBlockId > -1) {
      for (int j = 0; j < stridedBlockId; j++) {
        nStridedOffset += stridingInfo_[j];
      }

      nDofsPerNode      = stridingInfo_[stridedBlockId];
      numGlobalElements = numGlobalNodes * Teuchos::as<global_size_t>(nDofsPerNode);
    }
    size_t numLocalElements = numLocalNodes * Teuchos::as<size_t>(nDofsPerNode);

    std::vector<GlobalOrdinal> dofgids(numLocalElements);
    for (LocalOrdinal i = 0; i < Teuchos::as<LocalOrdinal>(numLocalNodes); i++) {
      GlobalOrdinal nodeGID = nodeMap->getGlobalElement(i);

      for (size_t j = 0; j < nDofsPerNode; j++) {
        dofgids[i * nDofsPerNode + j] = indexBase_ + offset_ + (nodeGID - indexBase_) * Teuchos::as<GlobalOrdinal>(blkSize) + Teuchos::as<GlobalOrdinal>(nStridedOffset + j);
      }
    }

    map_ = MapFactory_t::Build(xlib, numGlobalElements, dofgids, indexBase, comm);

    if (stridedBlockId == -1) {
      TEUCHOS_TEST_FOR_EXCEPTION(getLocalNumElements() != Teuchos::as<size_t>(nodeMap->getLocalNumElements() * nDofsPerNode),
                                 Exceptions::RuntimeError,
                                 "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");

      TEUCHOS_TEST_FOR_EXCEPTION(getGlobalNumElements() != Teuchos::as<size_t>(nodeMap->getGlobalNumElements() * nDofsPerNode),
                                 Exceptions::RuntimeError,
                                 "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");
    } else {
      size_t nDofsInStridedBlock = stridingInfo[stridedBlockId];
      TEUCHOS_TEST_FOR_EXCEPTION(getLocalNumElements() != Teuchos::as<size_t>(nodeMap->getLocalNumElements() * nDofsInStridedBlock),
                                 Exceptions::RuntimeError,
                                 "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");

      TEUCHOS_TEST_FOR_EXCEPTION(getGlobalNumElements() != Teuchos::as<size_t>(nodeMap->getGlobalNumElements() * nDofsInStridedBlock),
                                 Exceptions::RuntimeError,
                                 "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");
    }
  } else {
    map_ = MapFactory_t::Build(xlib, numGlobalElements, indexBase, comm, lg);
  }

  TEUCHOS_TEST_FOR_EXCEPTION(CheckConsistency() == false, Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: CheckConsistency() == false");
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    StridedMap(UnderlyingLib xlib,
               global_size_t numGlobalElements,
               size_t numLocalElements,
               GlobalOrdinal indexBase,
               std::vector<size_t>& stridingInfo,
               const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
               LocalOrdinal stridedBlockId,
               GlobalOrdinal offset)
  : stridingInfo_(stridingInfo)
  , stridedBlockId_(stridedBlockId)
  , offset_(offset)
  , indexBase_(indexBase) {
  using MapFactory_t = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>;

  size_t blkSize = getFixedBlockSize();
  TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo.size() == 0,
                             Exceptions::RuntimeError,
                             "StridedMap::StridedMap: stridingInfo not valid: stridingInfo.size() = 0?");
  if (numGlobalElements != Teuchos::OrdinalTraits<global_size_t>::invalid()) {
    TEUCHOS_TEST_FOR_EXCEPTION(numGlobalElements % blkSize != 0,
                               Exceptions::RuntimeError,
                               "StridedMap::StridedMap: stridingInfo not valid: getFixedBlockSize is not an integer "
                               "multiple of numGlobalElements.");
#ifdef HAVE_XPETRA_DEBUG
    // We have to do this check ourselves, as we don't necessarily construct the full Tpetra map
    global_size_t sumLocalElements;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, Teuchos::as<global_size_t>(numLocalElements), Teuchos::outArg(sumLocalElements));

    TEUCHOS_TEST_FOR_EXCEPTION(sumLocalElements != numGlobalElements,
                               std::invalid_argument,
                               "StridedMap::StridedMap: sum of numbers of local elements is different from the provided "
                               "number of global elements.");
#endif
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
      numLocalElements % blkSize != 0,
      Exceptions::RuntimeError,
      "StridedMap::StridedMap: stridingInfo not valid: getFixedBlockSize is not an integer multiple of numLocalElements.");

  if (stridedBlockId != -1) {
    TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo.size() < Teuchos::as<size_t>(stridedBlockId),
                               Exceptions::RuntimeError,
                               "StridedTpetraMap::StridedTpetraMap: stridedBlockId > stridingInfo.size()");
  }

  // Try to create a shortcut
  if (blkSize != 1 || offset_ != 0) {
    // check input data and reorganize map
    global_size_t numGlobalNodes = Teuchos::OrdinalTraits<global_size_t>::invalid();
    if (numGlobalElements != Teuchos::OrdinalTraits<global_size_t>::invalid()) {
      numGlobalNodes = numGlobalElements / blkSize;
    }
    global_size_t numLocalNodes = numLocalElements / blkSize;

    // build an equally distributed node map
    RCP<Map> nodeMap = MapFactory_t::Build(xlib, numGlobalNodes, numLocalNodes, indexBase, comm);

    // translate local node ids to local dofs
    size_t nStridedOffset = 0;
    size_t nDofsPerNode   = blkSize;  // dofs per node for local striding block
    if (stridedBlockId > -1) {
      for (int j = 0; j < stridedBlockId; j++) {
        nStridedOffset += stridingInfo_[j];
      }

      nDofsPerNode      = stridingInfo_[stridedBlockId];
      numGlobalElements = nodeMap->getGlobalNumElements() * Teuchos::as<global_size_t>(nDofsPerNode);
    }
    numLocalElements = numLocalNodes * Teuchos::as<size_t>(nDofsPerNode);

    std::vector<GlobalOrdinal> dofgids(numLocalElements);
    for (LocalOrdinal i = 0; i < Teuchos::as<LocalOrdinal>(numLocalNodes); i++) {
      GlobalOrdinal nodeGID = nodeMap->getGlobalElement(i);

      for (size_t j = 0; j < nDofsPerNode; j++) {
        dofgids[i * nDofsPerNode + j] = indexBase_ + offset_ + (nodeGID - indexBase_) * Teuchos::as<GlobalOrdinal>(blkSize) + Teuchos::as<GlobalOrdinal>(nStridedOffset + j);
      }
    }

    map_ = MapFactory_t::Build(xlib, numGlobalElements, dofgids, indexBase, comm);

    if (stridedBlockId == -1) {
      TEUCHOS_TEST_FOR_EXCEPTION(getLocalNumElements() != Teuchos::as<size_t>(nodeMap->getLocalNumElements() * nDofsPerNode),
                                 Exceptions::RuntimeError,
                                 "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");

      TEUCHOS_TEST_FOR_EXCEPTION(getGlobalNumElements() != Teuchos::as<size_t>(nodeMap->getGlobalNumElements() * nDofsPerNode),
                                 Exceptions::RuntimeError,
                                 "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");
    } else {
      int nDofsInStridedBlock = stridingInfo[stridedBlockId];

      TEUCHOS_TEST_FOR_EXCEPTION(getLocalNumElements() != Teuchos::as<size_t>(nodeMap->getLocalNumElements() * nDofsInStridedBlock),
                                 Exceptions::RuntimeError,
                                 "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");

      TEUCHOS_TEST_FOR_EXCEPTION(getGlobalNumElements() != Teuchos::as<size_t>(nodeMap->getGlobalNumElements() * nDofsInStridedBlock),
                                 Exceptions::RuntimeError,
                                 "StridedTpetraMap::StridedTpetraMap: wrong distribution of dofs among processors.");
    }
  } else {
    map_ = MapFactory_t::Build(xlib, numGlobalElements, numLocalElements, indexBase, comm);
  }

  TEUCHOS_TEST_FOR_EXCEPTION(CheckConsistency() == false, Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: CheckConsistency() == false");
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    StridedMap(UnderlyingLib xlib,
               global_size_t numGlobalElements,
               const Teuchos::ArrayView<const GlobalOrdinal>& elementList,
               GlobalOrdinal indexBase,
               std::vector<size_t>& stridingInfo,
               const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
               LocalOrdinal stridedBlockId)
  : stridingInfo_(stridingInfo)
  , stridedBlockId_(stridedBlockId)
  , indexBase_(indexBase) {
  using MapFactory_t = Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node>;

  size_t blkSize = getFixedBlockSize();

  TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo.size() == 0,
                             Exceptions::RuntimeError,
                             "StridedMap::StridedMap: stridingInfo not valid: stridingInfo.size() = 0?");
  if (stridedBlockId != -1)
    TEUCHOS_TEST_FOR_EXCEPTION(stridingInfo.size() < Teuchos::as<size_t>(stridedBlockId),
                               Exceptions::RuntimeError,
                               "StridedTpetraMap::StridedTpetraMap: stridedBlockId > stridingInfo.size()");
  if (numGlobalElements != Teuchos::OrdinalTraits<global_size_t>::invalid()) {
    TEUCHOS_TEST_FOR_EXCEPTION(numGlobalElements % blkSize != 0,
                               Exceptions::RuntimeError,
                               "StridedMap::StridedMap: stridingInfo not valid: getFixedBlockSize is not an integer "
                               "multiple of numGlobalElements.");
#ifdef HAVE_XPETRA_DEBUG
    // We have to do this check ourselves, as we don't necessarily construct the full Tpetra map
    global_size_t sumLocalElements, numLocalElements = elementList.size();
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, numLocalElements, Teuchos::outArg(sumLocalElements));
    TEUCHOS_TEST_FOR_EXCEPTION(sumLocalElements != numGlobalElements,
                               std::invalid_argument,
                               "StridedMap::StridedMap: sum of numbers of local elements is different from the provided "
                               "number of global elements.");
#endif
  }

  if (stridedBlockId == -1) {
    // numGlobalElements can be -1! FIXME
    // TEUCHOS_TEST_FOR_EXCEPTION(numGlobalElements  % blkSize != 0, Exceptions::RuntimeError,
    // "StridedMap::StridedMap: stridingInfo not valid: getFixedBlockSize is not an integer multiple of
    // numGlobalElements.");
    TEUCHOS_TEST_FOR_EXCEPTION(elementList.size() % blkSize != 0,
                               Exceptions::RuntimeError,
                               "StridedMap::StridedMap: stridingInfo not valid: getFixedBlockSize is not an integer "
                               "multiple of elementList.size().");
  } else {
    // numGlobalElements can be -1! FIXME
    // TEUCHOS_TEST_FOR_EXCEPTION(numGlobalElements  % stridingInfo[stridedBlockId] != 0, Exceptions::RuntimeError,
    // "StridedMap::StridedMap: stridingInfo not valid: stridingBlockInfo[stridedBlockId] is not an integer multiple of
    // numGlobalElements.");
    TEUCHOS_TEST_FOR_EXCEPTION(elementList.size() % stridingInfo[stridedBlockId] != 0,
                               Exceptions::RuntimeError,
                               "StridedMap::StridedMap: stridingInfo not valid: stridingBlockInfo[stridedBlockId] is not "
                               "an integer multiple of elementList.size().");
  }

  map_ = MapFactory_t::Build(xlib, numGlobalElements, elementList, indexBase, comm);

  // calculate offset_

  // find minimum GID over all procs
  GlobalOrdinal minGidOnCurProc = Teuchos::OrdinalTraits<GlobalOrdinal>::max();
  for (Teuchos_Ordinal k = 0; k < elementList.size(); k++)  // TODO fix occurence of Teuchos_Ordinal
  {
    if (elementList[k] < minGidOnCurProc) {
      minGidOnCurProc = elementList[k];
    }
  }

  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, minGidOnCurProc, Teuchos::outArg(offset_));

  // calculate striding index
  size_t nStridedOffset = 0;
  for (int j = 0; j < stridedBlockId; j++) {
    nStridedOffset += stridingInfo[j];
  }
  const GlobalOrdinal goStridedOffset = Teuchos::as<GlobalOrdinal>(nStridedOffset);

  // adapt offset_
  offset_ -= goStridedOffset + indexBase_;

  TEUCHOS_TEST_FOR_EXCEPTION(CheckConsistency() == false, Exceptions::RuntimeError, "StridedTpetraMap::StridedTpetraMap: CheckConsistency() == false");
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    StridedMap(const RCP<const Map>& map,
               std::vector<size_t>& stridingInfo,
               GlobalOrdinal /* indexBase */,
               LocalOrdinal stridedBlockId,
               GlobalOrdinal offset)
  : stridingInfo_(stridingInfo)
  , stridedBlockId_(stridedBlockId)
  , offset_(offset)
  , indexBase_(map->getIndexBase()) {
  // TAW: 11/24/15
  //      A strided map never can be built from a strided map. getMap always returns the underlying
  //      Xpetra::Map object which contains the data (either in a Xpetra::EpetraMapT or Xpetra::TpetraMap
  //      object)
  if (Teuchos::rcp_dynamic_cast<const StridedMap>(map) == Teuchos::null) {
    map_ = map;  // if map is not a strided map, just store it (standard case)
  } else {
    map_ = map->getMap();  // if map is also a strided map, store the underlying plain Epetra/Tpetra Xpetra map object
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    ~StridedMap() {
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
std::vector<size_t>
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getStridingData() const {
  return stridingInfo_;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    setStridingData(std::vector<size_t> stridingInfo) {
  stridingInfo_ = stridingInfo;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getFixedBlockSize() const {
  size_t blkSize = 0;
  for (std::vector<size_t>::const_iterator it = stridingInfo_.begin(); it != stridingInfo_.end(); ++it) {
    blkSize += *it;
  }
  return blkSize;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getStridedBlockId() const {
  return stridedBlockId_;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isStrided() const {
  return stridingInfo_.size() > 1 ? true : false;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isBlocked() const {
  return getFixedBlockSize() > 1 ? true : false;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getOffset() const {
  return offset_;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    setOffset(GlobalOrdinal offset) {
  offset_ = offset;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    GID2StridingBlockId(GlobalOrdinal gid) const {
  GlobalOrdinal tgid = gid - offset_ - indexBase_;
  tgid               = tgid % getFixedBlockSize();

  size_t nStridedOffset = 0;
  size_t stridedBlockId = 0;
  for (size_t j = 0; j < stridingInfo_.size(); j++) {
    nStridedOffset += stridingInfo_[j];
    if (Teuchos::as<size_t>(tgid) < nStridedOffset) {
      stridedBlockId = j;
      break;
    }
  }
  return stridedBlockId;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMap() const {
  return map_;
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    CheckConsistency() {
#ifndef HAVE_XPETRA_DEBUG
  return true;
#else
  if (getStridedBlockId() == -1) {
    // Strided map contains the full map
    if (getLocalNumElements() % getFixedBlockSize() != 0 ||  // number of local  elements is not a multiple of block size
        getGlobalNumElements() % getFixedBlockSize() != 0)   // number of global    -//-
      return false;
  } else {
    // Strided map contains only the partial map
    Teuchos::ArrayView<const GlobalOrdinal> dofGids = getLocalElementList();
    // std::sort(dofGids.begin(), dofGids.end());

    if (dofGids.size() == 0)  // special treatment for empty processors
    {
      return true;
    }

    if (dofGids.size() % stridingInfo_[stridedBlockId_] != 0) {
      return false;
    }

    // Calculate nStridedOffset
    size_t nStridedOffset = 0;
    for (int j = 0; j < stridedBlockId_; j++) {
      nStridedOffset += stridingInfo_[j];
    }

    const GlobalOrdinal goStridedOffset = Teuchos::as<GlobalOrdinal>(nStridedOffset);
    const GlobalOrdinal goZeroOffset    = (dofGids[0] - nStridedOffset - offset_ - indexBase_) / Teuchos::as<GlobalOrdinal>(getFixedBlockSize());

    GlobalOrdinal cnt = 0;
    for (size_t i = 0;
         i < Teuchos::as<size_t>(dofGids.size()) / stridingInfo_[stridedBlockId_];
         i += stridingInfo_[stridedBlockId_]) {
      const GlobalOrdinal first_gid = dofGids[i];

      // We expect this to be the same for all DOFs of the same node
      cnt = (first_gid - goStridedOffset - offset_ - indexBase_) / Teuchos::as<GlobalOrdinal>(getFixedBlockSize()) - goZeroOffset;

      // Loop over all DOFs that belong to current node
      for (size_t j = 0; j < stridingInfo_[stridedBlockId_]; j++) {
        const GlobalOrdinal gid = dofGids[i + j];
        const GlobalOrdinal r   = (gid - Teuchos::as<GlobalOrdinal>(j) - goStridedOffset - offset_ - indexBase_) / Teuchos::as<GlobalOrdinal>(getFixedBlockSize()) - goZeroOffset - cnt;
        // TAW 1/18/2016: We cannot use Teuchos::OrdinalTraits<GlobalOrdinal>::zero() ) here,
        //                If, e.g., GO=long long is disabled, OrdinalTraits<long long> is not available.
        //                But we instantiate stubs on GO=long long which might contain StridedMaps.
        //                These lead to compilation errors, then.
        if (0 != r) {
          std::cout << "goZeroOffset   : " << goZeroOffset << std::endl
                    << "dofGids[0]     : " << dofGids[0] << std::endl
                    << "stridedOffset  : " << nStridedOffset << std::endl
                    << "offset_        : " << offset_ << std::endl
                    << "goStridedOffset: " << goStridedOffset << std::endl
                    << "getFixedBlkSize: " << getFixedBlockSize() << std::endl
                    << "gid: " << gid << " GID: " << r << std::endl;

          return false;
        }
      }
    }
  }

  return true;
#endif
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalNumElements() const {
  return map_->getGlobalNumElements();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getLocalNumElements() const {
  return map_->getLocalNumElements();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getIndexBase() const {
  return map_->getIndexBase();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMinLocalIndex() const {
  return map_->getMinLocalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMaxLocalIndex() const {
  return map_->getMaxLocalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMinGlobalIndex() const {
  return map_->getMinGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMaxGlobalIndex() const {
  return map_->getMaxGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMinAllGlobalIndex() const {
  return map_->getMinAllGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMaxAllGlobalIndex() const {
  return map_->getMaxAllGlobalIndex();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getLocalElement(GlobalOrdinal globalIndex) const {
  return map_->getLocalElement(globalIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
GlobalOrdinal
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalElement(LocalOrdinal localIndex) const {
  return map_->getGlobalElement(localIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LookupStatus
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal>& GIDList,
                       const Teuchos::ArrayView<int>& nodeIDList,
                       const Teuchos::ArrayView<LocalOrdinal>& LIDList) const {
  return map_->getRemoteIndexList(GIDList, nodeIDList, LIDList);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
LookupStatus
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal>& GIDList,
                       const Teuchos::ArrayView<int>& nodeIDList) const {
  return map_->getRemoteIndexList(GIDList, nodeIDList);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayView<const GlobalOrdinal>
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getLocalElementList() const {
  return map_->getLocalElementList();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
typename Map<LocalOrdinal, GlobalOrdinal, Node>::global_indices_array_device_type
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getMyGlobalIndicesDevice() const {
  return map_->getMyGlobalIndicesDevice();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isNodeLocalElement(LocalOrdinal localIndex) const {
  return map_->isNodeLocalElement(localIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isNodeGlobalElement(GlobalOrdinal globalIndex) const {
  return map_->isNodeGlobalElement(globalIndex);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isContiguous() const {
  return map_->isContiguous();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isDistributed() const {
  return map_->isDistributed();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isCompatible(const Map& map) const {
  return map_->isCompatible(map);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
bool StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    isSameAs(const Map& map) const {
  return map_->isSameAs(map);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Teuchos::Comm<int>>
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    getComm() const {
  return map_->getComm();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    removeEmptyProcesses() const {
  return map_->removeEmptyProcesses();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    replaceCommWithSubset(const Teuchos::RCP<const Teuchos::Comm<int>>& newComm) const {
  return map_->replaceCommWithSubset(newComm);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
std::string
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    description() const {
  return map_->description();
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
void StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
  map_->describe(out, verbLevel);
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
UnderlyingLib
StridedMap<LocalOrdinal, GlobalOrdinal, Node>::
    lib() const {
  return map_->lib();
}

}  // namespace Xpetra

#endif  // XPETRA_STRIDEDMAP_DEF_HPP
