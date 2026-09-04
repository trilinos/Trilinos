// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_MAPEXTRACTOR_DEF_HPP_
#define XPETRA_MAPEXTRACTOR_DEF_HPP_

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_BlockedMultiVector.hpp>
#include <Xpetra_BlockedMap.hpp>
#include <Xpetra_TpetraMap.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraVector.hpp>

#include <Xpetra_MapExtractor_decl.hpp>

namespace Xpetra {

// Local helpers to bridge the Xpetra <-> Tpetra type gap.  These reuse the blocked-aware
// unwrap/wrap helpers from BlockedMultiVectorDetails so that blocked Xpetra vectors unwrap
// to Tpetra::BlockedMultiVector (and are re-wrapped preserving the blocked dynamic type),
// while plain vectors go through the ordinary TpetraMultiVector/TpetraVector conversions.
// Epetra-backed objects trigger Xpetra::Exceptions::BadCast.
namespace MapExtractorDetails {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::vector<Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>>
toTpetraMaps(const std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>>& maps) {
  std::vector<Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>> tmaps(maps.size());
  for (size_t i = 0; i < maps.size(); ++i)
    tmaps[i] = (maps[i] == Teuchos::null) ? Teuchos::null : toTpetra(maps[i]);
  return tmaps;
}

}  // namespace MapExtractorDetails

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    MapExtractor(const RCP<const Map>& fullmap, const std::vector<RCP<const Map>>& maps, bool bThyraMode) {
  // Build an Xpetra::BlockedMap wrapper first: it caches the original Xpetra
  // sub-maps and full map so their dynamic type (e.g. StridedMap) survives
  // getMap()/getFullMap().  The Tpetra core extractor is then built from that
  // wrapper's Tpetra blocked map -- identical to building it directly from the
  // unwrapped maps, but preserving the richer Xpetra identity for MueLu.
  blockedMapXpetra_ = Teuchos::rcp(new Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>(fullmap, maps, bThyraMode));
  te_               = Teuchos::rcp(new tpetra_mapextractor_type(blockedMapXpetra_->getTpetra_BlockedMap()));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    MapExtractor(const std::vector<RCP<const Map>>& maps, const std::vector<RCP<const Map>>& thyramaps) {
  // As above: route through an Xpetra::BlockedMap so the original sub-maps'
  // dynamic type is retained for getMap().
  blockedMapXpetra_ = Teuchos::rcp(new Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>(maps, thyramaps));
  te_               = Teuchos::rcp(new tpetra_mapextractor_type(blockedMapXpetra_->getTpetra_BlockedMap()));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    MapExtractor(const Teuchos::RCP<const BlockedMap>& blockedMap) {
  te_               = Teuchos::rcp(new tpetra_mapextractor_type(blockedMap->getTpetra_BlockedMap()));
  blockedMapXpetra_ = blockedMap;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    MapExtractor(const MapExtractor& input) {
  te_               = Teuchos::rcp(new tpetra_mapextractor_type(*input.getTpetra_MapExtractor()));
  blockedMapXpetra_ = input.blockedMapXpetra_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ~MapExtractor() {
  te_               = Teuchos::null;
  blockedMapXpetra_ = Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(const Vector& full, size_t block, Vector& partial) const {
  Teuchos::RCP<const MultiVector> fullMV    = Teuchos::rcpFromRef(full);
  Teuchos::RCP<MultiVector> partialMV       = Teuchos::rcpFromRef(partial);
  Teuchos::RCP<const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tfull =
      Teuchos::rcp_dynamic_cast<const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(
          BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fullMV), true);
  Teuchos::RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpartial =
      Teuchos::rcp_dynamic_cast<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(
          BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(partialMV), true);
  te_->ExtractVector(*tfull, block, *tpartial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(const MultiVector& full, size_t block, MultiVector& partial) const {
  Teuchos::RCP<const MultiVector> fullMV    = Teuchos::rcpFromRef(full);
  Teuchos::RCP<MultiVector> partialMV       = Teuchos::rcpFromRef(partial);
  Teuchos::RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tfull =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fullMV);
  Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpartial =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(partialMV);
  te_->ExtractVector(*tfull, block, *tpartial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<const Vector>& full, size_t block, RCP<Vector>& partial) const {
  ExtractVector(*full, block, *partial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<Vector>& full, size_t block, RCP<Vector>& partial) const {
  ExtractVector(*full, block, *partial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<const MultiVector>& full, size_t block, RCP<MultiVector>& partial) const {
  ExtractVector(*full, block, *partial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<MultiVector>& full, size_t block, RCP<MultiVector>& partial) const {
  ExtractVector(*full, block, *partial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& full, size_t block, bool bThyraMode) const {
  Teuchos::RCP<const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tfull = toTpetra(full);
  Teuchos::RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpartial     = te_->ExtractVector(tfull, block, bThyraMode);
  return toXpetra(tpartial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& full, size_t block, bool bThyraMode) const {
  Teuchos::RCP<const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tfull =
      toTpetra(Teuchos::rcp_const_cast<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(full));
  Teuchos::RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpartial = te_->ExtractVector(tfull, block, bThyraMode);
  return toXpetra(tpartial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& full, size_t block, bool bThyraMode) const {
  Teuchos::RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tfull =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(full);
  Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpartial = te_->ExtractVector(tfull, block, bThyraMode);
  return BlockedMultiVectorDetails::wrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tpartial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& full, size_t block, bool bThyraMode) const {
  Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tfull =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(full);
  Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpartial = te_->ExtractVector(tfull, block, bThyraMode);
  return BlockedMultiVectorDetails::wrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tpartial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<const Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& full, size_t block, bool bThyraMode) const {
  Teuchos::RCP<const Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tfull = full->getTpetra_BlockedMultiVector();
  Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpartial           = te_->ExtractVector(tfull, block, bThyraMode);
  return BlockedMultiVectorDetails::wrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tpartial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& full, size_t block, bool bThyraMode) const {
  Teuchos::RCP<Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tfull = full->getTpetra_BlockedMultiVector();
  Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpartial     = te_->ExtractVector(tfull, block, bThyraMode);
  return BlockedMultiVectorDetails::wrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(tpartial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& partial, size_t block, Vector& full, bool bThyraMode) const {
  Teuchos::RCP<const MultiVector> partialMV = Teuchos::rcpFromRef(partial);
  Teuchos::RCP<MultiVector> fullMV          = Teuchos::rcpFromRef(full);
  Teuchos::RCP<const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpartial =
      Teuchos::rcp_dynamic_cast<const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(
          BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(partialMV), true);
  Teuchos::RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tfull =
      Teuchos::rcp_dynamic_cast<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(
          BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fullMV), true);
  te_->InsertVector(*tpartial, block, *tfull, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& partial, size_t block, MultiVector& full, bool bThyraMode) const {
  Teuchos::RCP<const MultiVector> partialMV = Teuchos::rcpFromRef(partial);
  Teuchos::RCP<MultiVector> fullMV          = Teuchos::rcpFromRef(full);
  Teuchos::RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpartial =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(partialMV);
  Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tfull =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fullMV);
  te_->InsertVector(*tpartial, block, *tfull, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(RCP<const Vector> partial, size_t block, RCP<Vector> full, bool bThyraMode) const {
  InsertVector(*partial, block, *full, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(RCP<Vector> partial, size_t block, RCP<Vector> full, bool bThyraMode) const {
  InsertVector(*partial, block, *full, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> partial,
                 size_t block,
                 RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> full,
                 bool bThyraMode) const {
  Teuchos::RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpartial =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(partial);
  Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tfull =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(full);
  te_->InsertVector(tpartial, block, tfull, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> partial,
                 size_t block,
                 RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> full,
                 bool bThyraMode) const {
  Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpartial =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(partial);
  Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tfull =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(full);
  te_->InsertVector(tpartial, block, tfull, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> partial,
                 size_t block,
                 RCP<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> full,
                 bool bThyraMode) const {
  Teuchos::RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpartial =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(partial);
  Teuchos::RCP<Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tfull = full->getTpetra_BlockedMultiVector();
  te_->InsertVector(tpartial, block, tfull, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> partial,
                 size_t block,
                 RCP<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> full,
                 bool bThyraMode) const {
  Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tpartial =
      BlockedMultiVectorDetails::unwrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(partial);
  Teuchos::RCP<Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tfull = full->getTpetra_BlockedMultiVector();
  te_->InsertVector(tpartial, block, tfull, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVector(size_t i, bool bThyraMode, bool bZero) const {
  return toXpetra(te_->getVector(i, bThyraMode, bZero));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVector(size_t i, size_t numvec, bool bThyraMode, bool bZero) const {
  return BlockedMultiVectorDetails::wrapMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(
      te_->getVector(i, numvec, bThyraMode, bZero));
}

/// returns true, if sub maps are stored in Thyra-style numbering
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getThyraMode() const {
  return te_->getThyraMode();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    NumMaps() const {
  return te_->NumMaps();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMap(size_t i, bool bThyraMode) const {
  // Prefer the cached Xpetra sub-map so its dynamic type (e.g. StridedMap) is
  // preserved; the Tpetra core stores only flattened plain maps.
  if (!blockedMapXpetra_.is_null())
    return blockedMapXpetra_->getMap(i, bThyraMode);
  return toXpetra(te_->getMap(i, bThyraMode));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMap() const {
  return getBlockedMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getBlockedMap() const {
  if (blockedMapXpetra_.is_null())
    blockedMapXpetra_ = Teuchos::rcp(new Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>(te_->getBlockedMap()));
  return blockedMapXpetra_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getFullMap() const {
  // Prefer the cached Xpetra full map so its dynamic type (e.g. StridedMap) is
  // preserved for MueLu factories that rcp_dynamic_cast it.
  if (!blockedMapXpetra_.is_null()) {
    RCP<const Map> full = blockedMapXpetra_->getFullMap();
    if (!full.is_null())
      return full;
  }
  return toXpetra(te_->getFullMap());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMapIndexForGID(GlobalOrdinal gid) const {
  return te_->getMapIndexForGID(gid);
}

}  // namespace Xpetra

#endif /* XPETRA_MAPEXTRACTOR_DEF_HPP_ */
