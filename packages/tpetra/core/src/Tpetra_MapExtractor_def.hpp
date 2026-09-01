// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_MAPEXTRACTOR_DEF_HPP
#define TPETRA_MAPEXTRACTOR_DEF_HPP

#include "Tpetra_MapExtractor_decl.hpp"  // completes MapExtractor before cross-includes

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_BlockedMap_def.hpp"
#include "Tpetra_BlockedMultiVector_def.hpp"

namespace Tpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    MapExtractor(const Teuchos::RCP<const Map>& fullmap, const std::vector<Teuchos::RCP<const Map>>& maps, bool bThyraMode) {
  map_ = Teuchos::rcp(new BlockedMap(fullmap, maps, bThyraMode));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    MapExtractor(const std::vector<Teuchos::RCP<const Map>>& maps, const std::vector<Teuchos::RCP<const Map>>& thyramaps) {
  map_ = Teuchos::rcp(new BlockedMap(maps, thyramaps));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    MapExtractor(const Teuchos::RCP<const BlockedMap>& blockedMap)
  : map_(blockedMap) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    MapExtractor(const MapExtractor& input) {
  map_ = Teuchos::rcp(new BlockedMap(*(input.getBlockedMap())));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ~MapExtractor() {
  map_ = Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(const Vector& full, size_t block, Vector& partial) const {
  TEUCHOS_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range,
                             "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false) == Teuchos::null, std::runtime_error,
                             "ExtractVector: map_->getMap(" << block << ",false) is null");
  partial.doImport(full, *(map_->getImporter(block)), Tpetra::INSERT);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(const MultiVector& full, size_t block, MultiVector& partial) const {
  TEUCHOS_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range,
                             "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false) == Teuchos::null, std::runtime_error,
                             "ExtractVector: map_->getMap(" << block << ",false) is null");
  partial.doImport(full, *(map_->getImporter(block)), Tpetra::INSERT);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(Teuchos::RCP<const Vector>& full, size_t block, Teuchos::RCP<Vector>& partial) const {
  ExtractVector(*full, block, *partial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(Teuchos::RCP<Vector>& full, size_t block, Teuchos::RCP<Vector>& partial) const {
  ExtractVector(*full, block, *partial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(Teuchos::RCP<const MultiVector>& full, size_t block, Teuchos::RCP<MultiVector>& partial) const {
  ExtractVector(*full, block, *partial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(Teuchos::RCP<MultiVector>& full, size_t block, Teuchos::RCP<MultiVector>& partial) const {
  ExtractVector(*full, block, *partial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Vector>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(Teuchos::RCP<const Vector>& full, size_t block, bool bThyraMode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range,
                             "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false) == Teuchos::null, std::runtime_error,
                             "ExtractVector: map_->getMap(" << block << ",false) is null");
  // first extract partial vector from full vector (using xpetra style GIDs)
  const Teuchos::RCP<Vector> vv = Teuchos::rcp(new Vector(getMap(block, false), false));
  ExtractVector(*full, block, *vv);
  if (bThyraMode == false)
    return vv;
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, std::runtime_error,
                             "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been "
                             "created using Thyra-style numbered submaps.");
  vv->replaceMap(getMap(block, true));  // switch to Thyra-style map
  return vv;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Vector>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(Teuchos::RCP<Vector>& full, size_t block, bool bThyraMode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range,
                             "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false) == Teuchos::null, std::runtime_error,
                             "ExtractVector: map_->getMap(" << block << ",false) is null");
  const Teuchos::RCP<Vector> vv = Teuchos::rcp(new Vector(getMap(block, false), false));
  ExtractVector(*full, block, *vv);
  if (bThyraMode == false)
    return vv;
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, std::runtime_error,
                             "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been "
                             "created using Thyra-style numbered submaps.");
  vv->replaceMap(getMap(block, true));  // switch to Thyra-style map
  return vv;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiVector>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(Teuchos::RCP<const MultiVector>& full, size_t block, bool bThyraMode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range,
                             "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false) == Teuchos::null, std::runtime_error,
                             "ExtractVector: map_->getMap(" << block << ",false) is null");
  Teuchos::RCP<const BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(full);
  if (bfull.is_null() == true) {
    // standard case: full is not of type BlockedMultiVector
    const Teuchos::RCP<MultiVector> vv = Teuchos::rcp(new MultiVector(getMap(block, false), full->getNumVectors(), false));
    Teuchos::RCP<const Map> oldThyMapFull     = full->getMap();  // temporarily store map of full
    Teuchos::RCP<MultiVector> rcpNonConstFull = Teuchos::rcp_const_cast<MultiVector>(full);
    rcpNonConstFull->replaceMap(map_->getImporter(block)->getSourceMap());
    ExtractVector(*rcpNonConstFull, block, *vv);
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, std::runtime_error,
                               "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been "
                               "created using Thyra-style numbered submaps.");
    if (bThyraMode == true)
      vv->replaceMap(getMap(block, true));  // switch to Thyra-style map
    rcpNonConstFull->replaceMap(oldThyMapFull);
    // If sub-block is itself blocked, re-slice the extracted vector into a
    // nested BlockedMultiVector so nested identity is preserved (mirrors
    // Xpetra's MultiVectorFactory::Build(BlockedMap)).
    Teuchos::RCP<const BlockedMap> nested = map_->getBlockedSubMap(block);
    if (!nested.is_null())
      return Teuchos::rcp(new BlockedMultiVector(nested, Teuchos::rcp_implicit_cast<const MultiVector>(vv)));
    return vv;
  } else {
    // special case: full is of type BlockedMultiVector
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getNumMaps() != bfull->getBlockedMap()->getNumMaps(), std::runtime_error,
                               "ExtractVector: Number of blocks in map extractor is " << map_->getNumMaps() << " but should be "
                                                                                      << bfull->getBlockedMap()->getNumMaps()
                                                                                      << " (number of blocks in BlockedMultiVector)");
    return bfull->getMultiVector(block, bThyraMode);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiVector>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(Teuchos::RCP<MultiVector>& full, size_t block, bool bThyraMode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range,
                             "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false) == Teuchos::null, std::runtime_error,
                             "ExtractVector: map_->getMap(" << block << ",false) is null");
  Teuchos::RCP<BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(full);
  if (bfull.is_null() == true) {
    const Teuchos::RCP<MultiVector> vv = Teuchos::rcp(new MultiVector(getMap(block, false), full->getNumVectors(), false));
    Teuchos::RCP<const Map> oldThyMapFull = full->getMap();  // temporarily store map of full
    full->replaceMap(map_->getImporter(block)->getSourceMap());
    ExtractVector(*full, block, *vv);
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, std::runtime_error,
                               "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been "
                               "created using Thyra-style numbered submaps.");
    if (bThyraMode == true)
      vv->replaceMap(getMap(block, true));  // switch to Thyra-style map
    full->replaceMap(oldThyMapFull);
    // If sub-block is itself blocked, re-slice into a nested BlockedMultiVector.
    Teuchos::RCP<const BlockedMap> nested = map_->getBlockedSubMap(block);
    if (!nested.is_null())
      return Teuchos::rcp(new BlockedMultiVector(nested, Teuchos::rcp_implicit_cast<const MultiVector>(vv)));
    return vv;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getNumMaps() != bfull->getBlockedMap()->getNumMaps(), std::runtime_error,
                               "ExtractVector: Number of blocks in map extractor is " << map_->getNumMaps() << " but should be "
                                                                                      << bfull->getBlockedMap()->getNumMaps()
                                                                                      << " (number of blocks in BlockedMultiVector)");
    return bfull->getMultiVector(block, bThyraMode);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiVector>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(Teuchos::RCP<const BlockedMultiVector>& full, size_t block, bool bThyraMode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range,
                             "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false) == Teuchos::null, std::runtime_error,
                             "ExtractVector: map_->getMap(" << block << ",false) is null");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getNumMaps() != full->getBlockedMap()->getNumMaps(), std::runtime_error,
                             "ExtractVector: Number of blocks in map extractor is " << map_->getNumMaps() << " but should be "
                                                                                    << full->getBlockedMap()->getNumMaps()
                                                                                    << " (number of blocks in BlockedMultiVector)");
  Teuchos::RCP<MultiVector> vv = full->getMultiVector(block, bThyraMode);
  return vv;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiVector>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(Teuchos::RCP<BlockedMultiVector>& full, size_t block, bool bThyraMode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range,
                             "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false) == Teuchos::null, std::runtime_error,
                             "ExtractVector: map_->getMap(" << block << ",false) is null");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getNumMaps() != full->getBlockedMap()->getNumMaps(), std::runtime_error,
                             "ExtractVector: Number of blocks in map extractor is " << map_->getNumMaps() << " but should be "
                                                                                    << full->getBlockedMap()->getNumMaps()
                                                                                    << " (number of blocks in BlockedMultiVector)");
  Teuchos::RCP<MultiVector> vv = full->getMultiVector(block, bThyraMode);
  return vv;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(const Vector& partial, size_t block, Vector& full, bool bThyraMode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range,
                             "InsertVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false) == Teuchos::null, std::runtime_error,
                             "InsertVector: map_->getMap(" << block << ",false) is null");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, std::runtime_error,
                             "MapExtractor::InsertVector: InsertVector in Thyra-style numbering only possible if MapExtractor has been created "
                             "using Thyra-style numbered submaps.");
  if (bThyraMode) {
    Teuchos::RCP<const MultiVector> rcpPartial   = Teuchos::rcpFromRef(partial);
    Teuchos::RCP<MultiVector> rcpNonConstPartial = Teuchos::rcp_const_cast<MultiVector>(rcpPartial);
    Teuchos::RCP<const Map> oldThyMapPartial     = rcpNonConstPartial->getMap();  // temporarily store map of partial
    Teuchos::RCP<const Map> oldThyMapFull        = full.getMap();                 // temporarily store map of full

    TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false)->isSameAs(*(map_->getImporter(block)->getTargetMap())) == false, std::runtime_error,
                               "MapExtractor::InsertVector: InsertVector in Thyra-style mode: Xpetra GIDs of partial vector are not identical "
                               "to target Map of Importer. This should not be.");

    rcpNonConstPartial->replaceMap(getMap(block, false));       // temporarily switch to xpetra-style map
    full.replaceMap(map_->getImporter(block)->getSourceMap());  // temporarily switch to Xpetra GIDs

    full.doExport(*rcpNonConstPartial, *(map_->getImporter(block)), Tpetra::INSERT);

    full.replaceMap(oldThyMapFull);                    // reset original map (Thyra GIDs)
    rcpNonConstPartial->replaceMap(oldThyMapPartial);  // change map back to original map
  } else {
    full.doExport(partial, *(map_->getImporter(block)), Tpetra::INSERT);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(const MultiVector& partial, size_t block, MultiVector& full, bool bThyraMode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range,
                             "InsertVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false) == Teuchos::null, std::runtime_error,
                             "InsertVector: map_->getMap(" << block << ",false) is null");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, std::runtime_error,
                             "MapExtractor::InsertVector: InsertVector in Thyra-style numbering only possible if MapExtractor has been created "
                             "using Thyra-style numbered submaps.");
  if (bThyraMode) {
    Teuchos::RCP<const MultiVector> rcpPartial   = Teuchos::rcpFromRef(partial);
    Teuchos::RCP<MultiVector> rcpNonConstPartial = Teuchos::rcp_const_cast<MultiVector>(rcpPartial);
    Teuchos::RCP<const Map> oldThyMapPartial     = rcpNonConstPartial->getMap();  // temporarily store map of partial
    Teuchos::RCP<const Map> oldThyMapFull        = full.getMap();                 // temporarily store map of full

    TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false)->isSameAs(*(map_->getImporter(block)->getTargetMap())) == false, std::runtime_error,
                               "MapExtractor::InsertVector: InsertVector in Thyra-style mode: Xpetra GIDs of partial vector are not identical "
                               "to target Map of Importer. This should not be.");

    rcpNonConstPartial->replaceMap(getMap(block, false));       // temporarily switch to xpetra-style map
    full.replaceMap(map_->getImporter(block)->getSourceMap());  // temporarily switch to Xpetra GIDs

    full.doExport(*rcpNonConstPartial, *(map_->getImporter(block)), Tpetra::INSERT);

    full.replaceMap(oldThyMapFull);                    // reset original map (Thyra GIDs)
    rcpNonConstPartial->replaceMap(oldThyMapPartial);  // change map back to original map
  } else {
    full.doExport(partial, *(map_->getImporter(block)), Tpetra::INSERT);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(Teuchos::RCP<const Vector> partial, size_t block, Teuchos::RCP<Vector> full, bool bThyraMode) const {
  InsertVector(*partial, block, *full, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(Teuchos::RCP<Vector> partial, size_t block, Teuchos::RCP<Vector> full, bool bThyraMode) const {
  InsertVector(*partial, block, *full, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(Teuchos::RCP<const MultiVector> partial, size_t block, Teuchos::RCP<MultiVector> full, bool bThyraMode) const {
  Teuchos::RCP<BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(full);
  if (bfull.is_null() == true)
    InsertVector(*partial, block, *full, bThyraMode);
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false) == Teuchos::null, std::runtime_error,
                               "InsertVector: map_->getMap(" << block << ",false) is null");
    // As in Xpetra: a const partial cannot be handed to BlockedMultiVector::setMultiVector here.
    throw std::runtime_error("Tpetra::MapExtractor::InsertVector: cannot insert a const MultiVector into a BlockedMultiVector.");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(Teuchos::RCP<MultiVector> partial, size_t block, Teuchos::RCP<MultiVector> full, bool bThyraMode) const {
  Teuchos::RCP<BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(full);
  if (bfull.is_null() == true)
    InsertVector(*partial, block, *full, bThyraMode);
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false) == Teuchos::null, std::runtime_error,
                               "InsertVector: map_->getMap(" << block << ",false) is null");
    bfull->setMultiVector(block, partial, bThyraMode);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(Teuchos::RCP<const MultiVector> partial, size_t block, Teuchos::RCP<BlockedMultiVector> full, bool bThyraMode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false) == Teuchos::null, std::runtime_error,
                             "InsertVector: map_->getMap(" << block << ",false) is null");
  full->setMultiVector(block, partial, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(Teuchos::RCP<MultiVector> partial, size_t block, Teuchos::RCP<BlockedMultiVector> full, bool bThyraMode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block, false) == Teuchos::null, std::runtime_error,
                             "InsertVector: map_->getMap(" << block << ",false) is null");
  full->setMultiVector(block, partial, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Vector>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVector(size_t i, bool bThyraMode, bool bZero) const {
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, std::runtime_error,
                             "MapExtractor::getVector: getVector in Thyra-style numbering only possible if MapExtractor has been created using "
                             "Thyra-style numbered submaps.");
  return Teuchos::rcp(new Vector(getMap(i, bThyraMode), bZero));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiVector>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVector(size_t i, size_t numvec, bool bThyraMode, bool bZero) const {
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true, std::runtime_error,
                             "MapExtractor::getVector: getVector in Thyra-style numbering only possible if MapExtractor has been created using "
                             "Thyra-style numbered submaps.");
  // If sub-block i is itself blocked, build a nested BlockedMultiVector so the
  // nested blocked operators (e.g. a nested BGS smoother) see a properly nested
  // vector rather than a flattened plain one.
  Teuchos::RCP<const BlockedMap> nested = map_->getBlockedSubMap(i);
  if (!nested.is_null())
    return Teuchos::rcp(new BlockedMultiVector(nested, numvec, bZero));
  return Teuchos::rcp(new MultiVector(getMap(i, bThyraMode), numvec, bZero));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getThyraMode() const {
  return map_->getThyraMode();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    NumMaps() const {
  return map_->getNumMaps();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const typename MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMap(size_t i, bool bThyraMode) const {
  return map_->getMap(i, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const typename MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMap() const {
  // BlockedMap is not a Tpetra::Map; return the full map.
  return map_->getFullMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const typename MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockedMap>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getBlockedMap() const {
  return map_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const typename MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getFullMap() const {
  return map_->getFullMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMapIndexForGID(GlobalOrdinal gid) const {
  return map_->getMapIndexForGID(gid);
}

}  // namespace Tpetra

//
// Explicit instantiation macro
//
#define TPETRA_MAPEXTRACTOR_INSTANT(SCALAR, LO, GO, NODE) \
  template class MapExtractor<SCALAR, LO, GO, NODE>;

#endif  // TPETRA_MAPEXTRACTOR_DEF_HPP
