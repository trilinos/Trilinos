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

#include <Xpetra_MapExtractor_decl.hpp>

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    MapExtractor(const RCP<const Map>& fullmap, const std::vector<RCP<const Map>>& maps, bool bThyraMode) {
  map_ = Teuchos::rcp(new BlockedMap(fullmap, maps, bThyraMode));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    MapExtractor(const std::vector<RCP<const Map>>& maps, const std::vector<RCP<const Map>>& thyramaps) {
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
  XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(), std::out_of_range,
                            "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps() << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError,
                            "ExtractVector: map_->getMap(" << block << ",false) is null");

  partial.doImport(full, *(map_->getImporter(block)), Xpetra::INSERT);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(const MultiVector& full, size_t block, MultiVector& partial) const {
  XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                            std::out_of_range,
                            "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps()
                                                             << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "ExtractVector: map_->getMap(" << block << ",false) is null");

  partial.doImport(full, *(map_->getImporter(block)), Xpetra::INSERT);
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
  XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                            std::out_of_range,
                            "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps()
                                                             << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "ExtractVector: map_->getMap(" << block << ",false) is null");
  // first extract partial vector from full vector (using xpetra style GIDs)
  const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> vv = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(getMap(block, false), false);
  ExtractVector(*full, block, *vv);
  if (bThyraMode == false)
    return vv;
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true,
                             Xpetra::Exceptions::RuntimeError,
                             "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been "
                             "created using Thyra-style numbered submaps.");
  vv->replaceMap(getMap(block, true));  // switch to Thyra-style map
  return vv;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& full, size_t block, bool bThyraMode) const {
  XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                            std::out_of_range,
                            "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps()
                                                             << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "ExtractVector: map_->getmap(" << block << ",false) is null");
  // first extract partial vector from full vector (using xpetra style GIDs)
  const RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> vv = Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(getMap(block, false), false);

  ExtractVector(*full, block, *vv);
  if (bThyraMode == false)
    return vv;
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true,
                             Xpetra::Exceptions::RuntimeError,
                             "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been "
                             "created using Thyra-style numbered submaps.");
  vv->replaceMap(getMap(block, true));  // switch to Thyra-style map
  return vv;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& full, size_t block, bool bThyraMode) const {
  XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                            std::out_of_range,
                            "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps()
                                                             << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "ExtractVector: map_->getmap(" << block << ",false) is null");
  RCP<const BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(full);
  if (bfull.is_null() == true) {
    // standard case: full is not of type BlockedMultiVector
    // first extract partial vector from full vector (using xpetra style GIDs)
    const RCP<MultiVector> vv = MultiVectorFactory::Build(getMap(block, false), full->getNumVectors(), false);
    // if(bThyraMode == false) {
    //  ExtractVector(*full, block, *vv);
    //  return vv;
    //} else {
    RCP<const Map> oldThyMapFull     = full->getMap();  // temporarely store map of full
    RCP<MultiVector> rcpNonConstFull = Teuchos::rcp_const_cast<MultiVector>(full);
    rcpNonConstFull->replaceMap(map_->getImporter(block)->getSourceMap());
    ExtractVector(*rcpNonConstFull, block, *vv);
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true,
                               Xpetra::Exceptions::RuntimeError,
                               "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been "
                               "created using Thyra-style numbered submaps.");
    if (bThyraMode == true)
      vv->replaceMap(getMap(block, true));  // switch to Thyra-style map
    rcpNonConstFull->replaceMap(oldThyMapFull);
    return vv;
    //}
  } else {
    // special case: full is of type BlockedMultiVector
    XPETRA_TEST_FOR_EXCEPTION(map_->getNumMaps() != bfull->getBlockedMap()->getNumMaps(),
                              Xpetra::Exceptions::RuntimeError,
                              "ExtractVector: Number of blocks in map extractor is " << map_->getNumMaps() << " but should be "
                                                                                     << bfull->getBlockedMap()->getNumMaps()
                                                                                     << " (number of blocks in BlockedMultiVector)");
    return bfull->getMultiVector(block, bThyraMode);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& full, size_t block, bool bThyraMode) const {
  XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                            std::out_of_range,
                            "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps()
                                                             << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "ExtractVector: map_->getmap(" << block << ",false) is null");
  RCP<BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(full);
  if (bfull.is_null() == true) {
    // standard case: full is not of type BlockedMultiVector
    // first extract partial vector from full vector (using xpetra style GIDs)
    const RCP<MultiVector> vv = MultiVectorFactory::Build(getMap(block, false), full->getNumVectors(), false);
    // if(bThyraMode == false) {
    //  ExtractVector(*full, block, *vv);
    //  return vv;
    //} else {
    RCP<const Map> oldThyMapFull = full->getMap();  // temporarely store map of full
    full->replaceMap(map_->getImporter(block)->getSourceMap());
    ExtractVector(*full, block, *vv);
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true,
                               Xpetra::Exceptions::RuntimeError,
                               "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been "
                               "created using Thyra-style numbered submaps.");
    if (bThyraMode == true)
      vv->replaceMap(getMap(block, true));  // switch to Thyra-style map
    full->replaceMap(oldThyMapFull);
    return vv;
    //}
  } else {
    // special case: full is of type BlockedMultiVector
    XPETRA_TEST_FOR_EXCEPTION(map_->getNumMaps() != bfull->getBlockedMap()->getNumMaps(),
                              Xpetra::Exceptions::RuntimeError,
                              "ExtractVector: Number of blocks in map extractor is " << map_->getNumMaps() << " but should be "
                                                                                     << bfull->getBlockedMap()->getNumMaps()
                                                                                     << " (number of blocks in BlockedMultiVector)");
    return bfull->getMultiVector(block, bThyraMode);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<const Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& full, size_t block, bool bThyraMode) const {
  XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                            std::out_of_range,
                            "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps()
                                                             << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "ExtractVector: map_->getmap(" << block << ",false) is null");
  XPETRA_TEST_FOR_EXCEPTION(map_->getNumMaps() != full->getBlockedMap()->getNumMaps(),
                            Xpetra::Exceptions::RuntimeError,
                            "ExtractVector: Number of blocks in map extractor is " << map_->getNumMaps() << " but should be "
                                                                                   << full->getBlockedMap()->getNumMaps()
                                                                                   << " (number of blocks in BlockedMultiVector)");
  Teuchos::RCP<MultiVector> vv = full->getMultiVector(block, bThyraMode);
  return vv;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& full, size_t block, bool bThyraMode) const {
  XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                            std::out_of_range,
                            "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps()
                                                             << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "ExtractVector: map_->getmap(" << block << ",false) is null");
  XPETRA_TEST_FOR_EXCEPTION(map_->getNumMaps() != full->getBlockedMap()->getNumMaps(),
                            Xpetra::Exceptions::RuntimeError,
                            "ExtractVector: Number of blocks in map extractor is " << map_->getNumMaps() << " but should be "
                                                                                   << full->getBlockedMap()->getNumMaps()
                                                                                   << " (number of blocks in BlockedMultiVector)");
  Teuchos::RCP<MultiVector> vv = full->getMultiVector(block, bThyraMode);
  return vv;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& partial, size_t block, Vector& full, bool bThyraMode) const {
  XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                            std::out_of_range,
                            "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps()
                                                             << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "ExtractVector: map_->getmap(" << block << ",false) is null");
  XPETRA_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true,
                            Xpetra::Exceptions::RuntimeError,
                            "MapExtractor::InsertVector: InsertVector in Thyra-style numbering only possible if MapExtractor has been created "
                            "using Thyra-style numbered submaps.");
  if (bThyraMode) {
    // NOTE: the importer objects in the BlockedMap are always using Xpetra GIDs (or Thyra style Xpetra GIDs)
    // The source map corresponds to the full map (in Xpetra GIDs) starting with GIDs from zero. The GIDs are consecutive in Thyra mode
    // The target map is the partial map (in the corresponding Xpetra GIDs)

    // TODO can we skip the Export call in special cases (i.e. Src = Target map, same length, etc...)

    // store original GIDs (could be Thyra GIDs)
    RCP<const MultiVector> rcpPartial   = Teuchos::rcpFromRef(partial);
    RCP<MultiVector> rcpNonConstPartial = Teuchos::rcp_const_cast<MultiVector>(rcpPartial);
    RCP<const Map> oldThyMapPartial     = rcpNonConstPartial->getMap();  // temporarely store map of partial
    RCP<const Map> oldThyMapFull        = full.getMap();                 // temporarely store map of full

    // check whether getMap(block,false) is identical to target map of importer
    XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block, false)->isSameAs(*(map_->getImporter(block)->getTargetMap())) == false,
                              Xpetra::Exceptions::RuntimeError,
                              "MapExtractor::InsertVector: InsertVector in Thyra-style mode: Xpetra GIDs of partial vector are not identical "
                              "to target Map of Importer. This should not be.");

    // XPETRA_TEST_FOR_EXCEPTION(full.getMap()->isSameAs(*(map_->getImporter(block)->getSourceMap()))==false,
    // Xpetra::Exceptions::RuntimeError,
    //           "MapExtractor::InsertVector: InsertVector in Thyra-style mode: Xpetra GIDs of full vector are not identical to source Map of
    //           Importer. This should not be.");

    rcpNonConstPartial->replaceMap(getMap(block, false));       // temporarely switch to xpetra-style map
    full.replaceMap(map_->getImporter(block)->getSourceMap());  // temporarely switch to Xpetra GIDs

    // do the Export
    full.doExport(*rcpNonConstPartial, *(map_->getImporter(block)), Xpetra::INSERT);

    // switch back to original maps
    full.replaceMap(oldThyMapFull);                    // reset original map (Thyra GIDs)
    rcpNonConstPartial->replaceMap(oldThyMapPartial);  // change map back to original map
  } else {
    // Xpetra style numbering
    full.doExport(partial, *(map_->getImporter(block)), Xpetra::INSERT);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& partial, size_t block, MultiVector& full, bool bThyraMode) const {
  XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                            std::out_of_range,
                            "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps()
                                                             << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "ExtractVector: map_->getmap(" << block << ",false) is null");
  XPETRA_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true,
                            Xpetra::Exceptions::RuntimeError,
                            "MapExtractor::InsertVector: InsertVector in Thyra-style numbering only possible if MapExtractor has been created "
                            "using Thyra-style numbered submaps.");
  if (bThyraMode) {
    // NOTE: the importer objects in the BlockedMap are always using Xpetra GIDs (or Thyra style Xpetra GIDs)
    // The source map corresponds to the full map (in Xpetra GIDs) starting with GIDs from zero. The GIDs are consecutive in Thyra mode
    // The target map is the partial map (in the corresponding Xpetra GIDs)

    // TODO can we skip the Export call in special cases (i.e. Src = Target map, same length, etc...)

    // store original GIDs (could be Thyra GIDs)
    RCP<const MultiVector> rcpPartial   = Teuchos::rcpFromRef(partial);
    RCP<MultiVector> rcpNonConstPartial = Teuchos::rcp_const_cast<MultiVector>(rcpPartial);
    RCP<const Map> oldThyMapPartial     = rcpNonConstPartial->getMap();  // temporarely store map of partial
    RCP<const Map> oldThyMapFull        = full.getMap();                 // temporarely store map of full

    // check whether getMap(block,false) is identical to target map of importer
    XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block, false)->isSameAs(*(map_->getImporter(block)->getTargetMap())) == false,
                              Xpetra::Exceptions::RuntimeError,
                              "MapExtractor::InsertVector: InsertVector in Thyra-style mode: Xpetra GIDs of partial vector are not identical "
                              "to target Map of Importer. This should not be.");

    // XPETRA_TEST_FOR_EXCEPTION(full.getMap()->isSameAs(*(map_->getImporter(block)->getSourceMap()))==false,
    // Xpetra::Exceptions::RuntimeError,
    //           "MapExtractor::InsertVector: InsertVector in Thyra-style mode: Xpetra GIDs of full vector are not identical to source Map of
    //           Importer. This should not be.");

    rcpNonConstPartial->replaceMap(getMap(block, false));       // temporarely switch to xpetra-style map
    full.replaceMap(map_->getImporter(block)->getSourceMap());  // temporarely switch to Xpetra GIDs

    // do the Export
    full.doExport(*rcpNonConstPartial, *(map_->getImporter(block)), Xpetra::INSERT);

    // switch back to original maps
    full.replaceMap(oldThyMapFull);                    // reset original map (Thyra GIDs)
    rcpNonConstPartial->replaceMap(oldThyMapPartial);  // change map back to original map
  } else {
    // Xpetra style numbering
    full.doExport(partial, *(map_->getImporter(block)), Xpetra::INSERT);
  }
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
  RCP<BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(full);
  if (bfull.is_null() == true)
    InsertVector(*partial, block, *full, bThyraMode);
  else {
    XPETRA_TEST_FOR_EXCEPTION(
        map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "InsertVector: map_->getmap(" << block << ",false) is null");

#if 0
            // WCMCLEN - ETI: MultiVector::setMultiVector() doesn't exist.
            // WCMCLEN - ETI: but BlockedMultiVector::setMultiVector() does... should this be using bfull.
            full->setMultiVector(block, partial, bThyraMode);
#else
    throw std::runtime_error("Xpetra::MultiVector::setMultiVector() doesn't exist in " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
#endif
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> partial,
                 size_t block,
                 RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> full,
                 bool bThyraMode) const {
  RCP<BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(full);
  if (bfull.is_null() == true)
    InsertVector(*partial, block, *full, bThyraMode);
  else {
    XPETRA_TEST_FOR_EXCEPTION(
        map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "InsertVector: map_->getmap(" << block << ",false) is null");

    bfull->setMultiVector(block, partial, bThyraMode);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> partial,
                 size_t block,
                 RCP<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> full,
                 bool bThyraMode) const {
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "InsertVector: map_->getmap(" << block << ",false) is null");

  full->setMultiVector(block, partial, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> partial,
                 size_t block,
                 RCP<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> full,
                 bool bThyraMode) const {
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "InsertVector: map_->getmap(" << block << ",false) is null");
  full->setMultiVector(block, partial, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVector(size_t i, bool bThyraMode, bool bZero) const {
  XPETRA_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true,
                            Xpetra::Exceptions::RuntimeError,
                            "MapExtractor::getVector: getVector in Thyra-style numbering only possible if MapExtractor has been created using "
                            "Thyra-style numbered submaps.");
  // TODO check whether this can return a blocked multivector
  return Xpetra::VectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(getMap(i, bThyraMode), bZero);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVector(size_t i, size_t numvec, bool bThyraMode, bool bZero) const {
  XPETRA_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true,
                            Xpetra::Exceptions::RuntimeError,
                            "MapExtractor::getVector: getVector in Thyra-style numbering only possible if MapExtractor has been created using "
                            "Thyra-style numbered submaps.");
  // TODO check whether this can return a blocked multivector
  return MultiVectorFactory::Build(getMap(i, bThyraMode), numvec, bZero);
}

/// returns true, if sub maps are stored in Thyra-style numbering
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
const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMap(size_t i, bool bThyraMode) const {
  return map_->getMap(i, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMap() const {
  return map_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>
MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getBlockedMap() const {
  return map_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
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

}  // namespace Xpetra

#endif /* XPETRA_MAPEXTRACTOR_DEF_HPP_ */
