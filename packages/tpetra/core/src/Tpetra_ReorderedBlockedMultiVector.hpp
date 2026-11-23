// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_REORDEREDBLOCKEDMULTIVECTOR_HPP
#define TPETRA_REORDEREDBLOCKEDMULTIVECTOR_HPP

/// \file Tpetra_ReorderedBlockedMultiVector.hpp
/// \brief Declaration of the Tpetra::ReorderedBlockedMultiVector class.
///
/// Direct port of Xpetra::ReorderedBlockedMultiVector.  This is a header-only
/// class (and set of free functions) that applies a Tpetra::BlockReorderManager
/// tree to a Tpetra::BlockedMultiVector, producing a new (nested) blocked
/// multivector whose block structure follows the reorder manager.  The Xpetra
/// version relied on Xpetra::MapUtils::concatenateMaps; here we use the static
/// Tpetra::BlockedMap::concatenateMaps instead.  Leaf detection uses
/// rcp_dynamic_cast to Tpetra::BlockedMultiVector (there is no Xpetra
/// CrsMatrixWrap analog in Tpetra).

#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include "Tpetra_ConfigDefs.hpp"

#include "Tpetra_BlockReorderManager.hpp"
#include "Tpetra_BlockedMap_decl.hpp"
#include "Tpetra_BlockedMultiVector_decl.hpp"

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = ::Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class ReorderedBlockedMultiVector : public BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  typedef Scalar scalar_type;
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;

 private:
  typedef ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef ::Tpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node> BlockedMap;
  typedef ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;

 public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor
  /*!
   * \param rangeMap blocked map describing the (reordered) block structure
   * \param brm of type BlockReorderManager
   * \param bvec original full blocked multivector (we keep the RCP to make sure all subblocks are available)
   */
  ReorderedBlockedMultiVector(Teuchos::RCP<const BlockedMap>& rangeMap,
                              Teuchos::RCP<const ::Tpetra::BlockReorderManager> brm,
                              Teuchos::RCP<const ::Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bvec)
    : ::Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rangeMap, bvec->getNumVectors(), false) {
    brm_     = brm;
    fullVec_ = bvec;
  }

  // protected:

  //! Destructor
  virtual ~ReorderedBlockedMultiVector() {
    brm_     = Teuchos::null;
    fullVec_ = Teuchos::null;
  }

  //@}

 private:
  Teuchos::RCP<const Map> mergeSubBlockMaps(Teuchos::RCP<const ::Tpetra::BlockReorderManager> brm) {
    Teuchos::RCP<const BlockedMap> bmap = fullVec_->getBlockedMap();

    // number of sub blocks
    size_t numBlocks = brm->GetNumBlocks();

    Teuchos::RCP<const Map> map = Teuchos::null;

    if (numBlocks == 0) {
      // it is a leaf node
      Teuchos::RCP<const ::Tpetra::BlockReorderLeaf> leaf = Teuchos::rcp_dynamic_cast<const ::Tpetra::BlockReorderLeaf>(brm);

      // never extract Thyra style maps (since we have to merge them)
      map = bmap->getMap(Teuchos::as<size_t>(leaf->GetIndex()), false);
    } else {
      // initialize vector for sub maps
      std::vector<Teuchos::RCP<const Map>> subMaps(numBlocks, Teuchos::null);

      for (size_t i = 0; i < numBlocks; i++) {
        Teuchos::RCP<const ::Tpetra::BlockReorderManager> blkMgr = brm->GetBlock(Teuchos::as<int>(i));
        subMaps[i]                                               = mergeSubBlockMaps(blkMgr);
        TEUCHOS_ASSERT(subMaps[i].is_null() == false);
      }

      map = BlockedMap::concatenateMaps(subMaps);
    }
    TEUCHOS_ASSERT(map.is_null() == false);
    return map;
  }

 public:
  //! @name Overridden from Teuchos::Describable
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const { return "ReorderedBlockedMultiVector"; }

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {
    TEUCHOS_ASSERT(brm_ != Teuchos::null);
    out << description() << ": " << brm_->toString() << std::endl;
    fullVec_->describe(out, verbLevel);
  }

  //@}

 private:
  Teuchos::RCP<const ::Tpetra::BlockReorderManager> brm_;
  Teuchos::RCP<const ::Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> fullVec_;
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> mergeSubBlockMaps(Teuchos::RCP<const ::Tpetra::BlockReorderManager> brm, Teuchos::RCP<const ::Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bvec, bool bThyraMode) {
  typedef ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef ::Tpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node> BlockedMap;

  // TODO distinguish between range and domain map extractor! provide MapExtractor as parameter!
  Teuchos::RCP<const BlockedMap> bmap = bvec->getBlockedMap();

  // number of sub blocks
  size_t numBlocks = brm->GetNumBlocks();

  Teuchos::RCP<const Map> map = Teuchos::null;

  if (numBlocks == 0) {
    // it is a leaf node
    Teuchos::RCP<const ::Tpetra::BlockReorderLeaf> leaf = Teuchos::rcp_dynamic_cast<const ::Tpetra::BlockReorderLeaf>(brm);

    map = bmap->getMap(Teuchos::as<size_t>(leaf->GetIndex()), bThyraMode);
  } else {
    // initialize vector for sub maps
    std::vector<Teuchos::RCP<const Map>> subMaps(numBlocks, Teuchos::null);

    for (size_t i = 0; i < numBlocks; i++) {
      Teuchos::RCP<const ::Tpetra::BlockReorderManager> blkMgr = brm->GetBlock(Teuchos::as<int>(i));
      subMaps[i]                                               = mergeSubBlockMaps(blkMgr, bvec, bThyraMode);
      TEUCHOS_ASSERT(subMaps[i].is_null() == false);
    }

    // concatenate submaps
    // for Thyra mode this map isn't important
    // NOTE (Tpetra port): Xpetra wrapped the concatenated map in a BlockedMap here
    // (BlockedMap is-a Map in Xpetra).  Tpetra::BlockedMap is standalone (not a Map),
    // so a merged sub-block map is just the concatenated plain Map -- the nested block
    // structure is carried by the (Reordered)BlockedMultiVector, not by this map.
    map = BlockedMap::concatenateMaps(subMaps);
  }
  TEUCHOS_ASSERT(map.is_null() == false);
  return map;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mergeSubBlocks(Teuchos::RCP<const ::Tpetra::BlockReorderManager> rowMgr, Teuchos::RCP<const ::Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bvec) {
  typedef ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef ::Tpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node> BlockedMap;
  typedef ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
  typedef ::Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedMultiVector;
  typedef ::Tpetra::ReorderedBlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> ReorderedBlockedMultiVector;

  // number of sub blocks
  size_t rowSz = rowMgr->GetNumBlocks();

  Teuchos::RCP<BlockedMultiVector> rbvec = Teuchos::null;

  if (rowSz == 0) {
    // it is a leaf node
    Teuchos::RCP<const ::Tpetra::BlockReorderLeaf> rowleaf = Teuchos::rcp_dynamic_cast<const ::Tpetra::BlockReorderLeaf>(rowMgr);

    // extract leaf node
    Teuchos::RCP<MultiVector> vec = bvec->getMultiVector(rowleaf->GetIndex(), false);

    TEUCHOS_ASSERT(vec != Teuchos::null);

    // check, whether leaf node is of type Tpetra::BlockedMultiVector
    Teuchos::RCP<BlockedMultiVector> subBVec = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(vec);
    if (subBVec == Teuchos::null) {
      // If the leaf node is of type Tpetra::MultiVector. Wrap it into a ReorderedBlockMultiVector
      // with the corresponding MapExtractors for translating Thyra to Xpetra GIDs if necessary
      Teuchos::RCP<const BlockedMap> fullBlockedMap = bvec->getBlockedMap();
      Teuchos::RCP<const Map> submap                = fullBlockedMap->getMap(rowleaf->GetIndex(), false);
      std::vector<Teuchos::RCP<const Map>> rowSubMaps(1, submap);
      Teuchos::RCP<const BlockedMap> bbmap = Teuchos::rcp(new BlockedMap(submap, rowSubMaps, false));

      rbvec = Teuchos::rcp(new ReorderedBlockedMultiVector(bbmap, rowMgr, bvec));
      rbvec->setMultiVector(0, Teuchos::rcp_const_cast<const MultiVector>(vec), false);

    } else {
      // If leaf node is already wrapped into a blocked matrix do not wrap it again.
      rbvec = subBVec;
      TEUCHOS_ASSERT(rbvec != Teuchos::null);
    }
  } else {
    // create the map extractors
    // we cannot create block matrix in thyra mode since merged maps might not start with 0 GID
    Teuchos::RCP<const BlockedMap> rgBlockedMap = Teuchos::null;
    std::vector<Teuchos::RCP<const Map>> rowSubMaps(rowSz, Teuchos::null);
    for (size_t i = 0; i < rowSz; i++) {
      Teuchos::RCP<const ::Tpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
      rowSubMaps[i]                                               = mergeSubBlockMaps(rowSubMgr, bvec, false /*xpetra*/);
      TEUCHOS_ASSERT(rowSubMaps[i].is_null() == false);
    }
    Teuchos::RCP<const Map> rgMergedSubMaps = BlockedMap::concatenateMaps(rowSubMaps);
    rgBlockedMap                            = Teuchos::rcp(new BlockedMap(rgMergedSubMaps, rowSubMaps, false));
    rbvec                                   = Teuchos::rcp(new ReorderedBlockedMultiVector(rgBlockedMap, rowMgr, bvec));

    for (size_t i = 0; i < rowSz; i++) {
      Teuchos::RCP<const ::Tpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
      Teuchos::RCP<const MultiVector> subvec                      = mergeSubBlocks(rowSubMgr, bvec);
      rbvec->setMultiVector(i, subvec, false);
    }
  }
  return rbvec;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mergeSubBlocksThyra(Teuchos::RCP<const ::Tpetra::BlockReorderManager> rowMgr, Teuchos::RCP<const ::Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bvec) {
  typedef ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef ::Tpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node> BlockedMap;
  typedef ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
  typedef ::Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedMultiVector;
  typedef ::Tpetra::ReorderedBlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> ReorderedBlockedMultiVector;

  TEUCHOS_ASSERT(bvec->getBlockedMap()->getThyraMode() == true);

  // number of sub blocks
  size_t rowSz = rowMgr->GetNumBlocks();

  Teuchos::RCP<BlockedMultiVector> rbvec = Teuchos::null;

  if (rowSz == 0) {
    // it is a leaf node
    Teuchos::RCP<const ::Tpetra::BlockReorderLeaf> rowleaf = Teuchos::rcp_dynamic_cast<const ::Tpetra::BlockReorderLeaf>(rowMgr);

    // this MultiVector uses Thyra style GIDs as global row indices
    Teuchos::RCP<MultiVector> vec = bvec->getMultiVector(rowleaf->GetIndex(), true);

    TEUCHOS_ASSERT(vec.is_null() == false);

    // check, whether leaf node is of type Tpetra::BlockedMultiVector
    Teuchos::RCP<BlockedMultiVector> bbvec = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(vec);
    if (bbvec == Teuchos::null) {
      ///////////////////////////////////////////////////////////////////////////
      // build blocked map
      Teuchos::RCP<const BlockedMap> fullBlockedRangeMap = bvec->getBlockedMap();
      // extract Xpetra and Thyra based GIDs
      Teuchos::RCP<const Map> xpsubmap  = fullBlockedRangeMap->getMap(rowleaf->GetIndex(), false);
      Teuchos::RCP<const Map> thysubmap = fullBlockedRangeMap->getMap(rowleaf->GetIndex(), true);
      std::vector<Teuchos::RCP<const Map>> rowXpSubMaps(1, xpsubmap);
      std::vector<Teuchos::RCP<const Map>> rowTySubMaps(1, thysubmap);
      // use expert constructor
      Teuchos::RCP<const BlockedMap> rgBlockedMap = Teuchos::rcp(new BlockedMap(rowXpSubMaps, rowTySubMaps));

      ///////////////////////////////////////////////////////////////////////////
      // build reordered blocked multi vector
      rbvec = Teuchos::rcp(new ReorderedBlockedMultiVector(rgBlockedMap, rowMgr, bvec));
      rbvec->setMultiVector(0, vec, true);
    } else {
      // If leaf node is already wrapped into a blocked matrix do not wrap it again.
      rbvec = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(vec);
    }
  } else {
    // create the blocked map
    // we cannot create block multivector in thyra mode since merged maps might not start with 0 GID

    std::vector<Teuchos::RCP<const Map>> rowXpSubMaps(rowSz, Teuchos::null);
    std::vector<Teuchos::RCP<const Map>> rowTySubMaps(rowSz, Teuchos::null);
    for (size_t i = 0; i < rowSz; i++) {
      Teuchos::RCP<const ::Tpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
      // extract Xpetra and Thyra based merged GIDs
      rowXpSubMaps[i] = mergeSubBlockMaps(rowSubMgr, bvec, false);
      rowTySubMaps[i] = mergeSubBlockMaps(rowSubMgr, bvec, true);
      TEUCHOS_ASSERT(rowXpSubMaps[i].is_null() == false);
      TEUCHOS_ASSERT(rowTySubMaps[i].is_null() == false);
    }
    // use expert constructor
    Teuchos::RCP<const BlockedMap> rgBlockedMap = Teuchos::rcp(new BlockedMap(rowXpSubMaps, rowTySubMaps));

    rbvec = Teuchos::rcp(new ReorderedBlockedMultiVector(rgBlockedMap, rowMgr, bvec));

    for (size_t i = 0; i < rowSz; i++) {
      Teuchos::RCP<const ::Tpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
      Teuchos::RCP<const MultiVector> subvec                      = mergeSubBlocksThyra(rowSubMgr, bvec);
      rbvec->setMultiVector(i, subvec, true);
    }
  }
  return rbvec;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> buildReorderedBlockedMultiVector(Teuchos::RCP<const ::Tpetra::BlockReorderManager> brm, Teuchos::RCP<const ::Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bvec) {
  Teuchos::RCP<const ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> rbvec = Teuchos::null;
  if (bvec->getBlockedMap()->getThyraMode() == false) {
    rbvec = mergeSubBlocks(brm, bvec);
  } else {
    rbvec = mergeSubBlocksThyra(brm, bvec);
  }
  TEUCHOS_ASSERT(rbvec.is_null() == false);
  return rbvec;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> buildReorderedBlockedMultiVector(Teuchos::RCP<const ::Tpetra::BlockReorderManager> brm, Teuchos::RCP<::Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bvec) {
  typedef ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
  typedef ::Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedMultiVector;
  Teuchos::RCP<const MultiVector> rbvec = Teuchos::null;
  if (bvec->getBlockedMap()->getThyraMode() == false) {
    rbvec = mergeSubBlocks(brm, Teuchos::rcp_const_cast<const BlockedMultiVector>(bvec));
  } else {
    rbvec = mergeSubBlocksThyra(brm, Teuchos::rcp_const_cast<const BlockedMultiVector>(bvec));
  }
  TEUCHOS_ASSERT(rbvec.is_null() == false);
  return Teuchos::rcp_const_cast<MultiVector>(rbvec);
}

}  // namespace Tpetra

#endif /* TPETRA_REORDEREDBLOCKEDMULTIVECTOR_HPP */
