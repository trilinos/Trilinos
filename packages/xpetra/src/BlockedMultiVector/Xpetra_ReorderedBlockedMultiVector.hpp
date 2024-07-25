// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_REORDEREDBLOCKEDMULTIVECTOR_HPP
#define XPETRA_REORDEREDBLOCKEDMULTIVECTOR_HPP

#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_MapUtils.hpp"

#include "Xpetra_BlockReorderManager.hpp"
#include "Xpetra_BlockedMap.hpp"
#include "Xpetra_BlockedMultiVector.hpp"

/** \file Xpetra_ReorderedBlockedMultiVector.hpp

  Declarations for the class Xpetra::ReorderedBlockedMultiVector.
*/
namespace Xpetra {

typedef std::string viewLabel_t;

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class ReorderedBlockedMultiVector : public BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  typedef Scalar scalar_type;
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;

 private:
#undef XPETRA_REORDEREDBLOCKEDMULTIVECTOR_SHORT
#include "Xpetra_UseShortNames.hpp"

 public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor
  /*!
   * \param rangeMaps range maps for all blocks
   * \param domainMaps domain maps for all blocks
   * \param npr extimated number of entries per row in each block(!)
   * \param brm of type BlockReorderManager
   * \param bmat original full blocked operator (we keep the RCP to make sure all subblocks are available)
   */
  ReorderedBlockedMultiVector(Teuchos::RCP<const BlockedMap>& rangeMap,
                              Teuchos::RCP<const Xpetra::BlockReorderManager> brm,
                              Teuchos::RCP<const Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bvec)
    : Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rangeMap, bvec->getNumVectors(), false) {
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
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> mergeSubBlockMaps(Teuchos::RCP<const Xpetra::BlockReorderManager> brm) {
    RCP<const BlockedMap> bmap = fullVec_->getBlockedMap();

    // number of sub blocks
    size_t numBlocks = brm->GetNumBlocks();

    Teuchos::RCP<const Map> map = Teuchos::null;

    if (numBlocks == 0) {
      // it is a leaf node
      Teuchos::RCP<const Xpetra::BlockReorderLeaf> leaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(brm);

      // never extract Thyra style maps (since we have to merge them)
      map = bmap->getMap(Teuchos::as<size_t>(leaf->GetIndex()), false);
    } else {
      // initialize vector for sub maps
      std::vector<Teuchos::RCP<const Map>> subMaps(numBlocks, Teuchos::null);

      for (size_t i = 0; i < numBlocks; i++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> blkMgr = brm->GetBlock(Teuchos::as<int>(i));
        subMaps[i]                                             = mergeSubBlockMaps(blkMgr);
        TEUCHOS_ASSERT(subMaps[i].is_null() == false);
      }

      map = MapUtils::concatenateMaps(subMaps);
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
  Teuchos::RCP<const Xpetra::BlockReorderManager> brm_;
  Teuchos::RCP<const Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> fullVec_;
};

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> mergeSubBlockMaps(Teuchos::RCP<const Xpetra::BlockReorderManager> brm, Teuchos::RCP<const Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bvec, bool bThyraMode) {
  typedef Xpetra::MapUtils<LocalOrdinal, GlobalOrdinal, Node> MapUtils;

  // TODO distinguish between range and domain map extractor! provide MapExtractor as parameter!
  RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap = bvec->getBlockedMap();

  // number of sub blocks
  size_t numBlocks = brm->GetNumBlocks();

  Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> map = Teuchos::null;

  if (numBlocks == 0) {
    // it is a leaf node
    Teuchos::RCP<const Xpetra::BlockReorderLeaf> leaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(brm);

    map = bmap->getMap(Teuchos::as<size_t>(leaf->GetIndex()), bThyraMode);
  } else {
    // initialize vector for sub maps
    std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>> subMaps(numBlocks, Teuchos::null);

    for (size_t i = 0; i < numBlocks; i++) {
      Teuchos::RCP<const Xpetra::BlockReorderManager> blkMgr = brm->GetBlock(Teuchos::as<int>(i));
      subMaps[i]                                             = mergeSubBlockMaps(blkMgr, bvec, bThyraMode);
      TEUCHOS_ASSERT(subMaps[i].is_null() == false);
    }
#if 1
    // concatenate submaps
    // for Thyra mode this map isn't important
    Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> fullMap = MapUtils::concatenateMaps(subMaps);

    // create new BlockedMap (either in Thyra Mode or Xpetra mode)
    map = Teuchos::rcp(new Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>(fullMap, subMaps, bThyraMode));
#else
    // TAW: 11/27/16 we just concatenate the submaps to one monolithic Map object.
    // Alternatively, we could create a new BlockedMap using the concatenated map and the submaps
    // However, the block smoothers only need the concatenated map for creating MultiVectors...
    // But for the Thyra mode version concatenating would not be ok for the whole map
    map = MapUtils::concatenateMaps(subMaps);
#endif
  }
  TEUCHOS_ASSERT(map.is_null() == false);
  return map;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mergeSubBlocks(Teuchos::RCP<const Xpetra::BlockReorderManager> rowMgr, Teuchos::RCP<const Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bvec) {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node> BlockedMap;
  typedef Xpetra::MapUtils<LocalOrdinal, GlobalOrdinal, Node> MapUtils;
  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedMultiVector;
  typedef Xpetra::ReorderedBlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> ReorderedBlockedMultiVector;

  // number of sub blocks
  size_t rowSz = rowMgr->GetNumBlocks();

  Teuchos::RCP<BlockedMultiVector> rbvec = Teuchos::null;

  if (rowSz == 0) {
    // it is a leaf node
    Teuchos::RCP<const Xpetra::BlockReorderLeaf> rowleaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(rowMgr);

    // extract leaf node
    Teuchos::RCP<MultiVector> vec = bvec->getMultiVector(rowleaf->GetIndex(), false);

    TEUCHOS_ASSERT(vec != Teuchos::null);

    // check, whether leaf node is of type Xpetra::CrsMatrixWrap
    Teuchos::RCP<BlockedMultiVector> subBVec = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(vec);
    if (subBVec == Teuchos::null) {
      // DEBUG
      /*{
        std::cout << "MultiVector:" << std::endl;
        Teuchos::ArrayRCP<const Scalar> vData = vec->getData(0);
        for(size_t j=0; j< vec->getMap()->getLocalNumElements(); j++) {
            std::cout << j << ": " << vec->getMap()->getGlobalElement(j) << ": " << vData[j] << std::endl;
        }
      }*/
      // END DEBUG

      // If the leaf node is of type Xpetra::MultiVector. Wrap it into a ReorderedBlockMultiVector
      // with the corresponding MapExtractors for translating Thyra to Xpetra GIDs if necessary
      RCP<const BlockedMap> fullBlockedMap = bvec->getBlockedMap();
      Teuchos::RCP<const Map> submap       = fullBlockedMap->getMap(rowleaf->GetIndex(), false);
      std::vector<Teuchos::RCP<const Map>> rowSubMaps(1, submap);
      Teuchos::RCP<const BlockedMap> bbmap = Teuchos::rcp(new BlockedMap(submap, rowSubMaps, false));

      rbvec = Teuchos::rcp(new ReorderedBlockedMultiVector(bbmap, rowMgr, bvec));
      rbvec->setMultiVector(0, Teuchos::rcp_const_cast<MultiVector>(vec), false);

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
      Teuchos::RCP<const Xpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
      rowSubMaps[i]                                             = mergeSubBlockMaps(rowSubMgr, bvec, false /*xpetra*/);
      TEUCHOS_ASSERT(rowSubMaps[i].is_null() == false);
    }
    Teuchos::RCP<const Map> rgMergedSubMaps = MapUtils::concatenateMaps(rowSubMaps);
    rgBlockedMap                            = Teuchos::rcp(new BlockedMap(rgMergedSubMaps, rowSubMaps, false));
    rbvec                                   = Teuchos::rcp(new ReorderedBlockedMultiVector(rgBlockedMap, rowMgr, bvec));

    for (size_t i = 0; i < rowSz; i++) {
      Teuchos::RCP<const Xpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
      Teuchos::RCP<const MultiVector> subvec                    = mergeSubBlocks(rowSubMgr, bvec);
      rbvec->setMultiVector(i, Teuchos::rcp_const_cast<MultiVector>(subvec), false);
    }
  }
  return rbvec;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mergeSubBlocksThyra(Teuchos::RCP<const Xpetra::BlockReorderManager> rowMgr, Teuchos::RCP<const Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bvec) {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
  typedef Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node> BlockedMap;
  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedMultiVector;
  typedef Xpetra::ReorderedBlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> ReorderedBlockedMultiVector;

  TEUCHOS_ASSERT(bvec->getBlockedMap()->getThyraMode() == true);

  // number of sub blocks
  size_t rowSz = rowMgr->GetNumBlocks();

  Teuchos::RCP<BlockedMultiVector> rbvec = Teuchos::null;

  if (rowSz == 0) {
    // it is a leaf node
    Teuchos::RCP<const Xpetra::BlockReorderLeaf> rowleaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(rowMgr);

    // this MultiVector uses Thyra style GIDs as global row indices
    Teuchos::RCP<MultiVector> vec = bvec->getMultiVector(rowleaf->GetIndex(), true);

    TEUCHOS_ASSERT(vec.is_null() == false);

    // check, whether leaf node is of type Xpetra::CrsMatrixWrap
    Teuchos::RCP<BlockedMultiVector> bbvec = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(vec);
    if (bbvec == Teuchos::null) {
      ///////////////////////////////////////////////////////////////////////////
      // build blocked map
      RCP<const BlockedMap> fullBlockedRangeMap = bvec->getBlockedMap();
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
      Teuchos::RCP<const Xpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
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
      Teuchos::RCP<const Xpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
      Teuchos::RCP<const MultiVector> subvec                    = mergeSubBlocksThyra(rowSubMgr, bvec);
      rbvec->setMultiVector(i, Teuchos::rcp_const_cast<MultiVector>(subvec), true);
    }
  }
  return rbvec;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> buildReorderedBlockedMultiVector(Teuchos::RCP<const Xpetra::BlockReorderManager> brm, Teuchos::RCP<const Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bvec) {
  Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> rbvec = Teuchos::null;
  if (bvec->getBlockedMap()->getThyraMode() == false) {
    rbvec = mergeSubBlocks(brm, bvec);
  } else {
    rbvec = mergeSubBlocksThyra(brm, bvec);
  }
  TEUCHOS_ASSERT(rbvec.is_null() == false);
  return rbvec;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> buildReorderedBlockedMultiVector(Teuchos::RCP<const Xpetra::BlockReorderManager> brm, Teuchos::RCP<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bvec) {
  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> MultiVector;
  typedef Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedMultiVector;
  Teuchos::RCP<const MultiVector> rbvec = Teuchos::null;
  if (bvec->getBlockedMap()->getThyraMode() == false) {
    rbvec = mergeSubBlocks(brm, Teuchos::rcp_const_cast<const BlockedMultiVector>(bvec));
  } else {
    rbvec = mergeSubBlocksThyra(brm, Teuchos::rcp_const_cast<const BlockedMultiVector>(bvec));
  }
  TEUCHOS_ASSERT(rbvec.is_null() == false);
  return Teuchos::rcp_const_cast<MultiVector>(rbvec);
}

}  // namespace Xpetra

#define XPETRA_REORDEREDBLOCKEDMULTIVECTOR_SHORT
#endif /* XPETRA_REORDEREDBLOCKEDMULTIVECTOR_HPP */
