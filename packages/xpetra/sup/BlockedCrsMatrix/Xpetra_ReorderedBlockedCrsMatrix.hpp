// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_REORDEREDBLOCKEDCRSMATRIX_HPP
#define XPETRA_REORDEREDBLOCKEDCRSMATRIX_HPP

#include <Kokkos_DefaultNode.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_MapUtils.hpp"

#include "Xpetra_BlockedCrsMatrix.hpp"


/** \file Xpetra_ReorderedBlockedCrsMatrix.hpp

  Declarations for the class Xpetra::ReorderedBlockedCrsMatrix.
*/
namespace Xpetra {

  typedef std::string viewLabel_t;

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class ReorderedBlockedCrsMatrix :
    public BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  public:
    typedef Scalar scalar_type;
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;

  private:
#undef XPETRA_REORDEREDBLOCKEDCRSMATRIX_SHORT
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
     * \param pftype Xpetra profile type
     */
    ReorderedBlockedCrsMatrix
        (Teuchos::RCP<const MapExtractor>& rangeMaps,
         Teuchos::RCP<const MapExtractor>& domainMaps,
         size_t npr,
         Teuchos::RCP<const Xpetra::BlockReorderManager> brm,
         Teuchos::RCP<const Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > bmat,
         Xpetra::ProfileType pftype = Xpetra::DynamicProfile)
  : Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rangeMaps, domainMaps, npr,pftype) {
      brm_ = brm;
      fullOp_ = bmat;
    }

  //protected:

    //! Destructor
    virtual ~ReorderedBlockedCrsMatrix() {}

    //@}

  private:
    Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > mergeSubBlockMaps(Teuchos::RCP<const Xpetra::BlockReorderManager> brm) {
      RCP<const MapExtractor> fullRangeMapExtractor = fullOp_->getRangeMapExtractor();

      // number of sub blocks
      size_t numBlocks = brm->GetNumBlocks();

      Teuchos::RCP<const Map> map = Teuchos::null;

      if(numBlocks == 0) {
        // it is a leaf node
        Teuchos::RCP<const Xpetra::BlockReorderLeaf> leaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(brm);

        // never extract Thyra style maps (since we have to merge them)
        map = fullRangeMapExtractor->getMap(Teuchos::as<size_t>(leaf->GetIndex()), false);
      } else {
        // initialize vector for sub maps
        std::vector<Teuchos::RCP<const Map> > subMaps (numBlocks, Teuchos::null);

        for(size_t i = 0; i < numBlocks; i++) {
          Teuchos::RCP<const Xpetra::BlockReorderManager> blkMgr = brm->GetBlock(Teuchos::as<int>(i));
          subMaps[i] = mergeSubBlockMaps(blkMgr);
          TEUCHOS_ASSERT(subMaps[i].is_null()==false);
        }

        map = MapUtils::concatenateMaps(subMaps);
      }
      TEUCHOS_ASSERT(map.is_null()==false);
      return map;
    }

  public:
    //! @name Methods implementing Matrix
    //@{

    //! \brief Computes the sparse matrix-multivector multiplication.
    /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
      - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
      */
    virtual void apply(const MultiVector& X, MultiVector& Y,
                       Teuchos::ETransp mode = Teuchos::NO_TRANS,
                       Scalar alpha = ScalarTraits<Scalar>::one(),
                       Scalar beta  = ScalarTraits<Scalar>::zero()) const
    {
      Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply(X,Y,mode,alpha,beta);
    }

    // @}

    //! @name Overridden from Teuchos::Describable
    //@{

    /** \brief Return a simple one-line description of this object. */
    std::string description() const                       { return "ReorderedBlockedCrsMatrix"; }

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {
      typedef Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> BlockedCrsMatrix;

      out << "Xpetra::ReorderedBlockedCrsMatrix: " << BlockedCrsMatrix::Rows() << " x " << BlockedCrsMatrix::Cols() << std::endl;

      if (BlockedCrsMatrix::isFillComplete()) {
        out << "ReorderedBlockMatrix is fillComplete" << std::endl;

        out << "fullRowMap" << std::endl;
        BlockedCrsMatrix::getRangeMap(0,false)->describe(out,verbLevel);

        //out << "fullColMap" << std::endl;
        //fullcolmap_->describe(out,verbLevel);

      } else {
        out << "Xpetra::ReorderedBlockedCrsMatrix is NOT fillComplete" << std::endl;
      }

      for (size_t r = 0; r < BlockedCrsMatrix::Rows(); ++r)
        for (size_t c = 0; c < BlockedCrsMatrix::Cols(); ++c) {
          out << "Block(" << r << "," << c << ")" << std::endl;
          BlockedCrsMatrix::getMatrix(r,c)->describe(out,verbLevel);
        }
    }

    //@}

  private:
    Teuchos::RCP<const Xpetra::BlockReorderManager > brm_;
    Teuchos::RCP<const Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > fullOp_;


};

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > mergeSubBlockMaps(Teuchos::RCP<const Xpetra::BlockReorderManager> brm, Teuchos::RCP<const Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > bmat, bool bThyraMode) {
  typedef Xpetra::MapUtils<LocalOrdinal,GlobalOrdinal,Node> MapUtils;

  // TODO distinguish between range and domain map extractor! provide MapExtractor as parameter!
  RCP<const Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node> > fullRangeMapExtractor = bmat->getRangeMapExtractor();

  // number of sub blocks
  size_t numBlocks = brm->GetNumBlocks();

  Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>> map = Teuchos::null;

  if(numBlocks == 0) {
    // it is a leaf node
    Teuchos::RCP<const Xpetra::BlockReorderLeaf> leaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(brm);

    map = fullRangeMapExtractor->getMap(Teuchos::as<size_t>(leaf->GetIndex()), bThyraMode);
  } else {
    // initialize vector for sub maps
    std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > > subMaps (numBlocks, Teuchos::null);

    for(size_t i = 0; i < numBlocks; i++) {
      Teuchos::RCP<const Xpetra::BlockReorderManager> blkMgr = brm->GetBlock(Teuchos::as<int>(i));
      subMaps[i] = mergeSubBlockMaps(blkMgr,bmat,bThyraMode);
      TEUCHOS_ASSERT(subMaps[i].is_null()==false);
    }

    map = MapUtils::concatenateMaps(subMaps);
  }
  TEUCHOS_ASSERT(map.is_null()==false);
  return map;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mergeSubBlocks(Teuchos::RCP<const Xpetra::BlockReorderManager> rowMgr, Teuchos::RCP<const Xpetra::BlockReorderManager> colMgr, Teuchos::RCP<const Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > bmat) {

  typedef Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Map;
  typedef Xpetra::MapUtils<LocalOrdinal,GlobalOrdinal,Node> MapUtils;
  typedef Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node> MapExtractor;
  typedef Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Matrix;
  typedef Xpetra::ReorderedBlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> ReorderedBlockedCrsMatrix;

  // number of sub blocks
  size_t rowSz = rowMgr->GetNumBlocks();
  size_t colSz = colMgr->GetNumBlocks();

  Teuchos::RCP<ReorderedBlockedCrsMatrix> rbmat = Teuchos::null;

  if(rowSz == 0 && colSz == 0) {
    // it is a leaf node
    Teuchos::RCP<const Xpetra::BlockReorderLeaf> rowleaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(rowMgr);
    Teuchos::RCP<const Xpetra::BlockReorderLeaf> colleaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(colMgr);

    RCP<const MapExtractor> fullRangeMapExtractor = bmat->getRangeMapExtractor();
    Teuchos::RCP<const Map> submap = fullRangeMapExtractor->getMap(rowleaf->GetIndex(),false);
    std::vector<Teuchos::RCP<const Map> > rowSubMaps (1, submap);
    Teuchos::RCP<const MapExtractor> rgMapExtractor = Teuchos::rcp(new MapExtractor(submap, rowSubMaps, false));

    RCP<const MapExtractor> fullDomainMapExtractor = bmat->getDomainMapExtractor();
    Teuchos::RCP<const Map> submap2 = fullDomainMapExtractor->getMap(colleaf->GetIndex(),false);
    std::vector<Teuchos::RCP<const Map> > colSubMaps (1, submap2);
    Teuchos::RCP<const MapExtractor> doMapExtractor = Teuchos::rcp(new MapExtractor(submap2, colSubMaps, false));

    Teuchos::RCP<Matrix> mat = bmat->getMatrix(rowleaf->GetIndex(), colleaf->GetIndex());

    rbmat = Teuchos::rcp(new ReorderedBlockedCrsMatrix(rgMapExtractor,doMapExtractor, 33, rowMgr,bmat));
    rbmat->setMatrix(0,0,mat);
  } else {
    // create the map extractors
    // we cannot create block matrix in thyra mode since merged maps might not start with 0 GID
    Teuchos::RCP<const MapExtractor> rgMapExtractor = Teuchos::null;
    if(rowSz > 0) {
      std::vector<Teuchos::RCP<const Map> > rowSubMaps (rowSz, Teuchos::null);
      for(size_t i = 0; i < rowSz; i++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
        rowSubMaps[i] = mergeSubBlockMaps(rowSubMgr,bmat,false /*xpetra*/);
        TEUCHOS_ASSERT(rowSubMaps[i].is_null()==false);
      }
      Teuchos::RCP<const Map> rgMergedSubMaps = MapUtils::concatenateMaps(rowSubMaps);
      rgMapExtractor = Teuchos::rcp(new MapExtractor(rgMergedSubMaps, rowSubMaps, false));
    } else {
      Teuchos::RCP<const Xpetra::BlockReorderLeaf> rowleaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(rowMgr);
      RCP<const MapExtractor> fullRangeMapExtractor = bmat->getRangeMapExtractor();
      // TODO think about Thyra style maps: we cannot use thyra style maps when recombining several blocks!!!
      // The GIDs might not start with 0 and may not be consecutive!
      Teuchos::RCP<const Map> submap = fullRangeMapExtractor->getMap(rowleaf->GetIndex(),false);
      std::vector<Teuchos::RCP<const Map> > rowSubMaps (1, submap);
      rgMapExtractor = Teuchos::rcp(new MapExtractor(submap, rowSubMaps, false));
    }

    Teuchos::RCP<const MapExtractor> doMapExtractor = Teuchos::null;
    if(colSz > 0) {
      std::vector<Teuchos::RCP<const Map> > colSubMaps (colSz, Teuchos::null);
      for(size_t j = 0; j < colSz; j++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> colSubMgr = colMgr->GetBlock(Teuchos::as<int>(j));
        colSubMaps[j] = mergeSubBlockMaps(colSubMgr,bmat,false/*xpetra*/);
        TEUCHOS_ASSERT(colSubMaps[j].is_null()==false);
      }
      Teuchos::RCP<const Map> doMergedSubMaps = MapUtils::concatenateMaps(colSubMaps);
      doMapExtractor = Teuchos::rcp(new MapExtractor(doMergedSubMaps, colSubMaps, false));
    } else {
      Teuchos::RCP<const Xpetra::BlockReorderLeaf> colleaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(colMgr);
      RCP<const MapExtractor> fullDomainMapExtractor = bmat->getDomainMapExtractor();
      // TODO think about Thyra style maps: we cannot use thyra style maps when recombining several blocks!!!
      // The GIDs might not start with 0 and may not be consecutive!
      Teuchos::RCP<const Map> submap = fullDomainMapExtractor->getMap(colleaf->GetIndex(),false);
      std::vector<Teuchos::RCP<const Map> > colSubMaps (1, submap);
      doMapExtractor = Teuchos::rcp(new MapExtractor(submap, colSubMaps, false));
    }

    rbmat = Teuchos::rcp(new ReorderedBlockedCrsMatrix(rgMapExtractor,doMapExtractor, 33, rowMgr,bmat));

    if (rowSz == 0 && colSz > 0) {
      for(size_t j = 0; j < colSz; j++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> colSubMgr = colMgr->GetBlock(Teuchos::as<int>(j));
        Teuchos::RCP<const Matrix> submat = mergeSubBlocks(rowMgr, colSubMgr, bmat);
        rbmat->setMatrix(0,j,Teuchos::rcp_const_cast<Matrix>(submat));
      }
    } else if (rowSz > 0 && colSz == 0) {
      for(size_t i = 0; i < rowSz; i++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
        Teuchos::RCP<const Matrix> submat = mergeSubBlocks(rowSubMgr, colMgr, bmat);
        rbmat->setMatrix(i,0,Teuchos::rcp_const_cast<Matrix>(submat));
      }
    } else {
      for(size_t i = 0; i < rowSz; i++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
        for(size_t j = 0; j < colSz; j++) {
          Teuchos::RCP<const Xpetra::BlockReorderManager> colSubMgr = colMgr->GetBlock(Teuchos::as<int>(j));
          Teuchos::RCP<const Matrix> submat = mergeSubBlocks(rowSubMgr, colSubMgr, bmat);
          rbmat->setMatrix(i,j,Teuchos::rcp_const_cast<Matrix>(submat));
        }
      }
    }
  }
  rbmat->fillComplete();
  return rbmat;
}

  //MapExtractor(const std::vector<RCP<const Map> >& maps, const std::vector<RCP<const Map> >& thyramaps);

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mergeSubBlocksThyra(Teuchos::RCP<const Xpetra::BlockReorderManager> rowMgr, Teuchos::RCP<const Xpetra::BlockReorderManager> colMgr, Teuchos::RCP<const Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > bmat) {

  typedef Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> Map;
  typedef Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node> MapExtractor;
  typedef Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> Matrix;
  typedef Xpetra::ReorderedBlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> ReorderedBlockedCrsMatrix;

  TEUCHOS_ASSERT(bmat->getRangeMapExtractor()->getThyraMode() == true);
  TEUCHOS_ASSERT(bmat->getDomainMapExtractor()->getThyraMode() == true);

  // number of sub blocks
  size_t rowSz = rowMgr->GetNumBlocks();
  size_t colSz = colMgr->GetNumBlocks();

  Teuchos::RCP<ReorderedBlockedCrsMatrix> rbmat = Teuchos::null;

  if(rowSz == 0 && colSz == 0) {
    // it is a leaf node
    Teuchos::RCP<const Xpetra::BlockReorderLeaf> rowleaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(rowMgr);
    Teuchos::RCP<const Xpetra::BlockReorderLeaf> colleaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(colMgr);

    ///////////////////////////////////////////////////////////////////////////
    // build map extractors
    RCP<const MapExtractor> fullRangeMapExtractor = bmat->getRangeMapExtractor();
    // extract Xpetra and Thyra based GIDs
    Teuchos::RCP<const Map> xpsubmap  = fullRangeMapExtractor->getMap(rowleaf->GetIndex(),false);
    Teuchos::RCP<const Map> thysubmap = fullRangeMapExtractor->getMap(rowleaf->GetIndex(),true);
    std::vector<Teuchos::RCP<const Map> > rowXpSubMaps (1, xpsubmap);
    std::vector<Teuchos::RCP<const Map> > rowTySubMaps (1, thysubmap);
    // use expert constructor
    Teuchos::RCP<const MapExtractor> rgMapExtractor = Teuchos::rcp(new MapExtractor(rowXpSubMaps, rowTySubMaps));

    RCP<const MapExtractor> fullDomainMapExtractor = bmat->getDomainMapExtractor();
    // extract Xpetra and Thyra based GIDs
    Teuchos::RCP<const Map> xpsubmap2 = fullDomainMapExtractor->getMap(colleaf->GetIndex(),false);
    Teuchos::RCP<const Map> tysubmap2 = fullDomainMapExtractor->getMap(colleaf->GetIndex(),true);
    std::vector<Teuchos::RCP<const Map> > colXpSubMaps (1, xpsubmap2);
    std::vector<Teuchos::RCP<const Map> > colTySubMaps (1, tysubmap2);
    // use expert constructor
    Teuchos::RCP<const MapExtractor> doMapExtractor = Teuchos::rcp(new MapExtractor(colXpSubMaps, colTySubMaps));

    ///////////////////////////////////////////////////////////////////////////
    // build reordered block operator

    // this matrix uses Thyra style GIDs as global row, range, domain and column indices
    Teuchos::RCP<Matrix> mat = bmat->getMatrix(rowleaf->GetIndex(), colleaf->GetIndex());

    rbmat = Teuchos::rcp(new ReorderedBlockedCrsMatrix(rgMapExtractor,doMapExtractor, 33, rowMgr,bmat));
    rbmat->setMatrix(0,0,mat);
  } else {
    // create the map extractors
    // we cannot create block matrix in thyra mode since merged maps might not start with 0 GID
    Teuchos::RCP<const MapExtractor> rgMapExtractor = Teuchos::null;
    if(rowSz > 0) {
      std::vector<Teuchos::RCP<const Map> > rowXpSubMaps (rowSz, Teuchos::null);
      std::vector<Teuchos::RCP<const Map> > rowTySubMaps (rowSz, Teuchos::null);
      for(size_t i = 0; i < rowSz; i++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
        // extract Xpetra and Thyra based merged GIDs
        rowXpSubMaps[i] = mergeSubBlockMaps(rowSubMgr,bmat,false);
        rowTySubMaps[i] = mergeSubBlockMaps(rowSubMgr,bmat,true);
        TEUCHOS_ASSERT(rowXpSubMaps[i].is_null()==false);
        TEUCHOS_ASSERT(rowTySubMaps[i].is_null()==false);
      }
      // use expert constructor
      rgMapExtractor = Teuchos::rcp(new MapExtractor(rowXpSubMaps, rowTySubMaps));
    } else {
      Teuchos::RCP<const Xpetra::BlockReorderLeaf> rowleaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(rowMgr);
      RCP<const MapExtractor> fullRangeMapExtractor = bmat->getRangeMapExtractor();
      // extract Xpetra and Thyra based GIDs
      Teuchos::RCP<const Map> xpsubmap  = fullRangeMapExtractor->getMap(rowleaf->GetIndex(),false);
      Teuchos::RCP<const Map> thysubmap = fullRangeMapExtractor->getMap(rowleaf->GetIndex(),true);
      std::vector<Teuchos::RCP<const Map> > rowXpSubMaps (1, xpsubmap);
      std::vector<Teuchos::RCP<const Map> > rowTySubMaps (1, thysubmap);
      // use expert constructor
      rgMapExtractor = Teuchos::rcp(new MapExtractor(rowXpSubMaps, rowTySubMaps));
    }

    Teuchos::RCP<const MapExtractor> doMapExtractor = Teuchos::null;
    if(colSz > 0) {
      std::vector<Teuchos::RCP<const Map> > colXpSubMaps (colSz, Teuchos::null);
      std::vector<Teuchos::RCP<const Map> > colTySubMaps (colSz, Teuchos::null);
      for(size_t j = 0; j < colSz; j++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> colSubMgr = colMgr->GetBlock(Teuchos::as<int>(j));
        // extract Xpetra and Thyra based merged GIDs
        colXpSubMaps[j] = mergeSubBlockMaps(colSubMgr,bmat,false);
        colTySubMaps[j] = mergeSubBlockMaps(colSubMgr,bmat,true);
        TEUCHOS_ASSERT(colXpSubMaps[j].is_null()==false);
        TEUCHOS_ASSERT(colTySubMaps[j].is_null()==false);
      }
      // use expert constructor
      doMapExtractor = Teuchos::rcp(new MapExtractor(colXpSubMaps,colTySubMaps));
    } else {
      Teuchos::RCP<const Xpetra::BlockReorderLeaf> colleaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(colMgr);
      RCP<const MapExtractor> fullDomainMapExtractor = bmat->getDomainMapExtractor();
      // extract Xpetra and Thyra based GIDs
      Teuchos::RCP<const Map> xpsubmap = fullDomainMapExtractor->getMap(colleaf->GetIndex(),false);
      Teuchos::RCP<const Map> tysubmap = fullDomainMapExtractor->getMap(colleaf->GetIndex(),true);
      std::vector<Teuchos::RCP<const Map> > colXpSubMaps (1, xpsubmap);
      std::vector<Teuchos::RCP<const Map> > colTySubMaps (1, tysubmap);
      // use expert constructor
      doMapExtractor = Teuchos::rcp(new MapExtractor(colXpSubMaps, colTySubMaps));
    }

    // TODO matrix should have both rowMgr and colMgr??
    rbmat = Teuchos::rcp(new ReorderedBlockedCrsMatrix(rgMapExtractor,doMapExtractor, 33, rowMgr,bmat));

    if (rowSz == 0 && colSz > 0) {
      for(size_t j = 0; j < colSz; j++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> colSubMgr = colMgr->GetBlock(Teuchos::as<int>(j));
        Teuchos::RCP<const Matrix> submat = mergeSubBlocksThyra(rowMgr, colSubMgr, bmat);
        rbmat->setMatrix(0,j,Teuchos::rcp_const_cast<Matrix>(submat));
      }
    } else if (rowSz > 0 && colSz == 0) {
      for(size_t i = 0; i < rowSz; i++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
        Teuchos::RCP<const Matrix> submat = mergeSubBlocksThyra(rowSubMgr, colMgr, bmat);
        rbmat->setMatrix(i,0,Teuchos::rcp_const_cast<Matrix>(submat));
      }
    } else {
      for(size_t i = 0; i < rowSz; i++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
        for(size_t j = 0; j < colSz; j++) {
          Teuchos::RCP<const Xpetra::BlockReorderManager> colSubMgr = colMgr->GetBlock(Teuchos::as<int>(j));
          Teuchos::RCP<const Matrix> submat = mergeSubBlocksThyra(rowSubMgr, colSubMgr, bmat);
          rbmat->setMatrix(i,j,Teuchos::rcp_const_cast<Matrix>(submat));
        }
      }
    }
  }

  rbmat->fillComplete();
  return rbmat;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > buildReorderedBlockedCrsMatrix(Teuchos::RCP<const Xpetra::BlockReorderManager> brm, Teuchos::RCP<const Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > bmat) {
  TEUCHOS_ASSERT(bmat->getRangeMapExtractor()->getThyraMode() == bmat->getDomainMapExtractor()->getThyraMode());
  Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rbmat = Teuchos::null;
  if(bmat->getRangeMapExtractor()->getThyraMode() == false) {
    rbmat = mergeSubBlocks(brm, brm, bmat);
  } else {
    rbmat = mergeSubBlocksThyra(brm, brm, bmat);
  }

  TEUCHOS_ASSERT(rbmat != Teuchos::null);
  return rbmat;
}

} //namespace Xpetra

#define XPETRA_REORDEREDBLOCKEDCRSMATRIX_SHORT
#endif /* XPETRA_REORDEREDBLOCKEDCRSMATRIX_HPP */
