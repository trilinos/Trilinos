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

#if 0
      // extract map extractors of full block operator
      RCP<const MapExtractor> fullRangeMapExtractor = fullOp_->getRangeMapExtractor();
      RCP<const MapExtractor> fullDomainMapExtractor = fullOp_->getDomainMapExtractor();

      // create the map extractors from brm
      size_t numBlocks = brm->GetNumBlocks();
      std::vector<Teuchos::RCP<const Map> > subMaps (numBlocks, Teuchos::null);
      for(size_t i = 0; i < numBlocks; i++) {
        subMaps[i] = mergeSubBlockMaps(brm_);
        TEUCHOS_ASSERT(subMaps[i].is_null()==false);
      }

      // TODO: check thyra mode
      // we cannot create block matrix in thyra mode since merged maps might not start with 0 GID
      Teuchos::RCP<const Map> rgMergedSubMaps = MapUtils::concatenateMaps(subMaps);
      Teuchos::RCP<const MapExtractor> rgMapExtractor = Teuchos::rcp(new MapExtractor(rgMergedSubMaps, subMaps, false /*fullRangeMapExtractor->getThyraMode()*/));
      Teuchos::RCP<const Map> doMergedSubMaps = MapUtils::concatenateMaps(subMaps);
      Teuchos::RCP<const MapExtractor> doMapExtractor = Teuchos::rcp(new MapExtractor(doMergedSubMaps, subMaps, false /*fullDomainMapExtractor->getThyraMode()*/));

      // call constructor
      BlockedCrsMatrix::initialize(rgMapExtractor,
                       doMapExtractor,
                       33);
#endif

      // blocks must be filled externally...
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
      // TODO
      /*
      using Teuchos::RCP;

      TEUCHOS_TEST_FOR_EXCEPTION(mode != Teuchos::NO_TRANS && mode != Teuchos::TRANS, Xpetra::Exceptions::RuntimeError,
                                 "apply() only supports the following modes: NO_TRANS and TRANS." );

      RCP<const MultiVector> refX = rcpFromRef(X);
      RCP<MultiVector>       tmpY = MultiVectorFactory::Build(Y.getMap(), Y.getNumVectors());

      SC one = ScalarTraits<SC>::one();

      if (mode == Teuchos::NO_TRANS) {
        for (size_t row = 0; row < Rows(); row++) {
          RCP<MultiVector>    Yblock = rangemaps_->getVector(row, Y.getNumVectors(), bRangeThyraMode_);
          RCP<MultiVector> tmpYblock = rangemaps_->getVector(row, Y.getNumVectors(), bRangeThyraMode_);

          for (size_t col = 0; col < Cols(); col++) {
            RCP<const MultiVector> Xblock = domainmaps_->ExtractVector(refX, col, bDomainThyraMode_);
            RCP<Matrix>            Ablock = getMatrix(row, col);

            if (Ablock.is_null())
              continue;

            Ablock->apply(*Xblock, *tmpYblock);
            Yblock->update(one, *tmpYblock, one);
          }
          rangemaps_->InsertVector(Yblock, row, tmpY, bRangeThyraMode_);
        }

      } else if (mode == Teuchos::TRANS) {
        // TODO: test me!
        for (size_t col = 0; col < Cols(); col++) {
          RCP<MultiVector>    Yblock = domainmaps_->getVector(col, Y.getNumVectors(), bDomainThyraMode_);
          RCP<MultiVector> tmpYblock = domainmaps_->getVector(col, Y.getNumVectors(), bDomainThyraMode_);

          for (size_t row = 0; row<Rows(); row++) {
            RCP<const MultiVector> Xblock = rangemaps_->ExtractVector(refX, row, bRangeThyraMode_);
            RCP<Matrix>            Ablock = getMatrix(row, col);

            if (Ablock.is_null())
              continue;

            Ablock->apply(*Xblock, *tmpYblock, Teuchos::TRANS);

            Yblock->update(one, *tmpYblock, one);
          }
          domainmaps_->InsertVector(Yblock, col, tmpY, bDomainThyraMode_);
        }
      }

      Y.update(alpha, *tmpY, beta);
      */
    }

    //! \brief Returns the Map associated with the full domain of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    /*RCP<const Map > getDomainMap() const            { return domainmaps_->getFullMap(); }

    //! \brief Returns the Map associated with the i'th block domain of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    RCP<const Map > getDomainMap(size_t i, bool bThyraMode = false) const    { return domainmaps_->getMap(i, bDomainThyraMode_); }

    //! Returns the Map associated with the full range of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    RCP<const Map > getRangeMap() const             { return rangemaps_->getFullMap(); }

    //! Returns the Map associated with the i'th block range of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    RCP<const Map > getRangeMap(size_t i, bool bThyraMode = false) const     { return rangemaps_->getMap(i, bRangeThyraMode_); }

    //! Returns map extractor class for range map
    RCP<const MapExtractor> getRangeMapExtractor() const { return rangemaps_; }

    //! Returns map extractor for domain map
    RCP<const MapExtractor> getDomainMapExtractor() const { return domainmaps_; }*/

    //@}

    //! Implements DistObject interface
    //{@

    //! Access function for the Tpetra::Map this DistObject was constructed with.
    /*const Teuchos::RCP< const Map > getMap() const {
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::getMap(): operation not supported.");
    }

    //! Import.
    void doImport(const Matrix &source, const Import& importer, CombineMode CM) {
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doImport(): operation not supported.");
    }

    //! Export.
    void doExport(const Matrix& dest, const Import& importer, CombineMode CM) {
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doExport(): operation not supported.");
    }

    //! Import (using an Exporter).
    void doImport(const Matrix& source, const Export& exporter, CombineMode CM) {
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doImport(): operation not supported.");
    }

    //! Export (using an Importer).
    void doExport(const Matrix& dest, const Export& exporter, CombineMode CM) {
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doExport(): operation not supported.");
    }*/

    // @}

    //! @name Overridden from Teuchos::Describable
    //@{

    /** \brief Return a simple one-line description of this object. */
    std::string description() const                       { return "ReorderedBlockedCrsMatrix"; }

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {
      out << "Xpetra::ReorderedBlockedCrsMatrix: " << std::endl; // << Rows() << " x " << Cols() << std::endl;

      /*if (isFillComplete()) {
        out << "ReorderedBlockMatrix is fillComplete" << std::endl;

        //out << "fullRowMap" << std::endl;
        //fullrowmap_->describe(out,verbLevel);

        //out << "fullColMap" << std::endl;
        //fullcolmap_->describe(out,verbLevel);

      } else {
        out << "BlockMatrix is NOT fillComplete" << std::endl;
      }*/

      /*for (size_t r = 0; r < Rows(); ++r)
        for (size_t c = 0; c < Cols(); ++c) {
          out << "Block(" << r << "," << c << ")" << std::endl;
          getMatrix(r,c)->describe(out,verbLevel);
        }*/
    }

    //@}

  private:
    Teuchos::RCP<const Xpetra::BlockReorderManager > brm_;
    Teuchos::RCP<const Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > fullOp_;


};

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > mergeSubBlockMaps(Teuchos::RCP<const Xpetra::BlockReorderManager> brm, Teuchos::RCP<const Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > bmat) {
    typedef Xpetra::MapUtils<LocalOrdinal,GlobalOrdinal,Node> MapUtils;

    // TODO distinguish between range and domain map extractor! provide MapExtractor as parameter!
    RCP<const Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node> > fullRangeMapExtractor = bmat->getRangeMapExtractor();

    // number of sub blocks
    size_t numBlocks = brm->GetNumBlocks();

    Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node>> map = Teuchos::null;

    if(numBlocks == 0) {
      // it is a leaf node
      Teuchos::RCP<const Xpetra::BlockReorderLeaf> leaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(brm);

      // never extract Thyra style maps (since we have to merge them)
      map = fullRangeMapExtractor->getMap(Teuchos::as<size_t>(leaf->GetIndex()), false);
    } else {
      // initialize vector for sub maps
      std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > > subMaps (numBlocks, Teuchos::null);

      for(size_t i = 0; i < numBlocks; i++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> blkMgr = brm->GetBlock(Teuchos::as<int>(i));
        subMaps[i] = mergeSubBlockMaps(blkMgr,bmat);
        TEUCHOS_ASSERT(subMaps[i].is_null()==false);
      }

      map = MapUtils::concatenateMaps(subMaps);
    }
    TEUCHOS_ASSERT(map.is_null()==false);
    return map;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mergeSubBlocks(Teuchos::RCP<const Xpetra::BlockReorderManager> rowMgr, Teuchos::RCP<const Xpetra::BlockReorderManager> colMgr, Teuchos::RCP<const Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > bmat) {
  // number of sub blocks
  size_t rowSz = rowMgr->GetNumBlocks();
  size_t colSz = colMgr->GetNumBlocks();

  if(rowSz == 0 && colSz == 0) {
    // it is a leaf node
    Teuchos::RCP<const Xpetra::BlockReorderLeaf> rowleaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(rowMgr);
    Teuchos::RCP<const Xpetra::BlockReorderLeaf> colleaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(colMgr);

    Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mat = bmat->getMatrix(rowleaf->GetIndex(), colleaf->GetIndex());

    return mat;
  } else {
    // create the map extractors
    // we cannot create block matrix in thyra mode since merged maps might not start with 0 GID
    Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rgMapExtractor = Teuchos::null;
    if(rowSz > 0) {
      std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > > rowSubMaps (rowSz, Teuchos::null);
      for(size_t i = 0; i < rowSz; i++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
        rowSubMaps[i] = mergeSubBlockMaps(rowSubMgr,bmat);
        TEUCHOS_ASSERT(rowSubMaps[i].is_null()==false);
      }
      Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rgMergedSubMaps = Xpetra::MapUtils<LocalOrdinal,GlobalOrdinal,Node>::concatenateMaps(rowSubMaps);
      rgMapExtractor = Teuchos::rcp(new Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rgMergedSubMaps, rowSubMaps, false /*fullRangeMapExtractor->getThyraMode()*/));
    } else {
      Teuchos::RCP<const Xpetra::BlockReorderLeaf> rowleaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(rowMgr);
      RCP<const Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node> > fullRangeMapExtractor = bmat->getRangeMapExtractor();
      // TODO think about Thyra style maps: we cannot use thyra style maps when recombining several blocks!!!
      // The GIDs might not start with 0 and may not be consecutive!
      Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > submap = fullRangeMapExtractor->getMap(rowleaf->GetIndex(),false);
      std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > > rowSubMaps (1, submap);
      rgMapExtractor = Teuchos::rcp(new Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node>(submap, rowSubMaps, false /*fullRangeMapExtractor->getThyraMode()*/));
    }

    Teuchos::RCP<const Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node> > doMapExtractor = Teuchos::null;
    if(colSz > 0) {
      std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > > colSubMaps (colSz, Teuchos::null);
      for(size_t j = 0; j < colSz; j++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> colSubMgr = colMgr->GetBlock(Teuchos::as<int>(j));
        colSubMaps[j] = mergeSubBlockMaps(colSubMgr,bmat);
        TEUCHOS_ASSERT(colSubMaps[j].is_null()==false);
      }
      Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > doMergedSubMaps = Xpetra::MapUtils<LocalOrdinal,GlobalOrdinal,Node>::concatenateMaps(colSubMaps);
      doMapExtractor = Teuchos::rcp(new Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node>(doMergedSubMaps, colSubMaps, false /*fullDomainMapExtractor->getThyraMode()*/));
    } else {
      Teuchos::RCP<const Xpetra::BlockReorderLeaf> colleaf = Teuchos::rcp_dynamic_cast<const Xpetra::BlockReorderLeaf>(colMgr);
      RCP<const Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node> > fullDomainMapExtractor = bmat->getDomainMapExtractor();
      // TODO think about Thyra style maps: we cannot use thyra style maps when recombining several blocks!!!
      // The GIDs might not start with 0 and may not be consecutive!
      Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > submap = fullDomainMapExtractor->getMap(colleaf->GetIndex(),false);
      std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > > colSubMaps (1, submap);
      doMapExtractor = Teuchos::rcp(new Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node>(submap, colSubMaps, false /*fullRangeMapExtractor->getThyraMode()*/));
    }

    Teuchos::RCP<Xpetra::ReorderedBlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rbmat =
        Teuchos::rcp(new Xpetra::ReorderedBlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(rgMapExtractor,doMapExtractor, 33, rowMgr,bmat));

    if (rowSz == 0 && colSz > 0) {
      for(size_t j = 0; j < colSz; j++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> colSubMgr = colMgr->GetBlock(Teuchos::as<int>(j));
        Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > submat =
            mergeSubBlocks(rowMgr, colSubMgr, bmat);
        rbmat->setMatrix(0,j,Teuchos::rcp_const_cast<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(submat));
      }
    } else if (rowSz > 0 && colSz == 0) {
      for(size_t i = 0; i < rowSz; i++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
        Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > submat =
            mergeSubBlocks(rowSubMgr, colMgr, bmat);
        rbmat->setMatrix(i,0,Teuchos::rcp_const_cast<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(submat));
      }
    } else {
      for(size_t i = 0; i < rowSz; i++) {
        Teuchos::RCP<const Xpetra::BlockReorderManager> rowSubMgr = rowMgr->GetBlock(Teuchos::as<int>(i));
        for(size_t j = 0; j < colSz; j++) {
          Teuchos::RCP<const Xpetra::BlockReorderManager> colSubMgr = colMgr->GetBlock(Teuchos::as<int>(j));
          Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > submat =
              mergeSubBlocks(rowSubMgr, colSubMgr, bmat);
          rbmat->setMatrix(i,j,Teuchos::rcp_const_cast<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(submat));
        }
      }
    }

    rbmat->fillComplete();
    return rbmat;
  }

  return Teuchos::null;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > buildReorderedBlockedCrsMatrix(Teuchos::RCP<const Xpetra::BlockReorderManager> brm, Teuchos::RCP<const Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > bmat) {
  Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rbmat =
      mergeSubBlocks(brm, brm, bmat);
  TEUCHOS_ASSERT(rbmat != Teuchos::null);
  return rbmat;
}

} //namespace Xpetra

#define XPETRA_REORDEREDBLOCKEDCRSMATRIX_SHORT
#endif /* XPETRA_REORDEREDBLOCKEDCRSMATRIX_HPP */
