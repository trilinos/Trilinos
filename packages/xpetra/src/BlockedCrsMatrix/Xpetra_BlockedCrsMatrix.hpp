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
#ifndef XPETRA_BLOCKEDCRSMATRIX_HPP
#define XPETRA_BLOCKEDCRSMATRIX_HPP

#include <KokkosCompat_DefaultNode.hpp>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Hashtable.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_BlockedMultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_BlockedVector.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_CrsMatrixFactory.hpp"

#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_MapExtractorFactory.hpp"

#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixFactory.hpp"
#include "Xpetra_CrsMatrixWrap.hpp"

#ifdef HAVE_XPETRA_THYRA
#include <Thyra_ProductVectorSpaceBase.hpp>
#include <Thyra_VectorSpaceBase.hpp>
#include <Thyra_LinearOpBase.hpp>
#include <Thyra_BlockedLinearOpBase.hpp>
#include <Thyra_PhysicallyBlockedLinearOpBase.hpp>
#include "Xpetra_ThyraUtils.hpp"
#endif

#include "Xpetra_VectorFactory.hpp"


/** \file Xpetra_BlockedCrsMatrix.hpp

  Declarations for the class Xpetra::BlockedCrsMatrix.
*/
namespace Xpetra {

#ifdef HAVE_XPETRA_THYRA
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node> class ThyraUtils;
#endif

  typedef std::string viewLabel_t;

  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class BlockedCrsMatrix :
    public Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  public:
    typedef Scalar scalar_type;
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;

  private:
#undef XPETRA_BLOCKEDCRSMATRIX_SHORT
#include "Xpetra_UseShortNames.hpp"

  public:

    //! @name Constructor/Destructor Methods
    //@{

    //! Constructor
    /*!
     * \param rangeMaps range maps for all blocks
     * \param domainMaps domain maps for all blocks
     * \param numEntriesPerRow estimated number of entries per row in each block(!)
     */
    BlockedCrsMatrix(const Teuchos::RCP<const BlockedMap>& rangeMaps,
                     const Teuchos::RCP<const BlockedMap>& domainMaps,
                     size_t numEntriesPerRow)
      : is_diagonal_(true)
    {
      domainmaps_ = MapExtractorFactory::Build(domainMaps);
      rangemaps_ = MapExtractorFactory::Build(rangeMaps);
      bRangeThyraMode_ = rangeMaps->getThyraMode();
      bDomainThyraMode_ = domainMaps->getThyraMode();

      blocks_.reserve(Rows()*Cols());

      // add CrsMatrix objects in row,column order
      for (size_t r = 0; r < Rows(); ++r)
        for (size_t c = 0; c < Cols(); ++c)
          blocks_.push_back(MatrixFactory::Build(getRangeMap(r, bRangeThyraMode_), numEntriesPerRow));

      // Default view
      CreateDefaultView();
    }

    //! Constructor
    /*!
     * \param rangeMapExtractor range map extractor for all blocks
     * \param domainMapExtractor domain map extractor for all blocks
     * \param numEntriesPerRow estimated number of entries per row in each block(!)
     *
     * \note This constructor will be deprecated. Please use the constructor which takes BlockedMap objects instead.
     */
    BlockedCrsMatrix(Teuchos::RCP<const MapExtractor>& rangeMapExtractor,
                     Teuchos::RCP<const MapExtractor>& domainMapExtractor,
                     size_t numEntriesPerRow)
      : is_diagonal_(true), domainmaps_(domainMapExtractor), rangemaps_(rangeMapExtractor)
    {
      bRangeThyraMode_ = rangeMapExtractor->getThyraMode();
      bDomainThyraMode_ = domainMapExtractor->getThyraMode();

      blocks_.reserve(Rows()*Cols());

      // add CrsMatrix objects in row,column order
      for (size_t r = 0; r < Rows(); ++r)
        for (size_t c = 0; c < Cols(); ++c)
          blocks_.push_back(MatrixFactory::Build(getRangeMap(r, bRangeThyraMode_), numEntriesPerRow));

      // Default view
      CreateDefaultView();
    }

#ifdef HAVE_XPETRA_THYRA
    //! Constructor
    /*!
     * \param rangeMaps range maps for all blocks
     * \param domainMaps domain maps for all blocks
     * \param npr extimated number of entries per row in each block(!)
     */
    BlockedCrsMatrix(const Teuchos::RCP<const Thyra::BlockedLinearOpBase<Scalar> >& thyraOp, const Teuchos::RCP<const Teuchos::Comm<int> >& /* comm */)
      : is_diagonal_(true), thyraOp_(thyraOp)
    {
      // extract information from Thyra blocked operator and rebuilt information
      const Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > productRangeSpace  = thyraOp->productRange();
      const Teuchos::RCP<const Thyra::ProductVectorSpaceBase<Scalar> > productDomainSpace = thyraOp->productDomain();

      int numRangeBlocks  = productRangeSpace->numBlocks();
      int numDomainBlocks = productDomainSpace->numBlocks();

      // build range map extractor from Thyra::BlockedLinearOpBase object
      std::vector<Teuchos::RCP<const Map> > subRangeMaps(numRangeBlocks);
      for (size_t r=0; r<Teuchos::as<size_t>(numRangeBlocks); ++r) {
        for (size_t c=0; c<Teuchos::as<size_t>(numDomainBlocks); ++c) {
          if (thyraOp->blockExists(r,c)) {
            // we only need at least one block in each block row to extract the range map
            Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > const_op = thyraOp->getBlock(r,c); // nonConst access is not allowed.
            Teuchos::RCP<const Xpetra::Matrix<Scalar,LO,GO,Node> > xop =
                            Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toXpetra(const_op);
            subRangeMaps[r] = xop->getRangeMap();
            if(r!=c) is_diagonal_ = false;
            break;
          }
        }
      }
      Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > fullRangeMap = mergeMaps(subRangeMaps);

      // check whether the underlying Thyra operator uses Thyra-style numbering (default for most applications) or
      // Xpetra style numbering
      bool bRangeUseThyraStyleNumbering = false;
      size_t numAllElements = 0;
      for(size_t v = 0; v < subRangeMaps.size(); ++v) {
        numAllElements += subRangeMaps[v]->getGlobalNumElements();
      }
      if ( fullRangeMap->getGlobalNumElements() != numAllElements) bRangeUseThyraStyleNumbering = true;
      rangemaps_ = Teuchos::rcp(new Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node>(fullRangeMap, subRangeMaps, bRangeUseThyraStyleNumbering));

      // build domain map extractor from Thyra::BlockedLinearOpBase object
      std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > > subDomainMaps(numDomainBlocks);
      for (size_t c=0; c<Teuchos::as<size_t>(numDomainBlocks); ++c) {
        for (size_t r=0; r<Teuchos::as<size_t>(numRangeBlocks); ++r) {
          if (thyraOp->blockExists(r,c)) {
            // we only need at least one block in each block row to extract the range map
            Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > const_op = thyraOp->getBlock(r,c); // nonConst access is not allowed.
            Teuchos::RCP<const Xpetra::Matrix<Scalar,LO,GO,Node> > xop =
                            Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toXpetra(const_op);
            subDomainMaps[c] = xop->getDomainMap();
            break;
          }
        }
      }
      Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > fullDomainMap = mergeMaps(subDomainMaps);
      // plausibility check for numbering style (Xpetra or Thyra)
      bool bDomainUseThyraStyleNumbering = false;
      numAllElements = 0;
      for(size_t v = 0; v < subDomainMaps.size(); ++v) {
        numAllElements += subDomainMaps[v]->getGlobalNumElements();
      }
      if (fullDomainMap->getGlobalNumElements() != numAllElements) bDomainUseThyraStyleNumbering = true;
      domainmaps_ = Teuchos::rcp(new Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node>(fullDomainMap, subDomainMaps, bDomainUseThyraStyleNumbering));

      // store numbering mode
      bRangeThyraMode_  = bRangeUseThyraStyleNumbering;
      bDomainThyraMode_ = bDomainUseThyraStyleNumbering;

      // add CrsMatrix objects in row,column order
      blocks_.reserve(Rows()*Cols());
      for (size_t r = 0; r < Rows(); ++r) {
        for (size_t c = 0; c < Cols(); ++c) {
          if(thyraOp->blockExists(r,c)) {
            // TODO we do not support nested Thyra operators here!
            Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > const_op = thyraOp->getBlock(r,c); // nonConst access is not allowed.
            Teuchos::RCP<Thyra::LinearOpBase<Scalar> > op = Teuchos::rcp_const_cast<Thyra::LinearOpBase<Scalar> >(const_op); // cast away const
            Teuchos::RCP<Xpetra::Matrix<Scalar,LO,GO,Node> > xop =
                Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toXpetra(op);
            Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node> > xwrap =
              Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<Scalar,LO,GO,Node> >(xop, true);
            blocks_.push_back(xwrap);
          } else {
            // add empty block
            blocks_.push_back(MatrixFactory::Build(getRangeMap(r,bRangeThyraMode_), 0));
          }
        }
      }
      // Default view
      CreateDefaultView();
    }

  private:
    //! mergeMaps
    /*!
     * \param subMaps
     *
     * Merges all Xpetra::Map objects in std::vector subMaps and returns a new Xpetra::Map containing all entries.
     * Helper function only used in constructor of Xpetra_BlockedCrsMatrix for transforming a Thyra::BlockedLinearOp object
     * All GID entries are sorted and duplicates are eliminated.
     */
    Teuchos::RCP<const Map> mergeMaps(std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > > & subMaps) {
      // TODO merging for Thyra mode is missing (similar to what we do in constructor of MapExtractor

      // merge submaps to global map
      std::vector<GlobalOrdinal> gids;
      for(size_t tt = 0; tt<subMaps.size(); ++tt) {
        Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > subMap = subMaps[tt];
#if 1
        Teuchos::ArrayView< const GlobalOrdinal > subMapGids = subMap->getLocalElementList();
        gids.insert(gids.end(), subMapGids.begin(), subMapGids.end());
#else
        size_t myNumElements = subMap->getLocalNumElements();
        for(LocalOrdinal l = 0; l < Teuchos::as<LocalOrdinal>(myNumElements); ++l) {
          GlobalOrdinal gid = subMap->getGlobalElement(l);
          gids.push_back(gid);
        }
#endif
      }

      // we have to sort the matrix entries and get rid of the double entries
      // since we use this to detect Thyra-style numbering or Xpetra-style
      // numbering. In Thyra-style numbering mode, the Xpetra::MapExtractor builds
      // the correct row maps.
      const GO INVALID = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();
      std::sort(gids.begin(), gids.end());
      gids.erase(std::unique(gids.begin(), gids.end()), gids.end());
      Teuchos::ArrayView<GO> gidsView(&gids[0], gids.size());
      Teuchos::RCP<Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > fullMap = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(subMaps[0]->lib(), INVALID, gidsView, subMaps[0]->getIndexBase(), subMaps[0]->getComm());
      return fullMap;
    }

  public:
#endif

    //! Destructor
    virtual ~BlockedCrsMatrix() {}

    //@}


    //! @name Insertion/Removal Methods
    //@{

    //! Insert matrix entries, using global IDs.
    /**
      Note: this routine throws for Rows() > 1 and/or Cols() > 1

      All index values must be in the global space.
      \pre \c globalRow exists as an ID in the global row map
      \pre <tt>isLocallyIndexed() == false</tt>
      \pre <tt>isStorageOptimized() == false</tt>

      \post <tt>isGloballyIndexed() == true</tt>

      \note If \c globalRow does not belong to the matrix on this node, then it
      will be communicated to the appropriate node when globalAssemble() is
      called (which will, at the latest, occur during the next call to
      fillComplete().) Otherwise, the entries will be inserted in the local
      matrix.

      \note If the matrix row already contains values at the indices
      corresponding to values in \c cols, then the new values will be summed
      with the old values; this may happen at insertion or during the next call
      to fillComplete().

      \note If <tt>hasColMap() == true</tt>, only (cols[i],vals[i]) where
      cols[i] belongs to the column map on this node will be inserted into the
      matrix.
      */
    void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal>& cols, const ArrayView<const Scalar>& vals) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::insertGlobalValues");
      if (Rows() == 1 && Cols () == 1) {
        getMatrix(0,0)->insertGlobalValues(globalRow, cols, vals);
        return;
      }
      throw Xpetra::Exceptions::RuntimeError("insertGlobalValues not supported by BlockedCrsMatrix");
    }

    //! Insert matrix entries, using local IDs.
    /**
       Note: this routine throws if Rows() > 1 and/or Cols() > 1

       All index values must be in the local space.
      \pre \c localRow exists as an ID in the local row map
      \pre <tt>isGloballyIndexed() == false</tt>
      \pre <tt>isStorageOptimized() == false</tt>

      \post <tt>isLocallyIndexed() == true</tt>
      */
    void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal>& cols, const ArrayView<const Scalar>& vals) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::insertLocalValues");
      if (Rows() == 1 && Cols () == 1) {
        getMatrix(0,0)->insertLocalValues(localRow, cols, vals);
        return;
      }
      throw Xpetra::Exceptions::RuntimeError("insertLocalValues not supported by BlockedCrsMatrix");
    }

    void removeEmptyProcessesInPlace(const Teuchos::RCP<const Map>& newMap) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::removeEmptyProcessesInPlace");
      if (Rows() == 1 && Cols () == 1) {
        getMatrix(0,0)->removeEmptyProcessesInPlace(newMap);
        return;
      }
      throw Xpetra::Exceptions::RuntimeError("removeEmptyProcesses not supported by BlockedCrsMatrix");
    }

    //! \brief Replace matrix entries, using global IDs.
    /** All index values must be in the global space.

      \pre \c globalRow is a global row belonging to the matrix on this node.

      \note If (globalRow,cols[i]) corresponds to an entry that is duplicated in
      this matrix row (likely because it was inserted more than once and
      fillComplete() has not been called in the interim), the behavior of this
      function is not defined. */
    void replaceGlobalValues(GlobalOrdinal globalRow,
                             const ArrayView<const GlobalOrdinal> &cols,
                             const ArrayView<const Scalar>        &vals) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::replaceGlobalValues");
      if (Rows() == 1 && Cols () == 1) {
        getMatrix(0,0)->replaceGlobalValues(globalRow,cols,vals);
        return;
      }
      throw Xpetra::Exceptions::RuntimeError("replaceGlobalValues not supported by BlockedCrsMatrix");
    }

    //! Replace matrix entries, using local IDs.
    /** All index values must be in the local space.
      Note that if a value is not already present for the specified location in
      the matrix, the input value will be ignored silently.
      */
    void replaceLocalValues(LocalOrdinal localRow,
                            const ArrayView<const LocalOrdinal> &cols,
                            const ArrayView<const Scalar>       &vals) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::replaceLocalValues");
      if (Rows() == 1 && Cols () == 1) {
        getMatrix(0,0)->replaceLocalValues(localRow,cols,vals);
        return;
      }
      throw Xpetra::Exceptions::RuntimeError("replaceLocalValues not supported by BlockedCrsMatrix");
    }

    //! Set all matrix entries equal to scalar
    virtual void setAllToScalar(const Scalar& alpha) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::setAllToScalar");
      for (size_t row = 0; row < Rows(); row++) {
        for (size_t col = 0; col < Cols(); col++) {
          if (!getMatrix(row,col).is_null()) {
            getMatrix(row,col)->setAllToScalar(alpha);
          }
        }
      }
    }

    //! Scale the current values of a matrix, this = alpha*this.
    void scale(const Scalar& alpha) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::scale");
      for (size_t row = 0; row < Rows(); row++) {
        for (size_t col = 0; col < Cols(); col++) {
          if (!getMatrix(row,col).is_null()) {
            getMatrix(row,col)->scale(alpha);
          }
        }
      }
    }

    //@}

    //! @name Transformational Methods
    //@{

    /*! Resume fill operations.
      After calling fillComplete(), resumeFill() must be called before initiating any changes to the matrix.

      For BlockedCrsMatrix objects we call the routine iteratively for all sub-blocks.

      resumeFill() may be called repeatedly.
      */
    void resumeFill(const RCP< ParameterList >& params = null) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::resumeFill");
      for (size_t row = 0; row < Rows(); row++) {
        for (size_t col = 0; col < Cols(); col++) {
          if (!getMatrix(row,col).is_null()) {
            getMatrix(row,col)->resumeFill(params);
          }
        }
      }
    }

    /*! \brief Signal that data entry is complete.

      Note: for blocked operators the specified domain and range maps have no meaning.
            We just call fillComplete for all underlying blocks
      */
    void fillComplete(const RCP<const Map>& domainMap, const RCP<const Map>& rangeMap, const RCP<ParameterList>& params = null) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::fillComplete");
      if (Rows() == 1 && Cols () == 1) {
        getMatrix(0,0)->fillComplete(domainMap, rangeMap, params);
        return;
      }
      fillComplete(params);
    }

    /*! \brief Signal that data entry is complete.

      Off-node entries are distributed (via globalAssemble()), repeated entries are
      summed, and global indices are transformed to local indices.

      \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

      \pre  <tt>isFillActive() == true<tt>
      \pre <tt>isFillComplete()() == false<tt>

      \post <tt>isFillActive() == false<tt>
      \post <tt>isFillComplete() == true<tt>
      \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
      */
    void fillComplete(const RCP<ParameterList>& params = null) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::fillComplete");
      TEUCHOS_TEST_FOR_EXCEPTION(rangemaps_==Teuchos::null, Xpetra::Exceptions::RuntimeError,"BlockedCrsMatrix::fillComplete: rangemaps_ is not set. Error.");

      for (size_t r = 0; r < Rows(); ++r)
        for (size_t c = 0; c < Cols(); ++c) {
          if(getMatrix(r,c) != Teuchos::null) {
            Teuchos::RCP<Matrix> Ablock = getMatrix(r,c);
            if(r!=c) is_diagonal_ = false;
            if (Ablock != Teuchos::null && !Ablock->isFillComplete())
              Ablock->fillComplete(getDomainMap(c, bDomainThyraMode_), getRangeMap(r, bRangeThyraMode_), params);
          }
        }

#if 0
      // get full row map
      RCP<const Map> rangeMap = rangemaps_->getFullMap();
      fullrowmap_ = MapFactory::Build(rangeMap()->lib(), rangeMap()->getGlobalNumElements(), rangeMap()->getLocalElementList(), rangeMap()->getIndexBase(), rangeMap()->getComm());

      // build full col map
      fullcolmap_ = Teuchos::null; // delete old full column map

      std::vector<GO> colmapentries;
      for (size_t c = 0; c < Cols(); ++c) {
        // copy all local column lids of all block rows to colset
        std::set<GO> colset;
        for (size_t r = 0; r < Rows(); ++r) {
          Teuchos::RCP<CrsMatrix> Ablock = getMatrix(r,c);

          if (Ablock != Teuchos::null) {
            Teuchos::ArrayView<const GO> colElements = Ablock->getColMap()->getLocalElementList();
            Teuchos::RCP<const Map>      colmap      = Ablock->getColMap();
            copy(colElements.getRawPtr(), colElements.getRawPtr() + colElements.size(), inserter(colset, colset.begin()));
          }
        }

        // remove duplicates (entries which are in column maps of more than one block row)
        colmapentries.reserve(colmapentries.size() + colset.size());
        copy(colset.begin(), colset.end(), back_inserter(colmapentries));
        sort(colmapentries.begin(), colmapentries.end());
        typename std::vector<GO>::iterator gendLocation;
        gendLocation = std::unique(colmapentries.begin(), colmapentries.end());
        colmapentries.erase(gendLocation,colmapentries.end());
      }

      // sum up number of local elements
      size_t numGlobalElements = 0;
      Teuchos::reduceAll(*(rangeMap->getComm()), Teuchos::REDUCE_SUM, colmapentries.size(), Teuchos::outArg(numGlobalElements));

      // store global full column map
      const Teuchos::ArrayView<const GO> aView = Teuchos::ArrayView<const GO>(colmapentries);
      fullcolmap_ = MapFactory::Build(rangeMap->lib(), numGlobalElements, aView, 0, rangeMap->getComm());
#endif
    }

    //@}

    //! Returns the number of global rows.
    /** Undefined if isFillActive().
    */
    global_size_t getGlobalNumRows() const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getGlobalNumRows");
      global_size_t globalNumRows = 0;

      for (size_t row = 0; row < Rows(); row++)
        for (size_t col = 0; col < Cols(); col++)
          if (!getMatrix(row,col).is_null()) {
            globalNumRows += getMatrix(row,col)->getGlobalNumRows();
            break; // we need only one non-null matrix in a row
          }

      return globalNumRows;
    }

    //! \brief Returns the number of global columns in the matrix.
    /** Undefined if isFillActive().
    */
    global_size_t getGlobalNumCols() const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getGlobalNumCols");
      global_size_t globalNumCols = 0;

      for (size_t col = 0; col < Cols(); col++)
        for (size_t row = 0; row < Rows(); row++)
          if (!getMatrix(row,col).is_null()) {
            globalNumCols += getMatrix(row,col)->getGlobalNumCols();
            break; // we need only one non-null matrix in a col
          }

      return globalNumCols;
    }

    //! Returns the number of matrix rows owned on the calling node.
    size_t getLocalNumRows() const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getLocalNumRows");
      global_size_t nodeNumRows = 0;

      for (size_t row = 0; row < Rows(); ++row)
        for (size_t col = 0; col < Cols(); col++)
          if (!getMatrix(row,col).is_null()) {
            nodeNumRows += getMatrix(row,col)->getLocalNumRows();
            break; // we need only one non-null matrix in a row
          }

      return nodeNumRows;
    }

    //! Returns the global number of entries in this matrix.
    global_size_t getGlobalNumEntries() const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getGlobalNumEntries");
      global_size_t globalNumEntries = 0;

      for (size_t row = 0; row < Rows(); ++row)
        for (size_t col = 0; col < Cols(); ++col)
          if (!getMatrix(row,col).is_null())
            globalNumEntries += getMatrix(row,col)->getGlobalNumEntries();

      return globalNumEntries;
    }

    //! Returns the local number of entries in this matrix.
    size_t getLocalNumEntries() const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getLocalNumEntries");
      global_size_t nodeNumEntries = 0;

      for (size_t row = 0; row < Rows(); ++row)
        for (size_t col = 0; col < Cols(); ++col)
          if (!getMatrix(row,col).is_null())
            nodeNumEntries += getMatrix(row,col)->getLocalNumEntries();

      return nodeNumEntries;
    }

    //! Returns the current number of entries on this node in the specified local row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
    size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getNumEntriesInLocalRow");
      GlobalOrdinal gid = this->getRowMap()->getGlobalElement(localRow);
      size_t row = getBlockedRangeMap()->getMapIndexForGID(gid);
      LocalOrdinal lid = getBlockedRangeMap()->getMap(row)->getLocalElement(gid);
      size_t numEntriesInLocalRow = 0;
      for (size_t col = 0; col < Cols(); ++col)
        if (!getMatrix(row,col).is_null())
          numEntriesInLocalRow += getMatrix(row,col)->getNumEntriesInLocalRow(lid);
      return numEntriesInLocalRow;
    }

    //! Returns the current number of entries in the specified (locally owned) global row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
    size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getNumEntriesInGlobalRow");
      size_t row = getBlockedRangeMap()->getMapIndexForGID(globalRow);
      size_t numEntriesInGlobalRow = 0;
      for (size_t col = 0; col < Cols(); ++col)
        if (!getMatrix(row,col).is_null())
          numEntriesInGlobalRow += getMatrix(row,col)->getNumEntriesInGlobalRow(globalRow);
      return numEntriesInGlobalRow;
    }

    //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
    /** Undefined if isFillActive().
    */
    size_t getGlobalMaxNumRowEntries() const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getGlobalMaxNumRowEntries");
      global_size_t globalMaxEntries = 0;

      for (size_t row = 0; row < Rows(); row++) {
        global_size_t globalMaxEntriesBlockRows = 0;
        for (size_t col = 0; col < Cols(); col++) {
          if (!getMatrix(row,col).is_null()) {
            globalMaxEntriesBlockRows += getMatrix(row,col)->getGlobalMaxNumRowEntries();
          }
        }
        if(globalMaxEntriesBlockRows > globalMaxEntries)
          globalMaxEntries = globalMaxEntriesBlockRows;
      }
      return globalMaxEntries;
    }

    //! \brief Returns the maximum number of entries across all rows/columns on this node.
    /** Undefined if isFillActive().
    */
    size_t getLocalMaxNumRowEntries() const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getLocalMaxNumRowEntries");
      size_t localMaxEntries = 0;

      for (size_t row = 0; row < Rows(); row++) {
        size_t localMaxEntriesBlockRows = 0;
        for (size_t col = 0; col < Cols(); col++) {
          if (!getMatrix(row,col).is_null()) {
            localMaxEntriesBlockRows += getMatrix(row,col)->getLocalMaxNumRowEntries();
          }
        }
        if(localMaxEntriesBlockRows > localMaxEntries)
          localMaxEntries = localMaxEntriesBlockRows;
      }
      return localMaxEntries;
    }

    //! \brief If matrix indices of all matrix blocks are in the local range, this function returns true. Otherwise, this function returns false.
    /** if false, then this does not automatically mean that all blocks are globally indexed. The user has to make sure, that all matrix blocks
     * are indexed in the same way (locally or globally). Otherwise the block matrix is not valid...
     */
    bool isLocallyIndexed() const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::isLocallyIndexed");
      for (size_t i = 0; i < blocks_.size(); ++i)
        if (blocks_[i] != Teuchos::null && !blocks_[i]->isLocallyIndexed())
          return false;
      return true;
    }

    //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false.
    /** if false, then this does not automatically mean that all blocks are locally indexed. The user has to make sure, that all matrix blocks
     * are indexed in the same way (locally or globally). Otherwise the block matrix is not valid...
     */
    bool isGloballyIndexed() const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::isGloballyIndexed");
      for (size_t i = 0; i < blocks_.size(); i++)
        if (blocks_[i] != Teuchos::null && !blocks_[i]->isGloballyIndexed())
          return false;
      return true;
    }

    //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
    bool isFillComplete() const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::isFillComplete");
      for (size_t i = 0; i < blocks_.size(); i++)
        if (blocks_[i] != Teuchos::null && !blocks_[i]->isFillComplete())
          return false;
      return true;
    }

    //! Extract a list of entries in a specified local row of the matrix. Put into storage allocated by calling routine.
    /*!
      \param LocalRow - (In) Local row number for which indices are desired.
      \param Indices - (Out) Local column indices corresponding to values.
      \param Values - (Out) Matrix values.
      \param NumIndices - (Out) Number of indices.

      Note: A std::runtime_error exception is thrown if either \c Indices or \c
      Values is not large enough to hold the data associated with row \c
      LocalRow. If \c LocalRow is not valid for this node, then \c Indices and
      \c Values are unchanged and \c NumIndices is returned as
      OrdinalTraits<size_t>::invalid().

      \pre <tt>isLocallyIndexed()==true</tt> or <tt>hasColMap() == true</tt>
      */
    virtual void getLocalRowCopy(LocalOrdinal LocalRow,
                                 const ArrayView<LocalOrdinal>& Indices,
                                 const ArrayView<Scalar>& Values,
                                 size_t &NumEntries) const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getLocalRowCopy");
      if (Rows() == 1 && Cols () == 1) {
        getMatrix(0,0)->getLocalRowCopy(LocalRow, Indices, Values, NumEntries);
        return;
      }
      throw Xpetra::Exceptions::RuntimeError("getLocalRowCopy not supported by BlockedCrsMatrix" );
    }

    //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
    /*!
      \param GlobalRow - (In) Global row number for which indices are desired.
      \param Indices   - (Out) Global column indices corresponding to values.
      \param Values    - (Out) Row values
      \pre <tt>isLocallyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

      Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
      */
    void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal>& indices, ArrayView<const Scalar>& values) const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getGlobalRowView");
      if (Rows() == 1 && Cols () == 1) {
        getMatrix(0,0)->getGlobalRowView(GlobalRow, indices, values);
        return;
      }
      throw Xpetra::Exceptions::RuntimeError("getGlobalRowView not supported by BlockedCrsMatrix");
    }

    //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
    /*!
      \param LocalRow - (In) Local row number for which indices are desired.
      \param Indices  - (Out) Global column indices corresponding to values.
      \param Values   - (Out) Row values
      \pre <tt>isGloballyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

      Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
      */
    void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal>& indices, ArrayView<const Scalar>& values) const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getLocalRowView");
      if (Rows() == 1 && Cols () == 1) {
        getMatrix(0,0)->getLocalRowView(LocalRow, indices, values);
        return;
      }
      else if(is_diagonal_) {
        GlobalOrdinal gid = this->getRowMap()->getGlobalElement(LocalRow);
        size_t row = getBlockedRangeMap()->getMapIndexForGID(gid);
        getMatrix(row,row)->getLocalRowView(getMatrix(row,row)->getRowMap()->getLocalElement(gid),indices,values);
        return;
      }
      throw Xpetra::Exceptions::RuntimeError("getLocalRowView not supported by BlockedCrsMatrix");
    }

    //! \brief Get a copy of the diagonal entries owned by this node, with local row indices.
    /*! Returns a distributed Vector object partitioned according to this
      matrix's row map, containing the
      the zero and non-zero diagonals owned by this node. */
    void getLocalDiagCopy(Vector& diag) const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getLocalDiagCopy");

      //RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

      Teuchos::RCP<Vector> rcpdiag = Teuchos::rcpFromRef(diag);
      Teuchos::RCP<BlockedVector> bdiag = Teuchos::rcp_dynamic_cast<BlockedVector>(rcpdiag);

      // special treatment for 1x1 block matrices
      // ReorderedBlockedCrsMatrix object encapsulate single blocks in ReorderedBlocks
      // BlockedVectors have Vector objects as Leaf objects.
      if(Rows() == 1 && Cols() == 1 && bdiag.is_null() == true) {
        Teuchos::RCP<const Matrix> rm = getMatrix(0,0);
        rm->getLocalDiagCopy(diag);
        return;
      }

      TEUCHOS_TEST_FOR_EXCEPTION(bdiag.is_null() == true, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::getLocalDiagCopy: diag must be a Blocked(Multi)Vector.");
      TEUCHOS_TEST_FOR_EXCEPTION(bdiag->getNumVectors() != 1, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::getLocalDiagCopy: diag must be a Blocked(Multi)Vector with exactly one vector. However, the number of stored vectors is " << bdiag->getNumVectors());
      TEUCHOS_TEST_FOR_EXCEPTION(bdiag->getBlockedMap()->getNumMaps() != this->Rows(), Xpetra::Exceptions::RuntimeError,
        "BlockedCrsMatrix::getLocalDiagCopy(): the number of blocks in diag differ from the number of blocks in this operator." );
      //XPETRA_TEST_FOR_EXCEPTION(bdiag->getMap()->isSameAs(*(getMap())) == false, Xpetra::Exceptions::RuntimeError,
      //  "BlockedCrsMatrix::getLocalDiagCopy(): the map of the vector diag is not compatible with the map of the blocked operator." );

      for (size_t row = 0; row < Rows(); row++) {
        Teuchos::RCP<const Matrix> rm = getMatrix(row,row);
        if (!rm.is_null()) {
          Teuchos::RCP<Vector> rv = VectorFactory::Build(bdiag->getBlockedMap()->getMap(row,bdiag->getBlockedMap()->getThyraMode()));
          rm->getLocalDiagCopy(*rv);
          bdiag->setMultiVector(row,rv,bdiag->getBlockedMap()->getThyraMode());
        }
      }
    }

    //! Left scale matrix using the given vector entries
    void leftScale (const Vector& x) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::leftScale");

      Teuchos::RCP<const Vector> rcpx = Teuchos::rcpFromRef(x);
      Teuchos::RCP<const BlockedVector> bx = Teuchos::rcp_dynamic_cast<const BlockedVector>(rcpx);

      //RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

      // special treatment for 1xn block matrices
      // ReorderedBlockedCrsMatrix object encapsulate single blocks in ReorderedBlocks
      // BlockedVectors have Vector objects as Leaf objects.
      if(Rows() == 1 && bx.is_null() == true) {
        for (size_t col = 0; col < Cols(); ++col) {
          Teuchos::RCP<Matrix> rm = getMatrix(0,col);
          rm->leftScale(x);
        }
        return;
      }

      TEUCHOS_TEST_FOR_EXCEPTION(bx.is_null() == true, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::leftScale: x must be a Blocked(Multi)Vector.");
      TEUCHOS_TEST_FOR_EXCEPTION(bx->getNumVectors() != 1, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::leftScale: x must be a Blocked(Multi)Vector with exactly one vector. However, the number of stored vectors is " << bx->getNumVectors());
      TEUCHOS_TEST_FOR_EXCEPTION(bx->getBlockedMap()->getNumMaps() != this->Rows(), Xpetra::Exceptions::RuntimeError,
        "BlockedCrsMatrix::leftScale(): the number of blocks in diag differ from the number of blocks in this operator." );

      for (size_t row = 0; row < Rows(); row++) {
        Teuchos::RCP<const MultiVector> rmv = bx->getMultiVector(row);
        Teuchos::RCP<const Vector> rscale = rmv->getVector(0);
        XPETRA_TEST_FOR_EXCEPTION(rscale.is_null()==true, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::leftScale: x must be a Vector.");
        for (size_t col = 0; col < Cols(); ++col) {
          Teuchos::RCP<Matrix> rm = getMatrix(row,col);
          if (!rm.is_null()) {
            rm->leftScale(*rscale);
          }
        }
      }
    }

    //! Right scale matrix using the given vector entries
    void rightScale (const Vector& x) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::rightScale");

      Teuchos::RCP<const Vector> rcpx = Teuchos::rcpFromRef(x);
      Teuchos::RCP<const BlockedVector> bx = Teuchos::rcp_dynamic_cast<const BlockedVector>(rcpx);

      //RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

      // special treatment for nx1 block matrices
      // ReorderedBlockedCrsMatrix object encapsulate single blocks in ReorderedBlocks
      // BlockedVectors have Vector objects as Leaf objects.
      if(Cols() == 1 && bx.is_null() == true) {
        for (size_t row = 0; row < Rows(); ++row) {
          Teuchos::RCP<Matrix> rm = getMatrix(row,0);
          rm->rightScale(x);
        }
        return;
      }

      TEUCHOS_TEST_FOR_EXCEPTION(bx.is_null() == true, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::rightScale: x must be a Blocked(Multi)Vector.");
      TEUCHOS_TEST_FOR_EXCEPTION(bx->getNumVectors() != 1, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::rightScale: x must be a Blocked(Multi)Vector with exactly one vector. However, the number of stored vectors is " << bx->getNumVectors());
      TEUCHOS_TEST_FOR_EXCEPTION(bx->getBlockedMap()->getNumMaps() != this->Cols(), Xpetra::Exceptions::RuntimeError,
        "BlockedCrsMatrix::rightScale(): the number of blocks in diag differ from the number of blocks in this operator." );

      for (size_t col = 0; col < Cols(); ++col) {
        Teuchos::RCP<const MultiVector> rmv = bx->getMultiVector(col);
        Teuchos::RCP<const Vector> rscale = rmv->getVector(0);
        XPETRA_TEST_FOR_EXCEPTION(rscale.is_null()==true, Xpetra::Exceptions::RuntimeError, "BlockedCrsMatrix::leftScale: x must be a Vector.");
        for (size_t row = 0; row < Rows(); row++) {
          Teuchos::RCP<Matrix> rm = getMatrix(row,col);
          if (!rm.is_null()) {
            rm->rightScale(*rscale);
          }
        }
      }
    }


    //! Get Frobenius norm of the matrix
    virtual typename ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getFrobeniusNorm");
      typename ScalarTraits<Scalar>::magnitudeType ret = Teuchos::ScalarTraits<Scalar>::magnitude(Teuchos::ScalarTraits<Scalar>::zero());
      for (size_t col = 0; col < Cols(); ++col) {
        for (size_t row = 0; row < Rows(); ++row) {
          if(getMatrix(row,col)!=Teuchos::null) {
            typename ScalarTraits<Scalar>::magnitudeType n = getMatrix(row,col)->getFrobeniusNorm();
            ret += n * n;
          }
        }
      }
      return Teuchos::ScalarTraits< typename ScalarTraits<Scalar>::magnitudeType >::squareroot(ret);
    }


    //! Returns true if globalConstants have been computed; false otherwise
    virtual bool haveGlobalConstants() const {return true;}


    //@}

    //! @name Advanced Matrix-vector multiplication and solve methods
    //@{

    //! Multiplies this matrix by a MultiVector.
    /*! \c X is required to be post-imported, i.e., described by the column map
     * of the matrix. \c Y is required to be pre-exported, i.e., described by the
     * row map of the matrix.

     Both are required to have constant stride, and they are not permitted to
     ocupy overlapping space. No runtime checking will be performed in a
     non-debug build.

     This method is templated on the scalar type of MultiVector objects, allowing
     this method to be applied to MultiVector objects of arbitrary type. However,
     it is recommended that multiply() not be called directly; instead, use the
     CrsMatrixMultiplyOp, as it will handle the import/exprt operations required
     to apply a matrix with non-trivial communication needs.

     If \c beta is equal to zero, the operation will enjoy overwrite semantics
     (\c Y will be overwritten with the result of the multiplication). Otherwise,
     the result of the multiplication will be accumulated into \c Y.
     */

    //@}

    //! @name Methods implementing Matrix
    //@{

    //! sparse matrix-multivector multiplication for the region layout matrices (currently no blocked implementation)
    virtual void apply(const MultiVector &X, MultiVector &Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta, bool sumInterfaceValues,
                 const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> >& regionInterfaceImporter,
                 const Teuchos::ArrayRCP<LocalOrdinal>& regionInterfaceLIDs) const
    { }

    //! \brief Computes the sparse matrix-multivector multiplication.
    /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
      - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
      */
    virtual void apply(const MultiVector& X, MultiVector& Y,
                       Teuchos::ETransp mode = Teuchos::NO_TRANS,
                       Scalar alpha = ScalarTraits<Scalar>::one(),
                       Scalar beta  = ScalarTraits<Scalar>::zero()) const
    {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::apply");
      //using Teuchos::RCP;

      TEUCHOS_TEST_FOR_EXCEPTION(mode != Teuchos::NO_TRANS && mode != Teuchos::TRANS, Xpetra::Exceptions::RuntimeError,
                                 "apply() only supports the following modes: NO_TRANS and TRANS." );

      // check whether input parameters are blocked or not
      RCP<const MultiVector>         refX = rcpFromRef(X);
      RCP<const BlockedMultiVector> refbX = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(refX);
      //RCP<MultiVector>               tmpY = rcpFromRef(Y);
      //RCP<BlockedMultiVector>       tmpbY = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(tmpY);

      // TODO get rid of me: adapt MapExtractor
      bool bBlockedX = (refbX != Teuchos::null) ? true : false;

      // create (temporary) vectors for output
      // In the end we call Y.update(alpha, *tmpY, beta). Therefore we need a new vector storing the temporary results
      RCP<MultiVector> tmpY = MultiVectorFactory::Build(Y.getMap(), Y.getNumVectors(), true);

      //RCP<Teuchos::FancyOStream> out = rcp(new Teuchos::FancyOStream(rcp(&std::cout,false)));

      SC one = ScalarTraits<SC>::one();

      if (mode == Teuchos::NO_TRANS) {

        for (size_t row = 0; row < Rows(); row++) {
          RCP<MultiVector>    Yblock = rangemaps_->getVector(row, Y.getNumVectors(), bRangeThyraMode_, true);
          for (size_t col = 0; col < Cols(); col++) {

            // extract matrix block
            RCP<Matrix> Ablock = getMatrix(row, col);

            if (Ablock.is_null())
              continue;

            // check whether Ablock is itself a blocked operator
            // If it is a blocked operator we have to provide Xpetra style GIDs, i.e. we have to transform GIDs
            bool bBlockedSubMatrix = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(Ablock) == Teuchos::null ? false : true;

            // input/output vectors for local block operation
            RCP<const MultiVector> Xblock = Teuchos::null; // subpart of X vector to be applied to subblock of A

            // extract sub part of X using Xpetra or Thyra GIDs
            // if submatrix is again blocked, we extract it using Xpetra style gids. If it is a single
            // block matrix we use the Thyra or Xpetra style GIDs that are used to store the matrix
            if(bBlockedX) Xblock = domainmaps_->ExtractVector(refbX, col, bDomainThyraMode_);
            else          Xblock = domainmaps_->ExtractVector(refX,  col, bBlockedSubMatrix == true ? false : bDomainThyraMode_);

            RCP<MultiVector> tmpYblock = rangemaps_->getVector(row, Y.getNumVectors(), bRangeThyraMode_, false);  // subpart of Y vector containing part of solution of Xblock applied to Ablock
            Ablock->apply(*Xblock, *tmpYblock);

            Yblock->update(one, *tmpYblock, one);
          }
          rangemaps_->InsertVector(Yblock, row, tmpY, bRangeThyraMode_);
        }

      } else if (mode == Teuchos::TRANS) {
        // TODO: test me!
        for (size_t col = 0; col < Cols(); col++) {
          RCP<MultiVector>    Yblock = domainmaps_->getVector(col, Y.getNumVectors(), bDomainThyraMode_, true);

          for (size_t row = 0; row<Rows(); row++) {
            RCP<Matrix>            Ablock = getMatrix(row, col);

            if (Ablock.is_null())
              continue;

            // check whether Ablock is itself a blocked operator
            // If it is a blocked operator we have to provide Xpetra style GIDs, i.e. we have to transform GIDs
            bool bBlockedSubMatrix = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(Ablock) == Teuchos::null ? false : true;

            RCP<const MultiVector> Xblock = Teuchos::null;

            // extract sub part of X using Xpetra or Thyra GIDs
            if(bBlockedX) Xblock = rangemaps_->ExtractVector(refbX, row, bRangeThyraMode_);
            else          Xblock = rangemaps_->ExtractVector(refX,  row, bBlockedSubMatrix == true ? false : bRangeThyraMode_);
            RCP<MultiVector> tmpYblock = domainmaps_->getVector(col, Y.getNumVectors(), bDomainThyraMode_, false);
            Ablock->apply(*Xblock, *tmpYblock, Teuchos::TRANS);

            Yblock->update(one, *tmpYblock, one);
          }
          domainmaps_->InsertVector(Yblock, col, tmpY, bDomainThyraMode_);
        }
      }
      Y.update(alpha, *tmpY, beta);
    }

    //! \brief Returns the Map associated with the full domain of this operator.
    RCP<const Map > getFullDomainMap() const        { XPETRA_MONITOR("XpetraBlockedCrsMatrix::getFullDomainMap()"); return domainmaps_->getFullMap(); }

    //! \brief Returns the BlockedMap associated with the domain of this operator.
    RCP<const BlockedMap > getBlockedDomainMap() const     { XPETRA_MONITOR("XpetraBlockedCrsMatrix::getBlockedDomainMap()"); return domainmaps_->getBlockedMap(); }

    //! \brief Returns the Map associated with the domain of this operator.
    RCP<const Map > getDomainMap() const            { XPETRA_MONITOR("XpetraBlockedCrsMatrix::getDomainMap()"); return domainmaps_->getMap(); /*domainmaps_->getFullMap();*/ }

    //! \brief Returns the Map associated with the i'th block domain of this operator.
    RCP<const Map > getDomainMap(size_t i) const    { XPETRA_MONITOR("XpetraBlockedCrsMatrix::getDomainMap(size_t)"); return domainmaps_->getMap(i, bDomainThyraMode_); }

    //! \brief Returns the Map associated with the i'th block domain of this operator.
    RCP<const Map > getDomainMap(size_t i, bool bThyraMode) const    { XPETRA_MONITOR("XpetraBlockedCrsMatrix::getDomainMap(size_t,bool)"); return domainmaps_->getMap(i, bThyraMode); }

    //! Returns the Map associated with the full range of this operator.
    RCP<const Map > getFullRangeMap() const         { XPETRA_MONITOR("XpetraBlockedCrsMatrix::getRangeMap()"); return rangemaps_->getFullMap(); }

    //! \brief Returns the BlockedMap associated with the range of this operator.
    RCP<const BlockedMap > getBlockedRangeMap() const     { XPETRA_MONITOR("XpetraBlockedCrsMatrix::getBlockedRangeMap()"); return rangemaps_->getBlockedMap(); }

    //! Returns the Map associated with the range of this operator.
    RCP<const Map > getRangeMap() const             { XPETRA_MONITOR("XpetraBlockedCrsMatrix::getRangeMap()"); return rangemaps_->getMap(); /*rangemaps_->getFullMap();*/ }

    //! Returns the Map associated with the i'th block range of this operator.
    RCP<const Map > getRangeMap(size_t i) const     { XPETRA_MONITOR("XpetraBlockedCrsMatrix::getRangeMap(size_t)"); return rangemaps_->getMap(i, bRangeThyraMode_); }

    //! Returns the Map associated with the i'th block range of this operator.
    RCP<const Map > getRangeMap(size_t i, bool bThyraMode) const     { XPETRA_MONITOR("XpetraBlockedCrsMatrix::getRangeMap(size_t,bool)"); return rangemaps_->getMap(i, bThyraMode); }

    //! Returns map extractor class for range map
    RCP<const MapExtractor> getRangeMapExtractor() const { XPETRA_MONITOR("XpetraBlockedCrsMatrix::getRangeMapExtractor()"); return rangemaps_; }

    //! Returns map extractor for domain map
    RCP<const MapExtractor> getDomainMapExtractor() const { XPETRA_MONITOR("XpetraBlockedCrsMatrix::getDomainMapExtractor()"); return domainmaps_; }

    //@}

    //! Special multiplication routine (for BGS/Jacobi smoother)
    //{@

    /*! \brief Computes the sparse matrix-multivector multiplication (plus linear combination with input/result vector)
     *
     *  Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exception:
     *  - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
     *  - calculates result only for blocked row "row"
     *  - useful for BGS/Jacobi smoother in MueLu: there we have to calculate the residual for the current block row
     *    we can skip the MatVec calls in all other block rows
     */
    virtual void bgs_apply(
        const MultiVector& X, ///< Vector to be multiplied by matrix (input)
        MultiVector& Y, ///< result vector
        size_t row, ///< Index of block row to be treated
        Teuchos::ETransp mode = Teuchos::NO_TRANS, ///< Transpose mode
        Scalar alpha = ScalarTraits<Scalar>::one(), ///< scaling factor for result of matrix-vector product
        Scalar beta  = ScalarTraits<Scalar>::zero() ///< scaling factor for linear combination with result vector
        ) const
    {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::bgs_apply");
      //using Teuchos::RCP;

      TEUCHOS_TEST_FOR_EXCEPTION(mode != Teuchos::NO_TRANS && mode != Teuchos::TRANS, Xpetra::Exceptions::RuntimeError,
                                 "apply() only supports the following modes: NO_TRANS and TRANS." );

      // check whether input parameters are blocked or not
      RCP<const MultiVector>         refX = rcpFromRef(X);
      RCP<const BlockedMultiVector> refbX = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(refX);
      //RCP<MultiVector>               tmpY = rcpFromRef(Y);
      //RCP<BlockedMultiVector>       tmpbY = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(tmpY);

      bool bBlockedX = (refbX != Teuchos::null) ? true : false;

      // create (temporary) vectors for output
      // In the end we call Y.update(alpha, *tmpY, beta). Therefore we need a new vector storing the temporary results
      RCP<MultiVector> tmpY = MultiVectorFactory::Build(Y.getMap(), Y.getNumVectors(), true);

      SC one = ScalarTraits<SC>::one();

      if (mode == Teuchos::NO_TRANS) {
        RCP<MultiVector>    Yblock = rangemaps_->getVector(row, Y.getNumVectors(), bRangeThyraMode_, true);
        for (size_t col = 0; col < Cols(); col++) {

          // extract matrix block
          RCP<Matrix> Ablock = getMatrix(row, col);

          if (Ablock.is_null())
            continue;

          // check whether Ablock is itself a blocked operator
          // If it is a blocked operator we have to provide Xpetra style GIDs, i.e. we have to transform GIDs
          bool bBlockedSubMatrix = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(Ablock) == Teuchos::null ? false : true;

          // input/output vectors for local block operation
          RCP<const MultiVector> Xblock = Teuchos::null; // subpart of X vector to be applied to subblock of A

          // extract sub part of X using Xpetra or Thyra GIDs
          // if submatrix is again blocked, we extract it using Xpetra style gids. If it is a single
          // block matrix we use the Thyra or Xpetra style GIDs that are used to store the matrix
          if(bBlockedX) Xblock = domainmaps_->ExtractVector(refbX, col, bDomainThyraMode_);
          else          Xblock = domainmaps_->ExtractVector(refX,  col, bBlockedSubMatrix == true ? false : bDomainThyraMode_);

          RCP<MultiVector> tmpYblock = rangemaps_->getVector(row, Y.getNumVectors(), bRangeThyraMode_, false);  // subpart of Y vector containing part of solution of Xblock applied to Ablock
          Ablock->apply(*Xblock, *tmpYblock);

          Yblock->update(one, *tmpYblock, one);
        }
        rangemaps_->InsertVector(Yblock, row, tmpY, bRangeThyraMode_);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true,Xpetra::Exceptions::NotImplemented,"Xpetar::BlockedCrsMatrix::bgs_apply: not implemented for transpose case.");
      }
      Y.update(alpha, *tmpY, beta);
    }


    //@}


    //! Implements DistObject interface
    //{@

    //! Access function for the Tpetra::Map this DistObject was constructed with.
    const Teuchos::RCP< const Map > getMap() const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getMap");
      if (Rows() == 1 && Cols () == 1) {
        return getMatrix(0,0)->getMap();
      }
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::getMap(): operation not supported.");
    }

    //! Import.
    void doImport(const Matrix &source, const Import& importer, CombineMode CM) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::doImport");
      if (Rows() == 1 && Cols () == 1) {
        getMatrix(0,0)->doImport(source, importer, CM);
        return;
      }
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doImport(): operation not supported.");
    }

    //! Export.
    void doExport(const Matrix& dest, const Import& importer, CombineMode CM) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::doExport");
      if (Rows() == 1 && Cols () == 1) {
        getMatrix(0,0)->doExport(dest, importer, CM);
        return;
      }
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doExport(): operation not supported.");
    }

    //! Import (using an Exporter).
    void doImport(const Matrix& source, const Export& exporter, CombineMode CM) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::doImport");
      if (Rows() == 1 && Cols () == 1) {
        getMatrix(0,0)->doImport(source, exporter, CM);
        return;
      }
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doImport(): operation not supported.");
    }

    //! Export (using an Importer).
    void doExport(const Matrix& dest, const Export& exporter, CombineMode CM) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::doExport");
      if (Rows() == 1 && Cols () == 1) {
        getMatrix(0,0)->doExport(dest, exporter, CM);
        return;
      }
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::doExport(): operation not supported.");
    }

    // @}

    //! @name Overridden from Teuchos::Describable
    //@{

    /** \brief Return a simple one-line description of this object. */
    std::string description() const                       { return "Xpetra_BlockedCrsMatrix.description()"; }

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const {
      out << "Xpetra::BlockedCrsMatrix: " << Rows() << " x " << Cols() << std::endl;

      if (isFillComplete()) {
        out << "BlockMatrix is fillComplete" << std::endl;

        /*if(fullrowmap_ != Teuchos::null) {
          out << "fullRowMap" << std::endl;
          fullrowmap_->describe(out,verbLevel);
        } else {
          out << "fullRowMap not set. Check whether block matrix is properly fillCompleted!" << std::endl;
        }*/

        //out << "fullColMap" << std::endl;
        //fullcolmap_->describe(out,verbLevel);

      } else {
        out << "BlockMatrix is NOT fillComplete" << std::endl;
      }

      for (size_t r = 0; r < Rows(); ++r)
        for (size_t c = 0; c < Cols(); ++c) {
          if(getMatrix(r,c)!=Teuchos::null) {
            out << "Block(" << r << "," << c << ")" << std::endl;
            getMatrix(r,c)->describe(out,verbLevel);
          } else out << "Block(" << r << "," << c << ") = null" << std::endl;
        }
    }

    //! @name Overridden from Teuchos::LabeledObject
    //@{
    void setObjectLabel( const std::string &objectLabel ) {
      XPETRA_MONITOR("TpetraBlockedCrsMatrix::setObjectLabel");
      for (size_t r = 0; r < Rows(); ++r)
        for (size_t c = 0; c < Cols(); ++c) {
          if(getMatrix(r,c)!=Teuchos::null) {
            std::ostringstream oss; oss<< objectLabel << "(" << r << "," << c << ")";
            getMatrix(r,c)->setObjectLabel(oss.str());
          }
        }
    }
    //@}


    //! Supports the getCrsGraph() call
    bool hasCrsGraph() const {
      if (Rows() == 1 && Cols () == 1) return true;
      else return false;
    }

    //! Returns the CrsGraph associated with this matrix.
    RCP<const CrsGraph> getCrsGraph() const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getCrsGraph");
      if (Rows() == 1 && Cols () == 1) {
        return getMatrix(0,0)->getCrsGraph();
      }
      throw Xpetra::Exceptions::RuntimeError("getCrsGraph() not supported by BlockedCrsMatrix");
    }

    //@}

    //! @name Block matrix access
    //@{

    virtual bool isDiagonal() const {return is_diagonal_;}

    /// number of row blocks
    virtual size_t Rows() const                                       { XPETRA_MONITOR("XpetraBlockedCrsMatrix::Rows"); return rangemaps_->NumMaps(); }

    /// number of column blocks
    virtual size_t Cols() const                                       { XPETRA_MONITOR("XpetraBlockedCrsMatrix::Cols"); return domainmaps_->NumMaps(); }

    /// return unwrap 1x1 blocked operators
    Teuchos::RCP<Matrix> getCrsMatrix() const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getCrsMatrix");
      TEUCHOS_TEST_FOR_EXCEPTION(Rows()!=1, std::out_of_range, "Can only unwrap a 1x1 blocked matrix. The matrix has " << Rows() << " block rows, though.");
      TEUCHOS_TEST_FOR_EXCEPTION(Cols()!=1, std::out_of_range, "Can only unwrap a 1x1 blocked matrix. The matrix has " << Cols() << " block columns, though.");

      RCP<Matrix> mat = getMatrix(0,0);
      RCP<BlockedCrsMatrix> bmat = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(mat);
      if (bmat == Teuchos::null) return mat;
      return bmat->getCrsMatrix();
    }

    /// helper routine recursively returns the first inner-most non-null matrix block from a (nested) blocked operator
    Teuchos::RCP<Matrix> getInnermostCrsMatrix() {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getInnermostCrsMatrix");
      size_t row = Rows()+1, col = Cols()+1;
      for (size_t r = 0; r < Rows(); ++r)
        for(size_t c = 0; c < Cols(); ++c)
          if (getMatrix(r,c) != Teuchos::null) {
            row = r;
            col = c;
            break;
          }
      TEUCHOS_TEST_FOR_EXCEPTION(row == Rows()+1 || col == Cols()+1, Xpetra::Exceptions::Incompatible, "Xpetra::BlockedCrsMatrix::getInnermostCrsMatrix: Could not find a non-zero sub-block in blocked operator.")
      RCP<Matrix> mm = getMatrix(row,col);
      RCP<BlockedCrsMatrix> bmat = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(mm);
      if (bmat == Teuchos::null) return mm;
      return bmat->getInnermostCrsMatrix();
    }

    /// return block (r,c)
    Teuchos::RCP<Matrix> getMatrix(size_t r, size_t c) const       {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::getMatrix");
      TEUCHOS_TEST_FOR_EXCEPTION(r > Rows(), std::out_of_range, "Error, r = " << Rows() << " is too big");
      TEUCHOS_TEST_FOR_EXCEPTION(c > Cols(), std::out_of_range, "Error, c = " << Cols() << " is too big");

      // transfer strided/blocked map information
      /*      if (blocks_[r*Cols()+c] != Teuchos::null &&
          blocks_[r*Cols()+c]->IsView("stridedMaps") == false)
          blocks_[r*Cols()+c]->CreateView("stridedMaps", getRangeMap(r,bRangeThyraMode_), getDomainMap(c,bDomainThyraMode_));*/
      return blocks_[r*Cols()+c];
    }

    /// set matrix block
    //void setMatrix(size_t r, size_t c, Teuchos::RCP<CrsMatrix> mat) {
    void setMatrix(size_t r, size_t c, Teuchos::RCP<Matrix> mat) {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::setMatrix");
      // TODO: if filled -> return error

      TEUCHOS_TEST_FOR_EXCEPTION(r > Rows(), std::out_of_range, "Error, r = " << Rows() << " is too big");
      TEUCHOS_TEST_FOR_EXCEPTION(c > Cols(), std::out_of_range, "Error, c = " << Cols() << " is too big");
      if(!mat.is_null() && r!=c) is_diagonal_ = false;
      // set matrix
      blocks_[r*Cols() + c] = mat;
    }

    /// merge BlockedCrsMatrix blocks in a CrsMatrix
    // NOTE: This is a rather expensive operation, since all blocks are copied
    // into a new big CrsMatrix
    Teuchos::RCP<Matrix> Merge() const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::Merge");
      using Teuchos::RCP;
      using Teuchos::rcp_dynamic_cast;
      Scalar one = ScalarTraits<SC>::one();

      TEUCHOS_TEST_FOR_EXCEPTION(bRangeThyraMode_ != bDomainThyraMode_, Xpetra::Exceptions::RuntimeError,
                                 "BlockedCrsMatrix::Merge: only implemented for Xpetra-style or Thyra-style numbering. No mixup allowed!" );

      TEUCHOS_TEST_FOR_EXCEPTION(isFillComplete() == false, Xpetra::Exceptions::RuntimeError,
                                 "BlockedCrsMatrix::Merge: BlockMatrix must be fill-completed." );

      LocalOrdinal lclNumRows = getFullRangeMap()->getLocalNumElements();
      Teuchos::ArrayRCP<size_t> numEntPerRow (lclNumRows);
      for (LocalOrdinal lclRow = 0; lclRow < lclNumRows; ++lclRow)
        numEntPerRow[lclRow] = getNumEntriesInLocalRow(lclRow);

      RCP<Matrix> sparse = MatrixFactory::Build(getFullRangeMap(), numEntPerRow);

      if(bRangeThyraMode_ == false) {
        // Xpetra mode
        for (size_t i = 0; i < Rows(); i++) {
          for (size_t j = 0; j < Cols(); j++) {
            if (getMatrix(i,j) != Teuchos::null) {
              RCP<const Matrix> mat = getMatrix(i,j);

              // recursively call Merge routine
              RCP<const BlockedCrsMatrix> bMat = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(mat);
              if (bMat != Teuchos::null) mat = bMat->Merge();

              bMat = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(mat);
              TEUCHOS_TEST_FOR_EXCEPTION(bMat != Teuchos::null, Xpetra::Exceptions::RuntimeError,
                                         "BlockedCrsMatrix::Merge: Merging of blocked sub-operators failed?!" );

              // jump over empty blocks
              if(mat->getLocalNumEntries() == 0) continue;

              this->Add(*mat, one, *sparse, one);
            }
          }
        }
      } else {
        // Thyra mode
        for (size_t i = 0; i < Rows(); i++) {
          for (size_t j = 0; j < Cols(); j++) {
            if (getMatrix(i,j) != Teuchos::null) {
              RCP<const Matrix> mat = getMatrix(i,j);
              // recursively call Merge routine
              RCP<const BlockedCrsMatrix> bMat = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(mat);
              if (bMat != Teuchos::null) mat = bMat->Merge();

              bMat = Teuchos::rcp_dynamic_cast<const BlockedCrsMatrix>(mat);
              TEUCHOS_TEST_FOR_EXCEPTION(bMat != Teuchos::null, Xpetra::Exceptions::RuntimeError,
                                         "BlockedCrsMatrix::Merge: Merging of blocked sub-operators failed?!" );

              // check whether we have a CrsMatrix block (no blocked operator)
              RCP<const CrsMatrixWrap> crsMat = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(mat);
              TEUCHOS_ASSERT(crsMat != Teuchos::null);

              // these are the thyra style maps of the matrix
              RCP<const Map> trowMap = mat->getRowMap();
              RCP<const Map> tcolMap = mat->getColMap();
              RCP<const Map> tdomMap = mat->getDomainMap();

              // get Xpetra maps
              RCP<const Map> xrowMap = getRangeMapExtractor()->getMap(i,false);
              RCP<const Map> xdomMap = getDomainMapExtractor()->getMap(j,false);

              // generate column map with Xpetra GIDs
              // We have to do this separately for each block since the column
              // map of each block might be different in the same block column
              Teuchos::RCP<Map> xcolMap = MapUtils::transformThyra2XpetraGIDs(
                    *tcolMap,
                    *tdomMap,
                    *xdomMap);

              // jump over empty blocks
              if(mat->getLocalNumEntries() == 0) continue;

              size_t maxNumEntries = mat->getLocalMaxNumRowEntries();

              size_t    numEntries;
              Array<GO> inds (maxNumEntries);
              Array<GO> inds2(maxNumEntries);
              Array<SC> vals (maxNumEntries);

              // loop over all rows and add entries
              for (size_t k = 0; k < mat->getLocalNumRows(); k++) {
                GlobalOrdinal rowTGID = trowMap->getGlobalElement(k);
                crsMat->getCrsMatrix()->getGlobalRowCopy(rowTGID, inds(), vals(), numEntries);

                // create new indices array
                for (size_t l = 0; l < numEntries; ++l) {
                  LocalOrdinal lid = tcolMap->getLocalElement(inds[l]);
                  inds2[l] = xcolMap->getGlobalElement(lid);
                }

                GlobalOrdinal rowXGID = xrowMap->getGlobalElement(k);
                sparse->insertGlobalValues(
                    rowXGID, inds2(0, numEntries),
                    vals(0, numEntries));
              }
            }
          }
        }
      }

      sparse->fillComplete(getFullDomainMap(), getFullRangeMap());

      TEUCHOS_TEST_FOR_EXCEPTION(sparse->getLocalNumEntries() != getLocalNumEntries(), Xpetra::Exceptions::RuntimeError,
                                 "BlockedCrsMatrix::Merge: Local number of entries of merged matrix does not coincide with local number of entries of blocked operator." );

      TEUCHOS_TEST_FOR_EXCEPTION(sparse->getGlobalNumEntries() != getGlobalNumEntries(), Xpetra::Exceptions::RuntimeError,
                                 "BlockedCrsMatrix::Merge: Global number of entries of merged matrix does not coincide with global number of entries of blocked operator." );

      return sparse;
    }
    //@}

    typedef typename CrsMatrix::local_matrix_type local_matrix_type;
    /// \brief Access the underlying local Kokkos::CrsMatrix object
    local_matrix_type getLocalMatrixDevice () const {
      if (Rows() == 1 && Cols () == 1) {
        return getMatrix(0,0)->getLocalMatrixDevice();
      }
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::getLocalMatrix(): operation not supported.");
    }
    /// \brief Access the underlying local Kokkos::CrsMatrix object
    typename local_matrix_type::HostMirror getLocalMatrixHost () const {
      if (Rows() == 1 && Cols () == 1) {
        return getMatrix(0,0)->getLocalMatrixHost();
      }
      throw Xpetra::Exceptions::RuntimeError("BlockedCrsMatrix::getLocalMatrix(): operation not supported.");
    }


#ifdef HAVE_XPETRA_THYRA
    Teuchos::RCP<Thyra::BlockedLinearOpBase<Scalar> > getThyraOperator() {
      Teuchos::RCP<Thyra::LinearOpBase<Scalar> > thOp =
          Xpetra::ThyraUtils<Scalar,LO,GO,Node>::toThyra(Teuchos::rcpFromRef(*this));
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thOp));

      Teuchos::RCP<Thyra::BlockedLinearOpBase<Scalar> > thbOp =
          Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<Scalar> >(thOp);
      TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(thbOp));
      return thbOp;
    }
#endif
    //! Returns the block size of the storage mechanism
    LocalOrdinal GetStorageBlockSize() const {return 1;}


    //! Compute a residual R = B - (*this) * X
    void residual(const MultiVector & X,
                  const MultiVector & B,
                  MultiVector & R) const {
      using STS = Teuchos::ScalarTraits<Scalar>;
      R.update(STS::one(),B,STS::zero());
      this->apply (X, R, Teuchos::NO_TRANS, -STS::one(), STS::one());
    }

  private:

    /** \name helper functions */
    //@{

    /// Add a Xpetra::CrsMatrix to another: B = B*scalarB + A*scalarA
    /**
     * Note, that this routine works only correctly if A only has entries which are empty (zero) in B.
     * We use the insertGlobalValues routine for inserting the new values from A in B. The sumIntoGlobalValues
     * routine is not implemented in Xpetra (and would not extend the graph of B for new entries).
     * Here we need something to catch the exceptions of a future implementation of sumIntoGlobalValues that
     * then adds the remaining new entries with insertGlobal Values.
     *
     * This routine is private and used only by Merge. Since the blocks in BlockedCrsMatrix are seperated,
     * this routine works for merging a BlockedCrsMatrix.
     */
    void Add(const Matrix& A, const Scalar scalarA, Matrix& B, const Scalar scalarB) const {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::Add");
      TEUCHOS_TEST_FOR_EXCEPTION(!A.isFillComplete(), Xpetra::Exceptions::RuntimeError,
                                 "Matrix A is not completed");
      using Teuchos::Array;
      using Teuchos::ArrayView;

      B.scale(scalarB);

      Scalar one  = ScalarTraits<SC>::one();
      Scalar zero = ScalarTraits<SC>::zero();

      if (scalarA == zero)
        return;

      Teuchos::RCP<const Matrix> rcpA = Teuchos::rcpFromRef(A);
      Teuchos::RCP<const CrsMatrixWrap> rcpAwrap = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(rcpA);
      TEUCHOS_TEST_FOR_EXCEPTION(rcpAwrap == Teuchos::null, Xpetra::Exceptions::BadCast,
                                       "BlockedCrsMatrix::Add: matrix A must be of type CrsMatrixWrap.");
      Teuchos::RCP<const CrsMatrix> crsA = rcpAwrap->getCrsMatrix();

      size_t maxNumEntries = crsA->getLocalMaxNumRowEntries();

      size_t    numEntries;
      Array<GO> inds(maxNumEntries);
      Array<SC> vals(maxNumEntries);

      RCP<const Map> rowMap = crsA->getRowMap();
      RCP<const Map> colMap = crsA->getColMap();

      ArrayView<const GO> rowGIDs = crsA->getRowMap()->getLocalElementList();
      for (size_t i = 0; i < crsA->getLocalNumRows(); i++) {
        GO row = rowGIDs[i];
        crsA->getGlobalRowCopy(row, inds(), vals(), numEntries);

        if (scalarA != one)
          for (size_t j = 0; j < numEntries; ++j)
            vals[j] *= scalarA;

        B.insertGlobalValues(row, inds(0, numEntries), vals(0, numEntries)); // insert should be ok, since blocks in BlockedCrsOpeartor do not overlap!
      }
    }

    //@}

    // Default view is created after fillComplete()
    // Because ColMap might not be available before fillComplete().
    void CreateDefaultView() {
      XPETRA_MONITOR("XpetraBlockedCrsMatrix::CreateDefaultView");

      // Create default view
      this->defaultViewLabel_ = "point";
      this->CreateView(this->GetDefaultViewLabel(), getRangeMap(), getDomainMap());

      // Set current view
      this->currentViewLabel_ = this->GetDefaultViewLabel();
    }

  private:
    bool is_diagonal_; ///< If we're diagonal, a bunch of the extraction stuff should work
    Teuchos::RCP<const MapExtractor> domainmaps_; ///< full domain map together with all partial domain maps
    Teuchos::RCP<const MapExtractor> rangemaps_; ///< full range map together with all partial domain maps

    std::vector<Teuchos::RCP<Matrix> > blocks_; ///< row major matrix block storage
#ifdef HAVE_XPETRA_THYRA
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<Scalar> > thyraOp_; ///< underlying thyra operator
#endif
    bool bRangeThyraMode_; ///< boolean flag, which is true, if BlockedCrsMatrix has been created using Thyra-style numbering for sub blocks, i.e. all GIDs of submaps are contiguous and start from 0.
    bool bDomainThyraMode_; ///< boolean flag, which is true, if BlockedCrsMatrix has been created using Thyra-style numbering for sub blocks, i.e. all GIDs of submaps are contiguous and start from 0.

};

} //namespace Xpetra

#define XPETRA_BLOCKEDCRSMATRIX_SHORT
#endif /* XPETRA_BLOCKEDCRSMATRIX_HPP */
