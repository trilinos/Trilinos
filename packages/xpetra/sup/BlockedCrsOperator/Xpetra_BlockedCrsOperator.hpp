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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * Xpetra_BlockedCrsOperator.hpp
 *
 *  Created on: Aug 17, 2011
 *      Author: wiesner
 */

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_BLOCKEDCRSOPERATOR_HPP_
#define XPETRA_BLOCKEDCRSOPERATOR_HPP_

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Hashtable.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_MapFactory.hpp"
#include "Xpetra_MultiVector.hpp"
#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_CrsMatrixFactory.hpp"

#include "Xpetra_MapExtractor.hpp"

#include "Xpetra_Operator.hpp"

#define sumAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));

/** \file Xpetra_Operator.hpp

  Declarations for the class Xpetra::Operator.
*/
namespace Xpetra {

  typedef std::string viewLabel_t;

template <class Scalar,
          class LocalOrdinal  = int,
          class GlobalOrdinal = LocalOrdinal,
          class Node          = Kokkos::DefaultNode::DefaultNodeType,
          class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps > //TODO: or BlockSparseOp ?
class BlockedCrsOperator : public Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {

  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> MapClass;
  typedef Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> MapExtractorClass;
  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixClass;
  typedef Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Operator;
  typedef Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsGraph;
#ifdef HAVE_CTHULHU_TPETRA
  typedef Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsMatrix;
#endif
  typedef Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixFactory;
  typedef Xpetra::OperatorView<LocalOrdinal, GlobalOrdinal, Node> OperatorView;

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
  BlockedCrsOperator(Teuchos::RCP<const MapExtractorClass>& rangeMaps,
                     Teuchos::RCP<const MapExtractorClass>& domainMaps,
                     size_t npr,
                     Xpetra::ProfileType pftype = Xpetra::DynamicProfile)
  : domainmaps_(domainMaps),
    rangemaps_(rangeMaps)
  {
    // Set matrix data
    //matrixData_ = Teuchos::null; //CrsMatrixFactory::Build(rowMap, maxNumEntriesPerRow, pftype); // TODO remove me

    blocks_.reserve(Rows()*Cols());

    // add CrsMatrix objects in row,column order
    for (size_t r=0; r<Rows(); ++r)
    {
      for(size_t c=0; c<Cols(); ++c)
      {
        Teuchos::RCP<CrsMatrixClass> matblock = CrsMatrixFactory::Build(getRangeMap(r), npr, pftype);
        blocks_.push_back(matblock);
      }
    }

    // Default view
    CreateDefaultView();
  }

  //! Destructor
  virtual ~BlockedCrsOperator() {}

  //@}


  //! @name Insertion/Removal Methods
  //@{

  //! Insert matrix entries, using global IDs.
  /** All index values must be in the global space.
      \pre \c globalRow exists as an ID in the global row map
      \pre <tt>isLocallyIndexed() == false</tt>
      \pre <tt>isStorageOptimized() == false</tt>

      \post <tt>isGloballyIndexed() == true</tt>

      \note If \c globalRow does not belong to the matrix on this node, then it will be communicated to the appropriate node when globalAssemble() is called (which will, at the latest, occur during the next call to fillComplete().) Otherwise, the entries will be inserted in the local matrix.
      \note If the matrix row already contains values at the indices corresponding to values in \c cols, then the new values will be summed with the old values; this may happen at insertion or during the next call to fillComplete().
      \note If <tt>hasColMap() == true</tt>, only (cols[i],vals[i]) where cols[i] belongs to the column map on this node will be inserted into the matrix.
  */
  void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals)
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
        "insertGlobalValues not supported by BlockedCrsOperator!" );
  }

  //! Insert matrix entries, using local IDs.
  /** All index values must be in the local space.
      \pre \c localRow exists as an ID in the local row map
      \pre <tt>isGloballyIndexed() == false</tt>
      \pre <tt>isStorageOptimized() == false</tt>

      \post <tt>isLocallyIndexed() == true</tt>
  */
  void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals) {
    TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
        "insertLocalValues not supported by BlockedCrsOperator!" );
  }

  //@}

  //! @name Transformational Methods
  //@{

    /*! \brief Signal that data entry is complete, specifying domain and range maps.

    Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

    \pre  <tt>isFillActive() == true<tt>
    \pre <tt>isFillComplete()() == false<tt>

    \post <tt>isFillActive() == false<tt>
    \post <tt>isFillComplete() == true<tt>
    \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
    */
  void fillComplete(const RCP<const MapClass> &domainMap, const RCP<const MapClass> &rangeMap, const RCP<ParameterList> &params=null)
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
        "fillComplete with arguments not supported for block matrices!" );
  }

  /*! \brief Signal that data entry is complete.

  Off-node entries are distributed (via globalAssemble()), repeated entries are summed, and global indices are transformed to local indices.

  \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

  \pre  <tt>isFillActive() == true<tt>
  \pre <tt>isFillComplete()() == false<tt>

  \post <tt>isFillActive() == false<tt>
  \post <tt>isFillComplete() == true<tt>
  \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
  */
  //TODO : Get ride of "Tpetra"::OptimizeOption
  void fillComplete(const RCP<ParameterList> &params=null)
  {
    for (size_t r=0; r<Rows(); ++r)
    {
      for (size_t c=0; c<Cols(); ++c)
      {
        if(getMatrix(r,c)->isFillComplete() == false)
          getMatrix(r,c)->fillComplete(getDomainMap(c),getRangeMap(r),params);
      }
    }

    // get full row map
    //fullrowmap_ = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(rangemaps_->getFullMap()->lib(), rangemaps_->getFullMap()->getNodeNumElements(), rangemaps_->getFullMap()->getNodeElementList(), rangemaps_->getFullMap()->getIndexBase(), rangemaps_->getFullMap()->getComm());//rangemaps_->FullMap(); //->Clone();
    fullrowmap_ = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(rangemaps_->getFullMap()->lib(), rangemaps_->getFullMap()->getGlobalNumElements(), rangemaps_->getFullMap()->getNodeElementList(), rangemaps_->getFullMap()->getIndexBase(), rangemaps_->getFullMap()->getComm());//rangemaps_->FullMap(); //->Clone();

    // TODO: check me, clean up, use only ArrayView instead of std::vector
    // build full col map
    if (fullcolmap_ == Teuchos::null)
    {
      std::vector<GlobalOrdinal> colmapentries;
      for (size_t c=0; c<Cols(); ++c)
      {
        std::set<GlobalOrdinal> colset;
        for (size_t r=0; r<Rows(); ++r)
        {
          if(getMatrix(r,c) != Teuchos::null) {
            Teuchos::RCP<const MapClass> colmap = getMatrix(r,c)->getColMap();
            copy(colmap->getNodeElementList().getRawPtr(),
          colmap->getNodeElementList().getRawPtr()+colmap->getNodeNumElements(),
          inserter(colset,colset.begin()));
          }
        }
        colmapentries.reserve(colmapentries.size()+colset.size());
        copy(colset.begin(), colset.end(), back_inserter(colmapentries));
        sort(colmapentries.begin(), colmapentries.end());
        typename std::vector<GlobalOrdinal>::iterator gendLocation;
        gendLocation = std::unique(colmapentries.begin(), colmapentries.end());
        colmapentries.erase(gendLocation,colmapentries.end());
      }

      // sum up number of local elements
      size_t numGlobalElements = 0;
      sumAll(rangemaps_->getFullMap()->getComm(), colmapentries.size(), numGlobalElements) 

      const Teuchos::ArrayView<const GlobalOrdinal> aView = Teuchos::ArrayView<const GlobalOrdinal>(colmapentries);
      fullcolmap_ = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(rangemaps_->getFullMap()->lib(), numGlobalElements, aView, 0,rangemaps_->getFullMap()->getComm());
      //fullcolmap_ = Xpetra::MapFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(rangemaps_->getFullMap()->lib(), domainmaps_->getFullMap()->getGlobalNumElements(), aView, 0,domainmaps_->getFullMap()->getComm());
    }


  }

  //@}

  //! Returns the number of global rows.
  /** Undefined if isFillActive().
   */
  global_size_t getGlobalNumRows() const
  {
    global_size_t globalNumRows= 0;
    for (size_t r=0; r<Rows(); ++r)
    {
      globalNumRows += getMatrix(r,0)->getGlobalNumRows();
    }
    return globalNumRows;
  }

  //! \brief Returns the number of global columns in the matrix.
  /** Undefined if isFillActive().
   */
  global_size_t getGlobalNumCols() const
  {
    global_size_t globalNumCols= 0;
    for (size_t c=0; c<Cols(); ++c)
    {
      globalNumCols += getMatrix(0,c)->getGlobalNumCols();
    }
    return globalNumCols;
  }

  //! Returns the number of matrix rows owned on the calling node.
  size_t getNodeNumRows() const
  {
    global_size_t nodeNumRows= 0;
    for (size_t r=0; r<Rows(); ++r)
    {
      nodeNumRows += getMatrix(r,0)->getNodeNumRows();
    }
    return nodeNumRows;
  }

  //! Returns the global number of entries in this matrix.
  global_size_t getGlobalNumEntries() const
  {
    global_size_t globalNumEntries= 0;
    for (size_t r=0; r<Rows(); ++r)
    {
      for(size_t c=0; c<Cols(); ++c)
      {
        globalNumEntries += getMatrix(r,c)->getGlobalNumEntries();
      }
    }
    return globalNumEntries;
  }

  //! Returns the local number of entries in this matrix.
  size_t getNodeNumEntries() const
  {
    global_size_t nodeNumEntries= 0;
    for (size_t r=0; r<Rows(); ++r)
    {
      for(size_t c=0; c<Cols(); ++c)
      {
        nodeNumEntries += getMatrix(r,c)->getNodeNumEntries();
      }
    }
    return nodeNumEntries;
  }

  //! Returns the current number of entries on this node in the specified local row.
  /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
  size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
        "getNumEntriesInLocalRow not supported by BlockedCrsOperator!" );
    return 0;
  }

  //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons.
  /** Undefined if isFillActive().
   */
  global_size_t getGlobalNumDiags() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
        "getGlobalNumDiags() not supported by BlockedCrsOperator!" );
    return 0;
  }

  //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons.
  /** Undefined if isFillActive().
   */
  size_t getNodeNumDiags() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
        "getNodeNumDiags() not supported by BlockedCrsOperator!" );
    return 0;
  }

  //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
  /** Undefined if isFillActive().
   */
  size_t getGlobalMaxNumRowEntries() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
        "getGlobalMaxNumRowEntries() not supported by BlockedCrsOperator!" );
    return 0;
  }

  //! \brief Returns the maximum number of entries across all rows/columns on this node.
  /** Undefined if isFillActive().
   */
  size_t getNodeMaxNumRowEntries() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
        "getNodeMaxNumRowEntries() not supported by BlockedCrsOperator!" );
    return 0;
  }

  //! \brief If matrix indices of all matrix blocks are in the local range, this function returns true. Otherwise, this function returns false.
  /** if false, then this does not automatically mean that all blocks are globally indexed. The user has to make sure, that all matrix blocks
   * are indexed in the same way (locally or globally). Otherwise the block matrix is not valid...
  */
  bool isLocallyIndexed() const
  {
    for (size_t i=0; i<blocks_.size(); ++i)
      if (not blocks_[i]->isLocallyIndexed())
        return false;
    return true;
  }

  //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false.
  /** if false, then this does not automatically mean that all blocks are locally indexed. The user has to make sure, that all matrix blocks
   * are indexed in the same way (locally or globally). Otherwise the block matrix is not valid...
  */
  bool isGloballyIndexed() const
  {
    for (size_t i=0; i<blocks_.size(); ++i)
      if (not blocks_[i]->isGloballyIndexed())
        return false;
    return true;
  }

  //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
  bool isFillComplete() const
  {
    for (size_t i=0; i<blocks_.size(); ++i)
      if (not blocks_[i]->isFillComplete())
        return false;
    return true;
  }

  //! Extract a list of entries in a specified local row of the matrix. Put into storage allocated by calling routine.
  /*!
    \param LocalRow - (In) Local row number for which indices are desired.
    \param Indices - (Out) Local column indices corresponding to values.
    \param Values - (Out) Matrix values.
    \param NumIndices - (Out) Number of indices.

    Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
    with row \c LocalRow. If \c LocalRow is not valid for this node, then \c Indices and \c Values are unchanged and \c NumIndices is
    returned as OrdinalTraits<size_t>::invalid().

    \pre <tt>isLocallyIndexed()==true</tt> or <tt>hasColMap() == true</tt>
  */
    virtual void getLocalRowCopy(LocalOrdinal LocalRow,
                                 const ArrayView<LocalOrdinal> &Indices,
                                 const ArrayView<Scalar> &Values,
                                 size_t &NumEntries
                                 ) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
          "getLocalRowCopy not supported by BlockedCrsOperator!" );
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
  void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const Scalar> &values) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
        "getGlobalRowView not supported by BlockedCrsOperator!" );
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
  void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices, ArrayView<const Scalar> &values) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
        "getLocalRowView not supported by BlockedCrsOperator!" );
  }

  //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
  /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the
    the zero and non-zero diagonals owned by this node. */
  void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
        "getLocalDiagCopy not supported by BlockedCrsOperator!" );
  }

  //! Get Frobenius norm of the matrix
  virtual typename ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
        "getFrobeniusNorm() not supported by BlockedCrsOperator, yet!" );
  }

  //@}

  //! @name Advanced Matrix-vector multiplication and solve methods
  //@{

  //! Multiplies this matrix by a MultiVector.
  /*! \c X is required to be post-imported, i.e., described by the column map of the matrix. \c Y is required to be pre-exported, i.e., described by the row map of the matrix.

  Both are required to have constant stride, and they are not permitted to ocupy overlapping space. No runtime checking will be performed in a non-debug build.

  This method is templated on the scalar type of MultiVector objects, allowing this method to be applied to MultiVector objects of arbitrary type. However, it is recommended that multiply() not be called directly; instead, use the CrsMatrixMultiplyOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.

  If \c beta is equal to zero, the operation will enjoy overwrite semantics (\c Y will be overwritten with the result of the multiplication). Otherwise, the result of the multiplication
  will be accumulated into \c Y.
  */
  //TODO virtual=0 // TODO: Add default parameters ?
//   void multiply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, Scalar alpha, Scalar beta) const {
//      matrixData_->multiply(X, Y, trans, alpha, beta);
//   }

  //@}

  //! @name Methods implementing Operator
  //@{

  //! \brief Computes the sparse matrix-multivector multiplication.
  /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
    - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
  */
  virtual void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
                     Teuchos::ETransp mode = Teuchos::NO_TRANS,
                     Scalar alpha = ScalarTraits<Scalar>::one(),
                     Scalar beta = ScalarTraits<Scalar>::zero()) const
  {
    // TODO: check maps

    Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tX   = Teuchos::rcp(&X,false);
    Teuchos::RCP<      MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > tmpY = MultiVectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(Y.getMap(), Y.getNumVectors());

    if (mode == Teuchos::NO_TRANS)
    {
      for(size_t rblock=0; rblock<Rows(); ++rblock)
      {
        Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rowresult = rangemaps_->getVector(rblock,Y.getNumVectors()); // end result for block row
        Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rowy      = rangemaps_->getVector(rblock,Y.getNumVectors()); // helper vector
        for (size_t cblock=0; cblock<Cols(); ++cblock)
        {
          Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > colx = domainmaps_->ExtractVector(tX,cblock);
          Teuchos::RCP<CrsMatrixClass> bmat = getMatrix(rblock,cblock);
          bmat->apply(*colx,*rowy);
          rowresult->update(ScalarTraits< Scalar >::one(),*rowy,ScalarTraits< Scalar >::one());
        }
        rangemaps_->InsertVector(rowresult,rblock,tmpY);
      }
      Y.update(alpha,*tmpY,beta);
    }
    else if (mode == Teuchos::TRANS)
    {
      // TODO: test me!
      for (size_t cblock = 0; cblock<Cols(); ++cblock)
      {
        Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rowresult = domainmaps_->getVector(cblock,Y.getNumVectors()); // end result for block row
        Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > rowy      = domainmaps_->getVector(cblock,Y.getNumVectors()); // helper vector
        for (size_t rblock = 0; rblock<Rows(); ++rblock)
        {
          Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > colx = rangemaps_->ExtractVector(tX,rblock);
          Teuchos::RCP<CrsMatrixClass> bmat = getMatrix(rblock,cblock);
          bmat->apply(*colx,*rowy,Teuchos::TRANS);
          rowresult->update(ScalarTraits< Scalar >::one(),*rowy,ScalarTraits< Scalar >::one());
        }
        domainmaps_->InsertVector(rowresult,cblock,tmpY);
      }
      Y.update(alpha,*tmpY,beta);
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
          "apply() only supports NO_TRANS and TRANS." );

  }

  //! \brief Returns the Map associated with the full domain of this operator.
  //! This will be <tt>null</tt> until fillComplete() is called.
  const RCP<const MapClass > getDomainMap() const
  {
    return domainmaps_->getFullMap();
  }

  //! \brief Returns the Map associated with the i'th block domain of this operator.
  //! This will be <tt>null</tt> until fillComplete() is called.
  const RCP<const MapClass > getDomainMap(size_t i) const
  {
    return domainmaps_->getMap(i);
  }

  //! Returns the Map associated with the full range of this operator.
  //! This will be <tt>null</tt> until fillComplete() is called.
  const RCP<const MapClass > getRangeMap() const
  {
    return rangemaps_->getFullMap();
  }

  //! Returns the Map associated with the i'th block range of this operator.
  //! This will be <tt>null</tt> until fillComplete() is called.
  const RCP<const MapClass > getRangeMap(size_t i) const
  {
    return rangemaps_->getMap(i);
  }

  //! Returns map extractor class for range map
  const RCP<const MapExtractorClass> getRangeMapExtractor()
  {
    return rangemaps_;
  }

  //! Returns map extractor for domain map
  const RCP<const MapExtractorClass> getDomainMapExtractor()
  {
    return domainmaps_;
  }

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const {
    std::ostringstream oss;
    oss << "Xpetra_BlockedCrsOperator.description()" << std::endl;
    return oss.str();
  }

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const {
    //     Teuchos::EVerbosityLevel vl = verbLevel;
    //     if (vl == VERB_DEFAULT) vl = VERB_LOW;
    //     RCP<const Comm<int> > comm = this->getComm();
    //     const int myImageID = comm->getRank(),
    //       numImages = comm->getSize();

    //     if (myImageID == 0) out << this->description() << std::endl;

    out << "Xpetra::BlockedCrsOperator: " << Rows() << " x " << Cols() << std::endl;

    if(isFillComplete())
    {
      out << "BlockOperator is filled" << std::endl;
      out << "fullRowMap" << std::endl;
      fullrowmap_->describe(out,verbLevel);
      out << "fullColMap" << std::endl;
      fullcolmap_->describe(out,verbLevel);
    }
    else
      out << "BlockOperator is NOT filled" << std::endl;

    for (size_t r=0; r<Rows(); ++r)
    {
      for(size_t c=0; c<Cols(); ++c)
      {
        out << "Block(" << r << "," << c << ")" << std::endl;
        getMatrix(r,c)->describe(out,verbLevel);
      }
    }


    //matrixData_->describe(out,verbLevel);

    // Teuchos::OSTab tab(out);
  }

  // JG: Added:

  //! Returns the CrsGraph associated with this matrix.
  RCP<const CrsGraph> getCrsGraph() const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
        "getCrsGraph() not supported by BlockedCrsOperator!" );
    return Teuchos::null;
  }

  //@}

  //! @name Block matrix access
  //@{

  /// number of row blocks
  virtual size_t Rows() const { return rangemaps_->NumMaps(); }

  /// number of column blocks
  virtual size_t Cols() const { return domainmaps_->NumMaps(); }

  /// return block (r,c)
  Teuchos::RCP<CrsMatrixClass> getMatrix(size_t r, size_t c) const { return blocks_[r*Cols()+c]; }

  /// set matrix block
  void setMatrix(size_t r, size_t c, Teuchos::RCP<CrsMatrixClass>& mat)
  {
    // TODO: if filled -> return error

    TEUCHOS_TEST_FOR_EXCEPTION( r > Rows(), std::out_of_range, "Error, r = " << Rows() << " is too big" );
    TEUCHOS_TEST_FOR_EXCEPTION( c > Cols(), std::out_of_range, "Error, c = " << Cols() << " is too big" );

    // check row map
    //if (!rangemaps_->Map(r)->isSameAs(mat->getRowMap()))
    //  TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError, "Error. row maps do not fit." );

    // set matrix
    blocks_[r*Cols()+c] = mat;
  }

  /// merge BlockedCrsOperator blocks in a CrsMatrix
  /*
   * This is a rather expensive operation, since all blocks are copied into a new big CrsMatrix
   */
  Teuchos::RCP<CrsMatrixClass> Merge() const
  {
    Teuchos::RCP<CrsMatrixClass> sparse = Teuchos::rcp_dynamic_cast<CrsMatrixClass>(Xpetra::CrsMatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal>::Build(fullrowmap_,33));// Teuchos::rcp(new CrsMatrixClass(*fullrowmap_,33));
    for (size_t i=0; i<blocks_.size(); ++i)
    {
      Teuchos::RCP<CrsMatrixClass> block = Teuchos::rcp_dynamic_cast<CrsMatrixClass>(blocks_[i]);
      this->Add(block,ScalarTraits< Scalar >::one(),sparse,ScalarTraits< Scalar >::one());
    }
    sparse->fillComplete(getDomainMap(),getRangeMap());
    return sparse;
  }
  //@}


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
   * This routine is private and used only by Merge. Since the blocks in BlockedCrsOperator are seperated,
   * this routine works for merging a BlockedCrsOperator.
   */
  void Add(Teuchos::RCP<CrsMatrixClass>& A, const Scalar scalarA, Teuchos::RCP<CrsMatrixClass>& B, const Scalar scalarB) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION( !A->isFillComplete(), Xpetra::Exceptions::RuntimeError,
        "Matrix A is not completed" );

    B->scale(scalarB);

    size_t MaxNumEntries = std::max(A->getNodeMaxNumRowEntries(),B->getNodeMaxNumRowEntries());
    std::vector<GlobalOrdinal> vecIndices(MaxNumEntries);
    std::vector<Scalar>        vecValues (MaxNumEntries);

    const Teuchos::ArrayView<GlobalOrdinal> Indices(&vecIndices[0], vecIndices.size());
    const Teuchos::ArrayView<Scalar>        Values (&vecValues[0],  vecValues.size());
    size_t NumEntries;

    Teuchos::ArrayView<const GlobalOrdinal> MyGlobalRowIds = A->getRowMap()->getNodeElementList(); // global row ids

    if(scalarA)
    {
      for(size_t i=0; i<A->getNodeNumRows(); ++i)
      {
        GlobalOrdinal Row = MyGlobalRowIds[i];

        //A->getGlobalRowCopy(Row, Indices, Values, NumEntries);
        A->getLocalRowCopy(i, Indices, Values, NumEntries);

        if(scalarA != ScalarTraits< Scalar >::one())
          for (size_t j=0; j<NumEntries; ++j)
            Values[j] *= scalarA;


#if 0
        for (size_t j=0; j<NumEntries; ++j)
        {
          std::vector<GlobalOrdinal> tempVec; tempVec.push_back(Indices[j]);
          std::vector<Scalar> tempVal; tempVal.push_back(Values[j]);
          Teuchos::ArrayView<GlobalOrdinal> tempIndex(&tempVec[0], 1);
          Teuchos::ArrayView<Scalar>        tempValue(&tempVal[0],  1);
          B->insertGlobalValues(Row, tempIndex, tempValue); // insert should be ok, since blocks in BlockedCrsOpeartor do not overlap!
        }
#else

        std::vector<GlobalOrdinal> tempvecIndices(NumEntries);
        std::copy(Indices.getRawPtr(),
            Indices.getRawPtr()+NumEntries,
            std::inserter(tempvecIndices,tempvecIndices.begin()));
        std::vector<Scalar> tempvecValues(NumEntries);
        std::copy(Values.getRawPtr(),
            Values.getRawPtr()+NumEntries,
            std::inserter(tempvecValues,tempvecValues.begin()));
        Teuchos::ArrayView<GlobalOrdinal> tempIndex(&tempvecIndices[0], NumEntries);
        Teuchos::ArrayView<Scalar>        tempValue(&tempvecValues[0],  NumEntries);
        B->insertGlobalValues(Row, tempIndex, tempValue); // insert should be ok, since blocks in BlockedCrsOpeartor do not overlap!
#endif


      }
    }
  }

  //@}

  // Default view is created after fillComplete()
  // Because ColMap might not be available before fillComplete().
  void CreateDefaultView() {

    // Create default view
    this->defaultViewLabel_ = "point";
    this->CreateView(this->GetDefaultViewLabel(), getRangeMap(), getDomainMap());

    // Set current view
    this->currentViewLabel_ = this->GetDefaultViewLabel();
  }

private:
  // Teuchos::RCP<CrsMatrixClass> matrixData_;

  /// the full domain map together with all partial domain maps
  Teuchos::RCP<const MapExtractorClass> domainmaps_;

  /// the full range map together with all partial domain maps
  Teuchos::RCP<const MapExtractorClass> rangemaps_;

  /// row major matrix block storage
  std::vector<Teuchos::RCP<CrsMatrixClass> > blocks_;

  /// full matrix row map
  Teuchos::RCP<MapClass> fullrowmap_;

  /// full matrix column map
  Teuchos::RCP<MapClass> fullcolmap_;


}; //class BlockedCrsOperator

} //namespace Xpetra

#define XPETRA_BLOCKEDCRSOPERATOR_SHORT
#endif /* XPETRA_BLOCKEDCRSOPERATOR_HPP_ */
