//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef TPETRA_CRSMATRIX_HPP
#define TPETRA_CRSMATRIX_HPP

// TODO: add support for alpha,beta term coefficients: Y = alpha*A*X + beta*Y
// TODO: row-wise insertion of entries in globalAssemble() may be more efficient

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_CrsMatrix.hpp>
#include <Kokkos_DefaultSparseMultiply.hpp>
#include <Kokkos_DefaultSparseSolve.hpp>
#include <Kokkos_NodeHelpers.hpp>

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include "Tpetra_InverseOperator.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_NodeOps.hpp"

namespace Tpetra {
  // struct for i,j,v triplets
  template <class Ordinal, class Scalar>
  struct CrsIJV {
    CrsIJV() {}
    CrsIJV(Ordinal row, Ordinal col, const Scalar &val) {
      i = row;
      j = col;
      v = val;
    }
    Ordinal i,j;
    Scalar  v;
  };
}

namespace Teuchos {
  // SerializationTraits specialization for CrsIJV, using DirectSerialization
  template<typename Ordinal, typename Scalar>
  class SerializationTraits<int,Tpetra::CrsIJV<Ordinal,Scalar> >
  : public DirectSerializationTraits<int,Tpetra::CrsIJV<Ordinal,Scalar> >
  {};
}

namespace std {
  template <class Ordinal, class Scalar>
  bool operator<(const Tpetra::CrsIJV<Ordinal,Scalar> &ijv1, const Tpetra::CrsIJV<Ordinal,Scalar> &ijv2) {
    return ijv1.i < ijv2.i;
  }
}

namespace Tpetra 
{

  //! \brief A class for constructing and using sparse compressed matrices with row access.
  /*!
   * This class allows the construction of sparse matrices with row-access. 
   * Method insertGlobalValues() can be used to set both locally
   * owned and non-local elements; the shipping of data is done with hardcoded
   * MPI calls when fillComplete() is called.
   *
   * The nonzero elements of  locally owned row can be accessed by method
   * getLocalRowCopy() or getGlobalRowCopy(). The former returns the column
   * indices using local numbering, the latter using global numbering.
   *
   * This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
   * The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
   * type, if omitted, defaults to the \c LocalOrdinal type.
   * The class utilizes CrsGraph object which has the same local and global ordinal types.
   *
   */
  template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatVec = Kokkos::DefaultSparseMultiply<Scalar,LocalOrdinal,Node>, class LocalMatSolve = Kokkos::DefaultSparseSolve<Scalar,LocalOrdinal,Node> >
  class CrsMatrix : public RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
    public:
      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor specifying the number of non-zeros for all rows.
      CrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

      //! Constructor specifying the number of non-zeros for each row.
      CrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const Teuchos::ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile);

      //! Constructor specifying a column map and the number of non-zeros for all rows.
      CrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

      //! Constructor specifying a column map and the number of non-zeros for each row.
      CrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, const Teuchos::ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile);

      //! Constructor specifying a pre-constructed graph.
      explicit CrsMatrix(const Teuchos::RCP<CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > &graph);

      // !Destructor.
      virtual ~CrsMatrix();

      //@}

      //! @name Insertion/Removal Methods
      //@{ 

      //! Submit matrix entries, using global IDs.
      void insertGlobalValues(GlobalOrdinal globalRow, const Teuchos::ArrayView<const GlobalOrdinal> &cols, const Teuchos::ArrayView<const Scalar> &vals);

      //! Submit matrix entries, using local IDs.
      void insertLocalValues(LocalOrdinal localRow, const Teuchos::ArrayView<const LocalOrdinal> &cols, const Teuchos::ArrayView<const Scalar> &vals);

      //! Replace matrix entries, using global IDs.
      /*! All index values must be in the global space. 
         If (globalRow,cols[i]) corresponds to an entry 
         that is duplicated in this matrix (likely because it was inserted more than once and fillComplete() 
         has not been called), the behavior of this function is not defined. */
      void replaceGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &cols,
                         const Teuchos::ArrayView<const Scalar>        &vals);

      //! Sum into multiple entries, using global IDs.
      /*! All index values must be in the global space. */
      void sumIntoGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &cols,
                         const Teuchos::ArrayView<const Scalar>        &vals);


      //! Set all matrix entries equal to scalarThis.
      void setAllToScalar(const Scalar &alpha);

      //! Scale the current values of a matrix, this = alpha*this. 
      void scale(const Scalar &alpha);

      //@}

      //! @name Transformational Methods
      //@{ 

      //! \brief Communicate non-local contributions to other nodes.
      void globalAssemble();

      /*! \brief Signal that data entry is complete, specifying domain and range maps. 
          Off-node entries are distributed, repeated entries are summed, and global indices are transformed to local indices.
          If \c OptimizeStorage is true, then optimizeStorage() is called as well.
       */
      void fillComplete(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, OptimizeOption os = DoOptimizeStorage);

      /*! \brief Signal that data entry is complete. 
          Off-node entries are distributed, repeated entries are summed, and global indices are transformed to local indices.
          If \c OptimizeStorage is true, then optimizeStorage() is called as well.
          \note This method calls fillComplete( getRowMap(), getRowMap() ).
       */
      void fillComplete(OptimizeOption os = DoOptimizeStorage);

      //! \brief Re-allocate the data into contiguous storage.
      void optimizeStorage();

      //@}

      //! @name Methods implementing RowMatrix
      //@{ 

      //! Returns the communicator.
      const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const;

      //! Returns the underlying node.
      Teuchos::RCP<Node> getNode() const;

      //! Returns the Map that describes the row distribution in this matrix.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRowMap() const;

      //! \brief Returns the Map that describes the column distribution in this matrix.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getColMap() const;

      //! Returns the RowGraph associated with this matrix. 
      Teuchos::RCP<const RowGraph<LocalOrdinal,GlobalOrdinal,Node> > getGraph() const;

      //! Returns the CrsGraph associated with this matrix. 
      Teuchos::RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal> > getCrsGraph() const;

      //! Returns the number of global rows in this matrix.
      global_size_t getGlobalNumRows() const;

      //! \brief Returns the number of global columns in this matrix.
      global_size_t getGlobalNumCols() const;

      //! Returns the number of matrix rows owned by the calling image. 
      size_t getNodeNumRows() const;

      //! Returns the number of columns needed to apply the forward operator on this node, i.e., the number of elements listed in the column map.
      size_t getNodeNumCols() const;

      //! Returns the index base for global indices for this matrix. 
      GlobalOrdinal getIndexBase() const;

      //! Returns the global number of entries in this matrix.
      global_size_t getGlobalNumEntries() const;

      //! Returns the local number of entries in this matrix.
      size_t getNodeNumEntries() const;

      //! \brief Returns the current number of entries on this node in the specified global row.
      /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this graph. */
      size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const;

      //! Returns the current number of entries on this node in the specified local row.
      /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this graph. */
      size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const;

      //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
      global_size_t getGlobalNumDiags() const;

      //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
      size_t getNodeNumDiags() const;

      //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
      size_t getGlobalMaxNumRowEntries() const;

      //! \brief Returns the maximum number of entries across all rows/columns on this node.
      size_t getNodeMaxNumRowEntries() const;

      //! \brief Indicates whether this matrix has a well-defined column map. 
      bool hasColMap() const; 

      //! \brief Indicates whether this matrix is lower triangular.
      bool isLowerTriangular() const;

      //! \brief Indicates whether this matrix is upper triangular.
      bool isUpperTriangular() const;

      //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
      bool isLocallyIndexed() const;

      //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
      bool isGloballyIndexed() const;

      //! Returns \c true if fillComplete() has been called.
      bool isFillComplete() const;

      //! Extract a list of entries in a specified global row of this matrix. Put into pre-allocated storage.
      /*!
        \param LocalRow - (In) Global row number for which indices are desired.
        \param Indices - (Out) Global column indices corresponding to values.
        \param Values - (Out) Matrix values.
        \param NumEntries - (Out) Number of indices.

         Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
         with row \c GlobalRow. If \c GlobalRow does not belong to this node, then \c Indices and \c Values are unchanged and \c NumIndices is 
         returned as Teuchos::OrdinalTraits<size_t>::invalid().
       */
      void getGlobalRowCopy(GlobalOrdinal GlobalRow,
                            const Teuchos::ArrayView<GlobalOrdinal> &Indices,
                            const Teuchos::ArrayView<Scalar> &Values,
                            size_t &NumEntries) const;

      //! Extract a list of entries in a specified local row of the graph. Put into storage allocated by calling routine.
      /*!
        \param LocalRow - (In) Local row number for which indices are desired.
        \param Indices - (Out) Local column indices corresponding to values.
        \param Values - (Out) Matrix values.
        \param NumIndices - (Out) Number of indices.

         Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
         with row \c LocalRow. If \c LocalRow is not valid for this node, then \c Indices and \c Values are unchanged and \c NumIndices is 
         returned as Teuchos::OrdinalTraits<size_t>::invalid().
       */
      void getLocalRowCopy(LocalOrdinal LocalRow, 
                           const Teuchos::ArrayView<LocalOrdinal> &Indices, 
                           const Teuchos::ArrayView<Scalar> &Values,
                           size_t &NumEntries) const;

      //! Get a persisting const view of the entries in a specified global row of this matrix.
      /*!
        \param GlobalRow - (In) Global row from which to retrieve matrix entries.
        \param Indices - (Out) Indices for the global row.
        \param Values - (Out) Values for the global row.

         Note: If \c GlobalRow does not belong to this node, then \c Indices and \c Values are set to <tt>Teuchos::null</t>>.

        \pre isLocallyIndexed()==false
       */
      void getGlobalRowView(GlobalOrdinal GlobalRow, 
                            Teuchos::ArrayRCP<const GlobalOrdinal> &indices,
                            Teuchos::ArrayRCP<const Scalar>        &values) const;

      //! Get a persisting const view of the entries in a specified local row of this matrix.
      /*!
        \param LocalRow - (In) Local row from which to retrieve matrix entries.
        \param Indices - (Out) Indices for the local row.
        \param Values - (Out) Values for the local row.

         Note: If \c LocalRow is not valid for this node, then \c Indices and \c Values are set to <tt>Teuchos::null</tt>.

        \pre isGloballyIndexed()==false
       */
      void getLocalRowView(LocalOrdinal LocalRow,
                           Teuchos::ArrayRCP<const LocalOrdinal> &indices,
                           Teuchos::ArrayRCP<const Scalar>       &values) const;

      //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
      /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the 
          the zero and non-zero diagonals owned by this node. */
      void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal> &diag) const;

      //@}

      //! @name Matrix-vector multiplication and solve methods
      //@{

      //! Multiplies this matrix by a MultiVector.
      template <class DomainScalar, class RangeScalar>
      void multiply(const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans) const;

      //! Solves a linear system when the underlying matrix is triangular.
      template <class DomainScalar, class RangeScalar>
      void solve(const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> & Y, MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X, Teuchos::ETransp trans) const;
          
      //@}

      //! @name Methods implementing Operator
      //@{ 

      //! \brief Returns the Map associated with the domain of this operator.
      //! This will be equal to the row map until fillComplete() is called.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getOperatorDomainMap() const;

      //! Returns the Map associated with the domain of this operator.
      //! This will be equal to the row map until fillComplete() is called.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getOperatorRangeMap() const;

      //! Computes this matrix-vector multilication y = A x.
      void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

      //! Indicates whether this operator supports applying the adjoint operator.
      bool hasTransposeApply() const;

      //@}

      //! @name Methods implementing InverseOperator
      //@{ 

      //! \brief Returns the Map associated with the domain of this operator.
      //! This will be equal to the row map until fillComplete() is called.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getInverseOperatorDomainMap() const;

      //! Returns the Map associated with the domain of this operator.
      //! This will be equal to the row map until fillComplete() is called.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getInverseOperatorRangeMap() const;

      //! Computes this matrix-vector multilication y = A x.
      void applyInverse(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

      //! Indicates whether this operator supports inverting the adjoint operator.
      bool hasTransposeApplyInverse() const;

      //@}

      //! @name Miscellaneous Query Methods
      //@{

      //! Returns \c true if optimizeStorage() has been called.
      bool isStorageOptimized() const;

      //! Returns \c true if the graph data was allocated in static data structures.
      ProfileType getProfileType() const;

      //! Indicates that the graph is static, so that new entries cannot be added to this matrix. */
      bool isStaticGraph() const;

      //@} 

      //! @name Overridden from Teuchos::Describable 
      //@{

      /** \brief Return a simple one-line description of this object. */
      std::string description() const;

      /** \brief Print the object with some verbosity level to an FancyOStream object. */
      void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

      //@}

    private:
      // copy constructor disabled
      CrsMatrix(const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> &Source);
      // operator= disabled
      CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> & operator=(const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> &rhs);
    protected:

      // useful typedefs
      typedef Teuchos::OrdinalTraits<LocalOrdinal>    LOT;
      typedef Teuchos::OrdinalTraits<GlobalOrdinal>   GOT;
      typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

      void allocateValues(Scalar alpha = Teuchos::ScalarTraits<Scalar>::zero());
      void sortEntries();
      void mergeRedundantEntries();
      void checkInternalState() const;
      void updateAllocation(size_t lrow, size_t allocSize);
      void fillLocalMatrix();

      //! \brief Get a persisting const view of the elements in a specified local row of the matrix.
      /*! This protected method is used internally for almost all access to the matrix elements.
          No checking is done except in a debug build.

        \param myRow      - (In) \c size_t specifying the local row.
        \param sizeInfo   - (In) Size info from CrsGraph::getRowInfo()

        \returns values   - (Out) persisting, const view of the local values.
       */
      Teuchos::ArrayRCP<const Scalar> getFullView(size_t myRow, RowInfo sizeInfo) const;

      //! \brief Get a persisting non-const view of the elements in a specified local row of the matrix.
      /*! This protected method is used internally for almost all access to the matrix elements.
          No checking is done except in a debug build.

        \param myRow      - (In) \c size_t specifying the local row.
        \param sizeInfo   - (In) Size info from CrsGraph::getRowInfo()

        \returns values   - (Out) persisting, non-const view of the local values.
       */
      Teuchos::ArrayRCP<Scalar> getFullViewNonConst(size_t myRow, RowInfo sizeInfo);

      Teuchos::RCP<CrsGraph<LocalOrdinal,GlobalOrdinal> > graph_;

      Kokkos::CrsMatrix<Scalar,Node> lclMatrix_;
      LocalMatVec lclMatVec_;
      LocalMatSolve lclMatSolve_;

      bool valuesAreAllocated_,
           staticGraph_,
           constructedWithFilledGraph_,
           fillComplete_,
           storageOptimized_;

      // matrix values. before allocation, these are Teuchos::null.
      // after allocation, one is Teuchos::Null.
      // these are parallel compute buffers, not host memory. therefore, in general, they cannot be dereferenced in host cost.
      // 1D == StaticAllocation, 2D == DynamicAllocation
      // The allocation always matches that of graph_
      Teuchos::ArrayRCP<Scalar>                       pbuf_values1D_;
      Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >   pbuf_values2D_;

      // multivectors used for import/export dest/source in apply()
      mutable Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > importMV_, exportMV_;

      // a map between a (non-local) row and a list of (col,val)
      std::map<GlobalOrdinal, std::list<std::pair<GlobalOrdinal,Scalar> > > nonlocals_;

  }; // class CrsMatrix



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrix(
                                          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                          size_t maxNumEntriesPerRow, 
                                          ProfileType pftype)
  : lclMatrix_(rowMap->getNodeNumElements(), rowMap->getNode())
  , lclMatVec_(rowMap->getNode())
  , lclMatSolve_(rowMap->getNode())
  , valuesAreAllocated_(false)
  , staticGraph_(false)
  , constructedWithFilledGraph_(false)
  , fillComplete_(false)
  , storageOptimized_(false) {
    try {
      graph_ = Teuchos::rcp( new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowMap,maxNumEntriesPerRow,pftype) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          Teuchos::typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    checkInternalState();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrix(
                                          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                          const Teuchos::ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, 
                                          ProfileType pftype)
  : lclMatrix_(rowMap->getNodeNumElements(), rowMap->getNode())
  , lclMatVec_(rowMap->getNode())
  , lclMatSolve_(rowMap->getNode())
  , valuesAreAllocated_(false)
  , staticGraph_(false)
  , constructedWithFilledGraph_(false)
  , fillComplete_(false)
  , storageOptimized_(false) {
    try {
      graph_ = Teuchos::rcp( new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowMap,NumEntriesPerRowToAlloc,pftype) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          Teuchos::typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    checkInternalState();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrix(
                                          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, 
                                          size_t maxNumEntriesPerRow, 
                                          ProfileType pftype)
  : lclMatrix_(rowMap->getNodeNumElements(), rowMap->getNode())
  , lclMatVec_(rowMap->getNode())
  , lclMatSolve_(rowMap->getNode())
  , valuesAreAllocated_(false)
  , staticGraph_(false)
  , constructedWithFilledGraph_(false)
  , fillComplete_(false)
  , storageOptimized_(false) {
    try {
      graph_ = Teuchos::rcp( new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowMap,colMap,maxNumEntriesPerRow,pftype) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          Teuchos::typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    checkInternalState();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrix(
                                          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, 
                                          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, 
                                          const Teuchos::ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, 
                                          ProfileType pftype)
  : lclMatrix_(rowMap->getNodeNumElements(), rowMap->getNode())
  , lclMatVec_(rowMap->getNode())
  , lclMatSolve_(rowMap->getNode())
  , valuesAreAllocated_(false)
  , staticGraph_(false)
  , constructedWithFilledGraph_(false)
  , fillComplete_(false)
  , storageOptimized_(false) {
    try {
      graph_ = Teuchos::rcp( new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowMap,NumEntriesPerRowToAlloc,pftype) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          Teuchos::typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    checkInternalState();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrix(const Teuchos::RCP< CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > &graph)
  : graph_(graph)
  , lclMatrix_(graph->getRowMap()->getNodeNumElements(), graph->getRowMap()->getNode())
  , lclMatVec_(graph->getNode())
  , lclMatSolve_(graph->getNode())
  , valuesAreAllocated_(false)
  , staticGraph_(true)
  , fillComplete_(graph->isFillComplete())
  , storageOptimized_(graph->isStorageOptimized()) {
    TEST_FOR_EXCEPTION(graph_ == Teuchos::null, std::runtime_error,
        Teuchos::typeName(*this) << "::CrsMatrix(graph): specified pointer is null.");
    // we won't prohibit the case where the graph is not yet filled, but we will check below to ensure that the
    // graph isn't filled between now and when fillComplete() is called on this matrix.
    // the user is not allowed to modify the graph. this check is about the most that we can do to check whether they have.
    constructedWithFilledGraph_ = graph_->isFillComplete();
    checkInternalState();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::~CrsMatrix() {
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Teuchos::Comm<int> > &
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getComm() const {
    return graph_->getComm();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  Teuchos::RCP<Node>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNode() const {
    return lclMatVec_.getNode();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  ProfileType CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getProfileType() const {
    return graph_->getProfileType();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isFillComplete() const {
    return fillComplete_; 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isStorageOptimized() const {
    return storageOptimized_; 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isLocallyIndexed() const {
    return graph_->isLocallyIndexed();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isGloballyIndexed() const {
    return graph_->isGloballyIndexed();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::hasColMap() const {
    return graph_->hasColMap();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  global_size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalNumEntries() const {
    return graph_->getGlobalNumEntries();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNodeNumEntries() const {
    return graph_->getNodeNumEntries();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  global_size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalNumRows() const {
    return graph_->getGlobalNumRows(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  global_size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalNumCols() const { 
    return graph_->getGlobalNumCols(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNodeNumRows() const { 
    return graph_->getNodeNumRows(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNodeNumCols() const { 
    return graph_->getNodeNumCols(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  global_size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalNumDiags() const { 
    return graph_->getGlobalNumDiags(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNodeNumDiags() const { 
    return graph_->getNodeNumDiags(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const { 
    return graph_->getNumEntriesInGlobalRow(globalRow); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNumEntriesInLocalRow(LocalOrdinal localRow) const { 
    return graph_->getNumEntriesInLocalRow(localRow);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalMaxNumRowEntries() const { 
    return graph_->getGlobalMaxNumRowEntries(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNodeMaxNumRowEntries() const { 
    return graph_->getNodeMaxNumRowEntries(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getIndexBase() const { 
    return getRowMap()->getIndexBase(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getRowMap() const { 
    return graph_->getRowMap(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getColMap() const {
    return graph_->getColMap(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getOperatorDomainMap() const { 
    return graph_->getDomainMap(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getOperatorRangeMap() const { 
    return graph_->getRangeMap(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getInverseOperatorDomainMap() const { 
    return graph_->getRangeMap(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getInverseOperatorRangeMap() const { 
    return graph_->getDomainMap(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  Teuchos::RCP<const RowGraph<LocalOrdinal,GlobalOrdinal,Node> >
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGraph() const { 
    return graph_; 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  Teuchos::RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal> >
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getCrsGraph() const { 
    return graph_; 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isLowerTriangular() const { 
    return graph_->isLowerTriangular(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isUpperTriangular() const { 
    return graph_->isUpperTriangular(); 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isStaticGraph() const { 
    return staticGraph_; 
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::hasTransposeApply() const {
    return true;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::hasTransposeApplyInverse() const {
    return true;
  }
      

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                    Internal utility methods                             //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::allocateValues(Scalar alpha) {
#ifdef HAVE_TPETRA_DEBUG
    // need graph_->getNodeAllocationSize(), but it is only available after the graph indices are allocated. we will 
    // internally enforce that the graph is allocated before allocateValues() is called. this test serves to make sure
    // that we have done so.
    TEST_FOR_EXCEPTION( graph_->indicesAreAllocated() == false, std::logic_error,
        Teuchos::typeName(*this) << "::allocateValues(): Internal logic error. Please contact Tpetra team.");
#endif
    const size_t nlrs = getRowMap()->getNodeNumElements(),
                  nta = graph_->getNodeAllocationSize(),
           numEntries = graph_->getNodeNumEntries();
    if (valuesAreAllocated_) {
      return;
    }
    // do we even have anything to allocate?
    if (nta > 0) {
      Teuchos::RCP<Node> node = lclMatrix_.getNode();
      Kokkos::ReadyBufferHelper<Node> rbh(node);
      ////////////////////////////////////////
      if (getProfileType() == StaticProfile) {
        // allocate values
        pbuf_values1D_ = node->template allocBuffer<Scalar>(nta);
        // init values in parallel, if the graph already has valid entries
        if (numEntries > 0) {
          InitOp<Scalar> wdp;
          wdp.val = alpha;
          rbh.begin();
          wdp.x   = rbh.addNonConstBuffer(pbuf_values1D_);
          rbh.end();
          node->template parallel_for<InitOp<Scalar> >(0,nta,wdp);
        }
      }
      ///////////////////////////////////////////////
      else { // if getProfileType() == DynamicProfile
        InitOp<Scalar> wdp;
        wdp.val = alpha;
        // allocate array of buffers
        pbuf_values2D_ = Teuchos::arcp< Teuchos::ArrayRCP<Scalar> >(nlrs);
        for (size_t r=0; r<nlrs; ++r) {
          // this call to getNumAllocatedEntries() is cheap for the DynamicProfile case
          const size_t ntarow = graph_->getNumAllocatedEntriesInLocalRow(r),
                        nErow = graph_->getNumEntriesInLocalRow(r);
          if (ntarow > 0) {
            // allocate values for this row
            pbuf_values2D_[r] = node->template allocBuffer<Scalar>(ntarow);
            // initi values in parallel, if the graph already has valid entries
            if (nErow > 0) {
              rbh.begin();
              wdp.x   = rbh.addNonConstBuffer(pbuf_values2D_[r]);
              rbh.end();
              node->template parallel_for<InitOp<Scalar> >(0,ntarow,wdp);
            }
          }
        }
      }
    }
    valuesAreAllocated_ = true;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::checkInternalState() const {
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::RCP<Node> node = getNode();
    const global_size_t gsti = Teuchos::OrdinalTraits<global_size_t>::invalid();
    using Teuchos::null;
    std::string err = Teuchos::typeName(*this) + "::checkInternalState(): Likely internal logic error. Please contact Tpetra team.";
    // check the internal state of this data structure
    // this is called by numerous state-changing methods, in a debug build, to ensure that the object 
    // always remains in a valid state

    // constructedWithFilledGraph_ should only be true if matrix was constructed with a graph, in which case staticGraph_ should be true
    TEST_FOR_EXCEPTION( staticGraph_ == false && constructedWithFilledGraph_ == true,                  std::logic_error, err ); 
    // matrix values cannot be allocated without the graph indices being allocated
    TEST_FOR_EXCEPTION( valuesAreAllocated_ == true && graph_->indicesAreAllocated() == false,         std::logic_error, err );
    // if matrix is fill complete, then graph must be fill complete
    TEST_FOR_EXCEPTION( fillComplete_ == true && graph_->isFillComplete() == false,                    std::logic_error, err );
    // if matrix is storage optimized, then graph must be storage optimizied
    TEST_FOR_EXCEPTION( storageOptimized_   != graph_->isStorageOptimized(),                           std::logic_error, err );
    // if matrix is storage optimized, it must have been fill completed
    TEST_FOR_EXCEPTION( storageOptimized_ == true && fillComplete_ == false,                           std::logic_error, err );
    // if matrix is storage optimized, it should have a 1D allocation 
    TEST_FOR_EXCEPTION( storageOptimized_ == true && pbuf_values2D_ != Teuchos::null,                  std::logic_error, err );
    // if matrix/graph are static profile, then 2D allocation should not be present
    TEST_FOR_EXCEPTION( graph_->getProfileType() == StaticProfile  && pbuf_values2D_ != Teuchos::null, std::logic_error, err );
    // if matrix/graph are dynamic profile, then 1D allocation should not be present
    TEST_FOR_EXCEPTION( graph_->getProfileType() == DynamicProfile && pbuf_values1D_ != Teuchos::null, std::logic_error, err );
    // if values are allocated and they are non-zero in number, then one of the allocations should be present
    TEST_FOR_EXCEPTION( valuesAreAllocated_ == true && graph_->getNodeAllocationSize() && pbuf_values2D_ == Teuchos::null && pbuf_values1D_ == Teuchos::null,
                        std::logic_error, err );
    // we can nae have both a 1D and 2D allocation
    TEST_FOR_EXCEPTION( pbuf_values1D_ != Teuchos::null && pbuf_values2D_ != Teuchos::null, std::logic_error, err );
    // compare matrix allocations against graph allocations
    if (valuesAreAllocated_ && graph_->getNodeAllocationSize() > 0) {
      if (graph_->getProfileType() == StaticProfile) {
        if (graph_->isLocallyIndexed()) {
          TEST_FOR_EXCEPTION( pbuf_values1D_.size() != graph_->pbuf_lclInds1D_.size(),  std::logic_error, err );
        }
        else {
          TEST_FOR_EXCEPTION( pbuf_values1D_.size() != graph_->pbuf_gblInds1D_.size(),  std::logic_error, err );
        }
      } 
      else { // graph_->getProfileType() == DynamicProfile
        if (graph_->isLocallyIndexed()) {
          for (size_t r=0; r < getNodeNumRows(); ++r) {
            TEST_FOR_EXCEPTION( (pbuf_values2D_[r] == Teuchos::null) != (graph_->pbuf_lclInds2D_[r] == Teuchos::null), std::logic_error, err );
            if (pbuf_values2D_[r] != Teuchos::null && graph_->pbuf_lclInds2D_[r] != Teuchos::null) {
              TEST_FOR_EXCEPTION( pbuf_values2D_[r].size() != graph_->pbuf_lclInds2D_[r].size(), std::logic_error, err );
            }
          }
        }
        else {
          for (size_t r=0; r < getNodeNumRows(); ++r) {
            TEST_FOR_EXCEPTION( (pbuf_values2D_[r] == Teuchos::null) != (graph_->pbuf_gblInds2D_[r] == Teuchos::null), std::logic_error, err );
            if (pbuf_values2D_[r] != Teuchos::null && graph_->pbuf_gblInds2D_[r] != Teuchos::null) {
              TEST_FOR_EXCEPTION( pbuf_values2D_[r].size() != graph_->pbuf_gblInds2D_[r].size(), std::logic_error, err );
            }
          }
        }
      }
    }
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  Teuchos::ArrayRCP<const Scalar> 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getFullView(size_t myRow, RowInfo sizeInfo) const {
#ifdef HAVE_TPETRA_DEBUG
    std::string err = Teuchos::typeName(*this) + "::getFullView(): Internal logic error. Please contact Tpetra team.";
    TEST_FOR_EXCEPTION(getRowMap()->isNodeLocalElement(myRow) == false, std::logic_error, err);
#endif
    Teuchos::ArrayRCP<const Scalar> values = Teuchos::null;
    // sizeInfo indicates the allocation size for this row, whether it has actually been allocated or not
    if (sizeInfo.allocSize > 0 && valuesAreAllocated_) {
      Teuchos::RCP<Node> node = getNode();
      if (graph_->getProfileType() == StaticProfile) {
        values = node->template viewBuffer<Scalar>(sizeInfo.allocSize, pbuf_values1D_ + sizeInfo.offset1D);
      }
      else {  // dynamic profile
        values = node->template viewBuffer<Scalar>(sizeInfo.allocSize, pbuf_values2D_[myRow]);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(valuesAreAllocated_ == true && values != Teuchos::null && values.size() != sizeInfo.allocSize, std::logic_error, err);
#endif
    return values;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  Teuchos::ArrayRCP<Scalar> 
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getFullViewNonConst(size_t myRow, RowInfo sizeInfo) {
#ifdef HAVE_TPETRA_DEBUG
    std::string err = Teuchos::typeName(*this) + "::getFullViewNonConst(): Internal logic error. Please contact Tpetra team.";
    TEST_FOR_EXCEPTION(getRowMap()->isNodeLocalElement(myRow) == false, std::logic_error, err);
#endif
    Teuchos::ArrayRCP<Scalar> values = Teuchos::null;
    // sizeInfo indicates the allocation size for this row, whether it has actually been allocated or not
    if (sizeInfo.allocSize > 0 && valuesAreAllocated_) {
      Teuchos::RCP<Node> node = getNode();
      // if there are no valid entries, then this view can be constructed WriteOnly
      Kokkos::ReadWriteOption rw = (sizeInfo.numEntries == 0 ? Kokkos::WriteOnly : Kokkos::ReadWrite);
      if (graph_->getProfileType() == StaticProfile) {
        values = node->template viewBufferNonConst<Scalar>(rw, sizeInfo.allocSize, pbuf_values1D_ + sizeInfo.offset1D);
      }
      else {  // dynamic profile
        values = node->template viewBufferNonConst<Scalar>(rw, sizeInfo.allocSize, pbuf_values2D_[myRow]);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(valuesAreAllocated_ == true && values != Teuchos::null && values.size() != sizeInfo.allocSize, std::logic_error, err);
#endif
    return values;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::updateAllocation(size_t lrow, size_t allocSize) {
    using Teuchos::ArrayRCP;
    Teuchos::RCP<Node> node = getNode();
    RowInfo sizeInfo = graph_->getRowInfo(lrow);
    if (pbuf_values2D_ == Teuchos::null) {
      pbuf_values2D_ = Teuchos::arcp< ArrayRCP<Scalar> >(getNodeNumRows());
    }
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPT( getRowMap()->isNodeLocalElement(lrow) == false );
    TEST_FOR_EXCEPT( pbuf_values2D_[lrow] != Teuchos::null && allocSize < pbuf_values2D_[lrow].size() );
    TEST_FOR_EXCEPT( allocSize == 0 );
#endif
    ArrayRCP<Scalar> old_row, new_row;
    old_row = pbuf_values2D_[lrow];
    new_row = node->template allocBuffer<Scalar>(allocSize);
    if (sizeInfo.numEntries) {
      node->template copyBuffers<Scalar>(sizeInfo.numEntries,old_row,new_row);
    }
    old_row = Teuchos::null;
    pbuf_values2D_[lrow] = new_row;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::fillLocalMatrix() {
    lclMatrix_.clear();
    if (storageOptimized_) {
      // fill packed matrix; it is okay for pbuf_values1D_ to be null; the matrix will flag itself as empty
      lclMatrix_.setPackedValues(pbuf_values1D_);
    }
    else if (graph_->getProfileType() == StaticProfile) {
      if (pbuf_values1D_ != Teuchos::null) {
        const size_t nlrs = getNodeNumRows();
        for (size_t r=0; r < nlrs; ++r) {
          RowInfo sizeInfo = graph_->getRowInfo(r);
          Teuchos::ArrayRCP<const Scalar> rowvals;
          if (sizeInfo.numEntries > 0) {
            rowvals = pbuf_values1D_.persistingView(sizeInfo.offset1D, sizeInfo.numEntries);
            lclMatrix_.set2DValues(r,rowvals);
          }
        }
      }
    }
    else if (graph_->getProfileType() == DynamicProfile) {
      if (pbuf_values2D_ != Teuchos::null) {
        const size_t nlrs = getNodeNumRows();
        for (size_t r=0; r < nlrs; ++r) {
          RowInfo sizeInfo = graph_->getRowInfo(r);
          Teuchos::ArrayRCP<const Scalar> rowvals = pbuf_values2D_[r];
          if (sizeInfo.numEntries > 0) {
            rowvals = rowvals.persistingView(0,sizeInfo.numEntries);
            lclMatrix_.set2DValues(r,rowvals);
          }
        }
      }
    }

    // submit local matrix and local graph to lclMatVec_ and lclMatSolve_
    // lclMatVec_ and lclMatSolve_ are permitted to view, but we don't care whether they do or not
    lclMatVec_.clear();
    lclMatSolve_.clear();
    Teuchos::DataAccess ret;
    ret = lclMatVec_.initializeStructure( graph_->lclGraph_, Teuchos::View );
    ret = lclMatVec_.initializeValues( lclMatrix_, Teuchos::View );
    if (isLowerTriangular() || isUpperTriangular()) {
      ret = lclMatSolve_.initializeStructure( graph_->lclGraph_, Teuchos::View );
      ret = lclMatSolve_.initializeValues( lclMatrix_, Teuchos::View );
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                  User-visible class methods                             //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::insertLocalValues(
                                                         LocalOrdinal localRow, 
                         const Teuchos::ArrayView<const LocalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>       &values) {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(isStorageOptimized() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalValues(): cannot insert new values after optimizeStorage() has been called.");
    TEST_FOR_EXCEPTION(graph_->isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalValues(): graph indices are global; use insertGlobalValues().");
    TEST_FOR_EXCEPTION(hasColMap() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalValues(): cannot insert local indices without a column map; ");
    TEST_FOR_EXCEPTION(isStaticGraph() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalValues(): matrix was constructed with static graph; cannot insert new entries.");
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalValues(): values.size() must equal indices.size().");
    TEST_FOR_EXCEPTION(getRowMap()->isNodeLocalElement(localRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::insertLocalValues(): row does not belong to this node.");
    Teuchos::Array<LocalOrdinal> finds;
    Teuchos::Array<Scalar>       fvals;
    finds.reserve(indices.size());
    fvals.reserve(values.size());
    // use column map to filter the entries:
    const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &cmap = *getColMap();
    for (size_t i=0; i < indices.size(); ++i) {
      if (cmap.isNodeLocalElement(indices[i])) {
        finds.push_back(indices[i]);
        fvals.push_back(values[i]);
      }
    }
    if (finds.size() > 0) {
      if (valuesAreAllocated_ == false) {
        // allocate graph, so we can access views below
        if (graph_->indicesAreAllocated() == false || graph_->indicesAreAllocated() == false) {
          graph_->allocateIndices( CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::AllocateLocal );
        }
        // graph must be allocated before matrix, because matrix needs to know the total allocation of the graph
        allocateValues();
#ifdef HAVE_TPETRA_DEBUG
        TEST_FOR_EXCEPTION(valuesAreAllocated_ == false, std::logic_error, 
            Teuchos::typeName(*this) << "::insertLocalValues(): Internal logic error. Please contact Tpetra team.");
#endif
      }
      //
      ArrayRCP<LocalOrdinal> rowindsview;
      ArrayRCP<Scalar>       rowvalsview;
      RowInfo sizeInfo = graph_->getFullLocalViewNonConst(localRow, rowindsview);
      rowvalsview = getFullViewNonConst(localRow, sizeInfo);
      const size_t newSize = sizeInfo.numEntries + finds.size();
      if (newSize > sizeInfo.allocSize) {
        TEST_FOR_EXCEPTION(graph_->getProfileType() == StaticProfile, std::runtime_error,
            Teuchos::typeName(*this) << "::insertLocalValues(): new indices exceed statically allocated graph structure.");
        TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
            "::insertLocalValues(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
        // update allocation only as much as necessary
        graph_->updateLocalAllocation(localRow,newSize);
        updateAllocation(localRow,newSize);
        // get new views; inefficient, but acceptible in this already inefficient case
        sizeInfo = graph_->getFullLocalViewNonConst(localRow, rowindsview);
        rowvalsview = getFullViewNonConst(localRow, sizeInfo);
      }
      // add new indices, values to graph and matrix rows
      graph_->insertLocalIndicesViaView( localRow, finds, rowindsview + sizeInfo.numEntries );
      std::copy( fvals.begin(), fvals.end(), rowvalsview + sizeInfo.numEntries);
#ifdef HAVE_TPETRA_DEBUG
      {
        RowInfo sizeInfoNew = graph_->getRowInfo(localRow);
        TEST_FOR_EXCEPTION(sizeInfoNew.numEntries != newSize, std::logic_error,
            Teuchos::typeName(*this) << "::insertLocalValues(): Internal logic error. Please contact Tpetra team.");
      }
#endif
      // 
      rowindsview = Teuchos::null;
      rowvalsview = Teuchos::null;
      checkInternalState();
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::insertGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>        &values) {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(graph_->isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalValues(): graph indices are local; use insertLocalValues().");
    TEST_FOR_EXCEPTION(isStaticGraph() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalValues(): matrix was constructed with static graph. Cannot insert new entries.");
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalValues(): values.size() must equal indices.size().");
    const LocalOrdinal myRow = getRowMap()->getLocalElement(globalRow);
    if (myRow != LOT::invalid()) {
      // if we have a column map, use it to filter the entries.
      // only filter if this is our row.
      Teuchos::Array<GlobalOrdinal> finds_is_temporary;
      Teuchos::Array<Scalar>        fvals_is_temporary;
      Teuchos::ArrayView<const GlobalOrdinal> findices = indices;
      Teuchos::ArrayView<const Scalar       > fvalues  = values;
      if (hasColMap()) {
        // filter indices and values through the column map
        finds_is_temporary.reserve(indices.size());
        fvals_is_temporary.reserve(values.size());
        const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &cmap = *getColMap();
        for (size_t i=0; i<indices.size(); ++i) {
          if (cmap.isNodeGlobalElement(indices[i])) {
            finds_is_temporary.push_back(indices[i]);
            fvals_is_temporary.push_back(values[i]);
          }
        }
        findices = finds_is_temporary();
        fvalues  = fvals_is_temporary();
      }
      // add the new indices and values
      if (findices.size() > 0) {
        if (valuesAreAllocated_ == false) {
          // allocate graph, so we can access views below
          if (graph_->indicesAreAllocated() == false) {
            graph_->allocateIndices( CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::AllocateGlobal );
          }
          // graph must be allocated before matrix, because matrix needs to know the total allocation of the graph
          allocateValues();
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION(valuesAreAllocated_ == false || graph_->indicesAreAllocated() == false, std::logic_error, 
              Teuchos::typeName(*this) << "::insertGlobalValues(): Internal logic error. Please contact Tpetra team.");
#endif
        }
        // 
        ArrayRCP<GlobalOrdinal> rowindsview;
        ArrayRCP<Scalar>        rowvalsview;
        RowInfo sizeInfo = graph_->getFullGlobalViewNonConst(myRow, rowindsview);
        rowvalsview = getFullViewNonConst(myRow, sizeInfo);
        const size_t newSize = sizeInfo.numEntries + findices.size();
        if (newSize > sizeInfo.allocSize) {
          TEST_FOR_EXCEPTION(graph_->getProfileType() == StaticProfile, std::runtime_error,
              Teuchos::typeName(*this) << "::insertGlobalValues(): new indices exceed statically allocated graph structure.");
          TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
              "::insertGlobalValues(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
          // update allocation only as much as necessary
          graph_->updateGlobalAllocation(myRow,newSize);
          updateAllocation(myRow,newSize);
          // get new views; inefficient, but acceptible in this already inefficient case
          sizeInfo = graph_->getFullGlobalViewNonConst(myRow, rowindsview);
          rowvalsview = getFullViewNonConst(myRow, sizeInfo);
        }
        // add new indices, values to graph and matrix rows
        graph_->insertGlobalIndicesViaView( myRow, findices, rowindsview + sizeInfo.numEntries );
        std::copy(  fvalues.begin(),  fvalues.end(), rowvalsview + sizeInfo.numEntries );
#ifdef HAVE_TPETRA_DEBUG
        {
          RowInfo sizeInfoNew = graph_->getRowInfo(myRow);
          TEST_FOR_EXCEPTION(sizeInfoNew.numEntries != newSize, std::logic_error,
              Teuchos::typeName(*this) << "::insertGlobalValues(): Internal logic error. Please contact Tpetra team.");
        }
#endif
        // 
        rowindsview = Teuchos::null;
        rowvalsview = Teuchos::null;
        checkInternalState();
      }
    }
    else {
      typename Teuchos::ArrayView<const GlobalOrdinal>::iterator ind = indices.begin();
      typename Teuchos::ArrayView<const Scalar       >::iterator val =  values.begin();
      for (; val != values.end(); ++val, ++ind) {
        nonlocals_[globalRow].push_back(std::make_pair(*ind, *val));
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::replaceGlobalValues(      
                                        GlobalOrdinal globalRow, 
                                        const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                                        const Teuchos::ArrayView<const Scalar>        &values) {
    using Teuchos::ArrayRCP;
    // find the values for the specified indices
    // if the row is not ours, throw an exception
    // ignore values not in the matrix (indices not found)
    // operate whether indices are local or global
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid();
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::replaceGlobalValues(): values.size() must equal indices.size().");
    typename Teuchos::ArrayView<const GlobalOrdinal>::iterator ind = indices.begin();
    typename Teuchos::ArrayView<const        Scalar>::iterator val = values.begin();
    LocalOrdinal lrow = getRowMap()->getLocalElement(globalRow);
    TEST_FOR_EXCEPTION(lrow == LOT::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::replaceGlobalValues(): specified global row does not belong to this processor.");
    //
    if (valuesAreAllocated_ == false) {
      // allocate graph, so we can access views below
      if (graph_->indicesAreAllocated() == false) {
        graph_->allocateIndices( CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::AllocateGlobal );
      }
      // graph must be allocated before matrix, because matrix needs to know the total allocation of the graph
      allocateValues();
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(valuesAreAllocated_ == false || graph_->indicesAreAllocated() == false, std::logic_error, 
          Teuchos::typeName(*this) << "::insertLocalValues(): Internal logic error. Please contact Tpetra team.");
#endif
    }
    // 
    if (isLocallyIndexed() == true) {
      ArrayRCP<const LocalOrdinal> lindrowview;
      RowInfo sizeInfo = graph_->getFullLocalView(lrow, lindrowview);
      ArrayRCP<Scalar> valrowview = getFullViewNonConst(lrow, sizeInfo);
      while (ind != indices.end()) {
        LocalOrdinal lind = getColMap()->getLocalElement(*ind);
        size_t loc = graph_->findLocalIndex(lrow,lind,lindrowview);
        if (loc != STINV) {
          valrowview[loc] = (*val);
        }
        ++ind;
        ++val;
      }
      valrowview = Teuchos::null;
      lindrowview = Teuchos::null;
    }
    else if (isGloballyIndexed() == true) {
      ArrayRCP<const GlobalOrdinal> gindrowview;
      RowInfo sizeInfo = graph_->getFullGlobalView(lrow, gindrowview);
      ArrayRCP<Scalar> valrowview = getFullViewNonConst(lrow, sizeInfo);
      while (ind != indices.end()) {
        size_t loc = graph_->findGlobalIndex(lrow,*ind,gindrowview);
        if (loc != STINV) {
          valrowview[loc] = (*val);
        }
        ++ind;
        ++val;
      }
      valrowview = Teuchos::null;
      gindrowview = Teuchos::null;
    }
    //else {
    // graph indices are not allocated, i.e., are non-existant; nothing to do
    //}
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::sumIntoGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>        &values) {
    using Teuchos::ArrayRCP;
    // find the values for the specified indices
    // if the row is not ours, throw an exception
    // ignore values not in the matrix (indices not found)
    // operate whether indices are local or global
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid();
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::sumIntoGlobalValues(): values.size() must equal indices.size().");
    typename Teuchos::ArrayView<const GlobalOrdinal>::iterator ind = indices.begin();
    typename Teuchos::ArrayView<const        Scalar>::iterator val = values.begin();
    LocalOrdinal lrow = getRowMap()->getLocalElement(globalRow);
    TEST_FOR_EXCEPTION(lrow == LOT::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::sumIntoGlobalValues(): specified global row does not belong to this processor.");
    //
    if (valuesAreAllocated_ == false) {
      // allocate graph, so we can access views below
      if (graph_->indicesAreAllocated() == false) {
        graph_->allocateIndices( CrsGraph<LocalOrdinal,GlobalOrdinal,Node>::AllocateGlobal );
      }
      // graph must be allocated before matrix, because matrix needs to know the total allocation of the graph
      allocateValues();
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION(valuesAreAllocated_ == false || graph_->indicesAreAllocated() == false, std::logic_error, 
          Teuchos::typeName(*this) << "::insertLocalValues(): Internal logic error. Please contact Tpetra team.");
#endif
    }
    //
    if (isLocallyIndexed() == true) {
      ArrayRCP<const LocalOrdinal> lindrowview;
      RowInfo sizeInfo = graph_->getFullLocalView(lrow, lindrowview);
      ArrayRCP<Scalar> valrowview = getFullViewNonConst(lrow, sizeInfo);
      while (ind != indices.end()) {
        LocalOrdinal lind = getColMap()->getLocalElement(*ind);
        size_t loc = graph_->findLocalIndex(lrow,lind,lindrowview);
        if (loc != STINV) {
          valrowview[loc] += (*val);
        }
        ++ind;
        ++val;
      }
      valrowview = Teuchos::null;
      lindrowview = Teuchos::null;
    }
    else if (isGloballyIndexed() == true) {
      ArrayRCP<const GlobalOrdinal> gindrowview;
      RowInfo sizeInfo = graph_->getFullGlobalView(lrow, gindrowview);
      ArrayRCP<Scalar> valrowview = getFullViewNonConst(lrow, sizeInfo);
      while (ind != indices.end()) {
        size_t loc = graph_->findGlobalIndex(lrow,*ind,gindrowview);
        if (loc != STINV) {
          valrowview[loc] += (*val);
        }
        ++ind;
        ++val;
      }
      valrowview = Teuchos::null;
      gindrowview = Teuchos::null;
    }
    //else {
      // indices are not allocated; nothing to do
    //}
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getLocalRowCopy(
                                LocalOrdinal LocalRow, 
                                const Teuchos::ArrayView<LocalOrdinal> &indices, 
                                const Teuchos::ArrayView<Scalar>       &values,
                                size_t &numEntries) const {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(isGloballyIndexed()==true && hasColMap()==false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowCopy(): local indices cannot be produced.");
    TEST_FOR_EXCEPTION(getRowMap()->isNodeLocalElement(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowCopy(LocalRow,...): specified row (==" << LocalRow << ") is not valid on this node.");
    if (graph_->isLocallyIndexed()) {
      ArrayRCP<const LocalOrdinal> indrowview;
      ArrayRCP<const Scalar>       valrowview; 
      RowInfo sizeInfo = graph_->getFullLocalView(LocalRow, indrowview);
      valrowview = getFullView(LocalRow, sizeInfo);
      numEntries = sizeInfo.numEntries;
      TEST_FOR_EXCEPTION(indices.size() < numEntries || values.size() < numEntries, std::runtime_error, 
          Teuchos::typeName(*this) << "::getLocalRowCopy(LocalRow,indices,values): size of indices,values must be sufficient to store the specified row.");
      if (numEntries > 0) {
        std::copy( indrowview.begin(), indrowview.begin() + numEntries, indices.begin() );
        std::copy( valrowview.begin(), valrowview.begin() + numEntries, values.begin() );
      }
      valrowview = Teuchos::null;
      indrowview = Teuchos::null;
    }
    else if (graph_->isGloballyIndexed()) {
      ArrayRCP<const GlobalOrdinal> indrowview;
      ArrayRCP<const Scalar>        valrowview; 
      RowInfo sizeInfo = graph_->getFullGlobalView(LocalRow, indrowview);
      valrowview = getFullView(LocalRow, sizeInfo);
      numEntries = sizeInfo.numEntries;
      TEST_FOR_EXCEPTION(indices.size() < numEntries || values.size() < numEntries, std::runtime_error, 
          Teuchos::typeName(*this) << "::getLocalRowCopy(LocalRow,indices,values): size of indices,values must be sufficient to store the specified row.");
      if (numEntries > 0) {
        std::copy( valrowview.begin(), valrowview.begin() + numEntries, values.begin() );
      }
      for (size_t j=0; j < numEntries; ++j) {
        indices[j] = getColMap()->getLocalElement(indrowview[j]);
      }
      valrowview = Teuchos::null;
      indrowview = Teuchos::null;
    }
    else {
#ifdef HAVE_TPETRA_DEBUG
      // should have fallen in one of the above
      TEST_FOR_EXCEPTION( valuesAreAllocated_ == true, std::logic_error, 
          Teuchos::typeName(*this) << "::getLocalRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
      numEntries = 0;
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalRowCopy(
                                GlobalOrdinal globalRow, 
                                const Teuchos::ArrayView<GlobalOrdinal> &indices,
                                const Teuchos::ArrayView<Scalar>        &values,
                                size_t &numEntries) const {
    using Teuchos::ArrayRCP;
    // Only locally owned rows can be queried, otherwise complain
    size_t myRow = getRowMap()->getLocalElement(globalRow);
    TEST_FOR_EXCEPTION(myRow == LOT::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowCopy(globalRow,...): globalRow does not belong to this node.");
    if (graph_->isGloballyIndexed()) {
      ArrayRCP<const GlobalOrdinal> indrowview;
      ArrayRCP<const Scalar>        valrowview; 
      RowInfo sizeInfo = graph_->getFullGlobalView(myRow, indrowview);
      valrowview = getFullView(myRow, sizeInfo);
      numEntries = sizeInfo.numEntries;
      TEST_FOR_EXCEPTION(indices.size() < numEntries || values.size() < numEntries, std::runtime_error, 
          Teuchos::typeName(*this) << "::getGlobalRowCopy(GlobalRow,indices,values): size of indices,values must be sufficient to store the specified row.");
      if (numEntries > 0) {
        std::copy( indrowview.begin(), indrowview.begin() + numEntries, indices.begin() );
        std::copy( valrowview.begin(), valrowview.begin() + numEntries, values.begin() );
      }
      valrowview = Teuchos::null;
      indrowview = Teuchos::null;
    }
    else if (graph_->isLocallyIndexed()) {
      ArrayRCP<const LocalOrdinal> indrowview;
      ArrayRCP<const Scalar>       valrowview; 
      RowInfo sizeInfo = graph_->getFullLocalView(myRow, indrowview);
      valrowview = getFullView(myRow, sizeInfo);
      numEntries = sizeInfo.numEntries;
      TEST_FOR_EXCEPTION(indices.size() < numEntries || values.size() < numEntries, std::runtime_error, 
          Teuchos::typeName(*this) << "::getGlobalRowCopy(GlobalRow,indices,values): size of indices,values must be sufficient to store the specified row.");
      if (numEntries > 0) {
        std::copy( valrowview.begin(), valrowview.begin() + numEntries, values.begin() );
      }
      for (size_t j=0; j < numEntries; ++j) {
        indices[j] = getColMap()->getGlobalElement(indrowview[j]);
      }
      valrowview = Teuchos::null;
      indrowview = Teuchos::null;
    }
    else {
#ifdef HAVE_TPETRA_DEBUG
      // should have fallen in one of the above
      TEST_FOR_EXCEPTION( valuesAreAllocated_ == true, std::logic_error, 
          Teuchos::typeName(*this) << "::getGlobalRowCopy(): Internal logic error. Please contact Tpetra team.");
#endif
      numEntries = 0;
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalRowView(
                                GlobalOrdinal GlobalRow, 
                                Teuchos::ArrayRCP<const GlobalOrdinal> &indices,
                                Teuchos::ArrayRCP<const Scalar>        &values) const {
    TEST_FOR_EXCEPTION(isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowView(): global indices do not exist; call getLocalRowView().");
    LocalOrdinal lrow = getRowMap()->getLocalElement(GlobalRow);
    TEST_FOR_EXCEPTION(lrow == LOT::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowView(GlobalRow,...): GlobalRow (== " << GlobalRow << ") does not belong to this node.");
    RowInfo sizeInfo = graph_->getFullGlobalView(lrow,indices);
    values = getFullView(lrow,sizeInfo);
    if (sizeInfo.numEntries != sizeInfo.allocSize) {
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION( indices == Teuchos::null && values != Teuchos::null, std::logic_error, 
          Teuchos::typeName(*this) << "::getGlobalRowView(): Internal logic error. Please contact Tpetra team.");
#endif
      if (indices != Teuchos::null) {
        indices = indices.persistingView(0,sizeInfo.numEntries);
        if (values != Teuchos::null) {
          values  =  values.persistingView(0,sizeInfo.numEntries);
        }
      }
    }
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getLocalRowView(
                                LocalOrdinal LocalRow, 
                                Teuchos::ArrayRCP<const LocalOrdinal> &indices,
                                Teuchos::ArrayRCP<const Scalar>        &values) const {
    TEST_FOR_EXCEPTION(isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowView(): local indices do not exist; call getGlobalRowView().");
    TEST_FOR_EXCEPTION(getRowMap()->isNodeLocalElement(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowView(LocalRow,...): LocalRow (== " << LocalRow << ") is not valid on this node.");
    RowInfo sizeInfo = graph_->getFullLocalView(LocalRow,indices);
    values = getFullView(LocalRow,sizeInfo);
    if (sizeInfo.numEntries != sizeInfo.allocSize) {
#ifdef HAVE_TPETRA_DEBUG
      TEST_FOR_EXCEPTION( indices == Teuchos::null && values != Teuchos::null, std::logic_error, 
          Teuchos::typeName(*this) << "::getLocalRowView(): Internal logic error. Please contact Tpetra team.");
#endif
      if (indices != Teuchos::null) {
        indices = indices.persistingView(0,sizeInfo.numEntries);
        if (values != Teuchos::null) {
          values  =  values.persistingView(0,sizeInfo.numEntries);
        }
      }
    }
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::scale(const Scalar &alpha) {
    // scale all values in the matrix
    // it is easiest to scale all allocated values, instead of scaling only the ones with valid entries
    // however, if there are no valid entries, we can short-circuit
    // furthermore, if the values aren't allocated, we can short-circuit (unallocated values are zero, scaling to zero)
    const size_t     nlrs = graph_->getNodeNumRows(),
                 numAlloc = graph_->getNodeAllocationSize(),
               numEntries = graph_->getNodeNumEntries();
    if (valuesAreAllocated_ == false || numAlloc == 0 || numEntries == 0) {
      // do nothing
    }
    else {
      Teuchos::RCP<Node> node = lclMatrix_.getNode();
      Kokkos::ReadyBufferHelper<Node> rbh(node);
      ScaleOp<Scalar> wdp;
      wdp.val = alpha;
      if (graph_->getProfileType() == StaticProfile) {
        rbh.begin();
        wdp.x = rbh.addNonConstBuffer(pbuf_values1D_);
        rbh.end();
        node->template parallel_for<ScaleOp<Scalar> >(0,numAlloc,wdp);
      }
      else if (graph_->getProfileType() == DynamicProfile) {
        for (size_t row=0; row < nlrs; ++row) {
          if (pbuf_values2D_[row] != Teuchos::null) {
            rbh.begin();
            wdp.x = rbh.addNonConstBuffer(pbuf_values2D_[row]);
            rbh.end();
            node->template parallel_for<ScaleOp<Scalar> >(0,pbuf_values2D_[row].size(),wdp);
          }
        }
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::setAllToScalar(const Scalar &alpha) {
    // set all values in the matrix
    // it is easiest to set all allocated values, instead of setting only the ones with valid entries
    // however, if there are no valid entries, we can short-circuit
    // we must allocate, though we can set the values in the allocateValues() call.
    // this method is equivalent replacing all valid entries with alpha.
    const size_t     nlrs = graph_->getNodeNumRows(),
                 numAlloc = graph_->getNodeAllocationSize(),
               numEntries = graph_->getNodeNumEntries();
    if (numEntries == 0) {
      // do nothing
    }
    else if (valuesAreAllocated_) {
      // numEntries > 0, so graph must already be allocated; allocatedValues() will verify this
      allocateValues(alpha);
    }
    else {
      Teuchos::RCP<Node> node = lclMatrix_.getNode();
      Kokkos::ReadyBufferHelper<Node> rbh(node);
      InitOp<Scalar> wdp;
      wdp.val = alpha;
      if (graph_->getProfileType() == StaticProfile) {
        rbh.begin();
        wdp.x = rbh.addNonConstBuffer(pbuf_values1D_);
        rbh.end();
        node->template parallel_for<InitOp<Scalar> >(0,numAlloc,wdp);
      }
      else if (graph_->getProfileType() == DynamicProfile) {
        for (size_t row=0; row < nlrs; ++row) {
          if (pbuf_values2D_[row] != Teuchos::null) {
            rbh.begin();
            wdp.x = rbh.addNonConstBuffer(pbuf_values2D_[row]);
            rbh.end();
            node->template parallel_for<InitOp<Scalar> >(0,pbuf_values2D_[row].size(),wdp);
          }
        }
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal> &dvec) const {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << ": cannot call getLocalDiagCopy() until fillComplete() has been called.");
    TEST_FOR_EXCEPTION(dvec.getMap()->isSameAs(*getRowMap()) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalDiagCopy(dvec): dvec must have the same map as the CrsMatrix.");
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid();
#ifdef HAVE_TPETRA_DEBUG
    int numDiagFound = 0;
#endif
    const size_t nlrs = getNodeNumRows();
    Teuchos::ArrayRCP<Scalar> vecView = dvec.get1dViewNonConst();
    for (size_t r=0; r < nlrs; ++r) {
      vecView[r] = Teuchos::ScalarTraits<Scalar>::zero();
      GlobalOrdinal rgid = getRowMap()->getGlobalElement(r);
      if (getColMap()->isNodeGlobalElement(rgid)) {
        LocalOrdinal rlid = getColMap()->getLocalElement(rgid);
        Teuchos::ArrayRCP<const LocalOrdinal> inds;
        RowInfo sizeInfo = graph_->getFullLocalView(r,inds);
        if (sizeInfo.numEntries > 0) {
          const size_t j = graph_->findLocalIndex(r, rlid, inds);
          if (j != STINV) {
            Teuchos::ArrayRCP<const Scalar> vals;
            vals = getFullView(r,sizeInfo);
            vecView[r] = vals[j];
            vals = Teuchos::null;
#ifdef HAVE_TPETRA_DEBUG
            ++numDiagFound;
#endif
          }
        }
        inds = Teuchos::null;
      }
    }
    vecView = Teuchos::null;
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(numDiagFound != getNodeNumDiags(), std::logic_error, 
        "CrsMatrix::getLocalDiagCopy(): logic error. Please contact Tpetra team.");
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::globalAssemble() {
    using Teuchos::arcp;
    using Teuchos::OrdinalTraits;
    using Teuchos::Array;
    using Teuchos::SerialDenseMatrix;
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::rcp;
    using std::list;
    using std::pair;
    using std::make_pair;
    using Teuchos::tuple;
    typedef typename std::map<GlobalOrdinal,std::list<pair<GlobalOrdinal,Scalar> > >::const_iterator NLITER;
    int numImages = getComm()->getSize();
    int myImageID = getComm()->getRank();
    // Determine if any nodes have global entries to share
    size_t MyNonlocals = nonlocals_.size(), 
           MaxGlobalNonlocals;
    Teuchos::reduceAll<int,size_t>(*getComm(),Teuchos::REDUCE_MAX,MyNonlocals,&MaxGlobalNonlocals);
    if (MaxGlobalNonlocals == 0) return;  // no entries to share

    // compute a list of NLRs from nonlocals_ and use it to compute:
    //      IdsAndRows: a vector of (id,row) pairs
    //          NLR2Id: a map from NLR to the Id that owns it
    // globalNeighbors: a global graph of connectivity between images: globalNeighbors(i,j) indicates that j sends to i
    //         sendIDs: a list of all images I send to
    //         recvIDs: a list of all images I receive from (constructed later)
    Array<pair<int,GlobalOrdinal> > IdsAndRows;
    std::map<GlobalOrdinal,int> NLR2Id;
    SerialDenseMatrix<int,char> globalNeighbors;
    Array<int> sendIDs, recvIDs;
    {
      // nonlocals_ contains the entries we are holding for all non-local rows
      // we want a list of the rows for which we have data
      Array<GlobalOrdinal> NLRs;
      std::set<GlobalOrdinal> setOfRows;
      for (NLITER iter = nonlocals_.begin(); iter != nonlocals_.end(); ++iter)
      {
        setOfRows.insert(iter->first);
      }
      // copy the elements in the set into an Array
      NLRs.resize(setOfRows.size());
      std::copy(setOfRows.begin(), setOfRows.end(), NLRs.begin());

      // get a list of ImageIDs for the non-local rows (NLRs)
      Array<int> NLRIds(NLRs.size());
      {
        LookupStatus stat = getRowMap()->getRemoteIndexList(NLRs(),NLRIds());
        char lclerror = ( stat == IDNotPresent ? 1 : 0 );
        char gblerror;
        Teuchos::reduceAll(*getComm(),Teuchos::REDUCE_MAX,lclerror,&gblerror);
        TEST_FOR_EXCEPTION(gblerror, std::runtime_error,
            Teuchos::typeName(*this) << "::globalAssemble(): non-local entries correspond to invalid rows.");
      }

      // build up a list of neighbors, as well as a map between NLRs and Ids
      // localNeighbors[i] != 0 iff I have data to send to image i
      // put NLRs,Ids into an array of pairs
      IdsAndRows.reserve(NLRs.size());
      Array<char> localNeighbors(numImages,0);
      typename Array<GlobalOrdinal>::const_iterator nlr;
      typename Array<int>::const_iterator id;
      for (nlr = NLRs.begin(), id = NLRIds.begin();
           nlr != NLRs.end(); ++nlr, ++id) 
      {
        NLR2Id[*nlr] = *id;
        localNeighbors[*id] = 1;
        IdsAndRows.push_back(make_pair<int,GlobalOrdinal>(*id,*nlr));
      }
      for (int j=0; j<numImages; ++j)
      {
        if (localNeighbors[j]) {
          sendIDs.push_back(j);
        }
      }
      // sort IdsAndRows, by Ids first, then rows
      std::sort(IdsAndRows.begin(),IdsAndRows.end());
      // gather from other nodes to form the full graph
      globalNeighbors.shapeUninitialized(numImages,numImages);
      Teuchos::gatherAll(*getComm(),numImages,localNeighbors.getRawPtr(),numImages*numImages,globalNeighbors.values());
      // globalNeighbors at this point contains (on all images) the
      // connectivity between the images. 
      // globalNeighbors(i,j) != 0 means that j sends to i/that i receives from j
    }

    ////////////////////////////////////////////////////////////////////////////////////// 
    // FIGURE OUT WHO IS SENDING TO WHOM AND HOW MUCH
    // DO THIS IN THE PROCESS OF PACKING ALL OUTGOING DATA ACCORDING TO DESTINATION ID
    ////////////////////////////////////////////////////////////////////////////////////// 

    // loop over all columns to know from which images I can expect to receive something
    for (int j=0; j<numImages; ++j)
    {
      if (globalNeighbors(myImageID,j)) {
        recvIDs.push_back(j);
      }
    }
    size_t numRecvs = recvIDs.size();

    // we know how many we're sending to already
    // form a contiguous list of all data to be sent
    // track the number of entries for each ID
    Array<CrsIJV<GlobalOrdinal,Scalar> > IJVSendBuffer;
    Array<size_t> sendSizes(sendIDs.size(), 0);
    size_t numSends = 0;
    for (typename Array<pair<int,GlobalOrdinal> >::const_iterator IdAndRow = IdsAndRows.begin();
         IdAndRow != IdsAndRows.end(); ++IdAndRow) 
    {
      int            id = IdAndRow->first;
      GlobalOrdinal row = IdAndRow->second;
      // have we advanced to a new send?
      if (sendIDs[numSends] != id) {
        numSends++;
        TEST_FOR_EXCEPTION(sendIDs[numSends] != id, std::logic_error, Teuchos::typeName(*this) << "::globalAssemble(): internal logic error. Contact Tpetra team.");
      }
      // copy data for row into contiguous storage
      for (typename list<pair<GlobalOrdinal,Scalar> >::const_iterator jv = nonlocals_[row].begin(); jv != nonlocals_[row].end(); ++jv)
      {
        IJVSendBuffer.push_back( CrsIJV<GlobalOrdinal,Scalar>(row,jv->first,jv->second) );
        sendSizes[numSends]++;
      }
    }
    if (IdsAndRows.size() > 0) {
      numSends++; // one last increment, to make it a count instead of an index
    }
    TEST_FOR_EXCEPTION(Teuchos::as<typename Array<int>::size_type>(numSends) != sendIDs.size(), std::logic_error, Teuchos::typeName(*this) << "::globalAssemble(): internal logic error. Contact Tpetra team.");

    // don't need this data anymore
    nonlocals_.clear();

    ////////////////////////////////////////////////////////////////////////////////////// 
    // TRANSMIT SIZE INFO BETWEEN SENDERS AND RECEIVERS
    ////////////////////////////////////////////////////////////////////////////////////// 
    // perform non-blocking sends: send sizes to our recipients
    Array<Teuchos::RCP<Teuchos::CommRequest> > sendRequests;
    for (int s=0; s < numSends ; ++s) {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      sendRequests.push_back( Teuchos::isend<int,size_t>(*getComm(),Teuchos::rcpFromRef(sendSizes[s]),sendIDs[s]) );
    }
    // perform non-blocking receives: receive sizes from our senders
    Array<Teuchos::RCP<Teuchos::CommRequest> > recvRequests;
    Array<size_t> recvSizes(numRecvs);
    for (int r=0; r < numRecvs; ++r) {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      recvRequests.push_back( Teuchos::ireceive(*getComm(),rcp(&recvSizes[r],false),recvIDs[r]) );
    }
    // wait on all 
    if (!sendRequests.empty()) {
      Teuchos::waitAll(*getComm(),sendRequests());
    }
    if (!recvRequests.empty()) {
      Teuchos::waitAll(*getComm(),recvRequests());
    }
    Teuchos::barrier(*getComm());
    sendRequests.clear();
    recvRequests.clear();

    ////////////////////////////////////////////////////////////////////////////////////
    // NOW SEND/RECEIVE ALL ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // from the size info, build the ArrayViews into IJVSendBuffer
    Array<ArrayView<CrsIJV<GlobalOrdinal,Scalar> > > sendBuffers(numSends,Teuchos::null);
    {
      size_t cur = 0;
      for (size_t s=0; s<numSends; ++s) {
        sendBuffers[s] = IJVSendBuffer(cur,sendSizes[s]);
        cur += sendSizes[s];
      }
    }
    // perform non-blocking sends
    for (size_t s=0; s < numSends ; ++s)
    {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      ArrayRCP<CrsIJV<GlobalOrdinal,Scalar> > tmparcp = arcp(sendBuffers[s].getRawPtr(),0,sendBuffers[s].size(),false);
      sendRequests.push_back( Teuchos::isend<int,CrsIJV<GlobalOrdinal,Scalar> >(*getComm(),tmparcp,sendIDs[s]) );
    }
    // calculate amount of storage needed for receives
    // setup pointers for the receives as well
    size_t totalRecvSize = std::accumulate(recvSizes.begin(),recvSizes.end(),0);
    Array<CrsIJV<GlobalOrdinal,Scalar> > IJVRecvBuffer(totalRecvSize);
    // from the size info, build the ArrayViews into IJVRecvBuffer
    Array<ArrayView<CrsIJV<GlobalOrdinal,Scalar> > > recvBuffers(numRecvs,Teuchos::null);
    {
      size_t cur = 0;
      for (size_t r=0; r<numRecvs; ++r) {
        recvBuffers[r] = IJVRecvBuffer(cur,recvSizes[r]);
        cur += recvSizes[r];
      }
    }
    // perform non-blocking recvs
    for (size_t r=0; r < numRecvs ; ++r)
    {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      ArrayRCP<CrsIJV<GlobalOrdinal,Scalar> > tmparcp = arcp(recvBuffers[r].getRawPtr(),0,recvBuffers[r].size(),false);
      recvRequests.push_back( Teuchos::ireceive(*getComm(),tmparcp,recvIDs[r]) );
    }
    // perform waits
    if (!sendRequests.empty()) {
      Teuchos::waitAll(*getComm(),sendRequests());
    }
    if (!recvRequests.empty()) {
      Teuchos::waitAll(*getComm(),recvRequests());
    }
    Teuchos::barrier(*getComm());
    sendRequests.clear();
    recvRequests.clear();


    ////////////////////////////////////////////////////////////////////////////////////
    // NOW PROCESS THE RECEIVED ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // TODO: instead of adding one entry at a time, add one row at a time.
    //       this requires resorting; they arrived sorted by sending node, so that entries could be non-contiguous if we received
    //       multiple entries for a particular row from different processors.
    //       it also requires restoring the data, which may make it not worth the trouble.
    for (typename Array<CrsIJV<GlobalOrdinal,Scalar> >::const_iterator ijv = IJVRecvBuffer.begin(); ijv != IJVRecvBuffer.end(); ++ijv)
    {
      try {
        insertGlobalValues(ijv->i, tuple(ijv->j), tuple(ijv->v));
      }
      catch (std::runtime_error &e) {
        std::ostringstream omsg;
        omsg << e.what() << std::endl
          << "caught in globalAssemble() in " << __FILE__ << ":" << __LINE__ << std::endl ;
        throw std::runtime_error(omsg.str());
      }
    }

    // WHEW! THAT WAS TIRING!
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::fillComplete(OptimizeOption os) {
    fillComplete(getRowMap(),getRowMap(),os);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::fillComplete(
                                            const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, 
                                            const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, 
                                            OptimizeOption os) {
    if (getComm()->getSize() > 1) {
      globalAssemble();
    }
    else {
      TEST_FOR_EXCEPTION(nonlocals_.size() > 0, std::runtime_error,
          Teuchos::typeName(*this) << "::fillComplete(): cannot have non-local entries on a serial run. Invalid entry was submitted to the CrsMatrix.");
    }

    if (graph_->isFillComplete()) {
      // if the graph is filled, it should be because 
      // - it was given at construction, already filled
      // - it was filled on a previous call to fillComplete()
      // if neither, it means that the graph has been filled by the user.
      // this means that a graph passed at construction has been filled in the interim.
      // in this case, indices have been sorted and merged, so that we are no longer able
      // to sort or merge the associated values. nothing we can do here.
      if ((constructedWithFilledGraph_ || isFillComplete()) == false) {
        TPETRA_ABUSE_WARNING(true,std::runtime_error,
            "::fillComplete(): fillComplete() has been called on graph since matrix was constructed.");
        return;
      }
    }

    if (isStaticGraph() == false) {
      graph_->makeIndicesLocal(domainMap,rangeMap);
    }
    sortEntries();
    mergeRedundantEntries();
    if (isStaticGraph() == false) {
      // can't optimize graph storage before optimizing our storage
      graph_->fillComplete(domainMap,rangeMap,DoNotOptimizeStorage);
    }

    fillComplete_ = true;

    if (os == DoOptimizeStorage) {
      // this will call optimizeStorage() on the graph as well, if isStaticGraph() == false
      // these routines will fill the local graph and local matrix
      optimizeStorage();
    }
    else { 
      // local graph already filled.
      // fill the local matrix.
      fillLocalMatrix();
    }

    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::sortEntries() {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(graph_->isGloballyIndexed() == true, std::logic_error,
        Teuchos::typeName(*this) << "::sortEntries(): sortEntries() must be called after indices are transformed to local.\n"
        << "Likely internal logic error. Please contact Tpetra team.");
    if (graph_->isSorted()) return;
    if (graph_->getNodeAllocationSize()) {
      const size_t nlrs = getNodeNumRows();
      for (size_t r=0; r < nlrs; ++r) {
        // TODO: This is slightly inefficient, because it may query pbuf_rowOffsets_ repeatadly. 
        //       However, it is very simple code. Consider rewriting it.
        Teuchos::ArrayRCP<LocalOrdinal> inds;
        Teuchos::ArrayRCP<Scalar>       vals;
        RowInfo sizeInfo = graph_->getFullLocalViewNonConst(r, inds);
        vals = getFullViewNonConst(r, sizeInfo);
        sort2(inds.begin(), inds.begin() + sizeInfo.numEntries, vals);
      }
    }
    graph_->setSorted(true);  // we just sorted them
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::mergeRedundantEntries() {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(graph_->isSorted() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::mergeRedundantEntries() cannot be called before indices are sorted.\n"
        << "Likely interal logic error. Please contact Tpetra team.");
    if (graph_->notRedundant()) return;
    if (graph_->getNodeAllocationSize() > 0) {
      for (size_t r=0; r<getNodeNumRows(); ++r) {
        // TODO: This is slightly inefficient, because it may query pbuf_rowOffsets_ repeatadly. 
        //       However, it is very simple code. Consider rewriting it.
        Teuchos::ArrayRCP<LocalOrdinal> inds;
        Teuchos::ArrayRCP<Scalar>       vals;
        RowInfo sizeInfo = graph_->getFullLocalViewNonConst(r, inds);
        vals = getFullViewNonConst(r, sizeInfo);
        if (sizeInfo.numEntries > 0) {
          size_t curEntry = 0;
          Scalar curValue = vals[curEntry];
          for (size_t k=1; k < sizeInfo.numEntries; ++k) {
            if (inds[k] == inds[k-1]) {
              curValue += vals[k];
            }
            else {
              vals[curEntry++] = curValue;
              curValue = vals[k];
            }
          }
          vals[curEntry] = curValue;
        }
        vals = Teuchos::null;
        inds = Teuchos::null;
      }
      graph_->removeRedundantIndices();
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::optimizeStorage() {
    // optimizeStorage will perform two functions:
    // 1) create a single allocation of memory
    // 2) pack data in that allocation
    // if getProfileType() == StaticProfile, then 1) has already been done
#ifdef HAVE_TPETRA_DEBUG
    std::string err = Teuchos::typeName(*this) + "::optimizeStorage(): Internal logic error. Please contact Tpetra team.";
#endif
    if (isStorageOptimized() == true) return;
    TEST_FOR_EXCEPTION(isFillComplete() == false || graph_->isSorted() == false || graph_->notRedundant() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::optimizeStorage(): fillComplete() must be called before optimizeStorage().");
    TEST_FOR_EXCEPTION(graph_->isStorageOptimized(), std::logic_error,
        Teuchos::typeName(*this) << "::optimizeStorage(): optimizeStorage() already called on graph, but not on matrix.");
    // 
    Teuchos::RCP<Node> node = lclMatrix_.getNode();
    const size_t nlrs = getNodeNumRows(),
           numEntries = graph_->getNodeNumEntries();
    if (nlrs > 0) {
      if (valuesAreAllocated_ == false) {
        // haven't even allocated values yet. allocate to optimal size and fill with zero
        if (numEntries > 0) {
          // allocate
          pbuf_values1D_ = node->template allocBuffer<Scalar>(numEntries);        
          // init to zero
          Kokkos::ReadyBufferHelper<Node> rbh(node);
          InitOp<Scalar> wdp;
          wdp.val = Teuchos::ScalarTraits<Scalar>::zero();
          rbh.begin();
          wdp.x   = rbh.addNonConstBuffer(pbuf_values1D_);
          rbh.end();
          node->template parallel_for<InitOp<Scalar> >(0,numEntries,wdp);
        } 
        valuesAreAllocated_ = true;
      }
      else if (graph_->getProfileType() == DynamicProfile) {
        // allocate a single memory block
        if (numEntries > 0) {
          // numEntries > 0  =>  allocSize > 0  =>  pbuf_values2D_ != Teuchos::null
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION( pbuf_values2D_ == Teuchos::null, std::logic_error, err);
#endif
          pbuf_values1D_ = node->template allocBuffer<Scalar>(numEntries);
          size_t sofar = 0;
          for (size_t row=0; row<nlrs; ++row) {
            RowInfo sizeInfo = graph_->getRowInfo(row);
            if (sizeInfo.numEntries > 0) {
              node->template copyBuffers<Scalar>( sizeInfo.numEntries, pbuf_values2D_[row], pbuf_values1D_ + sofar );
              pbuf_values2D_[row] = Teuchos::null;
            }
            sofar += sizeInfo.numEntries;
          }
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION(sofar != numEntries, std::logic_error, err);
#endif
        }
        pbuf_values2D_ = Teuchos::null;
      }
      else {
        // storage is already allocated; just need to pack
        if (numEntries > 0) {
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION( pbuf_values1D_ == Teuchos::null, std::logic_error, err);
#endif
          size_t sofar = 0;
          for (size_t row=0; row<nlrs; ++row) {
            RowInfo sizeInfo = graph_->getRowInfo(row);
            if (sizeInfo.numEntries > 0 && sizeInfo.offset1D != sofar) {
              node->template copyBuffers<Scalar>( sizeInfo.numEntries, 
                                                  pbuf_values1D_ + sizeInfo.offset1D,
                                                  pbuf_values1D_ + sofar );
            }
            sofar += sizeInfo.numEntries;
          }
          pbuf_values1D_ = pbuf_values1D_.persistingView(0,sofar);
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION(sofar != numEntries, std::logic_error, err);
#endif
        }
      }
    }
    storageOptimized_ = true;

    if (isStaticGraph() == false) {
      graph_->optimizeStorage();
    }

    // local graph was filled during graph_->optimizeStorage()
    fillLocalMatrix(); 
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::applyInverse(
                                    const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                                          MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
                                          Teuchos::ETransp mode) const {
    solve<Scalar,Scalar>(Y,X,mode);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  template <class DomainScalar, class RangeScalar>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::solve(
                                    const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node>  &Y, 
                                          MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                                          Teuchos::ETransp trans) const {
//     // Solve U X = Y  or  L X = Y
//     // X belongs to domain map, while Y belongs to range map
//     typedef Teuchos::ScalarTraits<Scalar> ST;
//     using Teuchos::null;
//     using Teuchos::ArrayView;
//     TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error, 
//         Teuchos::typeName(*this) << ": cannot call applyInverse() until fillComplete() has been called.");
//     TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
//         Teuchos::typeName(*this) << "::applyInverse(X,Y): X and Y must have the same number of vectors.");
//     TEST_FOR_EXCEPTION(isLowerTriangular() == false && isUpperTriangular() == false, std::runtime_error,
//         Teuchos::typeName(*this) << "::applyInverser() requires either upper or lower triangular structure.");
// 
// #ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
//     int myImageID = Teuchos::rank(*getComm());
//     Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
//     if (myImageID == 0) {
//       *out << "Entering CrsMatrix::applyInverse()" << std::endl
//                 << "Column Map: " << std::endl;
//     }
//     *out << *this->getColMap() << std::endl;
//     if (myImageID == 0) {
//       *out << "Initial input: " << std::endl;
//     }
//     Y.describe(*out,Teuchos::VERB_EXTREME);
// #endif
// 
//     const size_t numVectors = X.getNumVectors();
//     Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal> > importer = graph_->getImporter();
//     Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal> > exporter = graph_->getExporter();
//     Teuchos::RCP<MV>        Xcopy;
//     Teuchos::RCP<const MV>  Ycopy;
//     Kokkos::MultiVector<Scalar,Node>        *lclX = &X.getLocalMVNonConst();
//     const Kokkos::MultiVector<Scalar,Node>  *lclY = &Y.getLocalMV();
// 
//     // cannot handle non-constant stride right now
//     if (X.isConstantStride() == false) {
//       // generate a copy of X 
//       Xcopy = Teuchos::rcp(new MV(X));
//       lclX = &Xcopy->getLocalMVNonConst();
// #ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
//       if (myImageID == 0) *out << "X is not constant stride, duplicating X results in a strided copy" << std::endl;
//       Xcopy->describe(*out,Teuchos::VERB_EXTREME);
// #endif
//     }
// 
//     // cannot handle non-constant stride right now
//     if (Y.isConstantStride() == false) {
//       // generate a copy of Y 
//       Ycopy = Teuchos::rcp(new MV(Y));
//       lclY = &Ycopy->getLocalMV();
// #ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
//       if (myImageID == 0) *out << "Y is not constant stride, duplicating Y results in a strided copy" << std::endl;
//       Ycopy->describe(*out,Teuchos::VERB_EXTREME);
// #endif
//     }
// 
//     // it is okay if X and Y reference the same data, because we can perform a triangular solve in-situ
//     // however, for simplicity, for now we will not support this use case
//     if (lclX==lclY && importer==null && exporter==null) {
//       TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
//           "::applyInverse(X,Y): If X and Y are the same, it necessitates a temporary copy of Y, which is inefficient.");
//       // generate a copy of Y
//       Ycopy = Teuchos::rcp(new MV(Y));
//       lclY = &Ycopy->getLocalMV();
// #ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
//       if (myImageID == 0) *out << "X and Y are co-located, duplicating Y results in a strided copy" << std::endl;
//       Ycopy->describe(*out,Teuchos::VERB_EXTREME);
// #endif
//     }
// 
//     if (importer != null) {
//       if (importMV_ != null && importMV_->getNumVectors() != numVectors) importMV_ = null;
//       if (importMV_ == null) {
//         importMV_ = Teuchos::rcp( new MV(getColMap(),numVectors) );
//       }
//     }
//     if (exporter != null) {
//       if (exportMV_ != null && exportMV_->getNumVectors() != numVectors) exportMV_ = null;
//       if (exportMV_ == null) {
//         exportMV_ = Teuchos::rcp( new MV(getRowMap(),numVectors) );
//       }
//     }
// 
//     // import/export patterns are different depending on transpose.
//     if (mode == Teuchos::NO_TRANS) {
//       // applyInverse(NO_TRANS): RangeMap -> DomainMap
//       // lclMatSolve_: RowMap -> ColMap
//       // importer: DomainMap -> ColMap
//       // exporter: RowMap -> RangeMap
//       // 
//       // applyInverse = reverse(exporter)  o   lclMatSolve_  o reverse(importer)
//       //                RangeMap   ->    RowMap     ->     ColMap         ->    DomainMap
//       // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
//       if (exporter != null) {
//         exportMV_->doImport(Y, *exporter, INSERT);
// #ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
//         if (myImageID == 0) {
//           *out << "Performed import of Y using exporter..." << std::endl;
//         }
//         exportMV_->describe(*out,Teuchos::VERB_EXTREME);
// #endif
//         lclY = &exportMV_->getLocalMV();
//       }
//       // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
//       // We will compute solution into the to-be-exported MV; get a view
//       if (importer != null) {
//         lclX = &importMV_->getLocalMVNonConst();
//       }
//       // Do actual computation
//       if (isUpperTriangular()) {
//         lclMatSolve_.template solve<Scalar,Scalar>(Teuchos::NO_TRANS, Teuchos::UPPER_TRI, Teuchos::NON_UNIT_DIAG, Teuchos::ScalarTraits<Scalar>::one(), *lclY, Teuchos::ScalarTraits<Scalar>::zero(), *lclX);
//       }
//       else {
//         lclMatSolve_.template solve<Scalar,Scalar>(Teuchos::NO_TRANS, Teuchos::LOWER_TRI, Teuchos::NON_UNIT_DIAG, Teuchos::ScalarTraits<Scalar>::one(), *lclY, Teuchos::ScalarTraits<Scalar>::zero(), *lclX);
//       }
// #ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
//       if (myImageID == 0) *out << "Matrix-MV solve..." << std::endl;
//       if (exporter != null) {
//         exportMV_->describe(*out,Teuchos::VERB_EXTREME);
//       } 
//       else {
//         X.describe(*out,Teuchos::VERB_EXTREME);
//       }
// #endif
//       // do the export
//       if (importer != null) {
//         // Make sure target is zero: necessary because we are adding. may need adjusting for alpha,beta apply()
//         X.putScalar(ST::zero());  
//         X.doExport(*importMV_, *importer, ADD); // Fill Y with Values from export vector
// #ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
//         if (myImageID == 0) *out << "Output vector after export() using exporter..." << std::endl;
//         X.describe(*out,Teuchos::VERB_EXTREME);
// #endif
//       }
//       // Handle case of rangemap being a local replicated map: in this case, sum contributions from each processor
//       if (X.isDistributed() == false) {
//         X.reduce();
// #ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
//         if (myImageID == 0) *out << "Output vector is local; result after reduce()..." << std::endl;
//         X.describe(*out,Teuchos::VERB_EXTREME);
// #endif
//       }
//     }
//     else {
//       // mode == CONJ_TRANS or TRANS
//       TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<Scalar>::isComplex && mode == Teuchos::TRANS, std::logic_error,
//           Teuchos::typeName(*this) << "::applyInverse() does not currently support transposed multiplications for complex scalar types.");
//       // applyInverse(TRANS): DomainMap -> RangeMap
//       // lclMatSolve_(TRANS): ColMap -> RowMap
//       // importer: DomainMap -> ColMap
//       // exporter: RowMap -> RangeMap
//       // 
//       // applyInverse =        importer o   lclMatSolve_  o  exporter
//       //                Domainmap -> ColMap     ->      RowMap -> RangeMap
//       // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
//       if (importer != null) {
//         importMV_->doImport(Y,*importer,INSERT);
//         lclY = &importMV_->getLocalMV();
// #ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
//         if (myImageID == 0) {
//           *out << "Performed import of Y using importer..." << std::endl;
//         }
//         importMV_->describe(*out,Teuchos::VERB_EXTREME);
// #endif
//       }
//       // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
//       // We will compute solution into the to-be-exported MV; get a view
//       if (exporter != null) {
//         lclX = &exportMV_->getLocalMVNonConst();
//       }
//       // Do actual computation
//       if (isUpperTriangular()) {
//         lclMatSolve_.template solve<Scalar,Scalar>(Teuchos::CONJ_TRANS, Teuchos::UPPER_TRI, Teuchos::NON_UNIT_DIAG, Teuchos::ScalarTraits<Scalar>::one(), *lclY, Teuchos::ScalarTraits<Scalar>::zero(), *lclX);
//       }
//       else {
//         lclMatSolve_.template solve<Scalar,Scalar>(Teuchos::CONJ_TRANS, Teuchos::LOWER_TRI, Teuchos::NON_UNIT_DIAG, Teuchos::ScalarTraits<Scalar>::one(), *lclY, Teuchos::ScalarTraits<Scalar>::zero(), *lclX);
//       }
// #ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
//       if (myImageID == 0) *out << "Matrix-MV solve..." << std::endl;
//       if (exporter != null) {
//         exportMV_->describe(*out,Teuchos::VERB_EXTREME);
//       } 
//       else {
//         X.describe(*out,Teuchos::VERB_EXTREME);
//       }
// #endif
//       // do the export
//       if (exporter != null) {
//         X.putScalar(ST::zero()); // Make sure target is zero: necessary because we are adding. may need adjusting for alpha,beta apply()
//         X.doExport(*exportMV_,*exporter,ADD);
// #ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
//         if (myImageID == 0) *out << "Output vector after export() using exporter..." << std::endl;
//         X.describe(*out,Teuchos::VERB_EXTREME);
// #endif
//       }
//       // Handle case of rangemap being a local replicated map: in this case, sum contributions from each processor
//       if (X.isDistributed() == false) {
//         X.reduce();
// #ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
//         if (myImageID == 0) *out << "Output vector is local; result after reduce()..." << std::endl;
//         X.describe(*out,Teuchos::VERB_EXTREME);
// #endif
//       }
//     }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::apply(
                                        const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
                                        MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                                        Teuchos::ETransp mode) const {
    multiply<Scalar,Scalar>(X,Y,mode);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  template <class DomainScalar, class RangeScalar>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::multiply(
                                        const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
                                        MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                                        Teuchos::ETransp mode) const {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::null;
    using Teuchos::ArrayView;
    TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error, 
        Teuchos::typeName(*this) << ": cannot call multiply() until fillComplete() has been called.");
    TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(X,Y): X and Y must have the same number of vectors.");

#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
    int myImageID = Teuchos::rank(*getComm());
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    if (myImageID == 0) {
      *out << "Entering CrsMatrix::multiply()" << std::endl
                << "Column Map: " << std::endl;
    }
    *out << *this->getColMap() << std::endl;
    if (myImageID == 0) {
      *out << "Initial input: " << std::endl;
    }
    X.describe(*out,Teuchos::VERB_EXTREME);
#endif

    const size_t numVectors = X.getNumVectors();
    // because of Views, it is difficult to determine if X and Y point to the same data. 
    // however, if they reference the exact same object, we will do the user the favor of copying X into new storage (with a warning)
    // we ony need to do this if we have trivial importers; otherwise, we don't actually apply the operator from X into Y
    Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal> > importer = graph_->getImporter();
    Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal> > exporter = graph_->getExporter();
    Teuchos::RCP<const MV> Xcopy;
    Teuchos::RCP<MV>       Ycopy;
    const Kokkos::MultiVector<Scalar,Node> *lclX = &X.getLocalMV();
    Kokkos::MultiVector<Scalar,Node>       *lclY = &Y.getLocalMVNonConst();

    // cannot handle non-constant stride right now
    if (X.isConstantStride() == false) {
      // generate a copy of X 
      Xcopy = Teuchos::rcp(new MV(X));
      lclX = &Xcopy->getLocalMV();
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "X is not constant stride, duplicating X results in a strided copy" << std::endl;
      Xcopy->describe(*out,Teuchos::VERB_EXTREME);
#endif
    }

    // cannot handle non-constant stride right now
    if (Y.isConstantStride() == false) {
      // generate a copy of Y 
      Ycopy = Teuchos::rcp(new MV(Y));
      lclY = &Ycopy->getLocalMVNonConst();
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "Y is not constant stride, duplicating Y results in a strided copy" << std::endl;
      Ycopy->describe(*out,Teuchos::VERB_EXTREME);
#endif
    }

    // cannot handle in-situ matvec
    if (lclX==lclY && importer==null && exporter==null) {
      TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
          "::multiply(X,Y): If X and Y are the same, it necessitates a temporary copy of X, which is inefficient.");
      // generate a copy of X 
      Xcopy = Teuchos::rcp(new MV(X));
      lclX = &Xcopy->getLocalMV();
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "X and Y are co-located, duplicating X results in a strided copy" << std::endl;
      Xcopy->describe(*out,Teuchos::VERB_EXTREME);
#endif
    }

    if (importer != null) {
      if (importMV_ != null && importMV_->getNumVectors() != numVectors) importMV_ = null;
      if (importMV_ == null) {
        importMV_ = Teuchos::rcp( new MV(getColMap(),numVectors) );
      }
    }
    if (exporter != null) {
      if (exportMV_ != null && exportMV_->getNumVectors() != numVectors) exportMV_ = null;
      if (exportMV_ == null) {
        exportMV_ = Teuchos::rcp( new MV(getRowMap(),numVectors) );
      }
    }

    // import/export patterns are different depending on transpose.
    if (mode == Teuchos::NO_TRANS) {
      // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
      if (importer != null) {
        importMV_->doImport(X, *importer, INSERT);
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) {
          *out << "Performed import of X using importer..." << std::endl;
        }
        importMV_->describe(*out,Teuchos::VERB_EXTREME);
#endif
        lclX = &importMV_->getLocalMV();
      }
      // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
      // We will compute solution into the to-be-exported MV; get a view
      if (exporter != null) {
        lclY = &exportMV_->getLocalMVNonConst();
      }
      // Do actual computation
      lclMatVec_.template multiply<Scalar,Scalar>(Teuchos::NO_TRANS, Teuchos::ScalarTraits<Scalar>::one(), *lclX, Teuchos::ScalarTraits<Scalar>::zero(), *lclY);
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "Matrix-MV product..." << std::endl;
      if (exportMV_ != null) {
        exportMV_->describe(*out,Teuchos::VERB_EXTREME);
      } 
      else {
        Y.describe(*out,Teuchos::VERB_EXTREME);
      }
#endif
      // do the export
      if (exporter != null) {
        // Make sure target is zero: necessary because we are adding. may need adjusting for alpha,beta apply()
        Y.putScalar(ST::zero());  
        Y.doExport(*exportMV_, *exporter, ADD); // Fill Y with Values from export vector
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) *out << "Output vector after export() using exporter..." << std::endl;
        Y.describe(*out,Teuchos::VERB_EXTREME);
#endif
      }
      // Handle case of rangemap being a local replicated map: in this case, sum contributions from each processor
      if (Y.isDistributed() == false) {
        Y.reduce();
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) *out << "Output vector is local; result after reduce()..." << std::endl;
        Y.describe(*out,Teuchos::VERB_EXTREME);
#endif
      }
    }
    else {
      // mode == CONJ_TRANS or TRANS
      TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<Scalar>::isComplex && mode == Teuchos::TRANS, std::logic_error,
          Teuchos::typeName(*this) << "::multiply() does not currently support transposed multiplications for complex scalar types.");
      // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
      if (exporter != null) {
        exportMV_->doImport(X,*exporter,INSERT);
        lclX = &exportMV_->getLocalMV();
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) {
          *out << "Performed import of X using exporter..." << std::endl;
        }
        exportMV_->describe(*out,Teuchos::VERB_EXTREME);
#endif
      }
      // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
      // We will compute colutioni into the to-be-exported MV; get a view
      if (importer != null) {
        lclY = &importMV_->getLocalMVNonConst();
      }
      // Do actual computation
      lclMatVec_.template multiply<Scalar,Scalar>(Teuchos::CONJ_TRANS, Teuchos::ScalarTraits<Scalar>::one(), *lclX, Teuchos::ScalarTraits<Scalar>::zero(), *lclY);
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "Matrix-MV product..." << std::endl;
      if (importMV_ != null) {
        importMV_->describe(*out,Teuchos::VERB_EXTREME);
      } 
      else {
        Y.describe(*out,Teuchos::VERB_EXTREME);
      }
#endif
      if (importer != null) {
        Y.putScalar(ST::zero()); // Make sure target is zero: necessary because we are adding. may need adjusting for alpha,beta multiply()
        Y.doExport(*importMV_,*importer,ADD);
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) *out << "Output vector after export() using importer..." << std::endl;
        Y.describe(*out,Teuchos::VERB_EXTREME);
#endif
      }
      // Handle case of rangemap being a local replicated map: in this case, sum contributions from each processor
      if (Y.isDistributed() == false) {
        Y.reduce();
#ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) *out << "Output vector is local; result after reduce()..." << std::endl;
        Y.describe(*out,Teuchos::VERB_EXTREME);
#endif
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  std::string CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::description() const {
    std::ostringstream oss;
    oss << Teuchos::Describable::description();
    if (isFillComplete()) {
      oss << "{status = fill complete"
          << ", global rows = " << getGlobalNumRows()
          << ", global cols = " << getGlobalNumCols()
          << ", global num entries = " << getGlobalNumEntries()
          << "}";
    }
    else {
      oss << "{status = fill not complete"
          << ", global rows = " << getGlobalNumRows()
          << "}";
    }
    return oss.str();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;
    Teuchos::RCP<const Teuchos::Comm<int> > comm = this->getComm();
    const int myImageID = comm->getRank(),
              numImages = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t>(width,11) + 2;
    Teuchos::OSTab tab(out);
    //    none: print nothing
    //     low: print O(1) info from node 0
    //  medium: print O(P) info, num entries per node
    //    high: print O(N) info, num entries per row
    // extreme: print O(NNZ) info: print indices and values
    // 
    // for medium and higher, print constituent objects at specified verbLevel
    if (vl != VERB_NONE) {
      if (myImageID == 0) out << this->description() << std::endl; 
      // O(1) globals, minus what was already printed by description()
      if (isFillComplete() && myImageID == 0) {
        out << "Global number of diagonals = " << getGlobalNumDiags() << std::endl;
        out << "Global max number of entries = " << getGlobalMaxNumRowEntries() << std::endl;
      }
      // constituent objects
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        if (myImageID == 0) out << "\nRow map: " << std::endl;
        getRowMap()->describe(out,vl);
        if (getColMap() != Teuchos::null) {
          if (myImageID == 0) out << "\nColumn map: " << std::endl;
          getColMap()->describe(out,vl);
        }
        if (getOperatorDomainMap() != Teuchos::null) {
          if (myImageID == 0) out << "\nDomain map: " << std::endl;
          getOperatorDomainMap()->describe(out,vl);
        }
        if (getOperatorRangeMap() != Teuchos::null) {
          if (myImageID == 0) out << "\nRange map: " << std::endl;
          getOperatorRangeMap()->describe(out,vl);
        }
      }
      // O(P) data
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << "Node ID = " << imageCtr << std::endl
                << "Node number of entries = " << getNodeNumEntries() << std::endl
                << "Node number of diagonals = " << getNodeNumDiags() << std::endl
                << "Node max number of entries = " << getNodeMaxNumRowEntries() << std::endl
                << "Node number of allocated entries = " << graph_->getNodeAllocationSize() << std::endl;
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
      // O(N) and O(NNZ) data
      if (vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << std::setw(width) << "Node ID"
                << std::setw(width) << "Global Row" 
                << std::setw(width) << "Num Entries";
            if (vl == VERB_EXTREME) {
              out << std::setw(width) << "(Index,Value)";
            }
            out << std::endl;
            for (size_t r=0; r < getNodeNumRows(); ++r) {
              const size_t nE = getNumEntriesInLocalRow(r);
              GlobalOrdinal gid = getRowMap()->getGlobalElement(r);
              out << std::setw(width) << myImageID 
                  << std::setw(width) << gid
                  << std::setw(width) << nE;
              if (vl == VERB_EXTREME) {
                if (isGloballyIndexed()) {
                  Teuchos::ArrayRCP<const GlobalOrdinal> rowinds;
                  Teuchos::ArrayRCP<const Scalar> rowvals;
                  getGlobalRowView(gid,rowinds,rowvals);
                  for (size_t j=0; j < nE; ++j) {
                    out << "(" << std::setw(width) << rowinds[j]
                        << "," << std::setw(width) << rowvals[j]
                        << ") ";
                  }
                }
                else if (isLocallyIndexed()) {
                  Teuchos::ArrayRCP<const LocalOrdinal> rowinds;
                  Teuchos::ArrayRCP<const Scalar> rowvals;
                  getLocalRowView(r,rowinds,rowvals);
                  for (size_t j=0; j < nE; ++j) {
                    out << "(" << std::setw(width) << getColMap()->getGlobalElement(rowinds[j]) 
                        << "," << std::setw(width) << rowvals[j]
                        << ") ";
                  }
                }
              }
              out << std::endl;
            }
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
    }
  }

} // namespace Tpetra

#endif
