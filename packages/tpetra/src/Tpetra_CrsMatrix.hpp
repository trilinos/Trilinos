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
// TODO: row-wise insertion of entries in globalAssemble()

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_CrsMatrix.hpp>
#include <Kokkos_DefaultSparseMultiply.hpp>
#include <Kokkos_DefaultSparseSolve.hpp>

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
  template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatVec = Kokkos::DefaultSparseMultiply<Scalar,Ordinal,Node>, class LocalMatSolve = Kokkos::DefaultSparseSolve<Scalar,Ordinal,Node> >
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
      CrsMatrix(const Teuchos::RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > &graph);

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
      Teuchos::RCP<Node> getNode() const = 0;

      //! Returns the Map that describes the row distribution in this matrix.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRowMap() const;

      //! \brief Returns the Map that describes the column distribution in this matrix.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getColMap() const;

      //! Returns the RowGraph associated with this matrix. 
      Teuchos::RCP<const RowGraph<LocalOrdinal,GlobalOrdinal> > getGraph() const;

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
                            Teuchos::ArrayRCP<const Scalar>        &values);

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
                           Teuchos::ArrayRCP<const Scalar>       &values);

      //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
      /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the 
          the zero and non-zero diagonals owned by this node. */
      void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal> &diag) const;

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

    protected:
      // copy constructor disabled
      CrsMatrix(const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> &Source);
      // operator= disabled
      CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> & operator=(const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> &rhs);

      enum AllocateLocalGlobal {
        AllocateLocal,
        AllocateGlobal
      };

      // useful typedefs
      typedef Teuchos::OrdinalTraits<LocalOrdinal>    LOT;
      typedef Teuchos::OrdinalTraits<GlobalOrdinal>   GOT;
      typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

      void allocateValues();
      void sortEntries();
      void mergeRedundantEntries();

      // multiplication routines
      inline typename Teuchos::ArrayRCP<const Scalar>::iterator getVptr(size_t row) const;
      inline typename Teuchos::ArrayRCP<Scalar>::iterator getVptr(size_t row);

      Teuchos::RCP<CrsGraph<LocalOrdinal,GlobalOrdinal> > graph_;
      Kokkos::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> lclMatrix_;
      LocalMatVec lclMatVec_;
      LocalMatSolve lclMatSolve_;
      bool staticGraph_;
      bool constructedWithFilledGraph_;
      bool fillComplete_;
      bool storageOptimized_;

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
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype);
  : lclMatrix_(rowMap->getNode()),
  , lclMatVec_(rowMap->getNode()),
  , lclMatSolve_(rowMap->getNode()),
  , staticGraph_(false)
  , constructedWithFilledGraph_(false)
  , fillComplete_(false)
  , storageOptimized_(false) {
    try {
      graph_ = Teuchos::rcp( new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowMap,maxNNZPerRow,staticProfile) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          Teuchos::typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    checkInternalState();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const Teuchos::ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile);
  : lclMatrix_(rowMap->getNode()),
  , lclMatVec_(rowMap->getNode()),
  , lclMatSolve_(rowMap->getNode()),
  , staticGraph_(false)
  , constructedWithFilledGraph_(false)
  , fillComplete_(false)
  , storageOptimized_(false) {
    try {
      graph_ = Teuchos::rcp( new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowMap,maxNNZPerRow,staticProfile) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          Teuchos::typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    checkInternalState();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype)
  : lclMatrix_(rowMap->getNode()),
  , lclMatVec_(rowMap->getNode()),
  , lclMatSolve_(rowMap->getNode()),
  , staticGraph_(false)
  , constructedWithFilledGraph_(false)
  , fillComplete_(false)
  , storageOptimized_(false) {
    try {
      graph_ = Teuchos::rcp( new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowMap,maxNNZPerRow,staticProfile) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          Teuchos::typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    checkInternalState();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrix(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, const Teuchos::ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype)
  : lclMatrix_(rowMap->getNode()),
  , lclMatVec_(rowMap->getNode()),
  , lclMatSolve_(rowMap->getNode()),
  , staticGraph_(false)
  , constructedWithFilledGraph_(false)
  , fillComplete_(false)
  , storageOptimized_(false) {
    try {
      graph_ = Teuchos::rcp( new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(rowMap,maxNNZPerRow,staticProfile) );
    }
    catch (std::exception &e) {
      TEST_FOR_EXCEPTION(true, std::runtime_error,
          Teuchos::typeName(*this) << "::CrsMatrix(): caught exception while allocating CrsGraph object: " 
          << std::endl << e.what() << std::endl);
    }
    checkInternalState();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::CrsMatrix(const Teuchos::RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > &graph)
  : graph_(graph)
  , staticGraph_(true)
  , fillComplete_(false)
  , storageOptimized_(false) {
    TEST_FOR_EXCEPTION(graph_ == Teuchos::null, std::runtime_error,
        Teuchos::typeName(*this) << "::CrsMatrix(graph): specified pointer is null.");
    // we won't prohibit the case where the graph is not yet filled, but we will check below to ensure that the
    // graph isn't filled between now and when fillComplete() is called on this matrix
    constructedWithFilledGraph_ = graph_->isFillComplete();
    checkInternalState();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::~CrsMatrix() {
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::allocateValues() {
    const size_t nlrs = getRowMap()->getNumLocalElements(),
                 nta = graph_->getNodeAllocationSize();
    if (nta > 0) {
      Teuchos::RCP<Node> node = lclMatrix_.getNode();
      if (getProfileType() == StaticProfile) {
        if (nta) {
          pbuf_values1D_ = node->template allocBuffer<Scalar>(nta);
        }
      }
      else { // if getProfileType() == DynamicProfile
        pbuf_values2D_ = Teuchos::arcp< Teuchos::ArrayRCP<Scalar> >(nlrs);
        for (size_t r=0; r<nlrs; ++r) {
          const size_t ntarow = graph_->getNumAllocatedEntriesInLocalRow(r);
          if (ntarow > 0) {
            pbuf_values2D_[r] = node->template allocBuffer<Scalar>(ntarow);
          }
        }
      }
    }
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
  const Teuchos::RCP<const Teuchos::Comm<int> > &
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getComm() const {
    return graph_->getComm();
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
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalNumCols() const { 
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
  global_size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalMaxNumRowEntries() const { 
    return graph_->getGlobalMaxNumRowEntries(); 
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  size_t CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNodeMaxNumRowEntries() const { 
    return graph_->getNodeMaxNumRowEntries(); 
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  GlobalOrdinal CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getIndexBase() const { 
    return getRowMap().getIndexBase(); 
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

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<const RowGraph<LocalOrdinal,GlobalOrdinal> >
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGraph() const { 
    return graph_; 
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  Teuchos::RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal> >
  CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getCrsGraph() const { 
    return graph_; 
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isLowerTriangular() const { 
    return graph_->isLowerTriangular(); 
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isUpperTriangular() const { 
    return graph_->isUpperTriangular(); 
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  bool CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::isStaticGraph() const { 
    return staticGraph_; 
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::insertMyValues(LocalOrdinal localRow, 
                         const Teuchos::ArrayView<const LocalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>       &values) {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(isStorageOptimized() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertMyValues(): cannot insert new values after optimizeStorage() has been called.");
    TEST_FOR_EXCEPTION(graph_->isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertMyValues(): graph indices are global; use insertGlobalValues().");
    TEST_FOR_EXCEPTION(hasColMap() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::insertMyValues(): cannot insert local indices without a column map; ");
    TEST_FOR_EXCEPTION(isStaticGraph() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertMyValues(): matrix was constructed with static graph; cannot insert new entries.");
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::insertMyValues(): values.size() must equal indices.size().");
    TEST_FOR_EXCEPTION(getRowMap().isNodeLocalElement(localRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::insertMyValues(): row does not belong to this node.");
    Teuchos::Array<LocalOrdinal> finds;
    Teuchos::Array<Scalar>       fvals;
    // use column map to filter the entries:
    const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &cmap = getColMap();
    for (size_t i=0; i<indices.size(); ++i) {
      if (cmap.isNodeLocalElement(indices[i])) {
        finds.push_back(indices[i]);
        fvals.push_back(values[i]);
      }
    }
    Teuchos::ArrayView<const LocalOrdinal> findices = finds();
    Teuchos::ArrayView<const Scalar      > fvalues  = fvals();
    size_t rowNNZ = getNumEntriesInLocalRow(localRow),
                     toAdd = findices.size(),
                  rowAlloc = graph_->numAllocatedEntriesForMyRow(localRow);
    if (rowNNZ+toAdd > rowAlloc) {
      TEST_FOR_EXCEPTION(graph_->isStaticProfile() == true, std::runtime_error,
          Teuchos::typeName(*this) << "::insertMyValues(): new indices exceed statically allocated graph structure.");
      TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
          "::insertMyValues(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
#ifdef HAVE_TPETRA_DEBUG
      // assumption: the number currently allocated in values_ is that stored in the graph.
      TEST_FOR_EXCEPTION( rowAlloc != values_[localRow].size(), std::logic_error,
          Teuchos::typeName(*this) << "::insertMyValues(): Internal logic error or unsupported use case. Please contact Tpetra team.");
#endif
      // increase the allocation, copy old entries to new storage
      rowAlloc = rowNNZ+toAdd;
      ArrayRCP<Scalar> newVals = Teuchos::arcp<Scalar>(rowAlloc);
      std::copy(values_[localRow].begin(), values_[localRow].begin()+rowNNZ, newVals.begin());
      values_[localRow] = newVals;
    }
    // insert indices and values
    graph_->insertMyIndices(localRow,findices);  // this will enlarge allocation and append new indices
    std::copy( fvalues.begin(), fvalues.end(),  getVptr(localRow)+rowNNZ);
#ifdef HAVE_TPETRA_DEBUG
    // the assumption is that graph_->numAllocatedEntriesForMyRow is the allocated size here
    TEST_FOR_EXCEPTION( rowAlloc != graph_->numAllocatedEntriesForMyRow(localRow) 
                        || rowNNZ+toAdd != getNumEntriesInLocalRow(localRow), std::logic_error,
                        Teuchos::typeName(*this) << "::insertMyValues(): Internal logic error or unsupported use case. Please contact Tpetra team.");
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::insertGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>  &values) {
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(graph_->isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalValues(): graph indices are local; use insertMyValues().");
    TEST_FOR_EXCEPTION(isStaticGraph() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalValues(): matrix was constructed with static graph. Cannot insert new entries.");
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::insertGlobalValues(): values.size() must equal indices.size().");

    LocalOrdinal myRow = getRowMap().getLocalIndex(globalRow);

    // if we have a column map, use it to filter the entries.
    // only filter if this is our row.
    Teuchos::Array<GlobalOrdinal> finds;
    Teuchos::Array<Scalar>        fvals;
    Teuchos::ArrayView<const GlobalOrdinal> findices = indices;
    Teuchos::ArrayView<const Scalar       > fvalues  = values;
    if (hasColMap() && myRow != LOT::invalid()) {
      const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &cmap = getColMap();
      for (size_t i=0; i<indices.size(); ++i) {
        if (cmap.isNodeGlobalElement(indices[i])) {
          finds.push_back(indices[i]);
          fvals.push_back(values[i]);
        }
      }
      findices = finds();
      fvalues  = fvals();
    }

    // add the new indices and values
    if (myRow != LOT::invalid()) {
      size_t rowNNZ = getNumEntriesInLocalRow(myRow),
                      toAdd = findices.size(),
                   rowAlloc = graph_->numAllocatedEntriesForMyRow(myRow);
      if (rowNNZ+toAdd > rowAlloc) {
        TEST_FOR_EXCEPTION(graph_->isStaticProfile() == true, std::runtime_error,
            Teuchos::typeName(*this) << "::insertGlobalValues(): new indices exceed statically allocated graph structure.");
        TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
            "::insertGlobalValues(): Pre-allocated space has been exceeded, requiring new allocation. To improve efficiency, suggest larger allocation.");
#ifdef HAVE_TPETRA_DEBUG
        // assumption: the number currently allocated in values_ is that stored in the graph.
        TEST_FOR_EXCEPTION( rowAlloc != values_[myRow].size(), std::logic_error,
            Teuchos::typeName(*this) << "::insertGlobalValues(): Internal logic error or unsupported use case. Please contact Tpetra team.");
#endif
        // increase the allocation, copy old entries to new storage
        rowAlloc = rowNNZ+toAdd;
        ArrayRCP<Scalar> newVals = Teuchos::arcp<Scalar>(rowAlloc);
        std::copy(values_[myRow].begin(), values_[myRow].begin()+rowNNZ, newVals.begin());
        values_[myRow] = newVals;
      }
      // insert indices and values
      graph_->insertGlobalIndices(globalRow,findices);  // this will enlarge allocation and append new indices
      std::copy( fvalues.begin(), fvalues.end(),  getVptr(myRow)+rowNNZ);
#ifdef HAVE_TPETRA_DEBUG
      // the assumption is that graph_->numAllocatedEntriesForMyRow is the allocated size here as well
      TEST_FOR_EXCEPTION( rowAlloc != graph_->numAllocatedEntriesForMyRow(myRow) 
                          || rowNNZ+toAdd != getNumEntriesInLocalRow(myRow), std::logic_error,
          Teuchos::typeName(*this) << "::insertGlobalValues(): Internal logic error or unsupported use case. Please contact Tpetra team.");
#endif
    }
    else {
      typename Teuchos::ArrayView<const GlobalOrdinal>::iterator ind = findices.begin();
      typename Teuchos::ArrayView<const Scalar       >::iterator val =  fvalues.begin();
      for (; val != fvalues.end(); ++val, ++ind) {
        nonlocals_[globalRow].push_back(std::make_pair(*ind, *val));
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::replaceGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>  &values) {
    // find the values for the specified indices
    // if the row is not ours, throw an exception
    // ignore values not in the matrix (indices not found)
    // operate whether indices are local or global
    const size_t TOINV = Teuchos::OrdinalTraits<size_t>::invalid();
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::replaceGlobalValues(): values.size() must equal indices.size().");
    typename Teuchos::ArrayView<const GlobalOrdinal>::iterator ind = indices.begin();
    typename Teuchos::ArrayView<const        Scalar>::iterator val = values.begin();
    LocalOrdinal lrow = getRowMap().getLocalIndex(globalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::replaceGlobalValues(): specified global row does not belong to this processor.");
    typename Teuchos::ArrayRCP<Scalar>::iterator vptr = getVptr(lrow);
    if (isLocallyIndexed() == true) {
      while (ind != indices.end()) {
        LocalOrdinal lind = getColMap().getLocalIndex(*ind);
        size_t loc = graph_->findMyIndex(lrow,lind);
        if (loc != TOINV) {
          vptr[loc] = (*val);
        }
        ++ind;
        ++val;
      }
    }
    else if (isGloballyIndexed() == true) {
      while (ind != indices.end()) {
        size_t loc = graph_->findGlobalIndex(lrow,*ind);
        if (loc != TOINV) {
          vptr[loc] = (*val);
        }
        ++ind;
        ++val;
      }
    }
    //else {
      // indices are not allocated; nothing to do
    //}
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::sumIntoGlobalValues(GlobalOrdinal globalRow, 
                         const Teuchos::ArrayView<const GlobalOrdinal> &indices,
                         const Teuchos::ArrayView<const Scalar>  &values) {
    // find the values for the specified indices
    // if the row is not ours, throw an exception
    // ignore values not in the matrix (indices not found)
    // operate whether indices are local or global
    const size_t TOINV = Teuchos::OrdinalTraits<size_t>::invalid();
    TEST_FOR_EXCEPTION(values.size() != indices.size(), std::runtime_error,
        Teuchos::typeName(*this) << "::sumIntoGlobalValues(): values.size() must equal indices.size().");
    typename Teuchos::ArrayView<const GlobalOrdinal>::iterator ind = indices.begin();
    typename Teuchos::ArrayView<const        Scalar>::iterator val = values.begin();
    LocalOrdinal lrow = getRowMap().getLocalIndex(globalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::sumIntoGlobalValues(): specified global row does not belong to this processor.");
    typename Teuchos::ArrayRCP<Scalar>::iterator vptr = getVptr(lrow);
    if (isLocallyIndexed() == true) {
      while (ind != indices.end()) {
        LocalOrdinal lind = getColMap().getLocalIndex(*ind);
        size_t loc = graph_->findMyIndex(lrow,lind);
        if (loc != TOINV) {
          vptr[loc] += (*val);
        }
        ++ind;
        ++val;
      }
    }
    else if (isGloballyIndexed() == true) {
      while (ind != indices.end()) {
        size_t loc = graph_->findGlobalIndex(lrow,*ind);
        if (loc != TOINV) {
          vptr[loc] += (*val);
        }
        ++ind;
        ++val;
      }
    }
    //else {
      // indices are not allocated; nothing to do
    //}
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getLocalRowCopy(LocalOrdinal LocalRow, 
                                                                     const Teuchos::ArrayView<LocalOrdinal> &indices, 
                                                                     const Teuchos::ArrayView<Scalar>       &values,
                                                                     size_t &numEntries) const {
    numEntries = getNumEntriesInLocalRow(LocalRow);
    TEST_FOR_EXCEPTION(getRowMap().isNodeLocalElement(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowCopy(LocalRow,...): specified row (==" << LocalRow << ") is not valid on this node.");
    TEST_FOR_EXCEPTION(indices.size() < numEntries || values.size() < numEntries, std::runtime_error, 
        Teuchos::typeName(*this) << "::getLocalRowCopy(LocalRow,indices,values): size of indices,values must be sufficient to store the specified row.");
#ifdef HAVE_TPETRA_DEBUG
    size_t nnzagain;
    graph_->getLocalRowCopy(LocalRow,indices,nnzagain);
    TEST_FOR_EXCEPTION(nnzagain != numEntries, std::logic_error, 
        Teuchos::typeName(*this) << "::getLocalRowCopy(): Internal logic error. Please contact Tpetra team.");
#else
    graph_->getLocalRowCopy(LocalRow,indices,numEntries);
#endif
    typename Teuchos::ArrayRCP<const Scalar>::iterator vptr = getVptr(LocalRow);
    std::copy( vptr, vptr+numEntries, values.begin() );
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalRowCopy(GlobalOrdinal globalRow, 
                                                                      const Teuchos::ArrayView<GlobalOrdinal> &indices,
                                                                      const Teuchos::ArrayView<Scalar>        &values,
                                                                      size_t &numEntries) const {
    // Only locally owned rows can be queried, otherwise complain
    size_t myRow = getRowMap().getLocalIndex(globalRow);
    TEST_FOR_EXCEPTION(myRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowCopy(globalRow,...): globalRow does not belong to this node.");
    numEntries = getNumEntriesInLocalRow(myRow);
    TEST_FOR_EXCEPTION(
        indices.size() < numEntries || values.size() < numEntries, std::runtime_error, 
        Teuchos::typeName(*this) << "::getGlobalRowCopy(globalRow,indices,values): size of indices,values must be sufficient to store the specified row.");
#ifdef HAVE_TPETRA_DEBUG
    size_t nnzagain;
    graph_->getGlobalRowCopy(globalRow,indices,nnzagain);
    TEST_FOR_EXCEPTION(nnzagain != numEntries, std::logic_error, 
        Teuchos::typeName(*this) << "::getLocalRowCopy(): Internal logic error. Please contact Tpetra team.");
#else
    graph_->getGlobalRowCopy(globalRow,indices,numEntries);
#endif
    typename Teuchos::ArrayRCP<const Scalar>::iterator vptr = getVptr(myRow);
    std::copy( vptr, vptr+numEntries, values.begin() );
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::extractGlobalRowView(GlobalOrdinal GlobalRow, 
                                Teuchos::ArrayView<GlobalOrdinal> &indices, 
                                Teuchos::ArrayView<Scalar> &values) {
    TEST_FOR_EXCEPTION(isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::extractGlobalRowView(): global indices do not exist; call getLocalRowViewNonConst().");
    size_t lrow = getRowMap().getLocalIndex(GlobalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::extractGlobalRowView(GlobalRow,...): GlobalRow (== " << GlobalRow << ") does not belong to this node.");
    size_t rnnz = getNumEntriesInLocalRow(lrow);
    graph_->extractGlobalRowView(GlobalRow,indices);
    if (rnnz == 0) {
      values = Teuchos::ArrayView<Scalar>(Teuchos::null);
    }
    else {
      if (isStorageOptimized() == true || graph_->isStaticProfile() == true) {
        values = Teuchos::arrayView<Scalar>(Teuchos::getRawPtr(valuesPtrs_[lrow]),rnnz);
      }
      else {
        values = values_[lrow](0,rnnz);
      }
    }
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getLocalRowViewNonConst(
                                LocalOrdinal LocalRow, 
                                Teuchos::ArrayRCP<LocalOrdinal> &indices, 
                                Teuchos::ArrayRCP<Scalar> &values) {
    TEST_FOR_EXCEPTION(isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowViewNonConst(): local indices do not exist; call extractGlobalRowVie().");
    TEST_FOR_EXCEPTION(getRowMap().isNodeLocalElement(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowViewNonConst(LocalRow,...): LocalRow (== " << LocalRow << ") is not valid on this node.");
    size_t rnnz = getNumEntriesInLocalRow(LocalRow);
    indicies = graph_->getLocalRowViewNonConst(LocalRow);
    if (rnnz == 0) {
      values = Teuchos::ArrayView<Scalar>(Teuchos::null);
    }
    else {
      if (isStorageOptimized() == true || graph_->isStaticProfile() == true) {
        values = Teuchos::arrayView<Scalar>(Teuchos::getRawPtr(valuesPtrs_[LocalRow]),rnnz);
      }
      else {
        values = values_[LocalRow](0,rnnz);
      }
    }
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalRowView(
                                GlobalOrdinal GlobalRow, 
                                Teuchos::ArrayRCP<const GlobalOrdinal> &indices,
                                Teuchos::ArrayrCP<const Scalar>        &values) const {
    TEST_FOR_EXCEPTION(isLocallyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowView(): global indices do not exist; call getLocalRowView().");
    size_t lrow = getRowMap().getLocalIndex(GlobalRow);
    TEST_FOR_EXCEPTION(lrow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        Teuchos::typeName(*this) << "::getGlobalRowView(GlobalRow,...): GlobalRow (== " << GlobalRow << ") does not belong to this node.");
    size_t rnnz = getNumEntriesInLocalRow(lrow);
    graph_->getGlobalRowView(GlobalRow,indices);
    if (rnnz == 0) {
      values = Teuchos::ArrayView<Scalar>(Teuchos::null);
    }
    else {
      if (isStorageOptimized() == true || graph_->isStaticProfile() == true) {
        values = Teuchos::arrayView<Scalar>(Teuchos::getRawPtr(valuesPtrs_[lrow]),rnnz);
      }
      else {
        values = values_[lrow](0,rnnz);
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
                                Teuchos::ArrayRCP<const Scalar> &values) const {
    TEST_FOR_EXCEPTION(isGloballyIndexed() == true, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowView(): local indices do not exist; call getGlobalRowView().");
    TEST_FOR_EXCEPTION(getRowMap().isNodeLocalElement(LocalRow) == false, std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalRowView(LocalRow,...): LocalRow (== " << LocalRow << ") is not valid on this node.");
    size_t rnnz = getNumEntriesInLocalRow(LocalRow);
    graph_->getLocalRowView(LocalRow,indices);
    if (rnnz == 0) {
      values = Teuchos::ArrayView<Scalar>(Teuchos::null);
    }
    else {
      if (isStorageOptimized() == true || graph_->isStaticProfile() == true) {
        values = Teuchos::arrayView<Scalar>(Teuchos::getRawPtr(valuesPtrs_[LocalRow]),rnnz);
      }
      else {
        values = values_[LocalRow](0,rnnz);
      }
    }
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
                                                                 MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                                                                 Teuchos::ETransp mode) const {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    using Teuchos::null;
    using Teuchos::ArrayView;
    TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error, 
        Teuchos::typeName(*this) << ": cannot call apply() until fillComplete() has been called.");
    TEST_FOR_EXCEPTION(X.numVectors() != Y.numVectors(), std::runtime_error,
        Teuchos::typeName(*this) << "::apply(X,Y): X and Y must have the same number of vectors.");

#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
    int myImageID = Teuchos::rank(*getComm());
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    if (myImageID == 0) {
      *out << "Entering CrsMatrix::apply()" << std::endl
                << "Column Map: " << std::endl;
    }
    *out << this->getColMap() << std::endl;
    if (myImageID == 0) {
      *out << "Initial input: " << std::endl;
    }
    X.print(*out); X.printValues(*out);
#   endif

    const size_t numVectors = X.numVectors();
    // because of Views, it is difficult to determine if X and Y point to the same data. 
    // however, if they reference the exact same object, we will do the user the favor of copying X into new storage (with a warning)
    // we ony need to do this if we have trivial importers; otherwise, we don't actually apply the operator from X into Y
    Teuchos::RCP<const Import<LocalOrdinal,GlobalOrdinal> > importer = graph_->getImporter();
    Teuchos::RCP<const Export<LocalOrdinal,GlobalOrdinal> > exporter = graph_->getExporter();
    Teuchos::RCP<const MV> Xcopy;
    typename MV::const_double_pointer Xdata = X.extractConstView2D();
    typename MV::double_pointer       Ydata = Y.extractView2D();
    if (&X==&Y && importer==null && exporter==null) {
      TPETRA_EFFICIENCY_WARNING(true,std::runtime_error,
          "::apply(X,Y): If X and Y are the same, it necessitates a temporary copy of X, which is inefficient.");
      // generate a copy of X 
      Xcopy = Teuchos::rcp(new MV(X));
      Xdata = Xcopy->extractConstView2D();
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "X and Y are co-located, duplicating X results in a stride copy" << std::endl;
      *out << this->getColMap() << std::endl;
      Xcopy->print(*out); Xcopy->printValues(*out);
#   endif
    }
    if (importer != null) {
      if (importMV_ != null && importMV_->numVectors() != numVectors) importMV_ = null;
      if (importMV_ == null) {
        importMV_ = Teuchos::rcp( new MV(getColMap(),numVectors) );
      }
    }
    if (exporter != null) {
      if (exportMV_ != null && exportMV_->numVectors() != numVectors) exportMV_ = null;
      if (exportMV_ == null) {
        exportMV_ = Teuchos::rcp( new MV(getRowMap(),numVectors) );
      }
    }
    if (mode == Teuchos::NO_TRANS) {
      // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
      if (importer != null) {
        importMV_->doImport(X, *importer, INSERT);
        Xdata = importMV_->extractConstView2D();
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) {
          *out << "Performed import of X using importer..." << std::endl;
        }
        importMV_->print(*out); importMV_->printValues(*out);
#   endif
      }
      // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
      // We will compute solution into the to-be-exported MV; get a view
      if (exporter != null) {
        Ydata = exportMV_->extractView2D();
      }
      // Do actual computation
      if (numVectors==1) {
        GeneralMV(Xdata[0],Ydata[0]);
      }
      else {
        GeneralMM(Xdata,Ydata,numVectors);
      }
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "Matrix-MV product..." << std::endl;
      if (exportMV_ != null) {
        exportMV_->print(*out); exportMV_->printValues(*out);
      } else {
        Y.print(*out); Y.printValues(*out);
      }
#   endif
      // do the export
      if (exporter != null) {
        Y.putScalar(ST::zero());  // Make sure target is zero: necessary because we are adding. may need adjusting for alpha,beta apply()
        Y.doExport(*exportMV_, *exporter, ADD); // Fill Y with Values from export vector
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) *out << "Output vector after export() using exporter..." << std::endl;
        Y.print(*out); Y.printValues(*out);
#   endif
      }
      // Handle case of rangemap being a local replicated map: in this case, sum contributions from each processor
      if (Y.isDistributed() == false) {
        Y.reduce();
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) *out << "Output vector is local; result after reduce()..." << std::endl;
        Y.print(*out); Y.printValues(*out);
#   endif
      }
    }
    else {
      // mode == CONJ_TRANS or TRANS
      TEST_FOR_EXCEPTION(Teuchos::ScalarTraits<Scalar>::isComplex && mode == Teuchos::TRANS, std::logic_error,
          Teuchos::typeName(*this) << "::apply() does not currently support transposed multiplications for complex scalar types.");
      // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors
      if (exporter != null) {
        exportMV_->doImport(X,*exporter,INSERT);
        Xdata = exportMV_->extractConstView2D();
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) {
          *out << "Performed import of X using exporter..." << std::endl;
        }
        exportMV_->print(*out); exportMV_->printValues(*out);
#   endif
      }
      // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
      // We will compute colutioni into the to-be-exported MV; get a view
      if (importer != null) {
        Ydata = importMV_->extractView2D();
      }
      // Do the actual transposed multiplication
      if (numVectors==1) {
        GeneralMhV(Xdata[0],Ydata[0]);
      }
      else {
        GeneralMhM(Xdata,Ydata,numVectors);
      }
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
      if (myImageID == 0) *out << "Matrix-MV product..." << std::endl;
      if (importMV_ != null) {
        importMV_->print(*out); importMV_->printValues(*out);
      } else {
        Y.print(*out); Y.printValues(*out);
      }
#   endif
      if (importer != null) {
        Y.putScalar(ST::zero()); // Make sure target is zero: necessary because we are adding. may need adjusting for alpha,beta apply()
        Y.doExport(*importMV_,*importer,ADD);
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) *out << "Output vector after export() using importer..." << std::endl;
        Y.print(*out); Y.printValues(*out);
#   endif
      }
      // Handle case of rangemap being a local replicated map: in this case, sum contributions from each processor
      if (Y.isDistributed() == false) {
        Y.reduce();
#   ifdef TPETRA_CRSMATRIX_MULTIPLY_DUMP
        if (myImageID == 0) *out << "Output vector is local; result after reduce()..." << std::endl;
        Y.print(*out); Y.printValues(*out);
#   endif
      }
    }
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

    if (!isStaticGraph()) {
      graph_->makeIndicesLocal(domainMap,rangeMap);
    }

    sortEntries();
    mergeRedundantEntries();

    if (!isStaticGraph()) {
      graph_->fillComplete(domainMap,rangeMap,os);
    }
  
    fillComplete_ = true;

    if (os == DoOptimizeStorage) optimizeStorage();
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
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::optimizeStorage() {
    // optimizeStorage will perform two functions:
    // 1) create a single allocation of memory
    // 2) pack data in that allocation
    // if isStaticProfile() == true, then 1) has already been done
    if (isStorageOptimized() == true) return;
    TEST_FOR_EXCEPTION(isFillComplete() == false || graph_->indicesAreSorted() == false || graph_->noRedundancies() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::optimizeStorage(): fillComplete() must be called before optimizeStorage().");
    
    // 1) allocate single memory block
    const size_t nlrs = getNodeNumRows();
    if (nlrs > 0) {
      if (graph_->isStaticProfile() == false) {
        // need to allocate storage, create pointers, copy data, and delete old data
        const size_t nE = getNodeNumEntries();
        valuesPtrs_ = Teuchos::arcp<typename Teuchos::ArrayRCP<Scalar>::iterator>(nlrs+1);
        if (nE > 0) {
          contigValues_ = Teuchos::arcp<Scalar>(nE);
          size_t sofar = 0;
          for (size_t r=0; r<nlrs; ++r) {
            size_t rne = graph_->RNNZ(r);
            valuesPtrs_[r] = contigValues_.persistingView(sofar,rne).begin();
            std::copy(values_[r].begin(),values_[r].begin()+rne,valuesPtrs_[r]);
            values_[r] = Teuchos::null;
            sofar += rne;
          }
          valuesPtrs_[nlrs] = contigValues_.end();
        }
        else {
          std::fill(valuesPtrs_.begin(),valuesPtrs_.end(),
                    Teuchos::NullIteratorTraits<typename Teuchos::ArrayRCP<Scalar>::iterator>::getNull());
        }
        values_ = Teuchos::null;
      }
      else {
        // storage is already allocated and pointers are set; just need to pack
        // need to change pointers and move data
        // remember, these aren't just pointers, but also bounds-checked 
        const size_t nE = getNodeNumEntries();
        if (nE > 0) {
          size_t sofar = 0;
          typename Teuchos::ArrayRCP<Scalar>::iterator newptr, oldptr;
          for (size_t r=0; r<nlrs; ++r) {
            size_t rne = graph_->RNNZ(r);
            newptr = contigValues_.persistingView(sofar,rne).begin();
            oldptr = valuesPtrs_[r];
            if (newptr != oldptr) {
              for (size_t j=0; j<rne; ++j) {
                newptr[j] = oldptr[j];
              }
              valuesPtrs_[r] = newptr;
            }
            sofar += rne;
          }
          valuesPtrs_[nlrs] = contigValues_.persistingView(sofar,0).begin();
        }
      }
    }
    storageOptimized_ = true;
#ifdef HAVE_TPETRA_DEBUG
    // check that we only have the data structures that we need
    TEST_FOR_EXCEPTION( values_ != Teuchos::null, std::logic_error,
        Teuchos::typeName(*this) << "::optimizeStorage(): Internal logic error. Please contact Tpetra team.");
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
    Teuchos::reduceAll<size_t>(*getComm(),Teuchos::REDUCE_MAX,MyNonlocals,&MaxGlobalNonlocals);
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
        LookupStatus stat = getRowMap().getRemoteIndexList(NLRs(),NLRIds());
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
    for (int s=0; s < numSends ; ++s)
    {
      // we'll fake the memory management, because all communication will be local to this method and the scope of our data
      sendRequests.push_back( Teuchos::isend<int,int>(*getComm(),rcp(&sendSizes[s],false),sendIDs[s]) );
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


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal> &dvec) const {
    TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
        Teuchos::typeName(*this) << ": cannot call getLocalDiagCopy() until fillComplete() has been called.");
    TEST_FOR_EXCEPTION(!dvec.getMap().isSameAs(getRowMap()), std::runtime_error,
        Teuchos::typeName(*this) << "::getLocalDiagCopy(dvec): dvec must have the same map as the CrsMatrix.");
#ifdef HAVE_TPETRA_DEBUG
    int numDiagFound = 0;
#endif
    Teuchos::ArrayView<Scalar> values;
    dvec.extractView1D(values);
    size_t nlrs = getNodeNumRows();
    typename Teuchos::ArrayRCP<const Scalar>::iterator v;
    typename Teuchos::ArrayView<Scalar>::iterator ov;
    Teuchos::ArrayView<const LocalOrdinal> cinds;
    typename Teuchos::ArrayView<const LocalOrdinal>::iterator i;
    ov = values.begin();
    for (size_t r=0; r < nlrs; ++r) {
      *ov = Teuchos::ScalarTraits<Scalar>::zero();
      GlobalOrdinal rgid = getRowMap().getGlobalIndex(r);
      if (getColMap().isNodeGlobalElement(rgid)) {
        LocalOrdinal rlid = getColMap().getLocalIndex(rgid);
        graph_->getLocalRowView(r,cinds);
        i = cinds.begin();
        v = getVptr(r);
        while (i != cinds.end()) {
          if (*i == rlid) {
            *ov = *v;
#ifdef HAVE_TPETRA_DEBUG
            ++numDiagFound;
#endif
            break;
          }
          else if (*i > rlid) {
            // indices are sorted; we have passed the diagonal, so quit
            break;
          }
          v++;
          i++;
        }
      }
      ++ov;
    }
#ifdef HAVE_TPETRA_DEBUG
    TEST_FOR_EXCEPTION(numDiagFound != graph_->getNodeNumDiags(), std::logic_error, 
        "CrsMatrix::getLocalDiagCopy(): logic error. Please contact Tpetra team.");
#endif
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::sortEntries() {
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(graph_->isGloballyIndexed() == true, std::logic_error,
        Teuchos::typeName(*this) << "::sortEntries(): sortEntries() must be called after indices are transformed to local.\n"
        << "Likely internal logic error. Please contact Tpetra team.");
    if (graph_->indicesAreSorted()) return;
    if (graph_->indicesAreAllocated()) {
      const size_t nlrs = getNodeNumRows();
      for (size_t r=0; r < nlrs; ++r)
      {
        ArrayView<LocalOrdinal> inds;
        typename ArrayRCP<Scalar>::iterator vals = getVptr(r);
        graph_->extractMyRowView(r,inds);
        sort2(inds.begin(), inds.end(), vals);
      }
    }
    graph_->indicesAreSorted(true);  // we just sorted them
    return;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::mergeRedundantEntries() {
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    TEST_FOR_EXCEPTION(graph_->indicesAreSorted() == false, std::runtime_error,
        Teuchos::typeName(*this) << "::mergeRedundantEntries() cannot be called before indices are sorted.\n"
        << "Likely interal logic error. Please contact Tpetra team.");
    if (graph_->noRedundancies()) return;
    for (size_t r=0; r<getNodeNumRows(); ++r) 
    {
      size_t rnnz = getNumEntriesInLocalRow(r);
      if (rnnz > 1) {
        ArrayView<const LocalOrdinal> inds;
        graph_->getLocalRowView(r,inds);
        typename ArrayRCP<Scalar>::iterator vals = getVptr(r);
        size_t curEntry = 0;
        Scalar curValue = vals[0];
        for (size_t k=1; k<rnnz; ++k) {
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
    }
    graph_->removeRedundantIndices();
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::scale(const Scalar &alpha) {
    TEST_FOR_EXCEPT(isStorageOptimized());
    for (size_t r=0; r<getNodeNumRows(); ++r) {
      typename Teuchos::ArrayRCP<Scalar>::iterator val, vend;
      val = getVptr(r);
      vend = val+getNumEntriesInLocalRow(r);
      while (val != vend) {
        (*val++) *= alpha;
      }
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
  void CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::print(std::ostream& os) const {
    // support this operation whether fillComplete() or not
    using std::endl;
    int myImageID = Teuchos::rank(*getComm());
    if (myImageID == 0)
    {
      os << "Tpetra::CrsMatrix, label = " << this->getObjectLabel() << endl;
      os << "Number of global rows    = " << getGlobalNumRows() << endl;
      if (isFillComplete())
      {
        os << "Number of global columns    = " << getGlobalNumCols() << endl;
        os << "Status = fill complete" << endl;
        os << "Number of global nonzeros   = " << getGlobalNumEntries() << endl;
        os << "Global max nonzeros per row = " << getGlobalMaxNumRowEntries() << endl;
      }
      else
      {
        os << "Status = fill not complete" << endl;
      }
    }
    if (isFillComplete())
    {
      for (int pid=0; pid < Teuchos::size(*getComm()); ++pid)
      {
        if (pid == myImageID)
        {
          Teuchos::Array<GlobalOrdinal> indices(getNodeMaxNumRowEntries());
          Teuchos::Array<Scalar>         values(getNodeMaxNumRowEntries());
          size_t rowSize;
          os << "% Number of rows on image " << myImageID << " = " << getNodeNumRows() << endl;
          for (size_t i=0; i < getNodeNumRows(); ++i)
          {
            GlobalOrdinal globalRow = getRowMap().getGlobalIndex(i);
            getGlobalRowCopy(globalRow, indices(), values(), rowSize);
            if (rowSize > Teuchos::OrdinalTraits<size_t>::zero()) {
              for (size_t j=0; j < rowSize; ++j) {
                os << "Matrix(" << globalRow << ", " << indices[j] << ") = " << values[j] << endl;
              }
            }
          }
        }
        Teuchos::barrier(*getComm());
        Teuchos::barrier(*getComm());
        Teuchos::barrier(*getComm());
      }
    }
  }


} // namespace Tpetra

#endif
