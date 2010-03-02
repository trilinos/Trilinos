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

#ifndef TPETRA_CRSMATRIX_DECL_HPP
#define TPETRA_CRSMATRIX_DECL_HPP

// TODO: row-wise insertion of entries in globalAssemble() may be more efficient

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_CrsMatrix.hpp>
#include <Kokkos_DefaultSparseMultiply.hpp>
#include <Kokkos_DefaultSparseSolve.hpp>

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>

#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrixMultiplyOp_decl.hpp"

namespace Tpetra {
  // struct for i,j,v triplets
  template <class Ordinal, class Scalar>
  struct CrsIJV {
    CrsIJV();
    CrsIJV(Ordinal row, Ordinal col, const Scalar &val);
    Ordinal i,j;
    Scalar  v;
  };
}

namespace Teuchos {
  // SerializationTraits specialization for CrsIJV, using DirectSerialization
  template <typename Ordinal, typename Scalar>
  class SerializationTraits<int,Tpetra::CrsIJV<Ordinal,Scalar> >
  : public DirectSerializationTraits<int,Tpetra::CrsIJV<Ordinal,Scalar> >
  {};
}

namespace std {
  template <class Ordinal, class Scalar>
  bool operator<(const Tpetra::CrsIJV<Ordinal,Scalar> &ijv1, const Tpetra::CrsIJV<Ordinal,Scalar> &ijv2);
}

namespace Tpetra {

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
      typedef Scalar        scalar_type;
      typedef LocalOrdinal  local_ordinal_type;
      typedef GlobalOrdinal global_ordinal_type;
      typedef Node          node_type;
      typedef LocalMatVec   mat_vec_type;
      typedef LocalMatSolve mat_solve_type;

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
      Teuchos::RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > getCrsGraph() const;

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
      void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const;

      //@}

      //! @name Advanced Matrix-vector multiplication and solve methods
      //@{

      //! Multiplies this matrix by a MultiVector.
      /*! \c X is required to be post-imported, i.e., described by the column map of the matrix. \c Y is required to be pre-exported, i.e., described by the row map of the matrix.

          Both are required to have constant stride, and they are not permitted to ocupy overlapping space. No runtime checking will be performed in a non-debug build.
          
          If \c beta is equal to zero, the operation will enjoy overwrite semantics (\c Y will be overwritten with the result of the multiplication). Otherwise, the result of the multiplication
          will be accumulated into \c Y.
       */
      template <class DomainScalar, class RangeScalar>
      void multiply(const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, RangeScalar alpha, RangeScalar beta) const;

      //! Solves a linear system when the underlying matrix is triangular.
      /*! \c X is required to be post-imported, i.e., described by the column map of the matrix. \c Y is required to be pre-exported, i.e., described by the row map of the matrix.

          Both are required to have constant stride. However, unlike multiply(), it is permissible for <tt>&X == &Y</tt>. No runtime checking will be performed in a non-debug build.

          If \c beta is equal to zero, the operation will enjoy overwrite semantics (\c Y will be overwritten with the result of the multiplication). Otherwise, the result of the multiplication
          will be accumulated into \c Y.
       */
      template <class DomainScalar, class RangeScalar>
      void solve(const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> & Y, MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X, Teuchos::ETransp trans) const;
          
      //@}

      //! @name Methods implementing Operator
      //@{ 

      //! \brief Computes the sparse matrix-multivector multiplication.
      /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
          - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
       */
      void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                 Teuchos::ETransp mode = Teuchos::NO_TRANS,
                 Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                 Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

      //! Indicates whether this operator supports applying the adjoint operator.
      bool hasTransposeApply() const;

      //! \brief Returns the Map associated with the domain of this operator.
      //! This will be <tt>Teuchos::null</tt> until fillComplete() is called.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const;

      //! Returns the Map associated with the domain of this operator.
      //! This will be <tt>Teuchos::null</tt> until fillComplete() is called.
      const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const;

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
      typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

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

      Teuchos::RCP<CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > graph_;

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
      Teuchos::ArrayRCP<Scalar>                       pbuf_values1D_, view_values1D_;
      Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >   pbuf_values2D_, view_values2D_;

      // a map between a (non-local) row and a list of (col,val)
      std::map<GlobalOrdinal, std::list<std::pair<GlobalOrdinal,Scalar> > > nonlocals_;

      // a wrapper around multiply, for use in apply; it contains a non-owning RCP to *this, therefore, it is not allowed 
      // to persist past the destruction of *this. therefore, we may not share it.
      Teuchos::RCP< const CrsMatrixMultiplyOp<Scalar,Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> > sameScalarMultiplyOp_;
  }; // class CrsMatrix

}

#endif
