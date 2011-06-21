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

// TODO: add typeglobs: CrsMatrix<Scalar,typeglob>
// TODO: add template (template) parameter for nonlocal container (this will be part of typeglob)

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
#include <Kokkos_CrsMatrix.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrixMultiplyOp_decl.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
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
#endif

namespace Tpetra {

  //! \brief A class for constructing and using sparse compressed matrices with row access.
  /*!
   \tparam Scalar        The scalar field describing the numerical entries of the matrix.
   \tparam LocalOrdinal  A ordinal type for lists of local indices. This specifies the \c LocalOrdinal type for Map objects used by this matrix.
   \tparam GlobalOrdinal A ordinal type for lists of global indices. This specifies the \c GlobalOrdinal type for Map objects used by this matrix.
   \tparam Node          A shared-memory node class, fulfilling the \ref kokkos_node_api "Kokkos Node API"
   \tparam LocalMatOps   A local sparse matrix operations class, fulfiling the \ref kokkos_crs_ops "Kokkos CRS Ops API".
   * This class allows the construction of sparse matrices with row-access. 
   * 
   * <b>Local vs. Global</b>
   * 
   * Matrix entries can be added using either local or global coordinates for the indices. The 
   * accessors isGloballyIndexed() and isLocallyIndexed() indicate whether the indices are currently
   * stored as global or local indices. Many of the class methods are divided into global and local 
   * versions, which differ only in whether they accept/return indices in the global or local coordinate
   * space. Some of these methods may only be used if the matrix coordinates are in the appropriate coordinates.
   * For example, getGlobalRowView() returns a View to the indices in global coordinates; if the indices are 
   * not in global coordinates, then no such View can be created.
   * 
   * The global/local distinction does distinguish between operation on the global/local matrix. Almost all methods 
   * operate on the local matrix, i.e., the rows of the matrix associated with the local node, per the distribution specified
   * by the row map. Access to non-local rows requires performing an explicit communication via the import/export capabilities of the
   * CrsMatrix object; see DistObject. However, the method insertGlobalValues() is an exception to this rule, as non-local rows are 
   * allowed to be added via the local matrix. These rows are stored in the local matrix and communicated to the appropriate node 
   * on the next call to globalAssemble() or fillComplete() (the latter calls the former).
   * 
   */
  template <class Scalar, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps >
  class CrsMatrix : public RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>,
                    public DistObject<char, LocalOrdinal,GlobalOrdinal,Node> {
    public:
      typedef Scalar        scalar_type;
      typedef LocalOrdinal  local_ordinal_type;
      typedef GlobalOrdinal global_ordinal_type;
      typedef Node          node_type;
      // backwards compatibility defines both of these
      typedef LocalMatOps   mat_vec_type;
      typedef LocalMatOps   mat_solve_type;

      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor specifying the number of non-zeros for all rows.
      CrsMatrix(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

      //! Constructor specifying the number of non-zeros for each row.
      CrsMatrix(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile);

      //! Constructor specifying a column map and the number of non-zeros for all rows.
      /** The column map will be used to filter any matrix entries inserted using insertLocalValues() or insertGlobalValues().
        */
      CrsMatrix(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, size_t maxNumEntriesPerRow, ProfileType pftype = DynamicProfile);

      //! Constructor specifying a column map and the number of non-zeros for each row.
      /** The column map will be used to filter any matrix entries inserted using insertLocalValues() or insertGlobalValues().
        */
      CrsMatrix(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &colMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile);

      //! Constructor specifying a pre-constructed graph.
      explicit CrsMatrix(const RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &graph);

      // !Destructor.
      virtual ~CrsMatrix();

      //@}

      //! @name Insertion/Removal Methods
      //@{ 

      //! Insert matrix entries, using global IDs.
      /** All index values must be in the global space. 
          \pre \c globalRow exists as an ID in the global row map
          \pre <tt>isStorageOptimized() == false</tt>

          \note If \c globalRow does not belong to the matrix on this node, then it will be communicated to the appropriate node when globalAssemble() is called (which will, at the latest, occur during the next call to fillComplete().) Otherwise, the entries will be inserted in the local matrix. 
          \note If the matrix row already contains values at the indices corresponding to values in \c cols, then the new values will be summed with the old values; this may happen at insertion or during the next call to fillComplete().
          \note If <tt>hasColMap() == true</tt>, only (cols[i],vals[i]) where cols[i] belongs to the column map on this node will be inserted into the matrix.
          \note If <tt>isLocallyIndexed() == true</tt>, then the global indices will be translated to local indices via the column map; indices not present in the column map will be discarded.
       */
      void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals);

      //! Insert matrix entries, using local IDs.
      /**
          \pre \c localRow is a local row belonging to the matrix on this node
          \pre <tt>isGloballyIndexed() == false</tt>
          \pre <tt>isStorageOptimized() == false</tt>
          \pre <tt>hasColMap() == true</tt>

          \post <tt>isLocallyIndexed() == true</tt>

          \note If the matrix row already contains entries at the indices corresponding to values in \c cols, then the new values will be summed with the old values; this may happen at insertion or during the next call to fillComplete().
          \note If <tt>hasColMap() == true</tt>, only (cols[i],vals[i]) where cols[i] belongs to the column map on this node will be inserted into the matrix.
        */
      void insertLocalValues(LocalOrdinal localRow, const ArrayView<const LocalOrdinal> &cols, const ArrayView<const Scalar> &vals);

      //! \brief Replace matrix entries, using global IDs.
      /** All index values must be in the global space. 

         \pre \c globalRow is a global row belonging to the matrix on this node.

         \note If (globalRow,cols[i]) corresponds to an entry that is duplicated in this matrix row (likely because it was inserted more than once and fillComplete() has not been called in the interim), the behavior of this function is not defined. */
      void replaceGlobalValues(GlobalOrdinal globalRow, 
                               const ArrayView<const GlobalOrdinal> &cols,
                               const ArrayView<const Scalar>        &vals);

      //! Replace matrix entries, using local IDs.
      /** All index values must be in the local space. 
        */
      void replaceLocalValues(LocalOrdinal localRow, 
                              const ArrayView<const LocalOrdinal> &cols,
                              const ArrayView<const Scalar>       &vals);

      //! Sum into multiple entries, using global IDs.
      /** All index values must be in the global space. 

         \pre \c globalRow is a global row belonging to the matrix on this node.

        */
      void sumIntoGlobalValues(GlobalOrdinal globalRow, 
                               const ArrayView<const GlobalOrdinal> &cols,
                               const ArrayView<const Scalar>        &vals);


      //! Sum into multiple entries, using local IDs.
      /** All index values must be in the local space. 

         \pre \c localRow is a local row belonging to the matrix on this node.

        */
      void sumIntoLocalValues(LocalOrdinal globalRow, 
                              const ArrayView<const LocalOrdinal>  &cols,
                              const ArrayView<const Scalar>        &vals); 

      //! Set all matrix entries equal to scalarThis.
      void setAllToScalar(const Scalar &alpha);

      //! Scale the current values of a matrix, this = alpha*this. 
      void scale(const Scalar &alpha);

      //@}
	  

      //! @name Transformational Methods
      //@{ 

      //! \brief Communicate non-local contributions to other nodes.
      void globalAssemble();

      /*! Resume fill operations.
          After calling fillComplete(), resumeFill() must be called before initiating any changes to the matrix.

          resumeFill() may be called repeatedly. 

          \post  <tt>isFillActive() == true<tt>
          \post  <tt>isFillComplete() == false<tt>
       */
      void resumeFill();

      /*! \brief Signal that data entry is complete, specifying domain and range maps.

          Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

          \pre  <tt>isFillActive() == true<tt>
          \pre <tt>isFillComplete()() == false<tt>

          \post <tt>isFillActive() == false<tt>
          \post <tt>isFillComplete() == true<tt>
          \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>. See isStorageOptimized() for consequences.
       */ 
      void fillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &domainMap, const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rangeMap, OptimizeOption os = DoOptimizeStorage);

      /*! \brief Signal that data entry is complete. 

          Off-node entries are distributed (via globalAssemble()), repeated entries are summed, and global indices are transformed to local indices.

          \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

          \pre  <tt>isFillActive() == true<tt>
          \pre <tt>isFillComplete()() == false<tt>

          \post <tt>isFillActive() == false<tt>
          \post <tt>isFillComplete() == true<tt>
          \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>. See isStorageOptimized() for consequences.
       */
      void fillComplete(OptimizeOption os = DoOptimizeStorage);

      //@}

      //! @name Methods implementing RowMatrix
      //@{ 

      //! Returns the communicator.
      const RCP<const Comm<int> > & getComm() const;

      //! Returns the underlying node.
      RCP<Node> getNode() const;

      //! Returns the Map that describes the row distribution in this matrix.
      const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRowMap() const;

      //! \brief Returns the Map that describes the column distribution in this matrix.
      const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getColMap() const;

      //! Returns the RowGraph associated with this matrix. 
      RCP<const RowGraph<LocalOrdinal,GlobalOrdinal,Node> > getGraph() const;

      //! Returns the CrsGraph associated with this matrix. 
      RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > getCrsGraph() const;

      /// \brief Number of global elements in the row map of this matrix.
      ///
      /// This is <it>not</it> the number of rows in the matrix as a
      /// mathematical object.  This method returns the global sum of
      /// the number of local elements in the row map on each
      /// processor, which is the row map's getGlobalNumElements().
      /// Since the row map is not one-to-one in general, that global
      /// sum could be different than the number of rows in the
      /// matrix.  If you want the number of rows in the matrix, ask
      /// the range map for its global number of elements, using the
      /// following code:
      /// <code>
      /// global_size_t globalNumRows = getRangeMap()->getGlobalNumElements();
      /// </code>
      /// This method retains the behavior of Epetra, which also asks
      /// the row map for the global number of rows, rather than
      /// asking the range map.
      ///
      /// \warning Undefined if isFillActive().
      /// 
      global_size_t getGlobalNumRows() const;

      /// \brief Number of global elements in the column map of this matrix.
      ///
      /// This is <it>not</it> the number of columns in the matrix as
      /// a mathematical object.  This method returns the global sum
      /// of the number of local elements in the column map on each
      /// processor, which is the column map's getGlobalNumElements().
      /// Since the column map is not one-to-one in general, that
      /// global sum could be different than the number of columns in
      /// the matrix.  If you want the number of rows in the matrix,
      /// ask the domain map for its global number of elements, using
      /// the following code:
      /// <code>
      /// global_size_t globalNumCols = getDomainMap()->getGlobalNumElements();
      /// </code>
      /// This method retains the behavior of Epetra, which also asks
      /// the column map for the global number of columns, rather than
      /// asking the domain map.
      ///
      /// \warning Undefined if isFillActive().
      ///
      global_size_t getGlobalNumCols() const;

      //! Returns the number of matrix rows owned on the calling node.
      size_t getNodeNumRows() const;

      //! Returns the number of columns connected to the locally owned rows of this matrix.
      /** Throws std::runtime_error if <tt>hasColMap() == false</tt>
        */
      size_t getNodeNumCols() const;

      //! Returns the index base for global indices for this matrix. 
      GlobalOrdinal getIndexBase() const;

      //! Returns the global number of entries in this matrix.
      global_size_t getGlobalNumEntries() const;

      //! Returns the local number of entries in this matrix.
      size_t getNodeNumEntries() const;

      //! \brief Returns the current number of entries on this node in the specified global row.
      /*! Returns OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this matrix. */
      size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const;

      //! Returns the current number of entries on this node in the specified local row.
      /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
      size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const;

      //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
      /** Undefined if isFillActive().
        */
      global_size_t getGlobalNumDiags() const;

      //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
      /** Undefined if isFillActive().
        */
      size_t getNodeNumDiags() const;

      //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
      /** Undefined if isFillActive().
        */
      size_t getGlobalMaxNumRowEntries() const;

      //! \brief Returns the maximum number of entries across all rows/columns on this node.
      /** Undefined if isFillActive().
        */
      size_t getNodeMaxNumRowEntries() const;

      //! \brief Indicates whether the matrix has a well-defined column map. 
      bool hasColMap() const; 

      //! \brief Indicates whether the matrix is lower triangular.
      /** Undefined if isFillActive().
        */
      bool isLowerTriangular() const;

      //! \brief Indicates whether the matrix is upper triangular.
      /** Undefined if isFillActive().
        */
      bool isUpperTriangular() const;

      //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
      bool isLocallyIndexed() const;

      //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
      bool isGloballyIndexed() const;

      //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
      bool isFillComplete() const;

      //! Returns \c true if resumeFill() has been called and the matrix is in edit mode.
      bool isFillActive() const;

      //! \brief Returns \c true if storage has been optimized.
      /**
        Optimized storage means that the allocation of each row is equal to the
        number of entries. The effect is that a pass through the matrix, i.e.,
        during a mat-vec, requires minimal memory traffic. One limitation of
        optimized storage is that no new indices can be added to the matrix.
        */
      bool isStorageOptimized() const;

      //! Returns \c true if the matrix was allocated with static data structures.
      ProfileType getProfileType() const;

      //! Indicates that the graph is static, so that new entries cannot be added to this matrix. */
      bool isStaticGraph() const;

      //! Extract a list of entries in a specified global row of this matrix. Put into pre-allocated storage.
      /*!
        \param LocalRow - (In) Global row number for which indices are desired.
        \param Indices - (Out) Global column indices corresponding to values.
        \param Values - (Out) Matrix values.
        \param NumEntries - (Out) Number of indices.

         Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
         with row \c GlobalRow. If \c GlobalRow does not belong to this node, then \c Indices and \c Values are unchanged and \c NumIndices is 
         returned as OrdinalTraits<size_t>::invalid().
       */
      void getGlobalRowCopy(GlobalOrdinal GlobalRow,
                            const ArrayView<GlobalOrdinal> &Indices,
                            const ArrayView<Scalar> &Values,
                            size_t &NumEntries
                            ) const;

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
      void getLocalRowCopy(LocalOrdinal LocalRow, 
                           const ArrayView<LocalOrdinal> &Indices, 
                           const ArrayView<Scalar> &Values,
                           size_t &NumEntries
                           ) const;

      //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
      /*!
        \param GlobalRow - (In) Global row number for which indices are desired.
        \param Indices   - (Out) Global column indices corresponding to values.
        \param Values    - (Out) Row values
        \pre <tt>isLocallyIndexed() == false</tt>
        \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

         Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
       */
      void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const Scalar> &values) const;

      //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
      /*!
        \param LocalRow - (In) Local row number for which indices are desired.
        \param Indices  - (Out) Global column indices corresponding to values.
        \param Values   - (Out) Row values
        \pre <tt>isGloballyIndexed() == false</tt>
        \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

         Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
       */
      void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices, ArrayView<const Scalar> &values) const;

      //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
      /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the 
          the zero and non-zero diagonals owned by this node. */
      void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const;

      //@}

      //! @name Advanced templated methods
      //@{

      //! Multiplies this matrix by a MultiVector.
      /*! \c X is required to be post-imported, i.e., described by the column map of the matrix. \c Y is required to be pre-exported, i.e., described by the row map of the matrix.

          Both are required to have constant stride, and they are not permitted to ocupy overlapping space. No runtime checking will be performed in a non-debug build.

          This method is templated on the scalar type of MultiVector objects, allowing this method to be applied to MultiVector objects of arbitrary type. However, it is recommended that multiply() not be called directly; instead, use the CrsMatrixMultiplyOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.
          
          If \c beta is equal to zero, the operation will enjoy overwrite semantics (\c Y will be overwritten with the result of the multiplication). Otherwise, the result of the multiplication
          will be accumulated into \c Y.
       */
      template <class DomainScalar, class RangeScalar>
      void multiply(const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, RangeScalar alpha, RangeScalar beta) const;

      //! Solves a linear system when the underlying matrix is triangular.
      /*! \c X is required to be post-imported, i.e., described by the column map of the matrix. \c Y is required to be pre-exported, i.e., described by the row map of the matrix.

          This method is templated on the scalar type of MultiVector objects, allowing this method to be applied to MultiVector objects of arbitrary type. However, it is recommended that solve() not be called directly; instead, use the CrsMatrixSolveOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.
          
          Both are required to have constant stride. However, unlike multiply(), it is permissible for <tt>&X == &Y</tt>. No runtime checking will be performed in a non-debug build.
       */
      template <class DomainScalar, class RangeScalar>
      void solve(const MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> & Y, MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> &X, Teuchos::ETransp trans) const;
          
      //@}

      //! Returns another CrsMatrix with the same entries, but represented as a different scalar type.
      template <class T>
      RCP<CrsMatrix<T,LocalOrdinal,GlobalOrdinal,Node> > convert() const;
          
      //@}

      //! @name Methods implementing Operator
      //@{ 

      //! \brief Computes the sparse matrix-multivector multiplication.
      /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
          - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
       */
      void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                 Teuchos::ETransp mode = Teuchos::NO_TRANS,
                 Scalar alpha = ScalarTraits<Scalar>::one(),
                 Scalar beta = ScalarTraits<Scalar>::zero()) const;

      //! Indicates whether this operator supports applying the adjoint operator.
      bool hasTransposeApply() const;

      //! \brief Returns the Map associated with the domain of this operator.
      //! This will be <tt>null</tt> until fillComplete() is called.
      const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const;

      //! Returns the Map associated with the domain of this operator.
      //! This will be <tt>null</tt> until fillComplete() is called.
      const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const;

      //@}

      //! @name Overridden from Teuchos::Describable 
      //@{

      /** \brief Return a simple one-line description of this object. */
      std::string description() const;

      /** \brief Print the object with some verbosity level to an FancyOStream object. */
      void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

      //@}

      //! @name Methods implementing Tpetra::DistObject
      //@{

      bool checkSizes(const DistObject<char, LocalOrdinal,GlobalOrdinal,Node>& source);

      void copyAndPermute(const DistObject<char, LocalOrdinal,GlobalOrdinal,Node>& source,
                          size_t numSameIDs,
                          const ArrayView<const LocalOrdinal> &permuteToLIDs,
                          const ArrayView<const LocalOrdinal> &permuteFromLIDs);

      void packAndPrepare(const DistObject<char, LocalOrdinal,GlobalOrdinal,Node>& source,
                          const ArrayView<const LocalOrdinal> &exportLIDs,
                          Array<char> &exports,
                          const ArrayView<size_t> & numPacketsPerLID,
                          size_t& constantNumPackets,
                          Distributor &distor);

      void unpackAndCombine(const ArrayView<const LocalOrdinal> &importLIDs,
                            const ArrayView<const char> &imports,
                            const ArrayView<size_t> &numPacketsPerLID,
                            size_t constantNumPackets,
                            Distributor &distor,
                            CombineMode CM);

      //@}

      //! \name Deprecated routines to be removed at some point in the future.
      //@{

      /** \brief Deprecated. Re-allocate the data into contiguous storage.

          This method is deprecated and will be removed in a future version of Tpetra, as 
          the implementation of storage optimization has been below Tpetra to Kokkos.

          Currently, the implementation simply calls resumeFill() and then fillComplete(OptimizeStorage). As such, it is 
          required to be called by all nodes that participate in the associated communicator.
       */
      TPETRA_DEPRECATED void optimizeStorage();

      //! Deprecated. Get a persisting const view of the entries in a specified global row of this matrix.
      TPETRA_DEPRECATED void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayRCP<const GlobalOrdinal> &indices, ArrayRCP<const Scalar> &values) const;

      //! Deprecated. Get a persisting const view of the entries in a specified local row of this matrix.
      TPETRA_DEPRECATED void getLocalRowView(LocalOrdinal LocalRow, ArrayRCP<const LocalOrdinal> &indices, ArrayRCP<const Scalar> &values) const;

      //@}

    private:
      // copy constructor disabled
      CrsMatrix(const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> &Source);
      // operator= disabled
      CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> & operator=(const CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> &rhs);
    protected:
      // useful typedefs
      typedef OrdinalTraits<LocalOrdinal>                     LOT;
      typedef OrdinalTraits<GlobalOrdinal>                    GOT;
      typedef ScalarTraits<Scalar>                             ST;
      typedef MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>       MV;
      typedef Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>             V;
      typedef CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>  Graph;
      // Enums
      enum GraphAllocationStatus {
        GraphAlreadyAllocated,
        GraphNotYetAllocated
      };
      // Allocation
      void allocateValues(ELocalGlobal lg, GraphAllocationStatus gas);
      // Sorting and merging
      void sortEntries();
      void mergeRedundantEntries();
      // global consts
      void clearGlobalConstants();
      void computeGlobalConstants();
      // matrix data accessors
      ArrayView<const Scalar>    getView(RowInfo rowinfo) const;
      ArrayView<      Scalar>    getViewNonConst(RowInfo rowinfo);
      // local Kokkos objects
      void pushToLocalMatrix();
      void pullFromLocalMatrix();
      void fillLocalMatrix(OptimizeOption os);
      void fillLocalSparseOps();
      // debugging
      void checkInternalState() const;

      // Two graph pointers needed in order to maintain const-correctness:
      // staticGraph_ is a graph passed to the constructor. We are not allowed to modify it. it is always a valid pointer.
      // myGraph_     is a graph created here. We are allowed to modify it. if myGraph_ != null, then staticGraph_ = myGraph_
      RCP<const Graph> staticGraph_;
      RCP<      Graph>     myGraph_;

      Kokkos::CrsMatrix<Scalar,LocalOrdinal,Node,LocalMatOps> lclMatrix_;
      typename LocalMatOps::template rebind<Scalar>::other    lclMatOps_;

      // matrix values. before allocation, both are null.
      // after allocation, one is null.
      // 1D == StaticAllocation, 2D == DynamicAllocation
      // The allocation always matches that of graph_, as the graph does the allocation for the matrix.
      ArrayRCP<Scalar>                       values1D_;
      ArrayRCP<ArrayRCP<Scalar> >            values2D_;
      // TODO: these could be allocated at resumeFill() and de-allocated at fillComplete() to make for very fast getView()/getViewNonConst()
      // ArrayRCP< typedef ArrayRCP<const Scalar>::iterator > rowPtrs_;
      // ArrayRCP< typedef ArrayRCP<      Scalar>::iterator > rowPtrsNC_;

      bool fillComplete_;

      // non-local data
      std::map<GlobalOrdinal, Array<std::pair<GlobalOrdinal,Scalar> > > nonlocals_;

      // a wrapper around multiply, for use in apply; it contains a non-owning RCP to *this, therefore, it is not allowed 
      // to persist past the destruction of *this. therefore, WE MAY NOT SHARE THIS POINTER.
      RCP< const CrsMatrixMultiplyOp<Scalar,Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > sameScalarMultiplyOp_;

  }; // class CrsMatrix

  /** \brief Non-member function to create an empty CrsMatrix given a row map and a non-zero profile.

      Returns a dynamically allocated (DynamicProfile) matrix with specified number of non-zeros per row (defaults to zero).

      \relatesalso CrsMatrix
   */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  createCrsMatrix(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, size_t maxNumEntriesPerRow = 0)
  {
    RCP<CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > ret;
    ret = rcp(new CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(map, maxNumEntriesPerRow, DynamicProfile) );
    return ret;
  }

} // namespace Tpetra

/**
  \example LocalMatOpExample.cpp
  An example using a different sparse mat-vec with Tpetra::CrsMatrix and Tpetra::CrsGraph.
 */

/** 
  \example CrsMatrix_NonlocalAfterResume.hpp
  An example for inserting non-local entries into a Tpetra::CrsMatrix using Tpetra::CrsMatrix::insertGlobalValues(), with multiple calls to Tpetra::CrsMatrix::fillComplete().
 */

#endif
