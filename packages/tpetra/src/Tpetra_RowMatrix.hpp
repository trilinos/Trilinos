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

#ifndef TPETRA_ROWMATRIX_HPP
#define TPETRA_ROWMATRIX_HPP

#include <Teuchos_Describable.hpp>
#include <Kokkos_DefaultNode.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_RowGraph.hpp"
#include "Tpetra_Vector.hpp"

namespace Tpetra {

  //! \brief A pure virtual interface for row-partitioned matrices.
  /*!
     This class is templated on \c Scalar, \c LocalOrdinal, \c GlobalOrdinal and \c Node.
     The \c LocalOrdinal type, if omitted, defaults to \c int. 
     The \c GlobalOrdinal type defaults to the \c LocalOrdinal type.
     The \c Node type defaults to the default node in Kokkos.
   */
  template <class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class RowMatrix : virtual public Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
    public:
      //! @name Destructor Method
      //@{ 

      //! Destructor.
      virtual ~RowMatrix();

      //@}

      //! @name Matrix Query Methods
      //@{ 

      //! Returns the communicator.
      virtual const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const = 0;

      //! Returns the underlying node.
      virtual Teuchos::RCP<Node> getNode() const = 0;

      //! Returns the Map that describes the row distribution in this matrix.
      virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRowMap() const = 0;

      //! \brief Returns the Map that describes the column distribution in this matrix.
      virtual const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getColMap() const = 0;

      //! Returns the RowGraph associated with this matrix. 
      virtual Teuchos::RCP<const RowGraph<LocalOrdinal,GlobalOrdinal,Node> > getGraph() const = 0;

      //! Returns the number of global rows in this matrix.
      virtual global_size_t getGlobalNumRows() const = 0;

      //! \brief Returns the number of global columns in this matrix.
      virtual global_size_t getGlobalNumCols() const = 0;

      //! Returns the number of rows owned on the calling node.
      virtual size_t getNodeNumRows() const = 0;

      //! Returns the number of columns needed to apply the forward operator on this node, i.e., the number of elements listed in the column map.
      virtual size_t getNodeNumCols() const = 0;

      //! Returns the index base for global indices for this matrix. 
      virtual GlobalOrdinal getIndexBase() const = 0;

      //! Returns the global number of entries in this matrix.
      virtual global_size_t getGlobalNumEntries() const = 0;

      //! Returns the local number of entries in this matrix.
      virtual size_t getNodeNumEntries() const = 0;

      //! \brief Returns the current number of entries on this node in the specified global row.
      /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this graph. */
      virtual size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const = 0;

      //! Returns the current number of entries on this node in the specified local row.
      /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this graph. */
      virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const = 0;

      //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
      virtual global_size_t getGlobalNumDiags() const = 0;

      //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
      virtual size_t getNodeNumDiags() const = 0;

      //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
      virtual size_t getGlobalMaxNumRowEntries() const = 0;

      //! \brief Returns the maximum number of entries across all rows/columns on this node.
      virtual size_t getNodeMaxNumRowEntries() const = 0;

      //! \brief Indicates whether this matrix has a well-defined column map. 
      virtual bool hasColMap() const = 0;

      //! \brief Indicates whether this matrix is lower triangular.
      virtual bool isLowerTriangular() const = 0;

      //! \brief Indicates whether this matrix is upper triangular.
      virtual bool isUpperTriangular() const = 0;

      //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
      virtual bool isLocallyIndexed() const = 0;

      //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
      virtual bool isGloballyIndexed() const = 0;

      //! Returns \c true if fillComplete() has been called.
      virtual bool isFillComplete() const = 0;


      //@}

      //! @name Extraction Methods
      //@{ 

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
      virtual void getGlobalRowCopy(GlobalOrdinal GlobalRow,
                                    const Teuchos::ArrayView<GlobalOrdinal> &Indices,
                                    const Teuchos::ArrayView<Scalar> &Values,
                                    size_t &NumEntries) const = 0;

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
      virtual void getLocalRowCopy(LocalOrdinal LocalRow, 
                                   const Teuchos::ArrayView<LocalOrdinal> &Indices, 
                                   const Teuchos::ArrayView<Scalar> &Values,
                                   size_t &NumEntries) const  = 0;

      //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
      /*!
        \param GlobalRow - (In) Global row number for which indices are desired.
        \param Indices   - (Out) Global column indices corresponding to values.
        \param Values    - (Out) Row values
        \pre <tt>isLocallyIndexed() == false</tt>
        \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

         Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
       */
      virtual void getGlobalRowView(GlobalOrdinal GlobalRow, 
                                    ArrayView<const GlobalOrdinal> &indices, 
                                    ArrayView<const Scalar> &values) const = 0;

      //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
      /*!
        \param LocalRow - (In) Local row number for which indices are desired.
        \param Indices  - (Out) Global column indices corresponding to values.
        \param Values   - (Out) Row values
        \pre <tt>isGloballyIndexed() == false</tt>
        \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

         Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
       */
      virtual void getLocalRowView(LocalOrdinal LocalRow, 
                                   ArrayView<const LocalOrdinal> &indices, 
                                   ArrayView<const Scalar> &values) const = 0;

      //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
      /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the 
          the zero and non-zero diagonals owned by this node. */
      virtual void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const = 0;

      //@}
      
      //! \name Mathematical Methods
      //@{

      /**
       * \brief Scales the RowMatrix on the left with the Vector x.
       *
       * This matrix will be scaled such that A(i,j) = x(i)*A(i,j) 
       * where i denoes the global row number of A and 
       * j denotes the global column number of A.
       *
       * \param x A vector to left scale this matrix.
       */
      virtual void leftScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) =0;

      /**
       * \brief Scales the RowMatrix on the right with the Vector x.
       *
       * This matrix will be scaled such that A(i,j) = x(j)*A(i,j) 
       * where i denoes the global row number of A and 
       * j denotes the global column number of A.
       *
       * \param x A vector to right scale this matrix.
       */
      virtual void rightScale(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x) =0;

      /** \brief Returns the Frobenius norm of the matrix */
      virtual Scalar getFrobeniusNorm() const = 0;

      //@}

      //! \name Deprecated routines to be removed at some point in the future.
      //@{

      //! Deprecated. Get a persisting const view of the entries in a specified global row of this matrix.
      /*!
        \param GlobalRow - (In) Global row from which to retrieve matrix entries.
        \param Indices - (Out) Indices for the global row.
        \param Values - (Out) Values for the global row.

         Note: If \c GlobalRow does not belong to this node, then \c Indices and \c Values are set to <tt>Teuchos::null</t>>.

        \pre isLocallyIndexed()==false
       */
      TPETRA_DEPRECATED virtual void getGlobalRowView(GlobalOrdinal GlobalRow, 
                                                      Teuchos::ArrayRCP<const GlobalOrdinal> &indices,
                                                      Teuchos::ArrayRCP<const Scalar>        &values) const = 0;

      //! Deprecated. Get a persisting const view of the entries in a specified local row of this matrix.
      /*!
        \param LocalRow - (In) Local row from which to retrieve matrix entries.
        \param Indices - (Out) Indices for the local row.
        \param Values - (Out) Values for the local row.

         Note: If \c LocalRow is not valid for this node, then \c Indices and \c Values are set to <tt>Teuchos::null</tt>.

        \pre isGloballyIndexed()==false
       */
      TPETRA_DEPRECATED virtual void getLocalRowView(LocalOrdinal LocalRow,
                                                    Teuchos::ArrayRCP<const LocalOrdinal> &indices,
                                                    Teuchos::ArrayRCP<const Scalar>       &values) const = 0;

      //@}

  }; // class RowMatrix

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~RowMatrix() {
  }

} // namespace Tpetra

#endif
