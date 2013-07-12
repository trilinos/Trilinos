/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_DETAILS_TPETRA_ROWGRAPH_DECL_HPP
#define IFPACK2_DETAILS_TPETRA_ROWGRAPH_DECL_HPP

/// \file Ifpack2_Details_Tpetra_RowGraph_decl.hpp
/// \brief Declaration of a minimalist Tpetra::RowGraph for use in AdditiveSchwarz
/// \author Tom Benson
///
/// This file is meant for Ifpack2 developers only, not for users.
/// It declares a Tpetra::RowGraph to be used with LocalFilter.

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_RowGraph.hpp"
#include "Tpetra_RowMatrix.hpp"

namespace Ifpack2 {

namespace Details {
/// \class Tpetra_RowGraph
/// \brief Minimal necessary functionality of Tpetra::RowGraph to satisfy a partitioner.
/// \tparam MatrixType A specialization of either Tpetra::RowGraph,
///   or Tpetra::CrsGraph (which is a subclass of Tpetra::RowGraph).
///
template<class MatrixType>
class Tpetra_RowGraph :
    virtual public Tpetra::RowGraph<typename MatrixType::local_ordinal_type,
                                    typename MatrixType::global_ordinal_type,
                                    typename MatrixType::node_type>
{
public:

  //! \name Typedefs
  //@{

  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type scalar_type;

  //! Preserved only for backwards compatibility.  Please use "scalar_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::scalar_type Scalar;


  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! Preserved only for backwards compatibility.  Please use "local_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::local_ordinal_type LocalOrdinal;


  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! Preserved only for backwards compatibility.  Please use "global_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::global_ordinal_type GlobalOrdinal;


  //! The type of the Kokkos Node used by the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! Preserved only for backwards compatibility.  Please use "node_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::node_type Node;


  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! Preserved only for backwards compatibility.  Please use "magnitude_type".
  TEUCHOS_DEPRECATED typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitudeType;

  //@}
  //! \name Constructor and destructor
  //@{

  /// \brief Constructor
  ///
  /// \param A [in] The sparse matrix to which to apply the local
  ///   filter, as a Tpetra::RowMatrix.  (Tpetra::CrsMatrix inherits
  ///   from this, so you may use a Tpetra::CrsMatrix here instead.)
  ///
  /// This class will <i>not</i> modify the input matrix.
  explicit Tpetra_RowGraph (const Teuchos::RCP<const Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> >& Matrix);

  //! Destructor
  virtual ~Tpetra_RowGraph();

  //@}
  //! \name Matrix Query Methods
  //@{

  //! Returns the communicator.
  virtual const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const;

  //! Returns the underlying Kokkos Node object.
  virtual Teuchos::RCP<node_type> getNode() const;

  //! Returns the Map that describes the row distribution in this matrix.
  virtual const Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > & getRowMap() const;

  //! Returns the Map that describes the column distribution in this matrix.
  virtual const Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > & getColMap() const;

  //! Returns the Map that describes the domain distribution in this matrix.
  virtual const Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > & getDomainMap() const;

  //! Returns the Map that describes the range distribution in this matrix.
  virtual const Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > & getRangeMap() const;

  //! The number of global rows in this matrix.
  virtual global_size_t getGlobalNumRows() const;

  //! The number of global columns in this matrix.
  virtual global_size_t getGlobalNumCols() const;

  //! The number of rows owned on the calling process.
  virtual size_t getNodeNumRows() const;

  //! The number of columns in the (locally filtered) matrix.
  virtual size_t getNodeNumCols() const;

  //! Returns the index base for global indices for this matrix.
  virtual global_ordinal_type getIndexBase() const;

  //! Returns the global number of entries in this matrix.
  virtual global_size_t getGlobalNumEntries() const;

  //! Returns the local number of entries in this matrix.
  virtual size_t getNodeNumEntries() const;

  //! Does nothing.
  virtual Teuchos::RCP<const Tpetra::Import<typename MatrixType::local_ordinal_type,global_ordinal_type,node_type> > getImporter() const;

  //! Does nothing.
  virtual Teuchos::RCP<const Tpetra::Export<local_ordinal_type,global_ordinal_type,node_type> > getExporter() const;


  /// \brief The current number of entries on this node in the specified global row.
  ///
  /// \return <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt> if
  ///   the specified row is not owned by this process, otherwise the
  ///   number of entries in that row on this process.
  virtual size_t getNumEntriesInGlobalRow (global_ordinal_type globalRow) const;

  /// \brief The current number of entries on this node in the specified local row.
  ///
  /// \return <tt>Teuchos::OrdinalTraits<size_t>::invalid()</tt> if
  ///   the specified local row is not valid on this process,
  ///   otherwise the number of entries in that row on this process.
  virtual size_t getNumEntriesInLocalRow (local_ordinal_type localRow) const;

  //! The number of global diagonal entries, based on global row/column index comparisons.
  virtual global_size_t getGlobalNumDiags() const;

  //! The number of local diagonal entries, based on global row/column index comparisons.
  virtual size_t getNodeNumDiags() const;

  //! The maximum number of entries across all rows/columns on all processes.
  virtual size_t getGlobalMaxNumRowEntries() const;

  //! The maximum number of entries across all rows/columns on this process.
  virtual size_t getNodeMaxNumRowEntries() const;

  //! Whether this matrix has a well-defined column Map.
  virtual bool hasColMap() const;

  //! Whether this matrix is lower triangular.
  virtual bool isLowerTriangular() const;

  //! Whether this matrix is upper triangular.
  virtual bool isUpperTriangular() const;

  //! Whether the underlying sparse matrix is locally (opposite of globally) indexed.
  virtual bool isLocallyIndexed() const;

  //! Whether the underlying sparse matrix is globally (opposite of locally) indexed.
  virtual bool isGloballyIndexed() const;

  //! Returns \c true if fillComplete() has been called.
  virtual bool isFillComplete() const;

  

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
  virtual void
  getGlobalRowCopy (global_ordinal_type GlobalRow,
                    const Teuchos::ArrayView<global_ordinal_type> &Indices,
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
  virtual void
  getLocalRowCopy (local_ordinal_type LocalRow,
                   const Teuchos::ArrayView<local_ordinal_type> &Indices,
                   size_t &NumEntries) const ;


  //@}


private:

  //! Pointer to the matrix to be preconditioned.
  Teuchos::RCP<const Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> > A_;


};// class Tpetra_RowGraph
}// namespace Details
}// namespace Ifpack2

#endif /* IFPACK2_TPETRA_ROWGRAPH_DECL_HPP */
