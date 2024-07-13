// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_OPERATORWITHDIAGONAL_HPP
#define XPETRA_OPERATORWITHDIAGONAL_HPP

#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_TpetraRowMatrix.hpp>

namespace Xpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class OperatorWithDiagonal : virtual public Xpetra::TpetraRowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  typedef Scalar scalar_type;
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;

  //! @name Constructor/Destructor Methods
  //@{

  //! Destructor.
  virtual ~OperatorWithDiagonal() = default;

  //@}

  //! @name Matrix Query Methods
  //@{

  //! Returns the Map that describes the row distribution in this matrix.
  const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getRowMap() const {
    throw std::runtime_error("Not implemented.");
  };

  //! Returns the Map that describes the column distribution in this matrix.
  const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > getColMap() const {
    throw std::runtime_error("Not implemented.");
  };

  //! Returns the number of global rows in this matrix.
  global_size_t getGlobalNumRows() const {
    return this->getRangeMap()->getGlobalNumElements();
  };

  //! Returns the number of global columns in this matrix.
  global_size_t getGlobalNumCols() const {
    return this->getDomainMap()->getGlobalNumElements();
  };

  //! Returns the number of rows owned on the calling node.
  size_t getLocalNumRows() const {
    return this->getRangeMap()->getLocalNumElements();
  };

  //! Returns the number of columns needed to apply the forward operator on this node, i.e., the number of elements listed in the column map.
  size_t getLocalNumCols() const {
    return this->getDomainMap()->getLocalNumElements();
  };

  //! Returns the global number of entries in this matrix.
  global_size_t getGlobalNumEntries() const {
    throw std::runtime_error("Not implemented.");
  };

  //! Returns the local number of entries in this matrix.
  size_t getLocalNumEntries() const {
    throw std::runtime_error("Not implemented.");
  };

  //! Returns the current number of entries on this node in the specified local row.
  size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const {
    throw std::runtime_error("Not implemented.");
  };

  //! Returns the maximum number of entries across all rows/columns on all nodes.
  size_t getGlobalMaxNumRowEntries() const {
    throw std::runtime_error("Not implemented.");
  };

  //! Returns the maximum number of entries across all rows/columns on this node.
  size_t getLocalMaxNumRowEntries() const {
    throw std::runtime_error("Not implemented.");
  };

  //! If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */.
  bool isLocallyIndexed() const {
    return true;
  };

  //! If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */.
  bool isGloballyIndexed() const {
    return true;
  };

  //! Returns true if fillComplete() has been called.
  bool isFillComplete() const {
    return true;
  };

  //! Returns true if getLocalRowView() and getGlobalRowView() are valid for this class.
  bool supportsRowViews() const {
    return false;
  };

  //@}

  //! @name Extraction Methods
  //@{

  //! Extract a list of entries in a specified local row of the graph. Put into storage allocated by calling routine.
  void getLocalRowCopy(LocalOrdinal LocalRow, const Teuchos::ArrayView<LocalOrdinal> &Indices, const Teuchos::ArrayView<Scalar> &Values, size_t &NumEntries) const {
    throw std::runtime_error("Not implemented.");
  };

  //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
  void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const Scalar> &values) const {
    throw std::runtime_error("Not implemented.");
  };

  //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
  void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices, ArrayView<const Scalar> &values) const {
    throw std::runtime_error("Not implemented.");
  };

  //! Get a copy of the diagonal entries owned by this node, with local row indices.
  void getLocalDiagCopy(Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &diag) const {
    throw std::runtime_error("Not implemented.");
  };

  //@}

  //! @name Mathematical Methods
  //@{

  //! Returns the Frobenius norm of the matrix.
  typename ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const {
    throw std::runtime_error("Not implemented.");
  };

  //@}

};  // RowMatrix class

}  // namespace Xpetra

#endif
