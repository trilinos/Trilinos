// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_MATRIXVIEW_HPP
#define XPETRA_MATRIXVIEW_HPP

#include <Teuchos_Describable.hpp>
#include <Tpetra_KokkosCompat_DefaultNode.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Map.hpp"

/** \file Xpetra_MatrixView.hpp

Declarations for the class Xpetra::MatrixView.
*/
namespace Xpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class MatrixView {  // TODO : virtual public Teuchos::Describable {
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;

 public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor
  MatrixView(const RCP<const Map> &rowMap, const RCP<const Map> &colMap)
    : rowMap_(rowMap)
    , colMap_(colMap)
    , maxEigValueEstimate_(-Teuchos::ScalarTraits<Scalar>::one()) {}

  //! Destructor
  virtual ~MatrixView() {}

  //@}

  //! @name Map access methods
  //@{
  //! Returns the Map that describes the row distribution in this matrix.
  const RCP<const Map> &GetRowMap() const { return rowMap_; }

  //! \brief Returns the Map that describes the column distribution in this matrix.
  const RCP<const Map> &GetColMap() const { return colMap_; }

  //! Returns the Map that describes the row distribution in this matrix.
  void SetRowMap(const RCP<const Map> &rowMap) { rowMap_ = rowMap; }

  //! \brief Set the Map that describes the column distribution in this matrix.
  void SetColMap(const RCP<const Map> &colMap) { colMap_ = colMap; }
  //@}

  //! \brief Set an maximum eigenvalue estimate for this matrix.
  void SetMaxEigenvalueEstimate(Scalar const &sigma) { maxEigValueEstimate_ = sigma; }

  //! \brief Return the maximum eigenvalue estimate for this matrix.
  Scalar GetMaxEigenvalueEstimate() const { return maxEigValueEstimate_; }

 private:
  RCP<const Map> rowMap_;
  RCP<const Map> colMap_;

  Scalar maxEigValueEstimate_;

};  // class MatrixView

}  // namespace Xpetra

#define XPETRA_MATRIXVIEW_SHORT
#endif  // XPETRA_MATRIX_VIEW_DECL_HPP
