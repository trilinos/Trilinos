// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_MULTIVECTORFACTORY_DECL_HPP
#define XPETRA_MULTIVECTORFACTORY_DECL_HPP

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_MultiVector_decl.hpp"

#include "Xpetra_TpetraMultiVector.hpp"

// #include "Xpetra_BlockedMap.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

/*!
  @class MultiVectorFactory
  @brief Factory for any type of Xpetra::MultiVector and its derived classes

  Creates instances of \c Xpetra::MulitVector and \c Xpetra::BlockedMultiVector ,
  depending on the type of map, i.e. \c Xpetra::Map vs. \c Xpetra::BlockedMap .

  @tparam Scalar
  @tparam LocalOrdinal
  @tparam GlobalOrdinal
  @tparam Node

  \note Although this class can gerenate \c Xpetra::BlockedMultiVector ,
  it always returns \c Xpetra::MultiVector .
  Don't forget to cast to \c Xpetra::BlockedMultiVector , if you need the blocked layout directly.
*/
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class MultiVectorFactory {
 private:
  //! Private constructor. This is a static class.
  MultiVectorFactory() {}

 public:
  //! Constructor specifying the number of non-zeros for all rows.
  static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &map, size_t NumVectors, bool zeroOut = true);

  //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy).
  static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Build(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>> &map,
        const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar>> &ArrayOfPtrs,
        size_t NumVectors);

  static Teuchos::RCP<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
  Build(const Teuchos::RCP<const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> &source,
        Teuchos::DataAccess copyOrView);
};

}  // namespace Xpetra

#define XPETRA_MULTIVECTORFACTORY_SHORT
#endif
