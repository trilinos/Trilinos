// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_REINDEX_MULTIVECTOR_DECL_HPP
#define TPETRA_REINDEX_MULTIVECTOR_DECL_HPP

/// \file Tpetra_Reindex_MultiVector_decl.hpp
/// \brief Declaration of the Tpetra::Reindex_MultiVector class

#include <Tpetra_Transform.hpp>
#include <Tpetra_MultiVector.hpp>

namespace Tpetra {

///
/** Given an input Tpetra MultiVector, a "reindexed" view is returned.
 */
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class Reindex_MultiVector : public ViewTransform<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > {
 public:
  using NewType      = typename ViewTransform<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::NewType;
  using OriginalType = typename ViewTransform<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::OriginalType;

  ///
  /** Constructor
   */
  Reindex_MultiVector(Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> const> newRowMap);

  ///
  /** Destructor
   */
  ~Reindex_MultiVector();

  ///
  /** Constructs a "reindexed" view of the original using the given NewRowMap.
   */
  NewType operator()(OriginalType const& origMultiVector);

 private:
  Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> const> newRowMap_;
};

}  // namespace Tpetra

#endif  // TPETRA_REINDEX_MULTIVECTOR_DECL_HPP
