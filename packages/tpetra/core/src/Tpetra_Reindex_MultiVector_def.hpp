// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_REINDEX_MULTIVECTOR_DEF_HPP
#define TPETRA_REINDEX_MULTIVECTOR_DEF_HPP

#include <Tpetra_Reindex_MultiVector_decl.hpp>

/// \file Tpetra_Reindex_MultiVector_def.hpp
/// \brief Definition of the Tpetra::Reindex_MultiVector class
///
/// If you want to use Tpetra::Reindex_MultiVector, include
/// "Tpetra_Reindex_MultiVector.hpp", a file which CMake generates
/// and installs for you.
///
/// If you only want the declaration of Tpetra::Reindex_MultiVector,
/// include "Tpetra_Reindex_MultiVector_decl.hpp".

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Reindex_MultiVector(Teuchos::RCP<Map<LocalOrdinal, GlobalOrdinal, Node> const> newRowMap)
  : ViewTransform<MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >()
  , newRowMap_(newRowMap) {
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~Reindex_MultiVector() {
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::operator()(OriginalType const& origMultiVector) {
  this->origObj_ = origMultiVector;
  assert(origMultiVector->getMap()->getLocalNumElements() == this->newRowMap_->getLocalNumElements());
  assert(origMultiVector->isConstantStride() == true);  // So that it is valid to call origMultiVector->getStride()

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > origValues = origMultiVector->get2dView();

  using size_type = typename Teuchos::ArrayRCP<Scalar>::size_type;

  size_type numEntries(origValues.size() * origValues[0].size());

  std::vector<Scalar> tmpVec(numEntries);
  size_t k(0);
  for (size_type v(0); v < origValues.size(); ++v) {
    for (size_type i(0); i < origValues[v].size(); ++i) {
      tmpVec[k++] = origValues[v][i];
    }
  }
  Teuchos::ArrayView<Scalar const> valuesToInsert(tmpVec.data(), numEntries);

  using mv_t    = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  this->newObj_ = Teuchos::RCP<mv_t>(new mv_t(newRowMap_  // const Teuchos::RCP<const map_type> & map
                                              ,
                                              valuesToInsert  // const Teuchos::ArrayView<const Scalar> & A
                                              ,
                                              origMultiVector->getStride()  // const size_t LDA
                                              ,
                                              origMultiVector->getNumVectors()  // const size_t NumVectors
                                              ));

  return this->newObj_;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_REINDEXMULTIVECTOR_INSTANT(SCALAR, LO, GO, NODE) \
  template class Reindex_MultiVector<SCALAR, LO, GO, NODE>;

}  // namespace Tpetra

#endif  // TPETRA_REINDEX_MULTIVECTOR_DEF_HPP
