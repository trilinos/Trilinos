// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_MULTIVECTOR_DEF_HPP
#define XPETRA_MULTIVECTOR_DEF_HPP

#include "Xpetra_MultiVector_decl.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ~MultiVector() {
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator=(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& rhs) {
  assign(rhs);  // dispatch to protected virtual method
  return *this;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Xpetra_randomize() {
  typedef Teuchos::ScalarTraits<Scalar> SCT;

  const size_t numVectors = getNumVectors();
  for (size_t i = 0; i < numVectors; i++) {
    Teuchos::ArrayRCP<Scalar> datai = getDataNonConst(i);

    const size_t myLength = getLocalLength();
    for (size_t j = 0; j < myLength; j++) {
      datai[j] = SCT::random();
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Xpetra_randomize(const Scalar& minVal, const Scalar& maxVal) {
  typedef Teuchos::ScalarTraits<Scalar> SCT;
  Scalar point5 = SCT::one() / (SCT::one() + SCT::one());

  const size_t numVectors = getNumVectors();
  for (size_t i = 0; i < numVectors; i++) {
    Teuchos::ArrayRCP<Scalar> datai = getDataNonConst(i);

    const size_t myLength = getLocalLength();
    for (size_t j = 0; j < myLength; j++) {
      datai[j] = point5 * (maxVal - minVal) * SCT::random() + point5 * (maxVal + minVal);
    }
  }
}

}  // namespace Xpetra

#endif  // XPETRA_MULTIVECTOR_DEF_HPP
