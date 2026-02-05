// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_TPETRAVECTOR_DEF_HPP
#define XPETRA_TPETRAVECTOR_DEF_HPP
#include "Xpetra_TpetraVector_decl.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraVector(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
                 bool zeroOut)
  : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, 1, zeroOut) {
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraVector(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node>>& map,
                 const Teuchos::ArrayView<const Scalar>& A)
  : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, A, map->getLocalNumElements(), 1) {
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ~TpetraVector() {
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceGlobalValue");
  getTpetra_Vector()->replaceGlobalValue(globalRow, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoGlobalValue");
  getTpetra_Vector()->sumIntoGlobalValue(globalRow, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceLocalValue(LocalOrdinal myRow, const Scalar& value) {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::replaceLocalValue");
  getTpetra_Vector()->replaceLocalValue(myRow, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoLocalValue(LocalOrdinal myRow, const Scalar& value) {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::sumIntoLocalValue");
  getTpetra_Vector()->sumIntoLocalValue(myRow, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm1() const {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm1");
  return getTpetra_Vector()->norm1();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm2() const {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::norm2");
  return getTpetra_Vector()->norm2();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    normInf() const {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::normInf");
  return getTpetra_Vector()->normInf();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    meanValue() const {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::meanValue");
  return getTpetra_Vector()->meanValue();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    description() const {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::description");
  return getTpetra_Vector()->description();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::describe");
  getTpetra_Vector()->describe(out, verbLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    dot(const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& a) const {
  XPETRA_MONITOR("TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dot");
  return getTpetra_Vector()->dot(*toTpetra(a));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraVector(const Teuchos::RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& vec)
  : TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(vec) {
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getTpetra_Vector() const {
  return this->TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getTpetra_MultiVector()->getVectorNonConst(0);
}

}  // namespace Xpetra

#endif
