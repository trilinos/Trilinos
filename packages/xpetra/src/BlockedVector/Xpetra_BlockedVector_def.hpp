// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_BLOCKEDVECTOR_DEF_HPP
#define XPETRA_BLOCKEDVECTOR_DEF_HPP

#include "Xpetra_BlockedVector_decl.hpp"

#include "Xpetra_BlockedMultiVector.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedVector(const Teuchos::RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>& map, bool zeroOut)
  : Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, 1, zeroOut) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedVector(Teuchos::RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap,
                  Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> v)
  : Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bmap, v) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedVector(Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mapExtractor,
                  Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> v)
  : Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(mapExtractor, v) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ~BlockedVector() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator=(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& rhs) {
  assign(rhs);  // dispatch to protected virtual method
  return *this;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar& value) {
  BlockedMultiVector::replaceGlobalValue(globalRow, vectorIndex, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar& value) {
  BlockedMultiVector::sumIntoGlobalValue(globalRow, vectorIndex, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar& value) {
  BlockedMultiVector::replaceLocalValue(myRow, vectorIndex, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar& value) {
  BlockedMultiVector::sumIntoLocalValue(myRow, vectorIndex, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {
  BlockedMultiVector::replaceGlobalValue(globalRow, 0, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {
  BlockedMultiVector::sumIntoGlobalValue(globalRow, 0, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceLocalValue(LocalOrdinal myRow, const Scalar& value) {
  BlockedMultiVector::replaceLocalValue(myRow, 0, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoLocalValue(LocalOrdinal myRow, const Scalar& value) {
  BlockedMultiVector::sumIntoLocalValue(myRow, 0, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    putScalar(const Scalar& value) {
  BlockedMultiVector::putScalar(value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVector(size_t j) const {
  return BlockedMultiVector::getVector(j);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVectorNonConst(size_t j) {
  return BlockedMultiVector::getVectorNonConst(j);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<const Scalar>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getData(size_t j) const {
  return BlockedMultiVector::getData(j);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<Scalar>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getDataNonConst(size_t j) {
  return BlockedMultiVector::getDataNonConst(j);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    dot(const MultiVector& A, const Teuchos::ArrayView<Scalar>& dots) const {
  BlockedMultiVector::dot(A, dots);
  return;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    dot(const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A) const {
  Teuchos::Array<Scalar> dots = Teuchos::Array<Scalar>(1);
  BlockedMultiVector::dot(A, dots);
  return dots[0];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    abs(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A) {
  BlockedMultiVector::abs(A);
  return;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    reciprocal(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A) {
  BlockedMultiVector::reciprocal(A);
  return;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    scale(const Scalar& alpha) {
  BlockedMultiVector::scale(alpha);
  return;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    scale(Teuchos::ArrayView<const Scalar> alpha) {
  BlockedMultiVector::scale(alpha);
  return;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    update(const Scalar& alpha,
           const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
           const Scalar& beta) {
  BlockedMultiVector::update(alpha, A, beta);
  return;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    update(const Scalar& alpha,
           const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
           const Scalar& beta,
           const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
           const Scalar& gamma) {
  BlockedMultiVector::update(alpha, A, beta, B, gamma);
  return;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm1() const {
  using Array = Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>;
  Array norm  = Array(1);
  this->norm1(norm);
  return norm[0];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm2() const {
  Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> norm =
      Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(1);
  this->norm2(norm);
  return norm[0];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    normInf() const {
  Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>
      norm = Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>(1);
  this->normInf(norm);
  return norm[0];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm1(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const {
  BlockedMultiVector::norm1(norms);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const {
  BlockedMultiVector::norm2(norms);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const {
  BlockedMultiVector::normInf(norms);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    meanValue(const Teuchos::ArrayView<Scalar>& /* means */) const {
  throw Xpetra::Exceptions::RuntimeError("BlockedVector::meanValue: Not (yet) supported by BlockedVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    meanValue() const {
  throw Xpetra::Exceptions::RuntimeError("BlockedVector::meanValue: Not (yet) supported by BlockedVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    multiply(Teuchos::ETransp /* transA */,
             Teuchos::ETransp /* transB */,
             const Scalar& /* alpha */,
             const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* A */,
             const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* B */,
             const Scalar& /* beta */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedVector::multiply: Not (yet) supported by BlockedVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    multiply(Teuchos::ETransp /* transA */,
             Teuchos::ETransp /* transB */,
             const Scalar& /* alpha */,
             const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* A */,
             const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* B */,
             const Scalar& /* beta */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedVector::multiply: Not (yet) supported by BlockedVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    elementWiseMultiply(Scalar /* scalarAB */,
                        const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* A */,
                        const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* B */,
                        Scalar /* scalarThis */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedVector::elementWiseMultiply: Not (yet) supported by BlockedVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    elementWiseMultiply(Scalar /* scalarAB */,
                        const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                        const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
                        Scalar /* scalarThis */) {
  XPETRA_TEST_FOR_EXCEPTION(B.getMap()->isSameAs(*(this->getMap())) == false,
                            Xpetra::Exceptions::RuntimeError,
                            "BlockedVector::elementWiseMultipy: B must have same blocked map than this.");
  TEUCHOS_TEST_FOR_EXCEPTION(A.getMap()->getLocalNumElements() != B.getMap()->getLocalNumElements(),
                             Xpetra::Exceptions::RuntimeError,
                             "BlockedVector::elementWiseMultipy: A has "
                                 << A.getMap()->getLocalNumElements() << " elements, B has " << B.getMap()->getLocalNumElements()
                                 << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(A.getMap()->getGlobalNumElements() != B.getMap()->getGlobalNumElements(),
                             Xpetra::Exceptions::RuntimeError,
                             "BlockedVector::elementWiseMultipy: A has " << A.getMap()->getGlobalNumElements()
                                                                         << " elements, B has "
                                                                         << B.getMap()->getGlobalNumElements() << ".");

  RCP<const BlockedMap> bmap                                                 = this->getBlockedMap();
  RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> rcpA  = Teuchos::rcpFromRef(A);
  RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bmvec = Teuchos::rcpFromRef(B);
  RCP<const BlockedVector> bbmvec                                            = Teuchos::rcp_dynamic_cast<const BlockedVector>(bmvec);
  TEUCHOS_TEST_FOR_EXCEPTION(bbmvec.is_null() == true,
                             Xpetra::Exceptions::RuntimeError,
                             "BlockedVector::elementWiseMultipy: B must be a BlockedVector.");

  // TODO implement me
  /*RCP<Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node> > me = Teuchos::rcp(new
  Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node>(bmap));

  for(size_t m = 0; m < bmap->getNumMaps(); m++) {
    // TODO introduce BlockedVector objects and "skip" this expensive ExtractVector call
    RCP<const Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > pd = me->ExtractVector(rcpA,m,bmap->getThyraMode());
    XPETRA_TEST_FOR_EXCEPTION(pd->getMap()->isSameAs(*(this->getBlockedMap()->getMap(m,bmap->getThyraMode())))==false,
  Xpetra::Exceptions::RuntimeError, "BlockedVector::elementWiseMultipy: sub map of B does not fit with sub map of this.");
    this->getMultiVector(m,bmap->getThyraMode())->elementWiseMultiply(scalarAB,*pd,*(bbmvec->getMultiVector(m,bmap->getThyraMode())),scalarThis);
  }*/
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getNumVectors() const {
  return 1;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalLength() const {
  throw Xpetra::Exceptions::RuntimeError(
      "BlockedVector::getLocalLength: routine not implemented. It has no value as one must iterate on the partial vectors.");
  TEUCHOS_UNREACHABLE_RETURN(0);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalLength() const {
  return this->getBlockedMap()->getFullMap()->getGlobalNumElements();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    isSameSize(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* vec */) const {
  throw Xpetra::Exceptions::RuntimeError(
      "BlockedVector::isSameSize: routine not implemented. It has no value as one must iterate on the partial vectors.");
  TEUCHOS_UNREACHABLE_RETURN(0);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    description() const {
  return std::string("BlockedVector");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
  out << description() << std::endl;
  for (size_t r = 0; r < this->getBlockedMap()->getNumMaps(); r++) {
    getMultiVector(r)->describe(out, verbLevel);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceMap(const RCP<const Map>& map) {
  BlockedMultiVector::replaceMap(map);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doImport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* source */,
             const Import& /* importer */,
             CombineMode /* CM */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedVector::doImport: Not supported by BlockedVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doExport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* dest */,
             const Import& /* importer */,
             CombineMode /* CM */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedVector::doExport: Not supported by BlockedVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doImport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* source */,
             const Export& /* exporter */,
             CombineMode /* CM */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedVector::doImport: Not supported by BlockedVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doExport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* dest */,
             const Export& /* exporter */,
             CombineMode /* CM */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedVector::doExport: Not supported by BlockedVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setSeed(unsigned int seed) {
  for (size_t r = 0; r < this->getBlockedMap()->getNumMaps(); ++r) {
    getMultiVector(r)->setSeed(seed);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    randomize(bool bUseXpetraImplementation) {
  for (size_t r = 0; r < this->getBlockedMap()->getNumMaps(); ++r) {
    getMultiVector(r)->randomize(bUseXpetraImplementation);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    randomize(const Scalar& minVal, const Scalar& maxVal, bool bUseXpetraImplementation) {
  for (size_t r = 0; r < this->getBlockedMap()->getNumMaps(); ++r) {
    getMultiVector(r)->randomize(minVal, maxVal, bUseXpetraImplementation);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Xpetra_randomize() {
  {
    Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Xpetra_randomize();
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Xpetra_randomize(const Scalar& minVal, const Scalar& maxVal) {
  {
    Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Xpetra_randomize(minVal, maxVal);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMap() const {
  XPETRA_MONITOR("BlockedVector::getMap");
  return this->getBlockedMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMultiVector(size_t r) const {
  return BlockedMultiVector::getMultiVector(r);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMultiVector(size_t r, bool bThyraMode) const {
  return BlockedMultiVector::getMultiVector(r, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setMultiVector(size_t r,
                   Teuchos::RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> v,
                   bool bThyraMode) {
  BlockedMultiVector::setMultiVector(r, v, bThyraMode);
  return;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Merge() const {
  return BlockedMultiVector::Merge();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    assign(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& rhs) {
  BlockedMultiVector::assign(rhs);
}

// template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
// virtual void BlockedVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
// assign (const XpetrA::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& rhs)
// {
//     throw Xpetra::Exceptions::RuntimeError("BlockedVector::assign: Not supported by BlockedVector.");
// }

}  // namespace Xpetra

#endif  // XPETRA_BLOCKEDVECTOR_DEF_HPP
