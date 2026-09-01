// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKEDVECTOR_DEF_HPP
#define TPETRA_BLOCKEDVECTOR_DEF_HPP

#include "Tpetra_BlockedVector_decl.hpp"

#include "Tpetra_BlockedMultiVector_def.hpp"
#include "Tpetra_Vector.hpp"

#include <Teuchos_ScalarTraits.hpp>

namespace Tpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedVector(const Teuchos::RCP<const BlockedMap>& map, bool zeroOut)
  : base_type() {
  bmv_ = Teuchos::rcp(new BlockedMultiVector(map, 1, zeroOut));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedVector(Teuchos::RCP<const BlockedMap> bmap,
                  Teuchos::RCP<Vector> v)
  : base_type() {
  bmv_ = Teuchos::rcp(new BlockedMultiVector(bmap, v));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedVector(Teuchos::RCP<const MapExtractor> mapExtractor,
                  Teuchos::RCP<Vector> v)
  : base_type() {
  bmv_ = Teuchos::rcp(new BlockedMultiVector(mapExtractor, v));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ~BlockedVector() {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator=(const MultiVector& rhs) {
  assign(rhs);  // dispatch to protected virtual method
  return *this;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar& value) {
  bmv_->replaceGlobalValue(globalRow, vectorIndex, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar& value) {
  bmv_->sumIntoGlobalValue(globalRow, vectorIndex, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar& value) {
  bmv_->replaceLocalValue(myRow, vectorIndex, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar& value) {
  bmv_->sumIntoLocalValue(myRow, vectorIndex, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {
  bmv_->replaceGlobalValue(globalRow, 0, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar& value) {
  bmv_->sumIntoGlobalValue(globalRow, 0, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceLocalValue(LocalOrdinal myRow, const Scalar& value) {
  bmv_->replaceLocalValue(myRow, 0, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoLocalValue(LocalOrdinal myRow, const Scalar& value) {
  bmv_->sumIntoLocalValue(myRow, 0, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    putScalar(const Scalar& value) {
  bmv_->putScalar(value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Vector>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVector(size_t j) const {
  return bmv_->getVector(j);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Vector>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVectorNonConst(size_t j) {
  return bmv_->getVectorNonConst(j);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<const Scalar>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getData(size_t j) const {
  return bmv_->getData(j);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<Scalar>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getDataNonConst(size_t j) {
  return bmv_->getDataNonConst(j);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    dot(const MultiVector& A, const Teuchos::ArrayView<dot_type>& dots) const {
  bmv_->dot(A, dots);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dot_type
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    dot(const Vector& A) const {
  Teuchos::Array<dot_type> dots = Teuchos::Array<dot_type>(1);
  bmv_->dot(A, dots);
  return dots[0];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    abs(const MultiVector& A) {
  bmv_->abs(A);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    reciprocal(const MultiVector& A) {
  bmv_->reciprocal(A);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    scale(const Scalar& alpha) {
  bmv_->scale(alpha);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    scale(Teuchos::ArrayView<const Scalar> alpha) {
  bmv_->scale(alpha);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    update(const Scalar& alpha, const MultiVector& A, const Scalar& beta) {
  bmv_->update(alpha, A, beta);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    update(const Scalar& alpha, const MultiVector& A, const Scalar& beta, const MultiVector& B, const Scalar& gamma) {
  bmv_->update(alpha, A, beta, B, gamma);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mag_type
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm1() const {
  Teuchos::Array<mag_type> norm = Teuchos::Array<mag_type>(1);
  this->norm1(norm);
  return norm[0];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mag_type
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm2() const {
  Teuchos::Array<mag_type> norm = Teuchos::Array<mag_type>(1);
  this->norm2(norm);
  return norm[0];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::mag_type
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    normInf() const {
  Teuchos::Array<mag_type> norm = Teuchos::Array<mag_type>(1);
  this->normInf(norm);
  return norm[0];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm1(const Teuchos::ArrayView<mag_type>& norms) const {
  bmv_->norm1(norms);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm2(const Teuchos::ArrayView<mag_type>& norms) const {
  bmv_->norm2(norms);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    normInf(const Teuchos::ArrayView<mag_type>& norms) const {
  bmv_->normInf(norms);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    meanValue(const Teuchos::ArrayView<impl_scalar_type>& /* means */) const {
  throw std::runtime_error("BlockedVector::meanValue: Not (yet) supported by BlockedVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::impl_scalar_type
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    meanValue() const {
  throw std::runtime_error("BlockedVector::meanValue: Not (yet) supported by BlockedVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    multiply(Teuchos::ETransp /* transA */,
             Teuchos::ETransp /* transB */,
             const Scalar& /* alpha */,
             const Vector& /* A */,
             const Vector& /* B */,
             const Scalar& /* beta */) {
  throw std::runtime_error("BlockedVector::multiply: Not (yet) supported by BlockedVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    multiply(Teuchos::ETransp /* transA */,
             Teuchos::ETransp /* transB */,
             const Scalar& /* alpha */,
             const MultiVector& /* A */,
             const MultiVector& /* B */,
             const Scalar& /* beta */) {
  throw std::runtime_error("BlockedVector::multiply: Not (yet) supported by BlockedVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    elementWiseMultiply(Scalar /* scalarAB */,
                        const Vector& /* A */,
                        const MultiVector& /* B */,
                        Scalar /* scalarThis */) {
  throw std::runtime_error("BlockedVector::elementWiseMultiply: Not (yet) supported by BlockedVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    elementWiseMultiply(Scalar /* scalarAB */,
                        const Vector& A,
                        const Vector& B,
                        Scalar /* scalarThis */) {
  TEUCHOS_TEST_FOR_EXCEPTION(A.getMap()->getLocalNumElements() != B.getMap()->getLocalNumElements(),
                             std::runtime_error,
                             "BlockedVector::elementWiseMultiply: A has " << A.getMap()->getLocalNumElements() << " elements, B has "
                                                                          << B.getMap()->getLocalNumElements() << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(A.getMap()->getGlobalNumElements() != B.getMap()->getGlobalNumElements(),
                             std::runtime_error,
                             "BlockedVector::elementWiseMultiply: A has " << A.getMap()->getGlobalNumElements() << " elements, B has "
                                                                          << B.getMap()->getGlobalNumElements() << ".");

  Teuchos::RCP<const Vector> bmvec = Teuchos::rcpFromRef(B);
  Teuchos::RCP<const BlockedVector> bbmvec = Teuchos::rcp_dynamic_cast<const BlockedVector>(bmvec);
  TEUCHOS_TEST_FOR_EXCEPTION(bbmvec.is_null() == true,
                             std::runtime_error,
                             "BlockedVector::elementWiseMultiply: B must be a BlockedVector.");

  // Note: as in Xpetra, the per-block elementWiseMultiply is left as future work
  // (the reference Xpetra implementation only performs the compatibility checks).
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
  throw std::runtime_error(
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
    isSameSize(const MultiVector& /* vec */) const {
  throw std::runtime_error(
      "BlockedVector::isSameSize: routine not implemented. It has no value as one must iterate on the partial vectors.");
  TEUCHOS_UNREACHABLE_RETURN(false);
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
  for (size_t r = 0; r < this->getBlockedMap()->getNumMaps(); r++)
    getMultiVector(r)->describe(out, verbLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceMap(const Teuchos::RCP<const Map>& map) {
  bmv_->replaceMap(map);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    randomize(bool bUseXpetraImplementation) {
  bmv_->randomize(bUseXpetraImplementation);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    randomize(const Scalar& minVal, const Scalar& maxVal, bool bUseXpetraImplementation) {
  bmv_->randomize(minVal, maxVal, bUseXpetraImplementation);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMap() const {
  return bmv_->getMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockedMap>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getBlockedMap() const {
  return bmv_->getBlockedMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiVector>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMultiVector(size_t r) const {
  return bmv_->getMultiVector(r);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiVector>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMultiVector(size_t r, bool bThyraMode) const {
  return bmv_->getMultiVector(r, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setMultiVector(size_t r, Teuchos::RCP<const Vector> v, bool bThyraMode) {
  bmv_->setMultiVector(r, v, bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiVector>
BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Merge() const {
  return bmv_->Merge();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    assign(const MultiVector& rhs) {
  *bmv_ = rhs;  // BlockedMultiVector::operator= dispatches to its (protected) assign
}

}  // namespace Tpetra

//
// Explicit instantiation macro
//
#define TPETRA_BLOCKEDVECTOR_INSTANT(SCALAR, LO, GO, NODE) \
  template class BlockedVector<SCALAR, LO, GO, NODE>;

#endif  // TPETRA_BLOCKEDVECTOR_DEF_HPP
