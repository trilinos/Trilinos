// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_TPETRAMULTIVECTOR_DEF_HPP
#define XPETRA_TPETRAMULTIVECTOR_DEF_HPP
#include "Xpetra_TpetraConfigDefs.hpp"

#include "Xpetra_TpetraMap.hpp"  //TMP
#include "Xpetra_Utils.hpp"
#include "Xpetra_TpetraImport.hpp"
#include "Xpetra_TpetraExport.hpp"

#include "Xpetra_TpetraMultiVector_decl.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Details_Random.hpp"

namespace Xpetra {

//! Basic constuctor.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::TpetraMultiVector(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, size_t NumVectors, bool zeroOut)
  : vec_(Teuchos::rcp(new Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(toTpetra(map), NumVectors, zeroOut))) {
  // TAW 1/30/2016: even though Tpetra allows numVecs == 0, Epetra does not. Introduce exception to keep behavior of Epetra and Tpetra consistent.
  TEUCHOS_TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument, "Xpetra::TpetraMultiVector(map,numVecs,zeroOut): numVecs = " << NumVectors << " < 1.");
}

//! Copy constructor (performs a deep copy).
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraMultiVector(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &source, const Teuchos::DataAccess copyOrView)
  : vec_(Teuchos::rcp(new Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(toTpetra(source), copyOrView))) {}

//! Create multivector by copying two-dimensional array of local data.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraMultiVector(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, const Teuchos::ArrayView<const Scalar> &A, size_t LDA, size_t NumVectors)
  : vec_(Teuchos::rcp(new Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(toTpetra(map), A, LDA, NumVectors))) {
  // TAW 1/30/2016: even though Tpetra allows numVecs == 0, Epetra does not. Introduce exception to keep behavior of Epetra and Tpetra consistent.
  TEUCHOS_TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument, "Xpetra::TpetraMultiVector(map,A,LDA,numVecs): numVecs = " << NumVectors << " < 1.");
}

//! Create multivector by copying array of views of local data.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraMultiVector(const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map, const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &ArrayOfPtrs, size_t NumVectors)
  : vec_(Teuchos::rcp(new Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(toTpetra(map), ArrayOfPtrs, NumVectors))) {
  // TAW 1/30/2016: even though Tpetra allows numVecs == 0, Epetra does not. Introduce exception to keep behavior of Epetra and Tpetra consistent.
  TEUCHOS_TEST_FOR_EXCEPTION(NumVectors < 1, std::invalid_argument, "Xpetra::TpetraMultiVector(map,ArrayOfPtrs,numVecs): numVecs = " << NumVectors << " < 1.");
}

//! Destructor (virtual for memory safety of derived classes).
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ~TpetraMultiVector() {}

//! Replace value, using global (row) index.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) {
  XPETRA_MONITOR("TpetraMultiVector::replaceGlobalValue");
  vec_->replaceGlobalValue(globalRow, vectorIndex, value);
}

//! Add value to existing value, using global (row) index.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) {
  XPETRA_MONITOR("TpetraMultiVector::sumIntoGlobalValue");
  vec_->sumIntoGlobalValue(globalRow, vectorIndex, value);
}

//! Replace value, using local (row) index.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) {
  XPETRA_MONITOR("TpetraMultiVector::replaceLocalValue");
  vec_->replaceLocalValue(myRow, vectorIndex, value);
}

//! Add value to existing value, using local (row) index
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) {
  XPETRA_MONITOR("TpetraMultiVector::sumIntoLocalValue");
  vec_->sumIntoLocalValue(myRow, vectorIndex, value);
}

//! Set all values in the multivector with the given value
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    putScalar(const Scalar &value) {
  XPETRA_MONITOR("TpetraMultiVector::putScalar");
  vec_->putScalar(value);
}

//! Sum values of a locally replicated multivector across all processes.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    reduce() {
  XPETRA_MONITOR("TpetraMultiVector::reduce");
  vec_->reduce();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVector(size_t j) const {
  XPETRA_MONITOR("TpetraMultiVector::getVector");
  return toXpetra(vec_->getVector(j));
}

//! Return a Vector which is a nonconst view of column j.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVectorNonConst(size_t j) {
  XPETRA_MONITOR("TpetraMultiVector::getVectorNonConst");
  return toXpetra(vec_->getVectorNonConst(j));
}

//! Const view of the local values in a particular vector of this multivector.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<const Scalar>
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getData(size_t j) const {
  XPETRA_MONITOR("TpetraMultiVector::getData");
  return vec_->getData(j);
}

//! View of the local values in a particular vector of this multivector.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<Scalar>
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getDataNonConst(size_t j) {
  XPETRA_MONITOR("TpetraMultiVector::getDataNonConst");
  return vec_->getDataNonConst(j);
}

//! Fill the given array with a copy of this multivector's local values.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    get1dCopy(Teuchos::ArrayView<Scalar> A, size_t LDA) const {
  XPETRA_MONITOR("TpetraMultiVector::get1dCopy");
  vec_->get1dCopy(A, LDA);
}

//! Fill the given array with a copy of this multivector's local values.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    get2dCopy(Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > ArrayOfPtrs) const {
  XPETRA_MONITOR("TpetraMultiVector::get2dCopy");
  vec_->get2dCopy(ArrayOfPtrs);
}

//! Const persisting (1-D) view of this multivector's local values.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<const Scalar>
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    get1dView() const {
  XPETRA_MONITOR("TpetraMultiVector::get1dView");
  return vec_->get1dView();
}

//! Return const persisting pointers to values.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> >
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    get2dView() const {
  XPETRA_MONITOR("TpetraMultiVector::get2dView");
  return vec_->get2dView();
}

//! Nonconst persisting (1-D) view of this multivector's local values.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<Scalar>
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    get1dViewNonConst() {
  XPETRA_MONITOR("TpetraMultiVector::get1dViewNonConst");
  return vec_->get1dViewNonConst();
}

//! Return non-const persisting pointers to values.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    get2dViewNonConst() {
  XPETRA_MONITOR("TpetraMultiVector::get2dViewNonConst");
  return vec_->get2dViewNonConst();
}

//! Compute dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    dot(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A, const Teuchos::ArrayView<Scalar> &dots) const {
  XPETRA_MONITOR("TpetraMultiVector::dot");
  vec_->dot(toTpetra(A), dots);
}

//! Put element-wise absolute values of input Multi-vector in target: A = abs(this).
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    abs(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A) {
  XPETRA_MONITOR("TpetraMultiVector::abs");
  vec_->abs(toTpetra(A));
}

//! Put element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    reciprocal(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A) {
  XPETRA_MONITOR("TpetraMultiVector::reciprocal");
  vec_->reciprocal(toTpetra(A));
}

//! Scale the current values of a multi-vector, this = alpha*this.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    scale(const Scalar &alpha) {
  XPETRA_MONITOR("TpetraMultiVector::scale");
  vec_->scale(alpha);
}

//! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    scale(Teuchos::ArrayView<const Scalar> alpha) {
  XPETRA_MONITOR("TpetraMultiVector::scale");
  vec_->scale(alpha);
}

//! Replace multi-vector values with scaled values of A, this = alpha*A.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    scale(const Scalar &alpha, const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A) {
  XPETRA_MONITOR("TpetraMultiVector::scale");
  vec_->scale(alpha, toTpetra(A));
}

//! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    update(const Scalar &alpha, const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A, const Scalar &beta) {
  XPETRA_MONITOR("TpetraMultiVector::update");
  vec_->update(alpha, toTpetra(A), beta);
}

//! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    update(const Scalar &alpha, const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A, const Scalar &beta, const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B, const Scalar &gamma) {
  XPETRA_MONITOR("TpetraMultiVector::update");
  vec_->update(alpha, toTpetra(A), beta, toTpetra(B), gamma);
}

//! Compute 1-norm of each vector in multi-vector.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm1(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const {
  XPETRA_MONITOR("TpetraMultiVector::norm1");
  vec_->norm1(norms);
}

//!
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const {
  XPETRA_MONITOR("TpetraMultiVector::norm2");
  vec_->norm2(norms);
}

//! Compute Inf-norm of each vector in multi-vector.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const {
  XPETRA_MONITOR("TpetraMultiVector::normInf");
  vec_->normInf(norms);
}

//! Compute mean (average) value of each vector in multi-vector. The outcome of this routine is undefined for non-floating point scalar types (e.g., int).
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    meanValue(const Teuchos::ArrayView<Scalar> &means) const {
  XPETRA_MONITOR("TpetraMultiVector::meanValue");
  vec_->meanValue(means);
}

//! Matrix-matrix multiplication: this = beta*this + alpha*op(A)*op(B).
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A, const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B, const Scalar &beta) {
  XPETRA_MONITOR("TpetraMultiVector::multiply");
  vec_->multiply(transA, transB, alpha, toTpetra(A), toTpetra(B), beta);
}

//! Number of columns in the multivector.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getNumVectors() const {
  XPETRA_MONITOR("TpetraMultiVector::getNumVectors");
  return vec_->getNumVectors();
}

//! Local number of rows on the calling process.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalLength() const {
  XPETRA_MONITOR("TpetraMultiVector::getLocalLength");
  return vec_->getLocalLength();
}

//! Global number of rows in the multivector.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalLength() const {
  XPETRA_MONITOR("TpetraMultiVector::getGlobalLength");
  return vec_->getGlobalLength();
}

// \brief Checks to see if the local length, number of vectors and size of Scalar type match
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    isSameSize(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &vec) const {
  XPETRA_MONITOR("TpetraMultiVector::isSameSize");
  return vec_->isSameSize(toTpetra(vec));
}

//! A simple one-line description of this object.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    description() const {
  XPETRA_MONITOR("TpetraMultiVector::description");
  return vec_->description();
}

//! Print the object with the given verbosity level to a FancyOStream.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
  XPETRA_MONITOR("TpetraMultiVector::describe");
  vec_->describe(out, verbLevel);
}

//! Set multi-vector values to random numbers.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    randomize(bool bUseXpetraImplementation) {
  XPETRA_MONITOR("TpetraMultiVector::randomize");

  if (bUseXpetraImplementation)
    MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Xpetra_randomize();
  else
    vec_->randomize();
}

//! Set multi-vector values to random numbers.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    randomize(const Scalar &minVal, const Scalar &maxVal, bool bUseXpetraImplementation) {
  XPETRA_MONITOR("TpetraMultiVector::randomize");

  if (bUseXpetraImplementation)
    MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Xpetra_randomize(minVal, maxVal);
  else
    vec_->randomize(minVal, maxVal);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMap() const {
  XPETRA_MONITOR("TpetraMultiVector::getMap");
  return toXpetra(vec_->getMap());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doImport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> &source, const Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) {
  XPETRA_MONITOR("TpetraMultiVector::doImport");

  XPETRA_DYNAMIC_CAST(const TpetraMultiVectorClass, source, tSource, "Xpetra::TpetraMultiVector::doImport only accept Xpetra::TpetraMultiVector as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tSource.getTpetra_MultiVector();
  this->getTpetra_MultiVector()->doImport(*v, toTpetra(importer), toTpetra(CM));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    beginImport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> &source, const Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) {
  XPETRA_MONITOR("TpetraMultiVector::beginImport");

  XPETRA_DYNAMIC_CAST(const TpetraMultiVectorClass, source, tSource, "Xpetra::TpetraMultiVector::doImport only accept Xpetra::TpetraMultiVector as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tSource.getTpetra_MultiVector();
  this->getTpetra_MultiVector()->beginImport(*v, toTpetra(importer), toTpetra(CM));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    endImport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> &source, const Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) {
  XPETRA_MONITOR("TpetraMultiVector::endImport");

  XPETRA_DYNAMIC_CAST(const TpetraMultiVectorClass, source, tSource, "Xpetra::TpetraMultiVector::doImport only accept Xpetra::TpetraMultiVector as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tSource.getTpetra_MultiVector();
  this->getTpetra_MultiVector()->endImport(*v, toTpetra(importer), toTpetra(CM));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doExport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> &dest, const Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) {
  XPETRA_MONITOR("TpetraMultiVector::doExport");

  XPETRA_DYNAMIC_CAST(const TpetraMultiVectorClass, dest, tDest, "Xpetra::TpetraMultiVector::doImport only accept Xpetra::TpetraMultiVector as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tDest.getTpetra_MultiVector();
  this->getTpetra_MultiVector()->doExport(*v, toTpetra(importer), toTpetra(CM));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    beginExport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> &dest, const Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) {
  XPETRA_MONITOR("TpetraMultiVector::beginExport");

  XPETRA_DYNAMIC_CAST(const TpetraMultiVectorClass, dest, tDest, "Xpetra::TpetraMultiVector::doImport only accept Xpetra::TpetraMultiVector as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tDest.getTpetra_MultiVector();
  this->getTpetra_MultiVector()->beginExport(*v, toTpetra(importer), toTpetra(CM));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    endExport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> &dest, const Import<LocalOrdinal, GlobalOrdinal, Node> &importer, CombineMode CM) {
  XPETRA_MONITOR("TpetraMultiVector::endExport");

  XPETRA_DYNAMIC_CAST(const TpetraMultiVectorClass, dest, tDest, "Xpetra::TpetraMultiVector::doImport only accept Xpetra::TpetraMultiVector as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tDest.getTpetra_MultiVector();
  this->getTpetra_MultiVector()->endExport(*v, toTpetra(importer), toTpetra(CM));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doImport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> &source, const Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM) {
  XPETRA_MONITOR("TpetraMultiVector::doImport");

  XPETRA_DYNAMIC_CAST(const TpetraMultiVectorClass, source, tSource, "Xpetra::TpetraMultiVector::doImport only accept Xpetra::TpetraMultiVector as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tSource.getTpetra_MultiVector();
  this->getTpetra_MultiVector()->doImport(*v, toTpetra(exporter), toTpetra(CM));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    beginImport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> &source, const Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM) {
  XPETRA_MONITOR("TpetraMultiVector::beginImport");

  XPETRA_DYNAMIC_CAST(const TpetraMultiVectorClass, source, tSource, "Xpetra::TpetraMultiVector::doImport only accept Xpetra::TpetraMultiVector as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tSource.getTpetra_MultiVector();
  this->getTpetra_MultiVector()->beginImport(*v, toTpetra(exporter), toTpetra(CM));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    endImport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> &source, const Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM) {
  XPETRA_MONITOR("TpetraMultiVector::endImport");

  XPETRA_DYNAMIC_CAST(const TpetraMultiVectorClass, source, tSource, "Xpetra::TpetraMultiVector::doImport only accept Xpetra::TpetraMultiVector as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tSource.getTpetra_MultiVector();
  this->getTpetra_MultiVector()->endImport(*v, toTpetra(exporter), toTpetra(CM));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doExport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> &dest, const Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM) {
  XPETRA_MONITOR("TpetraMultiVector::doExport");

  XPETRA_DYNAMIC_CAST(const TpetraMultiVectorClass, dest, tDest, "Xpetra::TpetraMultiVector::doImport only accept Xpetra::TpetraMultiVector as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tDest.getTpetra_MultiVector();
  this->getTpetra_MultiVector()->doExport(*v, toTpetra(exporter), toTpetra(CM));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    beginExport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> &dest, const Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM) {
  XPETRA_MONITOR("TpetraMultiVector::beginExport");

  XPETRA_DYNAMIC_CAST(const TpetraMultiVectorClass, dest, tDest, "Xpetra::TpetraMultiVector::doImport only accept Xpetra::TpetraMultiVector as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tDest.getTpetra_MultiVector();
  this->getTpetra_MultiVector()->beginExport(*v, toTpetra(exporter), toTpetra(CM));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    endExport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node> &dest, const Export<LocalOrdinal, GlobalOrdinal, Node> &exporter, CombineMode CM) {
  XPETRA_MONITOR("TpetraMultiVector::endExport");

  XPETRA_DYNAMIC_CAST(const TpetraMultiVectorClass, dest, tDest, "Xpetra::TpetraMultiVector::doImport only accept Xpetra::TpetraMultiVector as input arguments.");  // TODO: remove and use toTpetra()
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > v = tDest.getTpetra_MultiVector();
  this->getTpetra_MultiVector()->endExport(*v, toTpetra(exporter), toTpetra(CM));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceMap(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &map) {
  XPETRA_MONITOR("TpetraMultiVector::replaceMap");
  this->getTpetra_MultiVector()->replaceMap(toTpetra(map));
}

//! TpetraMultiVector constructor to wrap a Tpetra::MultiVector objecT
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    TpetraMultiVector(const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > &vec)
  : vec_(vec) {}  // TODO removed const

//! Get the underlying Tpetra multivector
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getTpetra_MultiVector() const { return vec_; }

//! Set seed for Random function.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setSeed(unsigned int seed) {
  XPETRA_MONITOR("TpetraMultiVector::seedrandom");
  Teuchos::ScalarTraits<Scalar>::seedrandom(seed);
  // Tell Tpetra to update its RNG pool for the new random seed
  Tpetra::Details::Static_Random_XorShift64_Pool<typename Node::device_type::execution_space>::resetPool(getMap()->getComm()->getRank());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type::t_host_const_um
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalViewHost(Tpetra::Access::ReadOnlyStruct) const {
  return subview(vec_->getLocalViewHost(Tpetra::Access::ReadOnly), Kokkos::ALL(), Kokkos::ALL());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type::t_dev_const_um
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalViewDevice(Tpetra::Access::ReadOnlyStruct) const {
  return subview(vec_->getLocalViewDevice(Tpetra::Access::ReadOnly), Kokkos::ALL(), Kokkos::ALL());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type::t_host_um
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalViewHost(Tpetra::Access::OverwriteAllStruct) const {
  return subview(vec_->getLocalViewHost(Tpetra::Access::OverwriteAll), Kokkos::ALL(), Kokkos::ALL());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type::t_dev_um
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalViewDevice(Tpetra::Access::OverwriteAllStruct) const {
  return subview(vec_->getLocalViewDevice(Tpetra::Access::OverwriteAll), Kokkos::ALL(), Kokkos::ALL());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type::t_host_um
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalViewHost(Tpetra::Access::ReadWriteStruct) const {
  return subview(vec_->getLocalViewHost(Tpetra::Access::ReadWrite), Kokkos::ALL(), Kokkos::ALL());
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::dual_view_type::t_dev_um
TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalViewDevice(Tpetra::Access::ReadWriteStruct) const {
  return subview(vec_->getLocalViewDevice(Tpetra::Access::ReadWrite), Kokkos::ALL(), Kokkos::ALL());
}

/// \brief Implementation of the assignment operator (operator=);
///   does a deep copy.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    assign(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &rhs) {
  typedef TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> this_type;
  const this_type *rhsPtr = dynamic_cast<const this_type *>(&rhs);
  TEUCHOS_TEST_FOR_EXCEPTION(
      rhsPtr == NULL, std::invalid_argument,
      "Xpetra::MultiVector::operator=:"
      " The left-hand side (LHS) of the assignment has a different type than "
      "the right-hand side (RHS).  The LHS has type Xpetra::TpetraMultiVector"
      " (which means it wraps a Tpetra::MultiVector), but the RHS has some "
      "other type.  This probably means that the RHS wraps an "
      "Epetra_MultiVector.  Xpetra::MultiVector does not currently implement "
      "assignment from an Epetra object to a Tpetra object, though this could"
      " be added with sufficient interest.");

  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TMV;
  RCP<const TMV> rhsImpl = rhsPtr->getTpetra_MultiVector();
  RCP<TMV> lhsImpl       = this->getTpetra_MultiVector();

  TEUCHOS_TEST_FOR_EXCEPTION(
      rhsImpl.is_null(), std::logic_error,
      "Xpetra::MultiVector::operator= "
      "(in Xpetra::TpetraMultiVector::assign): *this (the right-hand side of "
      "the assignment) has a null RCP<Tpetra::MultiVector> inside.  Please "
      "report this bug to the Xpetra developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      lhsImpl.is_null(), std::logic_error,
      "Xpetra::MultiVector::operator= "
      "(in Xpetra::TpetraMultiVector::assign): The left-hand side of the "
      "assignment has a null RCP<Tpetra::MultiVector> inside.  Please report "
      "this bug to the Xpetra developers.");

  Tpetra::deep_copy(*lhsImpl, *rhsImpl);
}

}  // namespace Xpetra

// Following header file inculsion is needed for the dynamic_cast to TpetraVector in
// elementWiseMultiply (because we cannot dynamic_cast if target is not a complete type)
// It is included here to avoid circular dependency between Vector and MultiVector
// TODO: there is certainly a more elegant solution...
#include "Xpetra_TpetraVector.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    elementWiseMultiply(Scalar scalarAB,
                        const Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &A,
                        const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> &B,
                        Scalar scalarThis) {
  XPETRA_MONITOR("TpetraMultiVector::elementWiseMultiply");

  // XPETRA_DYNAMIC_CAST won't take TpetraVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>
  // as an argument, hence the following typedef.
  typedef TpetraVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tpv;
  XPETRA_DYNAMIC_CAST(const tpv, A, tA, "Xpetra::TpetraMultiVectorMatrix->multiply() only accept Xpetra::TpetraMultiVector as input arguments.");
  XPETRA_DYNAMIC_CAST(const TpetraMultiVector, B, tB, "Xpetra::TpetraMultiVectorMatrix->multiply() only accept Xpetra::TpetraMultiVector as input arguments.");
  vec_->elementWiseMultiply(scalarAB, *tA.getTpetra_Vector(), *tB.getTpetra_MultiVector(), scalarThis);
}

}  // namespace Xpetra

#define XPETRA_TPETRAMULTIVECTOR_SHORT
#endif  // XPETRA_TPETRAMULTIVECTOR_HPP
