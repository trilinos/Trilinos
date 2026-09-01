// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_BLOCKEDMULTIVECTOR_DEF_HPP
#define XPETRA_BLOCKEDMULTIVECTOR_DEF_HPP

#include "Xpetra_BlockedMultiVector_decl.hpp"

#include "Xpetra_Exceptions.hpp"
#include "Xpetra_BlockedMap.hpp"
#include "Xpetra_MapExtractor.hpp"
#include "Xpetra_TpetraMap.hpp"
#include "Xpetra_TpetraMultiVector.hpp"
#include "Xpetra_TpetraVector.hpp"

namespace Xpetra {

// ---------------------------------------------------------------------------
// unwrap / wrap free-function helpers (blocked-aware)
// ---------------------------------------------------------------------------
namespace BlockedMultiVectorDetails {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
unwrapMultiVector(const Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& X) {
  typedef Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedMultiVectorClass;
  Teuchos::RCP<BlockedMultiVectorClass> bx = Teuchos::rcp_dynamic_cast<BlockedMultiVectorClass>(X);
  if (!bx.is_null())
    return bx->getTpetra_BlockedMultiVector();  // Tpetra::BlockedMultiVector upcast to Tpetra::MultiVector
  return toTpetra(X);                            // plain; throws Xpetra::Exceptions::BadCast for Epetra
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
unwrapMultiVector(const Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& X) {
  typedef Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> BlockedMultiVectorClass;
  Teuchos::RCP<const BlockedMultiVectorClass> bx = Teuchos::rcp_dynamic_cast<const BlockedMultiVectorClass>(X);
  if (!bx.is_null())
    return bx->getTpetra_BlockedMultiVector();
  return toTpetra(X);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
wrapMultiVector(const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& X) {
  typedef Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TpetraBlockedMultiVectorClass;
  if (X.is_null())
    return Teuchos::null;
  Teuchos::RCP<TpetraBlockedMultiVectorClass> bx = Teuchos::rcp_dynamic_cast<TpetraBlockedMultiVectorClass>(X);
  if (!bx.is_null())
    return Teuchos::rcp(new Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(bx));
  return toXpetra(X);
}

}  // namespace BlockedMultiVectorDetails

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMultiVector(const Teuchos::RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>& map,
                       size_t NumVectors,
                       bool zeroOut) {
  vec_ = Teuchos::rcp(new tpetra_blockedmultivector_type(map->getTpetra_BlockedMap(), NumVectors, zeroOut));
  // Retain the Xpetra BlockedMap so getBlockedMap() preserves nested BlockedMap
  // identity (the Tpetra core keeps only a flattened blocked map).
  blockedMapXpetra_ = map;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMultiVector(Teuchos::RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap,
                       Teuchos::RCP<const MultiVector> v) {
  vec_ = Teuchos::rcp(new tpetra_blockedmultivector_type(bmap->getTpetra_BlockedMap(),
                                                         BlockedMultiVectorDetails::unwrapMultiVector(v)));
  blockedMapXpetra_ = bmap;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMultiVector(Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mapExtractor,
                       Teuchos::RCP<const MultiVector> v) {
  // Build an (Xpetra, hence Tpetra-backed) BlockedMap out of the MapExtractor, then delegate.
  std::vector<RCP<const Map>> maps;
  maps.reserve(mapExtractor->NumMaps());
  for (size_t r = 0; r < mapExtractor->NumMaps(); ++r)
    maps.push_back(mapExtractor->getMap(r, mapExtractor->getThyraMode()));
  RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap =
      Teuchos::rcp(new Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>(mapExtractor->getFullMap(), maps, mapExtractor->getThyraMode()));
  vec_ = Teuchos::rcp(new tpetra_blockedmultivector_type(bmap->getTpetra_BlockedMap(),
                                                         BlockedMultiVectorDetails::unwrapMultiVector(v)));
  blockedMapXpetra_ = bmap;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMultiVector(const Teuchos::RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>& map,
                       std::vector<Teuchos::RCP<MultiVector>>& vin) {
  std::vector<Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>> tvin(vin.size());
  for (size_t i = 0; i < vin.size(); ++i)
    tvin[i] = BlockedMultiVectorDetails::unwrapMultiVector(vin[i]);
  vec_ = Teuchos::rcp(new tpetra_blockedmultivector_type(map->getTpetra_BlockedMap(), tvin));
  blockedMapXpetra_ = map;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMultiVector(const Teuchos::RCP<tpetra_blockedmultivector_type>& vec)
  : vec_(vec) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ~BlockedMultiVector() {
  vec_              = Teuchos::null;
  blockedMapXpetra_ = Teuchos::null;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator=(const MultiVector& rhs) {
  assign(rhs);  // dispatch to protected virtual method
  return *this;
}

// ---------------------------------------------------------------------------
// Post-construction modification
// ---------------------------------------------------------------------------
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar& value) {
  vec_->replaceGlobalValue(globalRow, vectorIndex, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar& value) {
  vec_->sumIntoGlobalValue(globalRow, vectorIndex, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar& value) {
  vec_->replaceLocalValue(myRow, vectorIndex, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar& value) {
  vec_->sumIntoLocalValue(myRow, vectorIndex, value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    putScalar(const Scalar& value) {
  vec_->putScalar(value);
}

// ---------------------------------------------------------------------------
// Data copy / view
// ---------------------------------------------------------------------------
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVector(size_t j) const {
  return toXpetra(vec_->getVector(j));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVectorNonConst(size_t j) {
  return toXpetra(vec_->getVectorNonConst(j));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<const Scalar>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getData(size_t j) const {
  return vec_->getData(j);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<Scalar>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getDataNonConst(size_t j) {
  return vec_->getDataNonConst(j);
}

// ---------------------------------------------------------------------------
// Mathematical methods
// ---------------------------------------------------------------------------
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    dot(const MultiVector& A, const Teuchos::ArrayView<Scalar>& dots) const {
  vec_->dot(*BlockedMultiVectorDetails::unwrapMultiVector(Teuchos::rcpFromRef(A)), dots);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    abs(const MultiVector& A) {
  vec_->abs(*BlockedMultiVectorDetails::unwrapMultiVector(Teuchos::rcpFromRef(A)));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    reciprocal(const MultiVector& A) {
  vec_->reciprocal(*BlockedMultiVectorDetails::unwrapMultiVector(Teuchos::rcpFromRef(A)));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    scale(const Scalar& alpha) {
  vec_->scale(alpha);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    scale(Teuchos::ArrayView<const Scalar> alpha) {
  vec_->scale(alpha);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    update(const Scalar& alpha, const MultiVector& A, const Scalar& beta) {
  vec_->update(alpha, *BlockedMultiVectorDetails::unwrapMultiVector(Teuchos::rcpFromRef(A)), beta);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    update(const Scalar& alpha, const MultiVector& A, const Scalar& beta, const MultiVector& B, const Scalar& gamma) {
  vec_->update(alpha, *BlockedMultiVectorDetails::unwrapMultiVector(Teuchos::rcpFromRef(A)), beta,
               *BlockedMultiVectorDetails::unwrapMultiVector(Teuchos::rcpFromRef(B)), gamma);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm1(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const {
  vec_->norm1(norms);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const {
  vec_->norm2(norms);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const {
  vec_->normInf(norms);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    meanValue(const Teuchos::ArrayView<Scalar>& means) const {
  vec_->meanValue(means);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar& alpha, const MultiVector& A, const MultiVector& B, const Scalar& beta) {
  vec_->multiply(transA, transB, alpha,
                 *BlockedMultiVectorDetails::unwrapMultiVector(Teuchos::rcpFromRef(A)),
                 *BlockedMultiVectorDetails::unwrapMultiVector(Teuchos::rcpFromRef(B)), beta);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    elementWiseMultiply(Scalar scalarAB, const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const MultiVector& B, Scalar scalarThis) {
  // A is a Vector; unwrap to a Tpetra::Vector via the multivector unwrap then downcast.
  Teuchos::RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tA =
      BlockedMultiVectorDetails::unwrapMultiVector(
          Teuchos::rcp_dynamic_cast<const MultiVector>(Teuchos::rcpFromRef(A)));
  Teuchos::RCP<const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> tAv =
      Teuchos::rcp_dynamic_cast<const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(tA, true);
  vec_->elementWiseMultiply(scalarAB, *tAv,
                            *BlockedMultiVectorDetails::unwrapMultiVector(Teuchos::rcpFromRef(B)), scalarThis);
}

// ---------------------------------------------------------------------------
// Attribute access
// ---------------------------------------------------------------------------
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getNumVectors() const {
  return vec_->getNumVectors();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalLength() const {
  return vec_->getLocalLength();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalLength() const {
  return vec_->getGlobalLength();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    isSameSize(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec) const {
  return vec_->isSameSize(*BlockedMultiVectorDetails::unwrapMultiVector(Teuchos::rcpFromRef(vec)));
}

// ---------------------------------------------------------------------------
// Describable
// ---------------------------------------------------------------------------
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    description() const {
  return vec_->description();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
  vec_->describe(out, verbLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceMap(const RCP<const Map>& map) {
  blockedMapXpetra_ = Teuchos::null;
  // A blocked map replaces the whole block structure; a plain map goes to the single-block overload.
  RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap =
      Teuchos::rcp_dynamic_cast<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>(map);
  if (!bmap.is_null()) {
    vec_->replaceMap(bmap->getTpetra_BlockedMap());
    // Retain the Xpetra blocked map so getBlockedMap() preserves nested identity.
    blockedMapXpetra_ = bmap;
  } else
    vec_->replaceMap(toTpetra(map));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doImport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* source */, const Import& /* importer */, CombineMode /* CM */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::doImport: Not supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doExport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* dest */, const Import& /* importer */, CombineMode /* CM */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::doExport: Not supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doImport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* source */, const Export& /* exporter */, CombineMode /* CM */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::doImport: Not supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doExport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* dest */, const Export& /* exporter */, CombineMode /* CM */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::doExport: Not supported by BlockedMultiVector.");
}

// ---------------------------------------------------------------------------
// Xpetra specific: random / seed
// ---------------------------------------------------------------------------
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setSeed(unsigned int seed) {
  for (size_t r = 0; r < vec_->getBlockedMap()->getNumMaps(); ++r)
    getMultiVector(r)->setSeed(seed);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    randomize(bool bUseXpetraImplementation) {
  vec_->randomize(bUseXpetraImplementation);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    randomize(const Scalar& minVal, const Scalar& maxVal, bool bUseXpetraImplementation) {
  vec_->randomize(minVal, maxVal, bUseXpetraImplementation);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Xpetra_randomize() {
  Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Xpetra_randomize();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Xpetra_randomize(const Scalar& minVal, const Scalar& maxVal) {
  Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Xpetra_randomize(minVal, maxVal);
}

// ---------------------------------------------------------------------------
// Map / block accessors
// ---------------------------------------------------------------------------
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMap() const {
  return getBlockedMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getBlockedMap() const {
  if (blockedMapXpetra_.is_null())
    blockedMapXpetra_ = Teuchos::rcp(new Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>(vec_->getBlockedMap()));
  return blockedMapXpetra_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMultiVector(size_t r) const {
  return BlockedMultiVectorDetails::wrapMultiVector(vec_->getMultiVector(r));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMultiVector(size_t r, bool bThyraMode) const {
  return BlockedMultiVectorDetails::wrapMultiVector(vec_->getMultiVector(r, bThyraMode));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setMultiVector(size_t r, Teuchos::RCP<const MultiVector> v, bool bThyraMode) {
  vec_->setMultiVector(r, BlockedMultiVectorDetails::unwrapMultiVector(v), bThyraMode);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Merge() const {
  return BlockedMultiVectorDetails::wrapMultiVector(vec_->Merge());
}

// ---------------------------------------------------------------------------
// Assignment
// ---------------------------------------------------------------------------
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    assign(const MultiVector& rhs) {
  // Tpetra::BlockedMultiVector::assign is protected; go through its operator= which
  // dispatches to assign() and performs the deep copy over the block structure.
  *vec_ = *BlockedMultiVectorDetails::unwrapMultiVector(Teuchos::rcpFromRef(rhs));
}

}  // namespace Xpetra

#endif  // XPETRA_BLOCKEDMULTIVECTOR_DEF_HPP
