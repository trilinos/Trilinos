// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKEDMULTIVECTOR_DEF_HPP
#define TPETRA_BLOCKEDMULTIVECTOR_DEF_HPP

#include "Tpetra_BlockedMultiVector_decl.hpp"  // completes BlockedMultiVector before cross-includes

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_BlockedMap_def.hpp"
#include "Tpetra_BlockedVector_def.hpp"   // getVector()/elementWiseMultiply construct/cast BlockedVector
#include "Tpetra_MapExtractor_def.hpp"    // ctor from MapExtractor calls its methods

#include <Teuchos_ScalarTraits.hpp>

namespace Tpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMultiVector(const Teuchos::RCP<const BlockedMap>& map,
                       size_t NumVectors,
                       bool zeroOut)
  : base_type()
  , map_(map) {
  numVectors_ = NumVectors;

  vv_.reserve(map->getNumMaps());

  // add sub-multivectors in row order
  for (size_t r = 0; r < map->getNumMaps(); ++r) {
    // If sub-block r is itself a BlockedMap (a nested blocked operator),
    // build a nested BlockedMultiVector so its dynamic type is preserved --
    // mirroring Xpetra, where MultiVectorFactory::Build(BlockedMap) returned a
    // nested BlockedMultiVector.  Leaf sub-blocks build a plain MultiVector.
    Teuchos::RCP<const BlockedMap> nested = map->getBlockedSubMap(r);
    if (!nested.is_null())
      vv_.push_back(Teuchos::rcp(new BlockedMultiVector(nested, NumVectors, zeroOut)));
    else
      vv_.push_back(Teuchos::rcp(new MultiVector(map->getMap(r, map_->getThyraMode()), NumVectors, zeroOut)));
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMultiVector(Teuchos::RCP<const BlockedMap> bmap,
                       Teuchos::RCP<const MultiVector> v)
  : base_type() {
  TEUCHOS_TEST_FOR_EXCEPTION(bmap->getFullMap()->getLocalNumElements() != v->getMap()->getLocalNumElements(),
                             std::runtime_error,
                             "BlockedMultiVector: inconsistent number of local elements of MultiVector and BlockedMap. The BlockedMap has "
                                 << bmap->getFullMap()->getLocalNumElements() << " local elements. The vector has "
                                 << v->getMap()->getLocalNumElements() << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(bmap->getFullMap()->getGlobalNumElements() != v->getMap()->getGlobalNumElements(),
                             std::runtime_error,
                             "BlockedMultiVector: inconsistent number of global elements of MultiVector and BlockedMap. The BlockedMap has "
                                 << bmap->getFullMap()->getGlobalNumElements() << " global elements. The vector has "
                                 << v->getMap()->getGlobalNumElements() << ".");

  numVectors_ = v->getNumVectors();
  map_        = bmap;

  vv_.reserve(bmap->getNumMaps());
  for (size_t r = 0; r < bmap->getNumMaps(); ++r)
    vv_.push_back(this->ExtractVector(v, r, bmap->getThyraMode()));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMultiVector(Teuchos::RCP<const MapExtractor> mapExtractor,
                       Teuchos::RCP<const MultiVector> v)
  : base_type() {
  numVectors_ = v->getNumVectors();

  // create blocked map out of MapExtractor
  std::vector<Teuchos::RCP<const Map>> maps;
  maps.reserve(mapExtractor->NumMaps());
  for (size_t r = 0; r < mapExtractor->NumMaps(); ++r)
    maps.push_back(mapExtractor->getMap(r, mapExtractor->getThyraMode()));
  map_ = Teuchos::rcp(new BlockedMap(mapExtractor->getFullMap(), maps, mapExtractor->getThyraMode()));

  vv_.reserve(mapExtractor->NumMaps());
  for (size_t r = 0; r < mapExtractor->NumMaps(); ++r)
    vv_.push_back(this->ExtractVector(v, r, mapExtractor->getThyraMode()));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMultiVector(const Teuchos::RCP<const BlockedMap>& map,
                       std::vector<Teuchos::RCP<MultiVector>>& vin)
  : base_type() {
  numVectors_ = vin[0]->getNumVectors();
  map_        = map;
  vv_.resize(vin.size());
  for (size_t i = 0; i < vv_.size(); i++)
    vv_[i] = vin[i];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ~BlockedMultiVector() {
  for (size_t r = 0; r < vv_.size(); ++r)
    vv_[r] = Teuchos::null;
  map_        = Teuchos::null;
  numVectors_ = 0;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator=(const MultiVector& rhs) {
  assign(rhs);  // dispatch to protected virtual method
  return *this;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceGlobalValue(GlobalOrdinal /* globalRow */, size_t /* vectorIndex */, const Scalar& /* value */) {
  throw std::runtime_error("BlockedMultiVector::replaceGlobalValue: Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoGlobalValue(GlobalOrdinal /* globalRow */, size_t /* vectorIndex */, const Scalar& /* value */) {
  throw std::runtime_error("BlockedMultiVector::sumIntoGlobalValue: Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceLocalValue(LocalOrdinal /* myRow */, size_t /* vectorIndex */, const Scalar& /* value */) {
  throw std::runtime_error("BlockedMultiVector::replaceLocalValue: Not supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoLocalValue(LocalOrdinal /* myRow */, size_t /* vectorIndex */, const Scalar& /* value */) {
  throw std::runtime_error("BlockedMultiVector::sumIntoLocalValue: Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    putScalar(const Scalar& value) {
  for (size_t r = 0; r < map_->getNumMaps(); r++)
    getMultiVector(r)->putScalar(value);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Vector>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVector(size_t j) const {
  using BlockedVector = ::Tpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  Teuchos::RCP<BlockedVector> ret = Teuchos::rcp(new BlockedVector(this->getBlockedMap(), false));

  for (size_t r = 0; r < getBlockedMap()->getNumMaps(); r++) {
    Teuchos::RCP<const Vector> subvec =
        this->getMultiVector(r, this->getBlockedMap()->getThyraMode())->getVector(j);
    ret->setMultiVector(r, subvec, this->getBlockedMap()->getThyraMode());
  }
  return ret;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Vector>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVectorNonConst(size_t j) {
  return Teuchos::rcp_const_cast<Vector>(getVector(j));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<const Scalar>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getData(size_t j) const {
  if (map_->getNumMaps() == 1)
    return vv_[0]->getData(j);
  throw std::runtime_error("BlockedMultiVector::getData: Not (yet) supported by BlockedMultiVector.");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<Scalar>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getDataNonConst(size_t j) {
  if (map_->getNumMaps() == 1)
    return vv_[0]->getDataNonConst(j);
  throw std::runtime_error("BlockedMultiVector::getDataNonConst: Not (yet) supported by BlockedMultiVector.");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    dot(const MultiVector& /* A */, const Teuchos::ArrayView<dot_type>& /* dots */) const {
  throw std::runtime_error("BlockedMultiVector::dot: Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    abs(const MultiVector& /* A */) {
  throw std::runtime_error("BlockedMultiVector::abs: Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    reciprocal(const MultiVector& /* A */) {
  throw std::runtime_error("BlockedMultiVector::reciprocal: Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    scale(const Scalar& alpha) {
  for (size_t r = 0; r < map_->getNumMaps(); ++r)
    if (getMultiVector(r) != Teuchos::null)
      getMultiVector(r)->scale(alpha);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    scale(const Teuchos::ArrayView<const Scalar>& alpha) {
  for (size_t r = 0; r < map_->getNumMaps(); ++r)
    if (getMultiVector(r) != Teuchos::null)
      getMultiVector(r)->scale(alpha);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    update(const Scalar& alpha, const MultiVector& A, const Scalar& beta) {
  Teuchos::RCP<const MultiVector> rcpA      = Teuchos::rcpFromRef(A);
  Teuchos::RCP<const BlockedMultiVector> bA = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpA);
  // NOTE (Tpetra port): Tpetra::MultiVector::getNumVectors() is NOT virtual (unlike
  // Xpetra's), so querying it through an RCP<const MultiVector> that actually refers to a
  // BlockedMultiVector reports the empty base object's column count.  Query the true count
  // through the blocked interface when A is itself a BlockedMultiVector.
  const size_t aNumVectors = bA.is_null() ? rcpA->getNumVectors() : bA->getNumVectors();
  TEUCHOS_TEST_FOR_EXCEPTION(numVectors_ != aNumVectors,
                             std::runtime_error,
                             "BlockedMultiVector::update: update with incompatible vector (different number of vectors in multivector).");
  if (bA != Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() != bA->getBlockedMap()->getThyraMode(),
                               std::runtime_error,
                               "BlockedMultiVector::update: update with incompatible vector (different thyra mode).");
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getNumMaps() != bA->getBlockedMap()->getNumMaps(),
                               std::runtime_error,
                               "BlockedMultiVector::update: update with incompatible vector (different number of partial vectors).");
    for (size_t r = 0; r < map_->getNumMaps(); r++) {
      Teuchos::RCP<MultiVector> lmv = getMultiVector(r);
      Teuchos::RCP<MultiVector> rmv = bA->getMultiVector(r);

      Teuchos::RCP<BlockedMultiVector> blmv = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(lmv);
      Teuchos::RCP<BlockedMultiVector> brmv = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(rmv);

      if (blmv.is_null() == true && brmv.is_null() == false) {
        // special case: lmv is standard MultiVector but rmv is BlockedMultiVector
        TEUCHOS_TEST_FOR_EXCEPTION(brmv->getBlockedMap()->getNumMaps() > 1,
                                   std::runtime_error,
                                   "BlockedMultiVector::update: Standard MultiVector object does not accept BlockedMultVector object as "
                                   "parameter in update call.");
        lmv->update(alpha, *(brmv->getMultiVector(0)), beta);
      } else
        lmv->update(alpha, *rmv, beta);
    }
  } else {
    if (getBlockedMap()->getNumMaps() == 1) {
      getMultiVector(0)->update(alpha, *rcpA, beta);
    } else {
      for (size_t r = 0; r < map_->getNumMaps(); r++) {
        Teuchos::RCP<const MultiVector> part = this->ExtractVector(rcpA, r, map_->getThyraMode());
        getMultiVector(r)->update(alpha, *part, beta);
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    update(const Scalar& alpha, const MultiVector& A, const Scalar& beta, const MultiVector& B, const Scalar& gamma) {
  Teuchos::RCP<const MultiVector> rcpA      = Teuchos::rcpFromRef(A);
  Teuchos::RCP<const MultiVector> rcpB      = Teuchos::rcpFromRef(B);
  Teuchos::RCP<const BlockedMultiVector> bA = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpA);
  Teuchos::RCP<const BlockedMultiVector> bB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpB);
  if (bA != Teuchos::null && bB != Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getNumMaps() != bA->getBlockedMap()->getNumMaps(),
                               std::runtime_error,
                               "BlockedMultiVector::update: update with incompatible vector (different number of partial vectors in vector A).");
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getNumMaps() != bB->getBlockedMap()->getNumMaps(),
                               std::runtime_error,
                               "BlockedMultiVector::update: update with incompatible vector (different number of partial vectors in vector B).");
    for (size_t r = 0; r < map_->getNumMaps(); r++)
      getMultiVector(r)->update(alpha, *(bA->getMultiVector(r)), beta, *(bB->getMultiVector(r)), gamma);
    return;
  }
  throw std::runtime_error("BlockedMultiVector::update: only supports update with other BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm1(const Teuchos::ArrayView<mag_type>& norms) const {
  Teuchos::Array<mag_type> temp_norms(getNumVectors());
  std::fill(norms.begin(), norms.end(), Teuchos::ScalarTraits<mag_type>::zero());
  std::fill(temp_norms.begin(), temp_norms.end(), Teuchos::ScalarTraits<mag_type>::zero());
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    if (getMultiVector(r) != Teuchos::null) {
      getMultiVector(r)->norm1(temp_norms());
      for (size_t c = 0; c < getNumVectors(); ++c)
        norms[c] += temp_norms[c];
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm2(const Teuchos::ArrayView<mag_type>& norms) const {
  Teuchos::Array<mag_type> results(getNumVectors());
  Teuchos::Array<mag_type> temp_norms(getNumVectors());
  std::fill(norms.begin(), norms.end(), Teuchos::ScalarTraits<mag_type>::zero());
  std::fill(results.begin(), results.end(), Teuchos::ScalarTraits<mag_type>::zero());
  std::fill(temp_norms.begin(), temp_norms.end(), Teuchos::ScalarTraits<mag_type>::zero());
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    if (getMultiVector(r) != Teuchos::null) {
      getMultiVector(r)->norm2(temp_norms());
      for (size_t c = 0; c < getNumVectors(); ++c)
        results[c] += temp_norms[c] * temp_norms[c];
    }
  }
  for (size_t c = 0; c < getNumVectors(); ++c)
    norms[c] = Teuchos::ScalarTraits<mag_type>::squareroot(results[c]);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    normInf(const Teuchos::ArrayView<mag_type>& norms) const {
  Teuchos::Array<mag_type> temp_norms(getNumVectors());
  std::fill(norms.begin(), norms.end(), Teuchos::ScalarTraits<mag_type>::zero());
  std::fill(temp_norms.begin(), temp_norms.end(), Teuchos::ScalarTraits<mag_type>::zero());
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    if (getMultiVector(r) != Teuchos::null) {
      getMultiVector(r)->normInf(temp_norms());
      for (size_t c = 0; c < getNumVectors(); ++c)
        norms[c] = std::max(norms[c], temp_norms[c]);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    meanValue(const Teuchos::ArrayView<impl_scalar_type>& /* means */) const {
  throw std::runtime_error("BlockedMultiVector::meanValue: Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    multiply(Teuchos::ETransp /* transA */,
             Teuchos::ETransp /* transB */,
             const Scalar& /* alpha */,
             const MultiVector& /* A */,
             const MultiVector& /* B */,
             const Scalar& /* beta */) {
  throw std::runtime_error("BlockedMultiVector::multiply: Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    elementWiseMultiply(Scalar scalarAB,
                        const Vector& A,
                        const MultiVector& B,
                        Scalar scalarThis) {
  TEUCHOS_TEST_FOR_EXCEPTION(A.getMap()->getLocalNumElements() != B.getMap()->getLocalNumElements(),
                             std::runtime_error,
                             "BlockedMultiVector::elementWiseMultiply: A has " << A.getMap()->getLocalNumElements() << " elements, B has "
                                                                               << B.getMap()->getLocalNumElements() << ".");

  Teuchos::RCP<const BlockedMap> bmap = getBlockedMap();
  Teuchos::RCP<const Vector> rcpA     = Teuchos::rcpFromRef(A);
  using BlockedVector                 = ::Tpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  Teuchos::RCP<const BlockedVector> bA = Teuchos::rcp_dynamic_cast<const BlockedVector>(rcpA);

  Teuchos::RCP<const MultiVector> rcpB      = Teuchos::rcpFromRef(B);
  Teuchos::RCP<const BlockedMultiVector> bB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpB);
  TEUCHOS_TEST_FOR_EXCEPTION(
      bB.is_null() == true, std::runtime_error, "BlockedMultiVector::elementWiseMultiply: B must be a BlockedMultiVector.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      bA.is_null() == true, std::runtime_error, "BlockedMultiVector::elementWiseMultiply: A must be a BlockedVector.");

  for (size_t m = 0; m < bmap->getNumMaps(); m++) {
    Teuchos::RCP<const Vector> partA      = bA->getMultiVector(m, bmap->getThyraMode())->getVector(0);
    Teuchos::RCP<const MultiVector> partB = bB->getMultiVector(m, bmap->getThyraMode());
    Teuchos::RCP<MultiVector> thisPart    = this->getMultiVector(m, bmap->getThyraMode());
    thisPart->elementWiseMultiply(scalarAB, *partA, *partB, scalarThis);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getNumVectors() const {
  return numVectors_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getLocalLength() const {
  return map_->getFullMap()->getLocalNumElements();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalLength() const {
  return map_->getFullMap()->getGlobalNumElements();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    isSameSize(const MultiVector& vec) const {
  const BlockedMultiVector* Vb = dynamic_cast<const BlockedMultiVector*>(&vec);
  if (!Vb)
    return false;
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    Teuchos::RCP<const MultiVector> a = getMultiVector(r);
    Teuchos::RCP<const MultiVector> b = Vb->getMultiVector(r);
    if ((a == Teuchos::null && b != Teuchos::null) || (a != Teuchos::null && b == Teuchos::null))
      return false;
    if (a != Teuchos::null && b != Teuchos::null && !a->isSameSize(*b))
      return false;
  }
  return true;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::string
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    description() const {
  return std::string("BlockedMultiVector");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const {
  out << description() << std::endl;
  for (size_t r = 0; r < map_->getNumMaps(); r++)
    getMultiVector(r)->describe(out, verbLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceMap(const Teuchos::RCP<const Map>& map) {
  // Non-blocked map: only valid if this multivector has a single block.
  if (this->getBlockedMap()->getNumMaps() > 1) {
    throw std::runtime_error(
        "BlockedMultiVector::replaceMap: map is not of type BlockedMap. General implementation not available, yet.");
    TEUCHOS_UNREACHABLE_RETURN();
  }
  std::vector<Teuchos::RCP<const Map>> subMaps(1, map);
  map_ = Teuchos::rcp(new BlockedMap(map, subMaps, this->getBlockedMap()->getThyraMode()));
  this->getMultiVector(0)->replaceMap(map);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceMap(const Teuchos::RCP<const BlockedMap>& bmap) {
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() != bmap->getThyraMode(),
                             std::runtime_error,
                             "BlockedMultiVector::replaceMap: inconsistent Thyra mode");
  map_ = bmap;
  for (size_t r = 0; r < map_->getNumMaps(); r++)
    getMultiVector(r)->replaceMap(bmap->getMap(r, map_->getThyraMode()));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    randomize(bool /* bUseXpetraImplementation */) {
  for (size_t r = 0; r < map_->getNumMaps(); ++r)
    getMultiVector(r)->randomize();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    randomize(const Scalar& minVal, const Scalar& maxVal, bool /* bUseXpetraImplementation */) {
  for (size_t r = 0; r < map_->getNumMaps(); ++r)
    getMultiVector(r)->randomize(minVal, maxVal);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Map>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMap() const {
  // BlockedMap is not a Tpetra::Map, so return the (real) full map.  Block-aware
  // callers should use getBlockedMap() instead.
  return map_->getFullMap();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const typename BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockedMap>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getBlockedMap() const {
  return map_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiVector>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMultiVector(size_t r) const {
  TEUCHOS_TEST_FOR_EXCEPTION(r > map_->getNumMaps(),
                             std::out_of_range,
                             "Error, r = " << r << " is too big. The BlockedMultiVector only contains " << map_->getNumMaps()
                                           << " partial blocks.");
  return vv_[r];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiVector>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMultiVector(size_t r, bool bThyraMode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(r > map_->getNumMaps(),
                             std::out_of_range,
                             "Error, r = " << r << " is too big. The BlockedMultiVector only contains " << map_->getNumMaps()
                                           << " partial blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() != bThyraMode,
                             std::runtime_error,
                             "BlockedMultiVector::getMultiVector: inconsistent Thyra mode");
  (void)bThyraMode;
  return vv_[r];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setMultiVector(size_t r,
                   Teuchos::RCP<const MultiVector> v,
                   bool bThyraMode) {
  TEUCHOS_TEST_FOR_EXCEPTION(r >= map_->getNumMaps(),
                             std::out_of_range,
                             "Error, r = " << r << " is too big. The BlockedMultiVector only contains " << map_->getNumMaps() << " partial blocks.");
  // NOTE (Tpetra port): Tpetra::MultiVector::getNumVectors() is NOT virtual (unlike
  // Xpetra's), so a blocked partial vector passed here through an RCP<const MultiVector>
  // would report the empty base object's column count.  Query the true count through the
  // blocked interface when the incoming vector is itself a BlockedMultiVector.
  Teuchos::RCP<const BlockedMultiVector> vblocked = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(v);
  const size_t vNumVectors                        = vblocked.is_null() ? v->getNumVectors() : vblocked->getNumVectors();
  TEUCHOS_TEST_FOR_EXCEPTION(numVectors_ != vNumVectors,
                             std::runtime_error,
                             "The BlockedMultiVectors expects " << getNumVectors() << " vectors. The provided partial multivector has "
                                                                << vNumVectors << " vectors.");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() != bThyraMode,
                             std::runtime_error,
                             "BlockedMultiVector::setMultiVector: inconsistent Thyra mode");
  (void)bThyraMode;
  Teuchos::RCP<MultiVector> vv = Teuchos::rcp_const_cast<MultiVector>(v);
  TEUCHOS_TEST_FOR_EXCEPTION(vv == Teuchos::null, std::runtime_error, "Partial vector must not be Teuchos::null");
  vv_[r] = vv;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiVector>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Merge() const {
  Teuchos::RCP<MultiVector> v = Teuchos::rcp(new MultiVector(map_->getFullMap(), getNumVectors()));
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    Teuchos::RCP<MultiVector> vi         = getMultiVector(r);
    Teuchos::RCP<BlockedMultiVector> bvi = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(vi);
    if (bvi.is_null() == true) {
      this->InsertVector(vi, r, v, map_->getThyraMode());
    } else {
      Teuchos::RCP<MultiVector> mvi = bvi->Merge();
      this->InsertVector(mvi, r, v, map_->getThyraMode());
    }
  }
  return v;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    assign(const MultiVector& rhs) {
  Teuchos::RCP<const MultiVector> rcpRhs      = Teuchos::rcpFromRef(rhs);
  Teuchos::RCP<const BlockedMultiVector> bRhs = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpRhs);
  if (bRhs == Teuchos::null)
    throw std::runtime_error("BlockedMultiVector::assign: argument is not a blocked multi vector.");

  map_ = Teuchos::rcp(new BlockedMap(*(bRhs->getBlockedMap())));
  vv_  = std::vector<Teuchos::RCP<MultiVector>>(map_->getNumMaps());
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    Teuchos::RCP<MultiVector> src = bRhs->getMultiVector(r, map_->getThyraMode());

    Teuchos::RCP<MultiVector> vv = Teuchos::rcp(new MultiVector(
        map_->getMap(r, bRhs->getBlockedMap()->getThyraMode()), rcpRhs->getNumVectors(), true));

    Teuchos::RCP<BlockedMultiVector> bsrc = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(src);
    Teuchos::RCP<BlockedMultiVector> bvv  = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(vv);
    if (bsrc.is_null() == true && bvv.is_null() == true) {
      *vv = *src;  // deep copy
    } else if (bsrc.is_null() == true && bvv.is_null() == false) {
      *vv = *src;  // deep copy
    } else if (bsrc.is_null() == false && bvv.is_null() == true) {
      if (bsrc->getBlockedMap()->getNumMaps() > 1) {
        throw std::runtime_error(
            "BlockedMultiVector::assign: source vector is of type BlockedMultiVector (with more than "
            "1 blocks) and target is a MultiVector.");
        TEUCHOS_UNREACHABLE_RETURN();
      }
      Teuchos::RCP<MultiVector> ssrc = bsrc->getMultiVector(0, map_->getThyraMode());
      *vv = *ssrc;
    } else {
      *vv = *src;  // deep copy
    }
    vv_[r] = vv;
  }
  numVectors_ = rcpRhs->getNumVectors();
}

// ---------------------------------------------------------------------------
// private ExtractVector/InsertVector helpers
// ---------------------------------------------------------------------------

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(Teuchos::RCP<const MultiVector>& full, size_t block, Teuchos::RCP<MultiVector>& partial) const {
  ExtractVector(*full, block, *partial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(Teuchos::RCP<MultiVector>& full, size_t block, Teuchos::RCP<MultiVector>& partial) const {
  ExtractVector(*full, block, *partial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiVector>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(Teuchos::RCP<const MultiVector>& full, size_t block, bool bThyraMode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                             std::out_of_range,
                             "ExtractVector: Error, block = " << block << " is too big. The BlockedMap only contains " << map_->getNumMaps()
                                                              << " partial blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == Teuchos::null, std::runtime_error, "ExtractVector: map_->getMap(" << block << ",false) is null");
  Teuchos::RCP<const BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(full);
  if (bfull.is_null() == true) {
    const Teuchos::RCP<MultiVector> vv =
        Teuchos::rcp(new MultiVector(map_->getMap(block, false), full->getNumVectors(), false));
    Teuchos::RCP<const Map> oldThyMapFull     = full->getMap();  // temporarily store map of full
    Teuchos::RCP<MultiVector> rcpNonConstFull = Teuchos::rcp_const_cast<MultiVector>(full);
    rcpNonConstFull->replaceMap(map_->getImporter(block)->getSourceMap());
    ExtractVector(*rcpNonConstFull, block, *vv);
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true,
                               std::runtime_error,
                               "BlockedMultiVector::ExtractVector: ExtractVector in Thyra-style numbering only possible if BlockedMap has been "
                               "created using Thyra-style numbered submaps.");
    if (bThyraMode == true)
      vv->replaceMap(map_->getMap(block, true));  // switch to Thyra-style map
    rcpNonConstFull->replaceMap(oldThyMapFull);
    return vv;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getNumMaps() != bfull->getBlockedMap()->getNumMaps(),
                               std::runtime_error,
                               "ExtractVector: Number of blocks in BlockedMap is " << map_->getNumMaps() << " but should be "
                                                                                   << bfull->getBlockedMap()->getNumMaps()
                                                                                   << " (number of blocks in BlockedMultiVector)");
    return bfull->getMultiVector(block, bThyraMode);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<typename BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::MultiVector>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(Teuchos::RCP<MultiVector>& full, size_t block, bool bThyraMode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                             std::out_of_range,
                             "ExtractVector: Error, block = " << block << " is too big. The BlockedMap only contains " << map_->getNumMaps()
                                                              << " partial blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == Teuchos::null, std::runtime_error, "ExtractVector: map_->getMap(" << block << ",false) is null");
  Teuchos::RCP<BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(full);
  if (bfull.is_null() == true) {
    const Teuchos::RCP<MultiVector> vv =
        Teuchos::rcp(new MultiVector(map_->getMap(block, false), full->getNumVectors(), false));
    Teuchos::RCP<const Map> oldThyMapFull = full->getMap();  // temporarily store map of full
    full->replaceMap(map_->getImporter(block)->getSourceMap());
    ExtractVector(*full, block, *vv);
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true,
                               std::runtime_error,
                               "BlockedMultiVector::ExtractVector: ExtractVector in Thyra-style numbering only possible if BlockedMap has been "
                               "created using Thyra-style numbered submaps.");
    if (bThyraMode == true)
      vv->replaceMap(map_->getMap(block, true));  // switch to Thyra-style map
    full->replaceMap(oldThyMapFull);
    return vv;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getNumMaps() != bfull->getBlockedMap()->getNumMaps(),
                               std::runtime_error,
                               "ExtractVector: Number of blocks in BlockedMap is " << map_->getNumMaps() << " but should be "
                                                                                   << bfull->getBlockedMap()->getNumMaps()
                                                                                   << " (number of blocks in BlockedMultiVector)");
    return bfull->getMultiVector(block, bThyraMode);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(const MultiVector& full, size_t block, MultiVector& partial) const {
  TEUCHOS_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                             std::out_of_range,
                             "ExtractVector: Error, block = " << block << " is too big. The BlockedMultiVector only contains " << map_->getNumMaps()
                                                              << " partial blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block) == Teuchos::null, std::runtime_error, "ExtractVector: maps_[" << block << "] is null");
  partial.doImport(full, *(map_->getImporter(block)), Tpetra::INSERT);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(const MultiVector& partial, size_t block, MultiVector& full, bool bThyraMode) const {
  TEUCHOS_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                             std::out_of_range,
                             "InsertVector: Error, block = " << block << " is too big. The BlockedMap only contains " << map_->getNumMaps()
                                                             << " partial blocks.");
  TEUCHOS_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == Teuchos::null, std::runtime_error, "InsertVector: map_->getMap(" << block << ",false) is null");
  TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true,
                             std::runtime_error,
                             "BlockedMultiVector::InsertVector: InsertVector in Thyra-style numbering only possible if BlockedMap has been created "
                             "using Thyra-style numbered submaps.");
  if (bThyraMode) {
    Teuchos::RCP<const MultiVector> rcpPartial   = Teuchos::rcpFromRef(partial);
    Teuchos::RCP<MultiVector> rcpNonConstPartial = Teuchos::rcp_const_cast<MultiVector>(rcpPartial);
    Teuchos::RCP<const Map> oldThyMapPartial     = rcpNonConstPartial->getMap();  // temporarily store map of partial
    Teuchos::RCP<const Map> oldThyMapFull        = full.getMap();                 // temporarily store map of full

    rcpNonConstPartial->replaceMap(map_->getMap(block, false));  // temporarily switch to xpetra-style map
    full.replaceMap(map_->getImporter(block)->getSourceMap());   // temporarily switch to Xpetra GIDs

    full.doExport(*rcpNonConstPartial, *(map_->getImporter(block)), Tpetra::INSERT);

    full.replaceMap(oldThyMapFull);                    // reset original map (Thyra GIDs)
    rcpNonConstPartial->replaceMap(oldThyMapPartial);  // change map back to original map
  } else {
    full.doExport(partial, *(map_->getImporter(block)), Tpetra::INSERT);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(Teuchos::RCP<const MultiVector> partial, size_t block, Teuchos::RCP<MultiVector> full, bool bThyraMode) const {
  Teuchos::RCP<BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(full);
  if (bfull.is_null() == true)
    InsertVector(*partial, block, *full, bThyraMode);
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block) == Teuchos::null, std::runtime_error, "InsertVector: maps_[" << block << "] is null");
    bfull->setMultiVector(block, partial, bThyraMode);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(Teuchos::RCP<MultiVector> partial, size_t block, Teuchos::RCP<MultiVector> full, bool bThyraMode) const {
  Teuchos::RCP<BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(full);
  if (bfull.is_null() == true)
    InsertVector(*partial, block, *full, bThyraMode);
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getMap(block) == Teuchos::null, std::runtime_error, "InsertVector: maps_[" << block << "] is null");
    bfull->setMultiVector(block, partial, bThyraMode);
  }
}

}  // namespace Tpetra

//
// Explicit instantiation macro
//
#define TPETRA_BLOCKEDMULTIVECTOR_INSTANT(SCALAR, LO, GO, NODE) \
  template class BlockedMultiVector<SCALAR, LO, GO, NODE>;

#endif  // TPETRA_BLOCKEDMULTIVECTOR_DEF_HPP
