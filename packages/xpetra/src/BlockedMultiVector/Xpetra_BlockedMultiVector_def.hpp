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

#include "Xpetra_MultiVectorFactory.hpp"
#include "Xpetra_BlockedVector.hpp"
#include "Xpetra_MapExtractor.hpp"

namespace Xpetra {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMultiVector(const Teuchos::RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>& map,
                       size_t NumVectors,
                       bool zeroOut)
  : map_(map) {
  numVectors_ = NumVectors;

  vv_.reserve(map->getNumMaps());

  // add CrsMatrix objects in row,column order
  for (size_t r = 0; r < map->getNumMaps(); ++r)
    vv_.push_back(
        Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(map->getMap(r, map_->getThyraMode()), NumVectors, zeroOut));
}

/*!
 * Const version of constructor which accepts a const version
 * of the multi-vector
 *
 * \note If you change the information in input vector v the data in the
 *       blocked multi-vector are not affected (and vice versa). Consider
 *       the blocked multivector to be a copy of the input multivector (not a view)
 *
 * \param bmap BlockedMap object containing information about the block splitting
 * \param v MultiVector that is to be splitted into a blocked multi vector
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMultiVector(Teuchos::RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap,
                       Teuchos::RCP<const MultiVector> v) {
  XPETRA_TEST_FOR_EXCEPTION(bmap->getMap()->getLocalNumElements() != v->getMap()->getLocalNumElements(),
                            Xpetra::Exceptions::RuntimeError,
                            "BlockedMultiVector: inconsistent number of local elements of MultiVector and BlockedMap. The BlockedMap has "
                                << bmap->getMap()->getLocalNumElements() << " local elements. The vector has " << v->getMap()->getLocalNumElements()
                                << ".");
  XPETRA_TEST_FOR_EXCEPTION(bmap->getMap()->getGlobalNumElements() != v->getMap()->getGlobalNumElements(),
                            Xpetra::Exceptions::RuntimeError,
                            "BlockedMultiVector: inconsistent number of global elements of MultiVector and BlockedMap. The BlockedMap has "
                                << bmap->getMap()->getGlobalNumElements() << " local elements. The vector has " << v->getMap()->getGlobalNumElements()
                                << ".");
  // TEUCHOS_TEST_FOR_EXCEPTION(bmap->getFullMap()->getLocalNumElements() != v->getMap()->getLocalNumElements(), Xpetra::Exceptions::RuntimeError,
  // "BlockedMultiVector: inconsistent number of local elements of MultiVector and BlockedMap. The BlockedMap has " <<
  // bmap->getFullMap()->getLocalNumElements() << " local elements. The vector has " << v->getMap()->getLocalNumElements() << ".");
  // TEUCHOS_TEST_FOR_EXCEPTION(bmap->getFullMap()->getGlobalNumElements() != v->getMap()->getGlobalNumElements(), Xpetra::Exceptions::RuntimeError,
  // "BlockedMultiVector: inconsistent number of global elements of MultiVector and BlockedMap. The BlockedMap has " <<
  // bmap->getFullMap()->getGlobalNumElements() << " local elements. The vector has " << v->getMap()->getGlobalNumElements() << ".");

  numVectors_ = v->getNumVectors();

  map_ = bmap;

  // Extract vector
  vv_.reserve(bmap->getNumMaps());

  // add CrsMatrix objects in row,column order
  for (size_t r = 0; r < bmap->getNumMaps(); ++r)
    vv_.push_back(this->ExtractVector(v, r, bmap->getThyraMode()));
}

/*!
 * Const version of constructor which accepts a const version
 * of the multi-vector
 *
 * \note If you change the information in input vector v the data in the
 *       blocked multi-vector are not affected (and vice versa). Consider
 *       the blocked multivector to be a copy of the input multivector (not a view)
 *
 * \param mapExtractor MapExtractor object containing information about the block splitting
 * \param v MultiVector that is to be splitted into a blocked multi vector
 */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMultiVector(Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mapExtractor,
                       Teuchos::RCP<const MultiVector> v) {
  numVectors_ = v->getNumVectors();

  // create blocked map out of MapExtractor
  std::vector<RCP<const Map>> maps;
  maps.reserve(mapExtractor->NumMaps());
  for (size_t r = 0; r < mapExtractor->NumMaps(); ++r)
    maps.push_back(mapExtractor->getMap(r, mapExtractor->getThyraMode()));
  map_ = Teuchos::rcp(new Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>(mapExtractor->getFullMap(), maps, mapExtractor->getThyraMode()));

  // Extract vector
  vv_.reserve(mapExtractor->NumMaps());

  // add CrsMatrix objects in row,column order
  for (size_t r = 0; r < mapExtractor->NumMaps(); ++r)
    vv_.push_back(this->ExtractVector(v, r, mapExtractor->getThyraMode()));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    BlockedMultiVector(const Teuchos::RCP<const BlockedMap>& map,
                       std::vector<Teuchos::RCP<MultiVector>>& vin) {
  numVectors_ = vin[0]->getNumVectors();
  map_        = map;
  vv_.resize(vin.size());
  for (size_t i = 0; i < vv_.size(); i++)
    vv_[i] = vin[i];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ~BlockedMultiVector() {
  for (size_t r = 0; r < vv_.size(); ++r) {
    vv_[r] = Teuchos::null;
  }

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
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::replaceGlobalValue: Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoGlobalValue(GlobalOrdinal /* globalRow */, size_t /* vectorIndex */, const Scalar& /* value */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::sumIntoGlobalValue: Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    replaceLocalValue(LocalOrdinal /* myRow */, size_t /* vectorIndex */, const Scalar& /* value */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::replaceLocalValue: Not supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    sumIntoLocalValue(LocalOrdinal /* myRow */, size_t /* vectorIndex */, const Scalar& /* value */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::sumIntoLocalValue:Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    putScalar(const Scalar& value) {
  XPETRA_MONITOR("BlockedMultiVector::putScalar");
  for (size_t r = 0; r < map_->getNumMaps(); r++) {
    getMultiVector(r)->putScalar(value);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVector(size_t j) const {
  RCP<Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> ret =
      Teuchos::rcp(new Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(this->getBlockedMap(), false));

  for (size_t r = 0; r < getBlockedMap()->getNumMaps(); r++) {
    RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> subvec =
        this->getMultiVector(r, this->getBlockedMap()->getThyraMode())->getVector(j);
    ret->setMultiVector(r, subvec, this->getBlockedMap()->getThyraMode());
  }
  return ret;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getVectorNonConst(size_t j) {
  RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> ret =
      Teuchos::rcp_const_cast<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(getVector(j));
  return ret;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<const Scalar>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getData(size_t j) const {
  if (map_->getNumMaps() == 1) {
    return vv_[0]->getData(j);
  }
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::getData: Not (yet) supported by BlockedMultiVector.");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<Scalar>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getDataNonConst(size_t j) {
  if (map_->getNumMaps() == 1) {
    return vv_[0]->getDataNonConst(j);
  }
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::getDataNonConst: Not (yet) supported by BlockedMultiVector.");
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    dot(const MultiVector& /* A */, const Teuchos::ArrayView<Scalar>& /* dots */) const {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::dot: Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    abs(const MultiVector& /* A */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::abs: Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    reciprocal(const MultiVector& /* A */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::reciprocal: Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    scale(const Scalar& alpha) {
  XPETRA_MONITOR("BlockedMultiVector::scale (Scalar)");
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    if (getMultiVector(r) != Teuchos::null) {
      getMultiVector(r)->scale(alpha);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    scale(Teuchos::ArrayView<const Scalar> alpha) {
  XPETRA_MONITOR("BlockedMultiVector::scale (Array)");
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    if (getMultiVector(r) != Teuchos::null) {
      getMultiVector(r)->scale(alpha);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    update(const Scalar& alpha, const MultiVector& A, const Scalar& beta) {
  XPETRA_MONITOR("BlockedMultiVector::update");
  Teuchos::RCP<const MultiVector> rcpA      = Teuchos::rcpFromRef(A);
  Teuchos::RCP<const BlockedMultiVector> bA = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpA);
  TEUCHOS_TEST_FOR_EXCEPTION(numVectors_ != rcpA->getNumVectors(),
                             Xpetra::Exceptions::RuntimeError,
                             "BlockedMultiVector::update: update with incompatible vector (different number of vectors in multivector).");
  if (bA != Teuchos::null) {
    // A is a BlockedMultiVector (and compatible with this)
    // Call update recursively on all sub vectors
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() != bA->getBlockedMap()->getThyraMode(),
                               Xpetra::Exceptions::RuntimeError,
                               "BlockedMultiVector::update: update with incompatible vector (different thyra mode).");
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getNumMaps() != bA->getBlockedMap()->getNumMaps(),
                               Xpetra::Exceptions::RuntimeError,
                               "BlockedMultiVector::update: update with incompatible vector (different number of partial vectors).");
    for (size_t r = 0; r < map_->getNumMaps(); r++) {
      XPETRA_TEST_FOR_EXCEPTION(getMultiVector(r)->getMap()->getLocalNumElements() != bA->getMultiVector(r)->getMap()->getLocalNumElements(),
                                Xpetra::Exceptions::RuntimeError,
                                "BlockedMultiVector::update: in subvector "
                                    << r << ": Cannot add a vector of (local) length " << bA->getMultiVector(r)->getMap()->getLocalNumElements()
                                    << " to the existing vector with " << getMultiVector(r)->getMap()->getLocalNumElements() << " entries.");
      XPETRA_TEST_FOR_EXCEPTION(getMultiVector(r)->getMap()->getGlobalNumElements() != bA->getMultiVector(r)->getMap()->getGlobalNumElements(),
                                Xpetra::Exceptions::RuntimeError,
                                "BlockedMultiVector::update: in subvector "
                                    << r << ": Cannot add a vector of length " << bA->getMultiVector(r)->getMap()->getGlobalNumElements()
                                    << " to the existing vector with " << getMultiVector(r)->getMap()->getGlobalNumElements() << " entries.");

      // TAW: 12/6 We basically want to do something like:
      //   getMultiVector(r)->update(alpha, *(bA->getMultiVector(r)), beta);
      // But if the left hand side is a standard MultiVector and bA->getMultiVector(r) is
      // a blocked MultiVector this will not work, as the update implementation of the standard
      // multivector cannot deal with blocked multivectors.
      // The only use case where this could happen is, if the rhs vector is a single block multivector
      Teuchos::RCP<MultiVector> lmv = getMultiVector(r);
      Teuchos::RCP<MultiVector> rmv = bA->getMultiVector(r);

      // check whether lmv/rmv are blocked or not
      Teuchos::RCP<BlockedMultiVector> blmv = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(lmv);
      Teuchos::RCP<BlockedMultiVector> brmv = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(rmv);

      if (blmv.is_null() == true && brmv.is_null() == false) {
        // special case: lmv is standard MultiVector but rmv is BlockedMultiVector
        TEUCHOS_TEST_FOR_EXCEPTION(brmv->getBlockedMap()->getNumMaps() > 1,
                                   Xpetra::Exceptions::RuntimeError,
                                   "BlockedMultiVector::update: Standard MultiVector object does not accept BlockedMultVector object as "
                                   "parameter in update call.");
        lmv->update(alpha, *(brmv->getMultiVector(0)), beta);
      } else
        lmv->update(alpha, *rmv, beta);
    }
  } else {
    // A is a MultiVector
    // If this is a BlockedMultiVector with only one sub-vector of same length we can just update
    // Otherwise, A is not compatible with this as BlockedMultiVector and we have to extract the vector data from A

    if (getBlockedMap()->getNumMaps() == 1) {
      // Actually call "update" on the underlying MultiVector (Epetra or Tpetra).
      // The maps have to be compatible.
      XPETRA_TEST_FOR_EXCEPTION(
          getMultiVector(0)->getMap()->isSameAs(*(rcpA->getMap())) == false,
          Xpetra::Exceptions::RuntimeError,
          "BlockedMultiVector::update: update with incompatible vector (maps of full vector do not match with map in MapExtractor).");
      getMultiVector(0)->update(alpha, *rcpA, beta);
    } else {
      // general case: A has to be splitted and subparts have to be extracted and stored in this BlockedMultiVector
      // XPETRA_TEST_FOR_EXCEPTION(map_->getFullMap()->isSameAs(*(rcpA->getMap()))==false, Xpetra::Exceptions::RuntimeError,
      // "BlockedMultiVector::update: update with incompatible vector (maps of full vector do not match with map in MapExtractor). - Note:
      // This test might be too strict and can probably be relaxed!");
      for (size_t r = 0; r < map_->getNumMaps(); r++) {
        // Call "update" on the subvector. Note, that getMultiVector(r) could return another BlockedMultiVector.
        // That is, in Thyra mode the maps could differ (local Xpetra versus Thyra style gids)
        Teuchos::RCP<const MultiVector> part = this->ExtractVector(rcpA, r, map_->getThyraMode());
        XPETRA_TEST_FOR_EXCEPTION(getMultiVector(r)->getMap()->getLocalNumElements() != part->getMap()->getLocalNumElements(),
                                  Xpetra::Exceptions::RuntimeError,
                                  "BlockedMultiVector::update: in subvector "
                                      << r << ": Cannot add a vector of (local) length " << part->getMap()->getLocalNumElements()
                                      << " to the existing vector with " << getMultiVector(r)->getMap()->getLocalNumElements() << " entries.");
        XPETRA_TEST_FOR_EXCEPTION(getMultiVector(r)->getMap()->getGlobalNumElements() != part->getMap()->getGlobalNumElements(),
                                  Xpetra::Exceptions::RuntimeError,
                                  "BlockedMultiVector::update: in subvector "
                                      << r << ": Cannot add a vector of length " << part->getMap()->getGlobalNumElements()
                                      << " to the existing vector with " << getMultiVector(r)->getMap()->getGlobalNumElements() << " entries.");
        getMultiVector(r)->update(alpha, *part, beta);
      }
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    update(const Scalar& alpha, const MultiVector& A, const Scalar& beta, const MultiVector& B, const Scalar& gamma) {
  XPETRA_MONITOR("BlockedMultiVector::update2");
  Teuchos::RCP<const MultiVector> rcpA      = Teuchos::rcpFromRef(A);
  Teuchos::RCP<const MultiVector> rcpB      = Teuchos::rcpFromRef(B);
  Teuchos::RCP<const BlockedMultiVector> bA = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpA);
  Teuchos::RCP<const BlockedMultiVector> bB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpB);
  if (bA != Teuchos::null && bB != Teuchos::null) {
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() != bA->getBlockedMap()->getThyraMode(),
                               Xpetra::Exceptions::RuntimeError,
                               "BlockedMultiVector::update: update with incompatible vector (different thyra mode in vector A).");
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getNumMaps() != bA->getBlockedMap()->getNumMaps(),
                               Xpetra::Exceptions::RuntimeError,
                               "BlockedMultiVector::update: update with incompatible vector (different number of partial vectors in vector A).");
    TEUCHOS_TEST_FOR_EXCEPTION(
        numVectors_ != bA->getNumVectors(),
        Xpetra::Exceptions::RuntimeError,
        "BlockedMultiVector::update: update with incompatible vector (different number of vectors in multivector in vector A).");
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() != bB->getBlockedMap()->getThyraMode(),
                               Xpetra::Exceptions::RuntimeError,
                               "BlockedMultiVector::update: update with incompatible vector (different thyra mode in vector B).");
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getNumMaps() != bB->getBlockedMap()->getNumMaps(),
                               Xpetra::Exceptions::RuntimeError,
                               "BlockedMultiVector::update: update with incompatible vector (different number of partial vectors in vector B).");
    TEUCHOS_TEST_FOR_EXCEPTION(
        numVectors_ != bB->getNumVectors(),
        Xpetra::Exceptions::RuntimeError,
        "BlockedMultiVector::update: update with incompatible vector (different number of vectors in multivector in vector B).");

    for (size_t r = 0; r < map_->getNumMaps(); r++) {
      XPETRA_TEST_FOR_EXCEPTION(getMultiVector(r)->getMap()->isSameAs(*(bA->getMultiVector(r)->getMap())) == false,
                                Xpetra::Exceptions::RuntimeError,
                                "BlockedMultiVector::update: update with incompatible vector (different maps in partial vector " << r << ").");
      getMultiVector(r)->update(alpha, *(bA->getMultiVector(r)), beta, *(bB->getMultiVector(r)), gamma);
    }
    return;
  }
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::update: only supports update with other BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm1(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const {
  XPETRA_MONITOR("BlockedMultiVector::norm1");
  typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
  Array<Magnitude> temp_norms(getNumVectors());
  std::fill(norms.begin(), norms.end(), ScalarTraits<Magnitude>::zero());
  std::fill(temp_norms.begin(), temp_norms.end(), ScalarTraits<Magnitude>::zero());
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    if (getMultiVector(r) != Teuchos::null) {
      getMultiVector(r)->norm1(temp_norms);
      for (size_t c = 0; c < getNumVectors(); ++c)
        norms[c] += temp_norms[c];
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const {
  XPETRA_MONITOR("BlockedMultiVector::norm2");
  typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
  Array<Magnitude> results(getNumVectors());
  Array<Magnitude> temp_norms(getNumVectors());
  std::fill(norms.begin(), norms.end(), ScalarTraits<Magnitude>::zero());
  std::fill(results.begin(), results.end(), ScalarTraits<Magnitude>::zero());
  std::fill(temp_norms.begin(), temp_norms.end(), ScalarTraits<Magnitude>::zero());
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    if (getMultiVector(r) != Teuchos::null) {
      getMultiVector(r)->norm2(temp_norms);
      for (size_t c = 0; c < getNumVectors(); ++c)
        results[c] += temp_norms[c] * temp_norms[c];
    }
  }
  for (size_t c = 0; c < getNumVectors(); ++c)
    norms[c] = Teuchos::ScalarTraits<Magnitude>::squareroot(results[c]);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const {
  XPETRA_MONITOR("BlockedMultiVector::normInf");
  typedef typename ScalarTraits<Scalar>::magnitudeType Magnitude;
  Array<Magnitude> temp_norms(getNumVectors());
  std::fill(norms.begin(), norms.end(), ScalarTraits<Magnitude>::zero());
  std::fill(temp_norms.begin(), temp_norms.end(), ScalarTraits<Magnitude>::zero());
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    if (getMultiVector(r) != Teuchos::null) {
      getMultiVector(r)->normInf(temp_norms);
      for (size_t c = 0; c < getNumVectors(); ++c)
        norms[c] = std::max(norms[c], temp_norms[c]);
    }
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    meanValue(const Teuchos::ArrayView<Scalar>& /* means */) const {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::meanValue: Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    multiply(Teuchos::ETransp /* transA */,
             Teuchos::ETransp /* transB */,
             const Scalar& /* alpha  */,
             const MultiVector& /* A      */,
             const MultiVector& /* B      */,
             const Scalar& /* beta   */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::multiply: Not (yet) supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    elementWiseMultiply(Scalar scalarAB,
                        const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A,
                        const MultiVector& B,
                        Scalar scalarThis) {
  XPETRA_MONITOR("BlockedMultiVector::elementWiseMultiply");
  XPETRA_TEST_FOR_EXCEPTION(B.getMap()->isSameAs(*(this->getMap())) == false,
                            Xpetra::Exceptions::RuntimeError,
                            "BlockedMultiVector::elementWiseMultipy: B must have same blocked map than this.");
  // XPETRA_TEST_FOR_EXCEPTION(A.getMap()->isSameAs(*(this->getMap()))==false, Xpetra::Exceptions::RuntimeError,
  // "BlockedMultiVector::elementWiseMultipy: A must have same blocked map than this.");
  TEUCHOS_TEST_FOR_EXCEPTION(A.getMap()->getLocalNumElements() != B.getMap()->getLocalNumElements(),
                             Xpetra::Exceptions::RuntimeError,
                             "BlockedMultiVector::elementWiseMultipy: A has " << A.getMap()->getLocalNumElements() << " elements, B has "
                                                                              << B.getMap()->getLocalNumElements() << ".");
  TEUCHOS_TEST_FOR_EXCEPTION(A.getMap()->getGlobalNumElements() != B.getMap()->getGlobalNumElements(),
                             Xpetra::Exceptions::RuntimeError,
                             "BlockedMultiVector::elementWiseMultipy: A has " << A.getMap()->getGlobalNumElements() << " elements, B has "
                                                                              << B.getMap()->getGlobalNumElements() << ".");

  RCP<const BlockedMap> bmap                                                = getBlockedMap();
  RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> rcpA = Teuchos::rcpFromRef(A);
  RCP<const Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bA =
      Teuchos::rcp_dynamic_cast<const Xpetra::BlockedVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(rcpA);

  RCP<const MultiVector> rcpB      = Teuchos::rcpFromRef(B);
  RCP<const BlockedMultiVector> bB = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpB);
  TEUCHOS_TEST_FOR_EXCEPTION(
      bB.is_null() == true, Xpetra::Exceptions::RuntimeError, "BlockedMultiVector::elementWiseMultipy: B must be a BlockedMultiVector.");

  // RCP<Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node> > me = Teuchos::rcp(new
  // Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node>(bmap));

  for (size_t m = 0; m < bmap->getNumMaps(); m++) {
    RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> partA      = bA->getMultiVector(m, bmap->getThyraMode())->getVector(0);
    RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> partB = bB->getMultiVector(m, bmap->getThyraMode());
    RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> thisPart    = this->getMultiVector(m, bmap->getThyraMode());

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
  XPETRA_MONITOR("BlockedMultiVector::getLocalLength()");
  return map_->getFullMap()->getLocalNumElements();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
global_size_t
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getGlobalLength() const {
  XPETRA_MONITOR("BlockedMultiVector::getGlobalLength()");
  return map_->getFullMap()->getGlobalNumElements();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    isSameSize(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec) const {
  const BlockedMultiVector* Vb = dynamic_cast<const BlockedMultiVector*>(&vec);
  if (!Vb)
    return false;
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    RCP<const MultiVector> a = getMultiVector(r);
    RCP<const MultiVector> b = Vb->getMultiVector(r);
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
    replaceMap(const RCP<const Map>& map) {
  XPETRA_MONITOR("BlockedMultiVector::replaceMap");
  RCP<const BlockedMap> bmap = Teuchos::rcp_dynamic_cast<const BlockedMap>(map);
  if (bmap.is_null() == true) {
    // if this has more than 1 sub blocks but "map" is not a blocked map, they are very likely not compatible
    if (this->getBlockedMap()->getNumMaps() > 1) {
      throw Xpetra::Exceptions::RuntimeError(
          "BlockedMultiVector::replaceMap: map is not of type BlockedMap. General implementation not available, yet.");
      TEUCHOS_UNREACHABLE_RETURN();
    }
    // special case: this is a blocked map with only one block
    // TODO add more debug check (especially for Thyra mode)
    std::vector<Teuchos::RCP<const Map>> subMaps(1, map);
    map_ = Teuchos::rcp(new BlockedMap(map, subMaps, this->getBlockedMap()->getThyraMode()));
    this->getMultiVector(0)->replaceMap(map);
    return;
  }
  RCP<const BlockedMap> mybmap = Teuchos::rcp_dynamic_cast<const BlockedMap>(map_);

  XPETRA_TEST_FOR_EXCEPTION(
      mybmap->getThyraMode() != bmap->getThyraMode(), Xpetra::Exceptions::RuntimeError, "BlockedMultiVector::replaceMap: inconsistent Thyra mode");
  map_ = bmap;
  for (size_t r = 0; r < map_->getNumMaps(); r++)
    getMultiVector(r)->replaceMap(bmap->getMap(r, map_->getThyraMode()));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doImport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* source */,
             const Import& /* importer */,
             CombineMode /* CM */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::doImport: Not supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doExport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* dest */,
             const Import& /* importer */,
             CombineMode /* CM */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::doExport: Not supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doImport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* source */,
             const Export& /* exporter */,
             CombineMode /* CM */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::doImport: Not supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    doExport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& /* dest */,
             const Export& /* exporter */,
             CombineMode /* CM */) {
  throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::doExport: Not supported by BlockedMultiVector.");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setSeed(unsigned int seed) {
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    getMultiVector(r)->setSeed(seed);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    randomize(bool bUseXpetraImplementation) {
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    getMultiVector(r)->randomize(bUseXpetraImplementation);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    randomize(const Scalar& minVal, const Scalar& maxVal, bool bUseXpetraImplementation) {
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    getMultiVector(r)->randomize(minVal, maxVal, bUseXpetraImplementation);
  }
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

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMap() const {
  XPETRA_MONITOR("BlockedMultiVector::getMap");
  return map_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getBlockedMap() const {
  XPETRA_MONITOR("BlockedMultiVector::getBlockedMap");
  return map_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMultiVector(size_t r) const {
  XPETRA_MONITOR("BlockedMultiVector::getMultiVector(r)");
  TEUCHOS_TEST_FOR_EXCEPTION(r > map_->getNumMaps(),
                             std::out_of_range,
                             "Error, r = " << r << " is too big. The BlockedMultiVector only contains " << map_->getNumMaps()
                                           << " partial blocks.");
  return vv_[r];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    getMultiVector(size_t r, bool bThyraMode) const {
  XPETRA_MONITOR("BlockedMultiVector::getMultiVector(r,bThyraMode)");
  TEUCHOS_TEST_FOR_EXCEPTION(r > map_->getNumMaps(),
                             std::out_of_range,
                             "Error, r = " << r << " is too big. The BlockedMultiVector only contains " << map_->getNumMaps()
                                           << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getThyraMode() != bThyraMode, Xpetra::Exceptions::RuntimeError, "BlockedMultiVector::getMultiVector: inconsistent Thyra mode");
  (void)bThyraMode;  // avoid unused parameter warning when HAVE_XPETRA_DEBUG isn't defined
  return vv_[r];
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    setMultiVector(size_t r,
                   Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> v,
                   bool bThyraMode) {
  // The map of the MultiVector should be the same as the stored submap
  // In thyra mode the vectors should live on the thyra maps
  // in xpetra mode the should live in the xpetra maps
  // that should be also ok in the nested case for thyra (if the vectors are distributed accordingly)
  XPETRA_MONITOR("BlockedMultiVector::setMultiVector");
  XPETRA_TEST_FOR_EXCEPTION(r >= map_->getNumMaps(),
                            std::out_of_range,
                            "Error, r = " << r << " is too big. The BlockedMultiVector only contains " << map_->getNumMaps() << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(numVectors_ != v->getNumVectors(),
                            Xpetra::Exceptions::RuntimeError,
                            "The BlockedMultiVectors expects " << getNumVectors() << " vectors. The provided partial multivector has "
                                                               << v->getNumVectors() << " vectors.");
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getThyraMode() != bThyraMode, Xpetra::Exceptions::RuntimeError, "BlockedMultiVector::setMultiVector: inconsistent Thyra mode");
  (void)bThyraMode;  // avoid unused parameter warning when HAVE_XPETRA_DEBUG isn't defined
  Teuchos::RCP<MultiVector> vv = Teuchos::rcp_const_cast<MultiVector>(v);
  TEUCHOS_TEST_FOR_EXCEPTION(vv == Teuchos::null, Xpetra::Exceptions::RuntimeError, "Partial vector must not be Teuchos::null");
  vv_[r] = vv;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Merge() const {
  XPETRA_MONITOR("BlockedMultiVector::Merge");
  using Teuchos::RCP;

  RCP<MultiVector> v = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(map_->getFullMap(), getNumVectors());
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    RCP<MultiVector> vi         = getMultiVector(r);
    RCP<BlockedMultiVector> bvi = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(vi);
    if (bvi.is_null() == true) {
      this->InsertVector(vi, r, v, map_->getThyraMode());
    } else {
      RCP<MultiVector> mvi = bvi->Merge();
      this->InsertVector(mvi, r, v, map_->getThyraMode());
    }
  }

  // TODO plausibility checks

  return v;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    assign(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& rhs) {
  Teuchos::RCP<const MultiVector> rcpRhs      = Teuchos::rcpFromRef(rhs);
  Teuchos::RCP<const BlockedMultiVector> bRhs = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(rcpRhs);
  if (bRhs == Teuchos::null)
    throw Xpetra::Exceptions::RuntimeError("BlockedMultiVector::assign: argument is not a blocked multi vector.");

  map_ = Teuchos::rcp(new BlockedMap(*(bRhs->getBlockedMap())));
  vv_  = std::vector<Teuchos::RCP<MultiVector>>(map_->getNumMaps());
  for (size_t r = 0; r < map_->getNumMaps(); ++r) {
    // extract source vector (is of type MultiVector or BlockedMultiVector)
    RCP<MultiVector> src = bRhs->getMultiVector(r, map_->getThyraMode());

    // create new (empty) multivector (is of type MultiVector or BlockedMultiVector)
    RCP<MultiVector> vv = Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
        map_->getMap(r, bRhs->getBlockedMap()->getThyraMode()), rcpRhs->getNumVectors(), true);

    // check type of source and target vector
    RCP<BlockedMultiVector> bsrc = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(src);
    RCP<BlockedMultiVector> bvv  = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(vv);
    if (bsrc.is_null() == true && bvv.is_null() == true) {
      *vv = *src;  // deep copy
    } else if (bsrc.is_null() == true && bvv.is_null() == false) {
      // special case: source vector is a merged MultiVector but target vector is blocked
      *vv = *src;  // deep copy (is this a problem???)
    } else if (bsrc.is_null() == false && bvv.is_null() == true) {
      // special case: source vector is blocked but target vector is a merged MultiVector
      // This is a problem and only works if bsrc has only one block
      if (bsrc->getBlockedMap()->getNumMaps() > 1) {
        throw Xpetra::Exceptions::RuntimeError(
            "BlockedMultiVector::assign: source vector is of type BlockedMultiVector (with more than "
            "1 blocks) and target is a MultiVector.");
        TEUCHOS_UNREACHABLE_RETURN();
      }
      RCP<MultiVector> ssrc = bsrc->getMultiVector(0, map_->getThyraMode());
      XPETRA_TEST_FOR_EXCEPTION(ssrc.is_null() == true, Xpetra::Exceptions::RuntimeError, "BlockedMultiVector::assign: cannot extract vector");
      XPETRA_TEST_FOR_EXCEPTION(Teuchos::rcp_dynamic_cast<BlockedMultiVector>(ssrc) != Teuchos::null,
                                Xpetra::Exceptions::RuntimeError,
                                "BlockedMultiVector::assign: sub block must not be of type BlockedMultiVector.");
      *vv = *ssrc;
    } else {
      // this should work (no exception necessary, i guess)
      XPETRA_TEST_FOR_EXCEPTION(bsrc->getBlockedMap()->getNumMaps() != bvv->getBlockedMap()->getNumMaps(),
                                Xpetra::Exceptions::RuntimeError,
                                "BlockedMultiVector::assign: source and target are BlockedMultiVectors with a different number of submaps.");
      *vv = *src;  // deep copy
    }
    vv_[r] = vv;
  }
  numVectors_ = rcpRhs->getNumVectors();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& full,
                  size_t block,
                  RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& partial) const {
  ExtractVector(*full, block, *partial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& full,
                  size_t block,
                  RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& partial) const {
  ExtractVector(*full, block, *partial);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& full,
                  size_t block,
                  bool bThyraMode) const {
  XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                            std::out_of_range,
                            "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps()
                                                             << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "ExtractVector: map_->getmap(" << block << ",false) is null");
  RCP<const BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<const BlockedMultiVector>(full);
  if (bfull.is_null() == true) {
    // standard case: full is not of type BlockedMultiVector
    // first extract partial vector from full vector (using xpetra style GIDs)
    const RCP<MultiVector> vv =
        Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(map_->getMap(block, false), full->getNumVectors(), false);
    // if(bThyraMode == false) {
    //  ExtractVector(*full, block, *vv);
    //  return vv;
    //} else {
    RCP<const Map> oldThyMapFull     = full->getMap();  // temporarely store map of full
    RCP<MultiVector> rcpNonConstFull = Teuchos::rcp_const_cast<MultiVector>(full);
    rcpNonConstFull->replaceMap(map_->getImporter(block)->getSourceMap());
    ExtractVector(*rcpNonConstFull, block, *vv);
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true,
                               Xpetra::Exceptions::RuntimeError,
                               "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been "
                               "created using Thyra-style numbered submaps.");
    if (bThyraMode == true)
      vv->replaceMap(map_->getMap(block, true));  // switch to Thyra-style map
    rcpNonConstFull->replaceMap(oldThyMapFull);
    return vv;
    //}
  } else {
    // special case: full is of type BlockedMultiVector
    XPETRA_TEST_FOR_EXCEPTION(map_->getNumMaps() != bfull->getBlockedMap()->getNumMaps(),
                              Xpetra::Exceptions::RuntimeError,
                              "ExtractVector: Number of blocks in map extractor is " << map_->getNumMaps() << " but should be "
                                                                                     << bfull->getBlockedMap()->getNumMaps()
                                                                                     << " (number of blocks in BlockedMultiVector)");
    return bfull->getMultiVector(block, bThyraMode);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& full,
                  size_t block,
                  bool bThyraMode) const {
  XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                            std::out_of_range,
                            "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps()
                                                             << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "ExtractVector: map_->getmap(" << block << ",false) is null");
  RCP<BlockedMultiVector> bfull = Teuchos::rcp_dynamic_cast<BlockedMultiVector>(full);
  if (bfull.is_null() == true) {
    // standard case: full is not of type BlockedMultiVector
    // first extract partial vector from full vector (using xpetra style GIDs)
    const RCP<MultiVector> vv =
        Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(map_->getMap(block, false), full->getNumVectors(), false);
    // if(bThyraMode == false) {
    //  ExtractVector(*full, block, *vv);
    //  return vv;
    //} else {
    RCP<const Map> oldThyMapFull = full->getMap();  // temporarely store map of full
    full->replaceMap(map_->getImporter(block)->getSourceMap());
    ExtractVector(*full, block, *vv);
    TEUCHOS_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true,
                               Xpetra::Exceptions::RuntimeError,
                               "MapExtractor::ExtractVector: ExtractVector in Thyra-style numbering only possible if MapExtractor has been "
                               "created using Thyra-style numbered submaps.");
    if (bThyraMode == true)
      vv->replaceMap(map_->getMap(block, true));  // switch to Thyra-style map
    full->replaceMap(oldThyMapFull);
    return vv;
    //}
  } else {
    // special case: full is of type BlockedMultiVector
    XPETRA_TEST_FOR_EXCEPTION(map_->getNumMaps() != bfull->getBlockedMap()->getNumMaps(),
                              Xpetra::Exceptions::RuntimeError,
                              "ExtractVector: Number of blocks in map extractor is " << map_->getNumMaps() << " but should be "
                                                                                     << bfull->getBlockedMap()->getNumMaps()
                                                                                     << " (number of blocks in BlockedMultiVector)");
    return bfull->getMultiVector(block, bThyraMode);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    ExtractVector(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& full,
                  size_t block,
                  Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& partial) const {
  XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                            std::out_of_range,
                            "ExtractVector: Error, block = " << block << " is too big. The BlockedMultiVector only contains " << map_->getNumMaps()
                                                             << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block) == null, Xpetra::Exceptions::RuntimeError, "ExtractVector: maps_[" << block << "] is null");
  partial.doImport(full, *(map_->getImporter(block)), Xpetra::INSERT);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& partial,
                 size_t block,
                 Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& full,
                 bool bThyraMode) const {
  XPETRA_TEST_FOR_EXCEPTION(block >= map_->getNumMaps(),
                            std::out_of_range,
                            "ExtractVector: Error, block = " << block << " is too big. The MapExtractor only contains " << map_->getNumMaps()
                                                             << " partial blocks.");
  XPETRA_TEST_FOR_EXCEPTION(
      map_->getMap(block, false) == null, Xpetra::Exceptions::RuntimeError, "ExtractVector: map_->getmap(" << block << ",false) is null");
  XPETRA_TEST_FOR_EXCEPTION(map_->getThyraMode() == false && bThyraMode == true,
                            Xpetra::Exceptions::RuntimeError,
                            "MapExtractor::InsertVector: InsertVector in Thyra-style numbering only possible if MapExtractor has been created "
                            "using Thyra-style numbered submaps.");
  if (bThyraMode) {
    // NOTE: the importer objects in the BlockedMap are always using Xpetra GIDs (or Thyra style Xpetra GIDs)
    // The source map corresponds to the full map (in Xpetra GIDs) starting with GIDs from zero. The GIDs are consecutive in Thyra mode
    // The target map is the partial map (in the corresponding Xpetra GIDs)

    // TODO can we skip the Export call in special cases (i.e. Src = Target map, same length, etc...)

    // store original GIDs (could be Thyra GIDs)
    RCP<const MultiVector> rcpPartial   = Teuchos::rcpFromRef(partial);
    RCP<MultiVector> rcpNonConstPartial = Teuchos::rcp_const_cast<MultiVector>(rcpPartial);
    RCP<const Map> oldThyMapPartial     = rcpNonConstPartial->getMap();  // temporarely store map of partial
    RCP<const Map> oldThyMapFull        = full.getMap();                 // temporarely store map of full

    // check whether getMap(block,false) is identical to target map of importer
    // XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block,false)->isSameAs(*(map_->getImporter(block)->getTargetMap()))==false,
    // Xpetra::Exceptions::RuntimeError,
    //           "MapExtractor::InsertVector: InsertVector in Thyra-style mode: Xpetra GIDs of partial vector are not identical to target Map
    //           of Importer. This should not be.");

    // XPETRA_TEST_FOR_EXCEPTION(full.getMap()->isSameAs(*(map_->getImporter(block)->getSourceMap()))==false,
    // Xpetra::Exceptions::RuntimeError,
    //           "MapExtractor::InsertVector: InsertVector in Thyra-style mode: Xpetra GIDs of full vector are not identical to source Map of
    //           Importer. This should not be.");

    rcpNonConstPartial->replaceMap(map_->getMap(block, false));  // temporarely switch to xpetra-style map
    full.replaceMap(map_->getImporter(block)->getSourceMap());   // temporarely switch to Xpetra GIDs

    // do the Export
    full.doExport(*rcpNonConstPartial, *(map_->getImporter(block)), Xpetra::INSERT);

    // switch back to original maps
    full.replaceMap(oldThyMapFull);                    // reset original map (Thyra GIDs)
    rcpNonConstPartial->replaceMap(oldThyMapPartial);  // change map back to original map
  } else {
    // Xpetra style numbering
    full.doExport(partial, *(map_->getImporter(block)), Xpetra::INSERT);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> partial,
                 size_t block,
                 RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> full,
                 bool bThyraMode) const {
  RCP<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bfull =
      Teuchos::rcp_dynamic_cast<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(full);
  if (bfull.is_null() == true)
    InsertVector(*partial, block, *full, bThyraMode);
  else {
    XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block) == null, Xpetra::Exceptions::RuntimeError, "InsertVector: maps_[" << block << "] is null");

#if 0
            // WCMCLEN - ETI: MultiVector doesn't have a setMultiVector method,
            // WCMCLEN - ETI: but BlockedMultiVector does... should this be bfull->...?
            full->setMultiVector(block, partial, bThyraMode);
#else
    throw Xpetra::Exceptions::RuntimeError("MultiVector::setMultiVector() is not implemented.");
#endif
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    InsertVector(RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> partial,
                 size_t block,
                 RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> full,
                 bool bThyraMode) const {
  RCP<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> bfull =
      Teuchos::rcp_dynamic_cast<Xpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(full);
  if (bfull.is_null() == true)
    InsertVector(*partial, block, *full, bThyraMode);
  else {
    XPETRA_TEST_FOR_EXCEPTION(map_->getMap(block) == null, Xpetra::Exceptions::RuntimeError, "InsertVector: maps_[" << block << "] is null");
    bfull->setMultiVector(block, partial, bThyraMode);
  }
}

}  // namespace Xpetra

#endif  // XPETRA_BLOCKEDMULTIVECTOR_DEF_HPP
