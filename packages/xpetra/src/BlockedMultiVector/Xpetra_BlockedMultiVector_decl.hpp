// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_BLOCKEDMULTIVECTOR_DECL_HPP
#define XPETRA_BLOCKEDMULTIVECTOR_DECL_HPP

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Map_decl.hpp"
#include "Xpetra_MultiVector_decl.hpp"

#include "Xpetra_BlockedMap_decl.hpp"

// Xpetra::BlockedMultiVector is now a thin wrapper around Tpetra::BlockedMultiVector.
// The real blocked-vector logic (per-block storage, block-wise math, Merge, thyra/xpetra
// numbering) lives in Tpetra; this class holds an RCP<Tpetra::BlockedMultiVector> and
// forwards every method to it.  Argument MultiVectors coming in through the Xpetra
// interface are unwrapped to their underlying Tpetra objects (a blocked Xpetra vector
// unwraps to a Tpetra::BlockedMultiVector, a plain one to a Tpetra::MultiVector); a
// non-Tpetra (e.g. Epetra) argument triggers Xpetra::Exceptions::BadCast, which is the
// intended behavior in the Tpetra-only world.
#include <Tpetra_BlockedMultiVector_decl.hpp>

namespace Xpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// forward declaration of Vector, needed to prevent circular inclusions
template <class S, class LO, class GO, class N>
class Vector;

// forward declaration of MapExtractor, we just need the class sig here.
template <class S, class LO, class GO, class N>
class MapExtractor;
#endif

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class BlockedMultiVector
  : public MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  typedef Scalar scalar_type;
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;

  //! The type of the underlying Tpetra::BlockedMultiVector.
  using tpetra_blockedmultivector_type = ::Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

 private:
#undef XPETRA_BLOCKEDMULTIVECTOR_SHORT
#include "Xpetra_UseShortNames.hpp"

 public:
  //! @name Constructor/Destructor Methods
  //@{

  BlockedMultiVector(const Teuchos::RCP<const BlockedMap>& map, size_t NumVectors, bool zeroOut = true);

  BlockedMultiVector(Teuchos::RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> bmap,
                     Teuchos::RCP<const MultiVector> v);

  BlockedMultiVector(Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mapExtractor,
                     Teuchos::RCP<const MultiVector> v);

  BlockedMultiVector(const Teuchos::RCP<const BlockedMap>& map, std::vector<Teuchos::RCP<MultiVector>>& vin);

  //! Xpetra-specific constructor: wrap an existing Tpetra::BlockedMultiVector object.
  BlockedMultiVector(const Teuchos::RCP<tpetra_blockedmultivector_type>& vec);

  //! Destructor.
  virtual ~BlockedMultiVector();

  BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& operator=(const MultiVector& rhs);

  //@}
  //! @name Post-construction modification routines
  //@{

  virtual void replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar& value);
  virtual void sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar& value);
  virtual void replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar& value);
  virtual void sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar& value);

  //! Set all values in the multivector with the given value.
  virtual void putScalar(const Scalar& value);

  //@}
  //! @name Data Copy and View get methods
  //@{

  virtual Teuchos::RCP<const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> getVector(size_t j) const;
  virtual Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> getVectorNonConst(size_t j);
  virtual Teuchos::ArrayRCP<const Scalar> getData(size_t j) const;
  virtual Teuchos::ArrayRCP<Scalar> getDataNonConst(size_t j);

  //@}
  //! @name Mathematical methods
  //@{

  virtual void dot(const MultiVector& A, const Teuchos::ArrayView<Scalar>& dots) const;
  virtual void abs(const MultiVector& A);
  virtual void reciprocal(const MultiVector& A);
  virtual void scale(const Scalar& alpha);
  virtual void scale(Teuchos::ArrayView<const Scalar> alpha);
  virtual void update(const Scalar& alpha, const MultiVector& A, const Scalar& beta);
  virtual void update(const Scalar& alpha, const MultiVector& A, const Scalar& beta, const MultiVector& B, const Scalar& gamma);
  virtual void norm1(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const;
  virtual void norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const;
  virtual void normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& norms) const;
  virtual void meanValue(const Teuchos::ArrayView<Scalar>& means) const;
  virtual void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar& alpha, const MultiVector& A, const MultiVector& B, const Scalar& beta);
  virtual void elementWiseMultiply(Scalar scalarAB, const Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const MultiVector& B, Scalar scalarThis);

  //@}
  //! @name Attribute access functions
  //@{

  virtual size_t getNumVectors() const;
  virtual size_t getLocalLength() const;
  virtual global_size_t getGlobalLength() const;
  virtual bool isSameSize(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& vec) const;

  //@}
  //! @name Overridden from Teuchos::Describable
  //@{

  virtual std::string description() const;
  virtual void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const;

  virtual void replaceMap(const RCP<const Map>& map);

  //! Import.
  virtual void doImport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source, const Import& importer, CombineMode CM);

  //! Export.
  virtual void doExport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& dest, const Import& importer, CombineMode CM);

  //! Import (using an Exporter).
  virtual void doImport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& source, const Export& exporter, CombineMode CM);

  //! Export (using an Importer).
  virtual void doExport(const DistObject<Scalar, LocalOrdinal, GlobalOrdinal, Node>& dest, const Export& exporter, CombineMode CM);

  //@}
  //! @name Xpetra specific
  //@{

  //! Set seed for Random function.
  virtual void setSeed(unsigned int seed);

  virtual void randomize(bool bUseXpetraImplementation = false);

  virtual void randomize(const Scalar& minVal, const Scalar& maxVal, bool bUseXpetraImplementation = false);

  //! Set multi-vector values to random numbers. XPetra implementation
  virtual void Xpetra_randomize();

  //! Set multi-vector values to random numbers. XPetra implementation
  virtual void Xpetra_randomize(const Scalar& minVal, const Scalar& maxVal);

  //@}

  //! Access function for the underlying Map this DistObject was constructed with.
  Teuchos::RCP<const Map> getMap() const;

  //! Access function for the underlying Map this DistObject was constructed with.
  Teuchos::RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> getBlockedMap() const;

  /// return partial multivector associated with block row r
  Teuchos::RCP<MultiVector> getMultiVector(size_t r) const;

  /// return partial multivector associated with block row r
  Teuchos::RCP<MultiVector> getMultiVector(size_t r, bool bThyraMode) const;

  /// set partial multivector associated with block row r
  void setMultiVector(size_t r, Teuchos::RCP<const MultiVector> v, bool bThyraMode);

  /// merge BlockedMultiVector blocks to a single MultiVector
  Teuchos::RCP<MultiVector> Merge() const;

  //! @name Xpetra-specific accessors for the wrapped Tpetra object
  //@{

  //! Get the underlying Tpetra::BlockedMultiVector object.
  RCP<tpetra_blockedmultivector_type> getTpetra_BlockedMultiVector() const { return vec_; }

  //@}

 protected:
  virtual void assign(const MultiVector& rhs);

 private:
  //! The wrapped Tpetra::BlockedMultiVector.
  RCP<tpetra_blockedmultivector_type> vec_;

  //! Lazy Xpetra wrapper around the wrapped BlockedMap (identity cache).
  mutable RCP<const Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>> blockedMapXpetra_;

};  // BlockedMultiVector class

// Unwrap an Xpetra MultiVector to its underlying Tpetra MultiVector.  A blocked Xpetra
// vector unwraps to its Tpetra::BlockedMultiVector (returned as a Tpetra::MultiVector so
// that Tpetra's own dynamic_cast<BlockedMultiVector> in ExtractVector still fires); a
// plain Xpetra vector unwraps via the TpetraMultiVector toTpetra.  Epetra-backed vectors
// throw Xpetra::Exceptions::BadCast.
namespace BlockedMultiVectorDetails {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
unwrapMultiVector(const Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& X);

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
unwrapMultiVector(const Teuchos::RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& X);

// Wrap a Tpetra MultiVector back into an Xpetra MultiVector, preserving the blocked
// dynamic type: a Tpetra::BlockedMultiVector is wrapped in an Xpetra::BlockedMultiVector
// (so MueLu's rcp_dynamic_cast<BlockedMultiVector> keeps working on nested blocks); a
// plain Tpetra::MultiVector is wrapped via the TpetraMultiVector toXpetra.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
wrapMultiVector(const Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& X);

}  // namespace BlockedMultiVectorDetails

}  // namespace Xpetra

#define XPETRA_BLOCKEDMULTIVECTOR_SHORT
#endif  // XPETRA_BLOCKEDMULTIVECTOR_DECL_HPP
