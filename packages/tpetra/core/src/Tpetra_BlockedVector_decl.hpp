// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKEDVECTOR_DECL_HPP
#define TPETRA_BLOCKEDVECTOR_DECL_HPP

/// \file Tpetra_BlockedVector_decl.hpp
/// \brief Declaration of the Tpetra::BlockedVector class.
///
/// Tpetra::BlockedVector is a direct port of Xpetra::BlockedVector.  Xpetra
/// realizes BlockedVector via virtual multiple inheritance from both
/// Xpetra::Vector and Xpetra::BlockedMultiVector.  Tpetra::Vector inherits
/// Tpetra::MultiVector *non-virtually*, so the analogous multiple inheritance
/// would create an ambiguous MultiVector base.  Instead, Tpetra::BlockedVector
/// derives (only) from Tpetra::Vector (an empty base sub-object) and *holds* a
/// Tpetra::BlockedMultiVector member (\c bmv_) that owns the block storage; all
/// blocked behavior is delegated to \c bmv_.  \c getBlockedMultiVector() exposes
/// that member so consumers that need to treat the vector as a
/// BlockedMultiVector (which a dynamic_cast cannot recover here) can do so.
///
/// The same non-virtual-shadowing caveat documented on Tpetra::BlockedMultiVector
/// applies: correct blocked behavior is obtained only through the BlockedVector
/// static type (or the Xpetra virtual wrapper), not through a bare
/// Tpetra::MultiVector& / Tpetra::Vector& handle.

#include "Tpetra_ConfigDefs.hpp"

#include "Tpetra_BlockedVector_fwd.hpp"
#include "Tpetra_BlockedMultiVector_fwd.hpp"
#include "Tpetra_BlockedMap_fwd.hpp"
#include "Tpetra_MapExtractor_fwd.hpp"

#include "Tpetra_Vector_decl.hpp"
#include "Tpetra_BlockedMultiVector_decl.hpp"
#include "Tpetra_BlockedMap_decl.hpp"

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class BlockedVector
  : public Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  using scalar_type         = Scalar;
  using local_ordinal_type  = LocalOrdinal;
  using global_ordinal_type = GlobalOrdinal;
  using node_type           = Node;

 private:
  using base_type         = ::Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Vector            = ::Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using MultiVector       = ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using BlockedMultiVector = ::Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Map               = ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using BlockedMap        = ::Tpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>;
  using MapExtractor      = ::Tpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

 public:
  using map_type         = Map;
  using impl_scalar_type = typename base_type::impl_scalar_type;
  using dot_type         = typename base_type::dot_type;
  using mag_type         = typename base_type::mag_type;

  //! @name Constructor/Destructor Methods
  //@{

  BlockedVector(const Teuchos::RCP<const BlockedMap>& map, bool zeroOut = true);

  BlockedVector(Teuchos::RCP<const BlockedMap> bmap,
                Teuchos::RCP<Vector> v);

  BlockedVector(Teuchos::RCP<const MapExtractor> mapExtractor,
                Teuchos::RCP<Vector> v);

  //! Destructor.
  virtual ~BlockedVector();

  /// \brief Assignment operator: does a deep copy.
  BlockedVector& operator=(const MultiVector& rhs);

  //@}
  //! @name Post-construction modification routines
  //@{

  virtual void replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar& value);
  virtual void sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar& value);
  virtual void replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar& value);
  virtual void sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar& value);

  virtual void replaceGlobalValue(GlobalOrdinal globalRow, const Scalar& value);
  virtual void sumIntoGlobalValue(GlobalOrdinal globalRow, const Scalar& value);
  virtual void replaceLocalValue(LocalOrdinal myRow, const Scalar& value);
  virtual void sumIntoLocalValue(LocalOrdinal myRow, const Scalar& value);

  //! Set all values in the vector with the given value.
  virtual void putScalar(const Scalar& value);

  //@}
  //! @name Data Copy and View get methods
  //@{

  virtual Teuchos::RCP<const Vector> getVector(size_t j) const;
  virtual Teuchos::RCP<Vector> getVectorNonConst(size_t j);
  virtual Teuchos::ArrayRCP<const Scalar> getData(size_t j) const;
  virtual Teuchos::ArrayRCP<Scalar> getDataNonConst(size_t j);

  //@}
  //! @name Mathematical methods
  //@{

  virtual void dot(const MultiVector& A, const Teuchos::ArrayView<dot_type>& dots) const;
  virtual dot_type dot(const Vector& A) const;
  virtual void abs(const MultiVector& A);
  virtual void reciprocal(const MultiVector& A);
  virtual void scale(const Scalar& alpha);
  virtual void scale(Teuchos::ArrayView<const Scalar> alpha);
  virtual void update(const Scalar& alpha, const MultiVector& A, const Scalar& beta);
  virtual void update(const Scalar& alpha, const MultiVector& A, const Scalar& beta, const MultiVector& B, const Scalar& gamma);
  virtual mag_type norm1() const;
  virtual mag_type norm2() const;
  virtual mag_type normInf() const;
  virtual void norm1(const Teuchos::ArrayView<mag_type>& norms) const;
  virtual void norm2(const Teuchos::ArrayView<mag_type>& norms) const;
  virtual void normInf(const Teuchos::ArrayView<mag_type>& norms) const;
  virtual void meanValue(const Teuchos::ArrayView<impl_scalar_type>& means) const;
  virtual impl_scalar_type meanValue() const;
  virtual void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar& alpha, const Vector& A, const Vector& B, const Scalar& beta);
  virtual void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar& alpha, const MultiVector& A, const MultiVector& B, const Scalar& beta);
  virtual void elementWiseMultiply(Scalar scalarAB, const Vector& A, const MultiVector& B, Scalar scalarThis);
  virtual void elementWiseMultiply(Scalar scalarAB, const Vector& A, const Vector& B, Scalar scalarThis);

  //@}
  //! @name Attribute access functions
  //@{

  virtual size_t getNumVectors() const;
  virtual size_t getLocalLength() const;
  virtual global_size_t getGlobalLength() const;
  virtual bool isSameSize(const MultiVector& vec) const;

  //@}
  //! @name Overridden from Teuchos::Describable
  //@{

  std::string description() const override;
  void describe(Teuchos::FancyOStream& out,
                const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const override;

  //@}

  virtual void replaceMap(const Teuchos::RCP<const Map>& map);

  virtual void randomize(bool bUseXpetraImplementation = false);
  virtual void randomize(const Scalar& minVal, const Scalar& maxVal, bool bUseXpetraImplementation = false);

  //! Access function for the (full) Map this object represents.
  Teuchos::RCP<const Map> getMap() const override;

  //! Access the underlying BlockedMap.
  Teuchos::RCP<const BlockedMap> getBlockedMap() const;

  /// return partial multivector associated with block row r
  Teuchos::RCP<MultiVector> getMultiVector(size_t r) const;

  /// return partial multivector associated with block row r
  Teuchos::RCP<MultiVector> getMultiVector(size_t r, bool bThyraMode) const;

  /// set partial Vector associated with block row r
  void setMultiVector(size_t r, Teuchos::RCP<const Vector> v, bool bThyraMode);

  /// merge BlockedVector blocks to a single Vector
  Teuchos::RCP<MultiVector> Merge() const;

  //! Expose the internal BlockedMultiVector (blocked identity is not recoverable via dynamic_cast here).
  Teuchos::RCP<BlockedMultiVector> getBlockedMultiVector() const { return bmv_; }

 protected:
  /// \brief Implementation of the assignment operator (operator=); does a deep copy.
  virtual void assign(const MultiVector& rhs);

 private:
  Teuchos::RCP<BlockedMultiVector> bmv_;  ///< the blocked multivector (single column) that owns the block storage

};  // BlockedVector class

}  // namespace Tpetra

#endif  // TPETRA_BLOCKEDVECTOR_DECL_HPP
