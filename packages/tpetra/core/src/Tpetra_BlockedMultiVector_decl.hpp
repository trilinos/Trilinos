// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKEDMULTIVECTOR_DECL_HPP
#define TPETRA_BLOCKEDMULTIVECTOR_DECL_HPP

/// \file Tpetra_BlockedMultiVector_decl.hpp
/// \brief Declaration of the Tpetra::BlockedMultiVector class.
///
/// Tpetra::BlockedMultiVector is a direct port of Xpetra::BlockedMultiVector.
/// It represents a MultiVector that is partitioned into an ordered list of
/// sub-MultiVectors according to a Tpetra::BlockedMap.  Unlike a plain
/// Tpetra::MultiVector it does NOT own a single monolithic block of storage:
/// the block sub-vectors (\c vv_) are the only storage.  It derives from
/// Tpetra::MultiVector purely so that it can be handled polymorphically by the
/// blocked-operator machinery (Tpetra::MapExtractor, Tpetra::BlockedCrsMatrix)
/// and by the thin Xpetra wrapper.
///
/// IMPORTANT: Tpetra::MultiVector's mathematical methods (dot/norm/scale/...)
/// are *not* virtual (only getMap(), description(), and the DistObject transfer
/// hooks are).  BlockedMultiVector therefore *shadows* the base math methods
/// (declaring its own \c virtual versions) rather than overriding them.  The
/// blocked behavior is only obtained when the object is manipulated through the
/// BlockedMultiVector static type (or through the Xpetra virtual interface that
/// wraps it).  A consumer that treats a BlockedMultiVector polymorphically as a
/// bare <tt>Tpetra::MultiVector&</tt> and calls a non-virtual math method will
/// hit the (empty) base implementation.  This faithfully preserves the behavior
/// of the current Xpetra-routed code paths, which always dispatch through the
/// blocked type.

#include "Tpetra_ConfigDefs.hpp"

#include <vector>

#include "Tpetra_BlockedMultiVector_fwd.hpp"
#include "Tpetra_BlockedMap_fwd.hpp"
#include "Tpetra_MapExtractor_fwd.hpp"
#include "Tpetra_Vector_fwd.hpp"
#include "Tpetra_Import_fwd.hpp"
#include "Tpetra_Export_fwd.hpp"

#include "Tpetra_MultiVector_decl.hpp"
#include "Tpetra_BlockedMap_decl.hpp"

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class BlockedMultiVector
  : public MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 public:
  using scalar_type         = Scalar;
  using local_ordinal_type  = LocalOrdinal;
  using global_ordinal_type = GlobalOrdinal;
  using node_type           = Node;

 private:
  using base_type      = ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using MultiVector    = ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Vector         = ::Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Map            = ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using BlockedMap     = ::Tpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>;
  using MapExtractor   = ::Tpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Import         = ::Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>;
  using Export         = ::Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>;

 public:
  using map_type       = Map;
  using impl_scalar_type = typename base_type::impl_scalar_type;
  using dot_type         = typename base_type::dot_type;
  using mag_type         = typename base_type::mag_type;

  //! @name Constructor/Destructor Methods
  //@{

  /*!
   * \param map        BlockedMap defining the block structure of the multi vector
   * \param NumVectors Number of vector columns in multi vector
   * \param zeroOut    If true initialize multivector with zeros
   */
  BlockedMultiVector(const Teuchos::RCP<const BlockedMap>& map, size_t NumVectors, bool zeroOut = true);

  /*!
   * Split the input (monolithic) MultiVector \c v into blocks according to \c bmap.
   * A copy is made; the blocked multivector is not a view of \c v.
   */
  BlockedMultiVector(Teuchos::RCP<const BlockedMap> bmap,
                     Teuchos::RCP<const MultiVector> v);

  /*!
   * Split the input (monolithic) MultiVector \c v into blocks according to a MapExtractor.
   * A copy is made; the blocked multivector is not a view of \c v.
   */
  BlockedMultiVector(Teuchos::RCP<const MapExtractor> mapExtractor,
                     Teuchos::RCP<const MultiVector> v);

  /*!
   * Build a blocked multivector from a blocked map and component vectors.
   * \note We do *NOT* check map compatibility between the BlockedMap and the array of multivectors
   */
  BlockedMultiVector(const Teuchos::RCP<const BlockedMap>& map, std::vector<Teuchos::RCP<MultiVector>>& vin);

  //! Destructor.
  virtual ~BlockedMultiVector();

  /// \brief Assignment operator: does a deep copy (dispatches to assign()).
  BlockedMultiVector& operator=(const MultiVector& rhs);

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

  //! Return a Vector which is a const view of column j.
  virtual Teuchos::RCP<const Vector> getVector(size_t j) const;

  //! Return a Vector which is a nonconst view of column j.
  virtual Teuchos::RCP<Vector> getVectorNonConst(size_t j);

  //! Const view of the local values in a particular vector of this multivector.
  virtual Teuchos::ArrayRCP<const Scalar> getData(size_t j) const;

  //! View of the local values in a particular vector of this multivector.
  virtual Teuchos::ArrayRCP<Scalar> getDataNonConst(size_t j);

  //@}
  //! @name Mathematical methods
  //@{

  virtual void dot(const MultiVector& A, const Teuchos::ArrayView<dot_type>& dots) const;
  virtual void abs(const MultiVector& A);
  virtual void reciprocal(const MultiVector& A);
  virtual void scale(const Scalar& alpha);
  virtual void scale(const Teuchos::ArrayView<const Scalar>& alpha);
  virtual void update(const Scalar& alpha, const MultiVector& A, const Scalar& beta);
  virtual void update(const Scalar& alpha, const MultiVector& A, const Scalar& beta, const MultiVector& B, const Scalar& gamma);
  virtual void norm1(const Teuchos::ArrayView<mag_type>& norms) const;
  virtual void norm2(const Teuchos::ArrayView<mag_type>& norms) const;
  virtual void normInf(const Teuchos::ArrayView<mag_type>& norms) const;
  virtual void meanValue(const Teuchos::ArrayView<impl_scalar_type>& means) const;
  virtual void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar& alpha, const MultiVector& A, const MultiVector& B, const Scalar& beta);
  virtual void elementWiseMultiply(Scalar scalarAB, const Vector& A, const MultiVector& B, Scalar scalarThis);

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

  //! Replace the underlying Map in place.  Accepts a BlockedMap (recommended) or a plain Map (single-block only).
  virtual void replaceMap(const Teuchos::RCP<const Map>& map);
  virtual void replaceMap(const Teuchos::RCP<const BlockedMap>& map);

  //! Set (multi)vector values to random numbers.
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

  /// set partial multivector associated with block row r
  void setMultiVector(size_t r, Teuchos::RCP<const MultiVector> v, bool bThyraMode);

  /// merge BlockedMultiVector blocks to a single (monolithic) MultiVector
  Teuchos::RCP<MultiVector> Merge() const;

 protected:
  /// \brief Implementation of the assignment operator (operator=); does a deep copy.
  virtual void assign(const MultiVector& rhs);

 private:
  // helper routines for interaction of MultiVector and BlockedMultiVectors
  void ExtractVector(Teuchos::RCP<const MultiVector>& full, size_t block, Teuchos::RCP<MultiVector>& partial) const;
  void ExtractVector(Teuchos::RCP<MultiVector>& full, size_t block, Teuchos::RCP<MultiVector>& partial) const;

  Teuchos::RCP<MultiVector> ExtractVector(Teuchos::RCP<const MultiVector>& full, size_t block, bool bThyraMode = false) const;
  Teuchos::RCP<MultiVector> ExtractVector(Teuchos::RCP<MultiVector>& full, size_t block, bool bThyraMode = false) const;

  void ExtractVector(const MultiVector& full, size_t block, MultiVector& partial) const;

  void InsertVector(const MultiVector& partial, size_t block, MultiVector& full, bool bThyraMode = false) const;
  void InsertVector(Teuchos::RCP<const MultiVector> partial, size_t block, Teuchos::RCP<MultiVector> full, bool bThyraMode = false) const;
  void InsertVector(Teuchos::RCP<MultiVector> partial, size_t block, Teuchos::RCP<MultiVector> full, bool bThyraMode = false) const;

 private:
  Teuchos::RCP<const BlockedMap> map_;                  ///< blocked map containing the sub block maps (either thyra or xpetra mode)
  std::vector<Teuchos::RCP<MultiVector>> vv_;           ///< array containing RCPs of the partial vectors
  size_t numVectors_;                                   ///< number of vectors (columns in multi vector)

};  // BlockedMultiVector class

}  // namespace Tpetra

#endif  // TPETRA_BLOCKEDMULTIVECTOR_DECL_HPP
