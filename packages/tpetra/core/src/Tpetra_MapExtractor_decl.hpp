// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_MAPEXTRACTOR_DECL_HPP
#define TPETRA_MAPEXTRACTOR_DECL_HPP

/// \file Tpetra_MapExtractor_decl.hpp
/// \brief Declaration of the Tpetra::MapExtractor class.
///
/// Direct port of Xpetra::MapExtractor.  A MapExtractor stores a
/// Tpetra::BlockedMap and provides Extract/Insert operations to move data
/// between a monolithic (full) (Multi)Vector and its block sub-vectors, in
/// either Xpetra-style or Thyra-style block numbering.

#include <map>
#include <iostream>
#include <vector>

#include <Teuchos_RCP.hpp>
#include <Teuchos_Describable.hpp>

#include "Tpetra_ConfigDefs.hpp"

#include "Tpetra_MapExtractor_fwd.hpp"
#include "Tpetra_BlockedMap_fwd.hpp"
#include "Tpetra_Import_fwd.hpp"

#include "Tpetra_Map_decl.hpp"
#include "Tpetra_MultiVector_decl.hpp"
#include "Tpetra_Vector_decl.hpp"
#include "Tpetra_BlockedMap_decl.hpp"
#include "Tpetra_BlockedMultiVector_decl.hpp"

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class MapExtractor : public Teuchos::Describable {
 public:
  using scalar_type         = Scalar;
  using local_ordinal_type  = LocalOrdinal;
  using global_ordinal_type = GlobalOrdinal;
  using node_type           = Node;

 private:
  using Map              = ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using MultiVector      = ::Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Vector           = ::Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using BlockedMap       = ::Tpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node>;
  using BlockedMultiVector = ::Tpetra::BlockedMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

 public:
  //! MapExtractor basic constructor
  MapExtractor(const Teuchos::RCP<const Map>& fullmap, const std::vector<Teuchos::RCP<const Map>>& maps, bool bThyraMode = false);

  //! Expert constructor for Thyra maps
  MapExtractor(const std::vector<Teuchos::RCP<const Map>>& maps, const std::vector<Teuchos::RCP<const Map>>& thyramaps);

  //! Constructor which accepts a const version of a blocked map
  MapExtractor(const Teuchos::RCP<const BlockedMap>& blockedMap);

  //! copy constructor
  MapExtractor(const MapExtractor& input);

  //! Destructor.
  virtual ~MapExtractor();

  /** \name Extract subblocks from full map */
  //@{
  void ExtractVector(const Vector& full, size_t block, Vector& partial) const;
  void ExtractVector(const MultiVector& full, size_t block, MultiVector& partial) const;
  void ExtractVector(Teuchos::RCP<const Vector>& full, size_t block, Teuchos::RCP<Vector>& partial) const;
  void ExtractVector(Teuchos::RCP<Vector>& full, size_t block, Teuchos::RCP<Vector>& partial) const;
  void ExtractVector(Teuchos::RCP<const MultiVector>& full, size_t block, Teuchos::RCP<MultiVector>& partial) const;
  void ExtractVector(Teuchos::RCP<MultiVector>& full, size_t block, Teuchos::RCP<MultiVector>& partial) const;

  Teuchos::RCP<Vector> ExtractVector(Teuchos::RCP<const Vector>& full, size_t block, bool bThyraMode = false) const;
  Teuchos::RCP<Vector> ExtractVector(Teuchos::RCP<Vector>& full, size_t block, bool bThyraMode = false) const;
  Teuchos::RCP<MultiVector> ExtractVector(Teuchos::RCP<const MultiVector>& full, size_t block, bool bThyraMode = false) const;
  Teuchos::RCP<MultiVector> ExtractVector(Teuchos::RCP<MultiVector>& full, size_t block, bool bThyraMode = false) const;

  Teuchos::RCP<MultiVector> ExtractVector(Teuchos::RCP<const BlockedMultiVector>& full, size_t block, bool bThyraMode = false) const;
  Teuchos::RCP<MultiVector> ExtractVector(Teuchos::RCP<BlockedMultiVector>& full, size_t block, bool bThyraMode = false) const;

  //@}

  /** \name Insert subblocks into full map */
  //@{
  void InsertVector(const Vector& partial, size_t block, Vector& full, bool bThyraMode = false) const;
  void InsertVector(const MultiVector& partial, size_t block, MultiVector& full, bool bThyraMode = false) const;
  void InsertVector(Teuchos::RCP<const Vector> partial, size_t block, Teuchos::RCP<Vector> full, bool bThyraMode = false) const;
  void InsertVector(Teuchos::RCP<Vector> partial, size_t block, Teuchos::RCP<Vector> full, bool bThyraMode = false) const;
  void InsertVector(Teuchos::RCP<const MultiVector> partial, size_t block, Teuchos::RCP<MultiVector> full, bool bThyraMode = false) const;
  void InsertVector(Teuchos::RCP<MultiVector> partial, size_t block, Teuchos::RCP<MultiVector> full, bool bThyraMode = false) const;
  void InsertVector(Teuchos::RCP<const MultiVector> partial, size_t block, Teuchos::RCP<BlockedMultiVector> full, bool bThyraMode = false) const;
  void InsertVector(Teuchos::RCP<MultiVector> partial, size_t block, Teuchos::RCP<BlockedMultiVector> full, bool bThyraMode = false) const;

  //@}

  Teuchos::RCP<Vector> getVector(size_t i, bool bThyraMode = false, bool bZero = true) const;
  Teuchos::RCP<MultiVector> getVector(size_t i, size_t numvec, bool bThyraMode = false, bool bZero = true) const;

  /// returns true, if sub maps are stored in Thyra-style numbering
  bool getThyraMode() const;

  /** \name Maps */
  //@{

  /// number of partial maps
  size_t NumMaps() const;

  /// returns the sub map i from list of sub maps
  const Teuchos::RCP<const Map> getMap(size_t i, bool bThyraMode = false) const;

  /// get the (full) Map (the BlockedMap is not a Tpetra::Map, so this returns the full map)
  const Teuchos::RCP<const Map> getMap() const;

  /// get the underlying BlockedMap object
  const Teuchos::RCP<const BlockedMap> getBlockedMap() const;

  /// the full map
  const Teuchos::RCP<const Map> getFullMap() const;

  /// returns map index in map extractor which contains GID
  size_t getMapIndexForGID(GlobalOrdinal gid) const;

  //@}

 private:
  //! blocked map containing the sub block maps (either thyra or xpetra mode)
  Teuchos::RCP<const BlockedMap> map_;
};

}  // namespace Tpetra

#endif  // TPETRA_MAPEXTRACTOR_DECL_HPP
