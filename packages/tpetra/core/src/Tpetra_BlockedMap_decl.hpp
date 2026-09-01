// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_BLOCKEDMAP_DECL_HPP
#define TPETRA_BLOCKEDMAP_DECL_HPP

/// \file Tpetra_BlockedMap_decl.hpp
/// \brief Declaration of the Tpetra::BlockedMap class.
///
/// A Tpetra::BlockedMap describes a partitioning of a "full" Tpetra::Map into
/// an ordered list of sub-maps.  It is the map-side infrastructure that backs
/// Tpetra::MapExtractor, Tpetra::BlockedMultiVector, and
/// Tpetra::BlockedCrsMatrix (the "matrix of matrices").  It is a direct port of
/// Xpetra::BlockedMap, but -- because Tpetra::Map is a value-like, non-derivable
/// class -- BlockedMap is a standalone Teuchos::Describable that *composes* a
/// full Tpetra::Map plus its sub-maps rather than deriving from Tpetra::Map.
/// Most of the Map-like query surface simply forwards to the full map, exactly
/// as the Xpetra implementation does.

#include "Tpetra_ConfigDefs.hpp"

#include <vector>
#include <utility>

#include <Teuchos_Describable.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_Comm.hpp>

#include "Tpetra_BlockedMap_fwd.hpp"
#include "Tpetra_Map_decl.hpp"
#include "Tpetra_Import_decl.hpp"

namespace Tpetra {

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class BlockedMap : public Teuchos::Describable {
 public:
  using local_ordinal_type  = LocalOrdinal;
  using global_ordinal_type = GlobalOrdinal;
  using node_type           = Node;

  using map_type    = ::Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;
  using import_type = ::Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>;

  using local_map_type = typename map_type::local_map_type;
  // Tpetra::Map declares this type privately, but exposes it through the public
  // getMyGlobalIndicesDevice() accessor; recover it via decltype.
  using global_indices_array_device_type =
      decltype(std::declval<const map_type>().getMyGlobalIndicesDevice());

  //! @name Constructor/Destructor Methods
  //@{

  //! Default constructor (creates an empty, Xpetra-mode BlockedMap).
  BlockedMap();

  //! BlockedMap basic constructor
  //!
  //! @param[in] fullmap Full map containing all GIDs throughout the full vector. This parameter is only
  //!                    important if bThyraMode == false (see below)
  //! @param[in] maps    Vector containing submaps. The set of all GIDs stored in the submaps should be
  //!                    the same than stored in fullmap, if bThyraMode == false. In Thyra mode, the
  //!                    submaps should contain consecutive GIDs starting with 0 in each submap.
  //! @param[in] bThyraMode Flag which allows to switch between generating a BlockedMap in Thyra mode
  //!                       or Xpetra mode
  //!
  //! In Thyra mode, fullmap is not important as a fullmap with unique blocked GIDs is automatically
  //! generated which map the GIDs of the submaps to uniquely defined GIDs in the fullmap.
  //!
  //! In Xpetra mode, the fullmap has to be the same as the union of the GIDs stored in the submaps in
  //! maps. The intersection of the GIDs of the sub-maps in maps must be empty.
  BlockedMap(const Teuchos::RCP<const map_type>& fullmap,
             const std::vector<Teuchos::RCP<const map_type>>& maps,
             bool bThyraMode = false);

  //! Expert constructor for Thyra maps
  BlockedMap(const std::vector<Teuchos::RCP<const map_type>>& maps,
             const std::vector<Teuchos::RCP<const map_type>>& thyramaps);

  //! copy constructor
  BlockedMap(const BlockedMap& input);

  //! Destructor.
  virtual ~BlockedMap();

  //@}

  //! @name Attributes
  //@{

  //! The number of elements in this Map.
  global_size_t getGlobalNumElements() const;

  //! The number of elements belonging to the calling process.
  size_t getLocalNumElements() const;

  //! The index base for this Map.
  GlobalOrdinal getIndexBase() const;

  //! The minimum local index.
  LocalOrdinal getMinLocalIndex() const;

  //! The maximum local index on the calling process.
  LocalOrdinal getMaxLocalIndex() const;

  //! The minimum global index owned by the calling process.
  GlobalOrdinal getMinGlobalIndex() const;

  //! The maximum global index owned by the calling process.
  GlobalOrdinal getMaxGlobalIndex() const;

  //! The minimum global index over all processes in the communicator.
  GlobalOrdinal getMinAllGlobalIndex() const;

  //! The maximum global index over all processes in the communicator.
  GlobalOrdinal getMaxAllGlobalIndex() const;

  //! The local index corresponding to the given global index.
  LocalOrdinal getLocalElement(GlobalOrdinal globalIndex) const;

  //! The global index corresponding to the given local index.
  GlobalOrdinal getGlobalElement(LocalOrdinal localIndex) const;

  //! Return the process ranks and corresponding local indices for the given global indices.
  LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal>& /* GIDList    */,
                                  const Teuchos::ArrayView<int>& /* nodeIDList */,
                                  const Teuchos::ArrayView<LocalOrdinal>& /* LIDList    */) const;

  //! Return the process ranks for the given global indices.
  LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal>& /* GIDList    */,
                                  const Teuchos::ArrayView<int>& /* nodeIDList */) const;

  //! Return a view of the global indices owned by this process.
  Teuchos::ArrayView<const GlobalOrdinal> getLocalElementList() const;

  //! Return a device view of the global indices owned by this process.
  global_indices_array_device_type getMyGlobalIndicesDevice() const;

  //@}

  //! @name Boolean tests
  //@{

  //! Whether the given local index is valid for this Map on this process.
  bool isNodeLocalElement(LocalOrdinal localIndex) const;

  //! Whether the given global index is valid for this Map on this process.
  bool isNodeGlobalElement(GlobalOrdinal globalIndex) const;

  //! True if this Map is distributed contiguously, else false.
  bool isContiguous() const;

  //! Whether this Map is globally distributed or locally replicated.
  bool isDistributed() const;

  //! True if and only if map is compatible with this Map.
  bool isCompatible(const map_type& map) const;

  //! True if and only if map is compatible with this Map.
  bool isCompatible(const BlockedMap& map) const;

  //! True if and only if map is identical to this Map.
  bool isSameAs(const map_type& map) const;

  //! True if and only if map is identical to this Map.
  bool isSameAs(const BlockedMap& map) const;

  //@}

  //! @name
  //@{

  //! Get this Map's Comm object.
  Teuchos::RCP<const Teuchos::Comm<int>> getComm() const;

  //@}

  /// \brief Assignment operator: Does a (shallow, immutable-map) copy.
  BlockedMap& operator=(const BlockedMap& rhs);

  //! @name Attribute access functions
  //@{

  //! returns true if internally stored sub maps are in Thyra mode (i.e. start all with GIDs=0)
  bool getThyraMode() const;

  //@}

  //! @name BlockedMap specific access
  //@{

  /// number of partial maps
  size_t getNumMaps() const;

  /// get the sub map i from the list of sub maps
  /// depending on the parameter bThyraMode the sub map that is returned uses Thyra or Xpetra numbering.
  /// Note: Thyra-numbering is only allowed if the BlockedMap is also constructed using Thyra numbering.
  Teuchos::RCP<const map_type> getMap(size_t i, bool bThyraMode = false) const;

  /// get the importer between full map and partial map
  Teuchos::RCP<import_type> getImporter(size_t i) const;

  /// the full map
  Teuchos::RCP<const map_type> getFullMap() const;

  /// returns map index in map extractor which contains GID
  size_t getMapIndexForGID(GlobalOrdinal gid) const;

  /// \brief Get the local Map for Kokkos kernels.
  local_map_type getLocalMap() const { return fullmap_->getLocalMap(); }

  /// \brief Record that sub-block \c i is itself a BlockedMap.
  ///
  /// Because Tpetra::Map is not derivable, \c maps_[i] can only ever be a
  /// flattened plain Map.  To preserve *nested* blocked identity (needed so
  /// that a BlockedMultiVector over sub-block \c i is itself a nested
  /// BlockedMultiVector, mirroring Xpetra where BlockedMap : public Map), the
  /// nested BlockedMap is recorded here alongside the flattened \c maps_[i].
  /// Only the reorder builders, which know the true nested structure, call
  /// this; leaf sub-blocks leave \c nestedMaps_[i] null.
  void setBlockedSubMap(size_t i, const Teuchos::RCP<const BlockedMap>& submap);

  /// \brief The nested BlockedMap for sub-block \c i, or null if sub-block \c i
  /// is a leaf (plain) map.  \sa setBlockedSubMap.
  Teuchos::RCP<const BlockedMap> getBlockedSubMap(size_t i) const;

  /// \brief Whether sub-block \c i is itself a BlockedMap.  \sa setBlockedSubMap.
  bool isBlocked(size_t i) const;

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! A simple one-line description of this object.
  std::string description() const override;

  //! Print the object with the given verbosity level to a FancyOStream.
  void describe(Teuchos::FancyOStream& out,
                const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const override;

  //@}

  /*! @brief Helper function to concatenate several maps

    @param  subMaps    vector of maps which are concatenated
    @return            concatenated map

    The routine builds a global map by concatenating all provided maps in the ordering defined by the vector.
    The GIDs are just appended in the same ordering as in the subMaps. No reordering or sorting is performed.
   */
  static Teuchos::RCP<const map_type>
  concatenateMaps(const std::vector<Teuchos::RCP<const map_type>>& subMaps);

 protected:
  /// \brief Implementation of the assignment operator (operator=).
  virtual void assign(const BlockedMap& input);

 private:
  bool CheckConsistency() const;

 private:
  Teuchos::RCP<const map_type> fullmap_;
  std::vector<Teuchos::RCP<const map_type>> maps_;
  std::vector<Teuchos::RCP<import_type>> importers_;
  bool bThyraMode_;                                        //< boolean flag: use Thyra numbering for local sub-block maps. default = false (Xpetra mode)
  std::vector<Teuchos::RCP<const map_type>> thyraMaps_;    //< Thyra-style numbering maps (empty in Xpetra mode).
  std::vector<Teuchos::RCP<const BlockedMap>> nestedMaps_;  //< nested BlockedMap per sub-block, parallel to maps_ (null entry = leaf). \sa setBlockedSubMap
};                                                          // BlockedMap class

}  // namespace Tpetra

#endif  // TPETRA_BLOCKEDMAP_DECL_HPP
