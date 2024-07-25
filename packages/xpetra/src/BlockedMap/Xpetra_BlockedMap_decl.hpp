// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PACKAGES_XPETRA_SUP_BLOCKEDMAP_XPETRA_BLOCKEDMAP_DECL_HPP_
#define PACKAGES_XPETRA_SUP_BLOCKEDMAP_XPETRA_BLOCKEDMAP_DECL_HPP_

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Import.hpp"
#include "Xpetra_Map_decl.hpp"

namespace Xpetra {

template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class BlockedMap : public Map<LocalOrdinal, GlobalOrdinal, Node> {
 public:
  typedef LocalOrdinal local_ordinal_type;
  typedef GlobalOrdinal global_ordinal_type;
  typedef Node node_type;
  typedef typename Map<LocalOrdinal, GlobalOrdinal, Node>::global_indices_array_device_type global_indices_array_device_type;

 private:
#undef XPETRA_BLOCKEDMAP_SHORT
#include "Xpetra_UseShortNamesOrdinal.hpp"

 public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor
  /*!
   */
  BlockedMap();

  //! BlockedMap basic constructor
  //!
  //! @param[in] fullmap Full map containing all GIDs throughout the full vector. This parameter is only
  //!                    important if bThyraMode == false (see below)
  //! @param[in] maps    Vector containing submaps. The set of all GIDs stored in the submaps should be
  //!                    the same than stored in fullmap, if bThyraMode == false. In Thyra mode, the
  //!                    submaps should contain consecutive GIDs starting with 0 in each submap.
  //! @param[in] bThyraMode Flag which allows to switch between generating a MapExtractor in Thyra mode
  //!                       or Xpetra mode
  //!
  //! In Thyra mode, fullmap is not important as a fullmap with unique blocked GIDs is automatically
  //! generated which map the GIDs of the submaps to uniquely defined GIDs in the fullmap. The user
  //! has to provide a fullmap in Thyra mode to specify the underlying linear algebra library
  //! (Epetra or Tpetra) and some other map information (e.g. indexBase). This could be fixed.
  //!
  //! In Xpetra mode, the fullmap has to be the same as the union of the GIDs stored in the submaps in
  //!  maps. The intersection of the GIDs of the sub- maps in maps must be empty.
  BlockedMap(const RCP<const Map>& fullmap, const std::vector<RCP<const Map>>& maps, bool bThyraMode = false);

  //! Expert constructor for Thyra maps
  BlockedMap(const std::vector<RCP<const Map>>& maps, const std::vector<RCP<const Map>>& thyramaps);

  //! copy constructor
  BlockedMap(const BlockedMap& input);

  //! Destructor.
  virtual ~BlockedMap();

  //! @name Attributes
  //@{

  //! The number of elements in this Map.
  virtual global_size_t getGlobalNumElements() const;

  //! The number of elements belonging to the calling process.
  virtual size_t getLocalNumElements() const;

  //! The index base for this Map.
  virtual GlobalOrdinal getIndexBase() const;

  //! The minimum local index.
  virtual LocalOrdinal getMinLocalIndex() const;

  //! The maximum local index on the calling process.
  virtual LocalOrdinal getMaxLocalIndex() const;

  //! The minimum global index owned by the calling process.
  virtual GlobalOrdinal getMinGlobalIndex() const;

  //! The maximum global index owned by the calling process.
  virtual GlobalOrdinal getMaxGlobalIndex() const;

  //! The minimum global index over all processes in the communicator.
  virtual GlobalOrdinal getMinAllGlobalIndex() const;

  //! The maximum global index over all processes in the communicator.
  virtual GlobalOrdinal getMaxAllGlobalIndex() const;

  //! The local index corresponding to the given global index.
  virtual LocalOrdinal getLocalElement(GlobalOrdinal globalIndex) const;

  //! The global index corresponding to the given local index.
  virtual GlobalOrdinal getGlobalElement(LocalOrdinal localIndex) const;

  //! Return the process ranks and corresponding local indices for the given global indices.
  virtual LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal>& /* GIDList    */,
                                          const Teuchos::ArrayView<int>& /* nodeIDList */,
                                          const Teuchos::ArrayView<LocalOrdinal>& /* LIDList    */) const;

  //! Return the process ranks for the given global indices.
  virtual LookupStatus getRemoteIndexList(const Teuchos::ArrayView<const GlobalOrdinal>& /* GIDList    */,
                                          const Teuchos::ArrayView<int>& /* nodeIDList */) const;

  //! Return a view of the global indices owned by this process.
  virtual Teuchos::ArrayView<const GlobalOrdinal> getLocalElementList() const;

  //! Return a view of the global indices owned by this process.
  virtual global_indices_array_device_type getMyGlobalIndicesDevice() const;

  //@}

  //! @name Boolean tests
  //@{

  //! Whether the given local index is valid for this Map on this process.
  virtual bool isNodeLocalElement(LocalOrdinal localIndex) const;

  //! Whether the given global index is valid for this Map on this process.
  virtual bool isNodeGlobalElement(GlobalOrdinal globalIndex) const;

  //! True if this Map is distributed contiguously, else false.
  virtual bool isContiguous() const;

  //! Whether this Map is globally distributed or locally replicated.
  virtual bool isDistributed() const;

  //! True if and only if map is compatible with this Map.
  virtual bool isCompatible(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map) const;

  //! True if and only if map is identical to this Map.
  virtual bool isSameAs(const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map) const;

  //@}

  //! @name
  //@{

  //! Get this Map's Comm object.
  virtual Teuchos::RCP<const Teuchos::Comm<int>> getComm() const;

  //@}

  /// \brief Assignment operator: Does a deep copy.
  ///
  /// The assignment operator does a deep copy, just like
  /// subclasses' copy constructors.
  ///
  /// \note This currently only works if both <tt>*this</tt> and the
  ///   input argument are instances of the same subclass.
  BlockedMap<LocalOrdinal, GlobalOrdinal, Node>& operator=(const BlockedMap& rhs);

  //@}

  //! @name Attribute access functions
  //@{
  //! Local number of rows on the calling process.
  /*virtual size_t getLocalLength() const {
    throw Xpetra::Exceptions::RuntimeError("BlockedMap::getLocalLength: routine not implemented.");
    return 0;
  }*/

  //! Global number of rows in the multivector.
  /*virtual global_size_t getGlobalLength() const {
    throw Xpetra::Exceptions::RuntimeError("BlockedMap::getGlobalLength: routine not implemented.");
    return 0;
  }*/

  //! returns true if internally stored sub maps are in Thyra mode (i.e. start all with GIDs=0)
  virtual bool getThyraMode() const;

  //@}

  //! \name Maps
  //@{

  //! Return a new Map with processes with zero elements removed.
  virtual RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> removeEmptyProcesses() const;

  //! Replace this Map's communicator with a subset communicator.
  virtual RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
  replaceCommWithSubset(const Teuchos::RCP<const Teuchos::Comm<int>>& /* newComm */) const;

  //@}

  //! @name Xpetra specific
  //@{

  //! Get the library used by this object (Tpetra or Epetra?)
  virtual UnderlyingLib lib() const;

  // TODO: find a better solution for this hack
  // The problem is that EpetraMap, TpetraMap and StridedMap all inherit Map. To have proper toEpetra() we
  // need to understand the type of underlying matrix. But in src/Map we have no knowledge of StridedMaps, so
  // we cannot check for it by casting. This function allows us to avoid the restriction, as StridedMap redefines
  // it to return the base map.
  virtual RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
  getMap() const;

  //@}

  /// number of partial maps
  size_t getNumMaps() const;

  /// get the map
  /// returns the sub map i from list of sub maps
  /// depending on the parameter bThyraMode the sub map that is returned uses Thyra or Xpetra numbering
  /// Note: Thyra-numbering is only allowed if the BlockedMap is also constructed using Thyra numbering
  const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
  getMap(size_t i, bool bThyraMode = false) const;

  /// get the importer between full map and partial map
  const RCP<Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node>>
  getImporter(size_t i) const;

  /// the full map
  const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
  getFullMap() const;

  /// returns map index in map extractor which contains GID
  size_t getMapIndexForGID(GlobalOrdinal gid) const;

#ifdef HAVE_XPETRA_TPETRA
  using local_map_type = typename Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>::local_map_type;

  /// \brief Get the local Map for Kokkos kernels.
  local_map_type getLocalMap() const { return fullmap_->getLocalMap(); }

#else  // HAVE_XPETRA_TPETRA
#ifdef __GNUC__
#warning "Xpetra Kokkos interface for CrsMatrix is enabled (HAVE_XPETRA_KOKKOS_REFACTOR) but Tpetra is disabled. The Kokkos interface needs Tpetra to be enabled, too."
#endif  // __GNUC__
#endif  // #else !HAVE_XPETRA_TPETRA

  //@}

  //! @name Overridden from Teuchos::Describable
  //@{

  //! A simple one-line description of this object.
  virtual std::string description() const;

  //! Print the object with the given verbosity level to a FancyOStream.
  virtual void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const;

  //@}

 protected:
  /// \brief Implementation of the assignment operator (operator=);
  ///   does a deep copy.
  virtual void assign(const BlockedMap& input);

  /*! @brief Helper function to concatenate several maps

    @param  subMaps    vector of maps which are concatenated
    @return            concatenated map

    The routine builds a global map by concatenating all provided maps in the ordering defined by the vector.
    The GIDs are just appended in the same ordering as in the subMaps. No reordering or sorting is performed.
    This routine is supposed to generate the full map in an Xpetra::MapExtractor for a block operator. Note, it
    should not be used for strided maps since the GIDs are not reordered.

    Example: subMap[0] = { 0, 1, 3, 4 };
             subMap[1] = { 2, 5 };
             concatenated map = { 0, 1, 3, 4, 2 ,5 };
    */
  static Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>
  concatenateMaps(const std::vector<Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>>>& subMaps);

 private:
  bool CheckConsistency() const;

 private:
  RCP<const Map> fullmap_;
  std::vector<RCP<const Map>> maps_;
  std::vector<RCP<Import>> importers_;
  bool bThyraMode_;                        //< boolean flag: use Thyra numbering for local sub-block maps. default = false (for Xpetra mode)
  std::vector<RCP<const Map>> thyraMaps_;  //< store Thyra-style numbering maps here in Thyra mode. In Xpetra mode this vector is empty.
};                                         // BlockedMap class

}  // namespace Xpetra

#define XPETRA_BLOCKEDMAP_SHORT

#endif /* PACKAGES_XPETRA_SUP_BLOCKEDMAP_XPETRA_BLOCKEDMAP_DECL_HPP_ */
