// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_STKConnManager_hpp__
#define __Panzer_STKConnManager_hpp__

#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Kokkos includes
#include "Kokkos_DynRankView.hpp"

// Panzer includes
#include "Panzer_ConnManager.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"

namespace panzer_stk {

/** \brief A panzer::ConnManager implementation backed by an STK mesh,
  * mapping local element IDs to global mesh-entity (node/edge/face/cell)
  * IDs according to a given field pattern.
  *
  * <b>Connectivity caching</b>
  *
  * When caching is enabled with cacheConnectivity(), the connectivity built
  * by buildConnectivity() is saved and reused by later calls on any
  * STKConnManager. A cached result is reused only when all of the following
  * match: the STK_Interface object (compared by address), the field pattern,
  * and the list of sidesets passed to associateElementsInSideset() (in the
  * same order).
  *
  * - The cache holds a reference to every mesh it has seen, so those meshes
  *   are not freed until their cache entries are removed with
  *   clearCachedConnectivityData(const Teuchos::RCP<const STK_Interface>&)
  *   or clearCachedConnectivityData().
  * - Cached connectivity data is shared, not copied, between the cache and
  *   every STKConnManager that reuses it. For this reason all accessors
  *   return read-only data.
  * - The cache cannot detect changes made to a mesh after its connectivity
  *   was cached (e.g. rebalancing or refinement). Reusing a cached entry for
  *   a modified mesh returns stale connectivity, so clear the cached data for
  *   a mesh whenever that mesh changes.
  */
class STKConnManager : public panzer::ConnManager {
public:
   typedef typename panzer::ConnManager::LocalOrdinal LocalOrdinal;
   typedef typename panzer::ConnManager::GlobalOrdinal GlobalOrdinal;
   typedef typename Kokkos::DynRankView<GlobalOrdinal,PHX::Device>::host_mirror_type::const_type GlobalOrdinalView;
   typedef typename Kokkos::DynRankView<LocalOrdinal, PHX::Device>::host_mirror_type::const_type LocalOrdinalView;

   /// \brief Construct from the STK mesh to build connectivity from. buildConnectivity() must be called before connectivity queries are valid.
   STKConnManager(const Teuchos::RCP<const STK_Interface> & stkMeshDB);

   /// \brief Copy constructor.
   STKConnManager(const panzer_stk::STKConnManager & cm);

   virtual ~STKConnManager() {}

   /** Tell the connection manager to build the connectivity assuming
     * a particular field pattern.
     *
     * \param[in] fp Field pattern to build connectivity for
     */
   virtual void buildConnectivity(const panzer::FieldPattern & fp);

   /** Build a clone of this connection manager, without any assumptions
     * about the required connectivity (e.g. <code>buildConnectivity</code>
     * has never been called).
     */
   virtual Teuchos::RCP<panzer::ConnManager> noConnectivityClone() const;

   /** Get ID connectivity for a particular element
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Pointer to beginning of indices, with total size
     *          equal to <code>getConnectivitySize(localElmtId)</code>
     */
   virtual const panzer::GlobalOrdinal * getConnectivity(LocalOrdinal localElmtId) const
   { return &(*connectivityPtr_)[(*elmtLidToConnPtr_)[localElmtId]]; }

   /** How many mesh IDs are associated with this element?
     *
     * \param[in] localElmtId Local element ID
     *
     * \returns Number of mesh IDs that are associated with this element.
     */
   virtual LocalOrdinal getConnectivitySize(LocalOrdinal localElmtId) const
   { return (*connSizePtr_)[localElmtId]; }

   /// \brief Returns a read-only view of the flat connectivity array (all elements' global IDs, concatenated).
   const GlobalOrdinalView getConnectivityView() const
   { return GlobalOrdinalView(connectivityPtr_->data(), connectivityPtr_->size()); }

   /// \brief Returns a read-only view of the per-element connectivity size array.
   const LocalOrdinalView getConnectivitySizeView() const
   { return LocalOrdinalView(connSizePtr_->data(), connSizePtr_->size()); }

   /// \brief Returns a read-only view of the per-element offset into the flat connectivity array returned by getConnectivityView().
   const LocalOrdinalView getElementLidToConnView() const
   { return LocalOrdinalView(elmtLidToConnPtr_->data(), elmtLidToConnPtr_->size()); }

   /** Get the block ID for a particular element.
     *
     * \param[in] localElmtId Local element ID
     */
   virtual std::string getBlockId(LocalOrdinal localElmtId) const;

   /** How many element blocks in this mesh?
     */
   virtual std::size_t numElementBlocks() const
   { return stkMeshDB_->getNumElementBlocks(); }

   /** Get block IDs from STK mesh object
     */
   virtual void getElementBlockIds(std::vector<std::string> & elementBlockIds) const
   { return stkMeshDB_->getElementBlockNames(elementBlockIds); }

   /** What are the cellTopologies linked to element blocks in this connection manager?
    */
   virtual void getElementBlockTopologies(std::vector<shards::CellTopology> & elementBlockTopologies) const{
     std::vector<std::string> elementBlockIds;
     getElementBlockIds(elementBlockIds);
     elementBlockTopologies.reserve(elementBlockIds.size());

     for (unsigned i=0; i<elementBlockIds.size(); ++i) {
       elementBlockTopologies.push_back(*(stkMeshDB_->getCellTopology(elementBlockIds[i])));
     }
   }
   /** Get the local element IDs for a paricular element
     * block. These are only the owned element ids.
     *
     * \param[in] blockIndex Block Index
     *
     * \returns Vector of local element IDs.
     */
   virtual const std::vector<LocalOrdinal> & getElementBlock(const std::string & blockId) const
   { return *(elementBlocks_.find(blockId)->second); }

   /** Get the local element IDs for a paricular element
     * block. These element ids are not owned, and the element
     * will live on another processor.
     *
     * \param[in] blockIndex Block Index
     *
     * \returns Vector of local element IDs.
     */
   virtual const std::vector<LocalOrdinal> & getNeighborElementBlock(const std::string & blockId) const
   { return *(neighborElementBlocks_.find(blockId)->second); }

   /** Get the coordinates (with local cell ids) for a specified element block and field pattern.
     *
     * \param[in] blockId Block containing the cells
     * \param[in] coordProvider Field pattern that builds the coordinates
     * \param[out] localCellIds Local cell Ids (indices)
     * \param[out] Resizable field container that contains the coordinates
     *             of the points on exit.
     */
   virtual void getDofCoords(const std::string & blockId,
                             const panzer::Intrepid2FieldPattern & coordProvider,
                             std::vector<std::size_t> & localCellIds,
                             Kokkos::DynRankView<double,PHX::Device> & points) const;

    /** Get STK interface that this connection manager is built on.
      */
    Teuchos::RCP<const STK_Interface> getSTKInterface() const
    { return stkMeshDB_; }

    /** How many elements are owned by this processor. Further,
      * the ordering of the local ids is suct that the first
      * <code>getOwnedElementCount()</code> elements are owned
      * by this processor. This is true only because of the
      * local element ids generated by the <code>STK_Interface</code>
      * object.
      */
    std::size_t getOwnedElementCount() const
    { return ownedElementCount_; }

    /** Before calling buildConnectivity, provide sideset IDs from which to
      * extract associated elements.
      */
    void associateElementsInSideset(const std::string sideset_id);

    /** After calling <code>buildConnectivity</code>, optionally check which
      * sidesets yielded no element associations in this communicator. This is a
      * parallel operation. In many applications, the outcome indicating
      * correctness is that the returned vector is empty.
      */
    std::vector<std::string> checkAssociateElementsInSidesets(const Teuchos::Comm<int>& comm) const;

    /** Get elements, if any, associated with <code>el</code>, excluding
      * <code>el</code> itself.
      */
    virtual const std::vector<LocalOrdinal>& getAssociatedNeighbors(const LocalOrdinal& el) const;

    /** Return whether getAssociatedNeighbors will return true for at least one
      * input. Default implementation returns false.
      */
    virtual bool hasAssociatedNeighbors() const;

  /** \brief Enables the caching of connectivity data.
    *
    * Cached data keeps its meshes alive until it is cleared. Be sure to call
    * clearCachedConnectivityData() before exiting your program, and see the
    * class documentation for when cached data must be cleared.
    */
  static void cacheConnectivity()
  { cache_connectivity_ = true; }

  /// If you enable caching with cacheConnectivity(), this must be called prior to exiting the code.
  static void clearCachedConnectivityData()
  {
    PANZER_FUNC_TIME_MONITOR("panzer::ConnectivityManager::clearCachedConnectivityData()");
    cached_conn_managers_.clear();
  }

  /** \brief Removes all cached connectivity data built on the given mesh.
    *
    * This releases the cache's references to that mesh. Call it whenever the
    * mesh is modified or is no longer needed.
    */
  static void clearCachedConnectivityData(const Teuchos::RCP<const STK_Interface>& mesh);

  /// This is purely for unit testing. Returns the number of times that buildConnectivity() was called, but a cached version was found to use instead.
  static int getCachedReuseCount()
  { return cached_reuse_count_; }

protected:
   /** Apply periodic boundary conditions associated with the mesh object.
     *
     * \note This function requires global All-2-All communication IFF
     *       periodic boundary conditions are required.
     */
   void applyPeriodicBCs( const panzer::FieldPattern & fp, GlobalOrdinal nodeOffset, GlobalOrdinal edgeOffset,
                                                           GlobalOrdinal faceOffset, GlobalOrdinal cellOffset);
   void applyInterfaceConditions();

   void buildLocalElementMapping();
   void clearLocalElementMapping();
   void buildOffsetsAndIdCounts(const panzer::FieldPattern & fp,
                                LocalOrdinal & nodeIdCnt, LocalOrdinal & edgeIdCnt,
                                LocalOrdinal & faceIdCnt, LocalOrdinal & cellIdCnt,
                                GlobalOrdinal & nodeOffset, GlobalOrdinal & edgeOffset,
                                GlobalOrdinal & faceOffset, GlobalOrdinal & cellOffset) const;

   /**
    * @brief Loops over relations of a given rank for a specified element and adds a unique ID to the connectivity vector
    *
    * @param[in] element Mesh element
    * @param[in] subcellRank Rank of the subcell entities to identify
    * @param[in] idCnt Number of IDs on the requested subcell type
    * @param[in] offset Offset for requested subcell type
    * @param[in,optional] maxIds If positive, maximum number of IDs to connect to this element. If 0 (default), add all IDs.
    *
    * @pre Should call buildOffsetsAndIdCounts() to obtain \p idCnt and \p offset.
    * @note The connectivity manager needs only the lowest order nodal information.
    *       Hence, maxIds should be set appropriately if the STK mesh is second order or higher.
   */
   LocalOrdinal addSubcellConnectivities(stk::mesh::Entity element,unsigned subcellRank,
                                         LocalOrdinal idCnt,GlobalOrdinal offset,const unsigned maxIds=0);

   void modifySubcellConnectivities(const panzer::FieldPattern & fp, stk::mesh::Entity element,
                                    unsigned subcellRank,unsigned subcellId,GlobalOrdinal newId,GlobalOrdinal offset);

   Teuchos::RCP<const STK_Interface> stkMeshDB_;

   Teuchos::RCP<std::vector<stk::mesh::Entity> > elements_;

   // element block information
   std::map<std::string,Teuchos::RCP<std::vector<LocalOrdinal> > > elementBlocks_;
   std::map<std::string,Teuchos::RCP<std::vector<LocalOrdinal> > > neighborElementBlocks_;

   // Large arrays are held by RCP so cached connectivity can be shared
   // instead of copied. A rebuild must allocate new arrays rather than
   // modify these in place, since they may be shared with the cache.
   Teuchos::RCP<std::vector<LocalOrdinal> > elmtLidToConnPtr_; // element LID to Connectivity map
   Teuchos::RCP<std::vector<LocalOrdinal> > connSizePtr_; // element LID to Connectivity size
   Teuchos::RCP<std::vector<GlobalOrdinal> > connectivityPtr_; // Connectivity

   std::size_t ownedElementCount_;

   std::vector<std::string> sidesetsToAssociate_;
   std::vector<bool> sidesetYieldedAssociations_;
   Teuchos::RCP<std::vector<std::vector<LocalOrdinal> > > elmtToAssociatedElmtsPtr_;

  using CachedEntry = std::pair<Teuchos::RCP<const panzer::FieldPattern>,Teuchos::RCP<panzer_stk::STKConnManager>>;
  /// Matches a cache entry on mesh, field pattern and sideset association list.
  struct CacheKeyCompare {
    const panzer::FieldPattern& fp_;
    const STK_Interface* mesh_;
    const std::vector<std::string>& sidesets_;
    CacheKeyCompare(const panzer::FieldPattern& fp,
                    const STK_Interface* mesh,
                    const std::vector<std::string>& sidesets)
      : fp_(fp), mesh_(mesh), sidesets_(sidesets) {}
    bool inline operator()(const CachedEntry& entry_to_compare) const
    {
      const STKConnManager& cm = *entry_to_compare.second;
      return (cm.stkMeshDB_.get() == mesh_) &&
             (cm.sidesetsToAssociate_ == sidesets_) &&
             fp_.equals(*entry_to_compare.first);
    }
  };

  static bool cache_connectivity_;
  static std::vector<CachedEntry> cached_conn_managers_;
  static int cached_reuse_count_;
};

}

#endif
