#ifndef XPETRA_REGION_MANAGER_DECL_HPP_
#define XPETRA_REGION_MANAGER_DECL_HPP_

// headers


// forward declarations
namespace Teuchos
{
  class Comm<std::size_t>;
}

namespace Xpetra {

/*! \class SplittingMapsInfo
 *
 *  This is an auxiliary class to store row maps for the composite matrix,
 *  region matrices and a regionToAll map to link region node indices with the
 *  composite ones.
 *
 *  \todo ToDo (mayr.mt) Add set/get access routines to control access to class members
 *
 *  ToDo (mayr.mt) Is this class obsolete? It only stores three variables,
 *  but everyting else (like printing the content of this variables) is implemented in Xpetra::RegionManager.
 *  Consider to move these variables to the Xpetra::RegionManager to drop this class entirely.
 *
 *  \author mayr.mt \date 09/2017
 */
template<class GO>
class SplittingMapsInfo
{

public:

  //! @name Construction/Destruction
  //@{

  //! Constructor
  SplittingMapsInfo();

  //! Destructor
  virtual ~SplittingMapsInfo(){};

  //@}

  //! @name Maps
  //@{

  //! used as a map for a RegionToAll node index
  Teuchos::Array<Teuchos::Array<std::tuple<GO,GO> > > regionToAll_;

  //! RowMap of composite matrix //ToDo (mayr.mt) Should this be replaced by an actual Xpetra::Map?
  Teuchos::Array<GO> compositeMap_;

  //! Collection of RowMaps for region matrices //ToDo (mayr.mt) Should this be replaced by an actual Xpetra::BlockedMap?
  Teuchos::Array<Teuchos::Array<GO> > regionMaps_;

  //@}

};

/*! \class RegionManager
 *
 *  \author mayr.mt \date 09/2017
 */
template<class SC, class LO, class GO, class NO>
class RegionManager
{
  public:

  //! @name Construction/Destruction
  //@{

  //! Constructor (empty)
  RegionManager();

  //! Constructor (read node-to-region mapping from file)
  RegionManager(const std::string& mappingFileName, ///< file name of file with nodes-to-region mapping information
      Teuchos::RCP<const Teuchos::Comm<int> > comm ///< communicator
      );

  //! Destructor
  virtual ~RegionManager(){};

  //@}

  //! @name Access routines
  //@{

  //! Get total number of nodes in the mesh
  GO getNumNodes() const { return numNodes_; }

  //! Get total number of regions in the mesh
  GO getNumRegions() const { return numRegions_; }

  //! Get number of nodes in region \c regionID
  GO getNumNodesPerRegion(GO regionID) const { return numNodesPerRegion_[regionID]; }

  //! Get row map of composite matrix
  Teuchos::Array<GO> getCompositeRowMap() const { return maps_->compositeMap_; }

  //! Get row map of region \c regionID
  Teuchos::Array<GO> getRegionRowMap(GO regionID) const;

  Teuchos::Array<Teuchos::Array<GO> > getRegionRowMaps() const { return maps_->regionMaps_; };

  //used as a map for a RegionToAll node index
  Teuchos::Array<Teuchos::Array< std::tuple<GO,GO> > > getRegionToAll() const { return maps_->regionToAll_; }

  // used as a map for a RegionToAll node index
  Teuchos::Array< std::tuple<GO,GO> > getRegionToAll(GO) const;

  //! Access mapping of interface nodes to regions
  Teuchos::Array<std::tuple<int,Teuchos::Array<GO> > > getInterfaceNodesToRegions() const { return interfaceNodesToRegions_; }

  //@}

  //! @name Print routines
  //@{

  //! Print all node-region pairs
  virtual void printNodeRegionPairs(Teuchos::FancyOStream& out ///< output stream
      ) const;

  //! Print mapping of nodes to regions
  virtual void printNodesToRegionMapping(Teuchos::FancyOStream& out ///< output stream
      ) const;

  //! Print mapping of interface nodes to regions
  virtual void printInterfaceNodesToRegionMapping(Teuchos::FancyOStream& out ///< output stream
      ) const;

  /*! \brief Print inactive processors
   *
   *  Inactive processors are those, that do not hold a row in the composite matrix
   */
  virtual void printInactiveProcs(Teuchos::FancyOStream& out ///< output stream
      ) const;

  //! \brief Print number of regions associated with each processor
  virtual void printNumRegionsPerProc(Teuchos::FancyOStream& out ///< output stream
      ) const;

  //! \brief Print all processors for each region
    virtual void printProcsPerRegion(Teuchos::FancyOStream& out ///< output stream
        ) const;

  //! \brief Print all regions that are (at least partially) owned by a processor
  virtual void printRegionsPerProc(Teuchos::FancyOStream& out ///< output stream
      ) const;

  //@}

  protected:

  private:

  //! @name Input
  //@{

  /*! \brief Read node-to-region mapping data from file
   *
   *  The information contained in the file is imported and stored in a
   *  Teuchos::Array of tuples<compostiteNodeIndex, regionIndex>.
   *  The first field of the tuple is the composite node index,
   *  the second field of the tuple is the region index.
   */
  virtual void readMappingFromFile(const std::string& mappingFileName ///< file name of file with node-to-region mapping data
      );

  //@}

  //! Communicator
  Teuchos::RCP<const Teuchos::Comm<int> > comm_;

  //! @name Mesh/Discretization data
  //@{

  /*! \brief Setup the mapping of nodes to regions
   *
   *  This routine associates a globally indexed node with the list of regions it belongs to
   *  This is helpful to spot which nodes lie on a interregion interface. In fact, these nodes must have
   *  the list of regions with more than one element.
   *
   *  First, create the mapping #nodesToRegions_, but also setup the interface
   *  mapping #interfaceNodesToRegions_.
   */
  virtual void setupNodeToRegionMapping();

  //! Count the number of nodes in each region
  virtual void countNodesPerRegion();

  //! Number of nodes in mesh
  GO numNodes_;

  //! Number of regions in mesh
  GO numRegions_;

  //! Number of nodes in each region
  Teuchos::Array<GO> numNodesPerRegion_;

  /*! \brief Node-Region pairings (read from file)
   *
   *  We collect each node-region pair in a std::tuple. If a node resides on an
   *  interface between two regions, it is part of two tuples (one for each
   *  region it is attached to).
   *
   *  The tuple is organized as follows:
   *  - first entry: global node ID
   *  - second entry: global region ID
   *
   *  After reading for file, we sort this list according to the region index.
   *
   *  \sa sortNodesByRegions()
   *  \sa nodesSortedByRegions_
   */
  Teuchos::Array<std::tuple<GO,GO> > nodeRegionPairs_;

  /*! \brief Mapping of nodes to regions
   *
   *  Interior nodes belong only to one region, while nodes on interfaces
   *  between regions belong to multiple regions.
   *
   *  The tuple is organized as follows:
   *  - first entry: global node ID
   *  - second entry: list of regions this node belongs to.
   *
   *  \note Size of #nodesToRegions_ always matches #numNodes_.
   */
  Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > > nodesToRegions_;

  /*! \brief Mapping of interface nodes to regions
   *
   *  Collect regions IDs of those regions, that an interface node belongs to.
   *
   *  The tuple is organized as follows:
   *  - first entry: global node ID
   *  - second entry: list of regions this node belongs to.
   */
  Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > > interfaceNodesToRegions_;

  //@}

  //! Ordering of node and region lists
  //@{

  /*! \brief Sort the list of nodes by ascending regions IDs
   *
   *  The list of node-region pairs #nodeRegionPairs_ is sorted by acending
   *  region IDs. After sorting, #nodesSortedByRegions_ is set to \c true.
   */
  virtual void sortNodesByRegions();

  /*! \brief Have the nodes been sorted by ascending region IDs?
   *
   *  This flag indicates whether the sorting has been performed.
   *
   *  \sa sortNodesByRegions()
   */
  bool nodesSortedByRegions_;

  //@}

  /*! @name Parallel distribuation and processes
   *
   *  We allow for cases, where the number of processors is smaller, equal, or
   *  larger than the number of regions.
   */
  //@{

  /*! \brief Setup the mapping of processors to regions
   *
   *  This routines computes the way regions are partitioned across processes
   *  The partitioning policies of course depends on whether the number of processes
   *  exceeds the number of regions or not.
   *
   *  ASSUMPTION: A PROCESS CANNOT OWN CHUNKS OF MULTIPLE REGIONS. EITHER A PROCESS IS CONFINED INSIDE A SINGLE REGION
   *  OR IT MUST POSSESS ENTIRE REGIONS.
   * The distribution of regions (or portions of them) across processes is conducted such that load balancing is guaranteed.
   */
  virtual void setupProcToRegionMapping();

  /*! \brief Regions per processor
   *
   *  Listing the processors that share a certain region
   *
   *  If number of processors is larger than #numRegions_, this lists the number
   *  of regions per processor.
   *
   *  Otherwise, #regionsPerProc_ is empty.
   */
  Teuchos::Array<GO> regionsPerProc_; // ToDo (mayr.mt) Do we really need this? Implementation seems to be wrong?

  Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > > regionsPerProc2_;

  //! Lists of processes instantiated for each region
  Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > > procsPerRegion_; // ToDo (mayr.mt) Do we really need this? Implementation seems to be wrong?

  //@}

  //! @name Maps
  //@{

  //! Create the row maps
  virtual void setupRowMaps(); // ToDo (mayr.mt) split this into subroutines

  //! Maps used for composite and region operators
  Teuchos::RCP<Xpetra::SplittingMapsInfo<GO> > maps_;

  //@}

};

} // namespace Xpetra

#endif /* XPETRA_REGION_MANAGER_DECL_HPP_ */
