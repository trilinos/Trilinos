#ifndef XPETRA_REGION_MANAGER_DECL_HPP_
#define XPETRA_REGION_MANAGER_DECL_HPP_

// headers

// Xpetra
#include "Xpetra_Map.hpp"

// forward declarations
namespace Teuchos
{
  template<>
  class Comm<std::size_t>;
}

namespace Xpetra {

/*! \class RegionNode
 *  \brief Data of a single node in an HHG-based/region-wise setting
 *
 *  A node itself is a unique object, although it association with a regions is
 *  not unique. Nodes on region interfaces belong to multiple regions.
 *
 *  \author mayr.mt \date 10/2017
 */
template<class GO>
class RegionNode
{

public:

  //! @name Construction/Destruction
  //@{

  //! Constructor
  RegionNode(const GO nodeID,
      const Teuchos::Array<GO> regions,
      const GO procID
      );

  //! Constructor
  RegionNode(const GO nodeID,
      const GO regionID,
      const GO procID
      );

  //! Destructor
  virtual ~RegionNode(){};

  //@}

  //! Modify existing node
  //@{

  //! Append a region ID to list of region IDs
  virtual void appendRegionID(const GO regionID);

  //@}


  //! Access routines
  //@{

  //! get the node ID
  virtual const GO getNodeID() const { return nodeID_; }

  //! get list of regions this node belongs to
  virtual const Teuchos::Array<GO> getRegions() const { return regions_; }

  //! get number of regions this node belongs to
  virtual const GO getNumRegions() const { return regions_.size(); }

  //! get ID of node-owning process
  virtual const GO getProc() const { return proc_; }

  //@}

  //! Check status of this node
  //@[

  //! Is this an interface node?
  virtual const bool isInterfaceNode() const;

  //@}

private:

  //! Node ID
  GO nodeID_;

  //! Regions this node belongs to
  Teuchos::Array<GO> regions_;

  //! Parallel process that owns this node
  GO proc_;
};

/*! \class RegionNodes
 *  \brief Collection of all \c RegionNode objects in the entire mesh
 *
 *  \author mayr.mt \date 10/2017
 */
template<class GO>
class RegionNodes
{

public:

  //! @name Construction/Destruction
  //@{

  //! Constructor
  RegionNodes(const GO numNodes);

  //! Destructor
  virtual ~RegionNodes(){};

  //@}

  //! @name Setup
  //@{

  //! Add a new node
  virtual void addNode(Teuchos::RCP<Xpetra::RegionNode<GO> > node);

  //! Modify an existing node
  virtual Xpetra::RegionNode<GO>& modifyNode(const GO nodeID);

  //! Setup mapping of nodes to regions
  virtual void setupMappingNodesPerRegion(const GO numRegions ///< total number of regions
      );

  //@}

  //! Access routines
  //@{

  //! Does node with ID \c nodeID already exist?
  virtual const bool hasNode(const GO nodeID) const;

  //! get the node (read only)
  virtual const Xpetra::RegionNode<GO>& getNode(const GO globalNodeID ///< global node ID
      ) const;

  //! get global number of nodes
  virtual const GO getNumNodes() const;

  //! get number of nodes in region \c regionID
  virtual const GO getNumNodesPerRegion(const GO regionID ///< ID of the region of interest
      ) const;

  //! Get list of node GIDs per processor
  virtual Teuchos::RCP<const Teuchos::Array<GO> > getNodeGIDsPerProc(const int myRank ///< MPI rank
      ) const;

  //! Get list of node GIDs per region
  virtual Teuchos::RCP<const Teuchos::Array<GO> > getNodeGIDsPerRegion(const GO regionID ///< Region ID
      ) const;

  //! Get list of node GIDs per region and per processor
  virtual Teuchos::RCP<const Teuchos::Array<GO> > getNodeGIDsPerRegionAndProc(const GO regionID, ///< Region ID
      const int myRank ///< MPI rank
      ) const;

  //! Access mapping of interface nodes to regions
  virtual Teuchos::Array<std::tuple<int,Teuchos::Array<GO> > > getMappingInterfaceNodesToRegions() const;

  //@}

  //! Print routines
  //@{

  //! Print list of nodes
  virtual void printAllNodes(Teuchos::FancyOStream& out ///< output stream
      ) const;

  //! Print region data
  virtual void printRegionData(Teuchos::FancyOStream& out ///< output stream
      ) const;

  //@}

private:

  //! \brief Collection of all nodes (redundant on all procs)
  Teuchos::Array<Teuchos::RCP<Xpetra::RegionNode<GO> > > nodes_;

  //! \brief Lists of nodes per region (redundant on all procs)
  Teuchos::Array<Teuchos::Array<GO> > nodesPerRegion_;

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
  virtual GO getNumNodes() const { return numNodes_; }

  //! Get total number of nodes in the mesh
  virtual const GO getNumNodesPerRegion(const GO regID ///< region ID
      ) const;

  //! Get total number of regions in the mesh
  virtual GO getNumRegions() const { return numRegions_; }

  //! Get composite map
  virtual Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > getCompositeMap() const { return compositeMap_; }

  //! Get map of region \c regID
  virtual Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > getRegionMap(const GO regID) const { return regionMaps_[regID]; }

  /*! \brief Get mapping of local node ID to global node ID for region \c regID
   *
   *  \warning This relies on std::tuples and is a relict from Max's implementation.
   *
   *  \todo Replace this once the splitting/collapsing routines are refactored.
   */
  virtual Teuchos::Array<std::tuple<GO,GO> > getRegionToAll(const GO regID) const;

  //! Access mapping of interface nodes to regions
  virtual Teuchos::Array<std::tuple<int,Teuchos::Array<GO> > > getInterfaceNodesToRegions() const;

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
   *
   *  We assume the file header to consists of 3 row, while the actual listing
   *  of all nodes starts in line 3 (when the first line is line 0). Line 1
   *  contains global information like the total number of nodes, regions, and
   *  parallel processes.
   */
  virtual void readMappingFromFile(const std::string& mappingFileName ///< file name of file with node-to-region mapping data
      );

  //! List of all node objects
  Teuchos::RCP<Xpetra::RegionNodes<GO> > nodes_;

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
  virtual void setupMappingNodesPerRegion();

//  //! Count the number of nodes in each region
//  virtual void countNodesPerRegion();

  //! Number of nodes in mesh
  GO numNodes_;

  //! Number of regions in mesh
  GO numRegions_;

  //@}

  //! @name Maps
  //@{

  /*! \brief Create the row maps
   *
   *  \sa setupCompositeRowMap()
   *  \sa setupRegionRowMaps()
   */
  virtual void setupRowMaps();

  /*! \brief Create the row map of the composite matrix
   *
   *  \todo ToDo (mayr.mt) Make the conversion of Teuchos::Array nicer.
   */
  virtual void setupCompositeRowMap();

  //! Create the row map of the composite matrix
  virtual void setupRegionRowMaps();

  //! Composite map
  Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > compositeMap_;

  //! Collection of region maps
  Teuchos::Array<Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > > regionMaps_;

  //@}

};

} // namespace Xpetra

#endif /* XPETRA_REGION_MANAGER_DECL_HPP_ */
