#ifndef XPETRA_REGION_MANAGER_IMPL_HPP_
#define XPETRA_REGION_MANAGER_IMPL_HPP_

// Teuchos
#include <Teuchos_Comm.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_RCP.hpp>

// Xpetra
#include "Xpetra_IO.hpp"
#include "Xpetra_RegionManager_decl.hpp"
#include "Xpetra_RegionUtils_decl.hpp"


template<class GO>
Xpetra::RegionNode<GO>::RegionNode(const GO nodeID,
    const Teuchos::Array<GO> regions, const GO procID)
    : nodeID_(nodeID),
      regions_(regions),
      proc_(procID)
{
  return;
}

template<class GO>
Xpetra::RegionNode<GO>::RegionNode(const GO nodeID,
    const GO regionID, const GO procID)
    : nodeID_(nodeID),
      proc_(procID)
{
  regions_.push_back(regionID);

  return;
}

template<class GO>
void Xpetra::RegionNode<GO>::appendRegionID(const GO regionID)
{
  regions_.push_back(regionID);

  return;
}

template<class GO>
const bool Xpetra::RegionNode<GO>::isInterfaceNode() const
{
  if (regions_.size() > 1)
    return true;
  else
    return false;
}

template<class GO>
Xpetra::RegionNodes<GO>::RegionNodes(const GO numNodes)
{
  // Initialize array to proper length and set Teuchos::null for every entry
  nodes_.clear();
  nodes_.resize(numNodes);
  typename Teuchos::Array<Teuchos::RCP<Xpetra::RegionNode<GO> > >::iterator it;
  for (it = nodes_.begin(); it < nodes_.end(); ++it)
    *it = Teuchos::null;

  nodesPerRegion_.clear();

  return;
}

template<class GO>
void Xpetra::RegionNodes<GO>::addNode(Teuchos::RCP<Xpetra::RegionNode<GO> > node)
{
  const GO nodeID = node->getNodeID();

  // Add node if not set yet
  if (nodes_[nodeID].is_null())
  {
    nodes_[nodeID] = node;
  }
  else
  {
    // modify existing node by adding region ID to its list of region IDs
//    nodes_[nodeID]->appendRegionID(node->getRegionID());
  }

  return;
}

template<class GO>
void Xpetra::RegionNodes<GO>::setupMappingNodesPerRegion(const GO numRegions)
{
  nodesPerRegion_.clear();
  nodesPerRegion_.resize(numRegions);

  // Loop over all nodes and extract sort them according to their region index.
  typename Teuchos::Array<Teuchos::RCP<Xpetra::RegionNode<GO> > >::const_iterator it;
  for (it = nodes_.begin(); it < nodes_.end(); ++it)
  {
    const Teuchos::Array<GO>& regIDs = (*it)->getRegions();

    typename Teuchos::Array<GO>::const_iterator regID;
    for (regID = regIDs.begin(); regID < regIDs.end(); ++ regID)
    {
      nodesPerRegion_[*regID].push_back((*it)->getNodeID());
    }
  }

  return;
}

template<class GO>
const bool Xpetra::RegionNodes<GO>::hasNode(const GO nodeID) const
{
  if (nodes_.size() <= nodeID)
    return false;

  if (nodes_[nodeID].is_null())
    return false;

  if (nodes_[nodeID].is_valid_ptr())
    return true;

  TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Logic failed. This is a case that hasn't been detected properly.");
  return false;
}

template<class GO>
const Xpetra::RegionNode<GO>& Xpetra::RegionNodes<GO>::getNode(const GO globalNodeID) const
{
  return *nodes_[globalNodeID];
}

template<class GO>
const GO Xpetra::RegionNodes<GO>::getNumNodes() const
{
  return nodes_.size();
}

template<class GO>
const GO Xpetra::RegionNodes<GO>::getNumNodesPerRegion(const GO regionID) const
{
  return nodesPerRegion_[regionID].size();
}

template<class GO>
Xpetra::RegionNode<GO>& Xpetra::RegionNodes<GO>::modifyNode(const GO nodeID)
{
  return *nodes_[nodeID];
}

template<class GO>
void Xpetra::RegionNodes<GO>::printAllNodes(Teuchos::FancyOStream& out) const
{
  out << std::endl << "nodes_:" << std::endl << "Node\tProc\tRegions" << std::endl;

  typename Teuchos::Array<Teuchos::RCP<Xpetra::RegionNode<GO> > >::const_iterator it;
  for (it = nodes_.begin(); it < nodes_.end(); ++it)
  {
    if (not (*it).is_null())
      out << (*it)->getNodeID() << "\t" << (*it)->getProc() << "\t" << (*it)->getRegions() << std::endl;
//    else
//      out << "Node " << it << " not defined, yet." << std::endl;
  }

  return;
}

template<class GO>
void Xpetra::RegionNodes<GO>::printRegionData(Teuchos::FancyOStream& out) const
{
  out << std::endl << "nodesPerRegion_: " << std::endl << "Region\tNodes" << std::endl;

  int regionIndex = 0;
  typename Teuchos::Array<Teuchos::Array<GO> >::const_iterator it;
  for (it = nodesPerRegion_.begin(); it < nodesPerRegion_.end(); ++it, ++regionIndex)
  {
      out << regionIndex << "\t" << *it << std::endl;
  }

  return;
}

template<class GO>
Teuchos::RCP<const Teuchos::Array<GO> > Xpetra::RegionNodes<GO>::getNodeGIDsPerProc(
    const int myRank) const
{
  Teuchos::RCP<Teuchos::Array<GO> > myNodeGIDs = Teuchos::rcp(new Teuchos::Array<GO>());
  myNodeGIDs->clear();

  typename Teuchos::Array<Teuchos::RCP<Xpetra::RegionNode<GO> > >::const_iterator it;
  for (it = nodes_.begin(); it < nodes_.end(); ++it)
  {
    if ((*it)->getProc() == myRank)
      myNodeGIDs->push_back((*it)->getNodeID());
  }

  return myNodeGIDs;
}

template<class GO>
Teuchos::RCP<const Teuchos::Array<GO> > Xpetra::RegionNodes<GO>::getNodeGIDsPerRegion(
    const GO regionID) const
{
  Teuchos::RCP<Teuchos::Array<GO> > myNodeGIDs = Teuchos::rcp(new Teuchos::Array<GO>());
  myNodeGIDs->clear();

  typename Teuchos::Array<GO>::const_iterator it;
  for (it = nodesPerRegion_[regionID].begin(); it < nodesPerRegion_[regionID].end(); ++it)
    myNodeGIDs->push_back(*it);

  std::cout << "Region " << regionID << ": " << *myNodeGIDs << std::endl;

  return myNodeGIDs;
}

template<class GO>
Teuchos::RCP<const Teuchos::Array<GO> > Xpetra::RegionNodes<GO>::getNodeGIDsPerRegionAndProc(
    const GO regionID, const int myRank) const
{
  Teuchos::RCP<Teuchos::Array<GO> > myNodeGIDs = Teuchos::rcp(new Teuchos::Array<GO>());
  myNodeGIDs->clear();

  typename Teuchos::Array<GO>::const_iterator it;
  for (it = nodesPerRegion_[regionID].begin(); it < nodesPerRegion_[regionID].end(); ++it)
  {
    const Xpetra::RegionNode<GO>& node = *nodes_[*it];

    if (node.getProc() == myRank)
      myNodeGIDs->push_back(*it);
  }

  std::cout << "Region/Proc " << regionID << "/" << myRank << ": " << *myNodeGIDs << std::endl;

  return myNodeGIDs;
}

template<class GO>
Teuchos::Array<std::tuple<int,Teuchos::Array<GO> > > Xpetra::RegionNodes<GO>::getMappingInterfaceNodesToRegions() const
{
  Teuchos::Array<std::tuple<int,Teuchos::Array<GO> > > intNodesToRegions;
  intNodesToRegions.clear();

  for (GO i = 0; i < nodes_.size(); ++i)
  {
    if (nodes_[i]->isInterfaceNode())
      intNodesToRegions.push_back(std::make_tuple(nodes_[i]->getNodeID(), nodes_[i]->getRegions()));
  }

  return intNodesToRegions;
}

template<class SC, class LO, class GO, class NO>
Xpetra::RegionManager<SC,LO,GO,NO>::RegionManager(
    const std::string& mappingFileName,
    Teuchos::RCP<const Teuchos::Comm<int> > comm)
    : nodes_(Teuchos::null),
      comm_(comm),
      compositeMap_(Teuchos::null)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  if (comm_->getRank() == 0)
    *out << "Constructing 'RegionManager'." << std::endl;

  // Read node-to-region mapping from file
  readMappingFromFile(mappingFileName);
  if (comm_->getRank() == 0)
    nodes_->printAllNodes(*out);

  setupMappingNodesPerRegion();
  if (comm_->getRank() == 0)
    nodes_->printRegionData(*out);

  setupRowMaps();

  setupMappingNodesPerRegion();

  if (comm_->getRank() == 0)
    nodes_->printRegionData(*out);

  setupRowMaps();

  setupMappingNodesPerRegion();

  if (comm_->getRank() == 0)
    nodes_->printRegionData(*out);

  setupRowMaps();

  return;
}

template<class SC, class LO, class GO, class NO>
void Xpetra::RegionManager<SC,LO,GO,NO>::readMappingFromFile(
    const std::string& mappingFileName)
{
  std::ifstream inputFile(mappingFileName, std::ifstream::in);
  std::string line;
  TEUCHOS_TEST_FOR_EXCEPTION(!inputFile.good(), Exceptions::RuntimeError,
      "Cannot read file '" << mappingFileName << "'.");

  GO lineIndex = 0;
  GO aux;
  Teuchos::Array<GO> lineContent;

  // Loop over all lines in file
  while (std::getline(inputFile, line))
  {
    std::istringstream is(line);
    lineContent.clear();

    switch (lineIndex)
    {
    case 0:
    case 2:
    {
      break;
    }
    case 1:
    {
      // parse the current line and extract the global number of nodes in the mesh
      while (is >> aux)
        lineContent.push_back(aux);

      numNodes_ = lineContent[0];
      numRegions_ = lineContent[1];

      // setup the nodes_ object with the appropriate number of nodes
      nodes_ = Teuchos::rcp(new Xpetra::RegionNodes<GO>(numNodes_));

      break;
    }
    default:
    {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(nodes_.is_null(), "'nodes_' has not been initialized, yet.");

      while (is >> aux)
        lineContent.push_back(aux);

      const GO nodeID = lineContent[0];
      const GO regionID = lineContent[1];
      const GO procID = lineContent[2];

      if (not nodes_->hasNode(nodeID))
      {
        // Node does not exist yet. Create it!
        Teuchos::RCP<Xpetra::RegionNode<GO> > node =
            Teuchos::rcp(new Xpetra::RegionNode<GO>(nodeID, regionID, procID));
        nodes_->addNode(node);
      }
      else
      {
        // Node already exists. Just add region ID to its list of region IDs.
        nodes_->modifyNode(nodeID).appendRegionID(regionID);
      }

      break;
    }
    }

    ++lineIndex;
  }

  // close the file
  inputFile.close();

  return;
}

template<class SC, class LO, class GO, class NO>
void Xpetra::RegionManager<SC,LO,GO,NO>::setupMappingNodesPerRegion()
{
  nodes_->setupMappingNodesPerRegion(numRegions_);

  for (GO regID = 0; regID < numRegions_; ++regID)
  {

  }

  return;
}

template<class SC, class LO, class GO, class NO>
void Xpetra::RegionManager<SC,LO,GO,NO>::setupRowMaps()
{
  setupCompositeRowMap();
  setupRegionRowMaps();

  return;
}

template<class SC, class LO, class GO, class NO>
void Xpetra::RegionManager<SC,LO,GO,NO>::setupCompositeRowMap()
{
  /* Access list of GIDs owned by each proc. Conversion of RCP<Array> to ArrayView
   * is not done in the most elegant way, but it works.
   *
   * Note: operator() creates and ArrayView of the Array myNodesGID.
   */
  Teuchos::RCP<const Teuchos::Array<GO> > myNodesGIDsRcp = nodes_->getNodeGIDsPerProc(comm_->getRank());
  const Teuchos::Array<GO>& myNodesGIDs = *myNodesGIDsRcp;

  compositeMap_ = Xpetra::MapFactory<LO,GO,NO>::Build(Xpetra::UseTpetra, nodes_->getNumNodes(), myNodesGIDs(), 0, comm_);

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  *out << std::endl << "compositeMap_:" << std::endl;
  compositeMap_->describe(*out, Teuchos::VERB_EXTREME);

  return;
}

template<class SC, class LO, class GO, class NO>
void Xpetra::RegionManager<SC,LO,GO,NO>::setupRegionRowMaps()
{
  comm_->barrier();

  regionMaps_.clear();
  regionMaps_.resize(numRegions_);

  // loop over all regions and create map for each region
  for (GO i = 0; i < numRegions_; ++i)
  {
    /* Access list of GIDs owned by each proc. Conversion of RCP<Array> to ArrayView
     * is not done in the most elegant way, but it works.
     *
     * Note: operator() creates and ArrayView of the Array myNodesGID.
     */
    Teuchos::RCP<const Teuchos::Array<GO> > myNodesGIDsRcp = nodes_->getNodeGIDsPerRegionAndProc(i, comm_->getRank());
    const Teuchos::Array<GO>& myNodesGIDs = *myNodesGIDsRcp;

    regionMaps_[i] = Xpetra::MapFactory<LO,GO,NO>::Build(Xpetra::UseTpetra, nodes_->getNumNodesPerRegion(i), myNodesGIDs(), 0, comm_);
  }

  comm_->barrier();

  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  for (GO i = 0; i < numRegions_; ++i)
  {
    *out << std::endl << "regionMaps[" << i << "]_:" << std::endl;
    regionMaps_[i]->describe(*out, Teuchos::VERB_EXTREME);
  }

  return;
}

template<class SC, class LO, class GO, class NO>
const GO Xpetra::RegionManager<SC,LO,GO,NO>::getNumNodesPerRegion(
    const GO regID) const
{
  return regionMaps_[regID]->getGlobalNumElements();
}

template<class SC, class LO, class GO, class NO>
Teuchos::Array<std::tuple<GO,GO> > Xpetra::RegionManager<SC,LO,GO,NO>::getRegionToAll(
    const GO regID) const
{
  const Teuchos::Array<GO> nodesPerRegion = *nodes_->getNodeGIDsPerRegion(regID);

  Teuchos::Array<std::tuple<GO,GO> > regionToAll;
  regionToAll.clear();
  regionToAll.resize(nodes_->getNumNodesPerRegion(regID));

  for (int i = 0; i < nodesPerRegion.size(); ++i)
    regionToAll[i] = std::make_tuple(i, nodesPerRegion[i]);

  return regionToAll;
}

template<class SC, class LO, class GO, class NO>
Teuchos::Array<std::tuple<int,Teuchos::Array<GO> > > Xpetra::RegionManager<SC,
    LO,GO,NO>::getInterfaceNodesToRegions() const
{
  return nodes_->getMappingInterfaceNodesToRegions();
}

#endif /* XPETRA_REGION_MANAGER_IMPL_HPP_ */
