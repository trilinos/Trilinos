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
Xpetra::SplittingMapsInfo<GO>::SplittingMapsInfo()
{
  return;
}

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
const bool Xpetra::RegionNode<GO>::IsInterfaceNode() const
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
const Xpetra::RegionNode<GO>& Xpetra::RegionNodes<GO>::getNode(const GO nodeID) const
{
  return *nodes_[nodeID];
}

template<class GO>
const GO Xpetra::RegionNodes<GO>::getNumNodes() const
{
  return nodes_.size();
}

template<class GO>
Xpetra::RegionNode<GO>& Xpetra::RegionNodes<GO>::modifyNode(const GO nodeID)
{
  return *nodes_[nodeID];
}

template<class GO>
void Xpetra::RegionNodes<GO>::printAllNodes(Teuchos::FancyOStream& out) const
{
  out << std::endl << "Node\tProc\tRegions" << std::endl;

  typename Teuchos::Array<Teuchos::RCP<Xpetra::RegionNode<GO> > >::const_iterator it;
  for (it = nodes_.begin(); it < nodes_.end(); ++it)
  {
    if (not (*it).is_null())
      out << (*it)->getNodeID() << "\t" << (*it)->getProc() << "\t" << (*it)->getRegions() << std::endl;
    else
      out << "Node " << it << " not defined, yet." << std::endl;
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

template<class SC, class LO, class GO, class NO>
Xpetra::RegionManager<SC,LO,GO,NO>::RegionManager(
    const std::string& mappingFileName,
    Teuchos::RCP<const Teuchos::Comm<int> > comm)
    : nodes_(Teuchos::null),
      comm_(comm),
      compositeMap_(Teuchos::null)
//      numNodes_(0),
//      numRegions_(0),
//      nodeRegionPairs_(0),
//      nodesToRegions_(0),
//      interfaceNodesToRegions_(0),
//      nodesSortedByRegions_(false),
//      maps_(Teuchos::null)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  if (comm_->getRank() == 0)
    *out << "Constructing 'RegionManager'." << std::endl;

  // Read node-to-region mapping from file
  readMappingFromFile(mappingFileName);
//  nodes_->printAllNodes(*out);

  setupMappingNodesPerRegion();
  setupRowMaps();

//  // Nodes are shuffled so that regions are sorted in ascending labeling order
//  sortNodesByRegions();
//
//  // create some more mappings
//  setupNodeToRegionMapping();
//  setupProcToRegionMapping();
//  setupRowMaps();
//
//  countNodesPerRegion();

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

//      std::cout << "Total number of nodes: " << numNodes << std::endl;

      // setup the nodes_ object with the appropriate number of nodes
      nodes_ = Teuchos::rcp(new Xpetra::RegionNodes<GO>(numNodes_));

      break;
    }
    default:
    {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(nodes_.is_null(), "'nodes_' has not been initialized, yet.");

//      std::cout << "Line " << lineIndex << ": " << line << std::endl;

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

//template<class SC, class LO, class GO, class NO>
//void Xpetra::RegionManager<SC,LO,GO,NO>::sortNodesByRegions()
//{
//  std::sort(nodeRegionPairs_.begin(), nodeRegionPairs_.end(), Xpetra::compareRegions<GO>);
//  nodesSortedByRegions_ = true;
//
//  return;
//}
//

template<class SC, class LO, class GO, class NO>
void Xpetra::RegionManager<SC,LO,GO,NO>::setupMappingNodesPerRegion()
{
  nodes_->setupMappingNodesPerRegion(numRegions_);

  return;
}
//template<class SC, class LO, class GO, class NO>
//void Xpetra::RegionManager<SC,LO,GO,NO>::setupNodeToRegionMapping()
//{
//  nodesToRegions_.clear();
//  interfaceNodesToRegions_.clear();
//  Teuchos::Array<std::tuple<GO, GO> > nodesReordered;
//  nodesReordered = nodeRegionPairs_;
//  std::sort(nodesReordered.begin(), nodesReordered.end(), Xpetra::compareNodes<GO>);
//
//  typename Teuchos::Array<std::tuple<GO,GO> >::iterator nodeIterator;
//  nodeIterator = nodesReordered.begin();
//
//  while (nodeIterator != nodesReordered.end())
//  {
//    GO currentNode = std::get<0>( *(nodeIterator) );
//    Teuchos::Array<GO> regions;
//    regions.clear();
//    regions.push_back(std::get<1>(*nodeIterator));
//
//    typename Array<std::tuple<GO, GO> >::iterator nextNodeIterator = nodeIterator + 1;
//
//    while (nextNodeIterator != nodesReordered.end())
//    {
//      GO nextNode = std::get<0>(*(nextNodeIterator));
//      if (currentNode == nextNode)
//      {
//        //As long as the information spanned regards the same node,
//        //the algorithm keeps on increasing the list of regions the given mesh node belong to
//        regions.push_back(std::get<1>(*(nextNodeIterator)));
//        ++nextNodeIterator;
//      }
//      else
//      {
//        //When the mesh node label changes, then the algorithm
//        //stops recording information about the previous node and it starts recording information for the new one
//        nodeIterator = nextNodeIterator;
//        break;
//      }
//    }
//    std::tuple<GO,Teuchos::Array<GO> > new_tuple;
//    std::sort(regions.begin(), regions.end());
//    new_tuple = std::make_tuple(currentNode, regions);
//    nodesToRegions_.push_back(new_tuple);
//
//    /* If a node belongs to multiple regions, it is located on an interface
//     * between regions. Hence, append it to the list of interface nodes.
//     */
//    if (regions.size() > 1)
//      interfaceNodesToRegions_.push_back(new_tuple);
//
//    if (nextNodeIterator == nodesReordered.end())
//      break;
//  }
//  TEUCHOS_TEST_FOR_EXCEPTION(!( nodesToRegions_.size() == numNodes_), Exceptions::RuntimeError,
//      "Number of nodes detected does not match with the initially declared one nodesToRegion tracks "
//      << nodesToRegions_.size() << " nodes whereas numNodes_ = " << numNodes_ << ".\n");
//
//  return;
//}
//
//template<class SC, class LO, class GO, class NO>
//void Xpetra::RegionManager<SC,LO,GO,NO>::setupProcToRegionMapping()
//{
//  int numProcs = comm_->getSize();
//  int myPID = comm_->getRank();
//
//  regionsPerProc_.clear();
//  regionsPerProc2_.clear();
//  procsPerRegion_.clear();
//
//  /* If the number of processes is smaller than the total number of regions,
//   * each process owns at least one entire regions.
//   *
//   * The number of regions per process is calculated to guarantee load balancing.
//   * After an initial distribution of regions, leftover regions that have not
//   * been assigned to any process yet are distributed in a round-robin fashion.
//   */
//  if (numProcs < numRegions_)
//  {
////    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not properly tested, yet.");
//
//    // compute minimum number of regions per processor
//    int minNumRegionsPerProc = std::floor(Teuchos::as<SC>(numRegions_) / Teuchos::as<SC>(numProcs));
//    std::cout << "Proc " << comm_->getRank() << " --- minNumRegionsPerProc = " << minNumRegionsPerProc << std::endl;
//
//    // compute number of 'leftover' regions
//    int numLeftoverRegions = numRegions_ % numProcs;
//    std::cout << "Proc " << comm_->getRank() << " --- numLeftoverRegions = " << numLeftoverRegions << std::endl;
//
//    for (int i = 0; i < minNumRegionsPerProc; ++i)
//    {
//      std::cout << "Proc " << comm_->getRank() << " --- i = " << i << std::endl;
//      regionsPerProc_.push_back(myPID * minNumRegionsPerProc + i + 1);
//    }
//
//    if (numLeftoverRegions >= myPID + 1 && numLeftoverRegions != 0)
//      regionsPerProc_.push_back(minNumRegionsPerProc * numProcs + (myPID + 1));
//
//    for (int i = 0; i < comm_->getSize(); ++i)
//    {
//      Teuchos::Array<GO> regions;
//      regions.clear();
//      for (int j = 0; j < minNumRegionsPerProc; ++j)
//      {
//        GO regionID =  i * minNumRegionsPerProc + j;
//        regions.push_back(regionID);
//      }
//      std::tuple<GO,Teuchos::Array<GO> > auxTuple = std::make_tuple(i, regions);
//      regionsPerProc2_.push_back(auxTuple);
//    }
//
//    for (int procID = 0; procID < comm_->getSize(); ++procID)
//    {
//      Teuchos::Array<GO> proc;
//      proc.clear();
//      proc.push_back(procID);
//      for (int i = 0; i < minNumRegionsPerProc; ++i)
//      {
//        GO region_index = (procID) * minNumRegionsPerProc + i;
//        std::tuple<GO,Teuchos::Array<GO> > tuple_aux = std::make_tuple(region_index, proc);
//        procsPerRegion_.push_back(tuple_aux);
//      }
//
//      if (numLeftoverRegions >= procID + 1 && numLeftoverRegions != 0)
//      {
//        GO region_index = minNumRegionsPerProc * numProcs + (procID + 1);
//        std::tuple<GO,Teuchos::Array<GO> > tuple_aux = std::make_tuple(region_index, proc);
//        procsPerRegion_.push_back(tuple_aux);
//      }
//    }
//    TEUCHOS_TEST_FOR_EXCEPTION(procsPerRegion_.size() != numRegions_, Exceptions::RuntimeError, "PID: "<<comm_->getRank()<<" - Number of regions detected does not match with the initially declared one \n procsPerRegion_ tracks "<<procsPerRegion_.size()<<" regions whereas numRegions_ = "<<numRegions_<<"\n");
//  }
//  //This is easy: if the number of regions coincides with the total number of processes instantiated, then
//  //a one-to-one relation between processes and regions is created
//  else if (numProcs == numRegions_)
//  {
//    regionsPerProc_.push_back(myPID + 1);
//
//    for (GO i = 0; i < numRegions_; ++i)
//    {
//      Teuchos::Array<GO> proc;
//      proc.clear();
//      proc.push_back(i);
//      std::tuple<GO,Teuchos::Array<GO> > tuple_aux = std::make_tuple(i, proc);
//      procsPerRegion_.push_back(tuple_aux);
//    }
//
//    for (int i = 0; i < comm_->getSize(); ++i)
//    {
//      Teuchos::Array<GO> regions;
//      regions.clear();
//      regions.push_back(i);
//      std::tuple<GO,Teuchos::Array<GO> > tuple_aux = std::make_tuple(i, regions);
//      regionsPerProc2_.push_back(std::make_tuple(i, regions));
//    }
//  }
//  //If the number of processes exceeds the number of regions in the domain,
//  //then each process is given a subset of a region.
//  //N.B.: A SINGLE PROCESS IS NOT ALLOWED TO OWN CHUNCKS OF MULTIPLE REGIONS.
//  //IN THIS CONFIGURATION EACH PROCESS IS CONFINED TO A SINGLE REGION
//  else if (numProcs > numRegions_)
//  {
//    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not properly tested, yet.");
//
//    int num_procs_region = std::ceil(
//        static_cast<double>(numProcs) / static_cast<double>(numRegions_));
//    int num_regions_extra_proc = numProcs % numRegions_;
//    int proc_count = 0;
//    std::tuple<int, Teuchos::Array<GO> > region_tuple;
//
//    for (int i = 0; i < numRegions_; ++i)
//    {
//      Teuchos::Array<GO> procs;
//      procs.clear();
//      if (i < num_regions_extra_proc || num_regions_extra_proc == 0)
//        for (int j = 0; j < num_procs_region; ++j)
//        {
//          procs.push_back(proc_count);
//          proc_count++;
//        }
//      else
//        for (int j = 0; j < num_procs_region - 1; ++j)
//        {
//          procs.push_back(proc_count);
//          proc_count++;
//        }
//      std::sort(procs.begin(), procs.end());
//      region_tuple = std::make_tuple(i, procs);
//      procsPerRegion_.push_back(region_tuple);
//    }
//    regionsPerProc_.clear();
//  }
//
//  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
//  printProcsPerRegion(*out);
//  printRegionsPerProc(*out);
//
//  return;
//}

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
  Teuchos::RCP<const Teuchos::Array<GO> > myNodesGIDsRcp = nodes_->getNodeGIDsPerProc(comm_->getRank());
  const Teuchos::Array<GO>& myNodesGIDs = *myNodesGIDsRcp;

  compositeMap_ = Xpetra::MapFactory<LO,GO,NO>::Build(Xpetra::UseTpetra, nodes_->getNumNodes(), myNodesGIDs(), 0, comm_);

//  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
//  compositeMap_->describe(*out, Teuchos::VERB_EXTREME);

  return;
}

template<class SC, class LO, class GO, class NO>
void Xpetra::RegionManager<SC,LO,GO,NO>::setupRegionRowMaps()
{


  return;
}

//
//template<class SC, class LO, class GO, class NO>
//void Xpetra::RegionManager<SC,LO,GO,NO>::setupRowMaps()
//{
//  TEUCHOS_TEST_FOR_EXCEPTION(
//      (procsPerRegion_.empty() && regionsPerProc_.empty()),
//      Exceptions::RuntimeError,
//      "Process ID: "<<comm_->getRank()<<" - Information about region partitioning across processors is not consistent: incorrect values for number of processors or number of regions \n");
//  Teuchos::Array<GO> elements;
//  Teuchos::Array<GO> region_elements;
//  Teuchos::Array<Teuchos::Array<GO> > elements_per_region;
//  Teuchos::Array<Teuchos::Array<std::tuple<GO, GO> > > regionToAll;
//  int myPID = comm_->getRank();
//
//  elements.clear();
//  region_elements.clear();
//  elements_per_region.clear();
//  elements_per_region.resize(numRegions_);
//  regionToAll.clear();
//  regionToAll.resize(numRegions_);
//
//  TEUCHOS_TEST_FOR_EXCEPTION(!nodesSortedByRegions_,
//      Exceptions::RuntimeError,
//      "Nodes are not sorted by regions in ascending order \n");
//  TEUCHOS_TEST_FOR_EXCEPTION(numNodes_ > nodeRegionPairs_.size(),
//      Exceptions::RuntimeError,
//      "Number of nodes declared in input file does not match with the effective number of nodes provided\n"<<"numNodes_ ="<<numNodes_<<" whereas nodeRegionPairs_ tracks "<<nodeRegionPairs_.size()<<" nodes \n");
//  TEUCHOS_TEST_FOR_EXCEPTION(
//      numRegions_ != std::get<1>(nodeRegionPairs_.back()) + 1,
//      Exceptions::RuntimeError,
//      "Number of regions declared in input file does not match with the effective number of regions provided.\n");
//
////  if (comm_->getSize() <= numRegions_)
////  {
////    // Loop over all processors
////    typename Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > >::const_iterator procIt;
////    for (procIt = regionsPerProc2_.begin(); procIt < regionsPerProc2_.end(); ++procIt)
////    {
////      const Teuchos::Array<GO>& myRegions = std::get<1>(*procIt);
////      Teuchos::Array<GO> regionElements;
////
////      // Loop over all regions on each processor
////      typename Teuchos::Array<GO>::const_iterator regionIt;
////      for (regionIt = myRegions.begin(); regionIt < myRegions.end(); ++regionIt)
////      {
////        regionElements.clear();
////        Xpetra::CheckerNode<GO> unaryPredicate(*regionIt);
////        typename Teuchos::Array<std::tuple<GO,GO> >::const_iterator nodeItLow;
////        typename Teuchos::Array<std::tuple<GO,GO> >::const_iterator nodeItHigh;
////
////        /* We position an iterator at the beginning of the information associated
////         * with the a region owned by the calling process and another iterator
////         * right at the end of the information
////         */
////        nodeItLow = std::find_if<typename Teuchos::Array<std::tuple<GO,GO> >::const_iterator,
////                Xpetra::CheckerNode<GO> >(nodeRegionPairs_.begin(), nodeRegionPairs_.end(), unaryPredicate);
////        nodeItHigh = std::find_if_not<typename Teuchos::Array<std::tuple<GO,GO> >::const_iterator,
////                Xpetra::CheckerNode<GO> >(nodeItLow, nodeRegionPairs_.end(), unaryPredicate);
////
////        // Count the number of mesh nodes inside a region
////        int numMyNodes = nodeItHigh - nodeItLow;
////        comm_->barrier();
////        std::cout << "Proc " << comm_->getRank() << ", procIt = " << procIt << "regionIt =  " << regionIt << " --- numMyNodes = " << numMyNodes << std::endl;
////
////        // Extract list of nodes of this region
////        GO region_node_label = 0;
////        typename Teuchos::Array<std::tuple<GO,GO> >::const_iterator nodes_iterator_aux;
////        for (nodes_iterator_aux = nodeItLow;
////            nodes_iterator_aux != nodeItHigh; ++nodes_iterator_aux)
////        {
////          GO nodeID = std::get<0>(*nodes_iterator_aux);
////          Xpetra::CheckerNodesToRegion<GO> unaryPredicateNode(nodeID);
////          typename Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > >::const_iterator nodes_to_region_iterator;
////          nodes_to_region_iterator = std::find_if<typename Teuchos::Array< std::tuple< GO,Teuchos::Array<GO> > >::const_iterator,
////              Xpetra::CheckerNodesToRegion<GO> >(nodesToRegions_.begin(), nodesToRegions_.end(), unaryPredicateNode);
////          Teuchos::Array<GO> nodal_regions = std::get<1>(*nodes_to_region_iterator);
////
////          //By default, I choose that a node is owned by the process associated with the region that shows up first in its list of beloning
////          //This guarantees that each row of the composite stiffness matrix is owned only by a single process, as Trilinos requires
////          if (*iter_array == nodal_regions[0])
////            elements.push_back(node);
////
////          //Nodes on the interface still belong to multiple regions, so
////          //it is important to keep track of this for the row maps of region matrices
////          regionElements.push_back(region_node_label);
////
////          //If a process owns a region (or even a portion of it), we provide to it a map
////          //from region indices to composite indices for all the nodes inside that region,
////          //even if a specific node is not owned by the calling process
////          regionToAll[*iter_array].push_back(std::make_tuple(region_node_label, std::get<0>(*nodes_iterator_aux)));
////          ++region_node_label;
////        }
////
////        elements_per_region[*regionIt] = regionElements;
////      }
////    }
////  }
//
//  TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Not implemented, yet.");
//
//  comm_->barrier();
//  if (comm_->getRank())
//  {
//    std::cout << "Region\tElements" << std::endl;
//    for (int i = 0; i < elements_per_region.size(); ++i)
//    {
//      std::cout << i << "\t" << elements_per_region[i] << std::endl;
//    }
//  }
//  comm_->barrier();
//
//  exit(0);
//
//  if (!(regionsPerProc_.empty()))
//  {
//    typename Teuchos::Array<GO>::iterator iter_array;
//    for (iter_array = regionsPerProc_.begin(); iter_array < regionsPerProc_.end(); ++iter_array)
//    {
//      region_elements.clear();
//      Xpetra::CheckerNode<GO> unaryPredicate(*iter_array);
//      typename Teuchos::Array<std::tuple<GO,GO> >::iterator nodes_iterator1;
//      typename Teuchos::Array<std::tuple<GO,GO> >::iterator nodes_iterator2;
//
//      /* We position an iterator at the beginning of the information associated
//       * with the a region owned by the calling process and another iterator
//       * right at the end of the information
//       */
//      nodes_iterator1 = std::find_if<typename Teuchos::Array<std::tuple<GO,GO> >::iterator,
//              Xpetra::CheckerNode<GO> >(nodeRegionPairs_.begin(), nodeRegionPairs_.end(), unaryPredicate);
//      nodes_iterator2 = std::find_if_not<typename Teuchos::Array<std::tuple<GO,GO> >::iterator,
//              Xpetra::CheckerNode<GO> >(nodes_iterator1, nodeRegionPairs_.end(), unaryPredicate);
//
//      // Count the number of mesh nodes inside a region
//      int num_region_nodes = nodes_iterator2 - nodes_iterator1;
//
//      //The policy assumes that in the input file the indexBase for the node label is 1 // (mayr.mt) changed this to 0
//      GO region_node_label = 0;
//      typename Teuchos::Array<std::tuple<GO,GO> >::iterator nodes_iterator_aux;
//      for (nodes_iterator_aux = nodes_iterator1;
//          nodes_iterator_aux != nodes_iterator2; ++nodes_iterator_aux)
//      {
//        GO node = std::get<0>(*nodes_iterator_aux);
//        Xpetra::CheckerNodesToRegion<GO> unaryPredicateNode(node);
//        typename Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > >::iterator nodes_to_region_iterator;
//        nodes_to_region_iterator = std::find_if<typename Teuchos::Array< std::tuple< GO,Teuchos::Array<GO> > >::iterator, Xpetra::CheckerNodesToRegion<GO> >(nodesToRegions_.begin(), nodesToRegions_.end(), unaryPredicateNode);
//        Teuchos::Array<GO> nodal_regions = std::get<1>(*nodes_to_region_iterator);
//
//        //By default, I choose that a node is owned by the process associated with the region that shows up first in its list of beloning
//        //This guarantees that each row of the composite stiffness matrix is owned only by a single process, as Trilinos requires
//        if (*iter_array == nodal_regions[0])
//          elements.push_back(node);
//
//        //Nodes on the interface still belong to multiple regions, so
//        //it is important to keep track of this for the row maps of region matrices
//        region_elements.push_back(region_node_label);
//
//        //If a process owns a region (or even a portion of it), we provide to it a map
//        //from region indices to composite indices for all the nodes inside that region,
//        //even if a specific node is not owned by the calling process
//        regionToAll[*iter_array].push_back(std::make_tuple(region_node_label, std::get<0>(*nodes_iterator_aux)));
//        ++region_node_label;
//      }
//
////      // C++ indexing starts from 0, so everything is shifted backward by one to make it consistent with programming language's policies
////      for (typename Teuchos::Array<GO>::iterator iter = region_elements.begin();
////          iter != region_elements.end(); ++iter)
////      {
////        *iter = *iter - 1;
////      }
//
//      elements_per_region[*iter_array-1] = region_elements;
//      TEUCHOS_TEST_FOR_EXCEPTION((num_region_nodes != regionToAll[*iter_array-1].size()), Exceptions::RuntimeError, "Process ID: "<<comm_->getRank()<<" - Number of region nodes does not match with number of nodes stored in regionToAll \n"<<"num_region_nodes= "<< num_region_nodes<<" whereas regionToAll["<<*iter_array-1<<"].size()= "<<regionToAll[*iter_array-1].size()<<"\n");
//    }
//
//    TEUCHOS_TEST_FOR_EXCEPTION((numRegions_ != regionToAll.size()), Exceptions::RuntimeError,
//        "regionToAll size has been corrupted\n"<<"numRegions_ = "<<numRegions_<<" whereas regionToAll.size() = "<<regionToAll.size()<<"\n");
//  }
//  else {
//    /* The code enters the scope of these curly brackets if the number of processes exceeds the number of regions in the domain.
//     * Therefore, each process owns only a subregion (or at most a single entire region).
//     * In this situation, the calling process must identify the region it is associated with
//     */
//    bool region_found = false;
//    GO myRegion = -1;
//    TEUCHOS_TEST_FOR_EXCEPTION(
//        !(procsPerRegion_.size() == numRegions_),
//        Exceptions::RuntimeError,
//        "Process ID: "<<comm_->getRank()<<" - Number of total regions does not match with regionHandler structures \n");
//
//    Teuchos::Array<GO> region_procs;
//    while (!region_found) {
//      typename Teuchos::Array<GO>::iterator iter_proc;
//      for (GO region_index = 0; region_index < procsPerRegion_.size(); ++region_index) {
//        region_procs = std::get<1>(procsPerRegion_[region_index]);
//        iter_proc = std::find(region_procs.begin(), region_procs.end(), myPID);
//        if (iter_proc != region_procs.end()) {
//          myRegion = region_index;
//          region_found = true;
//        }
//      }
//    }
//
//    TEUCHOS_TEST_FOR_EXCEPTION((myRegion == -1 || !region_found),
//        Exceptions::RuntimeError,
//        ("Region containing PROC ID: " + std::to_string(myPID) + " NOT FOUND \n"));
//    region_procs = std::get<1>(procsPerRegion_[myRegion]);
//
//    Xpetra::CheckerNode<GO> unaryPredicate(myRegion);
//    typename Teuchos::Array<std::tuple<GO,GO> >::iterator nodes_iterator1;
//    typename Teuchos::Array<std::tuple<GO,GO> >::iterator nodes_iterator2;
//    nodes_iterator1 = std::find_if<typename Teuchos::Array<std::tuple<GO,GO> >::iterator,
//        Xpetra::CheckerNode<GO> >(nodeRegionPairs_.begin(), nodeRegionPairs_.end(), unaryPredicate);
//    nodes_iterator2 = std::find_if_not<typename Teuchos::Array<std::tuple<GO,GO> >::iterator,
//        Xpetra::CheckerNode<GO> >(nodes_iterator1, nodeRegionPairs_.end(), unaryPredicate);
//
//    int num_region_nodes = nodes_iterator2 - nodes_iterator1;
//    int num_region_procs = region_procs.size();
//
//    if (num_region_nodes < num_region_procs)
//    {
//      Teuchos::Array<GO> region_procs_reduced;
//      region_procs_reduced.clear();
//      for (int i = 0; i < num_region_nodes; ++i)
//        region_procs_reduced.push_back(region_procs[i]);
//
//      typename Teuchos::Array<GO>::iterator proc_iterator;
//      proc_iterator = std::find<typename Teuchos::Array<GO>::iterator,GO>(
//          region_procs_reduced.begin(), region_procs_reduced.end(), myPID);
//
//      if (proc_iterator != region_procs_reduced.end()) //This reasoning works because the PROC ID for each region has been previously sorted in ascending order
//      {
//        GO node = std::get<0>(*(nodes_iterator1 + (proc_iterator - region_procs_reduced.begin() + 1)));
//        GO region_node_label = proc_iterator-region_procs_reduced.begin() + 1;
//        Xpetra::CheckerNodesToRegion<GO> unaryPredicateNode(node);
//        typename Teuchos::Array< std::tuple<GO, Teuchos::Array<GO> > >::iterator nodes_to_region_iterator;
//        nodes_to_region_iterator = std::find_if<typename Teuchos::Array< std::tuple< GO,Teuchos::Array<GO> > >::iterator,
//            Xpetra::CheckerNodesToRegion<GO> >(nodesToRegions_.begin(), nodesToRegions_.end(), unaryPredicateNode);
//        Teuchos::Array<GO> nodal_regions =  std::get<1>(*nodes_to_region_iterator);
//
//        // Guarantee uniqueness of RowMap
//        if (myRegion == nodal_regions[0])
//          elements.push_back(node);
//
//        /* Although a process does not own a row in the composite matrix,
//         * it still might own the row in the region matrix view.
//         */
//        region_elements.push_back(region_node_label);
//      }
//
//      // If a process owns a region (or even a portion of it), we provide to it a map
//      // from region indices to composite indices for all the nodes inside that region,
//      // even if a specific node is not owned by the calling process
//      // If a process owns something of a region, then the process has a global view of who owns what for that region
//      // Although this may seem more information than what actually needed, it is important for the computation of the collapsing.
//      // If the collapsing is not calculated, then this structure actually overestimates what a process needs to know.
//      for(proc_iterator = region_procs_reduced.begin(); proc_iterator!=region_procs_reduced.end(); ++proc_iterator)
//      {
//        GO node = std::get<0>(*(nodes_iterator1+(proc_iterator-region_procs_reduced.begin()+1)));
//        GO region_node_label = proc_iterator-region_procs_reduced.begin()+1;
//        Xpetra::CheckerNodesToRegion<GO> unaryPredicateNode(node);
//        typename Teuchos::Array< std::tuple<GO, Teuchos::Array<GO> > >::iterator nodes_to_region_iterator;
//        nodes_to_region_iterator = std::find_if<typename Teuchos::Array< std::tuple< GO,Teuchos::Array<GO> > >::iterator, Xpetra::CheckerNodesToRegion<GO> >(nodesToRegions_.begin(), nodesToRegions_.end(), unaryPredicateNode);
//        regionToAll[myRegion].push_back( std::make_tuple(region_node_label, node) );
//      }
//
//    }
//    else if (num_region_nodes == num_region_procs)
//    {
//      typename Teuchos::Array<GO>::iterator proc_iterator;
//      proc_iterator = std::find<typename Teuchos::Array<GO>::iterator, GO>(region_procs.begin(), region_procs.end(), myPID);
//
//      if( proc_iterator!=region_procs.end() )//This reasoning works because the PROC ID for each region has been previously sorted in ascending order
//      {
//        GO node = std::get<0>(*(nodes_iterator1+(proc_iterator-region_procs.begin()+1)));
//        GO region_node_label = proc_iterator-region_procs.begin()+1;
//        Xpetra::CheckerNodesToRegion<GO> unaryPredicateNode(node);
//        typename Teuchos::Array< std::tuple<GO, Teuchos::Array<GO> > >::iterator nodes_to_region_iterator;
//        nodes_to_region_iterator = std::find_if<typename Teuchos::Array< std::tuple< GO,Teuchos::Array<GO> > >::iterator,
//            Xpetra::CheckerNodesToRegion<GO> >(nodesToRegions_.begin(), nodesToRegions_.end(), unaryPredicateNode);
//        Teuchos::Array<GO> nodal_regions =  std::get<1>(*nodes_to_region_iterator);
//        if (myRegion == nodal_regions[0])
//          elements.push_back( node );
//
//        region_elements.push_back(region_node_label);
//      }
//
//      // If a process owns a region (or even a portion of it), we provide to it a map
//      // from region indices to composite indices for all the nodes inside that region,
//      // even if a specific node is not owned by the calling process
//      for (proc_iterator = region_procs.begin(); proc_iterator != region_procs.end(); ++proc_iterator)
//      {
//        GO node = std::get<0>(*(nodes_iterator1+(proc_iterator-region_procs.begin()+1)));
//        GO region_node_label = proc_iterator-region_procs.begin()+1;
//        Xpetra::CheckerNodesToRegion<GO> unaryPredicateNode(node);
//        typename Teuchos::Array< std::tuple<GO, Teuchos::Array<GO> > >::iterator nodes_to_region_iterator;
//        nodes_to_region_iterator = std::find_if<typename Teuchos::Array< std::tuple< GO,Teuchos::Array<GO> > >::iterator,
//            Xpetra::CheckerNodesToRegion<GO> >(nodesToRegions_.begin(), nodesToRegions_.end(), unaryPredicateNode);
//        regionToAll[myRegion].push_back(std::make_tuple(region_node_label, node));
//      }
//    }
//    else
//    {
//      typename Teuchos::Array<GO>::iterator proc_iterator;
//      proc_iterator = std::find<typename Teuchos::Array<GO>::iterator, GO>(region_procs.begin(), region_procs.end(), myPID);
//
//      int num_nodes_proc = std::ceil( static_cast<double>(num_region_nodes)/static_cast<double>(num_region_procs) );
//      int num_procs_extra_node = num_region_nodes % num_region_procs;
//
//      if (proc_iterator - region_procs.begin() < num_procs_extra_node + 1
//          or num_procs_extra_node == 0)
//      {
//        int init_node = num_nodes_proc * ( proc_iterator-region_procs.begin() );
//        for (int i = 0; i < num_nodes_proc; ++i)
//        {
//          GO node = std::get<0>( *( nodes_iterator1 + init_node + i ) );
//          GO region_node_label = init_node + i + 1;
//          Xpetra::CheckerNodesToRegion<GO> unaryPredicateNode(node);
//          typename Teuchos::Array< std::tuple<GO, Teuchos::Array<GO> > >::iterator nodes_to_region_iterator;
//          nodes_to_region_iterator = std::find_if<typename Teuchos::Array< std::tuple< GO,Teuchos::Array<GO> > >::iterator,
//              Xpetra::CheckerNodesToRegion<GO> >(nodesToRegions_.begin(), nodesToRegions_.end(), unaryPredicateNode);
//          Teuchos::Array<GO> nodal_regions =  std::get<1>(*nodes_to_region_iterator);
//          if (myRegion == nodal_regions[0])
//            elements.push_back(node);
//
//          region_elements.push_back(region_node_label);
//        }
//      }
//      else
//      {
//        int init_node = num_nodes_proc * num_procs_extra_node + (proc_iterator - region_procs.begin() - num_procs_extra_node) * (num_nodes_proc-1);
//        for (int i = 0; i < num_nodes_proc - 1; ++i)
//        {
//          GO node = std::get<0>( *( nodes_iterator1 + init_node + i ) );
//          GO region_node_label = init_node + i + 1;
//          Xpetra::CheckerNodesToRegion<GO> unaryPredicateNode(node);
//          typename Teuchos::Array< std::tuple<GO, Teuchos::Array<GO> > >::iterator nodes_to_region_iterator;
//          nodes_to_region_iterator = std::find_if<typename Teuchos::Array< std::tuple< GO,Teuchos::Array<GO> > >::iterator,
//              Xpetra::CheckerNodesToRegion<GO> >(nodesToRegions_.begin(), nodesToRegions_.end(), unaryPredicateNode);
//          Teuchos::Array<GO> nodal_regions = std::get<1>(*nodes_to_region_iterator);
//          if (myRegion == nodal_regions[0])
//            elements.push_back(node);
//
//          region_elements.push_back(region_node_label);
//        }
//      }
//
//      // If a process owns a region (or even a portion of it), we provide to it a map
//      // from region indices to composite indices for all the nodes inside that region,
//      // even if a specific node is not owned by the calling process
//      for (proc_iterator = region_procs.begin(); proc_iterator != region_procs.end(); ++proc_iterator)
//      {
//        if (proc_iterator - region_procs.begin() < num_procs_extra_node
//            or num_procs_extra_node == 0)
//        {
//          int init_node = num_nodes_proc * ( proc_iterator-region_procs.begin() );
//          for( int i=0; i<num_nodes_proc; ++i )
//          {
//            GO node = std::get<0>( *( nodes_iterator1 + init_node + i ) );
//            GO region_node_label = init_node + i + 1;
//            Xpetra::CheckerNodesToRegion<GO> unaryPredicateNode(node);
//            typename Teuchos::Array< std::tuple<GO, Teuchos::Array<GO> > >::iterator nodes_to_region_iterator;
//            nodes_to_region_iterator = std::find_if<typename Teuchos::Array< std::tuple< GO,Teuchos::Array<GO> > >::iterator, Xpetra::CheckerNodesToRegion<GO> >(nodesToRegions_.begin(), nodesToRegions_.end(), unaryPredicateNode);
//            regionToAll[myRegion].push_back( std::make_tuple(region_node_label, node) );
//          }
//        }
//        else
//        {
//          int init_node = num_nodes_proc * num_procs_extra_node + (proc_iterator - region_procs.begin() - num_procs_extra_node) * (num_nodes_proc-1);
//          for (int i = 0; i < num_nodes_proc - 1; ++i)
//          {
//            GO node = std::get<0>(*(nodes_iterator1 + init_node + i));
//            GO region_node_label = init_node + i + 1;
//            Xpetra::CheckerNodesToRegion<GO> unaryPredicateNode(node);
//            typename Teuchos::Array< std::tuple<GO, Teuchos::Array<GO> > >::iterator nodes_to_region_iterator;
//            nodes_to_region_iterator = std::find_if<typename Teuchos::Array< std::tuple< GO,Teuchos::Array<GO> > >::iterator,
//                Xpetra::CheckerNodesToRegion<GO> >(nodesToRegions_.begin(), nodesToRegions_.end(), unaryPredicateNode);
//            regionToAll[myRegion].push_back(std::make_tuple(region_node_label, node));
//          }
//        }
//      }
//
//    }
//    TEUCHOS_TEST_FOR_EXCEPTION( ( numRegions_!=regionToAll.size() ), Exceptions::RuntimeError, "Process ID: "<<comm_->getRank()<<" - regionToAll size has been corrupted\n"<<"numRegions_ = "<<numRegions_<<" whereas regionToAll.size()= "<<regionToAll.size()<<"\n");
//
////    //C++ indexing starts from 0, so everything is shifted backward by one to make it consistent with programming language's policies
////    for( typename Teuchos::Array<GO>::iterator iter = region_elements.begin(); iter!=region_elements.end(); ++iter )
////      *iter = *iter - 1;
//
//    elements_per_region[myRegion] = region_elements;
//  }
//
////  //C++ indexing starts from 0, so everything is shifted backward by one to make it consistent with programming language's policies
////  for( typename Teuchos::Array<GO>::iterator iter = elements.begin(); iter!=elements.end(); ++iter )
////    *iter = *iter - 1;
//
//  // Finally, create and fill the maps_ object
//  maps_ = Teuchos::rcp(new Xpetra::SplittingMapsInfo<GO>());
//  maps_->compositeMap_ = elements;
//  maps_->regionMaps_ = elements_per_region;
//  maps_->regionToAll_ = regionToAll;
//
//  return;
//}
//
//template<class SC, class LO, class GO, class NO>
//void Xpetra::RegionManager<SC,LO,GO,NO>::countNodesPerRegion()
//{
//  numNodesPerRegion_.clear();
//  numNodesPerRegion_.resize(numRegions_);
//
//  // Loop over all regions and count nodes
//  for (int i = 0; i < Teuchos::as<int>(numRegions_); ++i)
//  {
//    Xpetra::CheckerNode<GO> unaryPredicate(i);
//
//    typename Teuchos::Array<std::tuple<GO,GO> >::iterator nodes_iterator1;
//    typename Teuchos::Array<std::tuple<GO,GO> >::iterator nodes_iterator2;
//
//    nodes_iterator1 = std::find_if<typename Teuchos::Array<std::tuple<GO,GO> >::iterator,
//        Xpetra::CheckerNode<GO> >(nodeRegionPairs_.begin(), nodeRegionPairs_.end(), unaryPredicate);
//    nodes_iterator2 = std::find_if_not<typename Teuchos::Array<std::tuple<GO,GO> >::iterator,
//        Xpetra::CheckerNode<GO> >(nodes_iterator1, nodeRegionPairs_.end(), unaryPredicate);
//
//    numNodesPerRegion_[i] = nodes_iterator2 - nodes_iterator1;
//
//    std::cout << "Proc " << comm_->getRank() << " --- #nodes in region " << i << "/" << numRegions_<< ": " << numNodesPerRegion_[i] << std::endl;
//  }
//
//  return;
//}
//
//template<class SC, class LO, class GO, class NO>
//void Xpetra::RegionManager<SC,LO,GO,NO>::printNodeRegionPairs(
//    Teuchos::FancyOStream& out) const
//{
//  comm_->barrier();
//  if (comm_->getRank() == 0)
//  {
//    out << std::endl << "*** NodeRegionPairs:" << std::endl
//        <<              "    ----------------" << std::endl;
//
//    out << "Total number of mesh nodes: " << numNodes_ << std::endl
//        << "Total number of mesh regions: " << numRegions_ << std::endl
//        << "Number of rows in nodeRegionPairs_ structure: " << nodeRegionPairs_.size()
//        << std::endl;
//
//    for (GO i = 0; i < nodeRegionPairs_.size(); ++i)
//    {
//      out << std::get<0>(nodeRegionPairs_[i])
//          << "\t" << std::get<1>(nodeRegionPairs_[i])
//          << std::endl;
//    }
//  }
//
//  return;
//}
//
//template<class SC, class LO, class GO, class NO>
//void Xpetra::RegionManager<SC,LO,GO,NO>::printNodesToRegionMapping(
//    Teuchos::FancyOStream& out) const
//{
//  comm_->barrier();
//  if (comm_->getRank() == 0)
//  {
//    out << std::endl << "*** NodesToRegionMapping:" << std::endl
//        <<              "    ---------------------" << std::endl;
//
//    out << "Total number of mesh nodes: " << numNodes_ << std::endl
//        << "Total number of mesh regions: " << numRegions_ << std::endl
//        << "Number of rows in nodeRegionPairs_ structure: " << nodeRegionPairs_.size()
//        << std::endl;
//    for (GO i = 0; i < nodesToRegions_.size(); ++i)
//    {
//      out << "Node " << std::get<0>(nodesToRegions_[i])
//          << "\t belongs to regions: " << std::get<1>(nodesToRegions_[i])
//          << std::endl;
//    }
//  }
//
//  return;
//}
//
//template<class SC, class LO, class GO, class NO>
//void Xpetra::RegionManager<SC,LO,GO,NO>::printInterfaceNodesToRegionMapping(
//    Teuchos::FancyOStream& out) const
//{
//  comm_->barrier();
//  if (comm_->getRank() == 0)
//  {
//    out << std::endl << "*** InterfaceNodesToRegionMapping:" << std::endl
//        <<              "    ------------------------------" << std::endl;
//
//    out << "Total number of interface nodes: " << interfaceNodesToRegions_.size() << std::endl
//        << "Total number of mesh regions: " << numRegions_ << std::endl
//        << "Number of rows in interfaceNodesToRegions_ structure: " << interfaceNodesToRegions_.size()
//        << std::endl;
//    for (GO i = 0; i < interfaceNodesToRegions_.size(); ++i)
//    {
//      out << "Node " << std::get<0>(interfaceNodesToRegions_[i])
//          << "\t belongs to regions: " << std::get<1>(interfaceNodesToRegions_[i])
//          << std::endl;
//    }
//  }
//  comm_->barrier();
//
//  return;
//}
//
//template<class SC, class LO, class GO, class NO>
//void Xpetra::RegionManager<SC,LO,GO,NO>::printInactiveProcs(
//    Teuchos::FancyOStream& out) const
//{
//  comm_->barrier();
//  if (comm_->getRank() == 0)
//    out << std::endl << "*** InactiveProcs:" << std::endl
//        <<              "    --------------" << std::endl;
//
//  if (maps_->compositeMap_.empty())
//    out << "INACTIVE PROC ID: " << comm_->getRank() << std::endl;
//
//  comm_->barrier();
//
//  return;
//}
//
//template<class SC, class LO, class GO, class NO>
//void Xpetra::RegionManager<SC,LO,GO,NO>::printProcsPerRegion(
//    Teuchos::FancyOStream& out) const
//{
//  comm_->barrier();
//  if (comm_->getRank() == 0)
//  {
//    out << std::endl << "*** ProcsPerRegion:" << std::endl
//        <<              "    ---------------" << std::endl;
//
//    out << "Total number of procs: " << comm_->getSize() << std::endl
//        << "Total number of mesh regions: " << numRegions_ << std::endl
//        << "Number of rows in procsPerRegion_ structure: " << procsPerRegion_.size()
//        << std::endl;
//
//    out << std::endl << "Region\tProcs" << std::endl;
//
//    for (GO i = 0; i < procsPerRegion_.size(); ++i)
//    {
//      out << std::get<0>(procsPerRegion_[i])
//          << "\t" << std::get<1>(procsPerRegion_[i])
//          << std::endl;
//    }
//  }
//  comm_->barrier();
//
//  return;
//}
//
//template<class SC, class LO, class GO, class NO>
//void Xpetra::RegionManager<SC,LO,GO,NO>::printRegionsPerProc(
//    Teuchos::FancyOStream& out) const
//{
//  comm_->barrier();
//  if (comm_->getRank() == 0)
//  {
//    out << std::endl << "*** RegionsPerProc:" << std::endl
//        <<              "    ---------------" << std::endl;
//
//    out << "Total number of procs: " << comm_->getSize() << std::endl
//        << "Total number of mesh regions: " << numRegions_ << std::endl
//        << "Number of rows in regionsPerProc_ structure: " << regionsPerProc_.size()
//        << std::endl;
//
//    out << std::endl << "Proc\tRegions" << std::endl;
//
//    for (GO i = 0; i < regionsPerProc_.size(); ++i)
//    {
//      out << i
//          << "\t" << regionsPerProc_[i]
//          << std::endl;
//    }
//  }
//  comm_->barrier();
//
//  comm_->barrier();
//  if (comm_->getRank() == 0)
//  {
//    out << std::endl << "*** RegionsPerProc:" << std::endl
//        <<              "    ---------------" << std::endl;
//
//    out << "Total number of procs: " << comm_->getSize() << std::endl
//        << "Total number of mesh regions: " << numRegions_ << std::endl
//        << "Number of rows in regionsPerProc2_ structure: " << regionsPerProc2_.size()
//        << std::endl;
//
//    out << std::endl << "Proc\tRegions" << std::endl;
//
//    for (GO i = 0; i < regionsPerProc2_.size(); ++i)
//    {
//      out << std::get<0>(regionsPerProc2_[i])
//          << "\t" << std::get<1>(regionsPerProc2_[i])
//          << std::endl;
//    }
//  }
//  comm_->barrier();
//
//  return;
//}
//
//template<class SC, class LO, class GO, class NO>
//void Xpetra::RegionManager<SC,LO,GO,NO>::printNumRegionsPerProc(
//    Teuchos::FancyOStream& out) const
//{
//  comm_->barrier();
//  if (comm_->getRank() == 0)
//    out << std::endl << "*** NumRegionsPerProc:" << std::endl
//        <<              "    ------------------" << std::endl;
//
//  if (regionsPerProc_.empty())
//  {
//    // One-to-one mapping between procs and regions
//    out << "We have a one-to-one mapping between regions and MPI processes."
//        << std::endl;
//  }
//  else
//  {
//    if (comm_->getRank() == 0)
//    {
//      out << "Total number of procs: " << comm_->getSize() << std::endl
//          << "Total number of mesh regions: " << numRegions_ << std::endl
//          << "Number of rows in regionsPerProc_ structure: " << regionsPerProc_.size()
//          << std::endl;
//
//      out << "Number of regions per processor:" << std::endl;
//
//      for (GO i = 0; i < regionsPerProc_.size(); ++i)
//      {
//        out << "Proc " << i << ":\t" << regionsPerProc_[i] << " regions"
//            << std::endl;
//      }
//    }
//  }
//  comm_->barrier();
//
//  return;
//}
//
//template<class SC, class LO, class GO, class NO>
//Teuchos::Array<GO> Xpetra::RegionManager<SC,LO,GO,NO>::getRegionRowMap(
//    GO regionID) const
//{
//  TEUCHOS_TEST_FOR_EXCEPTION(regionID >= numRegions_, Exceptions::RuntimeError,
//      "Value of region index exceeds total number of stored regions.\n\n"<<"Trying to access information about region " << regionID << " while the total number of regions stored is " << numRegions_ << "\n");
//
//  return maps_->regionMaps_[regionID];
//}
//
//template<class SC, class LO, class GO, class NO>
//Teuchos::Array<std::tuple<GO,GO> > Xpetra::RegionManager<SC,LO,GO,NO>::getRegionToAll(
//    GO regionID) const
//{
//  TEUCHOS_TEST_FOR_EXCEPTION(regionID >= numRegions_, Exceptions::RuntimeError,
//      "Value of region index exceeds total number of regions stored \n" << "Trying to access informaiton about region " << regionID << " when the total number of regions stored is " <<numRegions_ << "\n");
//
//  return maps_->regionToAll_[regionID];
//}

#endif /* XPETRA_REGION_MANAGER_IMPL_HPP_ */
