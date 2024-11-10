// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This is a change for Eric.

#include "Panzer_IOSSConnManager.hpp"

#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <limits>

// MPI includes
#include <mpi.h>

// Teuchos includes
#include "Teuchos_Comm.hpp"
#include "Teuchos_EReductionType.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TestForException.hpp"

// Kokkos includes
#include "Kokkos_DynRankView.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// Phalanx includes
#include "Phalanx_KokkosDeviceTypes.hpp"

// Panzer includes
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
//#include "Panzer_STK_PeriodicBC_Matcher.hpp"
//#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_IOSSDatabaseTypeManager.hpp"

// Zoltan includes
#include "Zoltan2_findUniqueGids.hpp"

// Ioss includes
#include "Ioss_CodeTypes.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_Region.h"
#include "Ioss_EntityBlock.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_ElementBlock.h"

namespace panzer_ioss {

using Teuchos::RCP;
using Teuchos::rcp;

IOSSConnManager::IOSSConnManager(Ioss::DatabaseIO * iossMeshDB)
   : iossMeshDB_(iossMeshDB), iossRegion_(), iossNodeBlocks_(), iossElementBlocks_(),
	 iossConnectivity_(), iossToShardsTopology_(), iossElementBlockTopologies_(),
	 elementBlocks_(), neighborElementBlocks_(),
	 numUniqueEdges_(0), numUniqueFaces_(0), edgeNodeToEdgeMap_(), faceNodeToFaceMap_(),
	 elementEdgeNodes_(), elementFaceNodes_(),
	 elmtLidToGid_(), elmtLidToConn_(), connSize_(),
	 connectivity_(), ownedElementCount_(0), sidesetsToAssociate_(),
	 sidesetYieldedAssociations_(), elmtToAssociatedElmts_(), placeholder_()
{

  std::string error_message;
  int bad_count;

  // Check for a non-null database
  TEUCHOS_TEST_FOR_EXCEPTION(iossMeshDB_ == nullptr, std::logic_error, "Error, iossMeshDB is nullptr.");
  TEUCHOS_TEST_FOR_EXCEPTION(iossMeshDB_ == NULL, std::logic_error, "Error, iossMeshDB is NULL.");

  // Check that the database is OK and that the file can be opened.
  bool status_ok = iossMeshDB_->ok(true, &error_message, &bad_count);
  TEUCHOS_TEST_FOR_EXCEPTION(!status_ok, std::logic_error,
		  "Error, Ioss::DatabaseIO::ok() failed.\n error_message = "
		  << error_message << "\nbad_count = " << bad_count << "\n");

  // Get the metadata
  iossRegion_ = rcp(new Ioss::Region(iossMeshDB_, "iossRegion_"));
  iossNodeBlocks_ = iossRegion_->get_node_blocks();
  iossElementBlocks_ = iossRegion_->get_element_blocks();

  // Get element block topologies
  for (Ioss::ElementBlock * iossElementBlock : iossElementBlocks_) {
    iossElementBlockTopologies_[iossElementBlock->name()] = iossElementBlock->topology();
  }

  // Create iossToShardsTopology_ map;
  createTopologyMapping();

}

Teuchos::RCP<panzer::ConnManager> IOSSConnManager::noConnectivityClone(
  std::string & type, Ioss::PropertyManager & properties) const {
	   std::string filename = iossMeshDB_->get_filename();

	   // Verify that multiple open databases are supported
	   TEUCHOS_TEST_FOR_EXCEPTION(IossDatabaseTypeManager::supportsMultipleOpenDatabases(type), std::logic_error, "Error, " << type << "does not support multiple open databases.");

	   // temporarily try this with a _copy file
	   size_t dotPosition = filename.find(".");
	   std::string filenameBase(filename, 0, dotPosition);
	   std::string filenameExtension(filename, dotPosition+1);
	   std::string filenameCopy = filenameBase + "_copy." + filenameExtension;

	   // Make a copy of the property manager
	   Ioss::PropertyManager propertiesCopy;
	   Ioss::NameList propertyNames;
	   int numProperties = properties.describe(&propertyNames);
	   std::string propertyName;
	   for (int p = 0; p < numProperties; ++p) {
	     Ioss::Property propertyCopy(properties.get(propertyNames[p]));
	     propertiesCopy.add(propertyCopy);
	   }

	   Ioss::DatabaseUsage db_usage = iossMeshDB_->usage();
	   Ioss::ParallelUtils util = iossMeshDB_->util();
	   MPI_Comm communicator = util.communicator();
	   Ioss::DatabaseIO * iossMeshDB = Ioss::IOFactory::create(type, filenameCopy,
	   	       		db_usage, communicator, propertiesCopy);

	   return Teuchos::rcp(new IOSSConnManager(iossMeshDB));
}

std::string IOSSConnManager::getBlockId(LocalOrdinal localElmtId) const {
   std::string blockId = "";
   Teuchos::RCP<std::vector<LocalOrdinal>> elementBlock;
   bool found = false;
   for (auto it = elementBlocks_.begin(); it != elementBlocks_.end(); ++it) {
	 elementBlock = it->second;
	 if ((*elementBlock).size() == 0)
	   continue;
	 if (*(elementBlock->end()-1) < localElmtId)
       continue;
	 for (auto vec_it = elementBlock->begin(); vec_it != elementBlock->end(); ++vec_it) {
       if (*vec_it == localElmtId) {
         blockId = it->first;
         found = true;
         break;
       }
     if (found)
       break;
	 }
   }
   TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error,
     "Could not find blockId for local element" << localElmtId
	 << ", global element " << elmtLidToGid_[localElmtId] << ".");
   return blockId;
}

void IOSSConnManager::getDofCoords(const std::string & blockId,
                                       const panzer::Intrepid2FieldPattern & coordProvider,
                                       std::vector<LocalOrdinal> & localCellIds,
                                       Kokkos::DynRankView<double,PHX::Device> & points) const {

  // Check that the FieldPattern has a topology that is compatible with all elements;
  TEUCHOS_TEST_FOR_EXCEPTION(!compatibleTopology(coordProvider), std::logic_error,
		     		            "Error, coordProvider must have a topology compatible with all elements.");

  int dim = coordProvider.getDimension();
  int idsPerElement = coordProvider.numberIds();

  const CellTopologyData * baseTopologyData = coordProvider.getCellTopology().getBaseCellTopologyData();
  int cornerNodesPerElement = int(baseTopologyData->subcell_count[0]);

  // Get the local Cell Ids
  localCellIds = getElementBlock(blockId);
  int numElements = localCellIds.size();

  // Get the DOF Coordinates
  std::vector<int> iossNodeBlockNodeIds;
  std::vector<double> iossNodeBlockCoordinates;
  std::vector<double> iossNodeCoordinates(dim);
  std::map<GlobalOrdinal, std::vector<double>> iossNodeIdToCoordinateMap;
    for (Ioss::NodeBlock * NodeBlock : iossNodeBlocks_) {
      NodeBlock->get_field_data("ids", iossNodeBlockNodeIds);
      NodeBlock->get_field_data("mesh_model_coordinates", iossNodeBlockCoordinates);
      for (size_t node = 0; node < iossNodeBlockNodeIds.size(); ++node) {
        for (int k = 0; k < dim; ++k) {
          iossNodeCoordinates[k] = iossNodeBlockCoordinates[node*dim+k];
        }
        iossNodeIdToCoordinateMap[GlobalOrdinal(iossNodeBlockNodeIds[node])] = iossNodeCoordinates;
      }
    }

  Teuchos::RCP<std::vector<GlobalOrdinal>> blockConnectivity;
  Kokkos::DynRankView<double,PHX::Device> vertices("vertices",numElements,cornerNodesPerElement,dim);

  blockConnectivity = iossConnectivity_.find(blockId)->second;
  const Ioss::ElementTopology * iossElementBlockTopology = iossElementBlockTopologies_.find(blockId)->second;
  int iossNodesPerElement = iossElementBlockTopology->number_nodes();
  for (int elmtIdxInBlock = 0; elmtIdxInBlock < numElements; ++elmtIdxInBlock) {
    for (int nodeIdxInElmt = 0; nodeIdxInElmt < cornerNodesPerElement; ++nodeIdxInElmt) {
      iossNodeCoordinates = iossNodeIdToCoordinateMap[(*(blockConnectivity))[iossNodesPerElement*elmtIdxInBlock + iossElementBlockTopology->element_connectivity()[nodeIdxInElmt]]];
      for (LocalOrdinal k = 0; k < dim; ++k) {
        vertices(elmtIdxInBlock,nodeIdxInElmt,k) = iossNodeCoordinates[k];
      }
    }
  }

  // setup output array
  points = Kokkos::DynRankView<double,PHX::Device>("points",numElements,idsPerElement,dim);
  coordProvider.getInterpolatoryCoordinates(vertices,points);
}

void IOSSConnManager::clearLocalElementMapping()
{

   elementBlocks_.clear();
   elmtLidToGid_.clear();
   elmtLidToConn_.clear();
   connSize_.clear();
   //elmtToAssociatedElmts_.clear();
}

void IOSSConnManager::buildLocalElementMapping() {

  std::string name;
  std::vector<int> indices;
  std::size_t element_count = 0;

  // Initialize
  clearLocalElementMapping();
  ownedElementCount_ = 0;
  elmtLidToGid_.clear();

  for (Ioss::ElementBlock * iossElementBlock : iossElementBlocks_) {
    name = iossElementBlock->name();
	iossElementBlock->get_field_data("ids", indices);

	// Ioss::GroupingEntity::get_field_data() will resize indices
	// up or down as needed.
	ownedElementCount_ += indices.size();

    elementBlocks_[name] = rcp(new std::vector<LocalOrdinal>);
    for (GlobalOrdinal element : indices) {
      elementBlocks_[name]->push_back(element_count);
      elmtLidToGid_.push_back(element);
      ++element_count;
    }
  }

  TEUCHOS_TEST_FOR_EXCEPTION(ownedElementCount_ != element_count, std::logic_error,
		                     "Error, ownedElementCount_ != element_count.");

  // Allocate space for elmtLidToConn_ and initialize to zero.
  elmtLidToConn_.clear();
  elmtLidToConn_.resize(ownedElementCount_,0);

  // Allocate space for connSize_ and initialize to zero.
  connSize_.clear();
  connSize_.resize(ownedElementCount_,0);
}

void IOSSConnManager::buildEdgeFaceCornerNodeMapping() {

  // get the local maximum node, edge, and face numbers
  std::vector<std::string> elementBlockIds;
  getElementBlockIds(elementBlockIds);

  size_t dim;
  int nodesPerElement = 0;
  int edgesPerElement = 0;
  int facesPerElement = 0;
  int numCornerNodesThisSubcell;
  std::vector<int> numCornerNodesOnEdges, numCornerNodesOnFaces;
  std::vector<std::vector<int>> edgeConnectivities, faceConnectivities;
  std::vector<GlobalOrdinal> nodesThisElement;
  std::array<GlobalOrdinal,MAX_SUBCELL_CORNER_NODES_> cornerNodesThisSubcell;
  std::vector<std::array<GlobalOrdinal,MAX_SUBCELL_CORNER_NODES_>> edgesThisElement, facesThisElement;

  Teuchos::RCP<std::vector<GlobalOrdinal>> blockConnectivity;
  for (std::string elementBlockId : elementBlockIds) {
    blockConnectivity = iossConnectivity_.find(elementBlockId)->second;
	int blockConnectivitySize = blockConnectivity->size();
	if (blockConnectivitySize > 0) {
	  const Ioss::ElementTopology * top = iossElementBlockTopologies_.find(elementBlockId)->second;
	  dim = top->spatial_dimension();
	  nodesPerElement = top->number_nodes();
	  if (dim > 1) {
	    edgesPerElement = top->number_edges();
	    numCornerNodesOnEdges.resize(edgesPerElement);
	    edgeConnectivities.resize(edgesPerElement);
	    for (int edge = 0; edge < edgesPerElement; ++edge) {
	      numCornerNodesOnEdges[edge] = top->edge_type(edge)->number_corner_nodes();
	      edgeConnectivities[edge] = top->edge_connectivity(edge+1);
	    }
	    if (dim > 2) {
	      facesPerElement = top->number_faces();
	      numCornerNodesOnFaces.resize(facesPerElement);
	      faceConnectivities.resize(facesPerElement);
	      for (int face = 0; face < facesPerElement; ++face) {
	        numCornerNodesOnFaces[face] = top->face_type(face)->number_corner_nodes();
	        faceConnectivities[face] = top->face_connectivity(face+1);
	      }
	    }
	  }
	  elementEdgeNodes_[elementBlockId] = rcp(new std::vector<std::vector<std::array<GlobalOrdinal,MAX_SUBCELL_CORNER_NODES_>>>);
	  elementFaceNodes_[elementBlockId] = rcp(new std::vector<std::vector<std::array<GlobalOrdinal,MAX_SUBCELL_CORNER_NODES_>>>);
	  nodesThisElement.resize(nodesPerElement);
	  for(std::size_t elmtIdInBlock=0; elmtIdInBlock < getElementBlock(elementBlockId).size(); ++elmtIdInBlock) {
	    edgesThisElement.clear();
	    facesThisElement.clear();
		for (int node = 0; node < nodesPerElement; ++node) {
	      nodesThisElement[node] = (*blockConnectivity)[nodesPerElement*elmtIdInBlock + node];
	    }
	    if (dim > 1) {
	      for (int edge = 0; edge < edgesPerElement; ++edge) {
	        numCornerNodesThisSubcell = numCornerNodesOnEdges[edge];
	        TEUCHOS_TEST_FOR_EXCEPTION(numCornerNodesThisSubcell > MAX_SUBCELL_CORNER_NODES_, std::logic_error,
	                    "Error, currently IOSSConnManager only supports element edges with at most " << MAX_SUBCELL_CORNER_NODES_ << " corner nodes.");
	        //cornerNodesThisSubcell.resize(numCornerNodesThisSubcell);
	        for (int cornerNode = 0; cornerNode < numCornerNodesThisSubcell; ++cornerNode) {
	          cornerNodesThisSubcell[cornerNode] = nodesThisElement[ edgeConnectivities[edge][cornerNode] ];
	        }
	        for (int cornerNode = numCornerNodesThisSubcell; cornerNode < MAX_SUBCELL_CORNER_NODES_; ++cornerNode) {
	          cornerNodesThisSubcell[cornerNode] = -1;
	        }
	        // Sort is necessary to ensure edge identifiers are unique.
	        std::sort(cornerNodesThisSubcell.begin(), cornerNodesThisSubcell.begin()+numCornerNodesThisSubcell);
	        edgeNodeToEdgeMap_[cornerNodesThisSubcell] = -1;
	        edgesThisElement.push_back(cornerNodesThisSubcell);
	      }
	      if (dim > 2) {
	        for (int face = 0; face < facesPerElement; ++face) {
	    	  numCornerNodesThisSubcell = numCornerNodesOnFaces[face];
	    	  TEUCHOS_TEST_FOR_EXCEPTION(numCornerNodesThisSubcell > MAX_SUBCELL_CORNER_NODES_, std::logic_error,
	    	  	                    "Error, currently IOSSConnManager only supports element faces with at most " << MAX_SUBCELL_CORNER_NODES_ << " corner nodes.");
	    	  //cornerNodesThisSubcell.resize(numCornerNodesThisSubcell);
	    	  for (int cornerNode = 0; cornerNode < numCornerNodesThisSubcell; ++cornerNode) {
	    	    cornerNodesThisSubcell[cornerNode] = nodesThisElement[ faceConnectivities[face][cornerNode] ];
	    	  }
	    	  for (int cornerNode = numCornerNodesThisSubcell; cornerNode < MAX_SUBCELL_CORNER_NODES_; ++cornerNode) {
	    	    cornerNodesThisSubcell[cornerNode] = -1;
	    	  }
	    	  // Sort is necessary to ensure edge identifiers are unique.
	    	  std::sort(cornerNodesThisSubcell.begin(), cornerNodesThisSubcell.begin()+numCornerNodesThisSubcell);
	    	  faceNodeToFaceMap_[cornerNodesThisSubcell] = -1;
	    	  facesThisElement.push_back(cornerNodesThisSubcell);
	    	}
	      }
	    }
	    elementEdgeNodes_[elementBlockId]->push_back(edgesThisElement);
	    elementFaceNodes_[elementBlockId]->push_back(facesThisElement);
	  }
	}
  }

  dim = iossElementBlockTopologies_.find(elementBlockIds[0])->second->spatial_dimension();

  std::vector<std::array<GlobalOrdinal,MAX_SUBCELL_CORNER_NODES_>> keys;
  std::vector<GlobalOrdinal> gids;
  GlobalOrdinal numEdges = 0;
  GlobalOrdinal numFaces = 0;

  if (dim > 1) {
    numEdges = edgeNodeToEdgeMap_.size();
    if (dim > 2) {
      numFaces = faceNodeToFaceMap_.size();
    }
  }

#ifdef HAVE_MPI
  Ioss::ParallelUtils util = iossMeshDB_->util();
  MPI_Comm communicator = util.communicator();
  Teuchos::MpiComm<int> comm(communicator);

  if (dim > 1) {
	keys.clear();
	for (auto pair : edgeNodeToEdgeMap_) {
	  keys.push_back(pair.first);
	}
    gids.resize(numEdges);
    TEUCHOS_TEST_FOR_EXCEPTION(keys.size() != gids.size(), std::logic_error,
       "Error, keys() != gids.size().");
    numUniqueEdges_ = Zoltan2::findUniqueGids<std::array<GlobalOrdinal,MAX_SUBCELL_CORNER_NODES_>,GlobalOrdinal>(keys, gids, comm);
    for (GlobalOrdinal edge = 0; edge < numEdges; ++edge) {
      edgeNodeToEdgeMap_[keys[edge]] = gids[edge];
    }
    if (dim > 2) {
      keys.clear();
      for (auto pair : faceNodeToFaceMap_) {
        keys.push_back(pair.first);
      }
      gids.resize(numFaces);
      TEUCHOS_TEST_FOR_EXCEPTION(keys.size() != gids.size(), std::logic_error,
          		                     "Error, keys() != gids.size().");
      numUniqueFaces_ = Zoltan2::findUniqueGids<std::array<GlobalOrdinal,MAX_SUBCELL_CORNER_NODES_>,GlobalOrdinal>(keys, gids, comm);
      for (GlobalOrdinal face = 0; face < numFaces; ++face) {
        faceNodeToFaceMap_[keys[face]] = gids[face];
      }
    }
  }

#else
  if (dim > 1) {
	keys.clear();
	for (auto pair : edgeNodeToEdgeMap_) {
	  keys.push_back(pair.first);
	}
	for (GlobalOrdinal edge = 0; edge < numEdges; ++edge) {
	  edgeNodeToEdgeMap_[keys[edge]] = edge;
	}
	if (dim > 2) {
      keys.clear();
		for (auto pair : faceNodeToFaceMap_) {
		  keys.push_back(pair.first);
		}
	  for (GlobalOrdinal face = 0; face < numFaces; ++face) {
	    faceNodeToFaceMap_[keys[face]] = face;
	  }
	}
  }


#endif

}

void IOSSConnManager::buildConnectivity(const panzer::FieldPattern & fp) {

  //Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  //out.setShowProcRank(true);

  std::string blockName;
  int blockConnectivitySize = 0;
  std::vector<int> blockConnectivity;

  // Check that the FieldPattern has a topology that is compatible with all elements;
  TEUCHOS_TEST_FOR_EXCEPTION(!compatibleTopology(fp), std::logic_error,
		     		            "Error, fp must have a topology compatible with all elements.");

  // Get element info from IOSS database and build a local element mapping.
  buildLocalElementMapping();

  // build iossConnectivity_
  for (Ioss::ElementBlock * iossElementBlock : iossElementBlocks_) {

	// Ioss::GroupingEntity::get_field_data() will resize
	// blockConnectivity up or down as needed.
    iossElementBlock->get_field_data("connectivity", blockConnectivity);

    blockName = iossElementBlock->name();
    blockConnectivitySize = blockConnectivity.size();

    iossConnectivity_[blockName] = rcp(new std::vector<GlobalOrdinal>);
    iossConnectivity_[blockName]->resize(blockConnectivitySize);
    std::copy(blockConnectivity.begin(), blockConnectivity.end(), iossConnectivity_[blockName]->begin());

  }

  // Build mapping from edge and face corner node numbers to unique global edge and face numbers.
  buildEdgeFaceCornerNodeMapping();

  // Build sub cell ID counts and offsets
  //    ID counts = How many IDs belong on each subcell (number of mesh DOF used)
  //    Offset = What is starting index for subcell ID type?
  //             Global numbering goes like [node ids, edge ids, face ids, cell ids]
  LocalOrdinal nodeIdCnt=0, edgeIdCnt=0, faceIdCnt=0, cellIdCnt=0;
  GlobalOrdinal nodeOffset=0, edgeOffset=0, faceOffset=0, cellOffset=0;
  buildOffsetsAndIdCounts(fp, nodeIdCnt,  edgeIdCnt,  faceIdCnt,  cellIdCnt,
	                          nodeOffset, edgeOffset, faceOffset, cellOffset);

  //out << "node: count = " << nodeIdCnt << ", offset = " << nodeOffset << std::endl;
  //out << "edge: count = " << edgeIdCnt << ", offset = " << edgeOffset << std::endl;
  //out << "face: count = " << faceIdCnt << ", offset = " << faceOffset << std::endl;
  //out << "cell: count = " << cellIdCnt << ", offset = " << cellOffset << std::endl;

  // loop over elements and build global connectivity
  std::vector<std::string> elementBlockIds;
  getElementBlockIds(elementBlockIds);
  //RCP<std::vector<GlobalOrdinal>> elementBlock;
  std::size_t elmtLid = 0;

  for (std::string elementBlockId : elementBlockIds) {

	for(std::size_t elmtIdInBlock=0; elmtIdInBlock < getElementBlock(elementBlockId).size(); ++elmtIdInBlock) {

	  GlobalOrdinal numIds = 0;

      // get index into connectivity array
      elmtLidToConn_[elmtLid] = connectivity_.size();

      // add connectivities for sub cells
      numIds += addSubcellConnectivities(fp, elementBlockId, elmtIdInBlock, elmtLid, 0, nodeIdCnt, nodeOffset);
      numIds += addSubcellConnectivities(fp, elementBlockId, elmtIdInBlock, elmtLid, 1, edgeIdCnt, edgeOffset);
      numIds += addSubcellConnectivities(fp, elementBlockId, elmtIdInBlock, elmtLid, 2, faceIdCnt, faceOffset);
      numIds += addSubcellConnectivities(fp, elementBlockId, elmtIdInBlock, elmtLid, fp.getDimension(), cellIdCnt, cellOffset);

      connSize_[elmtLid] = numIds;
      ++elmtLid;
    }
  }

}

void IOSSConnManager::buildOffsetsAndIdCounts(const panzer::FieldPattern & fp,
                                  LocalOrdinal & nodeIdCnt, LocalOrdinal & edgeIdCnt,
                                  LocalOrdinal & faceIdCnt, LocalOrdinal & cellIdCnt,
                                  GlobalOrdinal & nodeOffset, GlobalOrdinal & edgeOffset,
                                  GlobalOrdinal & faceOffset, GlobalOrdinal & cellOffset) const {

  GlobalOrdinal nodeId;
  GlobalOrdinal blockmaxNodeId;
  GlobalOrdinal localmaxNodeId = -1;
  GlobalOrdinal maxNodeId = -1;

  GlobalOrdinal maxEdgeId = -1;
  GlobalOrdinal maxFaceId = -1;

  GlobalOrdinal blockmaxElementId;
  GlobalOrdinal localmaxElementId = std::numeric_limits<GlobalOrdinal>::min();
  GlobalOrdinal maxElementId;

  const CellTopologyData * baseTopologyData = fp.getCellTopology().getBaseCellTopologyData();

  int nodesPerElement = fp.getSubcellCount(0);
  int cornerNodesPerElement = int(baseTopologyData->subcell_count[0]);
  bool extended = nodesPerElement > cornerNodesPerElement;

  // get the local maximum node, edge, and face numbers
  std::vector<std::string> elementBlockIds;
  getElementBlockIds(elementBlockIds);

  int elementsInBlock;
  Teuchos::RCP<std::vector<GlobalOrdinal>> blockConnectivity;
  for (std::string elementBlockId : elementBlockIds) {
	std::vector<LocalOrdinal> elementBlock = *(elementBlocks_.find(elementBlockId)->second);
	blockmaxElementId = -1;
	for (auto element : elementBlock) {
	  if (element > blockmaxElementId) {
	    blockmaxElementId = elmtLidToGid_[element];
	  }
	}
	blockConnectivity = iossConnectivity_.find(elementBlockId)->second;
	int blockConnectivitySize = blockConnectivity->size();
	if (blockConnectivitySize > 0) {
	  const Ioss::ElementTopology * iossElementBlockTopology = iossElementBlockTopologies_.find(elementBlockId)->second;
	  int iossNodesPerElement = iossElementBlockTopology->number_nodes();
	  if (!extended && (cornerNodesPerElement < iossNodesPerElement)) {
	    // If the FieldPattern is not extended, but the Ioss::ElementTopology is, then
		// we must ignore the non-corner nodes.
		elementsInBlock = blockConnectivitySize / iossNodesPerElement;
		TEUCHOS_TEST_FOR_EXCEPTION(blockConnectivitySize % iossNodesPerElement != 0, std::logic_error,
		   "Error, Size of connectivity for this block not a multiple of the number of elements in the block.");
		blockmaxNodeId = -1;
		for (int element = 0; element < elementsInBlock; ++element) {
		  for (int node = 0; node < cornerNodesPerElement; ++node) {
		    nodeId = *(blockConnectivity->begin() + nodesPerElement*element + iossElementBlockTopology->element_connectivity()[node]);
		    if (nodeId > blockmaxNodeId) {
		      blockmaxNodeId = nodeId;
		    }
		  }
		}
	  }
	  else {
	    // We can use the entire blockConnectivity vector if the FieldPattern is extended
		// or if neither the FieldPattern nor the Ioss::ElementTopology is extended.
		blockmaxNodeId = *(std::max_element(blockConnectivity->begin(), blockConnectivity->end()));
	  }
          if (blockmaxNodeId > localmaxNodeId)
            localmaxNodeId = blockmaxNodeId;
	}
	if (blockmaxElementId > localmaxElementId)
          localmaxElementId = blockmaxElementId;
  }

#ifdef HAVE_MPI
  Ioss::ParallelUtils util = iossMeshDB_->util();
  MPI_Comm communicator = util.communicator();
  Teuchos::MpiComm<int> comm(communicator);

  Teuchos::reduceAll<int, GlobalOrdinal>(comm, Teuchos::REDUCE_MAX, 1, &localmaxNodeId, &maxNodeId);
  Teuchos::reduceAll<int, GlobalOrdinal>(comm, Teuchos::REDUCE_MAX, 1, &localmaxElementId, &maxElementId);
#else
  maxNodeId = localmaxNodeId;
  maxNodeId = localmaxElementId;
#endif

  int patternDim = fp.getDimension();
  switch(patternDim) {
    case 0:
      maxEdgeId = -1;
      maxFaceId = -1;
      break;
    case 1:
      maxEdgeId = maxElementId - 1;
      maxFaceId = -1;
      break;
    case 2:
      maxEdgeId = numUniqueEdges_ - 1;
      maxFaceId = maxElementId - 1;
      break;
    default:
      maxEdgeId = numUniqueEdges_ - 1;
      maxFaceId = numUniqueFaces_ - 1;
  };

  // compute ID counts for each sub cell type
  switch (patternDim)
  {
    case 3:
      faceIdCnt = fp.getSubcellIndices(2,0).size();
      // Intentional fall-through.
    case 2:
      edgeIdCnt = fp.getSubcellIndices(1,0).size();
      // Intentional fall-through.
    case 1:
      nodeIdCnt = fp.getSubcellIndices(0,0).size();
      cellIdCnt = fp.getSubcellIndices(patternDim,0).size();
      break;
    case 0:
    default:
      TEUCHOS_ASSERT(false);
  }

  // compute offsets for each sub cell type
  nodeOffset = 0;
  edgeOffset = nodeOffset+(maxNodeId+1)*nodeIdCnt;
  faceOffset = edgeOffset+(maxEdgeId+1)*edgeIdCnt;
  cellOffset = faceOffset+(maxFaceId+1)*faceIdCnt;

  // sanity check
  TEUCHOS_ASSERT(nodeOffset <= edgeOffset
	             && edgeOffset <= faceOffset
	             && faceOffset <= cellOffset);

}

typename IOSSConnManager::LocalOrdinal IOSSConnManager::addSubcellConnectivities(
		const panzer::FieldPattern & fp, std::string & blockId, std::size_t elmtIdInBlock, std::size_t elmtLid,
		unsigned subcellRank, LocalOrdinal idCnt, GlobalOrdinal offset) {

   if(idCnt<=0)
     return 0 ;

   int dim = fp.getDimension();
   TEUCHOS_TEST_FOR_EXCEPTION(int(subcellRank) > dim, std::logic_error,
          "Error, subcell rank (" << subcellRank << ") must not be greater than dimension (" << dim << ").");

   std::size_t subcellsPerElement = fp.getSubcellCount(subcellRank);

   const Ioss::ElementTopology * iossElementBlockTopology = iossElementBlockTopologies_.find(blockId)->second;
   int iossNodesPerElement = iossElementBlockTopology->number_nodes();


   LocalOrdinal numIds = 0;
   GlobalOrdinal subcellId = 0;

   // loop over all subcells
   for(std::size_t subcell=0; subcell < subcellsPerElement; ++subcell) {
	 if (subcellRank == 0) {
	   subcellId = (*(iossConnectivity_[blockId]))[iossNodesPerElement*elmtIdInBlock+iossElementBlockTopology->element_connectivity()[subcell]];
	 }
	 else {
	   if (int(subcellRank) == dim) {
	     subcellId = elmtLidToGid_[elmtLid];
	   }
	   else {
	     if (subcellRank == 1) {
	       subcellId = edgeNodeToEdgeMap_[(*(elementEdgeNodes_[blockId]))[elmtIdInBlock][subcell]];
	     }
	     else if (subcellRank == 2) {
	       subcellId = faceNodeToFaceMap_[(*(elementFaceNodes_[blockId]))[elmtIdInBlock][subcell]];
	     }
	     else {
	       TEUCHOS_TEST_FOR_EXCEPTION(false, std::logic_error,
	    	           "Error, unsupported combination of subcell rank (" << subcellRank << ") and dimension (" << dim << ").");
	     }
	   }
	 }

	 // loop over all ids for subcell
	 for(LocalOrdinal i=0; i < idCnt; i++) {
       connectivity_.push_back(offset+idCnt*subcellId+i);
	 }
     numIds += idCnt;
   }

   return numIds;
}

bool IOSSConnManager::compatibleTopology(const panzer::FieldPattern & fp) const {

   bool compatible = true;

   const CellTopologyData * baseTopologyData = fp.getCellTopology().getBaseCellTopologyData();

   int numNodes = fp.getSubcellCount(0);
   int numCornerNodes = int(baseTopologyData->subcell_count[0]);
   bool extended = numNodes > numCornerNodes;

   int numEdges = 0;
   int numFaces = 0;
   int numEdgeCornerNodes, numFaceCornerNodes;
   int uniqueNumEdgeCornerNodes = 1 + MAX_SUBCELL_CORNER_NODES_; // Needed for Zoltan2::findUniqueGids()
   int uniqueNumFaceCornerNodes = 1 + MAX_SUBCELL_CORNER_NODES_; // Needed for Zoltan2::findUniqueGids()

   std::vector<int> numCornerNodesOnEdges, numCornerNodesOnFaces, numEdgesOnFaces;

   if (baseTopologyData->dimension > 1) {

     numEdges = int(baseTopologyData->subcell_count[1]);

     numCornerNodesOnEdges.resize(numEdges);
     for (int edge = 0; edge < numEdges; ++edge) {
       numEdgeCornerNodes = int(baseTopologyData->subcell[1][edge].topology->subcell_count[0]);
       numCornerNodesOnEdges[edge] = numEdgeCornerNodes;
       if (edge == 0) {
         uniqueNumEdgeCornerNodes = numEdgeCornerNodes;
       }
       else {
         compatible &= numEdgeCornerNodes == uniqueNumEdgeCornerNodes;
       }
     }

     TEUCHOS_TEST_FOR_EXCEPTION(!compatible, std::logic_error,
       "Error, currently IOSSConnManager only supports elements for which all edges have the same number of corner nodes.");

     TEUCHOS_TEST_FOR_EXCEPTION(uniqueNumEdgeCornerNodes > MAX_SUBCELL_CORNER_NODES_, std::logic_error,
            "Error, currently IOSSConnManager only supports element edges with at most " << MAX_SUBCELL_CORNER_NODES_ << " corner nodes.");


     if (baseTopologyData->dimension > 2) {

       numFaces = int(baseTopologyData->subcell_count[2]);

       numCornerNodesOnFaces.resize(numFaces);
       for (int face = 0; face < numFaces; ++face) {
    	 numFaceCornerNodes = int(baseTopologyData->subcell[2][face].topology->subcell_count[0]);
         numCornerNodesOnFaces[face] = numFaceCornerNodes;
         if (face == 0) {
           uniqueNumFaceCornerNodes = numFaceCornerNodes;
         }
         else {
           compatible &= numFaceCornerNodes == uniqueNumFaceCornerNodes;
         }
       }

       TEUCHOS_TEST_FOR_EXCEPTION(!compatible, std::logic_error,
         "Error, currently IOSSConnManager only supports elements for which all faces have the same number of corner nodes.");

       TEUCHOS_TEST_FOR_EXCEPTION(uniqueNumEdgeCornerNodes > MAX_SUBCELL_CORNER_NODES_, std::logic_error,
              "Error, currently IOSSConnManager only supports element faces with at most " << MAX_SUBCELL_CORNER_NODES_ << " corner nodes.");

       numEdgesOnFaces.resize(numFaces);
       for (int face = 0; face < numFaces; ++face) {
         numEdgesOnFaces[face] = int(baseTopologyData->subcell[2][face].topology->subcell_count[1]);
       }

     }
   }

   for (Ioss::ElementBlock * iossElementBlock : iossElementBlocks_) {

	 const Ioss::ElementTopology * topology = iossElementBlockTopologies_.find(iossElementBlock->name())->second;

     // test same dimension
     std::size_t dim = topology->spatial_dimension();
     compatible &= (dim==(std::size_t) baseTopologyData->dimension);

     compatible &= topology->number_corner_nodes()==numCornerNodes;
     if (dim > 1) {
       compatible &= topology->number_edges()==numEdges;
       if (dim > 2) {
         compatible &= topology->number_faces()==numFaces;
       }
     }

     // If the field pattern has more nodes than the base cell topology,
     // make sure the Ioss element has the same number of nodes.
     if (extended) {
    	 compatible &= topology->number_nodes()==numNodes;
     }

     // Check that base topology of field pattern has same edge and face topologies as
     // base topology of the Ioss element.
     if (dim > 1) {
       for (int edge = 0; edge < numEdges; ++edge) {
         compatible &= topology->edge_type(edge)->number_corner_nodes()==numCornerNodesOnEdges[edge];
       }

       if (dim > 2) {
         for (int face = 0; face < numFaces; ++face) {
    	   compatible &= topology->face_type(face)->number_corner_nodes()==numCornerNodesOnFaces[face];
    	   compatible &= topology->face_type(face)->number_edges()==numEdgesOnFaces[face];
         }
       }
     }

   }

   return compatible;
}

void IOSSConnManager::createTopologyMapping() {

  // Miscellaneous
  // -------------

  //iossToShardsTopology_["unknown"] =
  //    shards::CellTopology(shards::getCellTopologyData<shards::>());

  //iossToShardsTopology_["super"] =
  //      shards::CellTopology(shards::getCellTopologyData<shards::>());

  // 0-D Continuum
  // -------------

  iossToShardsTopology_["node"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Node>());

  // 2-D Continuum
  // -------------

  iossToShardsTopology_["tri3"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3>>());

  iossToShardsTopology_["tri4"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Triangle<4>>());

  //iossToShardsTopology_["tri4a"] =
  //    shards::CellTopology(shards::getCellTopologyData<shards::>());

  iossToShardsTopology_["tri6"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Triangle<6>>());

  //iossToShardsTopology_["tri7"] =
  //    shards::CellTopology(shards::getCellTopologyData<shards::>());

  iossToShardsTopology_["quad4"] =
    shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4>>());

  iossToShardsTopology_["quad8"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<8>>());

  iossToShardsTopology_["quad9"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<9>>());

 // 3-D Continuum
 // -------------

  iossToShardsTopology_["tetra4"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4>>());

  //iossToShardsTopology_["tetra7"] =
  //    shards::CellTopology(shards::getCellTopologyData<shards::>());

  //iossToShardsTopology_["tetra8"] =
  //    shards::CellTopology(shards::getCellTopologyData<shards::>());

  iossToShardsTopology_["tetra10"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<10>>());

  iossToShardsTopology_["tetra11"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<11>>());

  //iossToShardsTopology_["tetra14"] =
  //    shards::CellTopology(shards::getCellTopologyData<shards::>());

  //iossToShardsTopology_["tetra15"] =
  //    shards::CellTopology(shards::getCellTopologyData<shards::>());

  iossToShardsTopology_["pyramid5"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<5>>());

  iossToShardsTopology_["pyramid13"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<13>>());

  iossToShardsTopology_["pyramid14"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Pyramid<14>>());

  iossToShardsTopology_["wedge6"] =
        shards::CellTopology(shards::getCellTopologyData<shards::Wedge<6>>());

  iossToShardsTopology_["wedge15"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Wedge<15>>());

  //iossToShardsTopology_["wedge16"] =
  //    shards::CellTopology(shards::getCellTopologyData<shards::>());

  iossToShardsTopology_["wedge18"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Wedge<18>>());

  //iossToShardsTopology_["wedge20"] =
  //    shards::CellTopology(shards::getCellTopologyData<shards::>());

  //iossToShardsTopology_["wedge21"] =
  //    shards::CellTopology(shards::getCellTopologyData<shards::>());

  iossToShardsTopology_["hex8"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8>>());

  iossToShardsTopology_["hex20"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<20>>());

  iossToShardsTopology_["hex27"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<27>>());

  // 1-D Structural
  // --------------


  iossToShardsTopology_["bar2"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Beam<2>>());

  iossToShardsTopology_["bar3"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Beam<3>>());

  // 2-D Structural
  // --------------

  //iossToShardsTopology_["sphere"] =
  //    shards::CellTopology(shards::getCellTopologyData<shards::>());

  iossToShardsTopology_["trishell3"] =
       shards::CellTopology(shards::getCellTopologyData<shards::ShellTriangle<3>>());

  //iossToShardsTopology_["trishell4"] =
  //    shards::CellTopology(shards::getCellTopologyData<shards::>());

  iossToShardsTopology_["trishell6"] =
      shards::CellTopology(shards::getCellTopologyData<shards::ShellTriangle<6>>());

  //iossToShardsTopology_["trishell7"] =
  //    shards::CellTopology(shards::getCellTopologyData<shards::>());

  iossToShardsTopology_["shell4"] =
      shards::CellTopology(shards::getCellTopologyData<shards::ShellQuadrilateral<4>>());

  iossToShardsTopology_["shell8"] =
      shards::CellTopology(shards::getCellTopologyData<shards::ShellQuadrilateral<8>>());

  iossToShardsTopology_["shell9"] =
      shards::CellTopology(shards::getCellTopologyData<shards::ShellQuadrilateral<9>>());

  // Unsure
  // ------

  iossToShardsTopology_["edge2"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Line<2>>());

  iossToShardsTopology_["edge3"] =
      shards::CellTopology(shards::getCellTopologyData<shards::Line<3>>());

  //iossToShardsTopology_["edge2d2"] =
  //    shards::CellTopology(shards::getCellTopologyData<shards::>());

  //iossToShardsTopology_["edge2d3"] =
  //    shards::CellTopology(shards::getCellTopologyData<shards::>());

  iossToShardsTopology_["shellline2d2"] =
      shards::CellTopology(shards::getCellTopologyData<shards::ShellLine<2>>());

  iossToShardsTopology_["shellline2d3"] =
      shards::CellTopology(shards::getCellTopologyData<shards::ShellLine<3>>());

}

} // end namespace panzer_ioss



