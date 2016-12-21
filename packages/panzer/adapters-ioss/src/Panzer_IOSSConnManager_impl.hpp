// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include "Panzer_IOSSConnManager.hpp"

#include <string>
#include <vector>
#include <algorithm>

// Teuchos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_TestForException.hpp"

// Kokkos includes
#include "Kokkos_DynRankView.hpp"

// Shards includes
#include "Shards_CellTopology.hpp"

// Panzer includes
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
//#include "Panzer_STK_PeriodicBC_Matcher.hpp"
//#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_NodalFieldPattern.hpp"

// Ioss includes
#include "Ioss_DatabaseIO.h"
#include "Ioss_GroupingEntity.h"
#include "Ioss_Region.h"
#include "Ioss_EntityBlock.h"
#include "Ioss_ElementTopology.h"
#include "Ioss_NodeBlock.h"
#include "Ioss_ElementBlock.h"

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_ioss {

template <typename GO>
IOSSConnManager<GO>::IOSSConnManager(Ioss::DatabaseIO * iossMeshDB)
   : iossMeshDB_(iossMeshDB), elementBlocks_(), neighborElementBlocks_(),
	 elmtLidToGid_(), elmtLidToConn_(), connSize_(),
	 connectivity_(), ownedElementCount_(0), sidesetsToAssociate_(),
	 sidesetYieldedAssociations_(), elmtToAssociatedElmts_()
{

  std::string error_message;
  int bad_count;

  // Check for a non-null database
  TEUCHOS_TEST_FOR_EXCEPTION(iossMeshDB_ == nullptr, std::logic_error, "Error, iossMeshDB is nullptr.");
  TEUCHOS_TEST_FOR_EXCEPTION(iossMeshDB_ == NULL, std::logic_error, "Error, iossMeshDB is NULL.");

  // Check that the database is OK and that the file can be opened.
  bool status_ok = iossMeshDB->ok(true, &error_message, &bad_count);
  TEUCHOS_TEST_FOR_EXCEPTION(!status_ok, std::logic_error,
		  "Error, Ioss::DatabaseIO::ok() failed.\n error_message = "
		  << error_message << "\nbad_count = " << bad_count << "\n");

  // Get the metadata
  iossRegion_ = rcp(new Ioss::Region(iossMeshDB_, "iossRegion_"));
  iossNodeBlocks_ = iossRegion_->get_node_blocks();
  iossElementBlocks_ = iossRegion_->get_element_blocks();

  // Create iossToShardsTopology_ map;
  createTopologyMapping();

}

template <typename GO>
std::string IOSSConnManager<GO>::getBlockId(LocalOrdinal localElmtId) const {
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

template <typename GO>
void IOSSConnManager<GO>::clearLocalElementMapping()
{

   elementBlocks_.clear();
   elmtLidToGid_.clear();
   elmtLidToConn_.clear();
   connSize_.clear();
   //elmtToAssociatedElmts_.clear();
}

template <typename GO>
void IOSSConnManager<GO>::buildLocalElementMapping() {

  std::string name;
  std::vector<LocalOrdinal> indices;
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
    for (LocalOrdinal element : indices) {
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

template <typename GO>
void IOSSConnManager<GO>::buildConnectivity(const panzer::FieldPattern & fp) {


  //Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  //out.setShowProcRank(true);

  std::string blockName;
  int blockConnectivitySize = 0;
  std::vector<GlobalOrdinal> blockConnectivity;

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
  RCP<std::vector<GlobalOrdinal>> elementBlock;
  std::size_t elmtLid = 0;

  for (std::string elementBlockId : elementBlockIds) {

	elementBlock = iossConnectivity_[elementBlockId];

	for(std::size_t elmtIdInBlock=0; elmtIdInBlock < getElementBlock(elementBlockId).size(); ++elmtIdInBlock) {

	  GlobalOrdinal numIds = 0;

      // get index into connectivity array
      elmtLidToConn_[elmtLid] = connectivity_.size();

      // add connectivities for sub cells
      numIds += addSubcellConnectivities(fp, elementBlockId, elmtIdInBlock, 0, nodeIdCnt, nodeOffset);
      numIds += addSubcellConnectivities(fp, elementBlockId, elmtIdInBlock, 1, edgeIdCnt, edgeOffset);
      numIds += addSubcellConnectivities(fp, elementBlockId, elmtIdInBlock, 2, faceIdCnt, faceOffset);
      numIds += addSubcellConnectivities(fp, elementBlockId, elmtIdInBlock, 3, faceIdCnt, faceOffset);

      connSize_[elmtLid] = numIds;
      ++elmtLid;
    }
  }

}

template <typename GO>
void IOSSConnManager<GO>::buildOffsetsAndIdCounts(const panzer::FieldPattern & fp,
                                  LocalOrdinal & nodeIdCnt, LocalOrdinal & edgeIdCnt,
                                  LocalOrdinal & faceIdCnt, LocalOrdinal & cellIdCnt,
                                  GlobalOrdinal & nodeOffset, GlobalOrdinal & edgeOffset,
                                  GlobalOrdinal & faceOffset, GlobalOrdinal & cellOffset) const {

  GlobalOrdinal blockmaxNodeId;

  GlobalOrdinal localmaxNodeId = -1;
  GlobalOrdinal localmaxEdgeId = -1;
  GlobalOrdinal localmaxFaceId = -1;

  GlobalOrdinal maxNodeId = -1;
  GlobalOrdinal maxEdgeId = -1;
  GlobalOrdinal maxFaceId = -1;

  // get the local maximum node, edge, and face numbers
  std::vector<std::string> elementBlockIds;
  getElementBlockIds(elementBlockIds);
  for (std::string elementBlockId : elementBlockIds) {
	if (iossConnectivity_.find(elementBlockId)->second->size() > 0) {
      blockmaxNodeId = *(std::max_element(iossConnectivity_.find(elementBlockId)->second->begin(), iossConnectivity_.find(elementBlockId)->second->end()));
      if (blockmaxNodeId > localmaxNodeId)
        localmaxNodeId = blockmaxNodeId;
	}
  }

  // Placeholder code. Need access to global communicator for reduceAll.
  // Requiring the communicator changes the class's interface.
  // get the global maximum node, edge, and face numbers
  //Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &localmaxNodeId, &maxNodeId); // Need to cast to char*
  //Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &localmaxEdgeId, &maxEdgeId); // Need to cast to char*
  //Teuchos::reduceAll(comm, Teuchos::REDUCE_MAX, 1, &localmaxFaceId, &maxFaceId); // Need to cast to char*
  maxNodeId = localmaxNodeId;
  maxEdgeId = localmaxEdgeId;
  maxFaceId = localmaxFaceId;


  // compute ID counts for each sub cell type
  int patternDim = fp.getDimension();
  switch(patternDim) {
    case 3:
	  faceIdCnt = fp.getSubcellIndices(2,0).size();
	case 2:
	  edgeIdCnt = fp.getSubcellIndices(1,0).size();
	case 1:
	  nodeIdCnt = fp.getSubcellIndices(0,0).size();
	  cellIdCnt = fp.getSubcellIndices(patternDim,0).size();
	  break;
	 case 0:
	 default:
	   TEUCHOS_ASSERT(false);
	 };

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

template <typename GO>
typename IOSSConnManager<GO>::LocalOrdinal IOSSConnManager<GO>::addSubcellConnectivities(
		const panzer::FieldPattern & fp, std::string & blockId, std::size_t elmtIdInBlock, unsigned subcellRank, LocalOrdinal idCnt, GlobalOrdinal offset) {

   if(idCnt<=0)
     return 0 ;

   if(subcellRank > 0)
     return 0; // Currently only set up for nodes, not edges, faces, or cells
               // since iossConnectivity_ only has node numbers.

   std::size_t subcellsPerElement = fp.getSubcellCount(subcellRank);

   // loop over all nodes
   LocalOrdinal numIds = 0;
   GlobalOrdinal subcellId;

   for(std::size_t subcell=0; subcell < subcellsPerElement; ++subcell) {
	 // Is this what we want, or do we want contiguous ids with a mapping back to application subcell numbers?
	 subcellId = (*(iossConnectivity_[blockId]))[subcellsPerElement*elmtIdInBlock+subcell];
	 // loop over all ids for subcell
	 for(LocalOrdinal i=0; i < idCnt; i++) {
       connectivity_.push_back(offset+idCnt*subcellId+i);
	 }
     numIds += idCnt;
   }

   return numIds;
}

template <typename GO>
void IOSSConnManager<GO>::createTopologyMapping() {

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



