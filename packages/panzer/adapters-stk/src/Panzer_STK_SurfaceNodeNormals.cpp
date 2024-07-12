// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_STK_SurfaceNodeNormals.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_IntegrationValues2.hpp"
#include "Panzer_STK_WorksetFactory.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/CreateAdjacentEntities.hpp>

#include "Shards_CellTopology.hpp"
//#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_CellTools_Serial.hpp"
#include "Teuchos_Assert.hpp"

namespace panzer_stk {

  void computeSidesetNodeNormals(std::unordered_map<unsigned,std::vector<double> >& normals,
				 const Teuchos::RCP<const panzer_stk::STK_Interface>& mesh,
				 const std::string& sidesetName,
				 const std::string& elementBlockName,
				 std::ostream* /* out */,
				 std::ostream* pout)
  {    
    using panzer::Cell;
    using panzer::NODE;
    using panzer::Dim;

    using Teuchos::RCP;

    panzer::MDFieldArrayFactory af("",true);
    
    RCP<stk::mesh::MetaData> metaData = mesh->getMetaData();
    RCP<stk::mesh::BulkData> bulkData = mesh->getBulkData();

    // Grab all nodes for a surface including ghosted to get correct contributions to normal average    
    stk::mesh::Part * sidePart = mesh->getSideset(sidesetName);
    stk::mesh::Part * elmtPart = mesh->getElementBlockPart(elementBlockName);
    stk::mesh::Selector sideSelector = *sidePart;
    stk::mesh::Selector blockSelector = *elmtPart;
    stk::mesh::Selector mySelector = metaData->universal_part() & blockSelector & sideSelector;
    std::vector<stk::mesh::Entity> sides;
    stk::mesh::get_selected_entities(mySelector,bulkData->buckets(metaData->side_rank()),sides);

    std::vector<std::size_t> missingElementIndices;
    std::vector<std::size_t> localSideTopoIDs;
    std::vector<stk::mesh::Entity> parentElements;
    panzer_stk::workset_utils::getUniversalSubcellElements(*mesh,elementBlockName,sides,localSideTopoIDs,parentElements,missingElementIndices);

    // Delete all sides whose neighboring element in elementBlockName is not in the current process
    std::vector<std::size_t>::reverse_iterator index;
    for(index=missingElementIndices.rbegin(); index != missingElementIndices.rend(); ++index) {
      sides.erase(sides.begin()+*index);
    }

    if (pout != NULL) {
      for (std::size_t i=0; i < localSideTopoIDs.size(); ++i) {
	*pout << "parent element: "
	      << " gid(" << bulkData->identifier(parentElements[i]) << ")"
	      << ", local_face(" << localSideTopoIDs[i] << ")"
	      << std::endl;
      }
    }

    // Do a single element at a time so that we don't have to batch up
    // into faces

    // maps a panzer local element id to a list of normals
    std::unordered_map<unsigned,std::vector<double> > nodeNormals;
    
    TEUCHOS_ASSERT(sides.size() == localSideTopoIDs.size());
    TEUCHOS_ASSERT(localSideTopoIDs.size() == parentElements.size());

    RCP<const shards::CellTopology> parentTopology = mesh->getCellTopology(elementBlockName);
    //Intrepid2::DefaultCubatureFactory cubFactory;
    int cubDegree = 1;

    std::vector<stk::mesh::Entity>::const_iterator side = sides.begin();
    std::vector<std::size_t>::const_iterator sideID = localSideTopoIDs.begin();
    std::vector<stk::mesh::Entity>::const_iterator parentElement = parentElements.begin();

    // KK: invoke serial interface; cubDegree is 1 and integration point is one 
    //     for debugging statement, use max dimension
    auto side_parametrization = Intrepid2::RefSubcellParametrization<Kokkos::HostSpace>::get(2,parentTopology->getKey());
    Kokkos::DynRankView<double,Kokkos::HostSpace> normal_at_point("normal",3); // parentTopology->getDimension());
    for ( ; sideID != localSideTopoIDs.end(); ++side,++sideID,++parentElement) {
    
      std::vector<stk::mesh::Entity> elementEntities;
      elementEntities.push_back(*parentElement); // notice this is size 1!
      PHX::MDField<double,panzer::Cell,panzer::NODE,panzer::Dim> nodes 
          = af.buildStaticArray<double,Cell,NODE,Dim>("",elementEntities.size(), parentTopology->getNodeCount(), mesh->getDimension());
      auto node_view = nodes.get_view();
      mesh->getElementNodesNoResize(elementEntities,elementBlockName,node_view);
      
      panzer::CellData sideCellData(1,*sideID,parentTopology); // this is size 1 because elementEntties is size 1!
      RCP<panzer::IntegrationRule> ir = Teuchos::rcp(new panzer::IntegrationRule(cubDegree,sideCellData));

      panzer::IntegrationValues2<double> iv("",true);
      iv.setupArrays(ir);
      iv.evaluateValues(nodes);
      
      // KK: use serial interface; jac_at_point (D,D) from (C,P,D,D)
      {
        auto jac_at_point = Kokkos::subview(iv.jac.get_view(), 0, 0, Kokkos::ALL(), Kokkos::ALL());
	auto jac_at_point_h = Kokkos::create_mirror_view(jac_at_point);
	Kokkos::deep_copy(jac_at_point_h, jac_at_point);
        Intrepid2::Impl::
          CellTools::Serial::getPhysicalSideNormal(normal_at_point, side_parametrization, jac_at_point_h, *sideID);
      }

      if (pout != NULL) {
      *pout << "element normals: "
	    << "gid(" << bulkData->identifier(*parentElement) << ")"
	    << ", normal(" << normal_at_point(0) << "," << normal_at_point(1) << "," << normal_at_point(2) << ")"
	    << std::endl;
      }

      // loop over nodes in nodes in side and add normal contribution for averaging
      const size_t numNodes = bulkData->num_nodes(*side);
      stk::mesh::Entity const* nodeRelations = bulkData->begin_nodes(*side);
      for (size_t n=0; n<numNodes; ++n) {
        stk::mesh::Entity node = nodeRelations[n];
	for (unsigned dim = 0; dim < parentTopology->getDimension(); ++dim) {
	  nodeNormals[bulkData->identifier(node)].push_back(normal_at_point(dim));
	}
      }

    }

    // Now do the averaging of contributions
    //std::unordered_map<unsigned,std::vector<double> > normals;
    for (std::unordered_map<unsigned,std::vector<double> >::const_iterator node = nodeNormals.begin(); node != nodeNormals.end(); ++node) {

      TEUCHOS_ASSERT( (node->second.size() % parentTopology->getDimension()) == 0);

      int numSurfaceContributions = node->second.size() / parentTopology->getDimension();
      std::vector<double> contribArea(numSurfaceContributions);
      double totalArea = 0.0;
      for (int surface = 0; surface < numSurfaceContributions; ++surface) {

	double sum = 0.0;
	for (unsigned dim = 0; dim < parentTopology->getDimension(); ++dim)
	  sum += (node->second[surface*parentTopology->getDimension() + dim]) * 
	    (node->second[surface*parentTopology->getDimension() + dim]);
	
	contribArea[surface] = std::sqrt(sum);

	totalArea += contribArea[surface];
      }

      // change the contribArea to the scale factor for each contribution
      for (std::size_t i = 0; i < contribArea.size(); ++i)
	contribArea[i] /= totalArea;

      // loop over contributions and compute final normal values
      normals[node->first].resize(parentTopology->getDimension());
      for (unsigned dim = 0; dim < parentTopology->getDimension(); ++dim) {
	normals[node->first][dim] = 0.0;
	for (int surface = 0; surface < numSurfaceContributions; ++surface) {
	  normals[node->first][dim] += node->second[surface*parentTopology->getDimension() + dim] * contribArea[surface] / totalArea;
	}
      }

      if (pout != NULL) {
	*pout << "surface normal before normalization: " 
	      << "gid(" << node->first << ")"
	      << ", normal(" << normals[node->first][0] << "," << normals[node->first][1] << "," << normals[node->first][2] << ")"
	      << std::endl;
      }
	    
      double sum = 0.0;
      for (unsigned dim = 0; dim < parentTopology->getDimension(); ++dim)
	sum += normals[node->first][dim] * normals[node->first][dim];

      for (unsigned dim = 0; dim < parentTopology->getDimension(); ++dim)
	normals[node->first][dim] /= std::sqrt(sum);
      
      if (pout != NULL) {
	*pout << "surface normal after normalization: " 
	      << "gid(" << node->first << ")"
	      << ", normal(" 
	      << normals[node->first][0] << "," 
	      << normals[node->first][1] << "," 
	      << normals[node->first][2] << ")"
	      << std::endl;
      }

    }
    
  }

  void computeSidesetNodeNormals(std::unordered_map<std::size_t,Kokkos::DynRankView<double,PHX::Device> >& normals,
				 const Teuchos::RCP<const panzer_stk::STK_Interface>& mesh,
				 const std::string& sidesetName,
				 const std::string& elementBlockName,
				 std::ostream* out,
				 std::ostream* pout)
  {    
    using Teuchos::RCP;
    
    std::unordered_map<unsigned,std::vector<double> > nodeEntityIdToNormals;
    
    computeSidesetNodeNormals(nodeEntityIdToNormals,mesh,sidesetName,elementBlockName,out,pout);

    RCP<stk::mesh::MetaData> metaData = mesh->getMetaData();
    RCP<stk::mesh::BulkData> bulkData = mesh->getBulkData();

    // Grab all nodes for a surface including ghosted to get correct contributions to normal average    
    stk::mesh::Part * sidePart = mesh->getSideset(sidesetName);
    stk::mesh::Part * elmtPart = mesh->getElementBlockPart(elementBlockName);
    stk::mesh::Selector sideSelector = *sidePart;
    stk::mesh::Selector blockSelector = *elmtPart;
    stk::mesh::Selector mySelector = metaData->universal_part() & blockSelector & sideSelector;
    std::vector<stk::mesh::Entity> sides;
    stk::mesh::get_selected_entities(mySelector,bulkData->buckets(metaData->side_rank()),sides);

    RCP<const shards::CellTopology> parentTopology = mesh->getCellTopology(elementBlockName);

    std::vector<std::size_t> missingElementIndices;
    std::vector<std::size_t> localSideTopoIDs;
    std::vector<stk::mesh::Entity> parentElements;
    panzer_stk::workset_utils::getUniversalSubcellElements(*mesh,elementBlockName,sides,localSideTopoIDs,parentElements,missingElementIndices);

    // Delete all sides whose neighboring element in elementBlockName is not in the current process
    std::vector<std::size_t>::reverse_iterator index;
    for(index=missingElementIndices.rbegin(); index != missingElementIndices.rend(); ++index) {
      sides.erase(sides.begin()+*index);
    }
    
    std::vector<stk::mesh::Entity>::const_iterator side = sides.begin();
    std::vector<std::size_t>::const_iterator sideID = localSideTopoIDs.begin();
    std::vector<stk::mesh::Entity>::const_iterator parentElement = parentElements.begin();
    for ( ; sideID != localSideTopoIDs.end(); ++side,++sideID,++parentElement) {
    
      // loop over nodes in nodes in side element
      const size_t numNodes = bulkData->num_nodes(*parentElement);
      stk::mesh::Entity const* nodeRelations = bulkData->begin_nodes(*parentElement);

      normals[mesh->elementLocalId(*parentElement)] = Kokkos::DynRankView<double,PHX::Device>("normals",numNodes,parentTopology->getDimension());
      auto normals_h = Kokkos::create_mirror_view(normals[mesh->elementLocalId(*parentElement)]);
      for (size_t nodeIndex=0; nodeIndex<numNodes; ++nodeIndex) {
        stk::mesh::Entity node = nodeRelations[nodeIndex];
	// if the node is on the sideset, insert, otherwise set normal
	// to zero (it is an interior node of the parent element).
	if (nodeEntityIdToNormals.find(bulkData->identifier(node)) != nodeEntityIdToNormals.end()) { 
	  for (unsigned dim = 0; dim < parentTopology->getDimension(); ++dim) {
	    normals_h(nodeIndex,dim) = (nodeEntityIdToNormals[bulkData->identifier(node)])[dim];
	  }
	}
	else {
	  for (unsigned dim = 0; dim < parentTopology->getDimension(); ++dim) {
	    normals_h(nodeIndex,dim) = 0.0;
	  }
	}
      }
      Kokkos::deep_copy(normals[mesh->elementLocalId(*parentElement)], normals_h);
    }

  }

}
