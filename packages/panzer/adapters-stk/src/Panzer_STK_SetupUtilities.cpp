// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Teuchos_Assert.hpp"

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace panzer_stk {

Teuchos::RCP<std::vector<panzer::Workset> >  
buildWorksets(const panzer_stk::STK_Interface & mesh,
              const std::string & eBlock,
              const panzer::WorksetNeeds & needs)
{
  using namespace workset_utils;

  std::vector<std::string> element_blocks;

  std::vector<std::size_t> local_cell_ids;
  Kokkos::DynRankView<double,PHX::Device> cell_node_coordinates;

  getIdsAndNodes(mesh, eBlock, local_cell_ids, cell_node_coordinates);

  // only build workset if there are elements to worry about
  // this may be processor dependent, so an element block
  // may not have elements and thus no contribution
  // on this processor
  return panzer::buildWorksets(needs, eBlock, local_cell_ids, cell_node_coordinates);
}

Teuchos::RCP<std::vector<panzer::Workset> >  
buildWorksets(const panzer_stk::STK_Interface & mesh,
              const panzer::WorksetNeeds & needs,
              const std::string & sideset,
              const std::string & eBlock,
              bool useCascade)
{
  using namespace workset_utils;
  using Teuchos::RCP;

  std::vector<stk::mesh::Entity> sideEntities; 

  try {
     // grab local entities on this side
     // ...catch any failure...primarily wrong side set and element block info
     if(!useCascade)
        mesh.getMySides(sideset,eBlock,sideEntities);
     else
        mesh.getAllSides(sideset,eBlock,sideEntities);
  } 
  catch(STK_Interface::SidesetException & e) {
     std::stringstream ss;
     std::vector<std::string> sideSets; 
     mesh.getSidesetNames(sideSets);
 
     // build an error message
     ss << e.what() << "\nChoose one of:\n";
     for(std::size_t i=0;i<sideSets.size();i++) 
        ss << "\"" << sideSets[i] << "\"\n";

     TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
  }
  catch(STK_Interface::ElementBlockException & e) {
     std::stringstream ss;
     std::vector<std::string> elementBlocks; 
     mesh.getElementBlockNames(elementBlocks);

     // build an error message
     ss << e.what() << "\nChoose one of:\n";
     for(std::size_t i=0;i<elementBlocks.size();i++) 
        ss << "\"" << elementBlocks[i] << "\"\n";

     TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
  }
  catch(std::logic_error & e) {
     std::stringstream ss;
     ss << e.what() << "\nUnrecognized logic error.\n";

     TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
  }
  
  std::vector<stk::mesh::Entity> elements;
  std::map<std::pair<unsigned,unsigned>,std::vector<std::size_t> > local_cell_ids;
  if(!useCascade) {
    unsigned subcell_dim = needs.cellData.baseCellDimension()-1;
    std::vector<std::size_t> local_side_ids;
    getSideElements(mesh, eBlock,
  		      sideEntities,local_side_ids,elements);

    // build local cell_ids, mapped by local side id
    for(std::size_t elm=0;elm<elements.size();++elm) {
      stk::mesh::Entity element = elements[elm];
	
      local_cell_ids[std::make_pair(subcell_dim,local_side_ids[elm])].push_back(mesh.elementLocalId(element));
    }
  }
  else {
    std::vector<std::size_t> local_subcell_ids, subcell_dim;
    getSideElementCascade(mesh, eBlock,
  		          sideEntities,subcell_dim,local_subcell_ids,elements);

    // build local cell_ids, mapped by local side id
    for(std::size_t elm=0;elm<elements.size();++elm) {
      stk::mesh::Entity element = elements[elm];
	
      local_cell_ids[std::make_pair(subcell_dim[elm],local_subcell_ids[elm])].push_back(mesh.elementLocalId(element));
    }
  }

  // only build workset if there are elements to worry about
  // this may be processor dependent, so a defined boundary
  // condition may have not elements and thus no contribution
  // on this processor
  if(elements.size()!=0) {
    Teuchos::RCP<const shards::CellTopology> topo = mesh.getCellTopology(eBlock);

    // worksets to be returned
    Teuchos::RCP<std::vector<panzer::Workset> > worksets = Teuchos::rcp(new std::vector<panzer::Workset>);

    // loop over each side
    for(std::map<std::pair<unsigned,unsigned>,std::vector<std::size_t> >::const_iterator itr=local_cell_ids.begin();
        itr!=local_cell_ids.end();++itr) {
 
      if(itr->second.size()==0)
        continue;

      Kokkos::DynRankView<double,PHX::Device> nodes;
      mesh.getElementNodes(itr->second,eBlock,nodes);
  
      Teuchos::RCP<std::vector<panzer::Workset> > current
         = panzer::buildWorksets(needs, eBlock, itr->second, nodes);

      // correct worksets so the sides are correct
      for(std::size_t w=0;w<current->size();w++) {
        (*current)[w].subcell_dim = itr->first.first;
        (*current)[w].subcell_index = itr->first.second;
      }

      // append new worksets
      worksets->insert(worksets->end(),current->begin(),current->end());
    }

    return worksets;
  }
  
  // return Teuchos::null;
  return Teuchos::rcp(new std::vector<panzer::Workset>());
}

Teuchos::RCP<std::map<unsigned,panzer::Workset> >
buildBCWorksets(const panzer_stk::STK_Interface & mesh,
                const panzer::WorksetNeeds & needs_a,
                const std::string & blockid_a,
                const panzer::WorksetNeeds & needs_b,
                const std::string & blockid_b,
                const std::string & sideset)
{
  using namespace workset_utils;
  using Teuchos::RCP;

  std::vector<stk::mesh::Entity> sideEntities; // we will reduce a_ and b_ to this vector

  try {
     // grab local entities on this side
     // ...catch any failure...primarily wrong side set and element block info
     
     // we can't use getMySides because it only returns locally owned sides
     // this gurantees all the sides are extracted (element ownership is considered
     // we we call getSideElements below)

     stk::mesh::Part * sidePart = mesh.getSideset(sideset);
     TEUCHOS_TEST_FOR_EXCEPTION(sidePart==0,std::logic_error,
                        "Unknown side set \"" << sideset << "\"");

     stk::mesh::Selector side = *sidePart;
     // stk::mesh::Selector ownedBlock = metaData_->locally_owned_part() & side;

     // grab elements
     stk::mesh::get_selected_entities(side,mesh.getBulkData()->buckets(mesh.getSideRank()),sideEntities);
  } 
  catch(STK_Interface::ElementBlockException & e) {
     std::stringstream ss;
     std::vector<std::string> elementBlocks; 
     mesh.getElementBlockNames(elementBlocks);

     // build an error message
     ss << e.what() << "\nChoose one of:\n";
     for(std::size_t i=0;i<elementBlocks.size();i++) 
        ss << "\"" << elementBlocks[i] << "\"\n";

     TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
  }
  catch(std::logic_error & e) {
     std::stringstream ss;
     ss << e.what() << "\nUnrecognized logic error.\n";

     TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
  }

  std::vector<stk::mesh::Entity> elements_a, elements_b;
  std::vector<std::size_t> local_cell_ids_a, local_cell_ids_b;
  std::vector<std::size_t> local_side_ids_a, local_side_ids_b;

  // this enforces that "a" elements must be owned.
  getSideElements(mesh, blockid_a,blockid_b, sideEntities,
                        local_side_ids_a,elements_a, 
                        local_side_ids_b,elements_b);

  TEUCHOS_TEST_FOR_EXCEPTION(elements_a.size()!=elements_b.size(),std::logic_error,
                             "For a DG type boundary, the number of elements on the \"left\" and \"right\" is not the same.");

  // only build workset if there are elements to worry about
  // this may be processor dependent, so a defined boundary
  // condition may have not elements and thus no contribution
  // on this processor
  if(elements_a.size()==0)
    return Teuchos::rcp(new std::map<unsigned,panzer::Workset>);

  // loop over elements of this block (note the assures that element_a and element_b
  // are the same size, the ordering is the same because the order of sideEntities is
  // the same
  for(std::size_t elm=0;elm<elements_a.size();++elm) {
    stk::mesh::Entity element_a = elements_a[elm];
    stk::mesh::Entity element_b = elements_b[elm];
	
    local_cell_ids_a.push_back(mesh.elementLocalId(element_a));
    local_cell_ids_b.push_back(mesh.elementLocalId(element_b));
  }

  Kokkos::DynRankView<double,PHX::Device> node_coordinates_a, node_coordinates_b;
  mesh.getElementNodes(local_cell_ids_a,blockid_a,node_coordinates_a);
  mesh.getElementNodes(local_cell_ids_b,blockid_b,node_coordinates_b);

  // worksets to be returned
  return buildBCWorkset(needs_a,blockid_a, local_cell_ids_a, local_side_ids_a, node_coordinates_a,
                        needs_b,blockid_b, local_cell_ids_b, local_side_ids_b, node_coordinates_b);
}

Teuchos::RCP<std::map<unsigned,panzer::Workset> >
buildBCWorksets(const panzer_stk::STK_Interface & mesh,
                const panzer::WorksetNeeds & needs,
                const std::string & eblockID,
                const std::string & sidesetID)
{
  using namespace workset_utils;
  using Teuchos::RCP;

  std::vector<stk::mesh::Entity> sideEntities; 

  try {
     // grab local entities on this side
     // ...catch any failure...primarily wrong side set and element block info
     mesh.getMySides(sidesetID,eblockID,sideEntities);
  } 
  catch(STK_Interface::SidesetException & e) {
     std::stringstream ss;
     std::vector<std::string> sideSets; 
     mesh.getSidesetNames(sideSets);
 
     // build an error message
     ss << e.what() << "\nChoose one of:\n";
     for(std::size_t i=0;i<sideSets.size();i++) 
        ss << "\"" << sideSets[i] << "\"\n";

     TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
  }
  catch(STK_Interface::ElementBlockException & e) {
     std::stringstream ss;
     std::vector<std::string> elementBlocks; 
     mesh.getElementBlockNames(elementBlocks);

     // build an error message
     ss << e.what() << "\nChoose one of:\n";
     for(std::size_t i=0;i<elementBlocks.size();i++) 
        ss << "\"" << elementBlocks[i] << "\"\n";

     TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
  }
  catch(std::logic_error & e) {
     std::stringstream ss;
     ss << e.what() << "\nUnrecognized logic error.\n";

     TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
  }
  
  std::vector<stk::mesh::Entity> elements;
  std::vector<std::size_t> local_cell_ids;
  std::vector<std::size_t> local_side_ids;
  getSideElements(mesh, eblockID,
		      sideEntities,local_side_ids,elements);

  // loop over elements of this block
  for(std::size_t elm=0;elm<elements.size();++elm) {
	stk::mesh::Entity element = elements[elm];
	
	local_cell_ids.push_back(mesh.elementLocalId(element));
  }

  // only build workset if there are elements to worry about
  // this may be processor dependent, so a defined boundary
  // condition may have not elements and thus no contribution
  // on this processor
  if(elements.size()!=0) {
      Teuchos::RCP<const shards::CellTopology> topo 
         = mesh.getCellTopology(eblockID);

      Kokkos::DynRankView<double,PHX::Device> nodes;
      mesh.getElementNodes(local_cell_ids,eblockID,nodes);
  
      return panzer::buildBCWorkset(needs, eblockID, local_cell_ids, local_side_ids, nodes);
  }
  
  return Teuchos::null;
}

namespace workset_utils { 

void getSubcellElements(const panzer_stk::STK_Interface & mesh,
	 	        const std::string & blockId, 
		        const std::vector<stk::mesh::Entity> & entities,
		        std::vector<std::size_t> & localEntityIds, 
		        std::vector<stk::mesh::Entity> & elements)
{
  // for verifying that an element is in specified block
  stk::mesh::Part * blockPart = mesh.getElementBlockPart(blockId);
  stk::mesh::Part * ownedPart = mesh.getOwnedPart();
  stk::mesh::BulkData& bulkData = *mesh.getBulkData();

  // loop over each entitiy extracting elements and local entity ID that
  // are containted in specified block.
  std::vector<stk::mesh::Entity>::const_iterator entityItr;
  for(entityItr=entities.begin();entityItr!=entities.end();++entityItr) {
    stk::mesh::Entity entity = *entityItr;

    const size_t num_rels = bulkData.num_elements(entity);
    stk::mesh::Entity const* relations = bulkData.begin_elements(entity);
    stk::mesh::ConnectivityOrdinal const* ordinals = bulkData.begin_element_ordinals(entity);
    for(std::size_t e=0; e<num_rels; ++e) {
      stk::mesh::Entity element = relations[e];
      std::size_t entityId = ordinals[e];

      // is this element in requested block
      stk::mesh::Bucket const& bucket = bulkData.bucket(element);
      bool inBlock = bucket.member(*blockPart);
      bool onProc = bucket.member(*ownedPart);
      if(inBlock && onProc) {
        // add element and Side ID to output vectors
        elements.push_back(element);
        localEntityIds.push_back(entityId);
      }
    }
  }
}

void getUniversalSubcellElements(const panzer_stk::STK_Interface & mesh,
				 const std::string & blockId, 
				 const std::vector<stk::mesh::Entity> & entities,
				 std::vector<std::size_t> & localEntityIds, 
				 std::vector<stk::mesh::Entity> & elements,
                                 std::vector<std::size_t> & missingElementIndices)
{
  // for verifying that an element is in specified block
  stk::mesh::Part * blockPart = mesh.getElementBlockPart(blockId);
  stk::mesh::Part * universalPart = &mesh.getMetaData()->universal_part();
  stk::mesh::BulkData& bulkData = *mesh.getBulkData();

  // loop over each entitiy extracting elements and local entity ID that
  // are containted in specified block.
  std::size_t entityIndex =-1;
  std::vector<stk::mesh::Entity>::const_iterator entityItr;
  for(entityItr=entities.begin();entityItr!=entities.end();++entityItr) {
    stk::mesh::Entity entity = *entityItr;
    entityIndex += 1;

    const size_t num_rels = bulkData.num_elements(entity);
    stk::mesh::Entity const* element_rels = bulkData.begin_elements(entity);
    stk::mesh::ConnectivityOrdinal const* ordinals = bulkData.begin_element_ordinals(entity);
    for(std::size_t e=0; e<num_rels; ++e) {
      stk::mesh::Entity element = element_rels[e];
      std::size_t entityId = ordinals[e];

      // is this element in requested block
      stk::mesh::Bucket const& bucket = bulkData.bucket(element);
      bool inBlock = bucket.member(*blockPart);
      bool onProc = bucket.member(*universalPart);
      if(inBlock && onProc) {
        // add element and Side ID to output vectors
        elements.push_back(element);
        localEntityIds.push_back(entityId);
      } else if(!inBlock && (num_rels == 1)) {
        // add index of side whose neighbor element in blockPart does not belong 
        // to the current processor
        missingElementIndices.push_back(entityIndex);
      }
    }
  }
}

void getSideElementCascade(const panzer_stk::STK_Interface & mesh,
                           const std::string & blockId, 
                           const std::vector<stk::mesh::Entity> & sides,
                           std::vector<std::size_t> & localSubcellDim, 
                           std::vector<std::size_t> & localSubcellIds, 
                           std::vector<stk::mesh::Entity> & elements)
{
  // This is the alogrithm, for computing the side element
  // cascade. The requirements are that for a particular set of sides
  // we compute all elements and subcells where they touch the side. Note
  // that elements can be and will be repeated within this list.

  std::vector<std::vector<stk::mesh::Entity> > subcells;
  getSubcellEntities(mesh,sides,subcells);
  subcells.push_back(sides);

  // subcells now contains a unique list of faces, edges and nodes that
  // intersect with the sides

  for(std::size_t d=0;d<subcells.size();d++) {
    std::vector<std::size_t> subcellIds;
    std::vector<stk::mesh::Entity> subcellElements;

    // find elements connected to the subcells and their local subcell information
    getSubcellElements(mesh,blockId,subcells[d],subcellIds,subcellElements);

    // sanity check
    TEUCHOS_ASSERT(subcellIds.size()==subcellElements.size());

    // concatenate with found elements
    localSubcellDim.insert(localSubcellDim.end(),subcellElements.size(),d);
    localSubcellIds.insert(localSubcellIds.end(),subcellIds.begin(),subcellIds.end());
    elements.insert(elements.end(),subcellElements.begin(),subcellElements.end());
  }
}

void getSideElements(const panzer_stk::STK_Interface & mesh,
                     const std::string & blockId, 
                     const std::vector<stk::mesh::Entity> & sides,
                     std::vector<std::size_t> & localSideIds, 
                     std::vector<stk::mesh::Entity> & elements)
{
   getSubcellElements(mesh,blockId,sides,localSideIds,elements);
}

void getSideElements(const panzer_stk::STK_Interface & mesh,
                     const std::string & blockId_a, 
                     const std::string & blockId_b, 
                     const std::vector<stk::mesh::Entity> & sides,
                     std::vector<std::size_t> & localSideIds_a, 
                     std::vector<stk::mesh::Entity> & elements_a,
                     std::vector<std::size_t> & localSideIds_b, 
                     std::vector<stk::mesh::Entity> & elements_b)
{
  // for verifying that an element is in specified block
  stk::mesh::Part * blockPart_a = mesh.getElementBlockPart(blockId_a);
  stk::mesh::Part * blockPart_b = mesh.getElementBlockPart(blockId_b);
  stk::mesh::Part * ownedPart = mesh.getOwnedPart();
  stk::mesh::Part * universalPart = &mesh.getMetaData()->universal_part();
  stk::mesh::BulkData& bulkData = *mesh.getBulkData();

  // loop over each entitiy extracting elements and local entity ID that
  // are containted in specified block.
  std::vector<stk::mesh::Entity>::const_iterator sidesItr;
  for(sidesItr=sides.begin();sidesItr!=sides.end();++sidesItr) {
    stk::mesh::Entity side = *sidesItr;

     // these are used below the loop to insert into the appropriate vectors
    stk::mesh::Entity element_a = stk::mesh::Entity(), element_b = stk::mesh::Entity();
    std::size_t entityId_a=0, entityId_b=0;

    const size_t num_rels = bulkData.num_elements(side);
    stk::mesh::Entity const* element_rels = bulkData.begin_elements(side);
    stk::mesh::ConnectivityOrdinal const* ordinals = bulkData.begin_element_ordinals(side);
    for(std::size_t e=0; e<num_rels; ++e) {
      stk::mesh::Entity element = element_rels[e];
      std::size_t entityId = ordinals[e];

      // is this element in requested block
      stk::mesh::Bucket const& bucket = bulkData.bucket(element);
      bool inBlock_a = bucket.member(*blockPart_a);
      bool inBlock_b = bucket.member(*blockPart_b);
      bool onProc = bucket.member(*ownedPart);
      bool unProc = bucket.member(*universalPart);

      if(inBlock_a && onProc) {
        TEUCHOS_ASSERT(element_a==stk::mesh::Entity()); // sanity check
        element_a = element;
        entityId_a = entityId;
      }
      if(inBlock_b && unProc) {
        TEUCHOS_ASSERT(element_b==stk::mesh::Entity()); // sanity check
        element_b = element;
        entityId_b = entityId;
      }
    }

    if(element_a!=stk::mesh::Entity() && element_b!=stk::mesh::Entity()) {      // add element and Side ID to output vectors
      elements_a.push_back(element_a);
      localSideIds_a.push_back(entityId_a);

      // add element and Side ID to output vectors
      elements_b.push_back(element_b);
      localSideIds_b.push_back(entityId_b);
    }
  }
}

void getNodeElements(const panzer_stk::STK_Interface & mesh,
                     const std::string & blockId, 
                     const std::vector<stk::mesh::Entity> & nodes,
                     std::vector<std::size_t> & localNodeIds, 
                     std::vector<stk::mesh::Entity> & elements)
{
   getSubcellElements(mesh,blockId,nodes,localNodeIds,elements);
}

void getSubcellEntities(const panzer_stk::STK_Interface & mesh,
		        const std::vector<stk::mesh::Entity> & entities,
	 	        std::vector<std::vector<stk::mesh::Entity> > & subcells)
{
  // exit if there is no work to do
  if(entities.size()==0) {
    subcells.clear();
    return;
  }

  stk::mesh::BulkData& bulkData = *mesh.getBulkData();
  stk::mesh::EntityRank master_rank = bulkData.entity_rank(entities[0]);

  std::vector<std::set<stk::mesh::Entity> > subcells_set(master_rank);

  // loop over each entitiy extracting elements and local entity ID that
  // are containted in specified block.
  std::vector<stk::mesh::Entity>::const_iterator entityItr;
  for(entityItr=entities.begin();entityItr!=entities.end();++entityItr) {
    stk::mesh::Entity entity = *entityItr;

    // sanity check, enforcing that there is only one rank
    TEUCHOS_ASSERT(bulkData.entity_rank(entity)==master_rank);

    for(int i=0; i<static_cast<int>(master_rank); i++) {
      stk::mesh::EntityRank const to_rank = static_cast<stk::mesh::EntityRank>(i);
      const size_t num_rels = bulkData.num_connectivity(entity, to_rank);
      stk::mesh::Entity const* relations = bulkData.begin(entity, to_rank);

      // for each relation insert the appropriate entity (into the set
      // which gurantees uniqueness
      for(std::size_t e=0; e<num_rels; ++e) {
        stk::mesh::Entity subcell = relations[e];

        subcells_set[i].insert(subcell);
      }
    }
  }

  // copy unique entities into vector
  subcells.clear();
  subcells.resize(subcells_set.size());
  for(std::size_t i=0;i<subcells_set.size();i++)
    subcells[i].assign(subcells_set[i].begin(),subcells_set[i].end());
}

}

}
