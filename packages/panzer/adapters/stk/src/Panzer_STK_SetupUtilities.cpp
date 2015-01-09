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

#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Teuchos_Assert.hpp"

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace panzer_stk_classic { 
Teuchos::RCP<std::vector<panzer::Workset> >  
buildWorksets(const panzer_stk_classic::STK_Interface & mesh,
              const panzer::PhysicsBlock & pb)
{
  using namespace workset_utils;

  std::vector<std::string> element_blocks;

  std::vector<std::size_t> local_cell_ids;
  Intrepid::FieldContainer<double> cell_vertex_coordinates;

  getIdsAndVertices(mesh, pb.elementBlockID(), local_cell_ids, cell_vertex_coordinates);

  // only build workset if there are elements to worry about
  // this may be processor dependent, so an element block
  // may not have elements and thus no contribution
  // on this processor
  return panzer::buildWorksets(pb, local_cell_ids, cell_vertex_coordinates);
}

Teuchos::RCP<std::vector<panzer::Workset> >  
buildWorksets(const panzer_stk_classic::STK_Interface & mesh,
              const panzer::PhysicsBlock & pb,
              const std::string & sideset,
              bool useCascade)
{
  using namespace workset_utils;
  using Teuchos::RCP;

  std::vector<stk_classic::mesh::Entity*> sideEntities; 

  try {
     // grab local entities on this side
     // ...catch any failure...primarily wrong side set and element block info
     mesh.getMySides(sideset,pb.elementBlockID(),sideEntities);
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
  
  std::vector<stk_classic::mesh::Entity*> elements;
  std::map<std::pair<unsigned,unsigned>,std::vector<std::size_t> > local_cell_ids;
  if(!useCascade) {
    unsigned subcell_dim = pb.cellData().baseCellDimension()-1;
    std::vector<std::size_t> local_side_ids;
    getSideElements(mesh, pb.elementBlockID(),
  		      sideEntities,local_side_ids,elements);

    // build local cell_ids, mapped by local side id
    for(std::size_t elm=0;elm<elements.size();++elm) {
      stk_classic::mesh::Entity * element = elements[elm];
	
      local_cell_ids[std::make_pair(subcell_dim,local_side_ids[elm])].push_back(mesh.elementLocalId(element));
    }
  }
  else {
    std::vector<std::size_t> local_subcell_ids, subcell_dim;
    getSideElementCascade(mesh, pb.elementBlockID(),
  		          sideEntities,subcell_dim,local_subcell_ids,elements);

    // build local cell_ids, mapped by local side id
    for(std::size_t elm=0;elm<elements.size();++elm) {
      stk_classic::mesh::Entity * element = elements[elm];
	
      local_cell_ids[std::make_pair(subcell_dim[elm],local_subcell_ids[elm])].push_back(mesh.elementLocalId(element));
    }
  }

  // only build workset if there are elements to worry about
  // this may be processor dependent, so a defined boundary
  // condition may have not elements and thus no contribution
  // on this processor
  if(elements.size()!=0) {
    Teuchos::RCP<const shards::CellTopology> topo = mesh.getCellTopology(pb.elementBlockID());

    // worksets to be returned
    Teuchos::RCP<std::vector<panzer::Workset> > worksets = Teuchos::rcp(new std::vector<panzer::Workset>);

    // loop over each side
    for(std::map<std::pair<unsigned,unsigned>,std::vector<std::size_t> >::const_iterator itr=local_cell_ids.begin();
        itr!=local_cell_ids.end();++itr) {
 
      if(itr->second.size()==0)
        continue;

      Intrepid::FieldContainer<double> vertices;
      mesh.getElementVertices(itr->second,pb.elementBlockID(),vertices);
  
      Teuchos::RCP<std::vector<panzer::Workset> > current
         = panzer::buildWorksets(pb, itr->second, vertices);

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

Teuchos::RCP<std::vector<panzer::Workset> >  
buildWorksets(const panzer_stk_classic::STK_Interface & mesh,
              const panzer::PhysicsBlock & pb_a,
              const panzer::PhysicsBlock & pb_b,
              const std::string & sideset)
{
  using namespace workset_utils;
  using Teuchos::RCP;

  std::vector<stk_classic::mesh::Entity*> sideEntities; // we will reduce a_ and b_ to this vector

  try {
     // grab local entities on this side
     // ...catch any failure...primarily wrong side set and element block info
     
     // we can't use getMySides because it only returns locally owned sides
     // this gurantees all the sides are extracted (element ownership is considered
     // we we call getSideElements below)

     stk_classic::mesh::Part * sidePart = mesh.getSideset(sideset);
     TEUCHOS_TEST_FOR_EXCEPTION(sidePart==0,std::logic_error,
                        "Unknown side set \"" << sideset << "\"");

     stk_classic::mesh::Selector side = *sidePart;
     // stk_classic::mesh::Selector ownedBlock = metaData_->locally_owned_part() & side;

     // grab elements
     stk_classic::mesh::get_selected_entities(side,mesh.getBulkData()->buckets(mesh.getSideRank()),sideEntities);
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

  std::vector<stk_classic::mesh::Entity*> elements_a, elements_b;
  std::vector<std::size_t> local_cell_ids_a, local_cell_ids_b;
  std::vector<std::size_t> local_side_ids_a, local_side_ids_b;

  // this enforces that "a" elements must be owned.
  getSideElements(mesh, pb_a.elementBlockID(), pb_b.elementBlockID(), sideEntities,
                        local_side_ids_a,elements_a, 
                        local_side_ids_b,elements_b);

  TEUCHOS_TEST_FOR_EXCEPTION(elements_a.size()!=elements_b.size(),std::logic_error,
                             "For a DG type boundary, the number of elements on the \"left\" and \"right\" is not the same.");

  // only build workset if there are elements to worry about
  // this may be processor dependent, so a defined boundary
  // condition may have not elements and thus no contribution
  // on this processor
  if(elements_a.size()==0)
    return Teuchos::rcp(new std::vector<panzer::Workset>());

  // loop over elements of this block (note the assures that element_a and element_b
  // are the same size, the ordering is the same because the order of sideEntities is
  // the same
  for(std::size_t elm=0;elm<elements_a.size();++elm) {
    stk_classic::mesh::Entity * element_a = elements_a[elm];
    stk_classic::mesh::Entity * element_b = elements_b[elm];
	
    local_cell_ids_a.push_back(mesh.elementLocalId(element_a));
    local_cell_ids_b.push_back(mesh.elementLocalId(element_b));
  }

  Intrepid::FieldContainer<double> vertex_coordinates_a, vertex_coordinates_b;
  mesh.getElementVertices(local_cell_ids_a,pb_a.elementBlockID(),vertex_coordinates_a);
  mesh.getElementVertices(local_cell_ids_b,pb_b.elementBlockID(),vertex_coordinates_b);

  // worksets to be returned
  Teuchos::RCP<std::vector<panzer::Workset> > worksets;
  worksets = buildEdgeWorksets(pb_a, local_cell_ids_a, local_side_ids_a, vertex_coordinates_a,
                               pb_b, local_cell_ids_b, local_side_ids_b, vertex_coordinates_b);

  return worksets;
}

Teuchos::RCP<std::map<unsigned,panzer::Workset> >
buildBCWorksets(const panzer_stk_classic::STK_Interface & mesh,
                const panzer::PhysicsBlock & pb,
                const std::string & sidesetID)
{
  using namespace workset_utils;
  using Teuchos::RCP;

  std::vector<stk_classic::mesh::Entity*> sideEntities; 

  try {
     // grab local entities on this side
     // ...catch any failure...primarily wrong side set and element block info
     mesh.getMySides(sidesetID,pb.elementBlockID(),sideEntities);
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
  
  std::vector<stk_classic::mesh::Entity*> elements;
  std::vector<std::size_t> local_cell_ids;
  std::vector<std::size_t> local_side_ids;
  getSideElements(mesh, pb.elementBlockID(),
		      sideEntities,local_side_ids,elements);

  // loop over elements of this block
  for(std::size_t elm=0;elm<elements.size();++elm) {
	stk_classic::mesh::Entity * element = elements[elm];
	
	local_cell_ids.push_back(mesh.elementLocalId(element));
  }

  // only build workset if there are elements to worry about
  // this may be processor dependent, so a defined boundary
  // condition may have not elements and thus no contribution
  // on this processor
  if(elements.size()!=0) {
      Teuchos::RCP<const shards::CellTopology> topo 
         = mesh.getCellTopology(pb.elementBlockID());

      Intrepid::FieldContainer<double> vertices;
      mesh.getElementVertices(local_cell_ids,pb.elementBlockID(),vertices);
  
      return panzer::buildBCWorkset(pb, local_cell_ids, local_side_ids, vertices);
  }
  
  return Teuchos::null;
}

namespace workset_utils { 

void getSubcellElements(const panzer_stk_classic::STK_Interface & mesh,
	 	        const std::string & blockId, 
		        const std::vector<stk_classic::mesh::Entity*> & entities,
		        std::vector<std::size_t> & localEntityIds, 
		        std::vector<stk_classic::mesh::Entity*> & elements)
{
  // for verifying that an element is in specified block
  stk_classic::mesh::Part * blockPart = mesh.getElementBlockPart(blockId);
  stk_classic::mesh::Part * ownedPart = mesh.getOwnedPart();
  stk_classic::mesh::EntityRank elementRank = mesh.getElementRank();
  
  // loop over each entitiy extracting elements and local entity ID that
  // are containted in specified block.
  std::vector<stk_classic::mesh::Entity*>::const_iterator entityItr;
  for(entityItr=entities.begin();entityItr!=entities.end();++entityItr) {
    stk_classic::mesh::Entity * entity = *entityItr;
    
    stk_classic::mesh::PairIterRelation relations = entity->relations(elementRank);

    for(std::size_t e=0;e<relations.size();++e) {
      stk_classic::mesh::Entity * element = relations[e].entity();
      std::size_t entityId = relations[e].identifier();
	
      // is this element in requested block
      bool inBlock = element->bucket().member(*blockPart);
      bool onProc = element->bucket().member(*ownedPart);
      if(inBlock && onProc) {
        // add element and Side ID to output vectors
        elements.push_back(element);
        localEntityIds.push_back(entityId);
      }
    }
  }
}

void getUniversalSubcellElements(const panzer_stk_classic::STK_Interface & mesh,
				 const std::string & blockId, 
				 const std::vector<stk_classic::mesh::Entity*> & entities,
				 std::vector<std::size_t> & localEntityIds, 
				 std::vector<stk_classic::mesh::Entity*> & elements)
{
  // for verifying that an element is in specified block
  stk_classic::mesh::Part * blockPart = mesh.getElementBlockPart(blockId);
  stk_classic::mesh::Part * universalPart = &mesh.getMetaData()->universal_part();
  stk_classic::mesh::EntityRank elementRank = mesh.getElementRank();
  
  // loop over each entitiy extracting elements and local entity ID that
  // are containted in specified block.
  std::vector<stk_classic::mesh::Entity*>::const_iterator entityItr;
  for(entityItr=entities.begin();entityItr!=entities.end();++entityItr) {
    stk_classic::mesh::Entity * entity = *entityItr;
    
    stk_classic::mesh::PairIterRelation relations = entity->relations(elementRank);

    for(std::size_t e=0;e<relations.size();++e) {
      stk_classic::mesh::Entity * element = relations[e].entity();
      std::size_t entityId = relations[e].identifier();
	
      // is this element in requested block
      bool inBlock = element->bucket().member(*blockPart);
      bool onProc = element->bucket().member(*universalPart);
      if(inBlock && onProc) {
        // add element and Side ID to output vectors
        elements.push_back(element);
        localEntityIds.push_back(entityId);
      }
    }
  }
}

void getSideElementCascade(const panzer_stk_classic::STK_Interface & mesh,
                           const std::string & blockId, 
                           const std::vector<stk_classic::mesh::Entity*> & sides,
                           std::vector<std::size_t> & localSubcellDim, 
                           std::vector<std::size_t> & localSubcellIds, 
                           std::vector<stk_classic::mesh::Entity*> & elements)
{
  // This is the alogrithm, for computing the side element
  // cascade. The requirements are that for a particular set of sides
  // we compute all elements and subcells where they touch the side. Note
  // that elements can be and will be repeated within this list.

  std::vector<std::vector<stk_classic::mesh::Entity*> > subcells;
  getSubcellEntities(mesh,sides,subcells);
  subcells.push_back(sides);

  // subcells now contains a unique list of faces, edges and nodes that
  // intersect with the sides

  for(std::size_t d=0;d<subcells.size();d++) {
    std::vector<std::size_t> subcellIds;
    std::vector<stk_classic::mesh::Entity*> subcellElements;

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

void getSideElements(const panzer_stk_classic::STK_Interface & mesh,
                     const std::string & blockId, 
                     const std::vector<stk_classic::mesh::Entity*> & sides,
                     std::vector<std::size_t> & localSideIds, 
                     std::vector<stk_classic::mesh::Entity*> & elements)
{
   getSubcellElements(mesh,blockId,sides,localSideIds,elements);
}

void getSideElements(const panzer_stk_classic::STK_Interface & mesh,
                     const std::string & blockId_a, 
                     const std::string & blockId_b, 
                     const std::vector<stk_classic::mesh::Entity*> & sides,
                     std::vector<std::size_t> & localSideIds_a, 
                     std::vector<stk_classic::mesh::Entity*> & elements_a,
                     std::vector<std::size_t> & localSideIds_b, 
                     std::vector<stk_classic::mesh::Entity*> & elements_b)
{
  // for verifying that an element is in specified block
  stk_classic::mesh::Part * blockPart_a = mesh.getElementBlockPart(blockId_a);
  stk_classic::mesh::Part * blockPart_b = mesh.getElementBlockPart(blockId_b);
  stk_classic::mesh::Part * ownedPart = mesh.getOwnedPart();
  stk_classic::mesh::Part * universalPart = &mesh.getMetaData()->universal_part();
  stk_classic::mesh::EntityRank elementRank = mesh.getElementRank();
  
  // loop over each entitiy extracting elements and local entity ID that
  // are containted in specified block.
  std::vector<stk_classic::mesh::Entity*>::const_iterator sidesItr;
  for(sidesItr=sides.begin();sidesItr!=sides.end();++sidesItr) {
    stk_classic::mesh::Entity * side = *sidesItr;
    
     // these are used below the loop to insert into the appropriate vectors
    stk_classic::mesh::Entity * element_a=0,* element_b=0;
    std::size_t entityId_a=0, entityId_b=0;

    stk_classic::mesh::PairIterRelation relations = side->relations(elementRank);
    for(std::size_t e=0;e<relations.size();++e) {
      stk_classic::mesh::Entity * element = relations[e].entity();
      std::size_t entityId = relations[e].identifier();
	
      // is this element in requested block
      bool inBlock_a = element->bucket().member(*blockPart_a);
      bool inBlock_b = element->bucket().member(*blockPart_b);
      bool onProc = element->bucket().member(*ownedPart);
      bool unProc = element->bucket().member(*universalPart);

      if(inBlock_a && onProc) {
        TEUCHOS_ASSERT(element_a==0); // sanity check
        element_a = element;
        entityId_a = entityId;
      }
      if(inBlock_b && unProc) {
        TEUCHOS_ASSERT(element_b==0); // sanity check
        element_b = element;
        entityId_b = entityId;
      }
    }

    if(element_a!=0 && element_b!=0) {
      // add element and Side ID to output vectors
      elements_a.push_back(element_a);
      localSideIds_a.push_back(entityId_a);

      // add element and Side ID to output vectors
      elements_b.push_back(element_b);
      localSideIds_b.push_back(entityId_b);
    }
  }
}

void getNodeElements(const panzer_stk_classic::STK_Interface & mesh,
                     const std::string & blockId, 
                     const std::vector<stk_classic::mesh::Entity*> & nodes,
                     std::vector<std::size_t> & localNodeIds, 
                     std::vector<stk_classic::mesh::Entity*> & elements)
{
   getSubcellElements(mesh,blockId,nodes,localNodeIds,elements);
}

void getSubcellEntities(const panzer_stk_classic::STK_Interface & mesh,
		        const std::vector<stk_classic::mesh::Entity*> & entities,
	 	        std::vector<std::vector<stk_classic::mesh::Entity*> > & subcells)
{
  // exit if there is no work to do
  if(entities.size()==0) {
    subcells.clear();
    return;
  }
 
  int maxRankIndex = mesh.getDimension()-1;
  stk_classic::mesh::EntityRank master_rank = entities[0]->entity_rank();
  std::vector<stk_classic::mesh::EntityRank> ranks(mesh.getDimension()+1);

  // build rank array, and compute maximum rank index (within rank array)
  // for these entities with "master_rank"
  switch(mesh.getDimension()) {
  case 3:
    ranks[2] = mesh.getFaceRank();
    maxRankIndex = (master_rank==mesh.getFaceRank() ? 1 : maxRankIndex);
  case 2:
    ranks[1] = mesh.getEdgeRank();
    maxRankIndex = (master_rank==mesh.getEdgeRank() ? 0 : maxRankIndex);
  case 1:
    ranks[0] = mesh.getNodeRank();
    maxRankIndex = (master_rank==mesh.getNodeRank() ? -1 : maxRankIndex);
    break;
  default:
    TEUCHOS_ASSERT(false);
    break;
  };
  ranks[mesh.getDimension()] = mesh.getElementRank();

  // make sure the rank index is ok
  TEUCHOS_ASSERT(maxRankIndex>-1);

  std::vector<std::set<stk_classic::mesh::Entity*> > subcells_set(maxRankIndex+1);

  // loop over each entitiy extracting elements and local entity ID that
  // are containted in specified block.
  std::vector<stk_classic::mesh::Entity*>::const_iterator entityItr;
  for(entityItr=entities.begin();entityItr!=entities.end();++entityItr) {
    stk_classic::mesh::Entity * entity = *entityItr;

    // sanity check, enforcing that there is only one rank
    TEUCHOS_ASSERT(entity->entity_rank()==master_rank); 
    
    for(int i=0;i<=maxRankIndex;i++) {
      stk_classic::mesh::PairIterRelation relations = entity->relations(ranks[i]);

      // for each relation insert the appropriate entity (into the set
      // which gurantees uniqueness
      for(std::size_t e=0;e<relations.size();++e) {
        stk_classic::mesh::Entity * subcell = relations[e].entity();

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
