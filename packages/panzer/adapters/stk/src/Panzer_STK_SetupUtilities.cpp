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

namespace panzer_stk { 

std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
buildWorksets(const panzer_stk::STK_Interface & mesh,
              const std::map<std::string,panzer::InputPhysicsBlock> & eb_to_ipb, 
              const std::size_t workset_size)
{
  using namespace workset_utils;

  std::vector<std::string> element_blocks;
  mesh.getElementBlockNames(element_blocks);

  std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > worksets;

  for (std::vector<std::string>::size_type i=0; i < element_blocks.size(); 
	 ++i) {

    std::map<std::string,panzer::InputPhysicsBlock>::const_iterator ipb_iterator = 
      eb_to_ipb.find(element_blocks[i]);
    
    // on error print ot all available worksets
    if(ipb_iterator==eb_to_ipb.end()) {
       std::stringstream ss;

       ss << "buildWorksets: Could not find input physics block corresponding to element block"
          << " \"" << element_blocks[i] << "\"\n\n Choose one of:\n";

       std::vector<std::string>::const_iterator str_iter;
       for(str_iter=element_blocks.begin();str_iter!=element_blocks.end();++str_iter)
          ss << "   \"" << *str_iter << "\"\n"; 

       TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true, std::logic_error,ss.str());

       // should never get here!
    }

    const panzer::InputPhysicsBlock& ipb = ipb_iterator->second;

    worksets.insert(std::make_pair(element_blocks[i],
                                   panzer_stk::buildWorksets(mesh,element_blocks[i],ipb,workset_size)));
  }
  
  return worksets;
}

Teuchos::RCP<std::vector<panzer::Workset> >  
buildWorksets(const panzer_stk::STK_Interface & mesh,
              const std::string & eBlock,
              const panzer::InputPhysicsBlock & ipb, 
              const std::size_t workset_size)
{
  using namespace workset_utils;

  std::vector<std::string> element_blocks;
  int base_cell_dimension = mesh.getDimension();

  std::vector<std::size_t> local_cell_ids;
  Intrepid::FieldContainer<double> cell_vertex_coordinates;

  getIdsAndVertices(mesh, eBlock, local_cell_ids, cell_vertex_coordinates);

  Teuchos::RCP<const shards::CellTopology> topo =  mesh.getCellTopology(eBlock);

  // only build workset if there are elements to worry about
  // this may be processor dependent, so an element block
  // may not have elements and thus no contribution
  // on this processor
  return panzer::buildWorksets(eBlock,topo, local_cell_ids, cell_vertex_coordinates,
                               ipb, workset_size, base_cell_dimension);
}

Teuchos::RCP<std::vector<panzer::Workset> >  
buildWorksets(const panzer_stk::STK_Interface & mesh,
              const panzer::PhysicsBlock & pb, 
              const std::size_t workset_size)
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
  return panzer::buildWorksets(pb, local_cell_ids, cell_vertex_coordinates,
                               workset_size);
}

const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC>
buildBCWorksets(const panzer_stk::STK_Interface & mesh,
                const std::map<std::string,panzer::InputPhysicsBlock> & eb_to_ipb,
                const std::vector<panzer::BC> & bcs) 
{
  using namespace workset_utils;
  using Teuchos::RCP;

  // int base_cell_dimension = mesh.getDimension();

  std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets;
  
  for (std::vector<panzer::BC>::const_iterator bc = bcs.begin();
	 bc != bcs.end(); ++bc) {

    std::map<std::string,panzer::InputPhysicsBlock>::const_iterator ipb_iterator = 
      eb_to_ipb.find(bc->elementBlockID());
    
    TEUCHOS_TEST_FOR_EXCEPTION(ipb_iterator == eb_to_ipb.end(), std::logic_error,
		       "Could not find input physics block corresponding to region");

    const panzer::InputPhysicsBlock& ipb = ipb_iterator->second;

    Teuchos::RCP<std::map<unsigned,panzer::Workset> > wMap = buildBCWorksets(mesh,ipb,*bc);
    if(wMap!=Teuchos::null)
       bc_worksets[*bc] = wMap;
  }
  
  return bc_worksets;
}

Teuchos::RCP<std::map<unsigned,panzer::Workset> >
buildBCWorksets(const panzer_stk::STK_Interface & mesh,
                const panzer::PhysicsBlock & pb,
                const panzer::BC & bc)
{
   return buildBCWorksets(mesh,pb.getInputPhysicsBlock(),bc);
}

Teuchos::RCP<std::map<unsigned,panzer::Workset> >
buildBCWorksets(const panzer_stk::STK_Interface & mesh,
                const panzer::InputPhysicsBlock & ipb,
                const panzer::BC & bc)
{
  using namespace workset_utils;
  using Teuchos::RCP;

  int base_cell_dimension = mesh.getDimension();

  std::vector<stk::mesh::Entity*> sideEntities; 

  try {
     // grab local entities on this side
     // ...catch any failure...primarily wrong side set and element block info
     mesh.getMySides(bc.sidesetID(),bc.elementBlockID(),sideEntities);
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
  
  std::vector<stk::mesh::Entity*> elements;
  std::vector<std::size_t> local_cell_ids;
  std::vector<std::size_t> local_side_ids;
  getSideElements(mesh, bc.elementBlockID(),
		      sideEntities,local_side_ids,elements);

  // loop over elements of this block
  for(std::size_t elm=0;elm<elements.size();++elm) {
	stk::mesh::Entity * element = elements[elm];
	
	local_cell_ids.push_back(mesh.elementLocalId(element));
  }

  // only build workset if there are elements to worry about
  // this may be processor dependent, so a defined boundary
  // condition may have not elements and thus no contribution
  // on this processor
  if(elements.size()!=0) {
      Teuchos::RCP<const shards::CellTopology> topo 
         = mesh.getCellTopology(bc.elementBlockID());

      Intrepid::FieldContainer<double> vertices;
      mesh.getElementVertices(local_cell_ids,vertices);
  
      return panzer::buildBCWorkset(bc,topo, local_cell_ids, local_side_ids,
	                            vertices, ipb, base_cell_dimension);
  }
  
  return Teuchos::null;
}

namespace workset_utils { 

void getSubcellElements(const panzer_stk::STK_Interface & mesh,
	 	        const std::string & blockId, 
		        const std::vector<stk::mesh::Entity*> & entities,
		        std::vector<std::size_t> & localEntityIds, 
		        std::vector<stk::mesh::Entity*> & elements)
{
  // for verifying that an element is in specified block
  stk::mesh::Part * blockPart = mesh.getElementBlockPart(blockId);
  stk::mesh::EntityRank elementRank = mesh.getElementRank();
  
  // loop over each entitiy extracting elements and local entity ID that
  // are containted in specified block.
  std::vector<stk::mesh::Entity*>::const_iterator entityItr;
  for(entityItr=entities.begin();entityItr!=entities.end();++entityItr) {
    stk::mesh::Entity * entity = *entityItr;
    
    stk::mesh::PairIterRelation relations = entity->relations(elementRank);

    for(std::size_t e=0;e<relations.size();++e) {
      stk::mesh::Entity * element = relations[e].entity();
      std::size_t entityId = relations[e].identifier();
	
      // is this element in requested block
      bool inBlock = element->bucket().member(*blockPart);
      if(inBlock) {
        // add element and Side ID to output vectors
        elements.push_back(element);
        localEntityIds.push_back(entityId);
      }
    }
  }
}

void getSideElements(const panzer_stk::STK_Interface & mesh,
                     const std::string & blockId, 
                     const std::vector<stk::mesh::Entity*> & sides,
                     std::vector<std::size_t> & localSideIds, 
                     std::vector<stk::mesh::Entity*> & elements)
{
   getSubcellElements(mesh,blockId,sides,localSideIds,elements);
}

void getNodeElements(const panzer_stk::STK_Interface & mesh,
                     const std::string & blockId, 
                     const std::vector<stk::mesh::Entity*> & nodes,
                     std::vector<std::size_t> & localNodeIds, 
                     std::vector<stk::mesh::Entity*> & elements)
{
   getSubcellElements(mesh,blockId,nodes,localNodeIds,elements);
}

}

}
