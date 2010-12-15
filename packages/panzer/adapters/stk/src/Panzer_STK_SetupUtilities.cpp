#include "Panzer_STK_SetupUtilities.hpp"

#include "Panzer_Workset_Builder.hpp"

namespace panzer_stk { 

std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
buildWorksets(const panzer_stk::STK_Interface & mesh,
              const panzer::InputPhysicsBlock & ipb, 
              const std::size_t workset_size)
{
  using namespace workset_utils;

  std::vector<std::string> element_blocks;
  mesh.getElementBlockNames(element_blocks);
  int base_cell_dimension = mesh.getDimension();

  std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > worksets;

  for (std::vector<std::string>::size_type i=0; i < element_blocks.size(); 
	 ++i) {

    std::vector<std::size_t> local_cell_ids;
    Intrepid::FieldContainer<double> cell_vertex_coordinates;

    getIdsAndVertices(mesh, element_blocks[i], local_cell_ids, 
				cell_vertex_coordinates);

    // only build workset if there are elements to worry about
    // this may be processor dependent, so an element block
    // may not have elements and thus no contribution
    // on this processor
    if(local_cell_ids.size()!=0) {
        worksets.insert(std::make_pair(element_blocks[i],panzer::buildWorksets(element_blocks[i],
      							                       local_cell_ids,
							                       cell_vertex_coordinates,
							                       ipb,
							                       workset_size,
							                       base_cell_dimension)));
    }
  }
  
  return worksets;
}

const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC>
buildBCWorksets(const panzer_stk::STK_Interface & mesh,
                const panzer::InputPhysicsBlock & ipb,
                const std::vector<panzer::BC> & bcs) 
{
  using namespace workset_utils;
  using Teuchos::RCP;

  int base_cell_dimension = mesh.getDimension();

  std::vector<std::string> sideSets; 
  std::vector<std::string> elementBlocks; 
  mesh.getSidesetNames(sideSets);
  mesh.getElementBlockNames(elementBlocks);
  
  std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets;
  
  for (std::vector<panzer::BC>::const_iterator bc = bcs.begin();
	 bc != bcs.end(); ++bc) {
    
    std::vector<stk::mesh::Entity*> sideEntities; 
    mesh.getMySides(bc->sidesetID(),bc->elementBlockID(),sideEntities);
    
    std::vector<stk::mesh::Entity*> elements;
    std::vector<std::size_t> local_cell_ids;
    std::vector<std::size_t> local_side_ids;
    getSideElements(mesh, bc->elementBlockID(),
		      sideEntities,local_side_ids,elements);

    // loop over elements of this block
    for(std::size_t elm=0;elm<elements.size();++elm) {
	std::vector<stk::mesh::EntityId> nodes;
	stk::mesh::Entity * element = elements[elm];
	
	local_cell_ids.push_back(mesh.elementLocalId(element));
    }

    // only build workset if there are elements to worry about
    // this may be processor dependent, so a defined boundary
    // condition may have not elements and thus no contribution
    // on this processor
    if(elements.size()!=0) {
        Intrepid::FieldContainer<double> vertices;
        mesh.getElementVertices(local_cell_ids,vertices);
    
        Teuchos::RCP<std::map<unsigned,panzer::Workset> > workset = 
              panzer::buildBCWorkset(*bc, local_cell_ids, local_side_ids,
  	                             vertices, ipb, base_cell_dimension);
   
       bc_worksets[*bc] = workset;
    }
  }
  
  return bc_worksets;
}

namespace workset_utils { 

void getSideElements(const panzer_stk::STK_Interface & mesh,
                     const std::string & blockId, 
                     const std::vector<stk::mesh::Entity*> & sides,
                     std::vector<std::size_t> & localSideIds, 
                     std::vector<stk::mesh::Entity*> & elements)
{
  // for verifying that an element is in specified block
  stk::mesh::Part * blockPart = mesh.getElementBlockPart(blockId);
  stk::mesh::EntityRank elementRank = mesh.getElementRank();
  
  // loop over each side extracting elements and local side ID that
  // are containted in specified block.
  std::vector<stk::mesh::Entity*>::const_iterator sideItr;
  for(sideItr=sides.begin();sideItr!=sides.end();++sideItr) {
    stk::mesh::Entity * side = *sideItr;
    
    stk::mesh::PairIterRelation relations = side->relations(elementRank);

    for(std::size_t e=0;e<relations.size();++e) {
      stk::mesh::Entity * element = relations[e].entity();
      std::size_t sideId = relations[e].identifier();
	
      // is this element in requested block
      bool inBlock = element->bucket().member(*blockPart);
      if(inBlock) {
        // add element and Side ID to output vectors
        elements.push_back(element);
        localSideIds.push_back(sideId);
      }
    }
  }
}

}

}
