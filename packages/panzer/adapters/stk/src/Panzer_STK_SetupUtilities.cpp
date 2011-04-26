#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Teuchos_TestForException.hpp"

namespace panzer_stk { 

std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > > 
buildWorksets(const panzer_stk::STK_Interface & mesh,
              const std::map<std::string,panzer::InputPhysicsBlock> & eb_to_ipb, 
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

       TEST_FOR_EXCEPTION_PURE_MSG(true, std::logic_error,ss.str());

       // should never get here!
    }

    const panzer::InputPhysicsBlock& ipb = ipb_iterator->second;

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
    else {
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
                const std::map<std::string,panzer::InputPhysicsBlock> & eb_to_ipb,
                const std::vector<panzer::BC> & bcs) 
{
  using namespace workset_utils;
  using Teuchos::RCP;

  int base_cell_dimension = mesh.getDimension();

  std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC> bc_worksets;
  
  for (std::vector<panzer::BC>::const_iterator bc = bcs.begin();
	 bc != bcs.end(); ++bc) {
    
    std::vector<stk::mesh::Entity*> sideEntities; 

    try {
       // grab local entities on this side
       // ...catch any failure...primarily wrong side set and element block info
       mesh.getMySides(bc->sidesetID(),bc->elementBlockID(),sideEntities);
    } 
    catch(STK_Interface::SidesetException & e) {
       std::stringstream ss;
       std::vector<std::string> sideSets; 
       mesh.getSidesetNames(sideSets);
 
       // build an error message
       ss << e.what() << "\nChoose one of:\n";
       for(std::size_t i=0;i<sideSets.size();i++) 
          ss << "\"" << sideSets[i] << "\"\n";

       TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
    }
    catch(STK_Interface::ElementBlockException & e) {
       std::stringstream ss;
       std::vector<std::string> elementBlocks; 
       mesh.getElementBlockNames(elementBlocks);
  
       // build an error message
       ss << e.what() << "\nChoose one of:\n";
       for(std::size_t i=0;i<elementBlocks.size();i++) 
          ss << "\"" << elementBlocks[i] << "\"\n";

       TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
    }
    catch(std::logic_error & e) {
       std::stringstream ss;
       ss << e.what() << "\nUnrecognized logic error.\n";

       TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
    }
    
    std::vector<stk::mesh::Entity*> elements;
    std::vector<std::size_t> local_cell_ids;
    std::vector<std::size_t> local_side_ids;
    getSideElements(mesh, bc->elementBlockID(),
		      sideEntities,local_side_ids,elements);

    // loop over elements of this block
    for(std::size_t elm=0;elm<elements.size();++elm) {
	stk::mesh::Entity * element = elements[elm];
	
	local_cell_ids.push_back(mesh.elementLocalId(element));
    }

    std::map<std::string,panzer::InputPhysicsBlock>::const_iterator ipb_iterator = 
      eb_to_ipb.find(bc->elementBlockID());
    
    TEST_FOR_EXCEPTION(ipb_iterator == eb_to_ipb.end(), std::logic_error,
		       "Could not find input physics block corresponding to region");

    const panzer::InputPhysicsBlock& ipb = ipb_iterator->second;

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
