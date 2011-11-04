#include "Panzer_ParameterList_ObjectBuilders.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_InputPhysicsBlock.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include <sstream>

namespace panzer {
  
  void buildInputPhysicsBlocks(std::map<std::string,panzer::InputPhysicsBlock>& ipb,
			       const Teuchos::ParameterList& p)
  {
    using std::string;
    using std::vector;
    using Teuchos::ParameterList;

    string block_names = p.get<string>("Physics Blocks");
    
    vector<string> block_names_vec;
    panzer::StringTokenizer(block_names_vec, block_names);
    
    // Validate that the sublists are physics block names
    ParameterList valid_params;
    valid_params.set("Physics Blocks", block_names);
    for (vector<string>::const_iterator block = block_names_vec.begin();
	 block != block_names_vec.end(); ++block)
      valid_params.sublist(*block);
    
    p.validateParameters(valid_params,0);
    
    
    for (vector<string>::const_iterator block = block_names_vec.begin();
	 block != block_names_vec.end(); ++block) {
      
      const ParameterList& block_list = p.sublist(*block); 
      int num_eqsets = block_list.get<int>("Number of Equation Sets");
      ipb[*block].physics_block_id = *block;

      for (int set=0; set < num_eqsets; ++set) {
	std::ostringstream set_name;
	set_name << "EQ " << set;
	const ParameterList& set_list = block_list.sublist(set_name.str());
	panzer::InputEquationSet ies(set_list);
	ipb[*block].eq_sets.push_back(ies);
      }
    }

  }
  
  void buildBCs(std::vector<panzer::BC>& bcs, 
		const Teuchos::ParameterList& p)
  {
    using Teuchos::ParameterList;

    int num_bcs = p.get<int>("Number of Boundary Conditions");
    for (int bc_id=0; bc_id < num_bcs; ++bc_id) {
      std::ostringstream bc_name;
      bc_name << "BC " << bc_id;
      const ParameterList& bc_list = p.sublist(bc_name.str());
      panzer::BC bc(bc_list);
      bcs.push_back(bc);
    }
  }

  void buildBlockIdToPhysicsIdMap(std::map<std::string,std::string>& b_to_p,
				  const Teuchos::ParameterList& p)
  {
    for (Teuchos::ParameterList::ConstIterator entry = p.begin();
	 entry != p.end(); ++entry) {
      std::string dummy_type;
      b_to_p[entry->first] = entry->second.getValue(&dummy_type);
    }
  }

  void buildPhysicsBlocks(const std::map<std::string,std::string>& block_ids_to_physics_ids,
                          const std::map<std::string,Teuchos::RCP<const shards::CellTopology> > & block_id_to_cell_topo,
                          const std::map<std::string,panzer::InputPhysicsBlock>& physics_id_to_input_physics_blocks,
                          const int base_cell_dimension,
                          const std::size_t workset_size,
                          const panzer::EquationSetFactory & eqset_factory,
                          const bool build_transient_support,
                          std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks)
  {
     using Teuchos::RCP;
     using Teuchos::rcp;
  
     // loop over all block id physics id pairs
     std::map<std::string,std::string>::const_iterator itr;
     for (itr = block_ids_to_physics_ids.begin(); itr!=block_ids_to_physics_ids.end();++itr) {
        std::string element_block_id = itr->first;
        std::string physics_block_id = itr->second;
   
        const panzer::CellData volume_cell_data(workset_size, base_cell_dimension);
        
        // find InputPhysicsBlock that corresponds to a paricular block ID
        std::map<std::string,panzer::InputPhysicsBlock>::const_iterator ipb_it = 
              physics_id_to_input_physics_blocks.find(physics_block_id);
  
        // sanity check: passes only if there is a paricular physics ID
        TEUCHOS_TEST_FOR_EXCEPTION(ipb_it == physics_id_to_input_physics_blocks.end(),
  			 std::runtime_error,
  			 "Falied to find InputPhysicsBlock for physics id: "
  			 << physics_block_id << "!");
  
        const panzer::InputPhysicsBlock& ipb = ipb_it->second;
        RCP<panzer::PhysicsBlock> pb = 
  	rcp(new panzer::PhysicsBlock(ipb, element_block_id, volume_cell_data, eqset_factory, build_transient_support));
        physicsBlocks.push_back(pb);
     }
  }
}
