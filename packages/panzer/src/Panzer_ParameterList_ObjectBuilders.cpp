#include "Panzer_ParameterList_ObjectBuilders.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_InputPhysicsBlock.hpp"
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


}
