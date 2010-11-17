
#ifndef PANZER_INPUT_PHYSICS_BLOCK
#define PANZER_INPUT_PHYSICS_BLOCK

#include <vector>
#include <string>
#include <iostream>
#include "Teuchos_ParameterList.hpp"
#include "Panzer_InputEquationSet.hpp"

namespace panzer {

  struct InputPhysicsBlock {

    std::string physics_block_id;
    
    std::vector<panzer::InputEquationSet> eq_sets;
 
  };

  std::ostream& operator<<(std::ostream& os, 
			   const panzer::InputPhysicsBlock& i);
}

#endif
