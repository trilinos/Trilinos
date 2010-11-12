
#ifndef PANZER_INPUT_PHYSICS_BLOCK
#define PANZER_INPUT_PHYSICS_BLOCK

#include <vector>
#include <string>
#include <iostream>
#include "Teuchos_ParameterList.hpp"

namespace panzer {

  struct InputPhysicsBlock {

    struct InputEquationSet {
      std::string name;
      std::string basis;
      int integration_order;
      int model_id;
      std::string model_factory;
      std::string prefix;
      Teuchos::ParameterList params;
    };

    //! ID from the input file
    std::string physics_block_id;
    
    std::vector<InputEquationSet> eq_sets;
 
  };

  std::ostream& operator<<(std::ostream& os, 
			   const panzer::InputPhysicsBlock& i);
}

#endif
