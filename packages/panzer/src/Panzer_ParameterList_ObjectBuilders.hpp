#ifndef PANZER_PARAMETER_LIST_OBJECT_BUILDERS_HPP
#define PANZER_PARAMETER_LIST_OBJECT_BUILDERS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include <vector>
#include <string>

namespace panzer {

  class InputPhysicsBlock;
  class BC;

  void buildInputPhysicsBlocks(std::vector<panzer::InputPhysicsBlock>& ipb,
			       const Teuchos::ParameterList& p);
  
  void buildBCs(std::vector<panzer::BC>& bcs, 
		const Teuchos::ParameterList& p);
  
  void StringTokenizer(std::vector<std::string>& tokens,
		       const std::string& str,
		       const std::string delimiter = ",");
  
}

#endif
