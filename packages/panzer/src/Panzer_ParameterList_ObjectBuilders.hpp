#ifndef PANZER_PARAMETER_LIST_OBJECT_BUILDERS_HPP
#define PANZER_PARAMETER_LIST_OBJECT_BUILDERS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include <map>
#include <vector>
#include <string>

namespace panzer {

  class InputPhysicsBlock;
  class BC;

  void buildInputPhysicsBlocks(std::map<std::string,panzer::InputPhysicsBlock>& ipb,
			       const Teuchos::ParameterList& p);
  
  void buildBCs(std::vector<panzer::BC>& bcs, 
		const Teuchos::ParameterList& p);

  void buildBlockIdToPhysicsIdMap(std::map<std::string,std::string>& b_to_p, 
				  const Teuchos::ParameterList& p);

}

#endif
