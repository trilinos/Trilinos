#ifndef PANZER_WORKSET_UTILITIES_HPP
#define PANZER_WORKSET_UTILITIES_HPP

#include "Panzer_Traits.hpp"
#include "Panzer_Workset.hpp"
#include <vector>

namespace panzer {

  std::vector<std::string>::size_type 
  getBasisIndex(std::string basis_name, panzer::Workset& workset);

  std::vector<int>::size_type
  getIntegrationRuleIndex(int ir_degree, panzer::Workset& workset);

}

#endif
