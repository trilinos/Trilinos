#ifndef PANZER_WORKSET_UTILITIES_HPP
#define PANZER_WORKSET_UTILITIES_HPP

#include "Panzer_Traits.hpp"
#include "Panzer_Workset.hpp"
#include "Teuchos_Assert.hpp"
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>

namespace panzer {

  std::vector<std::string>::size_type 
  getBasisIndex(std::string basis_name, panzer::Workset& workset)
  {
    std::vector<std::string>::iterator basis;

    basis = std::find(workset.basis_names->begin(),
		      workset.basis_names->end(),
		      basis_name);

    TEUCHOS_TEST_FOR_EXCEPTION(basis == workset.basis_names->end(),
			       std::logic_error,
			       "Could not find the basis named \"" 
                               << basis_name << "\" in the workset!");

    return std::distance(workset.basis_names->begin(), basis);
  }

  std::vector<std::string>::size_type
  getIntegrationRuleIndex(int ir_degree, panzer::Workset& workset)
  {
    std::vector<int>::iterator ir;

    ir = std::find(workset.ir_degrees->begin(),
		   workset.ir_degrees->end(),
		   ir_degree);
    
    TEUCHOS_TEST_FOR_EXCEPTION(ir == workset.ir_degrees->end(),
			       std::logic_error,
			       "Could not find the integration rule degree \"" 
                               << ir_degree << "\" in the workset!");

    return std::distance(workset.ir_degrees->begin(), ir);
  }

}

#endif
