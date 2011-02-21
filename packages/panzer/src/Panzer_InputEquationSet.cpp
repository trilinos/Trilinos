#include "Panzer_InputEquationSet.hpp"

namespace panzer {

  InputEquationSet::InputEquationSet()
  {
    
  }
  
  InputEquationSet::InputEquationSet(const Teuchos::ParameterList& p)
  {
    this->validateParameters(p);

    name = p.get<std::string>("Name");
    basis = p.get<std::string>("Basis");
    integration_order = p.get<int>("Integration Order");
    model_id = p.get<std::string>("Model ID");
    prefix = p.get<std::string>("Prefix");
    params = p.sublist("Options");
  }

  void InputEquationSet::
  validateParameters(const Teuchos::ParameterList& p) const
  {
    Teuchos::ParameterList valid_params;
    
    valid_params.set<std::string>("Name", "");
    valid_params.set<std::string>("Basis", "");
    valid_params.set<int>("Integration Order", -1);
    valid_params.set<std::string>("Model ID", "");
    valid_params.set<std::string>("Prefix", "");
    valid_params.sublist("Options");
    
    p.validateParameters(valid_params, 0);
  }

}
