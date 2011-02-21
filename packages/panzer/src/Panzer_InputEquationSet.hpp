
#ifndef PANZER_INPUT_EQUATION_SET
#define PANZER_INPUT_EQUATION_SET

#include <vector>
#include <string>
#include <iostream>
#include "Teuchos_ParameterList.hpp"

namespace panzer {

  struct InputEquationSet {
    std::string name;
    std::string basis;
    int integration_order;
    std::string model_id;
    std::string prefix;
    Teuchos::ParameterList params;

    InputEquationSet();

    InputEquationSet(const Teuchos::ParameterList& p);
    
  private:

    void validateParameters(const Teuchos::ParameterList& p) const;
  };

}

#endif
