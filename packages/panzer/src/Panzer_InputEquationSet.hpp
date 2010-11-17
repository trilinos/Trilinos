
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
    int model_id;
    std::string model_factory;
    std::string prefix;
    Teuchos::ParameterList params;
  };

}

#endif
