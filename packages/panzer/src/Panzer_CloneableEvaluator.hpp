#ifndef PANZER_CloneableEvaluator_H
#define PANZER_CloneableEvaluator_H

#include "Teuchos_RCP.hpp"

namespace panzer {

  //! Non-templated empty base class for template managers
  class CloneableEvaluator {
    
  public:
    
    CloneableEvaluator() {}
    
    virtual ~CloneableEvaluator() {}
    
    virtual Teuchos::RCP<CloneableEvaluator> clone(const Teuchos::ParameterList & pl) const = 0;
  };
  
}

#endif
