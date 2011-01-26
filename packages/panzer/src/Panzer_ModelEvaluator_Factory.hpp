#ifndef PANZER_MODEL_EVALUATOR_FACTORY_HPP
#define PANZER_MODEL_EVALUATOR_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"

namespace Thyra {
  template<typename Scalar> class ModelEvaluator;
}

namespace panzer {
  
  template <typename ScalarT>
  class ModelEvaluator_Factory : 
    public Teuchos::ParameterListAcceptorDefaultBase{

    void setParameterList(RCP<ParameterList> const& paramList);
    
    RCP<const ParameterList> getValidParameters() const;
    
    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > buildModelEvaluator() const;
    
  };

}

#endif
