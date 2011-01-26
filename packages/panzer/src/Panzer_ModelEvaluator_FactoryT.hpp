#ifndef PANZER_MODEL_EVALUATOR_FACTORY_T_HPP
#define PANZER_MODEL_EVALUATOR_FACTORY_T_HPP

#include "Thyra_ModelEvaluator.hpp"

namespace panzer {
  
  template <typename ScalarT>
  void setParameterList(RCP<ParameterList> const& paramList)
  {
    using Teuchos:RCP;
    using Teuchos::ParameterList;
    
    RCP<const ParameterList> valid_parameters = this->getValidParameters();
    
    paramList->validateParametersAndSetDefaults(valid_parameters);
    
    this->setMyParamList(paramList);
  }

  template <typename ScalarT>
  RCP<const ParameterList> getValidParameters() const
  {
    using Teuchos:RCP;
    using Teuchos::ParameterList;
   
    RCP<ParameterList> valid_parameters;
    
    Teuchos::setStringToIntegralParameter<int>(
      "Model Evaluator Type",
      "Epetra",
      "Determines the type of model evaluator to build.",
      Teuchos::tuple<std::string>("Epetra","Tpetra"),
      &valid_parameters
      );

  }

  template <typename ScalarT>
  Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > buildModelEvaluator() const
  {

  }

}

#endif
