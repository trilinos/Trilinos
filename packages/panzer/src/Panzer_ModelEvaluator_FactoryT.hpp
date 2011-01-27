#ifndef PANZER_MODEL_EVALUATOR_FACTORY_T_HPP
#define PANZER_MODEL_EVALUATOR_FACTORY_T_HPP

#include "Thyra_ModelEvaluator.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

namespace panzer {
  
  template <typename ScalarT>
  void ModelEvaluator_Factory<ScalarT>::setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList)
  {
    using Teuchos::RCP;
    using Teuchos::ParameterList;
    
    RCP<const ParameterList> valid_parameters = this->getValidParameters();
    
    paramList->validateParametersAndSetDefaults(*valid_parameters);
    
    this->setMyParamList(paramList);
  }

  template <typename ScalarT>
  Teuchos::RCP<const Teuchos::ParameterList> ModelEvaluator_Factory<ScalarT>::getValidParameters() const
  {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ParameterList;
   
    RCP<ParameterList> valid_parameters = rcp(new ParameterList);
    
    Teuchos::setStringToIntegralParameter<int>(
      "Model Evaluator Type",
      "Epetra",
      "Determines the type of model evaluator to build.",
      Teuchos::tuple<std::string>("Epetra","Tpetra"),
      valid_parameters.get()
      );

    return valid_parameters;
  }

  template <typename ScalarT>
  Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > ModelEvaluator_Factory<ScalarT>::buildModelEvaluator() const
  {
    std::string type = getMyParamList.get<std::string>("Model Evaluator Type");
    
    if (type == "Epetra") {
      
    }
    else if (type == "Tpetra") {

    }
    


  }

}

#endif
