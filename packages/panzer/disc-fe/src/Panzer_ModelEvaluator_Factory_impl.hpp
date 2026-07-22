// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_MODEL_EVALUATOR_FACTORY_IMPL_HPP
#define PANZER_MODEL_EVALUATOR_FACTORY_IMPL_HPP

#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Panzer_ModelEvaluator_Epetra.hpp"

namespace panzer {
  
  template <typename ScalarT, typename LO, typename GO>
  Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > 
  buildModelEvaluator() const
  { 
    Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > me;
    
    std::string type = getMyParamList.get<std::string>("Model Evaluator Type");
    
    if (type == "Epetra") {
      
      Teuchos::RCP<EpetraExt::ModelEvaluator> epetraModel = 
	Teuchos::rcp(new panzer::ModelEvaluator_Epetra(fmb, linObjFactory));
    
      Teuchos::RCP<Thyra::EpetraModelEvaluator>
	epetraThyraModel = rcp(new ::Thyra::EpetraModelEvaluator());
      epetraThyraModel->initialize(epetraModel,lowsFactory);
      Teuchos::RCP< ::Thyra::ModelEvaluator<double> > thyraModel = 
	epetraThyraModel;

    }
    else if (type == "Tpetra") {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			 "Tpetra version not supported yet, use Epetra!");
    }
    


  }

}

#endif
