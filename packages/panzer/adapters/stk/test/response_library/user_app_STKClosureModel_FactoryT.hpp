#ifndef USER_APP_STK_CLOSURE_MODEL_FACTORY_T_HPP
#define USER_APP_STK_CLOSURE_MODEL_FACTORY_T_HPP

#include <iostream>
#include <sstream>
#include "Panzer_InputEquationSet.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_TypeNameTraits.hpp"

// User application evaluators for this factory
#include "user_app_ConstantModel.hpp"
#include "TestEvaluators.hpp"

// ********************************************************************
// ********************************************************************
template<typename EvalT>
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
user_app::STKModelFactory<EvalT>::
buildClosureModels(const std::string& model_id,
		   const panzer::InputEquationSet& set,
		   const Teuchos::ParameterList& models, 
		   const Teuchos::ParameterList& default_params,
		   const Teuchos::ParameterList& user_data,
		   PHX::FieldManager<panzer::Traits>& fm) const
{

  using std::string;
  using std::vector;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using PHX::Evaluator;

  RCP< vector< RCP<Evaluator<panzer::Traits> > > > evaluators = 
    rcp(new vector< RCP<Evaluator<panzer::Traits> > > );

  {
     RCP<PHX::Evaluator<panzer::Traits> > e
        = rcp(new panzer::TestEvaluator<EvalT,panzer::Traits>(user_data));

     evaluators->push_back(e);
  }


  return evaluators;
}

#endif
