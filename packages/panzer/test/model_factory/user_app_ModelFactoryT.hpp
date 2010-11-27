#include <iostream>
#include <sstream>
#include "Panzer_InputEquationSet.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"

// User application evaluators for this factory
#include "user_app_ConstantModel.hpp"

// ********************************************************************
// ********************************************************************
template<typename EvalT>
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
user_app::MyModelFactory<EvalT>::
buildModels(const panzer::InputEquationSet& set,
	    const std::vector<Teuchos::ParameterList>& models, 
	    const Teuchos::ParameterList& default_params) const
{

  using std::string;
  using std::vector;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using PHX::Evaluator;

  RCP< vector< RCP<Evaluator<panzer::Traits> > > > evaluators = 
    rcp(new vector< RCP<Evaluator<panzer::Traits> > > );

  for (vector<ParameterList>::const_iterator model = models.begin(); 
       model != models.end(); ++model) {
    
    bool found = false;
    
    ParameterList plist = *model;

    if (model->get<string>("Object Type") == "Constant") {
      plist.set("Data Layout", default_params.get<RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
      RCP< Evaluator<panzer::Traits> > e = 
	rcp(new user_app::ConstantModel<EvalT,panzer::Traits>(plist));
      evaluators->push_back(e);
      found = true;
    }

    if (!found) {
      std::stringstream msg;
      msg << "Evaluator could not be created for \"Object Type\": " 
	  << model->get<string>("Object Type") << std::endl;
      TEST_FOR_EXCEPTION(true, std::logic_error, msg.str());
    }

  }

  return evaluators;
}
