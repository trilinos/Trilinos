#ifndef USER_APP_CLOSURE_MODEL_FACTORY_T_HPP
#define USER_APP_CLOSURE_MODEL_FACTORY_T_HPP

#include <iostream>
#include <sstream>
#include <typeinfo>
#include "Panzer_InputEquationSet.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Basis.hpp"
#include "Panzer_Integrator_Scalar.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_TypeNameTraits.hpp"

// User application evaluators for this factory
#include "user_app_ConstantModel.hpp"
#include "Panzer_GlobalStatistics.hpp"

// ********************************************************************
// ********************************************************************
template<typename EvalT>
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
user_app::MyModelFactory<EvalT>::
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

  if (!models.isSublist(model_id)) {
    std::stringstream msg;
    msg << "Falied to find requested model, \"" << model_id 
	<< "\", for equation set:\n" << std::endl;
    TEST_FOR_EXCEPTION(!models.isSublist(model_id), std::logic_error, msg.str());
  }

  const ParameterList& my_models = models.sublist(model_id);

  for (ParameterList::ConstIterator model_it = my_models.begin(); 
       model_it != my_models.end(); ++model_it) {
    
    bool found = false;
    
    const std::string key = model_it->first;
    ParameterList input;
    const Teuchos::ParameterEntry& entry = model_it->second;
    const ParameterList& plist = Teuchos::getValue<Teuchos::ParameterList>(entry);
    
    if (plist.isType<double>("Value")) {
      { // at IP
	input.set("Name", key);
	input.set("Value", plist.get<double>("Value"));
	input.set("Data Layout", default_params.get<RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
	RCP< Evaluator<panzer::Traits> > e = 
	  rcp(new user_app::ConstantModel<EvalT,panzer::Traits>(input));
	evaluators->push_back(e);
      }
      { // at BASIS
	input.set("Name", key);
	input.set("Value", plist.get<double>("Value"));
	input.set("Data Layout", default_params.get<RCP<panzer::Basis> >("Basis")->functional);
	RCP< Evaluator<panzer::Traits> > e = 
	  rcp(new user_app::ConstantModel<EvalT,panzer::Traits>(input));
	evaluators->push_back(e);
      }
      found = true;
    }

    if (plist.isType<std::string>("Value")) {
    
      const std::string value = plist.get<std::string>("Value");

      if (key == "Global Statistics") {
	if (typeid(EvalT) == typeid(panzer::Traits::Residual)) {
	  input.set("Comm", user_data.get<Teuchos::RCP<const Teuchos::Comm<int> > >("Comm"));
	  input.set("Names", value);
	  input.set("IR", default_params.get<RCP<panzer::IntegrationRule> >("IR"));
	  RCP< panzer::GlobalStatistics<EvalT,panzer::Traits> > e = 
	    rcp(new panzer::GlobalStatistics<EvalT,panzer::Traits>(input));
	  evaluators->push_back(e);
	  
	  // Require certain fields be evaluated
	  fm.template requireField<EvalT>(e->getRequiredFieldTag());
	}
	found = true;
      }

    }

    if (key == "Volume Integral") {

        {
	   input.set("Name", "Unit Value");
	   input.set("Value", 1.0);
	   input.set("Data Layout", default_params.get<RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
	   RCP< Evaluator<panzer::Traits> > e = 
   	     rcp(new user_app::ConstantModel<EvalT,panzer::Traits>(input));
   	   evaluators->push_back(e);
        }

        {
	   input.set("Integral Name", "Volume_Integral");
	   input.set("Integrand Name", "Unit Value");
	   input.set("IR", default_params.get<RCP<panzer::IntegrationRule> >("IR"));

	   RCP< Evaluator<panzer::Traits> > e = 
   	     rcp(new panzer::Integrator_Scalar<EvalT,panzer::Traits>(input));
   	   evaluators->push_back(e);
        }
    }

    if (!found) {
      std::stringstream msg;
      msg << "ClosureModelFactory failed to build evaluator for key \"" << key 
	  << "\"\nin model \"" << model_id 
	  << "\".  Please correct the type or add support to the \nfactory." <<std::endl;
      TEST_FOR_EXCEPTION(!found, std::logic_error, msg.str());
    }

  }

  return evaluators;
}

#endif
