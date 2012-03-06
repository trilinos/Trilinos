#ifndef PANZER_CLOSURE_MODEL_FACTORY_COMPOSITE_IMPL_HPP
#define PANZER_CLOSURE_MODEL_FACTORY_COMPOSITE_IMPL_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include "Panzer_InputEquationSet.hpp"
#include "Panzer_GlobalData.hpp"
#include "Teuchos_ParameterEntry.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Panzer_String_Utilities.hpp"

// ********************************************************************
template<typename EvalT>
panzer::ClosureModelFactoryComposite<EvalT>::
ClosureModelFactoryComposite(const std::vector<Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > >& factories) :
  m_factories(factories)
{ }

// ********************************************************************
template<typename EvalT>
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
panzer::ClosureModelFactoryComposite<EvalT>::
buildClosureModels(const std::string& model_id,
		   const panzer::InputEquationSet& set,
		   const Teuchos::ParameterList& models, 
		   const Teuchos::ParameterList& default_params, 
		   const Teuchos::ParameterList& user_data,
		   const Teuchos::RCP<panzer::GlobalData>& global_data,
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
	<< "\" for equation set:\n" << std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION(!models.isSublist(model_id), std::logic_error, msg.str());
  }

  const ParameterList& my_models = models.sublist(model_id);

  // for each model, try each factory until you get a nonnull
  // return, meaning you have built the evaluators for that model
  for (ParameterList::ConstIterator model_it = my_models.begin(); 
       model_it != my_models.end(); ++model_it) {
    
    bool found = false;

    std::string model_key = model_it->first;
    
    // Duplicate models sublist with just the particular model you
    // want to build
    ParameterList copy_of_models;
    Teuchos::ParameterList* tmp;
    copy_of_models.sublist(model_id).sublist(model_key) = model_it->second.getValue(tmp);
    
    // Loop over factories
    for (std::vector<Teuchos::RCP<panzer::ClosureModelFactory_TemplateManager<panzer::Traits> > >::const_iterator factory = m_factories.begin(); factory != m_factories.end(); ++factory) {

      RCP< vector< RCP<Evaluator<panzer::Traits> > > > tmp_evaluators =
	(*factory)->getAsObject<EvalT>()->buildClosureModels(model_id,set,copy_of_models,default_params,user_data,global_data,fm);

      if (tmp_evaluators->size() > 0) {
	
	for (vector< RCP<Evaluator<panzer::Traits> > >::const_iterator eval = tmp_evaluators->begin(); eval != tmp_evaluators->end(); ++eval)
	  evaluators->push_back(*eval);

	found = true;
	break;
      }
      
    }
    
    if (!found) {
      std::stringstream msg;
      msg << "ERROR: ClosureModelFactoryComposite failed to build model \"" << model_id << "\",\n"
	  << "with the sublist key \"" << model_key << "\".\n"
	  << "Please correct the model input or add support for this key to a\nclosure "
	  << "model factory." << std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(!found, std::logic_error, msg.str());
    }

  }

  return evaluators;
}

#endif
