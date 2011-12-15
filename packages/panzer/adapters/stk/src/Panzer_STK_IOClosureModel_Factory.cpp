#include "Panzer_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_STK_IOClosureModel_Factory.hpp"
#include "Panzer_STK_ScatterCellAvgQuantity.hpp"

template< >
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
panzer_stk::IOClosureModelFactory<panzer::Traits::Residual>::
buildClosureModels(const std::string& model_id,
		   const panzer::InputEquationSet& set,
		   const Teuchos::ParameterList& models, 
		   const Teuchos::ParameterList& default_params, 
		   const Teuchos::ParameterList& user_data,
		   PHX::FieldManager<panzer::Traits>& fm) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using PHX::Evaluator;

  // build user evaluators
  RCP< std::vector< RCP<Evaluator<panzer::Traits> > > > user_evals = 
     userCMF_->buildClosureModels(model_id,set,models,default_params,user_data,fm);

  // add user evaluators to evaluator list
  RCP< std::vector< RCP<Evaluator<panzer::Traits> > > > evaluators = 
    rcp(new std::vector< RCP<Evaluator<panzer::Traits> > > );

  // extract element block id
  std::string block_id = default_params.get<std::string>("Block ID");

  // if a requested field is found then add in cell avg quantity evaluator
  std::map<std::string,std::vector<std::string> >::const_iterator bi2FieldsItr  
     = blockIdToFields_.find(block_id);
  if(bi2FieldsItr!=blockIdToFields_.end()) {
     Teuchos::RCP<panzer::IntegrationRule> intRule = default_params.get<Teuchos::RCP<panzer::IntegrationRule> >("IR");
     Teuchos::RCP<std::vector<std::string> > fieldNames = Teuchos::rcp(new std::vector<std::string>(bi2FieldsItr->second));

     // setup averge cell fields
     Teuchos::ParameterList pl;
     pl.set("Mesh",mesh_);
     pl.set("IR",intRule);
     pl.set("Field Names",fieldNames);
     pl.set("Scatter Name", block_id+"_Cell_Avg_Fields");
     Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
         = Teuchos::rcp(new panzer_stk::ScatterCellAvgQuantity<panzer::Traits::Residual,panzer::Traits>(pl));
     fm.registerEvaluator<panzer::Traits::Residual>(eval);
     fm.requireField<panzer::Traits::Residual>(*eval->evaluatedFields()[0]);

     evaluators->push_back(eval);
  } 

  evaluators->insert(evaluators->end(),user_evals->begin(),user_evals->end()); 

  return evaluators;
}

#ifdef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_STK_IOClosureModel_FactoryT.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(panzer_stk::IOClosureModelFactory)

#endif
