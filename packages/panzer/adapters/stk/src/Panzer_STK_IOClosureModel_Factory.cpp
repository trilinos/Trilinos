#include "Panzer_config.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_STK_IOClosureModel_Factory.hpp"
#include "Panzer_STK_ScatterCellAvgQuantity.hpp"
#include "Panzer_STK_ScatterCellQuantity.hpp"

template< >
Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > > 
panzer_stk::IOClosureModelFactory<panzer::Traits::Residual>::
buildClosureModels(const std::string& model_id,
		   const panzer::InputEquationSet& set,
		   const Teuchos::ParameterList& models, 
		   const Teuchos::ParameterList& default_params, 
		   const Teuchos::ParameterList& user_data,
		   const Teuchos::RCP<panzer::GlobalData>& global_data,
		   PHX::FieldManager<panzer::Traits>& fm) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using PHX::Evaluator;
  static bool justOnce = false;

  // build user evaluators
  RCP< std::vector< RCP<Evaluator<panzer::Traits> > > > user_evals = 
    userCMF_->buildClosureModels(model_id,set,models,default_params,user_data,global_data,fm);

  // add user evaluators to evaluator list
  RCP< std::vector< RCP<Evaluator<panzer::Traits> > > > evaluators = 
    rcp(new std::vector< RCP<Evaluator<panzer::Traits> > > );

  // extract element block id
  std::string block_id = default_params.get<std::string>("Block ID");

  if(!blockIdEvaluated_[block_id]) {
     typedef std::map<std::string,std::vector<std::string> > BlockIdToFields;

     Teuchos::RCP<panzer::IntegrationRule> intRule = default_params.get<Teuchos::RCP<panzer::IntegrationRule> >("IR");
     int worksetsize = intRule->dl_scalar->dimension(0);

     // if a requested field is found then add in cell avg quantity evaluator
     BlockIdToFields::const_iterator cellAvgItr = blockIdToCellAvgFields_.find(block_id);
     if(cellAvgItr!=blockIdToCellAvgFields_.end() ) {
        justOnce = true;
        Teuchos::RCP<std::vector<std::string> > fieldNames = Teuchos::rcp(new std::vector<std::string>(cellAvgItr->second));
   
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
   
        blockIdEvaluated_[block_id] = true;
     } 

     // if a requested field is found then add in cell quantity evaluator
     BlockIdToFields::const_iterator cellItr = blockIdToCellFields_.find(block_id);
     if(cellItr!=blockIdToCellFields_.end() ) {
        justOnce = true;
        Teuchos::RCP<std::vector<std::string> > fieldNames = Teuchos::rcp(new std::vector<std::string>(cellItr->second));
   
        // setup averge cell fields
        Teuchos::ParameterList pl;
        pl.set("Mesh",mesh_);
        pl.set("Workset Size",worksetsize);
        pl.set("Field Names",fieldNames);
        pl.set("Scatter Name", block_id+"_Cell_Fields");
        Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval
            = Teuchos::rcp(new panzer_stk::ScatterCellQuantity<panzer::Traits::Residual,panzer::Traits>(pl));
        fm.registerEvaluator<panzer::Traits::Residual>(eval);
        fm.requireField<panzer::Traits::Residual>(*eval->evaluatedFields()[0]);
   
        evaluators->push_back(eval);
   
        blockIdEvaluated_[block_id] = true;
     } 
  }

  evaluators->insert(evaluators->end(),user_evals->begin(),user_evals->end()); 

  return evaluators;
}

#ifdef HAVE_PANZER_EXPLICIT_INSTANTIATION

#include "Panzer_ExplicitTemplateInstantiation.hpp"

#include "Panzer_STK_IOClosureModel_Factory_decl.hpp"
#include "Panzer_STK_IOClosureModel_Factory_impl.hpp"

PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(panzer_stk::IOClosureModelFactory)

#endif
