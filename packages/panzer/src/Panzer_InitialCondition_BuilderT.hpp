#ifndef PANZER_INITIAL_CONDITION_BUILDER_T_HPP
#define PANZER_INITIAL_CONDITION_BUILDER_T_HPP

#include "Teuchos_TestForException.hpp"

template<typename LO, typename GO>
void panzer::setupInitialConditionFieldManagers(const std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > >& volume_worksets,
						const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
						const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
						const Teuchos::ParameterList& ic_block_closure_models,
						const Teuchos::RCP<panzer::UniqueGlobalIndexer<LO,GO> >& dofManager,
						const panzer::LinearObjFactory<panzer::Traits>& lo_factory,
						const Teuchos::ParameterList& user_data,
						std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >& phx_ic_field_managers)

{
  
  std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator blkItr;
  for (blkItr=physicsBlocks.begin();blkItr!=physicsBlocks.end();++blkItr) {
    RCP<panzer::PhysicsBlock> pb = *blkItr;
    std::string blockId = pb->elementBlockID();

    // build a field manager object
    Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm 
          = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>);
    
    // Choose model sublist for this element block
    std::string closure_model_name = "";
    if (ic_block_closure_models.isSublist(blockId))
      closure_model_name = blockId;
    else if (ic_block_closure_models.isSublist("Default"))
      closure_model_name = "Default";
    else 
      TEST_FOR_EXCEPTION(true, std::logic_error, "Failed to find initial condition for element block \"" << blockId << "\".  You must provide an initial condition for each element block or set a default!");

    // use the physics block to register evaluators
    pb->buildAndRegisterInitialConditionEvaluators(*fm, cm_factory, closure_model_name, ic_block_closure_models, lo_factory, user_data);
    //pb->buildAndRegisterClosureModelEvaluators(*fm, cm_factory, closure_models);

    // build the setup data using passed in information
    Traits::SetupData setupData;
    setupData.worksets_ = volume_worksets.find(blockId)->second;

    fm->postRegistrationSetup(setupData);
    phx_ic_field_managers.push_back(fm);
  }
}

void panzer::evaluateInitialCondition(const std::vector< Teuchos::RCP<std::vector<panzer::Workset> > >& worksets,
				      const std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >& phx_ic_field_managers,
				      Teuchos::RCP<panzer::LinearObjContainer> loc,
				      const double time_stamp,
				      const bool write_graphviz_file)
{   
  for (std::size_t block = 0; block < worksets.size(); ++block) {
    
    std::vector<panzer::Workset>& w = *worksets[block]; 
    
    Teuchos::RCP< PHX::FieldManager<panzer::Traits> > fm = 
      phx_ic_field_managers[block];
    
    fm->writeGraphvizFile(std::string("IC_"+block));
    
    // Loop over worksets in this element block
    for (std::size_t i = 0; i < w.size(); ++i) {
      panzer::Workset& workset = w[i];
      
      workset.linContainer = loc;
      // Need to figure out how to get restart time from Rythmos.
      workset.time = time_stamp;
      
      fm->evaluateFields<panzer::Traits::Residual>(workset);
    }
  }
  
}

#endif
