#ifndef PANZER_INITIAL_CONDITION_BUILDER_T_HPP
#define PANZER_INITIAL_CONDITION_BUILDER_T_HPP

#include "Teuchos_TestForException.hpp"

template<typename LO, typename GO>
void panzer::setupInitialConditionFieldManagers(const std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > >& volume_worksets,
						const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
						const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
						const Teuchos::ParameterList& closure_models,
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
    
    // use the physics block to register evaluators
    pb->buildAndRegisterInitialConditionEvaluators(*fm, cm_factory, closure_models, lo_factory, user_data);
    //pb->buildAndRegisterClosureModelEvaluators(*fm, cm_factory, closure_models);

    // build the setup data using passed in information
    Traits::SetupData setupData;
    setupData.worksets_ = volume_worksets.find(blockId)->second;

    fm->postRegistrationSetup(setupData);
    phx_ic_field_managers.push_back(fm);
  }
  
}

#endif
