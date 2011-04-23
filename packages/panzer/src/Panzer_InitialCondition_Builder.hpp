#ifndef PANZER_INITIAL_CONDITION_BUILDER_HPP
#define PANZER_INITIAL_CONDITION_BUILDER_HPP

#include "Panzer_Traits.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_ClosureModel_Factory_TemplateManager.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include <vector>

namespace panzer {

  /** \brief Builds PHX::FieldManager objects for inital conditions and registers evaluators. 

  \param[in] volume_worksets Worksets for the corresponding evaluations.
  \param[in] physicsBlock PhysicsBlocks created by FieldManagerBuilder.
  \param[in] cm_factory Factory that provides all closure models required by the initial condition evaluation.
  \param[in] closure_models List of closure models for each element block required for initial conditions.
  \param[in] dofManager Degree of Freedom manager corresponding to the volume problem.
  \param[in] lo_factory LinearObjFactory corresponding to the problem.
  \param[in] user_data ParameterList with optional user data.
  \param[out] phx_ic_field_managers Allocated PHX::FieldManagers for each element block.

  */
  template <typename LO, typename GO>
  void setupInitialConditionFieldManagers(const std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > >& volume_worksets,
					  const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
					  const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
					  const Teuchos::ParameterList& ic_block_closure_models,
					  const Teuchos::RCP<panzer::UniqueGlobalIndexer<LO,GO> >& dofManager,
					  const panzer::LinearObjFactory<panzer::Traits>& lo_factory,
					  const Teuchos::ParameterList& user_data,
					  const bool write_graphviz_file,
					  const std::string& graphviz_file_prefix,
					  std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >& phx_ic_field_managers);

  void evaluateInitialCondition(const std::vector< Teuchos::RCP<std::vector<panzer::Workset> > >& worksets,
				const std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >& phx_ic_field_managers,
				Teuchos::RCP<panzer::LinearObjContainer> loc,
				const double time_stamp);
  
}

#include "Panzer_InitialCondition_BuilderT.hpp"

#endif
