#include <vector>
#include <string>
#include <sstream>
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_DOFManager.hpp"
#include "Panzer_ConnManager.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_Shards_Utilities.hpp"
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"
#include "Panzer_CellData.hpp"
#include "Shards_CellTopology.hpp"
#include "Panzer_InputPhysicsBlock.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Panzer_StlMap_Utilities.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"

//#include "EpetraExt_BlockMapOut.h"

//=======================================================================
//=======================================================================
template<typename LO, typename GO>
void panzer::FieldManagerBuilder<LO,GO>::print(std::ostream& os) const
{
  os << "panzer::FieldManagerBuilder<LO,GO> output:  Not implemented yet!";
}

//=======================================================================
//=======================================================================
template<typename LO, typename GO>
void panzer::FieldManagerBuilder<LO,GO>::setupVolumeFieldManagers( 
                                            const std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > >& volume_worksets, 
                                                                       // element block -> vector of worksets
                                            const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks, 
					    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
					    const Teuchos::ParameterList& closure_models,
                                            const panzer::LinearObjFactory<panzer::Traits> & lo_factory,
					    const Teuchos::ParameterList& user_data)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  worksets_.clear();
  phx_volume_field_managers_.clear();

  std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator blkItr;
  for (blkItr=physicsBlocks.begin();blkItr!=physicsBlocks.end();++blkItr) {
    RCP<panzer::PhysicsBlock> pb = *blkItr;
    std::string blockId = pb->elementBlockID();

    // build a field manager object
    Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm 
          = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>);
    
    // use the physics block to register evaluators
    pb->buildAndRegisterEquationSetEvaluators(*fm, user_data);
    pb->buildAndRegisterGatherScatterEvaluators(*fm,lo_factory, user_data);
    pb->buildAndRegisterClosureModelEvaluators(*fm, cm_factory, closure_models, user_data);

    // build the setup data using passed in information
    Traits::SetupData setupData;
    setupData.worksets_ = volume_worksets.find(blockId)->second;
    fm->postRegistrationSetup(setupData);

    // make sure to add the field manager & workset to the list 
    worksets_.push_back(setupData.worksets_);
    phx_volume_field_managers_.push_back(fm); 
  }
}
//=======================================================================
//=======================================================================
template<typename LO, typename GO>
void panzer::FieldManagerBuilder<LO,GO>::setupBCFieldManagers(
                           const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC>& bc_worksets,
                                                       // boundary condition -> map of (side_id,worksets)
                           const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
	                   const panzer::EquationSetFactory & eqset_factory,
			   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
                           const panzer::BCStrategyFactory& bc_factory,
			   const Teuchos::ParameterList& closure_models,
                           const panzer::LinearObjFactory<panzer::Traits> & lo_factory,
			   const Teuchos::ParameterList& user_data)
{

  // for convenience build a map (element block id => physics block)
  std::map<std::string,Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks_map;
  {
     std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator blkItr;
     for(blkItr=physicsBlocks.begin();blkItr!=physicsBlocks.end();++blkItr) {
        Teuchos::RCP<panzer::PhysicsBlock> pb = *blkItr;
        std::string blockId = pb->elementBlockID();

        // add block id, physics block pair to the map
        physicsBlocks_map.insert(std::make_pair(blockId,pb));
     }
  }

  // ***************************
  // BCs
  // ***************************
  std::map<panzer::BC,Teuchos::RCP<BCFaceWorksetMap>,panzer::LessBC>::const_iterator bc = 
    bc_worksets.begin();
  for (; bc != bc_worksets.end(); ++bc) {
    std::string element_block_id = bc->first.elementBlockID(); 
    Teuchos::RCP<const panzer::PhysicsBlock> volume_pb = physicsBlocks_map.find(element_block_id)->second;
    const shards::CellTopology volume_cell_topology = volume_pb->getBaseCellTopology();
    int base_cell_dimension = volume_pb->cellData().baseCellDimension();
    
    bc_worksets_[bc->first] = bc->second; 

    // Build one FieldManager for each local side workset for each dirichlet bc
    std::map<unsigned,PHX::FieldManager<panzer::Traits> >& field_managers = 
      bc_field_managers_[bc->first];

    // Loop over local face indices and setup each field manager
    for (std::map<unsigned,panzer::Workset>::const_iterator wkst = 
	   bc_worksets_[bc->first]->begin(); wkst != bc_worksets_[bc->first]->end();
	 ++wkst) {

      PHX::FieldManager<panzer::Traits>& fm = field_managers[wkst->first];
      
      // register evaluators from strategy
      Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs = 
	bc_factory.buildBCStrategy(bc->first);
      
      const panzer::CellData side_cell_data(wkst->second.num_cells,
					    base_cell_dimension,
					    wkst->first);      

      Teuchos::RCP<panzer::PhysicsBlock> side_pb 
            = volume_pb->copyWithCellData(side_cell_data, eqset_factory);
      
      // Iterate over evaluation types
      for (panzer::BCStrategy_TemplateManager<panzer::Traits>::iterator 
	     bcs_type = bcs->begin(); bcs_type != bcs->end(); ++bcs_type) {
	bcs_type->setup(*side_pb,user_data);
	bcs_type->buildAndRegisterEvaluators(fm,*side_pb,cm_factory,closure_models,user_data);
	bcs_type->buildAndRegisterGatherScatterEvaluators(fm,*side_pb,lo_factory,user_data);
      }

      // Setup the fieldmanager
      Traits::SetupData setupData;
      Teuchos::RCP<std::vector<panzer::Workset> > worksets = 
	Teuchos::rcp(new(std::vector<panzer::Workset>));
      worksets->push_back(wkst->second);
      setupData.worksets_ = worksets;
      fm.postRegistrationSetup(setupData);
    }
    
  }
}

//=======================================================================
//=======================================================================
template<typename LO, typename GO>
void panzer::FieldManagerBuilder<LO,GO>::
writeVolumeGraphvizDependencyFiles(std::string filename_prefix,
				   const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks) const
{  
  std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator blkItr;
  int index = 0;
  for (blkItr=physicsBlocks.begin();blkItr!=physicsBlocks.end();++blkItr,++index) {
    std::string blockId = (*blkItr)->elementBlockID();
    phx_volume_field_managers_[index]->writeGraphvizFile(filename_prefix+blockId);
  }

}

//=======================================================================
//=======================================================================
template<typename LO, typename GO>
std::ostream& panzer::operator<<(std::ostream& os, const panzer::FieldManagerBuilder<LO,GO>& rfd)
{
  rfd.print(os);
  return os;
}
