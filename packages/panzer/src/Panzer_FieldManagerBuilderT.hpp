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
void panzer::FieldManagerBuilder<LO,GO>::buildPhysicsBlocks(const std::map<std::string,std::string>& block_ids_to_physics_ids,
                                                            const std::map<std::string,panzer::InputPhysicsBlock>& physics_id_to_input_physics_blocks,
                                                            int base_cell_dimension, std::size_t workset_size,
	                                                    const panzer::EquationSetFactory & eqset_factory,
                                                            std::vector<Teuchos::RCP<panzer::PhysicsBlock> > & physicsBlocks
                                                            ) const 
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   // loop over all block id physics id pairs
   std::map<std::string,std::string>::const_iterator itr;
   for (itr = block_ids_to_physics_ids.begin(); itr!=block_ids_to_physics_ids.end();++itr) {
      std::string element_block_id = itr->first;
      std::string physics_block_id = itr->second;
 
      const panzer::CellData volume_cell_data(workset_size, base_cell_dimension);
      
      // find InputPhysicsBlock that corresponds to a paricular block ID
      std::map<std::string,panzer::InputPhysicsBlock>::const_iterator ipb_it = 
            physics_id_to_input_physics_blocks.find(physics_block_id);

      // sanity check: passes only if there is a paricular physics ID
      TEST_FOR_EXCEPTION(ipb_it == physics_id_to_input_physics_blocks.end(),
			 std::runtime_error,
			 "Falied to find InputPhysicsBlock for physics id: "
			 << physics_block_id << "!");

      const panzer::InputPhysicsBlock& ipb = ipb_it->second;
      RCP<panzer::PhysicsBlock> pb = 
	rcp(new panzer::PhysicsBlock(ipb, element_block_id, volume_cell_data, eqset_factory));
      physicsBlocks.push_back(pb);
   }
}


//=======================================================================
//=======================================================================
template<typename LO, typename GO>
void panzer::FieldManagerBuilder<LO,GO>::setupVolumeFieldManagers( 
                                            const std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > >& volume_worksets, 
                                                                       // element block -> vector of worksets
                                            const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks, 
                                            const Teuchos::RCP<panzer::UniqueGlobalIndexer<LO,GO> > & dofManager,
                                            const panzer::LinearObjFactory<panzer::Traits> & lo_factory,
                                            const Teuchos::RCP<const std::map<std::string, Teuchos::RCP<panzer::AuxiliaryEvaluator_TemplateManager<panzer::Traits> > > > 
                                               & auxManager)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  dofMngr_ = dofManager; // assign DOF manager for now

  worksets_.clear();
  phx_volume_field_managers_.clear();

  std::map<std::string, Teuchos::RCP<panzer::AuxiliaryEvaluator_TemplateManager<panzer::Traits> > >::const_iterator auxItr;
  
  std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator blkItr;
  for (blkItr=physicsBlocks.begin();blkItr!=physicsBlocks.end();++blkItr) {
    RCP<panzer::PhysicsBlock> pb = *blkItr;
    std::string blockId = pb->elementBlockID();

    // build a field manager object
    Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm 
          = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>);
    
    // use the physics block to register evaluators
    pb->buildAndRegisterEquationSetEvaluators(*fm);
    pb->buildAndRegisterGatherScatterEvaluators(*fm,lo_factory);
    //pb->buildAndRegisterModelEvaluators(fm, pb->getProvidedDOFs());

    if(auxManager!=Teuchos::null) {
       // find the set of auxilary factories for this block
       auxItr = auxManager->find(blockId);

       // add any auxiliary evaluators needed
       if(auxItr!=auxManager->end()) {
          typedef panzer::AuxiliaryEvaluator_TemplateManager<panzer::Traits> AuxEval_TM_Type;
          Teuchos::RCP<const AuxEval_TM_Type> auxEvalTM = auxItr->second; 

          // loop over all evaluator types stored in the tempalte manaager
          TEUCHOS_ASSERT(auxEvalTM!=Teuchos::null);
          std::cout << "Adding auxiliary evaluators in \"" << blockId << "\"" << std::endl;
          for(AuxEval_TM_Type::const_iterator auxEvalItr=auxEvalTM->begin();
                auxEvalItr!=auxEvalTM->end();++auxEvalItr) {
             const std::vector<Teuchos::RCP<panzer::AuxiliaryEvaluator_FactoryBase> > & auxFactory_vec = *auxEvalItr;

             // loop over vector of factories and have them build and register evaluators
             std::vector<Teuchos::RCP<panzer::AuxiliaryEvaluator_FactoryBase> >::const_iterator auxFactItr;
             for(auxFactItr=auxFactory_vec.begin();auxFactItr!=auxFactory_vec.end();++auxFactItr)
                (*auxFactItr)->buildAndRegisterEvaluators(*fm);
          }
       }
    }
     
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
                           const panzer::BCStrategyFactory& bc_factory,
                           const panzer::LinearObjFactory<panzer::Traits> & lo_factory)
{

  // for convenience build a map (element block id => physics block)
  std::map<std::string,Teuchos::RCP<panzer::PhysicsBlock> > physicsBlocks_map;
  {
     std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator blkItr;
     for(blkItr=physicsBlocks.begin();blkItr!=physicsBlocks.end();++blkItr) {
        RCP<panzer::PhysicsBlock> pb = *blkItr;
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
	bcs_type->buildAndRegisterEvaluators(fm,*side_pb);
	bcs_type->buildAndRegisterGatherScatterEvaluators(fm,*side_pb,lo_factory);
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
