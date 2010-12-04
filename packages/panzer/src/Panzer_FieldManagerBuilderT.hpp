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

typedef std::pair<std::string,Teuchos::RCP<panzer::Basis> > StrBasisPair;
struct StrBasisComp {
   bool operator() (const StrBasisPair & lhs, const StrBasisPair & rhs) const
   {return lhs.first<rhs.first;}
};

//=======================================================================
//=======================================================================
template<typename LO, typename GO>
void panzer::FieldManagerBuilder<LO,GO>::
setup(const Teuchos::RCP<panzer::ConnManager<LO,GO> >& conn_manager,
      MPI_Comm comm,
      const std::map<std::string,std::string>& block_ids_to_physics_ids,
      const std::map<std::string,panzer::InputPhysicsBlock>& physics_id_to_input_physics_blocks,
      const std::map<std::string,Teuchos::RCP<std::vector<panzer::Workset> > >& volume_worksets,
      const std::map<panzer::BC,Teuchos::RCP<std::map<unsigned,panzer::Workset> >,panzer::LessBC>& bc_worksets,
      int base_cell_dimension,
      const panzer::EquationSetFactory& eqset_factory,
      std::size_t workset_size)
{
  Teuchos::RCP<Teuchos::FancyOStream> pout = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  pout->setShowProcRank(true);
  pout->setOutputToRootOnly(0);

  // Build field managers for each element block
  phx_volume_field_managers_.resize(0);
  for (std::size_t i=0; i < block_ids_to_physics_ids.size(); ++i)
    phx_volume_field_managers_.push_back(Teuchos::rcp(new PHX::FieldManager<panzer::Traits>));
  
  // build the DOF manager for the problem
  dofMngr_ = Teuchos::rcp(new panzer::DOFManager<LO,GO>(conn_manager,comm));
  
  // Build the physics objects and register all variable providers
  std::vector<Teuchos::RCP<panzer::PhysicsBlock> > phxPhysicsBlocks;
  {
    using Teuchos::RCP;
    using Teuchos::rcp;

    // build a vector of physics blocks, get field names
    std::map<std::string,std::string>::const_iterator itr = block_ids_to_physics_ids.begin();
    for (;itr!=block_ids_to_physics_ids.end();++itr) {
      
      // build physics block
      const panzer::CellData volume_cell_data(workset_size, 
					      base_cell_dimension);
      
      
      std::map<std::string,panzer::InputPhysicsBlock>::const_iterator ipb_it = 
	physics_id_to_input_physics_blocks.find(itr->second);

      TEST_FOR_EXCEPTION(ipb_it == physics_id_to_input_physics_blocks.end(),
			 std::runtime_error,
			 "Falied to find InputPhysicsBlock for physics id: "
			 << itr->second << "!");

      const panzer::InputPhysicsBlock& ipb = ipb_it->second;

      RCP<panzer::PhysicsBlock> pb = 
	rcp(new panzer::PhysicsBlock(ipb, volume_cell_data, eqset_factory));
      
      phxPhysicsBlocks.push_back(pb);
      worksets_.push_back(volume_worksets.find(itr->first)->second);
      
      // save block's fields
      const std::vector<StrBasisPair> & blockFields = pb->getProvidedDOFs();

      // insert all fields into a set
      std::set<StrBasisPair,StrBasisComp> fieldNames;
      fieldNames.insert(blockFields.begin(),blockFields.end()); 

      // add basis to DOF manager: block specific
      std::set<StrBasisPair>::const_iterator fieldItr; 
      for (fieldItr=fieldNames.begin();fieldItr!=fieldNames.end();++fieldItr) {
         Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > intrepidBasis 
               = fieldItr->second->getIntrepidBasis();
         Teuchos::RCP<IntrepidFieldPattern> fp = Teuchos::rcp(new IntrepidFieldPattern(intrepidBasis));
         *pout << "\"" << fieldItr->first << "\" Field Pattern = \n";
         fp->print(*pout);
         dofMngr_->addField(itr->first,fieldItr->first,fp);
      }

    }

#ifdef PANZER_DEBUG_ROGER
    for (std::size_t block=0; block < worksets_.size(); ++block) {
      *pout << "Block Index = " << block << std::endl;
      *pout << "  Number of Worksets = " << worksets_[block]->size() 
		  << std::endl;
      for (std::size_t i=0; i < worksets_[block]->size(); ++i) {
	*pout << "  Workset[" << i << "] size = " 
		    << (*worksets_[block])[i].num_cells << endl;
	for (THashList::iterator cell = (*worksets_[block])[i].begin; 
	     cell != (*worksets_[block])[i].end; ++ cell) {
	  *pout << "    " << cell.element()->globalIndex() << endl;
	}
      }
    }
#endif
    
   

    // build field managers from physics blocks
    for (std::size_t block=0; block < phxPhysicsBlocks.size(); ++block) {
      RCP<panzer::PhysicsBlock> pb = phxPhysicsBlocks[block];

      // register the evaluators
      PHX::FieldManager<panzer::Traits>& fm = 
	*(phx_volume_field_managers_[block]);

      pb->buildAndRegisterEquationSetEvaluators(fm);
      pb->buildAndRegisterGatherScatterEvaluators(fm);
      //pb->buildAndRegisterModelEvaluators(fm, pb->getProvidedDOFs());

      Traits::SetupData setupData;
      setupData.dofManager_ = dofMngr_;
      setupData.worksets_ = worksets_[block];
      fm.postRegistrationSetup(setupData);

      int my_rank;
      MPI_Comm_rank(comm, &my_rank);
      if (my_rank == 0) {
	std::stringstream filename;
	filename << "panzer_dependency_graph_" << block;  
	fm.writeGraphvizFile(filename.str(), ".dot");
      }
    }
  }
  
  dofMngr_->buildGlobalUnknowns();
  dofMngr_->printFieldInformation(*pout);

  /*

  // ***************************
  // BCs
  // ***************************
  typedef std::vector<panzer::BoundaryCondition>::const_iterator bc_iter;
  for (bc_iter bc = bcs.begin(); bc != bcs.end(); ++bc) {

    SBC_Set* sideset = mesh.SBCbyId(bc->sidesetID());
    
    std::map<int,int>::const_iterator physics_id_iterator = 
      block_ids_to_physics_ids.find(bc->elementBlockID());
    
    TEST_FOR_EXCEPTION(physics_id_iterator == block_ids_to_physics_ids.end(), 
		       std::runtime_error,
		       "Error: the following boundary condition has invalid "
		       << "element block:\n" << *bc << endl);

    int physics_id = physics_id_iterator->second;

    const panzer::InputPhysicsBlock& ipb = 
      panzer::getEntry(physics_id_to_input_physics_blocks, physics_id);

    const panzer::CellData volume_cell_data(workset_size, base_cell_dimension);

    Teuchos::RCP<panzer::PhysicsBlock> bc_pb = 
      Teuchos::rcp(new panzer::PhysicsBlock(bc->elementBlockID(), ipb,
					    coordinate_index,
					    volume_cell_data));

    const shards::CellTopology volume_cell_topology = 
      bc_pb->getBaseCellTopology();
    
    bc_worksets_[*bc] = 
      panzer::buildBCWorkset(sideset, *bc, coordinate_index_, ipb, 
			     volume_cell_topology);
    
    // Build one FieldManager for each local side workset for each dirichlet bc
    std::map<unsigned,PHX::FieldManager<panzer::Traits> >& field_managers = 
      bc_field_managers_[*bc];

    // Loop over local face indices and setup each field manager
    for (std::map<unsigned,panzer::Workset>::const_iterator wkst = 
	   bc_worksets_[*bc]->begin(); wkst != bc_worksets_[*bc]->end();
	 ++wkst) {

      PHX::FieldManager<panzer::Traits>& fm = field_managers[wkst->first];
      
      // register evaluators from strategy
      panzer::BCStrategyFactory bcs_factory;
      Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs = 
	bcs_factory.buildBCStrategy(*bc);
      
      const panzer::CellData side_cell_data(wkst->second.num_cells,
					    base_cell_dimension,
					    wkst->first);      

      panzer::PhysicsBlock side_pb(bc->elementBlockID(), ipb,
				   coordinate_index,
				   side_cell_data);
      
      // Iterate over evaluation types
      for (panzer::BCStrategy_TemplateManager<panzer::Traits>::iterator 
	     bcs_type = bcs->begin(); bcs_type != bcs->end(); ++bcs_type) {
	bcs_type->buildAndRegisterEvaluators(fm,side_pb);
      }

      // Setup the fieldmanager
      Traits::SetupData setupData;
      setupData.dofManager_ = dofMngr_;
      Teuchos::RCP<std::vector<panzer::Workset> > worksets = 
	Teuchos::rcp(new(std::vector<panzer::Workset>));
      worksets->push_back(wkst->second);
      setupData.worksets_ = worksets;
      fm.postRegistrationSetup(setupData);
    }

    int my_rank;
    MPI_Comm_rank(NEVADA::comm.getComm(), &my_rank);
    // ROGER FIXME: This seg faults on "navierstokes_phx"
    // if (my_rank == 0) {
    //   std::stringstream filename;
    //   filename << "panzer_dependency_graph_bc_" << std::endl;  
    //   std::cout << "writing dependcy graph" << std::endl;
    //   field_managers.begin()->second.writeGraphvizFile(filename.str(), ".dot");
    // }
   
    */ 


    


}

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
std::ostream& panzer::operator<<(std::ostream& os, const panzer::FieldManagerBuilder<LO,GO>& rfd)
{
  rfd.print(os);
  return os;
}
