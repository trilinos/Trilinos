// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include "Panzer_FieldManagerBuilder.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Teuchos_FancyOStream.hpp"

#include "Shards_CellTopology.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_Shards_Utilities.hpp"
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_StlMap_Utilities.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_GlobalIndexer.hpp"

//#include "EpetraExt_BlockMapOut.h"

//=======================================================================
//=======================================================================
void panzer::FieldManagerBuilder::print(std::ostream& os) const
{
  os << "panzer::FieldManagerBuilder output:  Not implemented yet!";
}

//=======================================================================
panzer::FieldManagerBuilder::
FieldManagerBuilder(bool disablePhysicsBlockScatter,
                    bool disablePhysicsBlockGather)
  : disablePhysicsBlockScatter_(disablePhysicsBlockScatter)
  , disablePhysicsBlockGather_(disablePhysicsBlockGather)
  , active_evaluation_types_(Sacado::mpl::size<panzer::Traits::EvalTypes>::value, true)
{}

//=======================================================================
namespace {
  struct PostRegistrationFunctor {

    const std::vector<bool>& active_;
    PHX::FieldManager<panzer::Traits>& fm_;
    panzer::Traits::SD& setup_data_;

    PostRegistrationFunctor(const std::vector<bool>& active,
                            PHX::FieldManager<panzer::Traits>& fm,
                            panzer::Traits::SD& setup_data)
      : active_(active),fm_(fm),setup_data_(setup_data) {}

    template<typename T>
    void operator()(T) const {
      auto index = Sacado::mpl::find<panzer::Traits::EvalTypes,T>::value;
      if (active_[index])
        fm_.postRegistrationSetupForType<T>(setup_data_);
    }
  };
}

//=======================================================================
void panzer::FieldManagerBuilder::setupVolumeFieldManagers(
                                            const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                                            const std::vector<WorksetDescriptor> & wkstDesc,
					    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
					    const Teuchos::ParameterList& closure_models,
                                            const panzer::LinearObjFactory<panzer::Traits> & lo_factory,
					    const Teuchos::ParameterList& user_data,
                                            const GenericEvaluatorFactory & gEvalFact,
                                            bool closureModelByEBlock)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::FieldManagerBuilder::setupVolumeFieldManagers",setup_field_managers);

  TEUCHOS_TEST_FOR_EXCEPTION(getWorksetContainer()==Teuchos::null,std::logic_error,
                            "panzer::FMB::setupVolumeFieldManagers: method function getWorksetContainer() returns null. "
                            "Plase call setWorksetContainer() before calling this method");
  TEUCHOS_TEST_FOR_EXCEPTION(physicsBlocks.size()!=wkstDesc.size(),std::runtime_error,
                            "panzer::FMB::setupVolumeFieldManagers: physics block count must match workset descriptor count.");

  phx_volume_field_managers_.clear();

  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer = lo_factory.getRangeGlobalIndexer();

  for (std::size_t blkInd=0;blkInd<physicsBlocks.size();++blkInd) {
    RCP<panzer::PhysicsBlock> pb = physicsBlocks[blkInd];
    const WorksetDescriptor wd = wkstDesc[blkInd];

    Traits::SD setupData;
    {
      PANZER_FUNC_TIME_MONITOR_DIFF("getWorksets()",get_worksets);
      { setupData.worksets_ = getWorksetContainer()->getWorksets(wd); }
    }
    {
      PANZER_FUNC_TIME_MONITOR_DIFF("getOrientations()",get_orientations);
      { setupData.orientations_ = getWorksetContainer()->getOrientations(); }
    }
    if(setupData.worksets_->size()==0)
      continue;

    // sanity check
    TEUCHOS_ASSERT(wd.getElementBlock()==pb->elementBlockID());

    // build a field manager object
    Teuchos::RCP<PHX::FieldManager<panzer::Traits> > fm
          = Teuchos::rcp(new PHX::FieldManager<panzer::Traits>);

    // use the physics block to register active evaluators
    pb->setActiveEvaluationTypes(active_evaluation_types_);

    {
      PANZER_FUNC_TIME_MONITOR_DIFF("pb->buildAndRegisterEquationSetEvaluators()",build_and_reg_eq_set_eval);
      { pb->buildAndRegisterEquationSetEvaluators(*fm, user_data); }
    }

    if(!physicsBlockGatherDisabled()) {
      PANZER_FUNC_TIME_MONITOR_DIFF("pb->buildAndRegisterGatherAndOrientationEvaluators()",build_and_reg_gath_and_orient_eval);
      pb->buildAndRegisterGatherAndOrientationEvaluators(*fm,lo_factory,user_data);
    }

    {
      PANZER_FUNC_TIME_MONITOR_DIFF("pb->buildAndRegisterDOFProjectionsToIPEvaluators()",build_and_reg_dof_proj_eval);
      pb->buildAndRegisterDOFProjectionsToIPEvaluators(*fm,Teuchos::ptrFromRef(lo_factory),user_data);
    }

    if(!physicsBlockScatterDisabled()) {
      PANZER_FUNC_TIME_MONITOR_DIFF("pb->buildAndRegisterScatterEvaluators()",build_and_reg_scatter_eval);
      pb->buildAndRegisterScatterEvaluators(*fm,lo_factory,user_data);
    }

    if(closureModelByEBlock) {
      PANZER_FUNC_TIME_MONITOR_DIFF("pb->buildAndRegisterClosureModelEvaluators(): closureModelByEBlock==true",build_and_reg_closure_model_eval_if);
      pb->buildAndRegisterClosureModelEvaluators(*fm,cm_factory,pb->elementBlockID(),closure_models,user_data);
    }
    else {
      PANZER_FUNC_TIME_MONITOR_DIFF("pb->buildAndRegisterClosureModelEvaluators(): closureModelByEBlock==false",build_and_reg_closure_model_eval_else);
      pb->buildAndRegisterClosureModelEvaluators(*fm,cm_factory,closure_models,user_data);
    }

    // Reset active evaluation types
    pb->activateAllEvaluationTypes();

    // register additional model evaluator from the generic evaluator factory
    gEvalFact.registerEvaluators(*fm,wd,*pb);

    // setup derivative information
    setKokkosExtendedDataTypeDimensions(wd.getElementBlock(),*globalIndexer,user_data,*fm);

    // call postRegistrationSetup() for each active type
    Sacado::mpl::for_each_no_kokkos<panzer::Traits::EvalTypes>(PostRegistrationFunctor(active_evaluation_types_,*fm,setupData));

    // make sure to add the field manager & workset to the list
    volume_workset_desc_.push_back(wd);
    phx_volume_field_managers_.push_back(fm);
  }
}

//=======================================================================
void panzer::FieldManagerBuilder::setupVolumeFieldManagers(
                                            const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
					    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
					    const Teuchos::ParameterList& closure_models,
                                            const panzer::LinearObjFactory<panzer::Traits> & lo_factory,
					    const Teuchos::ParameterList& user_data)
{
   std::vector<WorksetDescriptor> wkstDesc;
   for(std::size_t i=0;i<physicsBlocks.size();i++)
     wkstDesc.push_back(blockDescriptor(physicsBlocks[i]->elementBlockID()));

   EmptyEvaluatorFactory eef;
   setupVolumeFieldManagers(physicsBlocks,wkstDesc,cm_factory,closure_models,lo_factory,user_data,eef);
}

//=======================================================================
//=======================================================================
void panzer::FieldManagerBuilder::
setupBCFieldManagers(const std::vector<panzer::BC> & bcs,
                     const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
                     const Teuchos::Ptr<const panzer::EquationSetFactory>& /* eqset_factory */,
                     const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
                     const panzer::BCStrategyFactory& bc_factory,
                     const Teuchos::ParameterList& closure_models,
                     const panzer::LinearObjFactory<panzer::Traits> & lo_factory,
                     const Teuchos::ParameterList& user_data)
{
  TEUCHOS_TEST_FOR_EXCEPTION(getWorksetContainer()==Teuchos::null,std::logic_error,
                            "panzer::FMB::setupBCFieldManagers: method function getWorksetContainer() returns null. "
                            "Plase call setWorksetContainer() before calling this method");

  Teuchos::RCP<const panzer::GlobalIndexer> globalIndexer = lo_factory.getRangeGlobalIndexer();

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
  std::vector<panzer::BC>::const_iterator bc;
  for (bc=bcs.begin(); bc != bcs.end(); ++bc) {
    WorksetDescriptor wd = panzer::bcDescriptor(*bc);
    const Teuchos::RCP<std::map<unsigned,panzer::Workset> >
      currentWkst = getWorksetContainer()->getSideWorksets(wd);
    if (currentWkst.is_null()) continue;

    BCType bc_type = bc->bcType();

    if (bc_type == BCT_Interface) {
      // Loop over local face indices and setup each field manager
      for (std::map<unsigned,panzer::Workset>::const_iterator wkst = currentWkst->begin();
           wkst != currentWkst->end(); ++wkst) {
        // Build one FieldManager for each local side workset for each bc
        std::map<unsigned,PHX::FieldManager<panzer::Traits> >& field_managers =
          bc_field_managers_[*bc];

        PHX::FieldManager<panzer::Traits>& fm = field_managers[wkst->first];

        int gid_count = 0;
        for (int block_id_index = 0; block_id_index < 2; ++block_id_index) {
          const std::string element_block_id = block_id_index == 0 ? bc->elementBlockID() : bc->elementBlockID2();

          std::map<std::string,Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator
            volume_pb_itr = physicsBlocks_map.find(element_block_id);

          TEUCHOS_TEST_FOR_EXCEPTION(volume_pb_itr == physicsBlocks_map.end(), std::logic_error,
            "panzer::FMB::setupBCFieldManagers: Cannot find physics block corresponding to element block \""
            << element_block_id << "\"");

          const Teuchos::RCP<const panzer::PhysicsBlock> volume_pb = physicsBlocks_map.find(element_block_id)->second;
          const Teuchos::RCP<const shards::CellTopology> volume_cell_topology = volume_pb->cellData().getCellTopology();

          // register evaluators from strategy
          const panzer::CellData side_cell_data(wkst->second.num_cells,
                                                wkst->second.details(block_id_index).subcell_index,
                                                volume_cell_topology);

          // Copy the physics block for side integrations
          Teuchos::RCP<panzer::PhysicsBlock> side_pb = volume_pb->copyWithCellData(side_cell_data);

          Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
            bcstm = bc_factory.buildBCStrategy(*bc, side_pb->globalData());

          // Iterate over evaluation types
          int i=0;
          for (panzer::BCStrategy_TemplateManager<panzer::Traits>::iterator
                 bcs_type = bcstm->begin(); bcs_type != bcstm->end(); ++bcs_type,++i) {
            if (active_evaluation_types_[i]) {
              bcs_type->setDetailsIndex(block_id_index);
              side_pb->setDetailsIndex(block_id_index);
              bcs_type->setup(*side_pb, user_data);
              bcs_type->buildAndRegisterEvaluators(fm, *side_pb, cm_factory, closure_models, user_data);
              bcs_type->buildAndRegisterGatherAndOrientationEvaluators(fm, *side_pb, lo_factory, user_data);
              if ( ! physicsBlockScatterDisabled())
                bcs_type->buildAndRegisterScatterEvaluators(fm, *side_pb, lo_factory, user_data);
            }
          }

          gid_count += globalIndexer->getElementBlockGIDCount(element_block_id);
        }

        { // Use gid_count to set up the derivative information.
          std::vector<PHX::index_size_type> derivative_dimensions;
          derivative_dimensions.push_back(gid_count);
          fm.setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);

          #ifdef Panzer_BUILD_HESSIAN_SUPPORT
            fm.setKokkosExtendedDataTypeDimensions<panzer::Traits::Hessian>(derivative_dimensions);
          #endif

          derivative_dimensions[0] = 1;
          if (user_data.isType<int>("Tangent Dimension"))
            derivative_dimensions[0] = user_data.get<int>("Tangent Dimension");
          fm.setKokkosExtendedDataTypeDimensions<panzer::Traits::Tangent>(derivative_dimensions);
        }

        // Set up the field manager
        Traits::SD setupData;
        Teuchos::RCP<std::vector<panzer::Workset> > worksets = Teuchos::rcp(new std::vector<panzer::Workset>);
        worksets->push_back(wkst->second);
        setupData.worksets_ = worksets;
        setupData.orientations_ = getWorksetContainer()->getOrientations();

        Sacado::mpl::for_each_no_kokkos<panzer::Traits::EvalTypes>(PostRegistrationFunctor(active_evaluation_types_,fm,setupData));

      }
    } else {
      const std::string element_block_id = bc->elementBlockID();

      std::map<std::string,Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator volume_pb_itr
	= physicsBlocks_map.find(element_block_id);

      TEUCHOS_TEST_FOR_EXCEPTION(volume_pb_itr==physicsBlocks_map.end(),std::logic_error,
				 "panzer::FMB::setupBCFieldManagers: Cannot find physics block corresponding to element block \"" << element_block_id << "\"");

      Teuchos::RCP<const panzer::PhysicsBlock> volume_pb = physicsBlocks_map.find(element_block_id)->second;
      Teuchos::RCP<const shards::CellTopology> volume_cell_topology = volume_pb->cellData().getCellTopology();

      // Build one FieldManager for each local side workset for each dirichlet bc
      std::map<unsigned,PHX::FieldManager<panzer::Traits> >& field_managers =
        bc_field_managers_[*bc];

      // Loop over local face indices and setup each field manager
      for (std::map<unsigned,panzer::Workset>::const_iterator wkst =
	     currentWkst->begin(); wkst != currentWkst->end();
	   ++wkst) {

        PHX::FieldManager<panzer::Traits>& fm = field_managers[wkst->first];

        // register evaluators from strategy
        const panzer::CellData side_cell_data(wkst->second.num_cells,
	                                      wkst->first,volume_cell_topology);

	// Copy the physics block for side integrations
	Teuchos::RCP<panzer::PhysicsBlock> side_pb = volume_pb->copyWithCellData(side_cell_data);

	Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcstm =
	  bc_factory.buildBCStrategy(*bc,side_pb->globalData());

	// Iterate over evaluation types
        int i=0;
	for (panzer::BCStrategy_TemplateManager<panzer::Traits>::iterator
	       bcs_type = bcstm->begin(); bcs_type != bcstm->end(); ++bcs_type,++i) {
          if (active_evaluation_types_[i]) {
            bcs_type->setup(*side_pb,user_data);
            bcs_type->buildAndRegisterEvaluators(fm,*side_pb,cm_factory,closure_models,user_data);
            bcs_type->buildAndRegisterGatherAndOrientationEvaluators(fm,*side_pb,lo_factory,user_data);
            if(!physicsBlockScatterDisabled())
              bcs_type->buildAndRegisterScatterEvaluators(fm,*side_pb,lo_factory,user_data);
          }
	}

	// Setup the fieldmanager
	Traits::SD setupData;
	Teuchos::RCP<std::vector<panzer::Workset> > worksets =
	  Teuchos::rcp(new(std::vector<panzer::Workset>));
	worksets->push_back(wkst->second);
	setupData.worksets_ = worksets;
        setupData.orientations_ = getWorksetContainer()->getOrientations();

	// setup derivative information
	setKokkosExtendedDataTypeDimensions(element_block_id,*globalIndexer,user_data,fm);

        Sacado::mpl::for_each_no_kokkos<panzer::Traits::EvalTypes>(PostRegistrationFunctor(active_evaluation_types_,fm,setupData));
      }
    }
  }
}

//=======================================================================
//=======================================================================
void panzer::FieldManagerBuilder::
writeVolumeGraphvizDependencyFiles(std::string filename_prefix,
				   const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks) const
{
  if(phx_volume_field_managers_.size()<1)
    return; // nothing to see here folks

  TEUCHOS_ASSERT(phx_volume_field_managers_.size()==physicsBlocks.size());

  std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator blkItr;
  int index = 0;
  for (blkItr=physicsBlocks.begin();blkItr!=physicsBlocks.end();++blkItr,++index) {
    std::string blockId = (*blkItr)->elementBlockID();
    phx_volume_field_managers_[index]->writeGraphvizFile(filename_prefix+"_VOLUME_"+blockId);
  }

}

//=======================================================================
//=======================================================================
void panzer::FieldManagerBuilder::
writeBCGraphvizDependencyFiles(std::string filename_prefix) const
{
  typedef std::map<panzer::BC,std::map<unsigned,PHX::FieldManager<panzer::Traits> >,panzer::LessBC> FMMap;

  FMMap::const_iterator blkItr;
  int bc_index = 0;
  for (blkItr=bc_field_managers_.begin();blkItr!=bc_field_managers_.end();++blkItr,++bc_index) {
    panzer::BC bc = blkItr->first;
    const PHX::FieldManager<panzer::Traits> & fm = blkItr->second.begin()->second; // get the first field manager

    BCType bc_type = bc.bcType();
    std::string type;
    if (bc_type == BCT_Dirichlet)
	type = "_Dirichlet_";
    else if (bc_type == BCT_Neumann)
        type = "_Neumann_";
    else if (bc_type == BCT_Interface)
        type = "_Interface_";
    else
        TEUCHOS_ASSERT(false);

    std::string blockId = bc.elementBlockID();
    std::string sideId = bc.sidesetID();
    fm.writeGraphvizFile(filename_prefix+"_BC_"+std::to_string(bc_index)+type+sideId+"_"+blockId);
  }

}

//=======================================================================
//=======================================================================
void panzer::FieldManagerBuilder::
writeVolumeTextDependencyFiles(std::string filename_prefix,
			       const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks) const
{
  if(phx_volume_field_managers_.size()<1)
    return; // nothing to see here folks

  TEUCHOS_ASSERT(phx_volume_field_managers_.size()==physicsBlocks.size());

  std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator blkItr;
  int index = 0;
  for (blkItr=physicsBlocks.begin();blkItr!=physicsBlocks.end();++blkItr,++index) {

    std::string blockId = (*blkItr)->elementBlockID();

    std::string filename = filename_prefix+"_VOLUME_"+blockId+".txt";
    std::ofstream ofs;
    ofs.open(filename.c_str());

    ofs << *(phx_volume_field_managers_[index]) << std::endl;

    ofs.close();
  }

}

//=======================================================================
//=======================================================================
void panzer::FieldManagerBuilder::
writeBCTextDependencyFiles(std::string filename_prefix) const
{
  typedef std::map<panzer::BC,std::map<unsigned,PHX::FieldManager<panzer::Traits> >,panzer::LessBC> FMMap;

  FMMap::const_iterator blkItr;
  int bc_index = 0;
  for (blkItr=bc_field_managers_.begin();blkItr!=bc_field_managers_.end();++blkItr,++bc_index) {
    panzer::BC bc = blkItr->first;
    const PHX::FieldManager<panzer::Traits> & fm = blkItr->second.begin()->second; // get the first field manager

    BCType bc_type = bc.bcType();
    std::string type;
    if (bc_type == BCT_Dirichlet)
	type = "_Dirichlet_";
    else if (bc_type == BCT_Neumann)
        type = "_Neumann_";
    else if (bc_type == BCT_Interface)
        type = "_Interface_";
    else
        TEUCHOS_ASSERT(false);

    std::string blockId = bc.elementBlockID();
    std::string sideId = bc.sidesetID();

    std::string filename = filename_prefix+"_BC_"+std::to_string(bc_index)+type+sideId+"_"+blockId+".txt";
    std::ofstream ofs;
    ofs.open(filename.c_str());

    ofs << fm << std::endl;

    ofs.close();
  }

}

//=======================================================================
//=======================================================================
void panzer::FieldManagerBuilder::
setKokkosExtendedDataTypeDimensions(const std::string & eblock,
                                    const panzer::GlobalIndexer & globalIndexer,
                                    const Teuchos::ParameterList& user_data,
                                    PHX::FieldManager<panzer::Traits> & fm) const
{
  // setup Jacobian derivative terms
  {
    std::vector<PHX::index_size_type> derivative_dimensions;
    derivative_dimensions.push_back(globalIndexer.getElementBlockGIDCount(eblock));

    fm.setKokkosExtendedDataTypeDimensions<panzer::Traits::Jacobian>(derivative_dimensions);

  }

  #ifdef Panzer_BUILD_HESSIAN_SUPPORT
  {
    std::vector<PHX::index_size_type> derivative_dimensions;
    derivative_dimensions.push_back(globalIndexer.getElementBlockGIDCount(eblock));

    fm.setKokkosExtendedDataTypeDimensions<panzer::Traits::Hessian>(derivative_dimensions);
  }
  #endif

  {
    std::vector<PHX::index_size_type> derivative_dimensions;
    derivative_dimensions.push_back(1);
    if (user_data.isType<int>("Tangent Dimension"))
      derivative_dimensions[0] = user_data.get<int>("Tangent Dimension");
    fm.setKokkosExtendedDataTypeDimensions<panzer::Traits::Tangent>(derivative_dimensions);
  }
}

void panzer::FieldManagerBuilder::setActiveEvaluationTypes(const std::vector<bool>& aet)
{active_evaluation_types_ = aet;}

//=======================================================================
//=======================================================================
void panzer::FieldManagerBuilder::clearVolumeFieldManagers(bool clearVolumeWorksets)
{
  phx_volume_field_managers_.clear();
  volume_workset_desc_.clear();
  if (clearVolumeWorksets)
    worksetContainer_->clearVolumeWorksets();
}

//=======================================================================
//=======================================================================
std::ostream& panzer::operator<<(std::ostream& os, const panzer::FieldManagerBuilder& rfd)
{
  rfd.print(os);
  return os;
}
