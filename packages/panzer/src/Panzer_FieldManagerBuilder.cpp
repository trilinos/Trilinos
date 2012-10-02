// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include <vector>
#include <string>
#include <sstream>

#include "Panzer_FieldManagerBuilder.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Teuchos_FancyOStream.hpp"

#include "Shards_CellTopology.hpp"

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
#include "Panzer_InputPhysicsBlock.hpp"
#include "Panzer_StlMap_Utilities.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"

//#include "EpetraExt_BlockMapOut.h"

//=======================================================================
//=======================================================================
void panzer::FieldManagerBuilder::print(std::ostream& os) const
{
  os << "panzer::FieldManagerBuilder output:  Not implemented yet!";
}

//=======================================================================
//=======================================================================
void panzer::FieldManagerBuilder::setupVolumeFieldManagers(
                                            const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks, 
					    const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
					    const Teuchos::ParameterList& closure_models,
                                            const panzer::LinearObjFactory<panzer::Traits> & lo_factory,
					    const Teuchos::ParameterList& user_data)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  TEUCHOS_TEST_FOR_EXCEPTION(getWorksetContainer()==Teuchos::null,std::logic_error,
                            "panzer::FMB::setupVolumeFieldManagers: method function getWorksetContainer() returns null. "
                            "Plase call setWorksetContainer() before calling this method");

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
    pb->buildAndRegisterGatherAndOrientationEvaluators(*fm,lo_factory,user_data);
    pb->buildAndRegisterDOFProjectionsToIPEvaluators(*fm,user_data);
    if(!physicsBlockScatterDisabled())
      pb->buildAndRegisterScatterEvaluators(*fm,lo_factory,user_data);
    pb->buildAndRegisterClosureModelEvaluators(*fm,cm_factory,closure_models,user_data);

    // build the setup data using passed in information
    Traits::SetupData setupData;
    setupData.worksets_ = getWorksetContainer()->getVolumeWorksets(blockId);
    fm->postRegistrationSetup(setupData);

    // make sure to add the field manager & workset to the list 
    element_block_names_.push_back(blockId);
    phx_volume_field_managers_.push_back(fm); 
  }
}

//=======================================================================
//=======================================================================
void panzer::FieldManagerBuilder::
setupBCFieldManagers(const std::vector<panzer::BC> & bcs,
                     const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
	             const panzer::EquationSetFactory & eqset_factory,
                     const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
                     const panzer::BCStrategyFactory& bc_factory,
                     const Teuchos::ParameterList& closure_models,
                     const panzer::LinearObjFactory<panzer::Traits> & lo_factory,
                     const Teuchos::ParameterList& user_data)
{
  TEUCHOS_TEST_FOR_EXCEPTION(getWorksetContainer()==Teuchos::null,std::logic_error,
                            "panzer::FMB::setupBCFieldManagers: method function getWorksetContainer() returns null. "
                            "Plase call setWorksetContainer() before calling this method");

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
    std::string element_block_id = bc->elementBlockID(); 
    Teuchos::RCP<const panzer::PhysicsBlock> volume_pb = physicsBlocks_map.find(element_block_id)->second;
    Teuchos::RCP<const shards::CellTopology> volume_cell_topology = volume_pb->cellData().getCellTopology();
    int base_cell_dimension = volume_pb->cellData().baseCellDimension();
    
    Teuchos::RCP<std::map<unsigned,panzer::Workset> > currentWkst = getWorksetContainer()->getSideWorksets(*bc);
    if(currentWkst==Teuchos::null) // if there is nothing to do...do nothing!
       continue;

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
					    base_cell_dimension,
					    wkst->first,volume_cell_topology);      

      Teuchos::RCP<panzer::PhysicsBlock> side_pb 
            = volume_pb->copyWithCellData(side_cell_data, eqset_factory);

      Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> > bcs = 
	bc_factory.buildBCStrategy(*bc,side_pb->globalData());

      // Iterate over evaluation types
      for (panzer::BCStrategy_TemplateManager<panzer::Traits>::iterator 
	     bcs_type = bcs->begin(); bcs_type != bcs->end(); ++bcs_type) {
	bcs_type->setup(*side_pb,user_data);
	bcs_type->buildAndRegisterEvaluators(fm,*side_pb,cm_factory,closure_models,user_data);
	bcs_type->buildAndRegisterGatherAndOrientationEvaluators(fm,*side_pb,lo_factory,user_data);
	bcs_type->buildAndRegisterScatterEvaluators(fm,*side_pb,lo_factory,user_data);
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
void panzer::FieldManagerBuilder::
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
std::ostream& panzer::operator<<(std::ostream& os, const panzer::FieldManagerBuilder& rfd)
{
  rfd.print(os);
  return os;
}
