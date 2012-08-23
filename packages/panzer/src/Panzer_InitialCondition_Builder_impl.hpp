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

#ifndef PANZER_INITIAL_CONDITION_BUILDER_IMPL_HPP
#define PANZER_INITIAL_CONDITION_BUILDER_IMPL_HPP

#include "Teuchos_Assert.hpp"

void panzer::setupInitialConditionFieldManagers(WorksetContainer & wkstContainer,
						const std::vector<Teuchos::RCP<panzer::PhysicsBlock> >& physicsBlocks,
						const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& cm_factory,
						const Teuchos::ParameterList& ic_block_closure_models,
						const panzer::LinearObjFactory<panzer::Traits>& lo_factory,
						const Teuchos::ParameterList& user_data,
						const bool write_graphviz_file,
						const std::string& graphviz_file_prefix,
						std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >& phx_ic_field_managers)
{
  std::vector<Teuchos::RCP<panzer::PhysicsBlock> >::const_iterator blkItr;
  for (blkItr=physicsBlocks.begin();blkItr!=physicsBlocks.end();++blkItr) {
    Teuchos::RCP<panzer::PhysicsBlock> pb = *blkItr;
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
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Failed to find initial condition for element block \"" << blockId << "\".  You must provide an initial condition for each element block or set a default!" << ic_block_closure_models);

    // use the physics block to register evaluators
    pb->buildAndRegisterInitialConditionEvaluators(*fm, cm_factory, closure_model_name, ic_block_closure_models, lo_factory, user_data);

    // build the setup data using passed in information
    Traits::SetupData setupData;
    setupData.worksets_ = wkstContainer.getVolumeWorksets(blockId);

    fm->postRegistrationSetup(setupData);
    phx_ic_field_managers.push_back(fm);
    
    if (write_graphviz_file)
      fm->writeGraphvizFile(graphviz_file_prefix+"IC_"+blockId);
  }
}

void panzer::evaluateInitialCondition(const std::vector< Teuchos::RCP<std::vector<panzer::Workset> > >& worksets,
				      const std::vector< Teuchos::RCP< PHX::FieldManager<panzer::Traits> > >& phx_ic_field_managers,
				      Teuchos::RCP<panzer::LinearObjContainer> loc,
				      const double time_stamp)
{   
  for (std::size_t block = 0; block < worksets.size(); ++block) {
    
    std::vector<panzer::Workset>& w = *worksets[block]; 
    
    Teuchos::RCP< PHX::FieldManager<panzer::Traits> > fm = 
      phx_ic_field_managers[block];

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
