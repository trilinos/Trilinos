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

#include <Panzer_SetupPartitionedWorksetUtilities.hpp>
#include "Panzer_STK_WorksetFactory.hpp"

#include "Panzer_LocalMeshInfo.hpp"

#include "Panzer_WorksetFactoryBase.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_LocalMeshUtilities.hpp"



namespace panzer_stk {

/** Set mesh
  */
void WorksetFactory::setMesh(const Teuchos::RCP<const panzer_stk::STK_Interface> & mesh)
{
   mesh_ = mesh;
}

/** Build sets of boundary condition worksets
  */
Teuchos::RCP<std::map<unsigned,panzer::Workset> > WorksetFactory::
getSideWorksets(const panzer::BC & bc,
              const panzer::PhysicsBlock & pb) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  panzer::WorksetNeeds needs;
  needs.cellData = pb.cellData();

  const std::map<int,RCP<panzer::IntegrationRule> >& int_rules = pb.getIntegrationRules();
  for(std::map<int,RCP<panzer::IntegrationRule> >::const_iterator ir_itr = int_rules.begin();
      ir_itr != int_rules.end(); ++ir_itr)
    needs.int_rules.push_back(ir_itr->second);

  const std::map<std::string,Teuchos::RCP<panzer::PureBasis> >& bases= pb.getBases();
  for(std::map<std::string,Teuchos::RCP<panzer::PureBasis> >::const_iterator b_itr = bases.begin();
      b_itr != bases.end(); ++b_itr)
    needs.bases.push_back(b_itr->second);

  return getSideWorksets(bc,needs);
}

Teuchos::RCP<std::map<unsigned,panzer::Workset> > WorksetFactory::
getSideWorksets(const panzer::BC & bc,
              const panzer::WorksetNeeds & needs) const
{
  return panzer_stk::buildBCWorksets(*mesh_,needs,bc.elementBlockID(),bc.sidesetID());
}

Teuchos::RCP<std::map<unsigned,panzer::Workset> > WorksetFactory::
getSideWorksets(const panzer::BC & bc,
                const panzer::PhysicsBlock & pb_a,
                const panzer::PhysicsBlock & pb_b) const
{
  TEUCHOS_ASSERT(bc.bcType() == panzer::BCT_Interface);
  return panzer_stk::buildBCWorksets(*mesh_, pb_a, pb_b, bc.sidesetID());
}

Teuchos::RCP<std::vector<panzer::Workset> > WorksetFactory::
getWorksets(const panzer::WorksetDescriptor & worksetDesc,
            const panzer::PhysicsBlock & pb) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  panzer::WorksetNeeds needs;
  needs.cellData = pb.cellData();

  const std::map<int,RCP<panzer::IntegrationRule> >& int_rules = pb.getIntegrationRules();
  for(std::map<int,RCP<panzer::IntegrationRule> >::const_iterator ir_itr = int_rules.begin();
      ir_itr != int_rules.end(); ++ir_itr)
    needs.int_rules.push_back(ir_itr->second);

  const std::map<std::string,Teuchos::RCP<panzer::PureBasis> >& bases= pb.getBases();
  for(std::map<std::string,Teuchos::RCP<panzer::PureBasis> >::const_iterator b_itr = bases.begin();
      b_itr != bases.end(); ++b_itr)
    needs.bases.push_back(b_itr->second);

  return getWorksets(worksetDesc,needs);
}

Teuchos::RCP<std::vector<panzer::Workset> > WorksetFactory::
getWorksets(const panzer::WorksetDescriptor & worksetDesc,
            const panzer::WorksetNeeds & needs) const
{

  if(worksetDesc.requiresPartitioning()){
    if(mesh_info_ == Teuchos::null){
      auto mesh_info = Teuchos::rcp(new panzer::LocalMeshInfo<int,int>());
      panzer_stk::generateLocalMeshInfo(*mesh_,*mesh_info);
      mesh_info_ = mesh_info;
    }
    return panzer::buildPartitionedWorksets(*mesh_info_,worksetDesc,needs);
  } else if(!worksetDesc.useSideset()) {
    return panzer_stk::buildWorksets(*mesh_,worksetDesc.getElementBlock(), needs);
  }
  else if(worksetDesc.useSideset() && worksetDesc.sideAssembly()) {
    // uses cascade by default, each subcell has its own workset
    return panzer_stk::buildWorksets(*mesh_,needs,worksetDesc.getSideset(),worksetDesc.getElementBlock(),true);
  }
  else {
    TEUCHOS_ASSERT(false);
    
    // The following code does not yet function in full generality, we need
    // to fix how the assembly process is handled for sidesets 
    /*
    Teuchos::RCP<std::map<unsigned,panzer::Workset> > workset_map =
      panzer_stk::buildBCWorksets(*mesh_,pb,worksetDesc.getSideset());

    // loop over worksets, adding them to vector
    Teuchos::RCP<std::vector<panzer::Workset> > worksets = Teuchos::rcp(new std::vector<panzer::Workset>);
    for(std::map<unsigned,panzer::Workset>::const_iterator itr=workset_map->begin();
        itr!=workset_map->end();++itr)
      worksets->push_back(itr->second);

    return worksets;
    */
  }
}

}
