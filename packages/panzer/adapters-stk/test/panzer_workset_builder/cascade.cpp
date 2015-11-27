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

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_Workset_Builder.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_BC.hpp"

#include "user_app_EquationSetFactory.hpp"

namespace panzer {

  void getNodeIds(stk_classic::mesh::EntityRank nodeRank,const stk_classic::mesh::Entity * element,
		  std::vector<stk_classic::mesh::EntityId> & nodeIds);

  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs);

  TEUCHOS_UNIT_TEST(cascade, side_element_cascade)
  {
    using Teuchos::RCP;

  #ifdef HAVE_MPI
     Teuchos::RCP<Teuchos::Comm<int> > tcomm = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));
  #else
     Teuchos::RCP<Teuchos::Comm<int> > tcomm = Teuchos::rcp(new Teuchos::SerialComm<int>);
  #endif
    
    // excercise subcell entities capability
    {
      RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
      pl->set("X Elements",4);
      pl->set("Y Elements",4);
      pl->set("X Procs",1);
      pl->set("Y Procs",2);
  
      panzer_stk_classic::SquareQuadMeshFactory factory;
      factory.setParameterList(pl);
      RCP<panzer_stk_classic::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
  
      std::vector<stk_classic::mesh::Entity*> sideEntities; 
      mesh->getMySides("left","eblock-0_0",sideEntities);

      std::vector<std::vector<stk_classic::mesh::Entity*> > subcells;
      panzer_stk_classic::workset_utils::getSubcellEntities(*mesh,sideEntities,subcells);

      if(tcomm->getRank()==0) {
        TEST_EQUALITY(subcells.size(),1);
        TEST_EQUALITY(subcells[0].size(),3);
      }
      else {
        TEST_EQUALITY(subcells.size(),1);
        TEST_EQUALITY(subcells[0].size(),3);
      }
    }

    {
      RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
      pl->set("X Elements",4);
      pl->set("Y Elements",4);
      pl->set("X Procs",1);
      pl->set("Y Procs",2);
  
      panzer_stk_classic::SquareQuadMeshFactory factory;
      factory.setParameterList(pl);
      RCP<panzer_stk_classic::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

      std::vector<stk_classic::mesh::Entity*> sideEntities; 
      mesh->getMySides("left","eblock-0_0",sideEntities);
      TEST_ASSERT(sideEntities.size()==2);

      std::vector<std::size_t> localSubcellDim,localSubcellIds;
      std::vector<stk_classic::mesh::Entity*> elements;

      panzer_stk_classic::workset_utils::getSideElementCascade(*mesh,"eblock-0_0",sideEntities,
                                                       localSubcellDim,localSubcellIds,elements);

      TEST_EQUALITY(localSubcellDim.size(),6);
      TEST_ASSERT(localSubcellDim.size()==localSubcellIds.size());
      TEST_ASSERT(localSubcellDim.size()==elements.size());
 
      int node_cnt = 0;
      int edge_cnt = 0;
      for(std::size_t i=0;i<localSubcellDim.size();i++) {
        node_cnt += (localSubcellDim[i]==0) ? 1 : 0;
        edge_cnt += (localSubcellDim[i]==1) ? 1 : 0;
      }

      TEST_EQUALITY(node_cnt,4);
      TEST_EQUALITY(edge_cnt,2);

      for(std::size_t i=0;i<localSubcellDim.size();i++)
        out << elements[i] << " " << localSubcellDim[i] << " " << localSubcellIds[i] << std::endl;
    }
  }

  void getNodeIds(stk_classic::mesh::EntityRank nodeRank,const stk_classic::mesh::Entity * element,
		  std::vector<stk_classic::mesh::EntityId> & nodeIds)
  {
    stk_classic::mesh::PairIterRelation nodeRel = element->relations(nodeRank);
    
    stk_classic::mesh::PairIterRelation::iterator itr;
    for(itr=nodeRel.begin();itr!=nodeRel.end();++itr) 
      nodeIds.push_back(itr->entity()->identifier());
  }
  
  void testInitialzation(const Teuchos::RCP<Teuchos::ParameterList>& ipb,
			 std::vector<panzer::BC>& bcs)
  {
    // Physics block
    Teuchos::ParameterList& physics_block = ipb->sublist("test physics");
    {
      Teuchos::ParameterList& p = physics_block.sublist("a");
      p.set("Type","Energy");
      p.set("Prefix","");
      p.set("Model ID","solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",2);
      p.set("Integration Order",1);
    }
    {
      Teuchos::ParameterList& p = physics_block.sublist("b");
      p.set("Type","Energy");
      p.set("Prefix","ION_");
      p.set("Model ID","ion solid");
      p.set("Basis Type","HGrad");
      p.set("Basis Order",1);
      p.set("Integration Order",1);
    }

    {
      std::size_t bc_id = 0;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "left";
      std::string element_block_id = "eblock-0_0";
      std::string dof_name = "UX";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }    
    {
      std::size_t bc_id = 0;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "right";
      std::string element_block_id = "eblock-1_0";
      std::string dof_name = "UX";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }   
    {
      std::size_t bc_id = 0;
      panzer::BCType neumann = BCT_Dirichlet;
      std::string sideset_id = "top";
      std::string element_block_id = "eblock-1_0";
      std::string dof_name = "UX";
      std::string strategy = "Constant";
      double value = 5.0;
      Teuchos::ParameterList p;
      p.set("Value",value);
      panzer::BC bc(bc_id, neumann, sideset_id, element_block_id, dof_name, 
		    strategy, p);
      bcs.push_back(bc);
    }
  }

}
