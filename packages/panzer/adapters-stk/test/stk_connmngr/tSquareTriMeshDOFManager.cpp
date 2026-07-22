// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;

namespace panzer_stk {

Teuchos::RCP<panzer::ConnManager> buildTriMesh(stk::ParallelMachine comm,int xelmts,int yelmts,int xblocks,int yblocks)
{
   Teuchos::ParameterList pl;
   pl.set<int>("X Elements",xelmts);
   pl.set<int>("Y Elements",yelmts);
   pl.set<int>("X Blocks",xblocks);
   pl.set<int>("Y Blocks",yblocks);

   panzer_stk::SquareTriMeshFactory meshFact;
   meshFact.setParameterList(Teuchos::rcpFromRef(pl));

   Teuchos::RCP<panzer_stk::STK_Interface> mesh = meshFact.buildMesh(comm);
   return Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
}

template <typename Intrepid2Type>
RCP<const panzer::Intrepid2FieldPattern> buildFieldPattern()
{
   // build a geometric pattern from a single basis
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2Type);
   RCP<const panzer::Intrepid2FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
   return pattern;
}

// quad tests
TEUCHOS_UNIT_TEST(tSquareTriMeshDOFManager, buildTest_tri)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      stk::ParallelMachine Comm = MPI_COMM_WORLD;
   #else
      stk::ParallelMachine Comm = WHAT_TO_DO_COMM;
   #endif

   int numProcs = stk::parallel_machine_size(Comm);
   int myRank = stk::parallel_machine_rank(Comm);

   TEUCHOS_ASSERT(numProcs==2);

   // build a geometric pattern from a single basis
   RCP<const panzer::FieldPattern> patternC1
         = buildFieldPattern<Intrepid2::Basis_HGRAD_TRI_C1_FEM<PHX::exec_space,double,double> >();

   RCP<panzer::ConnManager> connManager = buildTriMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager> dofManager = rcp(new panzer::DOFManager());

   TEST_EQUALITY(dofManager->getOrientationsRequired(),false);
   TEST_EQUALITY(dofManager->getConnManager(),Teuchos::null);

   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   TEST_EQUALITY(dofManager->getConnManager(),connManager);

   dofManager->addField("ux",patternC1);
   dofManager->addField("uy",patternC1);
   dofManager->addField("p",patternC1);

   std::vector<std::string> fieldOrder;
   fieldOrder.push_back("p");
   fieldOrder.push_back("ux");
   fieldOrder.push_back("uy");
   dofManager->setFieldOrder(fieldOrder);

   dofManager->buildGlobalUnknowns();
   dofManager->printFieldInformation(out);

   TEST_EQUALITY(dofManager->getFieldNum("p"),0);
   TEST_EQUALITY(dofManager->getFieldNum("ux"),1);
   TEST_EQUALITY(dofManager->getFieldNum("uy"),2);

   const std::vector<int> & uy_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("uy"));
   const std::vector<int> & p_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("p"));
   const std::vector<int> & ux_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("ux"));

   TEST_EQUALITY(uy_offsets.size(),p_offsets.size());
   TEST_EQUALITY(uy_offsets.size(),ux_offsets.size());

   if(myRank==0) {
      std::vector<panzer::GlobalOrdinal> gids;

      dofManager->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),9);
      TEST_EQUALITY(gids[0],0); TEST_EQUALITY(gids[1],1); TEST_EQUALITY(gids[2],2);
      TEST_EQUALITY(gids[3],3); TEST_EQUALITY(gids[4],4); TEST_EQUALITY(gids[5],5);
      TEST_EQUALITY(gids[6],9); TEST_EQUALITY(gids[7],10); TEST_EQUALITY(gids[8],11);

      for(std::size_t i=0;i<p_offsets.size();i++) {
         TEST_ASSERT(gids[p_offsets[i]]<gids[ux_offsets[i]]);
         TEST_ASSERT(gids[p_offsets[i]]<gids[uy_offsets[i]]);
         TEST_ASSERT(gids[ux_offsets[i]]<gids[uy_offsets[i]]);
      }

      dofManager->getElementGIDs(1,gids);
      TEST_EQUALITY(gids.size(),9);
      TEST_EQUALITY(gids[0],0); TEST_EQUALITY(gids[1],1); TEST_EQUALITY(gids[2],2);
      TEST_EQUALITY(gids[3],9); TEST_EQUALITY(gids[4],10); TEST_EQUALITY(gids[5],11);
      TEST_EQUALITY(gids[6],6); TEST_EQUALITY(gids[7],7); TEST_EQUALITY(gids[8],8);

      for(std::size_t i=0;i<p_offsets.size();i++) {
         TEST_ASSERT(gids[p_offsets[i]]<gids[ux_offsets[i]]);
         TEST_ASSERT(gids[p_offsets[i]]<gids[uy_offsets[i]]);
         TEST_ASSERT(gids[ux_offsets[i]]<gids[uy_offsets[i]]);
      }

      dofManager->getElementGIDs(3,gids);
      TEST_EQUALITY(gids.size(),9);
      TEST_EQUALITY(gids[0],6); TEST_EQUALITY(gids[1],7); TEST_EQUALITY(gids[2],8);
      TEST_EQUALITY(gids[3],15); TEST_EQUALITY(gids[4],16); TEST_EQUALITY(gids[5],17);
      TEST_EQUALITY(gids[6],12); TEST_EQUALITY(gids[7],13); TEST_EQUALITY(gids[8],14);

      for(std::size_t i=0;i<p_offsets.size();i++) {
         TEST_ASSERT(gids[p_offsets[i]]<gids[ux_offsets[i]]);
         TEST_ASSERT(gids[p_offsets[i]]<gids[uy_offsets[i]]);
         TEST_ASSERT(gids[ux_offsets[i]]<gids[uy_offsets[i]]);
      }
   }
   else if(myRank==1) {
      std::vector<panzer::GlobalOrdinal> gids;

      dofManager->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),9);
      TEST_EQUALITY(gids[0],3); TEST_EQUALITY(gids[1],4); TEST_EQUALITY(gids[2],5);
      TEST_EQUALITY(gids[3],18); TEST_EQUALITY(gids[4],19); TEST_EQUALITY(gids[5],20);
      TEST_EQUALITY(gids[6],21); TEST_EQUALITY(gids[7],22); TEST_EQUALITY(gids[8],23);

      for(std::size_t i=0;i<p_offsets.size();i++) {
         TEST_ASSERT(gids[p_offsets[i]]<gids[ux_offsets[i]]);
         TEST_ASSERT(gids[p_offsets[i]]<gids[uy_offsets[i]]);
         TEST_ASSERT(gids[ux_offsets[i]]<gids[uy_offsets[i]]);
      }

      dofManager->getElementGIDs(1,gids);
      TEST_EQUALITY(gids.size(),9);
      TEST_EQUALITY(gids[0],3); TEST_EQUALITY(gids[1],4); TEST_EQUALITY(gids[2],5);
      TEST_EQUALITY(gids[3],21); TEST_EQUALITY(gids[4],22); TEST_EQUALITY(gids[5],23);
      TEST_EQUALITY(gids[6],9); TEST_EQUALITY(gids[7],10); TEST_EQUALITY(gids[8],11);

      for(std::size_t i=0;i<p_offsets.size();i++) {
         TEST_ASSERT(gids[p_offsets[i]]<gids[ux_offsets[i]]);
         TEST_ASSERT(gids[p_offsets[i]]<gids[uy_offsets[i]]);
         TEST_ASSERT(gids[ux_offsets[i]]<gids[uy_offsets[i]]);
      }

      dofManager->getElementGIDs(3,gids);
      TEST_EQUALITY(gids.size(),9);
      TEST_EQUALITY(gids[0],9); TEST_EQUALITY(gids[1],10); TEST_EQUALITY(gids[2],11);
      TEST_EQUALITY(gids[3],24); TEST_EQUALITY(gids[4],25); TEST_EQUALITY(gids[5],26);
      TEST_EQUALITY(gids[6],15); TEST_EQUALITY(gids[7],16); TEST_EQUALITY(gids[8],17);

      for(std::size_t i=0;i<p_offsets.size();i++) {
         TEST_ASSERT(gids[p_offsets[i]]<gids[ux_offsets[i]]);
         TEST_ASSERT(gids[p_offsets[i]]<gids[uy_offsets[i]]);
         TEST_ASSERT(gids[ux_offsets[i]]<gids[uy_offsets[i]]);
      }
   }
}

/*
TEUCHOS_UNIT_TEST(tSquareTriMeshDOFManager, field_order)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      stk::ParallelMachine Comm = MPI_COMM_WORLD;
   #else
      stk::ParallelMachine Comm = WHAT_TO_DO_COMM;
   #endif

   int numProcs = stk::parallel_machine_size(Comm);
   int myRank = stk::parallel_machine_rank(Comm);

   TEUCHOS_ASSERT(numProcs==2);

   // build a geometric pattern from a single basis
   RCP<const panzer::FieldPattern> patternC1
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();

   RCP<panzer::ConnManager> connManager = buildTriMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager> dofManager = rcp(new panzer::DOFManager());

   TEST_EQUALITY(dofManager->getConnManager(),Teuchos::null);

   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   TEST_EQUALITY(dofManager->getConnManager(),connManager);

   dofManager->addField("ux",patternC1);
   dofManager->addField("uy",patternC1);
   dofManager->addField("p",patternC1);
   std::vector<std::string> fieldOrder;
   fieldOrder.push_back("uy");
   fieldOrder.push_back("p");
   fieldOrder.push_back("ux");

   dofManager->setFieldOrder(fieldOrder);

   dofManager->buildGlobalUnknowns();
   dofManager->printFieldInformation(out);

   TEST_EQUALITY(dofManager->getFieldNum("uy"),0);
   TEST_EQUALITY(dofManager->getFieldNum("p"),1);
   TEST_EQUALITY(dofManager->getFieldNum("ux"),2);

   const std::vector<int> & uy_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("uy"));
   const std::vector<int> & p_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("p"));
   const std::vector<int> & ux_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("ux"));

   TEST_EQUALITY(uy_offsets.size(),p_offsets.size());
   TEST_EQUALITY(uy_offsets.size(),ux_offsets.size());

   if(myRank==0) {
      std::vector<int> gids;

      dofManager->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),12);
      for(std::size_t i=0;i<uy_offsets.size();i++) {
         TEST_ASSERT(gids[uy_offsets[i]]<gids[p_offsets[i]]);
         TEST_ASSERT(gids[uy_offsets[i]]<gids[ux_offsets[i]]);
         TEST_ASSERT(gids[p_offsets[i]]<gids[ux_offsets[i]]);
      }

      dofManager->getElementGIDs(1,gids);
      TEST_EQUALITY(gids.size(),12);
      for(std::size_t i=0;i<uy_offsets.size();i++) {
         TEST_ASSERT(gids[uy_offsets[i]]<gids[p_offsets[i]]);
         TEST_ASSERT(gids[uy_offsets[i]]<gids[ux_offsets[i]]);
         TEST_ASSERT(gids[p_offsets[i]]<gids[ux_offsets[i]]);
      }
   }
   else if(myRank==1) {
      std::vector<int> gids;

      dofManager->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),12);
      for(std::size_t i=0;i<uy_offsets.size();i++) {
         TEST_ASSERT(gids[uy_offsets[i]]<gids[p_offsets[i]]);
         TEST_ASSERT(gids[uy_offsets[i]]<gids[ux_offsets[i]]);
         TEST_ASSERT(gids[p_offsets[i]]<gids[ux_offsets[i]]);
      }

      dofManager->getElementGIDs(1,gids);
      TEST_EQUALITY(gids.size(),12);
      for(std::size_t i=0;i<uy_offsets.size();i++) {
         TEST_ASSERT(gids[uy_offsets[i]]<gids[p_offsets[i]]);
         TEST_ASSERT(gids[uy_offsets[i]]<gids[ux_offsets[i]]);
         TEST_ASSERT(gids[p_offsets[i]]<gids[ux_offsets[i]]);
      }
   }
}

// quad tests
TEUCHOS_UNIT_TEST(tSquareTriMeshDOFManager, ghosted_owned_indices)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      stk::ParallelMachine Comm = MPI_COMM_WORLD;
   #else
      stk::ParallelMachine Comm = WHAT_TO_DO_COMM;
   #endif

   int numProcs = stk::parallel_machine_size(Comm);
   int myRank = stk::parallel_machine_rank(Comm);

   TEUCHOS_ASSERT(numProcs==2);

   // build a geometric pattern from a single basis
   RCP<const panzer::FieldPattern> patternC1
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();

   // build DOF manager
   RCP<panzer::ConnManager> connManager = buildTriMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager> dofManager = rcp(new panzer::DOFManager());
   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   dofManager->addField("u",patternC1);
   dofManager->buildGlobalUnknowns();

   // test GlobalIndexer
   RCP<panzer::GlobalIndexer> glbNum = dofManager;

   std::vector<int> owned, ownedAndGhosted;
   glbNum->getOwnedIndices(owned);
   glbNum->getOwnedAndGhostedIndices(ownedAndGhosted);

   if(myRank==0 && numProcs==2) {
      TEST_EQUALITY(owned.size(),6);
      TEST_EQUALITY(ownedAndGhosted.size(),6);
      bool ownedCorrect = true;
      bool ownedAndGhostedCorrect = true;

      std::sort(owned.begin(),owned.end());
      for(std::size_t i=0;i<owned.size();i++) {
         ownedCorrect &= (owned[i] == (int) i);
      }
      TEST_ASSERT(ownedCorrect);

      std::sort(ownedAndGhosted.begin(),ownedAndGhosted.end());
      for(std::size_t i=0;i<ownedAndGhosted.size();i++) {
         ownedAndGhostedCorrect &= (ownedAndGhosted[i] == (int) i);
      }
      TEST_ASSERT(ownedAndGhostedCorrect);
   }
   else if(myRank==1 && numProcs==2) {
      TEST_EQUALITY(owned.size(),3);
      TEST_EQUALITY(ownedAndGhosted.size(),6);
      bool ownedCorrect = true;
      bool ownedAndGhostedCorrect = true;

      std::sort(owned.begin(),owned.end());
      for(std::size_t i=0;i<owned.size();i++) {
         ownedCorrect &= (owned[i] == (int) i+6);
      }
      TEST_ASSERT(ownedCorrect);

      std::sort(ownedAndGhosted.begin(),ownedAndGhosted.end());
      for(std::size_t i=0;i<3;i++) {
         ownedAndGhostedCorrect &= (ownedAndGhosted[i] == (int) (2*i+1));
      }
      for(std::size_t i=0;i<3;i++) {
         ownedAndGhostedCorrect &= (ownedAndGhosted[i+3] == (int) i+6);
      }
      TEST_ASSERT(ownedAndGhostedCorrect);
   }
   else
      TEUCHOS_ASSERT(false);
}

// quad tests
TEUCHOS_UNIT_TEST(tSquareTriMeshDOFManager, multiple_dof_managers)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      stk::ParallelMachine Comm = MPI_COMM_WORLD;
   #else
      stk::ParallelMachine Comm = WHAT_TO_DO_COMM;
   #endif

   int numProcs = stk::parallel_machine_size(Comm);
   int myRank = stk::parallel_machine_rank(Comm);

   TEUCHOS_ASSERT(numProcs==2);

   // build a geometric pattern from a single basis
   RCP<const panzer::FieldPattern> patternC1
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();
   RCP<const panzer::FieldPattern> patternC2
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::exec_space,double,double> >();

   // build DOF manager
   RCP<panzer::ConnManager> connManager = buildTriMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager> dofManager_fluids = rcp(new panzer::DOFManager());
   dofManager_fluids->setConnManager(connManager,MPI_COMM_WORLD);
   dofManager_fluids->addField("ux",patternC2);
   dofManager_fluids->addField("uy",patternC2);
   dofManager_fluids->addField("p",patternC1);
   dofManager_fluids->buildGlobalUnknowns();

   RCP<panzer::DOFManager> dofManager_temp = rcp(new panzer::DOFManager());
   dofManager_temp->setConnManager(connManager,MPI_COMM_WORLD);
   dofManager_temp->addField("T",patternC1);
   dofManager_temp->buildGlobalUnknowns(dofManager_fluids->getGeometricFieldPattern());

   if(myRank==0) {
      std::vector<int> gids;

      dofManager_temp->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),4);
      TEST_EQUALITY(gids[0],0); TEST_EQUALITY(gids[1],1); TEST_EQUALITY(gids[2],3);
      TEST_EQUALITY(gids[3],2);

      dofManager_temp->getElementGIDs(1,gids);
      TEST_EQUALITY(gids.size(),4);
      TEST_EQUALITY(gids[0],2); TEST_EQUALITY(gids[1],3); TEST_EQUALITY(gids[2],5);
      TEST_EQUALITY(gids[3],4);
   }
   else if(myRank==1) {
      std::vector<int> gids;

      dofManager_temp->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),4);
      TEST_EQUALITY(gids[0],1); TEST_EQUALITY(gids[1],6); TEST_EQUALITY(gids[2],7);
      TEST_EQUALITY(gids[3],3);

      dofManager_temp->getElementGIDs(1,gids);
      TEST_EQUALITY(gids.size(),4);
      TEST_EQUALITY(gids[0],3); TEST_EQUALITY(gids[1],7); TEST_EQUALITY(gids[2],8);
      TEST_EQUALITY(gids[3],5);
   }
   else
      TEUCHOS_ASSERT(false);
}
*/

/*
TEUCHOS_UNIT_TEST(tSquareTriMeshDOFManager,getDofCoords)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      stk::ParallelMachine Comm = MPI_COMM_WORLD;
   #else
      stk::ParallelMachine Comm = WHAT_TO_DO_COMM;
   #endif

   int numProcs = stk::parallel_machine_size(Comm);
   int myRank = stk::parallel_machine_rank(Comm);

   TEUCHOS_ASSERT(numProcs==2);
   // build DOF manager
   RCP<panzer::ConnManager> connManager = buildTriMesh(Comm,2,2,2,1);
   RCP<const panzer_stk::STKConnManager> stkManager = rcp_dynamic_cast<panzer_stk::STKConnManager>(connManager);
   RCP<panzer_stk::STK_Interface> meshDB = stkManager->getSTKInterface();
   meshDB->print(out);

   // grab elements from mesh
   std::vector<stk::mesh::Entity> block00, block01;
   meshDB->getMyElements("eblock-0_0",block00);
   meshDB->getMyElements("eblock-1_0",block01);

   std::vector<std::size_t> localIds_00, localIds_01;
   FieldContainer coords00, coords01;
   RCP<const panzer::Intrepid2FieldPattern> patternC1_00
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();
   RCP<const panzer::Intrepid2FieldPattern> patternC1_01
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C2_FEM<PHX::exec_space,double,double> >();

   // get coordinates
   stkManager->getDofCoords("eblock-0_0",*patternC1_00,localIds_00,coords00);
   stkManager->getDofCoords("eblock-1_0",*patternC1_01,localIds_01,coords01);

   TEST_EQUALITY(localIds_00.size(),block00.size());
   TEST_EQUALITY(localIds_01.size(),block01.size());

   TEST_EQUALITY(coords00.extent(0),int(localIds_00.size()));
   TEST_EQUALITY(coords01.extent(0),int(localIds_01.size()));

   TEST_EQUALITY(coords00.extent(1),4); TEST_EQUALITY(coords00.extent(2),2);
   TEST_EQUALITY(coords01.extent(1),9); TEST_EQUALITY(coords01.extent(2),2);

   for(std::size_t i=0;i<block00.size();i++)
      TEST_EQUALITY(localIds_00[i],meshDB->elementLocalId(block00[i]));
   for(std::size_t i=0;i<block01.size();i++)
      TEST_EQUALITY(localIds_01[i],meshDB->elementLocalId(block01[i]));

   // for(std::size_t c=0;c<block00.size();c++) {
   //    stk::mesh::Entity element = block00[c];
   //    for(int i=0;i<4;i++) {
   //    }
   // }
}
*/

}
