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
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HCURL_QUAD_I1_FEM.hpp"

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;

namespace panzer_stk {

Teuchos::RCP<panzer::ConnManager> buildQuadMesh(stk::ParallelMachine comm,int xelmts,int yelmts,int xblocks,int yblocks)
{
   Teuchos::ParameterList pl;
   pl.set<int>("X Elements",xelmts);
   pl.set<int>("Y Elements",yelmts);
   pl.set<int>("X Blocks",xblocks);
   pl.set<int>("Y Blocks",yblocks);

   panzer_stk::SquareQuadMeshFactory meshFact;
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
TEUCHOS_UNIT_TEST(tSquareQuadMeshDOFManager, buildTest_quad)
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

   RCP<panzer::ConnManager> connManager = buildQuadMesh(Comm,2,2,1,1);
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
      TEST_EQUALITY(gids.size(),12);
      TEST_EQUALITY(gids[0],0); TEST_EQUALITY(gids[1],1); TEST_EQUALITY(gids[2],2);
      TEST_EQUALITY(gids[3],3); TEST_EQUALITY(gids[4],4); TEST_EQUALITY(gids[5],5);
      TEST_EQUALITY(gids[6],9); TEST_EQUALITY(gids[7],10); TEST_EQUALITY(gids[8],11);
      TEST_EQUALITY(gids[9],6); TEST_EQUALITY(gids[10],7); TEST_EQUALITY(gids[11],8);

      for(std::size_t i=0;i<p_offsets.size();i++) {
         TEST_ASSERT(gids[p_offsets[i]]<gids[ux_offsets[i]]);
         TEST_ASSERT(gids[p_offsets[i]]<gids[uy_offsets[i]]);
         TEST_ASSERT(gids[ux_offsets[i]]<gids[uy_offsets[i]]);
      }

      dofManager->getElementGIDs(1,gids);
      TEST_EQUALITY(gids.size(),12);
      TEST_EQUALITY(gids[0],6); TEST_EQUALITY(gids[1],7); TEST_EQUALITY(gids[2],8);
      TEST_EQUALITY(gids[3],9); TEST_EQUALITY(gids[4],10); TEST_EQUALITY(gids[5],11);
      TEST_EQUALITY(gids[6],15); TEST_EQUALITY(gids[7],16); TEST_EQUALITY(gids[8],17);
      TEST_EQUALITY(gids[9],12); TEST_EQUALITY(gids[10],13); TEST_EQUALITY(gids[11],14);

      for(std::size_t i=0;i<p_offsets.size();i++) {
         TEST_ASSERT(gids[p_offsets[i]]<gids[ux_offsets[i]]);
         TEST_ASSERT(gids[p_offsets[i]]<gids[uy_offsets[i]]);
         TEST_ASSERT(gids[ux_offsets[i]]<gids[uy_offsets[i]]);
      }
   }
   else if(myRank==1) {
      std::vector<panzer::GlobalOrdinal> gids;

      dofManager->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),12);
      TEST_EQUALITY(gids[0],3); TEST_EQUALITY(gids[1],4); TEST_EQUALITY(gids[2],5);
      TEST_EQUALITY(gids[3],18); TEST_EQUALITY(gids[4],19); TEST_EQUALITY(gids[5],20);
      TEST_EQUALITY(gids[6],21); TEST_EQUALITY(gids[7],22); TEST_EQUALITY(gids[8],23);
      TEST_EQUALITY(gids[9],9); TEST_EQUALITY(gids[10],10); TEST_EQUALITY(gids[11],11);

      for(std::size_t i=0;i<p_offsets.size();i++) {
         TEST_ASSERT(gids[p_offsets[i]]<gids[ux_offsets[i]]);
         TEST_ASSERT(gids[p_offsets[i]]<gids[uy_offsets[i]]);
         TEST_ASSERT(gids[ux_offsets[i]]<gids[uy_offsets[i]]);
      }

      dofManager->getElementGIDs(1,gids);
      TEST_EQUALITY(gids.size(),12);
      TEST_EQUALITY(gids[0],9); TEST_EQUALITY(gids[1],10); TEST_EQUALITY(gids[2],11);
      TEST_EQUALITY(gids[3],21); TEST_EQUALITY(gids[4],22); TEST_EQUALITY(gids[5],23);
      TEST_EQUALITY(gids[6],24); TEST_EQUALITY(gids[7],25); TEST_EQUALITY(gids[8],26);
      TEST_EQUALITY(gids[9],15); TEST_EQUALITY(gids[10],16); TEST_EQUALITY(gids[11],17);

      for(std::size_t i=0;i<p_offsets.size();i++) {
         TEST_ASSERT(gids[p_offsets[i]]<gids[ux_offsets[i]]);
         TEST_ASSERT(gids[p_offsets[i]]<gids[uy_offsets[i]]);
         TEST_ASSERT(gids[ux_offsets[i]]<gids[uy_offsets[i]]);
      }
   }
}

TEUCHOS_UNIT_TEST(tSquareQuadMeshDOFManager, field_order)
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

   RCP<panzer::ConnManager> connManager = buildQuadMesh(Comm,2,2,1,1);
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
      std::vector<panzer::GlobalOrdinal> gids;

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
      std::vector<panzer::GlobalOrdinal> gids;

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
TEUCHOS_UNIT_TEST(tSquareQuadMeshDOFManager, ghosted_owned_indices)
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
   RCP<panzer::ConnManager> connManager = buildQuadMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager> dofManager = rcp(new panzer::DOFManager());
   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   dofManager->addField("u",patternC1);
   dofManager->buildGlobalUnknowns();

   // test GlobalIndexer
   RCP<panzer::GlobalIndexer> glbNum = dofManager;

   std::vector<panzer::GlobalOrdinal> owned, ownedAndGhosted;
   glbNum->getOwnedIndices(owned);
   glbNum->getOwnedAndGhostedIndices(ownedAndGhosted);

   if(myRank==0 && numProcs==2) {
      TEST_EQUALITY(owned.size(),6);
      TEST_EQUALITY(ownedAndGhosted.size(),6);
      bool ownedCorrect = true;
      bool ownedAndGhostedCorrect = true;

      std::sort(owned.begin(),owned.end());
      ownedCorrect &= (owned[0] == (int) 0);
      ownedCorrect &= (owned[1] == (int) 1);
      ownedCorrect &= (owned[2] == (int) 2);
      ownedCorrect &= (owned[3] == (int) 3);
      ownedCorrect &= (owned[4] == (int) 4);
      ownedCorrect &= (owned[5] == (int) 5);
      TEST_ASSERT(ownedCorrect);

      std::sort(ownedAndGhosted.begin(),ownedAndGhosted.end());
      ownedAndGhostedCorrect &= (ownedAndGhosted[0] == (int) 0);
      ownedAndGhostedCorrect &= (ownedAndGhosted[1] == (int) 1);
      ownedAndGhostedCorrect &= (ownedAndGhosted[2] == (int) 2);
      ownedAndGhostedCorrect &= (ownedAndGhosted[3] == (int) 3);
      ownedAndGhostedCorrect &= (ownedAndGhosted[4] == (int) 4);
      ownedAndGhostedCorrect &= (ownedAndGhosted[5] == (int) 5);
      TEST_ASSERT(ownedAndGhostedCorrect);
   }
   else if(myRank==1 && numProcs==2) {
      TEST_EQUALITY(owned.size(),3);
      TEST_EQUALITY(ownedAndGhosted.size(),6);
      bool ownedCorrect = true;
      bool ownedAndGhostedCorrect = true;

      std::sort(owned.begin(),owned.end());
      ownedCorrect &= (owned[0] == (int) 6);
      ownedCorrect &= (owned[1] == (int) 7);
      ownedCorrect &= (owned[2] == (int) 8);
      TEST_ASSERT(ownedCorrect);

      std::sort(ownedAndGhosted.begin(),ownedAndGhosted.end());
      ownedAndGhostedCorrect &= (ownedAndGhosted[0] == (int) 1);
      ownedAndGhostedCorrect &= (ownedAndGhosted[1] == (int) 3);
      ownedAndGhostedCorrect &= (ownedAndGhosted[2] == (int) 5);
      ownedAndGhostedCorrect &= (ownedAndGhosted[3] == (int) 6);
      ownedAndGhostedCorrect &= (ownedAndGhosted[4] == (int) 7);
      ownedAndGhostedCorrect &= (ownedAndGhosted[5] == (int) 8);
      TEST_ASSERT(ownedAndGhostedCorrect);
   }
   else
      TEUCHOS_ASSERT(false);
}

// quad tests
TEUCHOS_UNIT_TEST(tSquareQuadMeshDOFManager, multiple_dof_managers)
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
   RCP<panzer::ConnManager> connManager = buildQuadMesh(Comm,2,2,1,1);
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
      std::vector<panzer::GlobalOrdinal> gids;

      dofManager_temp->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),4);
      TEST_EQUALITY(gids[0],0);
      TEST_EQUALITY(gids[1],1);
      TEST_EQUALITY(gids[2],3);
      TEST_EQUALITY(gids[3],2);

      dofManager_temp->getElementGIDs(1,gids);
      TEST_EQUALITY(gids.size(),4);
      TEST_EQUALITY(gids[0],2);
      TEST_EQUALITY(gids[1],3);
      TEST_EQUALITY(gids[2],5);
      TEST_EQUALITY(gids[3],4);
   }
   else if(myRank==1) {
      std::vector<panzer::GlobalOrdinal> gids;

      dofManager_temp->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),4);
      TEST_EQUALITY(gids[0],1);
      TEST_EQUALITY(gids[1],6);
      TEST_EQUALITY(gids[2],7);
      TEST_EQUALITY(gids[3],3);

      dofManager_temp->getElementGIDs(1,gids);
      TEST_EQUALITY(gids.size(),4);
      TEST_EQUALITY(gids[0],3);
      TEST_EQUALITY(gids[1],7);
      TEST_EQUALITY(gids[2],8);
      TEST_EQUALITY(gids[3],5);
   }
   else
      TEUCHOS_ASSERT(false);
}

TEUCHOS_UNIT_TEST(tSquareQuadMeshDOFManager,getDofCoords)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      stk::ParallelMachine Comm = MPI_COMM_WORLD;
   #else
      stk::ParallelMachine Comm = WHAT_TO_DO_COMM;
   #endif

   int numProcs = stk::parallel_machine_size(Comm);

   TEUCHOS_ASSERT(numProcs==2);
   // build DOF manager
   RCP<panzer::ConnManager> connManager = buildQuadMesh(Comm,2,2,2,1);
   RCP<const panzer_stk::STKConnManager> stkManager = rcp_dynamic_cast<panzer_stk::STKConnManager>(connManager);
   RCP<const panzer_stk::STK_Interface> meshDB = stkManager->getSTKInterface();
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

   TEST_EQUALITY(static_cast<int>(coords00.extent(0)),
     static_cast<int>(localIds_00.size()))
   TEST_EQUALITY(static_cast<int>(coords01.extent(0)),
     static_cast<int>(localIds_01.size()))

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

// quad tests
TEUCHOS_UNIT_TEST(tSquareQuadMeshDOFManager, buildTest_quad_edge_orientations)
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
   RCP<const panzer::FieldPattern> patternI1
         = buildFieldPattern<Intrepid2::Basis_HCURL_QUAD_I1_FEM<PHX::exec_space,double,double> >();

   RCP<panzer::ConnManager> connManager = buildQuadMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager> dofManager = rcp(new panzer::DOFManager());

   dofManager->setOrientationsRequired(true);
   TEST_EQUALITY(dofManager->getOrientationsRequired(),true);
   TEST_EQUALITY(dofManager->getConnManager(),Teuchos::null);

   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   TEST_EQUALITY(dofManager->getConnManager(),connManager);

   dofManager->addField("u",patternC1);
   dofManager->addField("b",patternI1);

   dofManager->buildGlobalUnknowns();
   dofManager->printFieldInformation(out);

   const std::vector<int> & u_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("u"));
   const std::vector<int> & b_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("b"));

   TEST_EQUALITY(u_offsets.size(),4);
   TEST_EQUALITY(b_offsets.size(),4);

   // unfortunatly this mesh is completly uniform
   double standardO[] = { 1.0, 1.0, -1.0, -1.0 };
   if(myRank==0) {
      std::vector<panzer::GlobalOrdinal> gids;
      std::vector<double> orientation;

      // element 0
      dofManager->getElementGIDs(0,gids);
      dofManager->getElementOrientation(0,orientation);
      TEST_EQUALITY(gids.size(),8);
      TEST_EQUALITY(orientation.size(),8);

      {
         double sum = 0;
         for(std::size_t i=0;i<orientation.size();i++)
            sum += orientation[i]*orientation[i];
         TEST_EQUALITY(sum,8);
      }

      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i],1.0);
      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i+4],standardO[i]);

      // element 1
      dofManager->getElementGIDs(1,gids);
      dofManager->getElementOrientation(1,orientation);
      TEST_EQUALITY(gids.size(),8);
      TEST_EQUALITY(orientation.size(),8);

      {
         double sum = 0;
         for(std::size_t i=0;i<orientation.size();i++)
            sum += orientation[i]*orientation[i];
         TEST_EQUALITY(sum,8);
      }

      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i],1.0);
      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i+4],standardO[i]);
   }
   else if(myRank==1) {
      std::vector<panzer::GlobalOrdinal> gids;
      std::vector<double> orientation;

      // element 0
      dofManager->getElementOrientation(0,orientation);
      dofManager->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),8);
      TEST_EQUALITY(orientation.size(),8);

      {
         double sum = 0;
         for(std::size_t i=0;i<orientation.size();i++)
            sum += orientation[i]*orientation[i];
         TEST_EQUALITY(sum,8);
      }

      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i],1.0);
      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i+4],standardO[i]);

      // element 1
      dofManager->getElementGIDs(1,gids);
      dofManager->getElementOrientation(1,orientation);
      TEST_EQUALITY(gids.size(),8);
      TEST_EQUALITY(orientation.size(),8);

      {
         double sum = 0;
         for(std::size_t i=0;i<orientation.size();i++)
            sum += orientation[i]*orientation[i];
         TEST_EQUALITY(sum,8);
      }

      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i],1.0);
      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i+4],standardO[i]);
   }
}

// quad tests
TEUCHOS_UNIT_TEST(tSquareQuadMeshDOFManager, buildTest_quad_edge_orientations2)
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
   RCP<const panzer::FieldPattern> patternI1
         = buildFieldPattern<Intrepid2::Basis_HCURL_QUAD_I1_FEM<PHX::exec_space,double,double> >();

   RCP<panzer::ConnManager> connManager = buildQuadMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager> dofManager = rcp(new panzer::DOFManager());

   dofManager->setOrientationsRequired(true);
   TEST_EQUALITY(dofManager->getOrientationsRequired(),true);
   TEST_EQUALITY(dofManager->getConnManager(),Teuchos::null);

   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   TEST_EQUALITY(dofManager->getConnManager(),connManager);

   dofManager->addField("b",patternI1);

   dofManager->buildGlobalUnknowns();
   dofManager->printFieldInformation(out);

   const std::vector<int> & b_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("b"));

   TEST_EQUALITY(b_offsets.size(),4);

   // unfortunatly this mesh is completly uniform
   double standardO[] = { 1.0, 1.0, -1.0, -1.0 };
   if(myRank==0) {
      std::vector<double> orientation;

      // element 0
      dofManager->getElementOrientation(0,orientation);
      TEST_EQUALITY(orientation.size(),4);

      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i],standardO[i]);

      // element 1
      dofManager->getElementOrientation(1,orientation);
      TEST_EQUALITY(orientation.size(),4);

      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i],standardO[i]);
   }
   else if(myRank==1) {
      std::vector<double> orientation;

      // element 0
      dofManager->getElementOrientation(0,orientation);
      TEST_EQUALITY(orientation.size(),4);

      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i],standardO[i]);

      // element 1
      dofManager->getElementOrientation(1,orientation);
      TEST_EQUALITY(orientation.size(),4);

      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i],standardO[i]);
   }
}

// quad tests

// This test verifies that the DOFManager can determine when you can't compute
// orientations because the field pattern is not correct (it must include nodes
// to work correctly). Thus its the ability of the DOFManager to protect itself
// from the insane.
TEUCHOS_UNIT_TEST(tSquareQuadMeshDOFManager, buildTest_quad_edge_orientations_fail)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      stk::ParallelMachine Comm = MPI_COMM_WORLD;
   #else
      stk::ParallelMachine Comm = WHAT_TO_DO_COMM;
   #endif

   int numProcs = stk::parallel_machine_size(Comm);

   TEUCHOS_ASSERT(numProcs==2);

   // build a geometric pattern from a single basis
   RCP<const panzer::FieldPattern> patternI1
         = buildFieldPattern<Intrepid2::Basis_HCURL_QUAD_I1_FEM<PHX::exec_space,double,double> >();

   RCP<panzer::ConnManager> connManager = buildQuadMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager> dofManager = rcp(new panzer::DOFManager());

   dofManager->setOrientationsRequired(true);
   TEST_EQUALITY(dofManager->getOrientationsRequired(),true);
   TEST_EQUALITY(dofManager->getConnManager(),Teuchos::null);

   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   TEST_EQUALITY(dofManager->getConnManager(),connManager);

   dofManager->addField("b",patternI1);

   TEST_THROW(dofManager->buildGlobalUnknowns(patternI1),std::logic_error);
}

TEUCHOS_UNIT_TEST(tSquareQuadMeshDOFManager, buildTest_q2q1)
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

   RCP<panzer::ConnManager> connManager = buildQuadMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager> dofManager = rcp(new panzer::DOFManager());

   TEST_EQUALITY(dofManager->getOrientationsRequired(),false);
   TEST_EQUALITY(dofManager->getConnManager(),Teuchos::null);

   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   TEST_EQUALITY(dofManager->getConnManager(),connManager);

   dofManager->addField("ux",patternC2);
   dofManager->addField("uy",patternC2);
   dofManager->addField("p",patternC1);

   std::vector<std::string> fieldOrder;
   fieldOrder.push_back("ux");
   fieldOrder.push_back("uy");
   fieldOrder.push_back("p");
   dofManager->setFieldOrder(fieldOrder);

   dofManager->buildGlobalUnknowns();
   dofManager->printFieldInformation(out);

   const std::vector<int> & ux_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("ux"));
   const std::vector<int> & uy_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("uy"));
   const std::vector<int> & p_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("p"));

   TEST_EQUALITY(ux_offsets.size(),9);
   TEST_EQUALITY(ux_offsets.size(),uy_offsets.size());
   TEST_EQUALITY(p_offsets.size(),4);

   TEST_EQUALITY(ux_offsets[0],0);
   TEST_EQUALITY(ux_offsets[1],3);
   TEST_EQUALITY(ux_offsets[2],6);
   TEST_EQUALITY(ux_offsets[3],9);
   TEST_EQUALITY(ux_offsets[4],12);
   TEST_EQUALITY(ux_offsets[5],14);
   TEST_EQUALITY(ux_offsets[6],16);
   TEST_EQUALITY(ux_offsets[7],18);
   TEST_EQUALITY(ux_offsets[8],20);

   TEST_EQUALITY(uy_offsets[0],1);
   TEST_EQUALITY(uy_offsets[1],4);
   TEST_EQUALITY(uy_offsets[2],7);
   TEST_EQUALITY(uy_offsets[3],10);
   TEST_EQUALITY(uy_offsets[4],13);
   TEST_EQUALITY(uy_offsets[5],15);
   TEST_EQUALITY(uy_offsets[6],17);
   TEST_EQUALITY(uy_offsets[7],19);
   TEST_EQUALITY(uy_offsets[8],21);

   TEST_EQUALITY(p_offsets[0],2);
   TEST_EQUALITY(p_offsets[1],5);
   TEST_EQUALITY(p_offsets[2],8);
   TEST_EQUALITY(p_offsets[3],11);

   if(myRank==0) {
      std::vector<panzer::GlobalOrdinal> gids;

      dofManager->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),9+9+4);

      // nodes
      TEST_EQUALITY(gids[0],0); TEST_EQUALITY(gids[1],1); TEST_EQUALITY(gids[2],2);
      TEST_EQUALITY(gids[3],3); TEST_EQUALITY(gids[4],4); TEST_EQUALITY(gids[5],5);
      TEST_EQUALITY(gids[6],9); TEST_EQUALITY(gids[7],10); TEST_EQUALITY(gids[8],11);
      TEST_EQUALITY(gids[9],6); TEST_EQUALITY(gids[10],7); TEST_EQUALITY(gids[11],8);

      // edges
      TEST_EQUALITY(gids[12],18); TEST_EQUALITY(gids[13],19);
      TEST_EQUALITY(gids[14],20); TEST_EQUALITY(gids[15],21);
      TEST_EQUALITY(gids[16],22); TEST_EQUALITY(gids[17],23);
      TEST_EQUALITY(gids[18],24); TEST_EQUALITY(gids[19],25);
      TEST_EQUALITY(gids[20],32); TEST_EQUALITY(gids[21],33);

      dofManager->getElementGIDs(1,gids);
      TEST_EQUALITY(gids.size(),9+9+4);

      // nodes
      TEST_EQUALITY(gids[0],6); TEST_EQUALITY(gids[1],7); TEST_EQUALITY(gids[2],8);
      TEST_EQUALITY(gids[3],9); TEST_EQUALITY(gids[4],10); TEST_EQUALITY(gids[5],11);
      TEST_EQUALITY(gids[6],15); TEST_EQUALITY(gids[7],16); TEST_EQUALITY(gids[8],17);
      TEST_EQUALITY(gids[9],12); TEST_EQUALITY(gids[10],13); TEST_EQUALITY(gids[11],14);

      // edges
      TEST_EQUALITY(gids[12],22); TEST_EQUALITY(gids[13],23);
      TEST_EQUALITY(gids[14],26); TEST_EQUALITY(gids[15],27);
      TEST_EQUALITY(gids[16],28); TEST_EQUALITY(gids[17],29);
      TEST_EQUALITY(gids[18],30); TEST_EQUALITY(gids[19],31);
      TEST_EQUALITY(gids[20],34); TEST_EQUALITY(gids[21],35);
   }
   else if(myRank==1) {
      std::vector<panzer::GlobalOrdinal> gids;

      dofManager->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),9+9+4);

      // nodes
      TEST_EQUALITY(gids[0],3); TEST_EQUALITY(gids[1],4); TEST_EQUALITY(gids[2],5);
      TEST_EQUALITY(gids[3],36); TEST_EQUALITY(gids[4],37); TEST_EQUALITY(gids[5],38);
      TEST_EQUALITY(gids[6],39); TEST_EQUALITY(gids[7],40); TEST_EQUALITY(gids[8],41);
      TEST_EQUALITY(gids[9],9); TEST_EQUALITY(gids[10],10); TEST_EQUALITY(gids[11],11);

      // edges
      TEST_EQUALITY(gids[12],45); TEST_EQUALITY(gids[13],46);
      TEST_EQUALITY(gids[14],47); TEST_EQUALITY(gids[15],48);
      TEST_EQUALITY(gids[16],49); TEST_EQUALITY(gids[17],50);
      TEST_EQUALITY(gids[18],20); TEST_EQUALITY(gids[19],21);
      TEST_EQUALITY(gids[20],55); TEST_EQUALITY(gids[21],56);

      dofManager->getElementGIDs(1,gids);
      TEST_EQUALITY(gids.size(),9+9+4);

      // nodes
      TEST_EQUALITY(gids[0],9); TEST_EQUALITY(gids[1],10); TEST_EQUALITY(gids[2],11);
      TEST_EQUALITY(gids[3],39); TEST_EQUALITY(gids[4],40); TEST_EQUALITY(gids[5],41);
      TEST_EQUALITY(gids[6],42); TEST_EQUALITY(gids[7],43); TEST_EQUALITY(gids[8],44);
      TEST_EQUALITY(gids[9],15); TEST_EQUALITY(gids[10],16); TEST_EQUALITY(gids[11],17);

      // edges
      TEST_EQUALITY(gids[12],49); TEST_EQUALITY(gids[13],50);
      TEST_EQUALITY(gids[14],51); TEST_EQUALITY(gids[15],52);
      TEST_EQUALITY(gids[16],53); TEST_EQUALITY(gids[17],54);
      TEST_EQUALITY(gids[18],26); TEST_EQUALITY(gids[19],27);
      TEST_EQUALITY(gids[20],57); TEST_EQUALITY(gids[21],58);
   }

   std::vector<panzer::GlobalOrdinal> owned, ownedAndGhosted;
   dofManager->getOwnedIndices(owned);
   dofManager->getOwnedAndGhostedIndices(ownedAndGhosted);

   if(myRank==0) {
      TEST_EQUALITY(owned.size(),36);
      TEST_EQUALITY(ownedAndGhosted.size(),36);
   }
   else if(myRank==1) {
      TEST_EQUALITY(owned.size(),59-36);
      TEST_EQUALITY(ownedAndGhosted.size(),36);
   }
   else
      TEUCHOS_ASSERT(false);
}

TEUCHOS_UNIT_TEST(tSquareQuadMeshDOFManager, buildTest_nabors)
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

   RCP<panzer::ConnManager> connManager = buildQuadMesh(Comm,4,2,1,1);
   RCP<panzer::DOFManager> dofManager = rcp(new panzer::DOFManager());
   RCP<panzer::DOFManager> dofManager_noNeighbors =
    rcp(new panzer::DOFManager());
   dofManager->useNeighbors(true);

   TEST_EQUALITY(dofManager->getOrientationsRequired(),false);
   TEST_EQUALITY(dofManager->getConnManager(),Teuchos::null);

   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   dofManager_noNeighbors->setConnManager(connManager,MPI_COMM_WORLD);
   TEST_EQUALITY(dofManager->getConnManager(),connManager);

   dofManager->addField("ux",patternC1);
   dofManager->addField("uy",patternC1);
   dofManager->addField("p",patternC1);
   dofManager_noNeighbors->addField("ux",patternC1);
   dofManager_noNeighbors->addField("uy",patternC1);
   dofManager_noNeighbors->addField("p",patternC1);

   dofManager->buildGlobalUnknowns();
   dofManager_noNeighbors->buildGlobalUnknowns();
   dofManager->printFieldInformation(out);

   TEST_EQUALITY(connManager->getElementBlock("eblock-0_0").size(),4);
   TEST_EQUALITY(connManager->getNeighborElementBlock("eblock-0_0").size(),2);

   TEST_EQUALITY(dofManager->getNumberElementGIDArrays(),6)

   if(myRank==0) {
      std::vector<panzer::GlobalOrdinal> gids;

      for(int i=0;i<6;i++) {
        dofManager->getElementGIDs(i,gids);

        out << "Element " << i << ": ";
        for(int j=0;j<12;j+=3)
          out << gids[j] << " ";
        out << std::endl;
      }

      dofManager->getElementGIDs(4,gids);
      TEST_EQUALITY(gids.size(),4+4+4);

      TEST_EQUALITY(gids[0],6); TEST_EQUALITY(gids[1],7); TEST_EQUALITY(gids[2],8);
      TEST_EQUALITY(gids[3],27); TEST_EQUALITY(gids[4],28); TEST_EQUALITY(gids[5],29);
      TEST_EQUALITY(gids[6],33); TEST_EQUALITY(gids[7],34); TEST_EQUALITY(gids[8],35);
      TEST_EQUALITY(gids[9],15); TEST_EQUALITY(gids[10],16); TEST_EQUALITY(gids[11],17);

      dofManager->getElementGIDs(5,gids);
      TEST_EQUALITY(gids.size(),4+4+4);

      TEST_EQUALITY(gids[0],15); TEST_EQUALITY(gids[1],16); TEST_EQUALITY(gids[2],17);
      TEST_EQUALITY(gids[3],33); TEST_EQUALITY(gids[4],34); TEST_EQUALITY(gids[5],35);
      TEST_EQUALITY(gids[6],39); TEST_EQUALITY(gids[7],40); TEST_EQUALITY(gids[8],41);
      TEST_EQUALITY(gids[9],24); TEST_EQUALITY(gids[10],25); TEST_EQUALITY(gids[11],26);
   }
   else {
      std::vector<panzer::GlobalOrdinal> gids;

      for(int i=0;i<6;i++) {
        dofManager->getElementGIDs(i,gids);

        out << "Element " << i << ": ";
        for(int j=0;j<12;j+=3)
          out << gids[j] << " ";
        out << std::endl;
      }

      dofManager->getElementGIDs(4,gids);
      TEST_EQUALITY(gids.size(),4+4+4);

      dofManager->getElementGIDs(5,gids);
      TEST_EQUALITY(gids.size(),4+4+4);
   }

   // owned vector
   {
     std::vector<panzer::GlobalOrdinal> owned, owned_noNeighbors;
     dofManager->getOwnedIndices(owned);
     dofManager_noNeighbors->getOwnedIndices(owned_noNeighbors);
     TEST_EQUALITY(owned.size(),owned_noNeighbors.size());

     bool owned_result = true;
     for(std::size_t j=0;j<owned.size();j++)
       owned_result &= (owned[j]==owned_noNeighbors[j]);
     TEST_ASSERT(owned_result);
   }

   // owned and ghosted vector
   {
     std::vector<panzer::GlobalOrdinal> ghosted;
     dofManager->getOwnedAndGhostedIndices(ghosted);

     std::set<panzer::GlobalOrdinal> ghosted_set;
     ghosted_set.insert(ghosted.begin(),ghosted.end());
     TEST_EQUALITY(ghosted_set.size(),ghosted.size()); // make sure there are no duplicated entries

     for(int e=0;e<6;e++) {
       std::vector<panzer::GlobalOrdinal> gids;

       dofManager->getElementGIDs(e,gids);
       bool allFound = true;
       for(std::size_t i=0;i<gids.size();i++)
         allFound &= (ghosted_set.find(gids[i])!=ghosted_set.end());

       TEST_ASSERT(allFound);
     }
   }
}

}
