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
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

namespace panzer_stk {

Teuchos::RCP<panzer_stk::STK_Interface> buildHexMesh(stk::ParallelMachine comm,int xelmts,int yelmts,int zelmts,
                                                                                    int xblocks,int yblocks,int zblocks)
{
   Teuchos::ParameterList pl;
   pl.set<int>("X Elements",xelmts);
   pl.set<int>("Y Elements",yelmts);
   pl.set<int>("Z Elements",zelmts);
   pl.set<int>("X Blocks",xblocks);
   pl.set<int>("Y Blocks",yblocks);
   pl.set<int>("Z Blocks",zblocks);

   panzer_stk::CubeHexMeshFactory meshFact;
   meshFact.setParameterList(Teuchos::rcpFromRef(pl));

   Teuchos::RCP<panzer_stk::STK_Interface> mesh = meshFact.buildMesh(comm);
   mesh->writeToExodus("whatish.exo");
   return mesh;
}

template <typename Intrepid2Type>
RCP<const panzer::FieldPattern> buildFieldPattern()
{
   // build a geometric pattern from a single basis
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2Type);
   RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
   return pattern;
}

// quad tests
TEUCHOS_UNIT_TEST(tCubeHexMeshDOFManager, buildTest_hex)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      stk::ParallelMachine Comm = MPI_COMM_WORLD;
   #else
      stk::ParallelMachine Comm = WHAT_TO_DO_COMM;
   #endif

   int numProcs = stk::parallel_machine_size(Comm);
   int myRank = stk::parallel_machine_rank(Comm);

   TEUCHOS_ASSERT(numProcs<=2);

   // build a geometric pattern from a single basis
   RCP<const panzer::FieldPattern> patternC1
         = buildFieldPattern<Intrepid2::Basis_HGRAD_HEX_C1_FEM<PHX::exec_space,double,double> >();

   Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildHexMesh(Comm,2,2,2,1,1,1);
   RCP<panzer::ConnManager> connManager
         = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
   RCP<panzer::DOFManager> dofManager = rcp(new panzer::DOFManager());

   TEST_EQUALITY(dofManager->getOrientationsRequired(),false);
   TEST_EQUALITY(dofManager->getConnManager(),Teuchos::null);

   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   TEST_EQUALITY(dofManager->getConnManager(),connManager);

   dofManager->addField("ux",patternC1);
   dofManager->addField("uy",patternC1);
   dofManager->addField("p",patternC1);

   std::vector<std::string> fieldOrder;
   fieldOrder.push_back("ux");
   fieldOrder.push_back("uy");
   fieldOrder.push_back("p");
   dofManager->setFieldOrder(fieldOrder);
   dofManager->buildGlobalUnknowns();
   dofManager->printFieldInformation(out);

   if(numProcs==1) {
      std::vector<panzer::GlobalOrdinal> gids_v;
      panzer::GlobalOrdinal * gids = nullptr;

      TEST_ASSERT(false);

      // element 0
      dofManager->getElementGIDs(0,gids_v);

      gids = &gids_v[0];
      TEST_EQUALITY(gids_v.size(),24);
      TEST_EQUALITY(gids[0],0);  TEST_EQUALITY(gids[1],1);   TEST_EQUALITY(gids[2],2);
      TEST_EQUALITY(gids[3],3);  TEST_EQUALITY(gids[4],4);   TEST_EQUALITY(gids[5],5);
      TEST_EQUALITY(gids[6],12); TEST_EQUALITY(gids[7],13);  TEST_EQUALITY(gids[8],14);
      TEST_EQUALITY(gids[9], 9); TEST_EQUALITY(gids[10],10); TEST_EQUALITY(gids[11],11);

      gids = &gids_v[12];
      TEST_EQUALITY(gids[0],27); TEST_EQUALITY(gids[1],28);  TEST_EQUALITY(gids[2],29);
      TEST_EQUALITY(gids[3],30); TEST_EQUALITY(gids[4],31);  TEST_EQUALITY(gids[5],32);
      TEST_EQUALITY(gids[6],39); TEST_EQUALITY(gids[7],40);  TEST_EQUALITY(gids[8],41);
      TEST_EQUALITY(gids[9],36); TEST_EQUALITY(gids[10],37); TEST_EQUALITY(gids[11],38);

      // element 6
      dofManager->getElementGIDs(mesh->elementLocalId(5),gids_v);

      gids = &gids_v[0];
      TEST_EQUALITY(gids_v.size(),24);
      TEST_EQUALITY(gids[0],27); TEST_EQUALITY(gids[1],28);  TEST_EQUALITY(gids[2],29);
      TEST_EQUALITY(gids[3],30); TEST_EQUALITY(gids[4],31);  TEST_EQUALITY(gids[5],32);
      TEST_EQUALITY(gids[6],39); TEST_EQUALITY(gids[7],40);  TEST_EQUALITY(gids[8],41);
      TEST_EQUALITY(gids[9],36); TEST_EQUALITY(gids[10],37); TEST_EQUALITY(gids[11],38);

      gids = &gids_v[12];
      TEST_EQUALITY(gids[0],54); TEST_EQUALITY(gids[1],55);  TEST_EQUALITY(gids[2],56);
      TEST_EQUALITY(gids[3],57); TEST_EQUALITY(gids[4],58);  TEST_EQUALITY(gids[5],59);
      TEST_EQUALITY(gids[6],66); TEST_EQUALITY(gids[7],67);  TEST_EQUALITY(gids[8],68);
      TEST_EQUALITY(gids[9],63); TEST_EQUALITY(gids[10],64); TEST_EQUALITY(gids[11],65);
   }
   else if(myRank==0) {
      // element 7
      const auto * gids = connManager->getConnectivity(mesh->elementLocalId(7));
      TEST_EQUALITY(connManager->getConnectivitySize(mesh->elementLocalId(7)),8);

      TEST_EQUALITY(gids[0],12); TEST_EQUALITY(gids[1],13);
      TEST_EQUALITY(gids[2],16); TEST_EQUALITY(gids[3],15);
      TEST_EQUALITY(gids[4],21); TEST_EQUALITY(gids[5],22);
      TEST_EQUALITY(gids[6],25); TEST_EQUALITY(gids[7],24);
   }
   else if(myRank==1) {
      // element 2
      const auto * gids = connManager->getConnectivity(mesh->elementLocalId(2));
      TEST_EQUALITY(connManager->getConnectivitySize(mesh->elementLocalId(2)),8);

      TEST_EQUALITY(gids[0],1); TEST_EQUALITY(gids[1],2);
      TEST_EQUALITY(gids[2],5); TEST_EQUALITY(gids[3],4);
      TEST_EQUALITY(gids[4],10); TEST_EQUALITY(gids[5],11);
      TEST_EQUALITY(gids[6],14); TEST_EQUALITY(gids[7],13);
   }

   // check that owned is_subset owned_and_ghosted
   //////////////////////////////////////////////////////////////////////////
   std::vector<panzer::GlobalOrdinal> owned, owned_and_ghosted;
   dofManager->getOwnedIndices(owned);
   dofManager->getOwnedAndGhostedIndices(owned_and_ghosted);

   if(numProcs==1) {
     TEST_EQUALITY(owned.size(),owned_and_ghosted.size());
   }
   else  {
     out << "owned size = " << owned.size() << std::endl;
     out << "owned_and_ghosted size = " << owned_and_ghosted.size() << std::endl;
     TEST_ASSERT(owned.size()<=owned_and_ghosted.size());
   }
   for(std::size_t i=0;i<owned.size();i++) {
     TEST_EQUALITY(owned[i],owned_and_ghosted[i]);
   }
   for(std::size_t i=owned.size();i<owned_and_ghosted.size();i++) {
     TEST_ASSERT(std::find(owned.begin(),owned.end(),owned_and_ghosted[i])==owned.end());
   }
}

// quad tests
TEUCHOS_UNIT_TEST(tCubeHexMeshDOFManager, buildTest_hex_face_orientations)
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
         = buildFieldPattern<Intrepid2::Basis_HDIV_HEX_I1_FEM<PHX::exec_space,double,double> >();

   RCP<panzer::ConnManager> connManager =
       Teuchos::rcp(new panzer_stk::STKConnManager(buildHexMesh(Comm,2,2,2,1,1,1)));
   RCP<panzer::DOFManager> dofManager = rcp(new panzer::DOFManager());

   dofManager->setOrientationsRequired(true);
   TEST_EQUALITY(dofManager->getOrientationsRequired(),true);
   TEST_EQUALITY(dofManager->getConnManager(),Teuchos::null);

   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   TEST_EQUALITY(dofManager->getConnManager(),connManager);

   dofManager->addField("b",patternI1);

   dofManager->buildGlobalUnknowns();
   dofManager->printFieldInformation(out);

   const std::vector<int> & b_offsets = dofManager->getGIDFieldOffsets("eblock-0_0_0",dofManager->getFieldNum("b"));

   TEST_EQUALITY(b_offsets.size(),6);

   // unfortunatly this mesh is completly uniform
   double standardO[] = { 1.0, 1.0, -1.0, -1.0, -1.0, 1.0 };
   if(myRank==0) {
      std::vector<double> orientation;

      // element 0
      dofManager->getElementOrientation(0,orientation);
      TEST_EQUALITY(orientation.size(),6);

      for(std::size_t i=0;i<6;i++)
         TEST_EQUALITY(orientation[i],standardO[i]);

      // element 1
      dofManager->getElementOrientation(1,orientation);
      TEST_EQUALITY(orientation.size(),6);

      for(std::size_t i=0;i<6;i++)
         TEST_EQUALITY(orientation[i],standardO[i]);
   }
   else if(myRank==1) {
      std::vector<double> orientation;

      // element 0
      dofManager->getElementOrientation(0,orientation);
      TEST_EQUALITY(orientation.size(),6);

      for(std::size_t i=0;i<6;i++)
         TEST_EQUALITY(orientation[i],standardO[i]);

      // element 1
      dofManager->getElementOrientation(1,orientation);
      TEST_EQUALITY(orientation.size(),6);

      for(std::size_t i=0;i<6;i++)
         TEST_EQUALITY(orientation[i],standardO[i]);
   }
}

}
