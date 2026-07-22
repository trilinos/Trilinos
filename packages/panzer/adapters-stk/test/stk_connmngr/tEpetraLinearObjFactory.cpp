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

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_Traits.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

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
RCP<const panzer::FieldPattern> buildFieldPattern()
{
   // build a geometric pattern from a single basis
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2Type);
   RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
   return pattern;
}

// quad tests
TEUCHOS_UNIT_TEST(tEpetraLinearObjFactory, buildTest_quad)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      stk::ParallelMachine Comm = MPI_COMM_WORLD;
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      stk::ParallelMachine Comm = WHAT_TO_DO_COMM;
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   Teuchos::RCP<const Teuchos::MpiComm<int> > tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

   int numProcs = stk::parallel_machine_size(Comm);
   int myRank = stk::parallel_machine_rank(Comm);

   TEUCHOS_ASSERT(numProcs<=2);

   // build a geometric pattern from a single basis
   RCP<const panzer::FieldPattern> patternC1 
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();

   RCP<panzer::ConnManager> connManager = buildQuadMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager> dofManager = rcp(new panzer::DOFManager());
   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   dofManager->addField("u",patternC1);
   dofManager->buildGlobalUnknowns();

   panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> laFactory(tComm.getConst(),dofManager);
   Teuchos::RCP<Epetra_Map> map = laFactory.getMap(0);
   Teuchos::RCP<Epetra_Map> gMap = laFactory.getGhostedMap(0);
   Teuchos::RCP<Epetra_CrsGraph> graph = laFactory.getGraph(0,0);
   Teuchos::RCP<Epetra_CrsGraph> gGraph = laFactory.getGhostedGraph(0,0);

   std::vector<panzer::GlobalOrdinal> owned,ownedAndGhosted;
   dofManager->getOwnedIndices(owned);
   dofManager->getOwnedAndGhostedIndices(ownedAndGhosted);
  
   // test maps
   {
      TEST_EQUALITY(map->NumMyElements(),(int) owned.size());
      TEST_EQUALITY(gMap->NumMyElements(),(int) ownedAndGhosted.size());

      // test indices
      for(std::size_t i=0;i<owned.size();i++) TEST_ASSERT(map->MyGID(owned[i]));
      for(std::size_t i=0;i<ownedAndGhosted.size();i++) TEST_ASSERT(gMap->MyGID(ownedAndGhosted[i]));
   }

   // test ograph
   {
      TEST_ASSERT(gGraph->Filled());

      TEST_EQUALITY(gGraph->NumMyRows(),(int) ownedAndGhosted.size());
      TEST_EQUALITY(gGraph->MaxNumIndices(),numProcs==2 ? 6 : 9);

      std::vector<int> indices(10);
      int numIndices = 0;

      // Take degree of freedom in middle of mesh: Then look at ghosted graph
      int err = gGraph->ExtractGlobalRowCopy(3 /* magic number */,10,numIndices,&indices[0]);
      TEST_EQUALITY(err,0);

      indices.resize(numIndices);
      std::sort(indices.begin(),indices.end());
      if(numProcs==2 && myRank==0) {
         TEST_EQUALITY(numIndices,6); 

         std::vector<int> compare(6);
         compare[0] = 0; compare[1] = 1; compare[2] = 2;
         compare[3] = 3; compare[4] = 4; compare[5] = 5;
         
         TEST_EQUALITY(compare.size(),indices.size());
         TEST_ASSERT(std::equal(compare.begin(),compare.end(),indices.begin()));
      }
      else if(numProcs==2 && myRank==1) {
         TEST_EQUALITY(numIndices,6); 

         std::vector<int> compare(6);
         compare[0] = 1; compare[1] = 3; compare[2] = 5;
         compare[3] = 6; compare[4] = 7; compare[5] = 8;
         
         TEST_EQUALITY(compare.size(),indices.size());
         TEST_ASSERT(std::equal(compare.begin(),compare.end(),indices.begin()));
      }
      
   }

   // test graph
   {
      TEST_ASSERT(graph->Filled());

      TEST_EQUALITY(graph->NumMyRows(),(int) owned.size());
      TEST_EQUALITY(graph->MaxNumIndices(),myRank==0 ? 9 : 6);
   }

}

}
