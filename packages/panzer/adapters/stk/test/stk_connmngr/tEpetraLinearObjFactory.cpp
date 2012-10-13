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
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_config.hpp"
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_Traits.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

typedef Intrepid::FieldContainer<double> FieldContainer;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

namespace panzer_stk {

Teuchos::RCP<panzer::ConnManager<int,int> > buildQuadMesh(stk::ParallelMachine comm,int xelmts,int yelmts,int xblocks,int yblocks)
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

template <typename IntrepidType>
RCP<const panzer::FieldPattern> buildFieldPattern()
{
   // build a geometric pattern from a single basis
   RCP<Intrepid::Basis<double,FieldContainer> > basis = rcp(new IntrepidType);
   RCP<const panzer::FieldPattern> pattern = rcp(new panzer::IntrepidFieldPattern(basis));
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

   int numProcs = stk::parallel_machine_size(Comm);
   int myRank = stk::parallel_machine_rank(Comm);

   TEUCHOS_ASSERT(numProcs<=2);

   // build a geometric pattern from a single basis
   RCP<const panzer::FieldPattern> patternC1 
         = buildFieldPattern<Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();

   RCP<panzer::ConnManager<int,int> > connManager = buildQuadMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager<int,int> > dofManager = rcp(new panzer::DOFManager<int,int>());
   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   dofManager->addField("u",patternC1);
   dofManager->buildGlobalUnknowns();

   panzer::EpetraLinearObjFactory<panzer::Traits,int> laFactory(eComm.getConst(),dofManager);
   Teuchos::RCP<Epetra_Map> map = laFactory.getMap();
   Teuchos::RCP<Epetra_Map> gMap = laFactory.getGhostedMap();
   Teuchos::RCP<Epetra_CrsGraph> graph = laFactory.getGraph();
   Teuchos::RCP<Epetra_CrsGraph> gGraph = laFactory.getGhostedGraph();

   std::vector<int> owned,ownedAndShared;
   dofManager->getOwnedIndices(owned);
   dofManager->getOwnedAndSharedIndices(ownedAndShared);
  
   // test maps
   {
      TEST_EQUALITY(map->NumMyElements(),(int) owned.size());
      TEST_EQUALITY(gMap->NumMyElements(),(int) ownedAndShared.size());

      // test indices
      for(std::size_t i=0;i<owned.size();i++) TEST_ASSERT(map->MyGID(owned[i]));
      for(std::size_t i=0;i<ownedAndShared.size();i++) TEST_ASSERT(gMap->MyGID(ownedAndShared[i]));
   }

   // test ograph
   {
      TEST_ASSERT(gGraph->Filled());

      TEST_EQUALITY(gGraph->NumMyRows(),(int) ownedAndShared.size());
      TEST_EQUALITY(gGraph->MaxNumIndices(),numProcs==2 ? 6 : 9);

      std::vector<int> indices(10);
      int numIndices = 0;

      // Take degree of freedom in middle of mesh: Then look at ghosted graph
      int err = gGraph->ExtractGlobalRowCopy(3,10,numIndices,&indices[0]);
      TEST_EQUALITY(err,0);

      indices.resize(numIndices);
      std::sort(indices.begin(),indices.end());
      if(numProcs==1) {
         TEST_EQUALITY(numIndices,9); 

         std::vector<int> compare(9);
         compare[0] = 0; compare[1] = 1; compare[2] = 2;
         compare[3] = 3; compare[4] = 4; compare[5] = 5;
         compare[3] = 6; compare[4] = 7; compare[5] = 8;
         
         TEST_EQUALITY(compare.size(),indices.size());
         TEST_ASSERT(std::equal(compare.begin(),compare.end(),indices.begin()));
      }
      else if(numProcs==2 && myRank==0) {
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
