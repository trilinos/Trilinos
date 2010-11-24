#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_EpetraLinearObjFactory.hpp"

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

template class panzer::EpetraLinearObjFactory<char>;

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

   panzer::EpetraLinearObjFactory<int> laFactory(eComm,dofManager);
   Teuchos::RCP<Epetra_Map> map = laFactory.getMap();
   Teuchos::RCP<Epetra_Map> oMap = laFactory.getOverlapMap();
   Teuchos::RCP<Epetra_CrsGraph> graph = laFactory.getGraph();
   Teuchos::RCP<Epetra_CrsGraph> oGraph = laFactory.getOverlapGraph();

   std::vector<int> owned,ownedAndShared;
   dofManager->getOwnedIndices(owned);
   dofManager->getOwnedAndSharedIndices(ownedAndShared);
  
   // test maps
   {
      TEST_EQUALITY(map->NumMyElements(),(int) owned.size());
      TEST_EQUALITY(oMap->NumMyElements(),(int) ownedAndShared.size());

      // test indices
      for(std::size_t i=0;i<owned.size();i++) TEST_ASSERT(map->MyGID(owned[i]));
      for(std::size_t i=0;i<ownedAndShared.size();i++) TEST_ASSERT(oMap->MyGID(ownedAndShared[i]));
   }

   // test ograph
   {
      TEST_ASSERT(oGraph->Filled());

      TEST_EQUALITY(oGraph->NumMyRows(),(int) ownedAndShared.size());
      TEST_EQUALITY(oGraph->MaxNumIndices(),numProcs==2 ? 6 : 9);

      std::vector<int> indices(10);
      int numIndices = 0;
      int err = oGraph->ExtractGlobalRowCopy(3,10,numIndices,&indices[0]);
      TEST_EQUALITY(err,0);
      TEST_EQUALITY(numIndices,6); 

      indices.resize(numIndices);
      std::sort(indices.begin(),indices.end());
      if(numProcs==1) {
         std::vector<int> compare(6);
         compare[0] = 0; compare[1] = 1; compare[2] = 3;
         compare[3] = 4; compare[4] = 6; compare[5] = 7;
         
         TEST_EQUALITY(compare.size(),indices.size());
         TEST_ASSERT(std::equal(compare.begin(),compare.end(),indices.begin()));
      }
      else if(numProcs==2 && myRank==0) {
         std::vector<int> compare(6);
         compare[0] = 0; compare[1] = 1; compare[2] = 2;
         compare[3] = 3; compare[4] = 4; compare[5] = 5;
         
         TEST_EQUALITY(compare.size(),indices.size());
         TEST_ASSERT(std::equal(compare.begin(),compare.end(),indices.begin()));
      }
      else if(numProcs==2 && myRank==1) {
         std::vector<int> compare(6);
         compare[0] = 1; compare[1] = 3; compare[2] = 5;
         compare[3] = 6; compare[4] = 7; compare[5] = 8;
         
         TEST_EQUALITY(compare.size(),indices.size());
         TEST_ASSERT(std::equal(compare.begin(),compare.end(),indices.begin()));
      }
      else TEUCHOS_ASSERT(false);
   }

   // test graph
   {
      TEST_ASSERT(graph->Filled());

      TEST_EQUALITY(graph->NumMyRows(),(int) owned.size());
      TEST_EQUALITY(graph->MaxNumIndices(),myRank==0 ? 9 : 6);
   }
}

}
