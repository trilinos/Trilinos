#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_AlbanySTKConnManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Albany_Quad2DSTKMeshStruct.hpp"
#include "Albany_STKDiscretization.hpp"

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

   TEUCHOS_ASSERT(numProcs<=2);

   // build a geometric pattern from a single basis
   RCP<const panzer::FieldPattern> patternC1 
         = buildFieldPattern<Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();

   RCP<panzer::ConnManager<int,int> > connManager = buildQuadMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager<int,int> > dofManager = rcp(new panzer::DOFManager<int,int>());

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
   dofManager->printFieldInformation(out);
   dofManager->buildGlobalUnknowns();

   if(numProcs==1) {
      std::vector<int> gids;

      dofManager->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),12);
      TEST_EQUALITY(gids[0],0); TEST_EQUALITY(gids[1],1); TEST_EQUALITY(gids[2],2);
      TEST_EQUALITY(gids[3],3); TEST_EQUALITY(gids[4],4); TEST_EQUALITY(gids[5],5);
      TEST_EQUALITY(gids[6],12); TEST_EQUALITY(gids[7],13); TEST_EQUALITY(gids[8],14);
      TEST_EQUALITY(gids[9], 9); TEST_EQUALITY(gids[10],10); TEST_EQUALITY(gids[11],11);
   
      dofManager->getElementGIDs(1,gids);
      TEST_EQUALITY(gids.size(),12);
      TEST_EQUALITY(gids[0],9); TEST_EQUALITY(gids[1],10); TEST_EQUALITY(gids[2],11);
      TEST_EQUALITY(gids[3],12); TEST_EQUALITY(gids[4],13); TEST_EQUALITY(gids[5],14);
      TEST_EQUALITY(gids[6],21); TEST_EQUALITY(gids[7],22); TEST_EQUALITY(gids[8],23);
      TEST_EQUALITY(gids[9],18); TEST_EQUALITY(gids[10],19); TEST_EQUALITY(gids[11],20);
   }
   else if(myRank==0) {
      std::vector<int> gids;

      dofManager->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),12);
      TEST_EQUALITY(gids[0],0); TEST_EQUALITY(gids[1],1); TEST_EQUALITY(gids[2],2);
      TEST_EQUALITY(gids[3],3); TEST_EQUALITY(gids[4],4); TEST_EQUALITY(gids[5],5);
      TEST_EQUALITY(gids[6],9); TEST_EQUALITY(gids[7],10); TEST_EQUALITY(gids[8],11);
      TEST_EQUALITY(gids[9],6); TEST_EQUALITY(gids[10],7); TEST_EQUALITY(gids[11],8);
   
      dofManager->getElementGIDs(1,gids);
      TEST_EQUALITY(gids.size(),12);
      TEST_EQUALITY(gids[0],6); TEST_EQUALITY(gids[1],7); TEST_EQUALITY(gids[2],8);
      TEST_EQUALITY(gids[3],9); TEST_EQUALITY(gids[4],10); TEST_EQUALITY(gids[5],11);
      TEST_EQUALITY(gids[6],15); TEST_EQUALITY(gids[7],16); TEST_EQUALITY(gids[8],17);
      TEST_EQUALITY(gids[9], 12); TEST_EQUALITY(gids[10],13); TEST_EQUALITY(gids[11],14);
   }
   else if(myRank==1) {
      std::vector<int> gids;

      dofManager->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),12);
      TEST_EQUALITY(gids[0],3); TEST_EQUALITY(gids[1],4); TEST_EQUALITY(gids[2],5);
      TEST_EQUALITY(gids[3],18); TEST_EQUALITY(gids[4],19); TEST_EQUALITY(gids[5],20);
      TEST_EQUALITY(gids[6],21); TEST_EQUALITY(gids[7],22); TEST_EQUALITY(gids[8],23);
      TEST_EQUALITY(gids[9],9); TEST_EQUALITY(gids[10],10); TEST_EQUALITY(gids[11],11);
   
      dofManager->getElementGIDs(1,gids);
      TEST_EQUALITY(gids.size(),12);
      TEST_EQUALITY(gids[0],9); TEST_EQUALITY(gids[1],10); TEST_EQUALITY(gids[2],11);
      TEST_EQUALITY(gids[3],21); TEST_EQUALITY(gids[4],22); TEST_EQUALITY(gids[5],23);
      TEST_EQUALITY(gids[6],24); TEST_EQUALITY(gids[7],25); TEST_EQUALITY(gids[8],26);
      TEST_EQUALITY(gids[9],15); TEST_EQUALITY(gids[10],16); TEST_EQUALITY(gids[11],17);
   }
}

// quad tests
TEUCHOS_UNIT_TEST(tSquareQuadMeshDOFManager, shared_owned_indices)
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
         = buildFieldPattern<Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();

   // build DOF manager
   RCP<panzer::ConnManager<int,int> > connManager = buildQuadMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager<int,int> > dofManager = rcp(new panzer::DOFManager<int,int>());
   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   dofManager->addField("u",patternC1);
   dofManager->buildGlobalUnknowns();

   // test GlobalNumProvider
   RCP<panzer::GlobalNumProvider<int,int> > glbNum = dofManager;

   std::vector<int> owned, ownedAndShared;
   glbNum->getOwnedIndices(owned);
   glbNum->getOwnedAndSharedIndices(ownedAndShared);

   if(numProcs==1) {
      TEST_EQUALITY(owned.size(),ownedAndShared.size());

      bool ownedCorrect = true;
      bool ownedAndSharedCorrect = true;
      std::sort(owned.begin(),owned.end());
      std::sort(ownedAndShared.begin(),ownedAndShared.end());
      for(std::size_t i=0;i<owned.size();i++) {
         ownedCorrect &= (owned[i] == (int) i);
         ownedAndSharedCorrect &= (ownedAndShared[i] == (int) i);
      }

      TEST_ASSERT(ownedCorrect);
      TEST_ASSERT(ownedAndSharedCorrect);
   }
   else if(myRank==0 && numProcs==2) {
      TEST_EQUALITY(owned.size(),6);
      TEST_EQUALITY(ownedAndShared.size(),6);
      bool ownedCorrect = true;
      bool ownedAndSharedCorrect = true;

      std::sort(owned.begin(),owned.end());
      for(std::size_t i=0;i<owned.size();i++) {
         ownedCorrect &= (owned[i] == (int) i);
      }
      TEST_ASSERT(ownedCorrect);

      std::sort(ownedAndShared.begin(),ownedAndShared.end());
      for(std::size_t i=0;i<ownedAndShared.size();i++) {
         ownedAndSharedCorrect &= (ownedAndShared[i] == (int) i);
      }
      TEST_ASSERT(ownedAndSharedCorrect);
   }
   else if(myRank==1 && numProcs==2) {
      TEST_EQUALITY(owned.size(),3);
      TEST_EQUALITY(ownedAndShared.size(),6);
      bool ownedCorrect = true;
      bool ownedAndSharedCorrect = true;

      std::sort(owned.begin(),owned.end());
      for(std::size_t i=0;i<owned.size();i++) {
         ownedCorrect &= (owned[i] == (int) i+6);
      }
      TEST_ASSERT(ownedCorrect);

      std::sort(ownedAndShared.begin(),ownedAndShared.end());
      for(std::size_t i=0;i<3;i++) {
         ownedAndSharedCorrect &= (ownedAndShared[i] == (int) (2*i+1));
      }
      for(std::size_t i=0;i<3;i++) {
         ownedAndSharedCorrect &= (ownedAndShared[i+3] == (int) i+6);
      }
      TEST_ASSERT(ownedAndSharedCorrect);
   }
   else 
      TEUCHOS_ASSERT(false);
}

}
