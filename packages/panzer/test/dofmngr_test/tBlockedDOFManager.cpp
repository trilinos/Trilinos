#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <string>
#include <iostream>
#include <vector>
#include <set>

#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_PauseToAttach.hpp"

#include "UnitTest_ConnManager.hpp"

// include some intrepid basis functions
// 2D basis 
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"

#include "Intrepid_FieldContainer.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

typedef Intrepid::FieldContainer<double> FieldContainer;

namespace panzer {

template <typename IntrepidType>
Teuchos::RCP<const panzer::FieldPattern> buildFieldPattern()
{
   // build a geometric pattern from a single basis
   Teuchos::RCP<Intrepid::Basis<double,FieldContainer> > basis = rcp(new IntrepidType);
   Teuchos::RCP<const panzer::FieldPattern> pattern = rcp(new panzer::IntrepidFieldPattern(basis));
   return pattern;
}

// this just excercises a bunch of functions
TEUCHOS_UNIT_TEST(tBlockedDOFManager_SimpleTests,assortedTests)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<ConnManager<short,int> > connManager = rcp(new unit_test::ConnManager(myRank,numProc));
   BlockedDOFManager<short,int> dofManager; 
   dofManager.setConnManager(connManager,MPI_COMM_WORLD);

   TEST_ASSERT(dofManager.getComm()!=Teuchos::null);
   TEST_EQUALITY(dofManager.getConnManager(),connManager);

   std::vector<std::string> eBlocks;
   dofManager.getElementBlockIds(eBlocks);
   std::sort(eBlocks.begin(),eBlocks.end());
   TEST_EQUALITY(eBlocks.size(),3);
   TEST_EQUALITY(eBlocks[0],"block_0");
   TEST_EQUALITY(eBlocks[1],"block_1");
   TEST_EQUALITY(eBlocks[2],"block_2");
   TEST_EQUALITY(dofManager.getNumFieldBlocks(),1); // if no field order is set this defaults to 1!

   std::vector<std::vector<std::string> > fieldOrder(3),fo_ut;
   fieldOrder[0].push_back("Ux");
   fieldOrder[0].push_back("Uy");
   fieldOrder[1].push_back("P");
   fieldOrder[2].push_back("rho");
   fieldOrder[2].push_back("T");
   dofManager.setFieldOrder(fieldOrder);
   dofManager.getFieldOrder(fo_ut);
   TEST_ASSERT(fieldOrder==fo_ut);
   TEST_EQUALITY(dofManager.getNumFieldBlocks(),3); // 

   TEST_ASSERT(dofManager.getElementBlock("block_0")==connManager->getElementBlock("block_0"));
   TEST_ASSERT(dofManager.getElementBlock("block_1")==connManager->getElementBlock("block_1"));
   TEST_ASSERT(dofManager.getElementBlock("block_2")==connManager->getElementBlock("block_2"));
}

TEUCHOS_UNIT_TEST(tBlockedDOFManager_SimpleTests,registerFields)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<ConnManager<short,int> > connManger = rcp(new unit_test::ConnManager(myRank,numProc));
   BlockedDOFManager<short,int> dofManager; 
   dofManager.setConnManager(connManger,MPI_COMM_WORLD);

   TEST_EQUALITY(dofManager.getMaxSubFieldNumber(),-1);

   RCP<const panzer::FieldPattern> patternC1 
         = buildFieldPattern<Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();

   dofManager.addField("T",patternC1); // add it to all three blocks

   dofManager.addField("block_0","Ux",patternC1);
   dofManager.addField("block_0","Uy",patternC1);
   dofManager.addField("block_0","P",patternC1);

   dofManager.addField("block_2","rho",patternC1);

   TEST_ASSERT(!dofManager.fieldsRegistered());
   TEST_EQUALITY(dofManager.getNumFields(),5);

   TEST_EQUALITY(dofManager.getFieldPattern("block_0","T"),patternC1);
   TEST_EQUALITY(dofManager.getFieldPattern("block_1","T"),patternC1);
   TEST_EQUALITY(dofManager.getFieldPattern("block_2","T"),patternC1);
   TEST_EQUALITY(dofManager.getFieldPattern("block_0","Ux"),patternC1);
   TEST_EQUALITY(dofManager.getFieldPattern("block_1","Ux"),Teuchos::null);
   TEST_EQUALITY(dofManager.getFieldPattern("block_2","Ux"),Teuchos::null);

   TEST_ASSERT(dofManager.fieldInBlock("T","block_0"));
   TEST_ASSERT(dofManager.fieldInBlock("T","block_1"));
   TEST_ASSERT(dofManager.fieldInBlock("T","block_2"));
   TEST_ASSERT(dofManager.fieldInBlock("Ux","block_0"));
   TEST_ASSERT(!dofManager.fieldInBlock("Ux","block_2"));
   TEST_ASSERT(dofManager.fieldInBlock("Uy","block_0"));
   TEST_ASSERT(!dofManager.fieldInBlock("Uy","block_1"));
   TEST_ASSERT(dofManager.fieldInBlock("P","block_0"));
   TEST_ASSERT(!dofManager.fieldInBlock("P","block_1"));
   TEST_ASSERT(!dofManager.fieldInBlock("rho","block_1"));
   TEST_ASSERT(dofManager.fieldInBlock("rho","block_2"));

   // set up a blocking structure
   std::vector<std::vector<std::string> > fieldOrder(3);
   fieldOrder[0].push_back("Ux");
   fieldOrder[0].push_back("Uy");
   fieldOrder[1].push_back("P");
   fieldOrder[2].push_back("rho");
   fieldOrder[2].push_back("T");
   dofManager.setFieldOrder(fieldOrder);

   dofManager.registerFields();
   TEST_ASSERT(dofManager.fieldsRegistered());
   const std::vector<RCP<panzer::DOFManager<short,int> > > & subManagers = 
         dofManager.getFieldDOFManagers();
   TEST_EQUALITY(subManagers.size(),fieldOrder.size());

   TEST_EQUALITY(subManagers[0]->getNumFields(),2);
   TEST_EQUALITY(subManagers[0]->getFieldPattern("block_0","Ux"),patternC1);
   TEST_EQUALITY(subManagers[0]->getFieldPattern("block_0","Uy"),patternC1);
   TEST_EQUALITY(subManagers[0]->getFieldPattern("block_1","Uy"),Teuchos::null);
   TEST_EQUALITY(subManagers[0]->getFieldPattern("block_1","T"),Teuchos::null);

   TEST_EQUALITY(subManagers[1]->getNumFields(),1);
   TEST_EQUALITY(subManagers[1]->getFieldPattern("block_0","P"),patternC1);
   TEST_EQUALITY(subManagers[1]->getFieldPattern("block_1","T"),Teuchos::null);

   TEST_EQUALITY(subManagers[2]->getNumFields(),2);
   TEST_EQUALITY(subManagers[2]->getFieldPattern("block_0","T"),patternC1);
   TEST_EQUALITY(subManagers[2]->getFieldPattern("block_1","T"),patternC1);
   TEST_EQUALITY(subManagers[2]->getFieldPattern("block_2","T"),patternC1);
   TEST_EQUALITY(subManagers[2]->getFieldPattern("block_0","rho"),Teuchos::null);
   TEST_EQUALITY(subManagers[2]->getFieldPattern("block_1","rho"),Teuchos::null);
   TEST_EQUALITY(subManagers[2]->getFieldPattern("block_2","rho"),patternC1);

   // test field numbers, should be based on a field block index * largest field
   // number+1 (in this case the largest field number is 1...hence this size for
   // the field blocks is 2
   TEST_EQUALITY(dofManager.getMaxSubFieldNumber(),1);
   TEST_EQUALITY(subManagers[0]->getFieldNum("Ux")+0*2,  dofManager.getFieldNum("Ux"));
   TEST_EQUALITY(subManagers[0]->getFieldNum("Uy")+0*2,  dofManager.getFieldNum("Uy"));
   TEST_EQUALITY(subManagers[1]->getFieldNum("P")+1*2,   dofManager.getFieldNum("P"));
   TEST_EQUALITY(subManagers[2]->getFieldNum("rho")+2*2, dofManager.getFieldNum("rho"));
   TEST_EQUALITY(subManagers[2]->getFieldNum("T")+2*2,   dofManager.getFieldNum("T"));

   // check field blocks, notice we are copying here
   std::vector<int> blk0fn = dofManager.getBlockFieldNumbers("block_0"); std::sort(blk0fn.begin(),blk0fn.end());
   std::vector<int> blk1fn = dofManager.getBlockFieldNumbers("block_1"); std::sort(blk1fn.begin(),blk1fn.end());
   std::vector<int> blk2fn = dofManager.getBlockFieldNumbers("block_2"); std::sort(blk2fn.begin(),blk2fn.end());

   TEST_EQUALITY(blk0fn.size(),4);
   TEST_EQUALITY(blk1fn.size(),1);
   TEST_EQUALITY(blk2fn.size(),2);

   TEST_EQUALITY(blk0fn[0],0);
   TEST_EQUALITY(blk0fn[1],1);
   TEST_EQUALITY(blk0fn[2],2);
   TEST_EQUALITY(blk0fn[3],dofManager.getFieldNum("T"));

   TEST_EQUALITY(blk1fn[0],dofManager.getFieldNum("T"));

   TEST_EQUALITY(blk2fn[0],4);
   TEST_EQUALITY(blk2fn[1],5);
}

TEUCHOS_UNIT_TEST(tBlockedDOFManager_SimpleTests,buildGlobalUnknowns)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<ConnManager<short,int> > connManger = rcp(new unit_test::ConnManager(myRank,numProc));
   BlockedDOFManager<short,int> dofManager; 
   dofManager.setConnManager(connManger,MPI_COMM_WORLD);

   TEST_EQUALITY(dofManager.getMaxSubFieldNumber(),-1);

   RCP<const panzer::FieldPattern> patternC1 
         = buildFieldPattern<Intrepid::Basis_HGRAD_QUAD_C1_FEM<double,FieldContainer> >();

   dofManager.addField("T",patternC1); // add it to all three blocks
   dofManager.addField("block_0","Ux", patternC1);
   dofManager.addField("block_0","Uy", patternC1);
   dofManager.addField("block_0","P",  patternC1);
   dofManager.addField("block_2","rho",patternC1);

   // set up a blocking structure
   std::vector<std::vector<std::string> > fieldOrder(3);
   fieldOrder[0].push_back("Ux");
   fieldOrder[0].push_back("Uy");
   fieldOrder[1].push_back("P");
   fieldOrder[2].push_back("rho");
   fieldOrder[2].push_back("T");
   dofManager.setFieldOrder(fieldOrder);

   dofManager.buildGlobalUnknowns();

   std::vector<BlockedDOFManager<short,int>::GlobalOrdinal> ownedAndShared, owned;
   dofManager.getOwnedAndSharedIndices(ownedAndShared);
   dofManager.getOwnedIndices(owned);

   for(std::size_t i=0;i<ownedAndShared.size();i++)
      out << "(" << ownedAndShared[i].first << ","
                 << ownedAndShared[i].second << "), ";
   out << std::endl;

   const std::vector<RCP<panzer::DOFManager<short,int> > > & subManagers = 
         dofManager.getFieldDOFManagers();
   for(std::size_t i=0;i<subManagers.size();i++) {
      std::stringstream ss;
 
      subManagers[i]->printFieldInformation(ss);

      out << "******************************************************************" 
          << ss.str() 
          << "******************************************************************" 
          << std::endl;
   }

   if(myRank==0)
   {  TEST_EQUALITY(ownedAndShared.size(),39); }
   else
   {  TEST_EQUALITY(ownedAndShared.size(),30); }
}

TEUCHOS_UNIT_TEST(tBlockedDOFManager_SimpleTests,validFieldOrder)
{
   BlockedDOFManager<int,int> dofManager; 

   std::set<std::string> validFields;
   validFields.insert("horse");
   validFields.insert("cat");
   validFields.insert("monkey");
   validFields.insert("dog");
    
   {
      std::vector<std::vector<std::string> > order(1);
      order[0].push_back("cat");
      order[0].push_back("dog");
      order[0].push_back("horse");
      order[0].push_back("monkey");
   
      dofManager.validFieldOrder(order,validFields);
      TEST_ASSERT(dofManager.validFieldOrder(order,validFields));
   }

   {
      std::vector<std::vector<std::string> > order(1);
      order[0].push_back("cat");
      order[0].push_back("horse");
      order[0].push_back("monkey");
 
      TEST_ASSERT(!dofManager.validFieldOrder(order,validFields));
   }

   {
      std::vector<std::vector<std::string> > order(1);
      order[0].push_back("cat");
      order[0].push_back("dog");
      order[0].push_back("horse");
      order[0].push_back("monkey");
      order[0].push_back("monkey");
 
      TEST_ASSERT(!dofManager.validFieldOrder(order,validFields));
   }

   {
      std::vector<std::vector<std::string> > order(1);
      order[0].push_back("cat");
      order[0].push_back("dog");
      order[0].push_back("horse");
      order[0].push_back("tank");
 
      TEST_ASSERT(!dofManager.validFieldOrder(order,validFields));
   }

   {
      std::vector<std::vector<std::string> > order(2);
      order[0].push_back("cat");
      order[0].push_back("dog");
      order[0].push_back("horse");
      order[1].push_back("monkey");
 
      TEST_ASSERT(dofManager.validFieldOrder(order,validFields));
   }

   {
      std::vector<std::vector<std::string> > order(2);
      order[0].push_back("cat");
      order[0].push_back("dog");
      order[0].push_back("horse");
      order[1].push_back("cat");
      order[1].push_back("monkey");
 
      TEST_ASSERT(!dofManager.validFieldOrder(order,validFields));
   }

   {
      std::vector<std::vector<std::string> > order(2);
      order[0].push_back("dog");
      order[0].push_back("horse");
      order[1].push_back("monkey");
 
      TEST_ASSERT(!dofManager.validFieldOrder(order,validFields));
   }

   {
      std::vector<std::vector<std::string> > order(2);
      order[0].push_back("dog");
      order[0].push_back("horse");
      order[1].push_back("monkey");
      order[1].push_back("Zebra");
 
      TEST_ASSERT(!dofManager.validFieldOrder(order,validFields));
   }
}

}
