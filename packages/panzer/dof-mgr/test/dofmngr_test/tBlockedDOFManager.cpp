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
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"

#include "Kokkos_DynRankView.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

namespace panzer {

template <typename Intrepid2Type>
Teuchos::RCP<const panzer::FieldPattern> buildFieldPattern()
{
   // build a geometric pattern from a single basis
   Teuchos::RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2Type);
   Teuchos::RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
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

   // panzer::pauseToAttach();

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<ConnManager<int,int> > connManager = rcp(new unit_test::ConnManager<int>(myRank,numProc));
   BlockedDOFManager<int,int> dofManager; 
   dofManager.setUseDOFManagerFEI(false);
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

   RCP<ConnManager<int,int> > connManger = rcp(new unit_test::ConnManager<int>(myRank,numProc));
   BlockedDOFManager<int,int> dofManager; 
   dofManager.setUseDOFManagerFEI(false);
   dofManager.setConnManager(connManger,MPI_COMM_WORLD);

   TEST_EQUALITY(dofManager.getMaxSubFieldNumber(),-1);

   RCP<const panzer::FieldPattern> patternC1 
     = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();

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

   dofManager.registerFields(true);
   TEST_ASSERT(dofManager.fieldsRegistered());
   const std::vector<RCP<panzer::UniqueGlobalIndexer<int,int> > > & subManagers = 
         dofManager.getFieldDOFManagers();
   TEST_EQUALITY(subManagers.size(),fieldOrder.size());

   typedef panzer::DOFManager<int,int> DOFManager;

   TEST_EQUALITY(subManagers[0]->getNumFields(),2);
   TEST_EQUALITY(rcp_dynamic_cast<DOFManager>(subManagers[0])->getFieldPattern("block_0","Ux"),patternC1);
   TEST_EQUALITY(rcp_dynamic_cast<DOFManager>(subManagers[0])->getFieldPattern("block_0","Uy"),patternC1);
   TEST_EQUALITY(rcp_dynamic_cast<DOFManager>(subManagers[0])->getFieldPattern("block_1","Uy"),Teuchos::null);
   TEST_EQUALITY(rcp_dynamic_cast<DOFManager>(subManagers[0])->getFieldPattern("block_1","T"),Teuchos::null);

   TEST_EQUALITY(subManagers[1]->getNumFields(),1);
   TEST_EQUALITY(rcp_dynamic_cast<DOFManager>(subManagers[1])->getFieldPattern("block_0","P"),patternC1);
   TEST_EQUALITY(rcp_dynamic_cast<DOFManager>(subManagers[1])->getFieldPattern("block_1","T"),Teuchos::null);

   TEST_EQUALITY(subManagers[2]->getNumFields(),2);
   TEST_EQUALITY(rcp_dynamic_cast<DOFManager>(subManagers[2])->getFieldPattern("block_0","T"),patternC1);
   TEST_EQUALITY(rcp_dynamic_cast<DOFManager>(subManagers[2])->getFieldPattern("block_1","T"),patternC1);
   TEST_EQUALITY(rcp_dynamic_cast<DOFManager>(subManagers[2])->getFieldPattern("block_2","T"),patternC1);
   TEST_EQUALITY(rcp_dynamic_cast<DOFManager>(subManagers[2])->getFieldPattern("block_0","rho"),Teuchos::null);
   TEST_EQUALITY(rcp_dynamic_cast<DOFManager>(subManagers[2])->getFieldPattern("block_1","rho"),Teuchos::null);
   TEST_EQUALITY(rcp_dynamic_cast<DOFManager>(subManagers[2])->getFieldPattern("block_2","rho"),patternC1);

   // test field numbers, should be based on a field block index * largest field
   // number+1 (in this case the largest field number is 1...hence this size for
   // the field blocks is 2
   TEST_EQUALITY(dofManager.getMaxSubFieldNumber(),1);
   TEST_EQUALITY(subManagers[0]->getFieldNum("Ux")+0*2,  dofManager.getFieldNum("Ux"));
   TEST_EQUALITY(subManagers[0]->getFieldNum("Uy")+0*2,  dofManager.getFieldNum("Uy"));
   TEST_EQUALITY(subManagers[1]->getFieldNum("P")+1*2,   dofManager.getFieldNum("P"));
   TEST_EQUALITY(subManagers[2]->getFieldNum("rho")+2*2, dofManager.getFieldNum("rho"));
   TEST_EQUALITY(subManagers[2]->getFieldNum("T")+2*2,   dofManager.getFieldNum("T"));

   // check field order of sub managers
   for(int i=0;i<3;i++) { 
      std::vector<std::string> subFieldOrder; 
      subManagers[i]->getFieldOrder(subFieldOrder); 
      TEST_ASSERT(subFieldOrder==fieldOrder[i]); 
   }

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

   // panzer::pauseToAttach();

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<ConnManager<int,int> > connManger = rcp(new unit_test::ConnManager<int>(myRank,numProc));
   BlockedDOFManager<int,int> dofManager; 
   dofManager.setUseDOFManagerFEI(false);
   dofManager.setConnManager(connManger,MPI_COMM_WORLD);

   TEST_EQUALITY(dofManager.getMaxSubFieldNumber(),-1);

   RCP<const panzer::FieldPattern> patternC1 
     = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();

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
   if(myRank==0)
      dofManager.printFieldInformation(out);

   TEST_ASSERT(dofManager.getGeometricFieldPattern()!=Teuchos::null);

   std::vector<BlockedDOFManager<int,int>::GlobalOrdinal> ownedAndGhosted, owned;
   std::vector<bool> ownedAndGhosted_bool, owned_bool;
   dofManager.getOwnedAndGhostedIndices(ownedAndGhosted);
   dofManager.getOwnedIndices(owned);
/*
   if(myRank==0)
   { TEST_EQUALITY(ownedAndGhosted.size(),39); }
   else
   { TEST_EQUALITY(ownedAndGhosted.size(),30); }
*/

   int sum = 0;
   int mySize = (int) owned.size();
   eComm->SumAll(&mySize,&sum,1);
   TEST_EQUALITY(sum,51);

   // give it a shuffle to make it interesting
   std::random_shuffle(owned.begin(),owned.end());
   std::random_shuffle(ownedAndGhosted.begin(),ownedAndGhosted.end());
   dofManager.ownedIndices(owned,owned_bool);
   dofManager.ownedIndices(ownedAndGhosted,ownedAndGhosted_bool);

   bool ownedCheck = true;
   for(std::size_t i=0;i<owned_bool.size();i++) 
      ownedCheck &= owned_bool[i];
   TEST_ASSERT(ownedCheck);

   ownedCheck = true;
   for(std::size_t i=0;i<ownedAndGhosted_bool.size();i++) {
      bool isOwned = std::find(owned.begin(),owned.end(),ownedAndGhosted[i])!=owned.end();

      ownedCheck &= (isOwned==ownedAndGhosted_bool[i]);
   }
   TEST_ASSERT(ownedCheck);

   // WARNING: at this point we assume unknowns are unique! WRONG THING TO DO!

   {
      int fb0 = dofManager.getBlockGIDOffset("block_0",0);
      int fb1 = dofManager.getBlockGIDOffset("block_0",1);
      int fb2 = dofManager.getBlockGIDOffset("block_0",2);

      TEST_EQUALITY(fb0,0);
      TEST_EQUALITY(fb1,8);
      TEST_EQUALITY(fb2,12);
   }
   {
      int fb0 = dofManager.getBlockGIDOffset("block_1",0);
      int fb1 = dofManager.getBlockGIDOffset("block_1",1);
      int fb2 = dofManager.getBlockGIDOffset("block_1",2);

      TEST_EQUALITY(fb0,0);
      TEST_EQUALITY(fb1,0);
      TEST_EQUALITY(fb2,0);
   }
   {
      int fb0 = dofManager.getBlockGIDOffset("block_2",0);
      int fb1 = dofManager.getBlockGIDOffset("block_2",1);
      int fb2 = dofManager.getBlockGIDOffset("block_2",2);

      TEST_EQUALITY(fb0,0);
      TEST_EQUALITY(fb1,0);
      TEST_EQUALITY(fb2,0);
   }

}

TEUCHOS_UNIT_TEST(tBlockedDOFManager_SimpleTests,getElement_gids_fieldoffsets)
{

   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   // panzer::pauseToAttach();

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<ConnManager<int,int> > connManger = rcp(new unit_test::ConnManager<int>(myRank,numProc));
   BlockedDOFManager<int,int> dofManager; 
   dofManager.setUseDOFManagerFEI(false);
   dofManager.setConnManager(connManger,MPI_COMM_WORLD);

   TEST_EQUALITY(dofManager.getMaxSubFieldNumber(),-1);

   RCP<const panzer::FieldPattern> patternC1 
     = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();

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

   // check from element block 0
   std::vector<std::pair<int,int> > gids;
   if(myRank==0) {
      // std::vector<std::pair<int,int> > gids;
      dofManager.getElementGIDs(1,gids);

      TEST_EQUALITY(gids.size(),16);
      for(std::size_t i=0; i<8;i++) TEST_EQUALITY(gids[i].first,0);
      for(std::size_t i=8; i<12;i++) TEST_EQUALITY(gids[i].first,1);
      for(std::size_t i=12;i<16;i++) TEST_EQUALITY(gids[i].first,2);
   }
   else if(myRank==1) {
      // std::vector<std::pair<int,int> > gids;
      dofManager.getElementGIDs(1,gids);

      TEST_EQUALITY(gids.size(),16);
      for(std::size_t i=0; i<8;i++) TEST_EQUALITY(gids[i].first,0);
      for(std::size_t i=8; i<12;i++) TEST_EQUALITY(gids[i].first,1);
      for(std::size_t i=12;i<16;i++) TEST_EQUALITY(gids[i].first,2);
   }

   // check from element block 1 
   if(myRank==0) {
      // std::vector<std::pair<int,int> > gids;
      dofManager.getElementGIDs(4,gids);

      TEST_EQUALITY(gids.size(),4);
      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(gids[i+0].first,2);
   }
   else if(myRank==1) {
      // std::vector<std::pair<int,int> > gids;
      dofManager.getElementGIDs(3,gids);

      TEST_EQUALITY(gids.size(),4);
      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(gids[i+0].first,2);
   }

   // check from element block 2 
   if(myRank==0) {
      // std::vector<std::pair<int,int> > gids;
      dofManager.getElementGIDs(3,gids);

      TEST_EQUALITY(gids.size(),8);
      for(std::size_t i=0;i<4;i++) {
         TEST_EQUALITY(gids[i+0].first,2);
         TEST_EQUALITY(gids[i+1].first,2);
      }
   }

   // WARNING: More full tests of actutal DOFs are required
   //          however, since its built on the DOFManager we
   //          should be OK.

   {
      const std::vector<int> *p1=0,*p2=0;

      // test that pointers are the same
      p1 = &dofManager.getGIDFieldOffsets("block_0",dofManager.getFieldNum("Ux"));
      p2 = &dofManager.getGIDFieldOffsets("block_0",dofManager.getFieldNum("Ux"));

      TEST_EQUALITY(p1,p2);
   }

   // test getGIDFieldOffsets
   {
      const std::vector<int> * vec = 0;
   
      // block 0
      vec = &dofManager.getGIDFieldOffsets("block_0",dofManager.getFieldNum("Ux"));
      TEST_EQUALITY(vec->size(),4);
      TEST_EQUALITY((*vec)[0],0); TEST_EQUALITY((*vec)[1],2); TEST_EQUALITY((*vec)[2],4); TEST_EQUALITY((*vec)[3],6);
      vec = &dofManager.getGIDFieldOffsets("block_0",dofManager.getFieldNum("P"));
      TEST_EQUALITY(vec->size(),4);
      TEST_EQUALITY((*vec)[0],8+0); TEST_EQUALITY((*vec)[1],8+1); TEST_EQUALITY((*vec)[2],8+2); TEST_EQUALITY((*vec)[3],8+3);
      vec = &dofManager.getGIDFieldOffsets("block_0",dofManager.getFieldNum("T"));
      TEST_EQUALITY(vec->size(),4);
      TEST_EQUALITY((*vec)[0],12+0); TEST_EQUALITY((*vec)[1],12+1); TEST_EQUALITY((*vec)[2],12+2); TEST_EQUALITY((*vec)[3],12+3);
   
      // block 1
      vec = &dofManager.getGIDFieldOffsets("block_1",dofManager.getFieldNum("T"));
      TEST_EQUALITY(vec->size(),4);
      TEST_EQUALITY((*vec)[0],0); TEST_EQUALITY((*vec)[1],1); TEST_EQUALITY((*vec)[2],2); TEST_EQUALITY((*vec)[3],3);
   
      // block 2
      vec = &dofManager.getGIDFieldOffsets("block_2",dofManager.getFieldNum("rho"));
      TEST_EQUALITY(vec->size(),4);
      TEST_EQUALITY((*vec)[0],0); TEST_EQUALITY((*vec)[1],2); TEST_EQUALITY((*vec)[2],4); TEST_EQUALITY((*vec)[3],6);
      vec = &dofManager.getGIDFieldOffsets("block_2",dofManager.getFieldNum("T"));
      TEST_EQUALITY(vec->size(),4);
      TEST_EQUALITY((*vec)[0],1); TEST_EQUALITY((*vec)[1],3); TEST_EQUALITY((*vec)[2],5); TEST_EQUALITY((*vec)[3],7);
   }

   // test getGIDFieldOffsets_closure
   {
      const std::pair<std::vector<int>,std::vector<int> > *p1=0,*p2=0;

      // test that pointers are the same
      p1 = &dofManager.getGIDFieldOffsets_closure("block_0",dofManager.getFieldNum("Ux"),1,0);
      p2 = &dofManager.getGIDFieldOffsets_closure("block_0",dofManager.getFieldNum("Ux"),1,0);
      TEST_EQUALITY(p1,p2);

      p1 = &dofManager.getGIDFieldOffsets_closure("block_0",dofManager.getFieldNum("Ux"),1,3);
      TEST_ASSERT(p1!=p2);

      p1 = &dofManager.getGIDFieldOffsets_closure("block_0",dofManager.getFieldNum("Ux"),1,3);
      p2 = &dofManager.getGIDFieldOffsets_closure("block_0",dofManager.getFieldNum("Ux"),1,3);
      TEST_EQUALITY(p1,p2);
   }

   { 
      const std::pair<std::vector<int>,std::vector<int> > * vec = 0;
      const std::pair<std::vector<int>,std::vector<int> > * sub_vec = 0;

      Teuchos::RCP<const UniqueGlobalIndexer<int,int> > subManager;
   
      // block 0
      subManager = dofManager.getFieldDOFManagers()[2];
      vec = &dofManager.getGIDFieldOffsets_closure("block_2",dofManager.getFieldNum("T"),1,0);
      sub_vec = &subManager->getGIDFieldOffsets_closure("block_2",subManager->getFieldNum("T"),1,0);
      TEST_EQUALITY(vec->first.size(),sub_vec->first.size());
      TEST_EQUALITY(vec->second.size(),sub_vec->second.size());
      for(std::size_t i=0;i<vec->second.size();i++) 
         TEST_EQUALITY(vec->second[i],sub_vec->second[i]);
      for(std::size_t i=0;i<vec->first.size();i++) 
         TEST_EQUALITY(vec->first[i],dofManager.getBlockGIDOffset("block_2",2)+sub_vec->first[i]);
   }

}

TEUCHOS_UNIT_TEST(tBlockedDOFManager_SimpleTests,validFieldOrder)
{

   BlockedDOFManager<int,int> dofManager; 
   dofManager.setUseDOFManagerFEI(false);

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

TEUCHOS_UNIT_TEST(tBlockedDOFManager,mergetests)
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

   RCP<ConnManager<int,int> > connManager = rcp(new unit_test::ConnManager<int>(myRank,numProc));

   RCP<const panzer::FieldPattern> patternC1 
     = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();

   // Setup two DOF managers that will correspond to the blocks of the 
   // system.
   /////////////////////////////////////////////////////////////////////////
   DOFManager<int,int> dofManager[2]; 

   dofManager[0].setConnManager(connManager,MPI_COMM_WORLD);
   dofManager[0].addField("T",patternC1); // add it to all three blocks
   dofManager[0].addField("block_0","Ux",patternC1);
   dofManager[0].addField("block_0","Uy",patternC1);
   dofManager[0].addField("block_0","P",patternC1);
   dofManager[0].addField("block_2","rho",patternC1);
   dofManager[0].buildGlobalUnknowns();

   dofManager[1].setConnManager(connManager,MPI_COMM_WORLD);
   dofManager[1].addField("x",patternC1);
   dofManager[1].buildGlobalUnknowns(dofManager[0].getGeometricFieldPattern());

   std::vector<std::vector<std::string> > fieldOrder(2);
   dofManager[0].getFieldOrder(fieldOrder[0]);
   dofManager[1].getFieldOrder(fieldOrder[1]);

   std::vector<RCP<panzer::UniqueGlobalIndexer<int,int> > > ugi_vector;
   ugi_vector.push_back(Teuchos::rcpFromRef(dofManager[0]));
   ugi_vector.push_back(Teuchos::rcpFromRef(dofManager[1]));

   // Setup the BlockedDOFManager using the previously constructed
   // DOF managers
   ////////////////////////////////////////////////////////
   BlockedDOFManager<int,int> blkDofManager; 
   blkDofManager.setUseDOFManagerFEI(false);
   blkDofManager.setConnManager(connManager,MPI_COMM_WORLD);
   blkDofManager.setFieldOrder(fieldOrder);
   blkDofManager.buildGlobalUnknowns(ugi_vector);

   // Test to make sure the BlockedDOFManager understands all the fields
   // it has.
   ////////////////////////////////////////////////////////

   TEST_ASSERT(blkDofManager.fieldsRegistered());
   TEST_EQUALITY(blkDofManager.getNumFields(),6);

   TEST_EQUALITY(blkDofManager.getFieldPattern("block_0","T"),patternC1);
   TEST_EQUALITY(blkDofManager.getFieldPattern("block_1","T"),patternC1);
   TEST_EQUALITY(blkDofManager.getFieldPattern("block_2","T"),patternC1);
   TEST_EQUALITY(blkDofManager.getFieldPattern("block_0","Ux"),patternC1);
   TEST_EQUALITY(blkDofManager.getFieldPattern("block_1","Ux"),Teuchos::null);
   TEST_EQUALITY(blkDofManager.getFieldPattern("block_2","Ux"),Teuchos::null);
   TEST_EQUALITY(blkDofManager.getFieldPattern("block_0","x"),patternC1);
   TEST_EQUALITY(blkDofManager.getFieldPattern("block_1","x"),patternC1);
   TEST_EQUALITY(blkDofManager.getFieldPattern("block_2","x"),patternC1);

   TEST_ASSERT(blkDofManager.fieldInBlock("T","block_0"));
   TEST_ASSERT(blkDofManager.fieldInBlock("T","block_1"));
   TEST_ASSERT(blkDofManager.fieldInBlock("T","block_2"));
   TEST_ASSERT(blkDofManager.fieldInBlock("Ux","block_0"));
   TEST_ASSERT(!blkDofManager.fieldInBlock("Ux","block_2"));
   TEST_ASSERT(blkDofManager.fieldInBlock("Uy","block_0"));
   TEST_ASSERT(!blkDofManager.fieldInBlock("Uy","block_1"));
   TEST_ASSERT(blkDofManager.fieldInBlock("P","block_0"));
   TEST_ASSERT(!blkDofManager.fieldInBlock("P","block_1"));
   TEST_ASSERT(!blkDofManager.fieldInBlock("rho","block_1"));
   TEST_ASSERT(blkDofManager.fieldInBlock("rho","block_2"));
   TEST_ASSERT(blkDofManager.fieldInBlock("x","block_0"));
   TEST_ASSERT(blkDofManager.fieldInBlock("x","block_1"));
   TEST_ASSERT(blkDofManager.fieldInBlock("x","block_2"));

   // verify the elements are correct
   /////////////////////////////////////////////////////////

   std::vector<std::string> eBlocks;
   blkDofManager.getElementBlockIds(eBlocks);
   for(std::size_t eb=0;eb<eBlocks.size();eb++) {
     std::string block_name = eBlocks[eb];

     const std::vector<int> & blk_elements = blkDofManager.getElementBlock(block_name);
     const std::vector<int> & elements_0   = dofManager[0].getElementBlock(block_name);
     const std::vector<int> & elements_1   = dofManager[1].getElementBlock(block_name);
     
     // test the size for equality
     TEST_EQUALITY(blk_elements.size(),elements_0.size());
     TEST_EQUALITY(blk_elements.size(),elements_1.size());

     // test each element value (a future change could break this, you might want to sort them 
     // all instead)
     for(std::size_t i=0;i<blk_elements.size();i++) {
       TEST_EQUALITY(blk_elements[i],elements_0[i]);
       TEST_EQUALITY(blk_elements[i],elements_1[i]);
     }
   }

   // verify the GID offsets are correct
   /////////////////////////////////////////////////////////

   for(std::size_t eb=0;eb<eBlocks.size();eb++) {
     std::string block_name = eBlocks[eb];
 
     const std::vector<int> & offsets_blk_T = blkDofManager.getGIDFieldOffsets(block_name,blkDofManager.getFieldNum("T"));
     const std::vector<int> & offsets_0_T   = dofManager[0].getGIDFieldOffsets(block_name,dofManager[0].getFieldNum("T"));

     const std::vector<int> & offsets_blk_x = blkDofManager.getGIDFieldOffsets(block_name,blkDofManager.getFieldNum("x"));
     const std::vector<int> & offsets_1_x   = dofManager[1].getGIDFieldOffsets(block_name,dofManager[1].getFieldNum("x"));

     int offset0 = blkDofManager.getBlockGIDOffset(block_name,0);
     int offset1 = blkDofManager.getBlockGIDOffset(block_name,1);

     TEST_EQUALITY(offsets_blk_T.size(),offsets_0_T.size());
     for(std::size_t i=0; i<offsets_blk_T.size(); i++)
        TEST_EQUALITY(offsets_blk_T[i],offsets_0_T[i]+offset0);

     TEST_EQUALITY(offsets_blk_x.size(),offsets_1_x.size());
      for(std::size_t i=0; i<offsets_blk_x.size(); i++)
        TEST_EQUALITY(offsets_blk_x[i],offsets_1_x[i]+offset1);
   }

   // verify the GIDs are correct
   /////////////////////////////////////////////////////////
   typedef BlockedDOFManager<int,int>::GlobalOrdinal GIDType;

   for(std::size_t eb=0;eb<eBlocks.size();eb++) {
     std::string block_name = eBlocks[eb];

     int offset0 = blkDofManager.getBlockGIDOffset(block_name,0);
     int offset1 = blkDofManager.getBlockGIDOffset(block_name,1);

     const std::vector<int> & elements = dofManager[0].getElementBlock(block_name);
     for(std::size_t i=0;i<elements.size();i++) {
       std::vector<GIDType> blkgids;
       std::vector<int> gids0, gids1;
       dofManager[0].getElementGIDs(elements[i],gids0);
       dofManager[1].getElementGIDs(elements[i],gids1);
       blkDofManager.getElementGIDs(elements[i],blkgids);

       for(std::size_t j=0;j<gids0.size();j++) {
         TEST_EQUALITY(0,blkgids[j+offset0].first);
         TEST_EQUALITY(gids0[j],blkgids[j+offset0].second);
       }
       for(std::size_t j=0;j<gids1.size();j++) {
         TEST_EQUALITY(1,blkgids[j+offset1].first);
         TEST_EQUALITY(gids1[j],blkgids[j+offset1].second);
       }
     }
   }
}

}
