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

#include "Panzer_EpetraLinearObjFactory.hpp"
#include "Panzer_Traits.hpp"

// for testing gather/scatter construction
#include "Panzer_PureBasis.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_EpetraLinearObjContainer.hpp"
#include "Panzer_BlockedEpetraLinearObjContainer.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_PauseToAttach.hpp"

#include "UnitTest_UniqueGlobalIndexer.hpp"

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

Teuchos::RCP<Epetra_CrsMatrix> getSubBlock(int i,int j,Thyra::LinearOpBase<double> & lo)
{
   Thyra::BlockedLinearOpBase<double> & blo = Teuchos::dyn_cast<Thyra::BlockedLinearOpBase<double> >(lo);
   Teuchos::RCP<Epetra_Operator> e_blo = Thyra::get_Epetra_Operator(*blo.getNonconstBlock(i,j));

   return rcp_dynamic_cast<Epetra_CrsMatrix>(e_blo);
}

Teuchos::RCP<const Epetra_CrsMatrix> getSubBlock(int i,int j,const Thyra::LinearOpBase<double> & lo)
{
   const Thyra::BlockedLinearOpBase<double> & blo = Teuchos::dyn_cast<const Thyra::BlockedLinearOpBase<double> >(lo);
   Teuchos::RCP<const Epetra_Operator> e_blo = Thyra::get_Epetra_Operator(*blo.getBlock(i,j));

   return rcp_dynamic_cast<const Epetra_CrsMatrix>(e_blo);
}

TEUCHOS_UNIT_TEST(tBlockedLinearObjFactory, intializeContainer_epetra)
{
   panzer::BlockedEpetraLinearObjContainer container;

   TEST_ASSERT(container.checkCompatibility());
}

TEUCHOS_UNIT_TEST(tBlockedEpetraLinearObjFactory, epetra_factory_tests)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   // pauseToAttach();

   typedef LinearObjContainer LOC;
   typedef BlockedEpetraLinearObjContainer BLOC;

   int numBlocks = 3;
   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<panzer::UniqueGlobalIndexer<int,int> > indexer 
         = rcp(new panzer::unit_test::UniqueGlobalIndexer<int>(myRank,numProc));
   RCP<const panzer::UniqueGlobalIndexer<int,std::pair<int,int> > > blkIndexer 
         = rcp(new panzer::unit_test::BlockUniqueGlobalIndexer<int>(numBlocks,myRank,numProc));

   std::vector<int> ownedIndices, ownedAndSharedIndices;
   indexer->getOwnedIndices(ownedIndices);
   indexer->getOwnedAndSharedIndices(ownedAndSharedIndices);

   std::vector<RCP<const panzer::UniqueGlobalIndexer<int,int> > > indexers;
   for(int i=0;i<numBlocks;i++)
      indexers.push_back(indexer); // 3x3 square blocks

   BlockedEpetraLinearObjFactory<panzer::Traits,int> factory(eComm,blkIndexer,indexers);

   RCP<LinearObjContainer> container = factory.buildLinearObjContainer();
   RCP<LinearObjContainer> ghosted = factory.buildGhostedLinearObjContainer();
   TEST_ASSERT(container!=Teuchos::null);
   TEST_ASSERT(ghosted!=Teuchos::null);


   RCP<BLOC> bContainer = rcp_dynamic_cast<BLOC>(container);
   RCP<BLOC> b_ghosted = rcp_dynamic_cast<BLOC>(ghosted);
   TEST_ASSERT(bContainer!=Teuchos::null);
   TEST_ASSERT(b_ghosted!=Teuchos::null);

   // tests global initialize
   {
      // Generic code
      /////////////////////////////////////////////////////////////
   
      // test individial initializers
      factory.initializeContainer(LOC::X,*container);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::DxDt,*container);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::F,*container);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::Mat,*container);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   
      // jacobian and residual vector output
      factory.initializeContainer(LOC::F | LOC::Mat,*container);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   
      // x and time dertivative input
      factory.initializeContainer(LOC::X | LOC::DxDt,*container);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      // everything
      factory.initializeContainer(LOC::X | LOC::DxDt | LOC::F | LOC::Mat,*container);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   
      // Epetra specific code
      /////////////////////////////////////////////////////////////
   
      // test individial initializers
      factory.initializeContainer(LOC::X,*bContainer);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::DxDt,*bContainer);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::F,*bContainer);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      factory.initializeContainer(LOC::Mat,*bContainer);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   
      // jacobian and residual vector output
      factory.initializeContainer(LOC::F | LOC::Mat,*bContainer);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   
      // x and time dertivative input
      factory.initializeContainer(LOC::X | LOC::DxDt,*bContainer);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_A(),    Teuchos::null)
   
      // everything
      factory.initializeContainer(LOC::X | LOC::DxDt | LOC::F | LOC::Mat,*bContainer);
      TEST_ASSERT(bContainer->get_x()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_dxdt()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_f()!=Teuchos::null);
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);
   }
}

TEUCHOS_UNIT_TEST(tBlockedEpetraLinearObjFactory, ghostToGlobal)
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
   using Thyra::ProductVectorBase;
   using Thyra::SpmdVectorBase;

   // pauseToAttach();

   int numBlocks = 2;
   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();
 
   typedef BlockedEpetraLinearObjContainer BLOC;

   RCP<panzer::UniqueGlobalIndexer<int,int> > indexer 
         = rcp(new panzer::unit_test::UniqueGlobalIndexer<int>(myRank,numProc));
   RCP<const panzer::UniqueGlobalIndexer<int,std::pair<int,int> > > blkIndexer 
         = rcp(new panzer::unit_test::BlockUniqueGlobalIndexer<int>(numBlocks,myRank,numProc));

   std::vector<RCP<const panzer::UniqueGlobalIndexer<int,int> > > indexers;
   for(int i=0;i<numBlocks;i++)
      indexers.push_back(indexer); // 2x2 square blocks

   Teuchos::RCP<BlockedEpetraLinearObjFactory<panzer::Traits,int> > la_factory
         = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(eComm,blkIndexer,indexers));

   Teuchos::RCP<LinearObjContainer> global  = la_factory->buildLinearObjContainer();
   Teuchos::RCP<LinearObjContainer> ghosted = la_factory->buildGhostedLinearObjContainer();

   la_factory->initializeContainer(LinearObjContainer::Mat,*global);
   la_factory->initializeGhostedContainer(LinearObjContainer::Mat,*ghosted);

   Teuchos::rcp_dynamic_cast<BlockedEpetraLinearObjContainer>(ghosted)->initializeMatrix(1.0);

   la_factory->ghostToGlobalContainer(*ghosted,*global,LinearObjContainer::Mat);

   Teuchos::RCP<Thyra::LinearOpBase<double> > th_A = Teuchos::rcp_dynamic_cast<BlockedEpetraLinearObjContainer>(global)->get_A();
   Teuchos::RCP<Thyra::BlockedLinearOpBase<double> > blk_A = Teuchos::rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(th_A);
   
  
   Teuchos::RCP<Epetra_Operator> cA_00 = Thyra::get_Epetra_Operator(*blk_A->getNonconstBlock(0,0));
   Teuchos::RCP<Epetra_Operator> cA_01 = Thyra::get_Epetra_Operator(*blk_A->getNonconstBlock(0,1));
   Teuchos::RCP<Epetra_Operator> cA_10 = Thyra::get_Epetra_Operator(*blk_A->getNonconstBlock(1,0));
   Teuchos::RCP<Epetra_Operator> cA_11 = Thyra::get_Epetra_Operator(*blk_A->getNonconstBlock(1,1));

   Teuchos::RCP<Epetra_CrsMatrix> A_00 = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(cA_00);
   Teuchos::RCP<Epetra_CrsMatrix> A_01 = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(cA_01);

   Teuchos::RCP<Epetra_CrsMatrix> A_10 = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(cA_10);
   Teuchos::RCP<Epetra_CrsMatrix> A_11 = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(cA_11);

   std::cout << "A_00 = \n";
   A_00->Print(out);

   std::cout << "A_01 = \n";
   A_01->Print(out);

   std::cout << "A_10 = \n";
   A_10->Print(out);

   std::cout << "A_11 = \n";
   A_11->Print(out);
}

TEUCHOS_UNIT_TEST(tBlockedEpetraLinearObjFactory, graph_constr)
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
   using Thyra::ProductVectorBase;
   using Thyra::SpmdVectorBase;

   // pauseToAttach();

   int numBlocks = 2;
   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();
 
   typedef BlockedEpetraLinearObjContainer BLOC;

   RCP<panzer::UniqueGlobalIndexer<int,int> > indexer 
         = rcp(new panzer::unit_test::UniqueGlobalIndexer<int>(myRank,numProc));
   RCP<const panzer::UniqueGlobalIndexer<int,std::pair<int,int> > > blkIndexer 
         = rcp(new panzer::unit_test::BlockUniqueGlobalIndexer<int>(numBlocks,myRank,numProc));

   std::vector<RCP<const panzer::UniqueGlobalIndexer<int,int> > > indexers;
   for(int i=0;i<numBlocks;i++)
      indexers.push_back(indexer); // 2x2 square blocks

   Teuchos::RCP<BlockedEpetraLinearObjFactory<panzer::Traits,int> > la_factory
         = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(eComm,blkIndexer,indexers));

   Teuchos::RCP<Epetra_CrsMatrix> A_00 = la_factory->getGhostedEpetraMatrix(0,0); A_00->PutScalar(1.0);
   Teuchos::RCP<Epetra_CrsMatrix> A_01 = la_factory->getGhostedEpetraMatrix(0,1); A_01->PutScalar(1.0);

   Teuchos::RCP<Epetra_CrsMatrix> A_10 = la_factory->getGhostedEpetraMatrix(1,0); A_10->PutScalar(1.0);
   Teuchos::RCP<Epetra_CrsMatrix> A_11 = la_factory->getGhostedEpetraMatrix(1,1); A_11->PutScalar(1.0);

   std::cout << "A_00 = \n";
   A_00->Print(out);

   std::cout << "A_01 = \n";
   A_01->Print(out);

   std::cout << "A_10 = \n";
   A_10->Print(out);

   std::cout << "A_11 = \n";
   A_11->Print(out);
}

TEUCHOS_UNIT_TEST(tBlockedEpetraLinearObjFactory, adjustDirichlet)
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
   using Thyra::ProductVectorBase;
   using Thyra::SpmdVectorBase;

   // pauseToAttach();

   int numBlocks = 3;
   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();
 
   typedef BlockedEpetraLinearObjContainer BLOC;

   RCP<panzer::UniqueGlobalIndexer<int,int> > indexer 
         = rcp(new panzer::unit_test::UniqueGlobalIndexer<int>(myRank,numProc));
   RCP<const panzer::UniqueGlobalIndexer<int,std::pair<int,int> > > blkIndexer 
         = rcp(new panzer::unit_test::BlockUniqueGlobalIndexer<int>(numBlocks,myRank,numProc));

   std::vector<int> ownedIndices, ownedAndSharedIndices;
   indexer->getOwnedIndices(ownedIndices);
   indexer->getOwnedAndSharedIndices(ownedAndSharedIndices);

   std::vector<RCP<const panzer::UniqueGlobalIndexer<int,int> > > indexers;
   for(int i=0;i<numBlocks;i++)
      indexers.push_back(indexer); // 3x3 square blocks

   Teuchos::RCP<BlockedEpetraLinearObjFactory<panzer::Traits,int> > la_factory
         = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(eComm,blkIndexer,indexers));

   RCP<LinearObjContainer> ghosted_0   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_1   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_sys = la_factory->buildGhostedLinearObjContainer();

   la_factory->initializeGhostedContainer(LinearObjContainer::X,*ghosted_0);
   la_factory->initializeGhostedContainer(LinearObjContainer::X,*ghosted_1);
   la_factory->initializeGhostedContainer(LinearObjContainer::F | LinearObjContainer::Mat,*ghosted_sys);

   RCP<BLOC> b_0   = rcp_dynamic_cast<BLOC>(ghosted_0);
   RCP<BLOC> b_1   = rcp_dynamic_cast<BLOC>(ghosted_1);
   RCP<BLOC> b_sys = rcp_dynamic_cast<BLOC>(ghosted_sys);

   TEST_ASSERT(!Teuchos::is_null(b_0->get_x()));
   TEST_ASSERT(!Teuchos::is_null(b_1->get_x()));
   TEST_ASSERT(!Teuchos::is_null(b_sys->get_f()));
   TEST_ASSERT(!Teuchos::is_null(b_sys->get_A()));

   Thyra::assign(b_0->get_x().ptr(),0.0); // put some garbage in the systems
   Thyra::assign(b_1->get_x().ptr(),0.0); // put some garbage in the systems
   Thyra::assign(b_sys->get_f().ptr(),-3.0); // put some garbage in the systems

   // b_sys->get_A()->PutScalar(-3.0);
   for(int i=0;i<numBlocks;i++)
      for(int j=0;j<numBlocks;j++)
         getSubBlock(i,j,*b_sys->get_A())->PutScalar(-3.0);

   // there are 3 cases for adjustDirichlet
   //   1. Local set only for GID
   //   2. Set on multiple processors
   //   3. Set remotely

   if(myRank==0) {   
      for(int i=0;i<numBlocks;i++) {
         RCP<Thyra::VectorBase<double> > x_0 = rcp_dynamic_cast<ProductVectorBase<double> >(b_0->get_x())->getNonconstVectorBlock(i);
         RCP<Thyra::VectorBase<double> > x_1 = rcp_dynamic_cast<ProductVectorBase<double> >(b_1->get_x())->getNonconstVectorBlock(i);

         Teuchos::ArrayRCP<double> data_0,data_1;
         rcp_dynamic_cast<SpmdVectorBase<double> >(x_0)->getNonconstLocalData(Teuchos::ptrFromRef(data_0)); 
         rcp_dynamic_cast<SpmdVectorBase<double> >(x_1)->getNonconstLocalData(Teuchos::ptrFromRef(data_1)); 

         // case 0
         data_0[0] = 1.0; // GID = 0
         data_1[0] = 1.0;
   
         // case 1
         data_0[2] = 1.0; // GID = 2
         data_1[2] = 2.0;
   
         // case 2
         data_1[5] = 2.0; // GID = 5
      }
   }
   else if(myRank==1) {
      for(int i=0;i<numBlocks;i++) {
         RCP<Thyra::VectorBase<double> > x_0 = rcp_dynamic_cast<ProductVectorBase<double> >(b_0->get_x())->getNonconstVectorBlock(i);
         RCP<Thyra::VectorBase<double> > x_1 = rcp_dynamic_cast<ProductVectorBase<double> >(b_1->get_x())->getNonconstVectorBlock(i);

         Teuchos::ArrayRCP<double> data_0,data_1;
         rcp_dynamic_cast<SpmdVectorBase<double> >(x_0)->getNonconstLocalData(Teuchos::ptrFromRef(data_0)); 
         rcp_dynamic_cast<SpmdVectorBase<double> >(x_1)->getNonconstLocalData(Teuchos::ptrFromRef(data_1)); 

         // case 0
         data_0[3] = 1.0; // GID = 9
         data_1[3] = 1.0;
   
         // case 1
         data_0[0] = 1.0; // GID =2
         data_1[0] = 2.0;
   
         // case 2
         data_1[6] = 2.0; // GID = 4
      }
   }
   else 
      TEUCHOS_ASSERT(false);

   out << "LOCAL " << std::endl;
   b_0->get_x()->describe(out,Teuchos::VERB_HIGH);
   out << std::endl;
   out << "GLOBAL " << std::endl;
   b_1->get_x()->describe(out,Teuchos::VERB_HIGH);
   out << std::endl;

   // run test for conditions
   la_factory->adjustForDirichletConditions(*ghosted_0,*ghosted_1,*ghosted_sys);

   int numEntries = 0;
   double * values = 0;
   int * indices = 0;

   if(myRank==0) {   
      RCP<const Thyra::LinearOpBase<double> > A = b_sys->get_A();

      for(int i=0;i<numBlocks;i++) {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[0],-3.0);     // case 0
         TEST_EQUALITY(data[2],-3.0/2.0); // case 1
         TEST_EQUALITY(data[5],0.0);      // case 2

         for(int j=0;j<numBlocks;j++) {
            RCP<const Epetra_CrsMatrix> subA = getSubBlock(i,j,*A);

            subA->ExtractMyRowView(0,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0);

            subA->ExtractMyRowView(2,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0/2.0);

            subA->ExtractMyRowView(5,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],0.0);
         }
      }
   }
   else if(myRank==1) {
      RCP<const Thyra::LinearOpBase<double> > A = b_sys->get_A();

      for(int i=0;i<numBlocks;i++) {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[3],-3.0);     // case 0
         TEST_EQUALITY(data[0],-3.0/2.0); // case 1
         TEST_EQUALITY(data[6],0.0);     // case 2

         for(int j=0;j<numBlocks;j++) {
            RCP<const Epetra_CrsMatrix> subA = getSubBlock(i,j,*A);

            subA->ExtractMyRowView(3,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0);
   
            subA->ExtractMyRowView(0,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0/2.0);
   
            subA->ExtractMyRowView(6,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],0.0);
         }
      }
   }
   else 
      TEUCHOS_ASSERT(false);
}

TEUCHOS_UNIT_TEST(tBlockedEpetraLinearObjFactory, node_cell)
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
   using Thyra::ProductVectorBase;
   using Thyra::SpmdVectorBase;

   // pauseToAttach();

   int numBlocks = 2;
   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();
 
   typedef BlockedEpetraLinearObjContainer BLOC;

   RCP<panzer::UniqueGlobalIndexer<int,int> > indexer_node
         = rcp(new panzer::unit_test::UniqueGlobalIndexer<int>(myRank,numProc));
   RCP<panzer::UniqueGlobalIndexer<int,int> > indexer_cell
         = rcp(new panzer::unit_test::UniqueGlobalIndexer_Element<int>(myRank,numProc));
   RCP<const panzer::UniqueGlobalIndexer<int,std::pair<int,int> > > blkIndexer 
         = rcp(new panzer::unit_test::BlockUniqueGlobalIndexer<int>(numBlocks,myRank,numProc));

   std::vector<RCP<const panzer::UniqueGlobalIndexer<int,int> > > indexers;
   indexers.push_back(indexer_node);
   indexers.push_back(indexer_cell);

   Teuchos::RCP<BlockedEpetraLinearObjFactory<panzer::Traits,int> > la_factory
         = Teuchos::rcp(new panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int>(eComm,blkIndexer,indexers));

   RCP<LinearObjContainer> ghosted_0   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_1   = la_factory->buildGhostedLinearObjContainer();
   RCP<LinearObjContainer> ghosted_sys = la_factory->buildGhostedLinearObjContainer();

   la_factory->initializeGhostedContainer(LinearObjContainer::X,*ghosted_0);
   la_factory->initializeGhostedContainer(LinearObjContainer::X,*ghosted_1);
   la_factory->initializeGhostedContainer(LinearObjContainer::F | LinearObjContainer::Mat,*ghosted_sys);

   RCP<BLOC> b_0   = rcp_dynamic_cast<BLOC>(ghosted_0);
   RCP<BLOC> b_1   = rcp_dynamic_cast<BLOC>(ghosted_1);
   RCP<BLOC> b_sys = rcp_dynamic_cast<BLOC>(ghosted_sys);

   TEST_ASSERT(!Teuchos::is_null(b_0->get_x()));
   TEST_ASSERT(!Teuchos::is_null(b_1->get_x()));
   TEST_ASSERT(!Teuchos::is_null(b_sys->get_f()));
   TEST_ASSERT(!Teuchos::is_null(b_sys->get_A()));

   Thyra::assign(b_0->get_x().ptr(),0.0); // put some garbage in the systems
   Thyra::assign(b_1->get_x().ptr(),0.0); // put some garbage in the systems
   Thyra::assign(b_sys->get_f().ptr(),-3.0); // put some garbage in the systems

   // b_sys->get_A()->PutScalar(-3.0);
   for(int i=0;i<numBlocks;i++)
      for(int j=0;j<numBlocks;j++)
         getSubBlock(i,j,*b_sys->get_A())->PutScalar(-3.0);

   // there are 3 cases for adjustDirichlet
   //   1. Local set only for GID
   //   2. Set on multiple processors
   //   3. Set remotely

   if(myRank==0) {   
      RCP<Thyra::VectorBase<double> > x_0 = rcp_dynamic_cast<ProductVectorBase<double> >(b_0->get_x())->getNonconstVectorBlock(0);
      RCP<Thyra::VectorBase<double> > x_1 = rcp_dynamic_cast<ProductVectorBase<double> >(b_1->get_x())->getNonconstVectorBlock(0);

      Teuchos::ArrayRCP<double> data_0,data_1;
      rcp_dynamic_cast<SpmdVectorBase<double> >(x_0)->getNonconstLocalData(Teuchos::ptrFromRef(data_0)); 
      rcp_dynamic_cast<SpmdVectorBase<double> >(x_1)->getNonconstLocalData(Teuchos::ptrFromRef(data_1)); 

      // case 0
      data_0[0] = 1.0; // GID = 0
      data_1[0] = 1.0;

      // case 1
      data_0[2] = 1.0; // GID = 2
      data_1[2] = 2.0;

      // case 2
      data_1[5] = 2.0; // GID = 5

      {
         x_0 = rcp_dynamic_cast<ProductVectorBase<double> >(b_0->get_x())->getNonconstVectorBlock(1);
         x_1 = rcp_dynamic_cast<ProductVectorBase<double> >(b_1->get_x())->getNonconstVectorBlock(1);

         rcp_dynamic_cast<SpmdVectorBase<double> >(x_0)->getNonconstLocalData(Teuchos::ptrFromRef(data_0)); 
         rcp_dynamic_cast<SpmdVectorBase<double> >(x_1)->getNonconstLocalData(Teuchos::ptrFromRef(data_1)); 

         data_1[0] = 2.0;
      }
   }
   else if(myRank==1) {
      RCP<Thyra::VectorBase<double> > x_0 = rcp_dynamic_cast<ProductVectorBase<double> >(b_0->get_x())->getNonconstVectorBlock(0);
      RCP<Thyra::VectorBase<double> > x_1 = rcp_dynamic_cast<ProductVectorBase<double> >(b_1->get_x())->getNonconstVectorBlock(0);

      Teuchos::ArrayRCP<double> data_0,data_1;
      rcp_dynamic_cast<SpmdVectorBase<double> >(x_0)->getNonconstLocalData(Teuchos::ptrFromRef(data_0)); 
      rcp_dynamic_cast<SpmdVectorBase<double> >(x_1)->getNonconstLocalData(Teuchos::ptrFromRef(data_1)); 

      // case 0
      data_0[3] = 1.0; // GID = 9
      data_1[3] = 1.0;

      // case 1
      data_0[0] = 1.0; // GID =2
      data_1[0] = 2.0;

      // case 2
      data_1[6] = 2.0; // GID = 4
   }
   else 
      TEUCHOS_ASSERT(false);

   out << "LOCAL " << std::endl;
   b_0->get_x()->describe(out,Teuchos::VERB_HIGH);
   out << std::endl;
   out << "GLOBAL " << std::endl;
   b_1->get_x()->describe(out,Teuchos::VERB_HIGH);
   out << std::endl;

   // run test for conditions
   la_factory->adjustForDirichletConditions(*ghosted_0,*ghosted_1,*ghosted_sys);

   int numEntries = 0;
   double * values = 0;
   int * indices = 0;

   if(myRank==0) {   
      RCP<const Thyra::LinearOpBase<double> > A = b_sys->get_A();

      int i = 0;
      {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[0],-3.0);     // case 0
         TEST_EQUALITY(data[2],-3.0/2.0); // case 1
         TEST_EQUALITY(data[5],0.0);     // case 2

         for(int j=0;j<numBlocks;j++) {
            RCP<const Epetra_CrsMatrix> subA = getSubBlock(i,j,*A);

            subA->ExtractMyRowView(0,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0);

            subA->ExtractMyRowView(2,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0/2.0);

            subA->ExtractMyRowView(5,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],0.0);
         }
      }

      i = 1;
      {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[0],0.0);

         for(int j=0;j<numBlocks;j++) {
            RCP<const Epetra_CrsMatrix> subA = getSubBlock(i,j,*A);

            subA->ExtractMyRowView(0,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],0.0);
         }
      }
   }
   else if(myRank==1) {
      RCP<const Thyra::LinearOpBase<double> > A = b_sys->get_A();

      int i = 0;
      {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[3],-3.0);     // case 0
         TEST_EQUALITY(data[0],-3.0/2.0); // case 1
         TEST_EQUALITY(data[6],0.0);     // case 2

         for(int j=0;j<numBlocks;j++) {
            RCP<const Epetra_CrsMatrix> subA = getSubBlock(i,j,*A);

            subA->ExtractMyRowView(3,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0);
   
            subA->ExtractMyRowView(0,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0/2.0);
   
            subA->ExtractMyRowView(6,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],0.0);
         }
      }

      i = 1;
      {
         Teuchos::ArrayRCP<const double> data;
         RCP<const Thyra::VectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_sys->get_f())->getVectorBlock(i);
         rcp_dynamic_cast<const SpmdVectorBase<double> >(f)->getLocalData(Teuchos::ptrFromRef(data)); 
   
         TEST_EQUALITY(data[0],-3.0);

         for(int j=0;j<numBlocks;j++) {
            RCP<const Epetra_CrsMatrix> subA = getSubBlock(i,j,*A);

            subA->ExtractMyRowView(0,numEntries,values,indices);
            for(int k=0;k<numEntries;k++) TEST_EQUALITY(values[k],-3.0);
         }
      }
   }
   else 
      TEUCHOS_ASSERT(false);
}

TEUCHOS_UNIT_TEST(tBlockedEpetraLinearObjFactory, exclusion)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else
      Teuchos::RCP<Epetra_Comm> eComm = Teuchos::rcp(new Epetra_SerialComm());
   #endif

   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   // pauseToAttach();

   typedef LinearObjContainer LOC;
   typedef BlockedEpetraLinearObjContainer BLOC;

   int numBlocks = 3;
   int myRank = eComm->MyPID();
   int numProc = eComm->NumProc();

   RCP<panzer::UniqueGlobalIndexer<int,int> > indexer 
         = rcp(new panzer::unit_test::UniqueGlobalIndexer<int>(myRank,numProc));
   RCP<const panzer::UniqueGlobalIndexer<int,std::pair<int,int> > > blkIndexer 
         = rcp(new panzer::unit_test::BlockUniqueGlobalIndexer<int>(numBlocks,myRank,numProc));

   std::vector<int> ownedIndices, ownedAndSharedIndices;
   indexer->getOwnedIndices(ownedIndices);
   indexer->getOwnedAndSharedIndices(ownedAndSharedIndices);

   std::vector<RCP<const panzer::UniqueGlobalIndexer<int,int> > > indexers;
   for(int i=0;i<numBlocks;i++)
      indexers.push_back(indexer); // 3x3 square blocks

   BlockedEpetraLinearObjFactory<panzer::Traits,int> factory(eComm,blkIndexer,indexers);
 
   // exclude some pairs
   std::vector<std::pair<int,int> > exPairs;
   exPairs.push_back(std::make_pair(0,2));
   exPairs.push_back(std::make_pair(2,1));
   factory.addExcludedPairs(exPairs);
   factory.addExcludedPair(1,2);

   RCP<LinearObjContainer> container = factory.buildLinearObjContainer();
   RCP<LinearObjContainer> ghosted = factory.buildGhostedLinearObjContainer();
   TEST_ASSERT(container!=Teuchos::null);
   TEST_ASSERT(ghosted!=Teuchos::null);

   RCP<BLOC> bContainer = rcp_dynamic_cast<BLOC>(container);
   RCP<BLOC> b_ghosted = rcp_dynamic_cast<BLOC>(ghosted);
   TEST_ASSERT(bContainer!=Teuchos::null);
   TEST_ASSERT(b_ghosted!=Teuchos::null);

   // tests global initialize
   {
      // Generic code
      /////////////////////////////////////////////////////////////
      factory.initializeContainer(LOC::Mat,*container);
      TEST_EQUALITY(bContainer->get_x(),    Teuchos::null)
      TEST_EQUALITY(bContainer->get_dxdt(), Teuchos::null)
      TEST_EQUALITY(bContainer->get_f(),    Teuchos::null)
      TEST_ASSERT(bContainer->get_A()!=Teuchos::null);

      RCP<Thyra::BlockedLinearOpBase<double> > blo 
         = rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(bContainer->get_A());

      TEST_ASSERT(!blo->getNonconstBlock(0,0).is_null());
      TEST_ASSERT(!blo->getNonconstBlock(0,1).is_null());
      TEST_ASSERT(blo->getNonconstBlock(0,2).is_null());

      TEST_ASSERT(!blo->getNonconstBlock(1,0).is_null());
      TEST_ASSERT(!blo->getNonconstBlock(1,1).is_null());
      TEST_ASSERT(blo->getNonconstBlock(1,2).is_null());

      TEST_ASSERT(!blo->getNonconstBlock(2,0).is_null());
      TEST_ASSERT(blo->getNonconstBlock(2,1).is_null());
      TEST_ASSERT(!blo->getNonconstBlock(2,2).is_null());
   }
}

}
