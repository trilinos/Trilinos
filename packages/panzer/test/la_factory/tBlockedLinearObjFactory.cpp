// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
#include "Panzer_BlockedLinearObjContainer.hpp"
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#include "Panzer_PauseToAttach.hpp"

#include "UnitTest_UniqueGlobalIndexer.hpp"

#include "Thyra_EpetraThyraWrappers.hpp"

#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {

TEUCHOS_UNIT_TEST(tBlockedLinearObjFactory, intializeContainer_epetra)
{
   panzer::BlockedLinearObjContainer<EpetraLinearObjContainer> container;

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
   typedef BlockedLinearObjContainer<EpetraLinearObjContainer> BLOC;

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

}
