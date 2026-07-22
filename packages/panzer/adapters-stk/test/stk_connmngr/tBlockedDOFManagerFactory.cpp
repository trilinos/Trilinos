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

#include "Panzer_BlockedDOFManagerFactory.hpp"
#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

namespace panzer_stk {

Teuchos::RCP<panzer_stk::STK_Interface> buildQuadMesh(stk::ParallelMachine comm,int xelmts,int yelmts,
                                                                                    int xblocks,int yblocks)
{
   Teuchos::ParameterList pl;
   pl.set<int>("X Elements",xelmts);
   pl.set<int>("Y Elements",yelmts);
   pl.set<int>("X Blocks",xblocks);
   pl.set<int>("Y Blocks",yblocks);

   panzer_stk::SquareQuadMeshFactory meshFact;
   meshFact.setParameterList(Teuchos::rcpFromRef(pl));
   
   return meshFact.buildMesh(comm);
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
TEUCHOS_UNIT_TEST(tBlockedDOFManagerFactory, basic_test)
{
//    // build global (or serial communicator)
//    #ifdef HAVE_MPI
//       stk::ParallelMachine Comm = MPI_COMM_WORLD;
//    #else
//       stk::ParallelMachine Comm = WHAT_TO_DO_COMM;
//    #endif
// 
   // int numProcs = stk::parallel_machine_size(Comm);
   // int myRank = stk::parallel_machine_rank(Comm);

   typedef panzer::BlockedDOFManagerFactory BDFii;

   bool result = false;
   result = BDFii::requiresBlocking("");                     TEST_ASSERT(!result);
   result = BDFii::requiresBlocking("UX UY P");              TEST_ASSERT(!result);
   result = BDFii::requiresBlocking("blocked: UX - UY - P"); TEST_ASSERT(result);
   result = BDFii::requiresBlocking("blocked: UX UY - P");   TEST_ASSERT(result);
   result = BDFii::requiresBlocking("blocked: UX - UY P");   TEST_ASSERT(result);
   result = BDFii::requiresBlocking("blocked: UX UY P");     TEST_ASSERT(result);
   TEST_THROW(BDFii::requiresBlocking("blocked: - UX"),std::logic_error);
   TEST_THROW(BDFii::requiresBlocking("blocked: UX - - P"),std::logic_error);
   
   
   {
      std::vector<std::vector<std::string> > blocks;
      BDFii::buildBlocking("blocked: UX - UY - P",blocks); 

      TEST_EQUALITY(blocks.size(),1);
      TEST_EQUALITY(blocks[0].size(),3);
      TEST_EQUALITY(blocks[0][0],"UX");
      TEST_EQUALITY(blocks[0][1],"UY");
      TEST_EQUALITY(blocks[0][2],"P");
   }

   {
      std::vector<std::vector<std::string> > blocks;
      BDFii::buildBlocking("blocked: UX UY - P",blocks); 

      TEST_EQUALITY(blocks.size(),2);
      TEST_EQUALITY(blocks[0].size(),1);
      TEST_EQUALITY(blocks[0][0],"UX");

      TEST_EQUALITY(blocks[1].size(),2);
      TEST_EQUALITY(blocks[1][0],"UY");
      TEST_EQUALITY(blocks[1][1],"P");
   }

   {
      std::vector<std::vector<std::string> > blocks;
      BDFii::buildBlocking("blocked: UX - UY P",blocks); 

      TEST_EQUALITY(blocks.size(),2);
      TEST_EQUALITY(blocks[0].size(),2);
      TEST_EQUALITY(blocks[0][0],"UX");
      TEST_EQUALITY(blocks[0][1],"UY");

      TEST_EQUALITY(blocks[1].size(),1);
      TEST_EQUALITY(blocks[1][0],"P");
   }

   {
      std::vector<std::vector<std::string> > blocks;
      BDFii::buildBlocking("blocked: UX P UY",blocks); 

      TEST_EQUALITY(blocks.size(),3);

      TEST_EQUALITY(blocks[0].size(),1);
      TEST_EQUALITY(blocks[0][0],"UX");

      TEST_EQUALITY(blocks[1].size(),1);
      TEST_EQUALITY(blocks[1][0],"P");

      TEST_EQUALITY(blocks[2].size(),1);
      TEST_EQUALITY(blocks[2][0],"UY");
   }
}

}
