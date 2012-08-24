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

#include "Panzer_BlockedDOFManagerFactory.hpp"
#include "Panzer_STK_Version.hpp"
#include "Panzer_STK_config.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"

typedef Intrepid::FieldContainer<double> FieldContainer;

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

template <typename IntrepidType>
RCP<const panzer::FieldPattern> buildFieldPattern()
{
   // build a geometric pattern from a single basis
   RCP<Intrepid::Basis<double,FieldContainer> > basis = rcp(new IntrepidType);
   RCP<const panzer::FieldPattern> pattern = rcp(new panzer::IntrepidFieldPattern(basis));
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

   typedef panzer::BlockedDOFManagerFactory<int,int> BDFii;

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
