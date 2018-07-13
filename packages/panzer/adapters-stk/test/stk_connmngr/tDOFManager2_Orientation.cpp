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

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_GeometricAggFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HCURL_QUAD_I1_FEM.hpp"
#include "Intrepid2_HDIV_QUAD_I1_FEM.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

typedef Kokkos::DynRankView<double,PHX::Device> FieldContainer;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::rcp_dynamic_cast;

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
   return Teuchos::rcp(new panzer_stk::STKConnManager<int>(mesh));
}

template <typename Intrepid2Type>
RCP<const panzer::Intrepid2FieldPattern> buildFieldPattern()
{
   // build a geometric pattern from a single basis
   RCP<Intrepid2::Basis<PHX::exec_space,double,double> > basis = rcp(new Intrepid2Type);
   RCP<const panzer::Intrepid2FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
   return pattern;
}


// quad tests
TEUCHOS_UNIT_TEST(tDOFManager_Orientation, buildTest_quad_edge_orientations)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      stk::ParallelMachine Comm = MPI_COMM_WORLD;
   #else
      stk::ParallelMachine Comm = WHAT_TO_DO_COMM;
   #endif

   int numProcs = stk::parallel_machine_size(Comm);
   int myRank = stk::parallel_machine_rank(Comm);

   TEUCHOS_ASSERT(numProcs==2);

   // build a geometric pattern from a single basis
   RCP<const panzer::FieldPattern> patternC1 
         = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space,double,double> >();
   RCP<const panzer::FieldPattern> patternI1 
         = buildFieldPattern<Intrepid2::Basis_HCURL_QUAD_I1_FEM<PHX::exec_space,double,double> >();
   RCP<const panzer::FieldPattern> patternDIV
         = buildFieldPattern<Intrepid2::Basis_HDIV_QUAD_I1_FEM<PHX::exec_space,double,double> >();

   RCP<panzer::ConnManager<int,int> > connManager = buildQuadMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager<int,int> > dofManager = rcp(new panzer::DOFManager<int,int>());
   // build global (or serial communicator)
   dofManager->setOrientationsRequired(true);
   TEST_EQUALITY(dofManager->getOrientationsRequired(),true);
   TEST_EQUALITY(dofManager->getConnManager(),Teuchos::null);

   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   TEST_EQUALITY(dofManager->getConnManager(),connManager);

   dofManager->addField("u",patternC1);
   dofManager->addField("b",patternI1);
   dofManager->addField("e",patternDIV);

   dofManager->buildGlobalUnknowns();
   
   const std::vector<int> & u_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("u"));
   const std::vector<int> & b_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("b"));
   const std::vector<int> & e_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("e"));

   TEST_EQUALITY(u_offsets.size(),4);
   TEST_EQUALITY(b_offsets.size(),4);
   TEST_EQUALITY(e_offsets.size(),4);

   for (std::size_t i = 0; i < 4; ++i)
      out << u_offsets[i] << " " << std::endl;
   for (std::size_t i = 0; i < 4; ++i)
      out << b_offsets[i] << " " << std::endl;
   for (std::size_t i = 0; i < 4; ++i)
      out << e_offsets[i] << " " << std::endl;

   // unfortunatly this mesh is completly uniform
   double standardO[] = { 1.0, 1.0, -1.0, -1.0 };
   if(myRank==0) {
      std::vector<int> gids;
      std::vector<double> orientation;

      // element 0
      dofManager->getElementGIDs(0,gids);
      dofManager->getElementOrientation(0,orientation);
      TEST_EQUALITY(gids.size(),12);
      TEST_EQUALITY(orientation.size(),12);

      {
         double sum = 0;
         for(std::size_t i=0;i<orientation.size();i++)
            sum += orientation[i]*orientation[i];
         TEST_EQUALITY(sum,12);
      }

      for(std::size_t i=0;i<4;i++) {
         TEST_EQUALITY(orientation[i],1.0);
      }
      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[2*i+4],standardO[i]);
      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[2*i+5],standardO[i]);

      // element 1
      dofManager->getElementGIDs(1,gids);
      dofManager->getElementOrientation(1,orientation);
      TEST_EQUALITY(gids.size(),12);
      TEST_EQUALITY(orientation.size(),12);

      {
         double sum = 0;
         for(std::size_t i=0;i<orientation.size();i++)
            sum += orientation[i]*orientation[i];
         TEST_EQUALITY(sum,12);
      }

      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i],1.0);
      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[2*i+4],standardO[i]);
      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[2*i+5],standardO[i]);
   }
   else if(myRank==1) {
      std::vector<int> gids;
      std::vector<double> orientation;

      // element 0
      dofManager->getElementOrientation(0,orientation);
      dofManager->getElementGIDs(0,gids);
      TEST_EQUALITY(gids.size(),12);
      TEST_EQUALITY(orientation.size(),12);

      {
         double sum = 0;
         for(std::size_t i=0;i<orientation.size();i++)
            sum += orientation[i]*orientation[i];
         TEST_EQUALITY(sum,12);
      }

      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i],1.0);
      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[2*i+4],standardO[i]);
      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[2*i+5],standardO[i]);

      // element 1
      dofManager->getElementGIDs(1,gids);
      dofManager->getElementOrientation(1,orientation);
      TEST_EQUALITY(gids.size(),12);
      TEST_EQUALITY(orientation.size(),12);

      {
         double sum = 0;
         for(std::size_t i=0;i<orientation.size();i++)
            sum += orientation[i]*orientation[i];
         TEST_EQUALITY(sum,12);
      }

      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i],1.0);
      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[2*i+4],standardO[i]);
      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[2*i+5],standardO[i]);
   }
}

// quad tests
TEUCHOS_UNIT_TEST(tDOFManager_Orientation, buildTest_quad_edge_orientations2)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      stk::ParallelMachine Comm = MPI_COMM_WORLD;
   #else
      stk::ParallelMachine Comm = WHAT_TO_DO_COMM;
   #endif

   int numProcs = stk::parallel_machine_size(Comm);
   int myRank = stk::parallel_machine_rank(Comm);

   TEUCHOS_ASSERT(numProcs==2);

   // build a geometric pattern from a single basis
   RCP<const panzer::FieldPattern> patternI1 
         = buildFieldPattern<Intrepid2::Basis_HCURL_QUAD_I1_FEM<PHX::exec_space,double,double> >();

   RCP<panzer::ConnManager<int,int> > connManager = buildQuadMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager<int,int> > dofManager = rcp(new panzer::DOFManager<int,int>());

   dofManager->setOrientationsRequired(true);
   TEST_EQUALITY(dofManager->getOrientationsRequired(),true);
   TEST_EQUALITY(dofManager->getConnManager(),Teuchos::null);

   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   TEST_EQUALITY(dofManager->getConnManager(),connManager);

   dofManager->addField("b",patternI1);
   dofManager->buildGlobalUnknowns();

   const std::vector<int> & b_offsets = dofManager->getGIDFieldOffsets("eblock-0_0",dofManager->getFieldNum("b"));

   TEST_EQUALITY(b_offsets.size(),4);

   // unfortunatly this mesh is completly uniform
   double standardO[] = { 1.0, 1.0, -1.0, -1.0 };
   if(myRank==0) {
      std::vector<double> orientation;

      // element 0
      dofManager->getElementOrientation(0,orientation);
      TEST_EQUALITY(orientation.size(),4);

      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i],standardO[i]);

      // element 1
      dofManager->getElementOrientation(1,orientation);
      TEST_EQUALITY(orientation.size(),4);

      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i],standardO[i]);
   }
   else if(myRank==1) {
      std::vector<double> orientation;

      // element 0
      dofManager->getElementOrientation(0,orientation);
      TEST_EQUALITY(orientation.size(),4);

      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i],standardO[i]);

      // element 1
      dofManager->getElementOrientation(1,orientation);
      TEST_EQUALITY(orientation.size(),4);

      for(std::size_t i=0;i<4;i++)
         TEST_EQUALITY(orientation[i],standardO[i]);
   }
}

// quad tests
TEUCHOS_UNIT_TEST(tDOFManager_Orientation, buildTest_quad_edge_orientations_fail)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      stk::ParallelMachine Comm = MPI_COMM_WORLD;
   #else
      stk::ParallelMachine Comm = WHAT_TO_DO_COMM;
   #endif

   int numProcs = stk::parallel_machine_size(Comm);

   TEUCHOS_ASSERT(numProcs==2);

   // build a geometric pattern from a single basis
   RCP<const panzer::FieldPattern> patternI1 
         = buildFieldPattern<Intrepid2::Basis_HCURL_QUAD_I1_FEM<PHX::exec_space,double,double> >();

   RCP<panzer::ConnManager<int,int> > connManager = buildQuadMesh(Comm,2,2,1,1);
   RCP<panzer::DOFManager<int,int> > dofManager = rcp(new panzer::DOFManager<int,int>());

   dofManager->setOrientationsRequired(true);
   TEST_EQUALITY(dofManager->getOrientationsRequired(),true);
   TEST_EQUALITY(dofManager->getConnManager(),Teuchos::null);

   dofManager->setConnManager(connManager,MPI_COMM_WORLD);
   TEST_EQUALITY(dofManager->getConnManager(),connManager);

   dofManager->addField("b",patternI1);

   //TEST_THROW(dofManager->buildGlobalUnknowns(patternI1),std::logic_error);
}

}
