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

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Panzer_Traits.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"

typedef Kokkos::DynRankView<double, PHX::Device> FieldContainer;
typedef panzer::BlockedTpetraLinearObjFactory<panzer::Traits, double, panzer::LocalOrdinal, panzer::GlobalOrdinal> BTLOF;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

namespace panzer_stk
{

   Teuchos::RCP<panzer::ConnManager> buildQuadMesh(stk::ParallelMachine comm, int xelmts, int yelmts, int xblocks, int yblocks)
   {
      Teuchos::ParameterList pl;
      pl.set<int>("X Elements", xelmts);
      pl.set<int>("Y Elements", yelmts);
      pl.set<int>("X Blocks", xblocks);
      pl.set<int>("Y Blocks", yblocks);

      panzer_stk::SquareQuadMeshFactory meshFact;
      meshFact.setParameterList(Teuchos::rcpFromRef(pl));

      Teuchos::RCP<panzer_stk::STK_Interface> mesh = meshFact.buildMesh(comm);
      return Teuchos::rcp(new panzer_stk::STKConnManager(mesh));
   }

   template <typename Intrepid2Type>
   RCP<const panzer::FieldPattern> buildFieldPattern()
   {
      // build a geometric pattern from a single basis
      RCP<Intrepid2::Basis<PHX::exec_space, double, double>> basis = rcp(new Intrepid2Type);
      RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
      return pattern;
   }

   // quad tests
   TEUCHOS_UNIT_TEST(tTpetraLinearObjFactory, buildTest_quad)
   {

// build global (or serial communicator)
#ifdef HAVE_MPI
      stk::ParallelMachine Comm = MPI_COMM_WORLD;
#else
      stk::ParallelMachine Comm = WHAT_TO_DO_COMM;
#endif

      Teuchos::RCP<const Teuchos::MpiComm<int>> tComm = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

      int numProcs = stk::parallel_machine_size(Comm);
      int myRank = stk::parallel_machine_rank(Comm);

      TEUCHOS_ASSERT(numProcs <= 2);

      // build a geometric pattern from a single basis
      RCP<const panzer::FieldPattern> patternC1 = buildFieldPattern<Intrepid2::Basis_HGRAD_QUAD_C1_FEM<PHX::exec_space, double, double>>();

      RCP<panzer::ConnManager> connManager = buildQuadMesh(Comm, 2, 2, 1, 1);
      RCP<panzer::DOFManager> dofManager = rcp(new panzer::DOFManager());
      dofManager->setConnManager(connManager, MPI_COMM_WORLD);
      dofManager->addField("u", patternC1);
      dofManager->buildGlobalUnknowns();
      std::vector<Teuchos::RCP<const panzer::GlobalIndexer>> providers{dofManager};

      BTLOF laFactory(tComm.getConst(), providers);
      auto map = laFactory.getMap(0);
      auto gMap = laFactory.getGhostedMap(0);
      auto graph = laFactory.getGraph(0, 0);
      auto gGraph = laFactory.getGhostedGraph(0, 0);

      std::vector<panzer::GlobalOrdinal> owned, ownedAndGhosted;
      dofManager->getOwnedIndices(owned);
      dofManager->getOwnedAndGhostedIndices(ownedAndGhosted);

      // test maps
      {
         TEST_EQUALITY(map->getLocalNumElements(), (int)owned.size());
         TEST_EQUALITY(gMap->getLocalNumElements(), (int)ownedAndGhosted.size());

         // test indices
         TEST_EQUALITY(owned[0], 0);
         TEST_EQUALITY(ownedAndGhosted[0], 0);
         for (std::size_t i = 1; i < owned.size(); i++)
            TEST_ASSERT(map->getGlobalElement(owned[i]));
         for (std::size_t i = 1; i < ownedAndGhosted.size(); i++)
            TEST_ASSERT(gMap->getGlobalElement(ownedAndGhosted[i]));
      }

      // test ograph
      {
         TEST_ASSERT(gGraph->isFillComplete());

         TEST_EQUALITY(gGraph->getLocalNumRows(), (int)ownedAndGhosted.size());
         TEST_EQUALITY(gGraph->getLocalMaxNumRowEntries(), numProcs == 2 ? 6 : 9);

         BTLOF::CrsGraphType::nonconst_global_inds_host_view_type indicesView("indices", 10);
         size_t numIndices;

         // Take degree of freedom in middle of mesh: Then look at ghosted graph
         gGraph->getGlobalRowCopy(map->getGlobalElement(3), indicesView, numIndices);

         std::vector<int> indices;
         for (int i = 0; i < numIndices; i++)
         {
            indices.push_back(indicesView(i));
         }

         indices.resize(numIndices);
         std::sort(indices.begin(), indices.end());
         if (numProcs == 2 && myRank == 0)
         {
            TEST_EQUALITY(numIndices, 6);

            std::vector<int> compare(6);
            compare[0] = 0;
            compare[1] = 1;
            compare[2] = 2;
            compare[3] = 3;
            compare[4] = 4;
            compare[5] = 5;

            TEST_EQUALITY(compare.size(), indices.size());
            for (int i = 0; i < numIndices; i++)
            {
               TEST_EQUALITY(compare[i], indices[i]);
            }
         }
         else if (numProcs == 2 && myRank == 1)
         {
            TEST_EQUALITY(numIndices, 6);

            std::vector<int> compare(6);
            compare[0] = 1;
            compare[1] = 3;
            compare[2] = 5;
            compare[3] = 6;
            compare[4] = 7;
            compare[5] = 8;

            TEST_EQUALITY(compare.size(), indices.size());
            for (int i = 0; i < numIndices; i++)
            {
               TEST_EQUALITY(compare[i], indices[i]);
            }
         }
      }

      // test graph
      {
         TEST_ASSERT(graph->isFillComplete());

         TEST_EQUALITY(graph->getLocalNumRows(), (int)owned.size());
         TEST_EQUALITY(graph->getLocalMaxNumRowEntries(), myRank == 0 ? 9 : 6);
      }
   }

}
