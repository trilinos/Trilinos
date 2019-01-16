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

using Teuchos::RCP;
using Teuchos::rcp;

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_STK_Version.hpp"
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_WorksetUtilities.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_LocalMeshUtilities.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_LocalMeshInfo.hpp"
#include "Panzer_WorksetDescriptor.hpp"

#include "user_app_EquationSetFactory.hpp"

namespace panzer {

  void testIpMatch(const panzer::Workset& d0,
                   const panzer::Workset& d1,
                   Teuchos::FancyOStream& out,
                   bool& success);

  TEUCHOS_UNIT_TEST(workset_builder, stk_edge)
  {

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
    pl->set("X Blocks",2);
    pl->set("Y Blocks",1);
    pl->set("X Elements",2);  // in each block
    pl->set("Y Elements",2);  // in each block

    // This mesh will have element IDs
    //    4 5 6 7
    //    0 1 2 3

    // There are 4 cases to test for
    // Test    Rank    Interface       Result
    //  0       0        a -> b       empty
    //  1       1        a -> b       not empty
    //  2       0        b -> a       not empty
    //  3       1        b -> a       empty

    // Here: a = element_blocks[0], b = element_blocks[1]

    int myRank=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    panzer_stk::SquareQuadMeshFactory factory;
    factory.setParameterList(pl);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);
    mesh->writeToExodus("test.exo");

    auto mesh_info_rcp = panzer_stk::generateLocalMeshInfo(*mesh);
    auto & mesh_info = *mesh_info_rcp;

    std::vector<std::string> element_blocks;
    mesh->getElementBlockNames(element_blocks);

    // Tests 0 and 1 (a->b)
    {
      std::string sideset = "vertical_0";
      auto worksets = panzer::buildWorksets(mesh_info, WorksetDescriptor(element_blocks[0],element_blocks[1], sideset, sideset));

      if(myRank==0) {

        // Test 0

        TEST_EQUALITY(worksets->size(), 0);


      } else {

        // Test 1

        TEST_EQUALITY(worksets->size(), 1);

        const panzer::Workset& workset = (*worksets)[0];

        TEST_EQUALITY(workset.numCells(),2);
        TEST_EQUALITY(workset.getSubcellDimension(),1);
        TEST_EQUALITY(workset.size(),2);

        // this is identical to workset(0)
        TEST_EQUALITY(workset.getSubcellIndex(), 1);
        TEST_EQUALITY(workset.getElementBlock(), "eblock-0_0");
        TEST_EQUALITY(workset.getLocalCellIDs().size(),2);
        TEST_EQUALITY(workset.getLocalCellIDs()(0),2);
        TEST_EQUALITY(workset.getLocalCellIDs()(1),0);

        TEST_EQUALITY(workset(0).getSubcellIndex(), 1);
        TEST_EQUALITY(workset(0).getElementBlock(), "eblock-0_0");
        TEST_EQUALITY(workset(0).getLocalCellIDs().size(),2);
        TEST_EQUALITY(workset(0).getLocalCellIDs()(0),2);
        TEST_EQUALITY(workset(0).getLocalCellIDs()(1),0);

        TEST_EQUALITY(workset(1).getSubcellIndex(), 3);
        TEST_EQUALITY(workset(1).getElementBlock(), "eblock-1_0");
        TEST_EQUALITY(workset(1).getLocalCellIDs().size(),2);
        TEST_EQUALITY(workset(1).getLocalCellIDs()(0),7);
        TEST_EQUALITY(workset(1).getLocalCellIDs()(1),5);

        testIpMatch(workset(0), workset(1), out, success);

      }

    }

    // Tests 2 and 3 (b->a)

    {
      std::string sideset = "vertical_0";
      auto worksets = panzer::buildWorksets(mesh_info, WorksetDescriptor(element_blocks[1],element_blocks[0], sideset, sideset));

      if(myRank==0) {

        // Test 2

        TEST_EQUALITY(worksets->size(),1);

        panzer::Workset& workset = *worksets->begin();

        TEST_EQUALITY(workset.numCells(),2);
        TEST_EQUALITY(workset.getSubcellDimension(),1);
        TEST_EQUALITY(workset.size(),2);

        // this is identical to details[0]
        TEST_EQUALITY(workset.getSubcellIndex(), 3);
        TEST_EQUALITY(workset.getElementBlock(), "eblock-1_0");
        TEST_EQUALITY(workset.getLocalCellIDs().size(),2);
        TEST_EQUALITY(workset.getLocalCellIDs()(0),3);
        TEST_EQUALITY(workset.getLocalCellIDs()(1),1);

        TEST_EQUALITY(workset(0).getSubcellIndex(), 3);
        TEST_EQUALITY(workset(0).getElementBlock(), "eblock-1_0");
        TEST_EQUALITY(workset(0).getLocalCellIDs().size(),2);
        TEST_EQUALITY(workset(0).getLocalCellIDs()(0),3);
        TEST_EQUALITY(workset(0).getLocalCellIDs()(1),1);

        TEST_EQUALITY(workset(1).getSubcellIndex(), 1);
        TEST_EQUALITY(workset(1).getElementBlock(), "eblock-0_0");
        TEST_EQUALITY(workset(1).getLocalCellIDs()(0),6);
        TEST_EQUALITY(workset(1).getLocalCellIDs()(1),4);

        testIpMatch(workset(0), workset(1), out, success);

      } else {

        // Test 3

        TEST_EQUALITY(worksets->size(),0);

      }
    }
    
  }

  void testIpMatch(const panzer::Workset& d0,
                   const panzer::Workset& d1,
                   Teuchos::FancyOStream& out,
                   bool& success)
  {
    const auto & coords0 = d0.getIntegrationValues(IntegrationDescriptor(4, IntegrationDescriptor::SIDE, d0.getSubcellIndex())).ip_coordinates;
    const auto & coords1 = d1.getIntegrationValues(IntegrationDescriptor(4, IntegrationDescriptor::SIDE, d1.getSubcellIndex())).ip_coordinates;

    TEST_EQUALITY(d0.numCells(), d1.numCells());
    {
      const std::size_t num_ip = coords0.extent(1);
      const std::size_t num_dim = coords0.extent(2);
      for (index_t cell = 0; cell < d0.numCells(); ++cell)
        for (std::size_t ip = 0; ip < num_ip; ++ip)
          for (std::size_t dim = 0; dim < num_dim; ++dim)
            TEST_EQUALITY(coords0(cell, ip, dim), coords1(cell, ip, dim));
    }
  }

}
