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
#include "Teuchos_ParameterList.hpp"

#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_LocalMeshUtilities.hpp"

#include "Panzer_LocalMeshInfo.hpp"
#include "PanzerSTK_UnitTest_BuildMesh.hpp"

#include <string>

namespace panzer_stk {

Teuchos::RCP<panzer::LocalMeshInfo>
buildLocalMeshInfo(const std::vector<int> & N,   // Cells per dimension
                   const std::vector<int> & B,   // Blocks per dimension
                   const std::vector<double> &L) // Domain length per block
{
  Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildMesh(N,B,L);
  Teuchos::RCP<panzer::LocalMeshInfo> mesh_info(new panzer::LocalMeshInfo);
  generateLocalMeshInfo(*mesh,*mesh_info);
  return mesh_info;
}


TEUCHOS_UNIT_TEST(localMeshUtilities, basic)
{

  // This is not an all encompassing test.
  // It just makes sure that the local mesh infos are being generated and sized correctly.

  // Make sure if fails when you pass in a null mesh
  {
    panzer_stk::STK_Interface mesh;
    panzer::LocalMeshInfo mesh_info;

    TEST_THROW((generateLocalMeshInfo(mesh, mesh_info)),std::logic_error);
  }

  // 1D Mesh test
  {
    auto mesh_info = buildLocalMeshInfo({3},{2},{2.});

    // Make sure there are two blocks (eblock-0, and eblock-1)
    TEST_EQUALITY(mesh_info->element_blocks.size(), 2);

    // Make sure there are two sidesets (left and right)
    TEST_EQUALITY(mesh_info->sidesets.size(), 2);

    // Test block 0
    {
      const auto & block = mesh_info->element_blocks.at("eblock-0");

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 3);
      TEST_EQUALITY(block.num_ghstd_cells, 1);
      TEST_EQUALITY(block.num_virtual_cells, 1);
    }
  }


  // 2D Mesh test
  {
    auto mesh_info = buildLocalMeshInfo({3,5},{4,3},{2.,6.});

    // Make sure there are blocks
    TEST_EQUALITY(mesh_info->element_blocks.size(), 4*3);

    // Make sure there are sidesets
    TEST_EQUALITY(mesh_info->sidesets.size(), 4*3);
    TEST_EQUALITY(mesh_info->sidesets["eblock-0_0"].size(), (4+1) + (3+1));

    // Test block 0,0
    {
      const auto & block = mesh_info->element_blocks.at("eblock-0_0");

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 3*5);
      TEST_EQUALITY(block.num_ghstd_cells, 5+3);
      TEST_EQUALITY(block.num_virtual_cells, 5+3);
    }

    // Test block 1,1
    {
      const auto & block = mesh_info->element_blocks.at("eblock-1_1");

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 3*5);
      TEST_EQUALITY(block.num_ghstd_cells, 2*5+2*3);
      TEST_EQUALITY(block.num_virtual_cells, 0);
    }
  }

  // 3D Mesh test
  {
    auto mesh_info = buildLocalMeshInfo({3,5,2},{4,3,5},{2.,6.,1.});

    // Make sure there are blocks
    TEST_EQUALITY(mesh_info->element_blocks.size(), 4*3*5);

    // Make sure there are sidesets
    TEST_EQUALITY(mesh_info->sidesets.size(), 4*3*5);

    // Note the naming scheme for 3D is different from 2D - there are no internal sidesets
    //TEST_EQUALITY(mesh_info->sidesets["eblock-0_0_0"].size(), (4+1) + (3+1) + (5+1));
    TEST_EQUALITY(mesh_info->sidesets["eblock-0_0_0"].size(), 6);

    // Test block 0,0,0
    {
      const auto & block = mesh_info->element_blocks.at("eblock-0_0_0");

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 3*5*2);
      TEST_EQUALITY(block.num_ghstd_cells, 3*5+2*3+2*5);
      TEST_EQUALITY(block.num_virtual_cells, 3*5+2*3+2*5);
    }

    // Test block 1,1,1
    {
      const auto & block = mesh_info->element_blocks.at("eblock-1_1_1");

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 3*5*2);
      TEST_EQUALITY(block.num_ghstd_cells, 2*3*5+2*2*3+2*2*5);
      TEST_EQUALITY(block.num_virtual_cells, 0);
    }
  }

}

}
