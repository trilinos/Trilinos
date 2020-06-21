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
#include "Teuchos_DefaultComm.hpp"
#include "Tpetra_Core.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_LocalMeshUtilities.hpp"

#include "Panzer_LocalMeshInfo.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"

#include <string>

namespace panzer_stk {

Teuchos::RCP<panzer::LocalMeshInfo>
buildLocalMeshInfo(const std::string & filename,
                   Teuchos::RCP<const Teuchos::Comm<int> > comm)
{
  Teuchos::RCP<const Teuchos::MpiComm<int> > mpi_comm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm,true);
  MPI_Comm raw_comm = *(mpi_comm->getRawMpiComm());
  STK_ExodusReaderFactory factory(filename);
  Teuchos::RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(raw_comm);
  Teuchos::RCP<panzer::LocalMeshInfo> mesh_info(new panzer::LocalMeshInfo);
  generateLocalMeshInfo(*mesh,*mesh_info);
  return mesh_info;
}

TEUCHOS_UNIT_TEST(parallelLocalMeshUtilities, basic)
{

  // The goal of this test is to make sure that when a mesh is divided on processor boundaries, both chunks of the mesh see the sideset properly

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  TEST_ASSERT(comm->getSize() == 2);

  auto mesh_info = buildLocalMeshInfo("test_mesh.exo", comm);

  const std::set<std::string> block_names {"left_block", "right_block"};
  const std::set<std::string> sideset_names {"top", "bottom", "left_wall", "right_wall", "center_wall", "front", "back"};

  // All ranks will see the same names for blocks and sidesets - even if they don't exist locally
  TEST_EQUALITY(mesh_info->element_blocks.size(), block_names.size());
  for(const auto & pr : mesh_info->element_blocks){
    TEST_ASSERT(block_names.find(pr.first) != block_names.end());
  }
  TEST_EQUALITY(mesh_info->sidesets.size(), block_names.size());
  for(const auto & pr : mesh_info->sidesets){
    TEST_ASSERT(block_names.find(pr.first) != block_names.end());
    const auto & sidesets = pr.second;

    TEST_EQUALITY(sidesets.size(), sideset_names.size());
    for(const auto & pr2 : pr.second){
      TEST_ASSERT(sideset_names.find(pr2.first) != sideset_names.end());
    }
  }

#define TEST_BLOCK(BLOCK, NUM_OWNED_CELLS, NUM_GHOST_CELLS, NUM_VIRTUAL_CELLS) \
    { \
      const auto & block = mesh_info->element_blocks[BLOCK]; \
      TEST_EQUALITY(block.num_owned_cells, NUM_OWNED_CELLS); \
      TEST_EQUALITY(block.num_ghstd_cells, NUM_GHOST_CELLS); \
      TEST_EQUALITY(block.num_virtual_cells, NUM_VIRTUAL_CELLS); \
    }

#define TEST_SIDESET(BLOCK, SIDESET, NUM_OWNED_CELLS, NUM_GHOST_CELLS, NUM_VIRTUAL_CELLS) \
    { \
      const auto & block = mesh_info->sidesets[BLOCK][SIDESET]; \
      TEST_EQUALITY(block.num_owned_cells, NUM_OWNED_CELLS); \
      TEST_EQUALITY(block.num_ghstd_cells, NUM_GHOST_CELLS); \
      TEST_EQUALITY(block.num_virtual_cells, NUM_VIRTUAL_CELLS); \
    }

  if(comm->getRank() == 0){

    // Test Blocks
    TEST_BLOCK("left_block",0,0,0);
    TEST_BLOCK("right_block",8,4,20);

    // Test Sidesets
    TEST_SIDESET("left_block","top",0,0,0);
    TEST_SIDESET("left_block","bottom",0,0,0);
    TEST_SIDESET("left_block","front",0,0,0);
    TEST_SIDESET("left_block","back",0,0,0);
    TEST_SIDESET("left_block","left_wall",0,0,0);
    TEST_SIDESET("left_block","center_wall",0,0,0);
    TEST_SIDESET("left_block","right_wall",0,0,0);

    TEST_SIDESET("right_block","top",4,0,4);
    TEST_SIDESET("right_block","bottom",4,0,4);
    TEST_SIDESET("right_block","front",4,0,4);
    TEST_SIDESET("right_block","back",4,0,4);
    TEST_SIDESET("right_block","left_wall",0,0,0);
    TEST_SIDESET("right_block","center_wall",4,4,0);
    TEST_SIDESET("right_block","right_wall",4,0,4);
  } else {

    // Test Blocks
    TEST_BLOCK("left_block",8,4,20);
    TEST_BLOCK("right_block",0,0,0);

    // Test Sidesets
    TEST_SIDESET("left_block","top",4,0,4);
    TEST_SIDESET("left_block","bottom",4,0,4);
    TEST_SIDESET("left_block","front",4,0,4);
    TEST_SIDESET("left_block","back",4,0,4);
    TEST_SIDESET("left_block","left_wall",4,0,4);
    TEST_SIDESET("left_block","center_wall",4,4,0);
    TEST_SIDESET("left_block","right_wall",0,0,0);

    TEST_SIDESET("right_block","top",0,0,0);
    TEST_SIDESET("right_block","bottom",0,0,0);
    TEST_SIDESET("right_block","front",0,0,0);
    TEST_SIDESET("right_block","back",0,0,0);
    TEST_SIDESET("right_block","left_wall",0,0,0);
    TEST_SIDESET("right_block","center_wall",0,0,0);
    TEST_SIDESET("right_block","right_wall",0,0,0);
  }

#undef TEST_BLOCK
#undef TEST_SIDESET

}

}
