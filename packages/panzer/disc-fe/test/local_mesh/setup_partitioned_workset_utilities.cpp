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

#include "PanzerCore_config.hpp"

#include "Panzer_Workset.hpp"

#include "Panzer_WorksetNeeds.hpp"
#include "Panzer_WorksetDescriptor.hpp"
#include "Panzer_BasisDescriptor.hpp"
#include "Panzer_IntegrationDescriptor.hpp"

#include "Panzer_LocalMeshInfo.hpp"
#include "Panzer_WorksetUtilities.hpp"

#include "Panzer_UnitTest_LocalMeshUtilities.hpp"

#include <vector>
#include <string>
#include <iostream>

namespace panzer {

TEUCHOS_UNIT_TEST(setupPartitionedWorksetUtilities, basic)
{

  // Useful things
  const panzer::BasisDescriptor basis(1,"HGrad");
  const panzer::IntegrationDescriptor integrator(1,panzer::IntegrationDescriptor::VOLUME);

  // Test empty mesh
  {
    panzer::LocalMeshInfo empty_mesh;

    // It shouldn't return any worksets
    auto worksets = buildWorksets(empty_mesh, panzer::blockGhostedDescriptor("block"));

    TEST_ASSERT(not worksets.is_null());
    TEST_EQUALITY(worksets->size(),0);

  }

  // Generate a local mesh info object

  Teuchos::RCP<panzer::LocalMeshInfo> mesh_info = generateLocalMeshInfo();
  {
    // Should get a single workset
    auto worksets = buildWorksets(*mesh_info, panzer::blockDescriptor("block0"));

    TEST_ASSERT(not worksets.is_null());
    TEST_EQUALITY(worksets->size(), 1);
  }

  // Test the worksets

  // Non-existant block
  {
    // It shouldn't return any worksets
    auto worksets = buildWorksets(*mesh_info, panzer::blockGhostedDescriptor("block32"));

    TEST_ASSERT(not worksets.is_null());
    TEST_EQUALITY(worksets->size(), 0);

  }

  // Non-existant sideset
  {
    // It shouldn't return any worksets
    auto worksets = buildWorksets(*mesh_info, panzer::sidesetGhostedDescriptor("block32","sideset0"));

    TEST_ASSERT(not worksets.is_null());

    // Nothing should have been found
    TEST_EQUALITY(worksets->size(), 0);

  }
  {
    // It shouldn't return any worksets
    auto worksets = buildWorksets(*mesh_info, panzer::sidesetGhostedDescriptor("block0","sideset32"));

    TEST_ASSERT(not worksets.is_null());

    // Nothing should have been found
    TEST_EQUALITY(worksets->size(), 0);

  }
  {
    // It shouldn't return any worksets
    auto worksets = buildWorksets(*mesh_info, panzer::sidesetGhostedDescriptor("block0","sideset1"));

    TEST_ASSERT(not worksets.is_null());

    // Nothing should have been found
    TEST_EQUALITY(worksets->size(), 0);
  }

  // Existing block partition
  {

    // It shouldn't return any worksets
    auto worksets = buildWorksets(*mesh_info, panzer::blockGhostedDescriptor("block0"));

    TEST_ASSERT(not worksets.is_null());

    // Only one partition since we requested a full partition
    TEST_EQUALITY(worksets->size(), 1);

    const auto & workset = (*worksets)[0];
    TEST_EQUALITY(workset.numOwnedCells(),2);
    TEST_EQUALITY(workset.numGhostCells(),1);
    TEST_EQUALITY(workset.numVirtualCells(),1);
  }

  // Existing block, but don't partition
  {
    // Should only get a single workset with all elements
    auto worksets = buildWorksets(*mesh_info, panzer::blockDescriptor("block0"));

    // Should have found something
    TEST_ASSERT(not worksets.is_null());
    TEST_EQUALITY(worksets->size(), 1);

    // No partition means no ghost or virtual cells
    const auto & workset = (*worksets)[0];
    TEST_EQUALITY(workset.numOwnedCells(),2);
    TEST_EQUALITY(workset.numGhostCells(),0);
    TEST_EQUALITY(workset.numVirtualCells(),0);
  }

  // Existing sideset partition
  {

    // It shouldn't return any worksets
    auto worksets = buildWorksets(*mesh_info, panzer::sidesetGhostedDescriptor("block0","sideset0"));

    TEST_INEQUALITY(worksets,Teuchos::null);

    // Only one partition since we requested a full partition
    TEST_EQUALITY(worksets->size(), 1);

    const auto & workset = (*worksets)[0];
    TEST_EQUALITY(workset.numOwnedCells(),1);
    TEST_EQUALITY(workset.numGhostCells(),0);
    TEST_EQUALITY(workset.numVirtualCells(),1);
  }
  {
    // It shouldn't return any worksets
    auto worksets = buildWorksets(*mesh_info, panzer::sidesetGhostedDescriptor("block1","sideset2"));

    TEST_INEQUALITY(worksets,Teuchos::null);

    // Only one partition since we requested a full partition
    TEST_EQUALITY(worksets->size(), 1);

    const auto & workset = (*worksets)[0];
    TEST_EQUALITY(workset.numOwnedCells(),1);
    TEST_EQUALITY(workset.numGhostCells(),1);
    TEST_EQUALITY(workset.numVirtualCells(),0);
  }

  // Existing block double partition
  {
    // We want two partitions, each with a size of 1 cell
    auto worksets = buildWorksets(*mesh_info, panzer::blockGhostedDescriptor("block1",1));

    TEST_INEQUALITY(worksets,Teuchos::null);

    // Two partitions
    TEST_EQUALITY(worksets->size(), 2);

    {
      const auto & workset = (*worksets)[0];
      TEST_EQUALITY(workset.numOwnedCells(),1);

      // Either both extra cells are ghosts, or one is virtual and one is ghost
      TEST_EQUALITY(workset.numGhostCells()+workset.numVirtualCells(),2);
    }

    {
      const auto & workset = (*worksets)[1];
      TEST_EQUALITY(workset.numOwnedCells(),1);

      // Either both extra cells are ghosts, or one is virtual and one is ghost
      TEST_EQUALITY(workset.numGhostCells()+workset.numVirtualCells(),2);
    }
  }

  // Existing sideset double partition, but only one partition is available
  {
    // We want two partitions, each with a size of 1 cell
    auto worksets = buildWorksets(*mesh_info, panzer::sidesetGhostedDescriptor("block1","sideset1",1));

    TEST_INEQUALITY(worksets,Teuchos::null);

    // Only one partition should have been built
    TEST_EQUALITY(worksets->size(), 1);

    {
      const auto & workset = (*worksets)[0];
      TEST_EQUALITY(workset.numOwnedCells(),1);
      TEST_EQUALITY(workset.numGhostCells(),0);
      TEST_EQUALITY(workset.numVirtualCells(),1);
    }
  }


}

} // end namespace panzer
