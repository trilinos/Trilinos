// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <math.h>                                     // for sqrt
#include <cstddef>                                    // for size_t
#include <cstdint>                                    // for int64_t, uint64_t
#include <limits>                                     // for numeric_limits
#include <stdexcept>                                  // for logic_error
#include <string>                                     // for string, basic_s...
#include <typeinfo>                                   // for type_info
#include <utility>                                    // for move, pair
#include <vector>                                     // for vector, swap
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include "mpi.h"                                      // for MPI_COMM_WORLD

#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/Field.hpp>
#include "stk_mesh/base/Entity.hpp"                   // for Entity
#include "stk_mesh/base/FieldBase.hpp"                // for field_data, Fie...
#include "stk_mesh/base/FieldState.hpp"               // for FieldState
#include "stk_mesh/base/Part.hpp"                     // for Part
#include "stk_mesh/base/Selector.hpp"                 // for Selector, opera...
#include "stk_mesh/base/SkinBoundary.hpp"             // for create_all_sides
#include "stk_mesh/base/Types.hpp"                    // for EntityRank, Ent...
#include "stk_topology/topology.hpp"                  // for topology, topol...
#include "stk_unit_test_utils/MeshFixture.hpp"        // for MeshFixtureNoTest
#include "stk_unit_test_utils/TextMesh.hpp"           // for setup_text_mesh
#include "stk_unit_test_utils/getOption.h"            // for get_command_lin...
#include "stk_util/diag/String.hpp"                   // for String
#include "stk_util/parallel/Parallel.hpp"             // for parallel_machin...
#include "stk_util/util/ReportHandler.hpp"            // for ThrowRequireMsg
#include <stk_transfer_util/Patch.hpp>
#include <stk_transfer_util/EntityCentroidRecoverField.hpp>

namespace {

class PatchTester : public stk::unit_test_util::MeshFixtureNoTest, public ::testing::Test {
 public:
  PatchTester()
    : stk::unit_test_util::MeshFixtureNoTest(3)
  {
    setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
  }

  ~PatchTester()
  {
    delete m_patch;
    delete m_patchFilter;
    delete m_selector;
  }

 protected:
  stk::transfer::EntityPatchFilter* m_patchFilter = nullptr;
  stk::transfer::Patch<stk::transfer::EntityPatchFilter>* m_patch = nullptr;

  stk::mesh::Part* m_part = nullptr;
  stk::mesh::Selector *m_selector = nullptr;

  void create_mesh(unsigned numElemsPerAxis = 3)
  {
    std::string meshSpec("generated:" + std::to_string(numElemsPerAxis) + "x" + std::to_string(numElemsPerAxis) + "x" +
                         std::to_string(numElemsPerAxis));
    stk::io::fill_mesh(meshSpec, get_bulk());
  }

  void create_patch(stk::mesh::EntityId seedId, stk::mesh::EntityRank rank)
  {
    delete m_patchFilter;
    delete m_patch;
    delete m_selector;

    m_part = get_meta().get_part("block_1");
    m_selector = new stk::mesh::Selector(get_meta().universal_part());

    stk::mesh::Entity patchSeed = get_bulk().get_entity(rank, seedId);
    ASSERT_TRUE(get_bulk().is_valid(patchSeed));
    m_patchFilter = new stk::transfer::EntityPatchFilter(m_part, m_selector);
    m_patch = new stk::transfer::LinearPatch<stk::transfer::EntityPatchFilter>(get_bulk(), patchSeed, *m_patchFilter, *m_selector);
  }

  void create_cubic_patch(stk::mesh::EntityId seedId, stk::mesh::EntityRank rank)
  {
    delete m_patchFilter;
    delete m_patch;
    delete m_selector;

    m_part = get_meta().get_part("block_1");
    m_selector = new stk::mesh::Selector(get_meta().universal_part());

    stk::mesh::Entity patchSeed = get_bulk().get_entity(rank, seedId);
    ASSERT_TRUE(get_bulk().is_valid(patchSeed));
    m_patchFilter = new stk::transfer::EntityPatchFilter(m_part, m_selector);
    m_patch = new stk::transfer::CubicPatch<stk::transfer::EntityPatchFilter>(get_bulk(), patchSeed, *m_patchFilter, *m_selector);
  }
};

TEST_F(PatchTester, nodePatch_middleOf2x2x2Mesh)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 1) {
    return;
  }

  create_mesh(2);
  create_patch(stk::mesh::EntityId(14u), stk::topology::NODE_RANK);

  EXPECT_EQ(27u, m_patch->get_patch_entities().size());
  EXPECT_EQ( 8u, m_patch->get_patch_elements().size());
  EXPECT_EQ(27u, m_patch->get_patch_nodes().size());
}

TEST_F(PatchTester, nodePatch_cornerOf2x2x2Mesh)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 1) {
    return;
  }

  create_mesh(2);
  create_patch(stk::mesh::EntityId(3u), stk::topology::NODE_RANK);

  EXPECT_EQ( 8u, m_patch->get_patch_entities().size());
  EXPECT_EQ( 1u, m_patch->get_patch_elements().size());
  EXPECT_EQ( 8u, m_patch->get_patch_nodes().size());
}

TEST_F(PatchTester, elementPatch_middleOf3x3x3Mesh)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 1) {
    return;
  }

  create_mesh(3);
  create_patch(stk::mesh::EntityId(14u), stk::topology::ELEM_RANK);

  EXPECT_EQ(27u, m_patch->get_patch_entities().size());
  EXPECT_EQ(27u, m_patch->get_patch_elements().size());
  EXPECT_EQ(64u, m_patch->get_patch_nodes().size());
}

TEST_F(PatchTester, elementPatch_cornerOf3x3x3Mesh)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 1) {
    return;
  }

  create_mesh(3);
  create_patch(stk::mesh::EntityId(3u), stk::topology::ELEM_RANK);

  EXPECT_EQ( 8u, m_patch->get_patch_entities().size());
  EXPECT_EQ( 8u, m_patch->get_patch_elements().size());
  EXPECT_EQ(27u, m_patch->get_patch_nodes().size());
}

TEST_F(PatchTester, elementPatch_edgeOf3x3x3Mesh)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 1) {
    return;
  }

  create_mesh(3);
  create_patch(stk::mesh::EntityId(2u), stk::topology::ELEM_RANK);

  EXPECT_EQ(12u, m_patch->get_patch_entities().size());
  EXPECT_EQ(12u, m_patch->get_patch_elements().size());
  EXPECT_EQ(36u, m_patch->get_patch_nodes().size());
}

TEST_F(PatchTester, cubicElementPatch_middleOf5x5x5Mesh)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  if(stk::parallel_machine_size(comm) != 1) {
    return;
  }

  create_mesh(5);
  create_cubic_patch(stk::mesh::EntityId(63u), stk::topology::ELEM_RANK);

  EXPECT_EQ(125u, m_patch->get_patch_entities().size());
  EXPECT_EQ(125u, m_patch->get_patch_elements().size());
  EXPECT_EQ(216u, m_patch->get_patch_nodes().size());
}

}
