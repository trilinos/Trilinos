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
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/Zoltan2ParallelGraph.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_balance/internal/Balancer.hpp>
#include <stk_balance/mesh/BalanceMesh.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <vector>
#include <string>

class StkBalanceDecomposition : public stk::unit_test_util::MeshFixture
{
protected:
  StkBalanceDecomposition()
    : MeshFixture()
  { }

  ~StkBalanceDecomposition() = default;

  void setup_initial_mesh(const std::string & inputMeshFile)
  {
    setup_mesh(inputMeshFile, stk::mesh::BulkData::AUTO_AURA);
  }

  void unbalance_mesh()
  {
    stk::mesh::EntityProcVec entitiesToChange;
    stk::mesh::EntityVector entities;
    get_bulk().get_entities(stk::topology::ELEM_RANK, get_meta().locally_owned_part(), entities);
    for (const stk::mesh::Entity & element : entities) {
      entitiesToChange.push_back(stk::mesh::EntityProc(element, 0));
    }

    entities.clear();
    get_bulk().get_entities(stk::topology::NODE_RANK, get_meta().locally_owned_part(), entities);
    for (const stk::mesh::Entity & node : entities) {
      entitiesToChange.push_back(stk::mesh::EntityProc(node, 0));
    }

    get_bulk().change_entity_owner(entitiesToChange);
  }

  void balance_mesh(const std::vector<stk::mesh::Selector> & selectors)
  {
    stk::balance::GraphCreationSettings balanceSettings;
    stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;
    stk::balance::Balancer balancer(balanceSettings);
    stk::balance::BalanceMesh balanceMesh(get_bulk());
    balancer.balance(balanceMesh, selectors);
    stk::EnvData::instance().m_outputP0 = &std::cout;
  }

  void test_mesh_element_distribution(const std::vector<int> & expectedElemsPerProc)
  {
    const int numElems = stk::mesh::count_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part());

    EXPECT_EQ(numElems, expectedElemsPerProc[get_parallel_rank()]);
  }
};


TEST_F(StkBalanceDecomposition, 4Elem1ProcMesh_EntireDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  setup_initial_mesh("generated:1x1x4");
  balance_mesh({get_meta().universal_part()});

  test_mesh_element_distribution({4});
}

TEST_F(StkBalanceDecomposition, 6Elem2ProcMesh_EntireDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  unbalance_mesh();
  balance_mesh({get_meta().universal_part()});

  test_mesh_element_distribution({3, 3});
}

TEST_F(StkBalanceDecomposition, 6Elem2ProcMesh_OneElem)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}});
  balance_mesh({*get_meta().get_part("partA")});

  test_mesh_element_distribution({2, 4});  // Moved the one elem from p0 to p1
}

TEST_F(StkBalanceDecomposition, 6Elem2ProcMesh_TwoElems)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {2, "partA"}});
  balance_mesh({*get_meta().get_part("partA")});

  test_mesh_element_distribution({2, 4});  // Moved one of the two elems from p0 to p1
}

TEST_F(StkBalanceDecomposition, 6Elem2ProcMesh_TwoElemsEachProc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {2, "partA"}, {4, "partA"}, {5, "partA"}});
  balance_mesh({*get_meta().get_part("partA")});

  test_mesh_element_distribution({3, 3});
}

TEST_F(StkBalanceDecomposition, 6Elem2ProcMesh_TwoElemsAcrossProcBoundary)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{3, "partA"}, {4, "partA"}});
  balance_mesh({*get_meta().get_part("partA")});

  test_mesh_element_distribution({3, 3});
}

TEST_F(StkBalanceDecomposition, 6Elem2ProcMesh_TwoElemsEachProcWithOneAdjacentToProcBoundary)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {3, "partA"}, {4, "partA"}, {6, "partA"}});
  balance_mesh({*get_meta().get_part("partA")});

  test_mesh_element_distribution({3, 3});
}
