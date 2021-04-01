#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/Zoltan2ParallelGraph.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <vector>
#include <string>

struct ElementAndPart {
  stk::mesh::EntityId id;
  std::string partName;
};

class ZoltanGraphGeneration : public stk::unit_test_util::MeshFixture
{
protected:
  ZoltanGraphGeneration()
    : MeshFixture()
  { }

  ~ZoltanGraphGeneration() override {
  }

  void setup_initial_mesh(const std::string & inputMeshFile)
  {
    setup_mesh(inputMeshFile, stk::mesh::BulkData::AUTO_AURA);
  }

  void fill_zoltan_graph(const stk::mesh::Selector & selector, const stk::balance::BalanceSettings & balanceSettings)
  {
    testing::internal::CaptureStdout();
    stk::balance::internal::createZoltanParallelGraph(get_bulk(), selector, get_comm(), balanceSettings, m_graph);
    testing::internal::GetCapturedStdout();
  }

  void fill_zoltan_graph_for_decomp(const stk::mesh::Selector & selector)
  {
    stk::balance::GraphCreationSettings balanceSettings;
    fill_zoltan_graph(selector, balanceSettings);
  }

  void fill_zoltan_graph_for_coloring(const stk::mesh::Selector & selector)
  {
    stk::balance::BasicColoringSettings balanceSettings;
    fill_zoltan_graph(selector, balanceSettings);
  }

  struct ElemConnectivity {
    stk::mesh::EntityId id;
    stk::mesh::EntityIdVector connectedIds;
  };

  void test_connectivity_for_each_owned_element(const std::vector<ElemConnectivity> & expectedElemConnectivity)
  {
    std::vector<BalanceLocalNumber> expectedOffsets;
    std::vector<BalanceGlobalNumber> expectedAdjacency;

    int adjacencyIndex = 0;
    for (const ElemConnectivity & elemConnectivity : expectedElemConnectivity) {
      expectedOffsets.push_back(adjacencyIndex);

      for (const stk::mesh::EntityId & connectedId : elemConnectivity.connectedIds) {
        expectedAdjacency.push_back(connectedId);
        ++adjacencyIndex;
      }
    }
    expectedOffsets.push_back(adjacencyIndex);

    const std::vector<BalanceLocalNumber> & offsets = m_graph.get_offsets();
    const std::vector<BalanceGlobalNumber> & adjacency = m_graph.get_adjacency();

    EXPECT_EQ(offsets, expectedOffsets);
    EXPECT_EQ(adjacency, expectedAdjacency);
  }

  void print_offsets_and_adjacency() {
    const std::vector<BalanceLocalNumber> & offsets = m_graph.get_offsets();
    const std::vector<BalanceGlobalNumber> & adjacency = m_graph.get_adjacency();

    for (int proc = 0; proc < get_parallel_size(); ++proc) {
      if (proc == get_parallel_rank()) {
        std::cout << "### parallel_rank = " << proc << std::endl << std::endl;
        for (unsigned i = 0; i < offsets.size(); ++i) {
          std::cout << "offset[" << i << "] = " << offsets[i] << std::endl;
        }
        std::cout << std::endl;

        for (unsigned i = 0; i < adjacency.size(); ++i) {
          std::cout << "adjacency[" << i << "] = " << adjacency[i] << std::endl;
        }
        std::cout << std::endl;
      }
      MPI_Barrier(get_comm());
    }
  }

  Zoltan2ParallelGraph m_graph;
};


TEST_F(ZoltanGraphGeneration, 6Elem2ProcMesh_Coloring_EntireDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  fill_zoltan_graph_for_coloring(get_meta().universal_part());

  if (get_parallel_rank() == 0) {
    //  |    p0     |
    //  | 0 | 1 | 2 |
    test_connectivity_for_each_owned_element({ {0, {1}}, {1, {0, 2}}, {2, {1}} });
  }
  else if (get_parallel_rank() == 1) {
    //  |    p1     |
    //  | 0 | 1 | 2 |
    test_connectivity_for_each_owned_element({ {0, {1}}, {1, {0, 2}}, {2, {1}} });
  }
}

TEST_F(ZoltanGraphGeneration, 6Elem2ProcMesh_Coloring_OneElem)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}});
  fill_zoltan_graph_for_decomp(*get_meta().get_part("partA"));

  if (get_parallel_rank() == 0) {
    //  |    p0     |
    //  | 0 | 1 | 2 |
    test_connectivity_for_each_owned_element({ {0, {}} });
  }
  else if (get_parallel_rank() == 1) {
    //  |    p1     |
    //  | 0 | 1 | 2 |
    test_connectivity_for_each_owned_element({ });
  }
}

TEST_F(ZoltanGraphGeneration, 6Elem2ProcMesh_Coloring_TwoElems)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {2, "partA"}});
  fill_zoltan_graph_for_coloring(*get_meta().get_part("partA"));

  if (get_parallel_rank() == 0) {
    //  |    p0     |
    //  | 0 | 1 | 2 |
    test_connectivity_for_each_owned_element({ {0, {1}}, {1, {0}} });
  }
  else if (get_parallel_rank() == 1) {
    //  |    p1     |
    //  | 0 | 1 | 2 |
    test_connectivity_for_each_owned_element({ });
  }
}

TEST_F(ZoltanGraphGeneration, 6Elem2ProcMesh_Coloring_TwoElemsEachProc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {2, "partA"}, {4, "partA"}, {5, "partA"}});
  fill_zoltan_graph_for_coloring(*get_meta().get_part("partA"));

  if (get_parallel_rank() == 0) {
    //  |    p0     |
    //  | 0 | 1 | 2 |
    test_connectivity_for_each_owned_element({ {0, {1}}, {1, {0}} });
  }
  else if (get_parallel_rank() == 1) {
    //  |    p1     |
    //  | 0 | 1 | 2 |
    test_connectivity_for_each_owned_element({ {0, {1}}, {1, {0}} });
  }
}

TEST_F(ZoltanGraphGeneration, 6Elem2ProcMesh_Decomp_EntireDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  fill_zoltan_graph_for_decomp(get_meta().universal_part());

  if (get_parallel_rank() == 0) {
    //  |    p0     ||p1 |
    //  | 1 | 2 | 3 || 4 |
    test_connectivity_for_each_owned_element({ {1, {2}}, {2, {1, 3}}, {3, {2, 4}} });
  }
  else if (get_parallel_rank() == 1) {
    //  |p0 ||    p1     |
    //  | 3 || 4 | 5 | 6 |
    test_connectivity_for_each_owned_element({ {4, {3, 5}}, {5, {4, 6}}, {6, {5}} });
  }
}

TEST_F(ZoltanGraphGeneration, 6Elem2ProcMesh_Decomp_OneElem)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}});
  fill_zoltan_graph_for_decomp(*get_meta().get_part("partA"));

  if (get_parallel_rank() == 0) {
    //  |    p0     ||p1 |
    //  | 1 | 2 | 3 || 4 |
    test_connectivity_for_each_owned_element({ {1, {}} });
  }
  else if (get_parallel_rank() == 1) {
    //  |p0 ||    p1     |
    //  | 3 || 4 | 5 | 6 |
    test_connectivity_for_each_owned_element({ });
  }
}

TEST_F(ZoltanGraphGeneration, 6Elem2ProcMesh_Decomp_TwoElems)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {2, "partA"}});
  fill_zoltan_graph_for_decomp(*get_meta().get_part("partA"));

  if (get_parallel_rank() == 0) {
    //  |    p0     ||p1 |
    //  | 1 | 2 | 3 || 4 |
    test_connectivity_for_each_owned_element({ {1, {2}}, {2, {1}} });
  }
  else if (get_parallel_rank() == 1) {
    //  |p0 ||    p1     |
    //  | 3 || 4 | 5 | 6 |
    test_connectivity_for_each_owned_element({ });
  }
}

TEST_F(ZoltanGraphGeneration, 6Elem2ProcMesh_Decomp_TwoElemsEachProc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {2, "partA"}, {4, "partA"}, {5, "partA"}});
  fill_zoltan_graph_for_decomp(*get_meta().get_part("partA"));

  if (get_parallel_rank() == 0) {
    //  |    p0     ||p1 |
    //  | 1 | 2 | 3 || 4 |
    test_connectivity_for_each_owned_element({ {1, {2}}, {2, {1}} });
  }
  else if (get_parallel_rank() == 1) {
    //  |p0 ||    p1     |
    //  | 3 || 4 | 5 | 6 |
    test_connectivity_for_each_owned_element({ {4, {5}}, {5, {4}} });
  }
}

TEST_F(ZoltanGraphGeneration, 6Elem2ProcMesh_Decomp_TwoElemsAcrossProcBoundary)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{3, "partA"}, {4, "partA"}});
  fill_zoltan_graph_for_decomp(*get_meta().get_part("partA"));

  if (get_parallel_rank() == 0) {
    //  |    p0     ||p1 |
    //  | 1 | 2 | 3 || 4 |
    test_connectivity_for_each_owned_element({ {3, {4}} });
  }
  else if (get_parallel_rank() == 1) {
    //  |p0 ||    p1     |
    //  | 3 || 4 | 5 | 6 |
    test_connectivity_for_each_owned_element({ {4, {3}} });
  }
}

TEST_F(ZoltanGraphGeneration, 6Elem2ProcMesh_Decomp_TwoElemsEachProcWithOneAdjacentToProcBoundary)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {3, "partA"}, {4, "partA"}, {6, "partA"}});
  fill_zoltan_graph_for_decomp(*get_meta().get_part("partA"));

  if (get_parallel_rank() == 0) {
    //  |    p0     ||p1 |
    //  | 1 | 2 | 3 || 4 |
    test_connectivity_for_each_owned_element({ {1, {}}, {3, {4}} });
  }
  else if (get_parallel_rank() == 1) {
    //  |p0 ||    p1     |
    //  | 3 || 4 | 5 | 6 |
    test_connectivity_for_each_owned_element({ {4, {3}}, {6, {}} });
  }
}
