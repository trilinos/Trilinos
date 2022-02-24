#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/Zoltan2ParallelGraph.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/EnvData.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <vector>
#include <string>

class StkBalancePartitioning : public stk::unit_test_util::MeshFixture
{
protected:
  StkBalancePartitioning()
    : MeshFixture(),
      m_numFinalSubdomains(0)
  { }

  ~StkBalancePartitioning() = default;

  void setup_initial_mesh(const std::string & inputMeshFile)
  {
    setup_mesh(inputMeshFile, stk::mesh::BulkData::AUTO_AURA);
  }

  void setup_4hex_contact_perpendicular_to_proc_boundary()
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);

    const std::string meshDesc = "0,1,HEX_8,1,2,5,4,13,14,17,16\n"
                                 "1,2,HEX_8,2,3,6,5,14,15,18,17\n"
                                 "0,3,HEX_8,7,8,11,10,19,20,23,22\n"
                                 "1,4,HEX_8,8,9,12,11,20,21,24,23";

    std::vector<double> coordinates = {
      0,0,0, 1,0,0, 2,0,0, 0,1,0, 1,1,0, 2,1,0,
      0,1,0, 1,1,0, 2,1,0, 0,2,0, 1,2,0, 2,2,0,
      0,0,1, 1,0,1, 2,0,1, 0,1,1, 1,1,1, 2,1,1,
      0,1,1, 1,1,1, 2,1,1, 0,2,1, 1,2,1, 2,2,1
    };

    stk::unit_test_util::setup_text_mesh(
        get_bulk(), stk::unit_test_util::get_full_text_mesh_desc(meshDesc, coordinates));
  }

  void balance_mesh(const stk::ParallelMachine & decompCommunicator,
                    int numFinalSubdomains,
                    const std::vector<stk::mesh::Selector> & selectors,
                    stk::balance::BalanceSettings & balanceSettings)
  {
    m_numFinalSubdomains = numFinalSubdomains;
    stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;
    stk::balance::internal::calculateGeometricOrGraphBasedDecomp(get_bulk(), selectors,
                                                                 decompCommunicator, numFinalSubdomains,
                                                                 balanceSettings, m_decomp);
    stk::EnvData::instance().m_outputP0 = &std::cout;

  }

  void balance_mesh(const stk::ParallelMachine & decompCommunicator,
                    int numFinalSubdomains,
                    const std::vector<stk::mesh::Selector> & selectors)
  {
    stk::balance::GraphCreationSettings balanceSettings;
    balance_mesh(decompCommunicator, numFinalSubdomains, selectors, balanceSettings);
  }

  void balance_mesh_scotch(const stk::ParallelMachine & decompCommunicator,
                           int numFinalSubdomains,
                           const std::vector<stk::mesh::Selector> & selectors)
  {
    stk::balance::GraphCreationSettings balanceSettings;
    balanceSettings.setDecompMethod("scotch");
    balance_mesh(decompCommunicator, numFinalSubdomains, selectors, balanceSettings);
  }

  void test_partition_element_distribution(const std::vector<int> & expectedElemsPerProc)
  {
    ASSERT_EQ(static_cast<unsigned>(m_numFinalSubdomains), expectedElemsPerProc.size());

    for (int subdomain = 0; subdomain < m_numFinalSubdomains; ++subdomain) {
      int numLocalElementsForSubdomain = 0;
      for (const stk::mesh::EntityProc & entityProc : m_decomp) {
        if (get_bulk().is_valid(entityProc.first) && entityProc.second == subdomain) {
          ++numLocalElementsForSubdomain;
        }
      }

      const int numGlobalElementsForSubdomain = stk::get_global_sum(get_comm(), numLocalElementsForSubdomain);
      if (get_parallel_rank() == 0) {
        EXPECT_EQ(numGlobalElementsForSubdomain, expectedElemsPerProc[subdomain]);
      }
    }
  }

  int m_numFinalSubdomains;
  stk::mesh::EntityProcVec m_decomp;
};


TEST_F(StkBalancePartitioning, 6Elem1ProcMesh_EntireDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  setup_initial_mesh("generated:1x1x6");
  balance_mesh(get_bulk().parallel(), 2, {get_meta().universal_part()});

  test_partition_element_distribution({3, 3});
}


TEST_F(StkBalancePartitioning, 6Elem2ProcMesh_EntireDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  balance_mesh(get_bulk().parallel(), 2, {get_meta().universal_part()});

  test_partition_element_distribution({3, 3});
}

TEST_F(StkBalancePartitioning, 6Elem2ProcMesh_OneElem)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}});
  balance_mesh(get_bulk().parallel(), 2, {*get_meta().get_part("partA")});

  test_partition_element_distribution({0, 1});
}

TEST_F(StkBalancePartitioning, 6Elem2ProcMesh_TwoElems)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {2, "partA"}});
  balance_mesh(get_bulk().parallel(), 2, {*get_meta().get_part("partA")});

  test_partition_element_distribution({1, 1});
}

TEST_F(StkBalancePartitioning, 6Elem2ProcMesh_TwoElemsEachProc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {2, "partA"}, {4, "partA"}, {5, "partA"}});
  balance_mesh(get_bulk().parallel(), 2, {*get_meta().get_part("partA")});

  test_partition_element_distribution({2, 2});
}

TEST_F(StkBalancePartitioning, 6Elem2ProcMesh_TwoElemsAcrossProcBoundary)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{3, "partA"}, {4, "partA"}});
  balance_mesh(get_bulk().parallel(), 2, {*get_meta().get_part("partA")});

  test_partition_element_distribution({1, 1});
}

TEST_F(StkBalancePartitioning, 6Elem2ProcMesh_TwoElemsEachProcWithOneAdjacentToProcBoundary)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {3, "partA"}, {4, "partA"}, {6, "partA"}});
  balance_mesh(get_bulk().parallel(), 2, {*get_meta().get_part("partA")});

  test_partition_element_distribution({2, 2});
}

TEST_F(StkBalancePartitioning, 6Elem2to3ProcMesh_EntireDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  balance_mesh(get_bulk().parallel(), 3, {get_meta().universal_part()});

  test_partition_element_distribution({2, 2, 2});
}

TEST_F(StkBalancePartitioning, 6Elem1to1ProcMesh_EntireDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 1) return;

  setup_initial_mesh("generated:1x1x6");
  balance_mesh(get_bulk().parallel(), 1, {get_meta().universal_part()});

  test_partition_element_distribution({6});
}

TEST_F(StkBalancePartitioning, 6Elem2to1ProcMesh_EntireDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  balance_mesh(get_bulk().parallel(), 1, {get_meta().universal_part()});

  test_partition_element_distribution({6});
}

TEST_F(StkBalancePartitioning, 6Elem3to1ProcMesh_EntireDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 3) return;

  setup_initial_mesh("generated:1x1x6");
  balance_mesh(get_bulk().parallel(), 1, {get_meta().universal_part()});

  test_partition_element_distribution({6});
}

TEST_F(StkBalancePartitioning, 6Elem2to1ProcMesh_HalfDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {2, "partA"}, {3, "partA"}});
  balance_mesh(get_bulk().parallel(), 1, {*get_meta().get_part("partA")});

  test_partition_element_distribution({3});
}

TEST_F(StkBalancePartitioning, 6Elem3to2ProcMesh_EntireDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 3) return;

  setup_initial_mesh("generated:1x1x6");
  balance_mesh(get_bulk().parallel(), 2, {get_meta().universal_part()});

  test_partition_element_distribution({3, 3});
}

TEST_F(StkBalancePartitioning, 4Elem2ProcMeshWithContact_EmptyOnOneProc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_4hex_contact_perpendicular_to_proc_boundary();
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {3, "partA"}});
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{2, "partB"}, {4, "partB"}});
  balance_mesh(get_bulk().parallel(), 2, {*get_meta().get_part("partB")});

  test_partition_element_distribution({1, 1});
}

TEST_F(StkBalancePartitioning, 4Elem2ProcMesh_SeparateCommunicator)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x4");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {2, "partA"}, {3, "partB"}, {4, "partB"}});
  if (get_parallel_rank() == 0) {
    balance_mesh(MPI_COMM_SELF, 2, {*get_meta().get_part("partA")});
  }
  else {
    balance_mesh(MPI_COMM_SELF, 2, {*get_meta().get_part("partB")});
  }

  test_partition_element_distribution({2, 2});
}

TEST_F(StkBalancePartitioning, 4Elem2to1ProcMesh_SeparateCommunicator)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x4");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {2, "partA"}, {3, "partB"}, {4, "partB"}});
  if (get_parallel_rank() == 0) {
    balance_mesh(MPI_COMM_SELF, 1, {*get_meta().get_part("partA")});
  }
  else {
    balance_mesh(MPI_COMM_SELF, 1, {*get_meta().get_part("partB")});
  }

  test_partition_element_distribution({4});
}

TEST_F(StkBalancePartitioning, 4Elem2ProcMesh_SeparateCommunicator_EmptyOnOneProc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x4");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {2, "partA"}});
  balance_mesh(MPI_COMM_SELF, 2, {*get_meta().get_part("partA")});

  test_partition_element_distribution({1, 1});
}

TEST_F(StkBalancePartitioning, 4Elem2ProcMesh_Geometric_SeparateCommunicator_EmptyOnOneProc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x4");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {2, "partA"}});
  stk::balance::BasicGeometricSettings balanceSettings;
  balance_mesh(MPI_COMM_SELF, 2, {*get_meta().get_part("partA")}, balanceSettings);

  test_partition_element_distribution({1, 1});
}

TEST_F(StkBalancePartitioning, 6Elem1ProcMesh_EntireDomain_Scotch)
{
  if (stk::parallel_machine_size(get_comm()) != 1) GTEST_SKIP();

  setup_initial_mesh("generated:1x1x6");
  balance_mesh_scotch(get_bulk().parallel(), 2, {get_meta().universal_part()});

  test_partition_element_distribution({3, 3});
}


TEST_F(StkBalancePartitioning, 6Elem2ProcMesh_EntireDomain_Scotch)
{
  if (stk::parallel_machine_size(get_comm()) != 2) GTEST_SKIP();

  setup_initial_mesh("generated:1x1x6");
  balance_mesh_scotch(get_bulk().parallel(), 2, {get_meta().universal_part()});

  test_partition_element_distribution({3, 3});
}

