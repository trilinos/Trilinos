#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/Zoltan2ParallelGraph.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
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
    setup_mesh(inputMeshFile, stk::mesh::BulkData::NO_AUTO_AURA);
  }

  void balance_mesh(int numFinalSubdomains, const std::vector<stk::mesh::Selector> & selectors)
  {
    m_numFinalSubdomains = numFinalSubdomains;
    stk::balance::GraphCreationSettings balanceSettings;
    testing::internal::CaptureStdout();
    stk::balance::internal::calculateGeometricOrGraphBasedDecomp(balanceSettings, numFinalSubdomains,
                                                                 m_decomp, get_bulk(), selectors);
    testing::internal::GetCapturedStdout();
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
  balance_mesh(2, {get_meta().universal_part()});

  test_partition_element_distribution({3, 3});
}


TEST_F(StkBalancePartitioning, 6Elem2ProcMesh_EntireDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  balance_mesh(2, {get_meta().universal_part()});

  test_partition_element_distribution({3, 3});
}

TEST_F(StkBalancePartitioning, 6Elem2ProcMesh_OneElem)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}});
  balance_mesh(2, {*get_meta().get_part("partA")});

  test_partition_element_distribution({0, 1});
}

TEST_F(StkBalancePartitioning, 6Elem2ProcMesh_TwoElems)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {2, "partA"}});
  balance_mesh(2, {*get_meta().get_part("partA")});

  test_partition_element_distribution({1, 1});
}

TEST_F(StkBalancePartitioning, 6Elem2ProcMesh_TwoElemsEachProc)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {2, "partA"}, {4, "partA"}, {5, "partA"}});
  balance_mesh(2, {*get_meta().get_part("partA")});

  test_partition_element_distribution({2, 2});
}

TEST_F(StkBalancePartitioning, 6Elem2ProcMesh_TwoElemsAcrossProcBoundary)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{3, "partA"}, {4, "partA"}});
  balance_mesh(2, {*get_meta().get_part("partA")});

  test_partition_element_distribution({1, 1});
}

TEST_F(StkBalancePartitioning, 6Elem2ProcMesh_TwoElemsEachProcWithOneAdjacentToProcBoundary)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  stk::unit_test_util::put_elements_into_part(get_bulk(), {{1, "partA"}, {3, "partA"}, {4, "partA"}, {6, "partA"}});
  balance_mesh(2, {*get_meta().get_part("partA")});

  test_partition_element_distribution({2, 2});
}

TEST_F(StkBalancePartitioning, 6Elem2to3ProcMesh_EntireDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  balance_mesh(3, {get_meta().universal_part()});

  test_partition_element_distribution({2, 2, 2});
}

TEST_F(StkBalancePartitioning, 6Elem2to1ProcMesh_EntireDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh("generated:1x1x6");
  balance_mesh(1, {get_meta().universal_part()});

  test_partition_element_distribution({6});
}

TEST_F(StkBalancePartitioning, 6Elem3to2ProcMesh_EntireDomain)
{
  if (stk::parallel_machine_size(get_comm()) != 3) return;

  setup_initial_mesh("generated:1x1x6");
  balance_mesh(2, {get_meta().universal_part()});

  test_partition_element_distribution({3, 3});
}

