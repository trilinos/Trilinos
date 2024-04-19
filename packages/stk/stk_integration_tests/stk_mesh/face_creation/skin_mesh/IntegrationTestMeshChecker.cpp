/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t, nullptr
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include <stk_unit_test_utils/FaceCreationTestUtils.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include "mpi.h"
#include <map>
#include <string>
#include <stk_mesh/baseImpl/EquivalentEntityBlocks.hpp>
#include <stk_mesh/base/MeshDiagnostics.hpp>

namespace
{
using stk::unit_test_util::build_mesh;

struct SplitElementInfo
{
  stk::mesh::EntityId localElementId;
  stk::mesh::EntityId remoteElementId;
  int neighboringProc;
};

struct TestCase
{
  std::string filename;
  int maxNumProcs;
  std::vector<std::vector<SplitElementInfo>> expectedSplitElementsPerProc;
};

typedef std::vector<TestCase> TestCaseData;

const TestCaseData badDecomps =
{
  /* filename, #procs,      {local element id, remoted element id, remote proc} */
  {"Aef.e",     2,        { {{3u, 2u, 1}}, {{2u, 3u, 0}} }},
  {"ef.e",      2,        { {{1u, 2u, 1}}, {{2u, 1u, 0}} }},
  {"AefB.e",    3,        { {}, {{2u, 3u, 2}}, {{3u, 2u, 1}} }},
  {"AP.e",      2,        { {{1u, 2u, 1}}, {{2u, 1u, 0}} }}
};

void expect_split_coincidents(const stk::mesh::BulkData &bulkData,
                              const std::vector<std::vector<SplitElementInfo>> &expectedSplitElementsPerProc,
                              const stk::mesh::SplitCoincidentInfo &splitCoincidentElements)
{
  const std::vector<SplitElementInfo> &expectedSplits = expectedSplitElementsPerProc[bulkData.parallel_rank()];
  ASSERT_EQ(expectedSplits.size(), splitCoincidentElements.size());
  auto foundSplitCoincident = splitCoincidentElements.begin();
  for(size_t i = 0; i < expectedSplits.size(); i++)
  {
    EXPECT_EQ(expectedSplits[i].localElementId,  foundSplitCoincident->first);
    EXPECT_EQ(expectedSplits[i].remoteElementId, foundSplitCoincident->second[0].first);
    EXPECT_EQ(expectedSplits[i].neighboringProc, foundSplitCoincident->second[0].second);
    foundSplitCoincident++;
  }
}

class MeshChecker : public ::testing::Test
{
public:
  MeshChecker()
  {
    comm = MPI_COMM_WORLD;
  }

  void run_all_test_cases(const TestCaseData &testCases, stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    for(const TestCase& testCase : testCases)
    {
      if(stk::parallel_machine_size(get_comm()) == testCase.maxNumProcs)
      {
        EXPECT_THROW(test_one_case(testCase, auraOption), std::logic_error);
      }
    }
  }

  void test_one_case(const TestCase &testCase, stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(get_comm(), auraOption);
    stk::mesh::BulkData& bulkData = *bulkPtr;
    SideTestUtil::read_and_decompose_mesh(testCase.filename, bulkData);

    stk::mesh::SplitCoincidentInfo splitCoincidentElements = stk::mesh::get_split_coincident_elements(bulkData);

    expect_split_coincidents(bulkData, testCase.expectedSplitElementsPerProc, splitCoincidentElements);
    stk::mesh::throw_if_any_proc_has_false(bulkData.parallel(), splitCoincidentElements.empty());
  }

  MPI_Comm get_comm() const
  {
    return comm;
  }

private:
  MPI_Comm comm;
};


TEST_F(MeshChecker, diagnose_bad_meshes)
{
  run_all_test_cases(badDecomps, stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST(MeshCheckerIncremental, createSplitCoincident)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  if(stk::parallel_machine_size(comm) == 2)
  {
    std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(comm, stk::mesh::BulkData::NO_AUTO_AURA);
    stk::mesh::BulkData& bulkData = *bulkPtr;
    stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
    SideTestUtil::read_and_decompose_mesh("Ae.e", bulkData);

    stk::mesh::SplitCoincidentInfo splitCoincidentElements = stk::mesh::get_split_coincident_elements(bulkData);
    ASSERT_EQ(0u, splitCoincidentElements.size());

    bulkData.modification_begin();
    if(bulkData.parallel_rank() == 0)
    {
      stk::mesh::Part *block2 = metaData.get_part("block_2");
      stk::mesh::declare_element(bulkData, *block2, 3, {5, 6, 7, 8});
    }
    bulkData.modification_end();

    splitCoincidentElements = stk::mesh::get_split_coincident_elements(bulkData);

    std::vector<std::vector<SplitElementInfo>> expectedSplitElementsPerProc = { {{3u, 2u, 1}}, {{2u, 3u, 0}} };
    expect_split_coincidents(bulkData, expectedSplitElementsPerProc, splitCoincidentElements);
  }
}

}
