/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <Ionit_Initializer.h>     // for Initializer
#include <stddef.h>                     // for size_t, nullptr
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <string>                       // for string
#include "mpi.h"                        // for MPI_COMM_WORLD
#include "stk_io/DatabasePurpose.hpp"
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include <stk_unit_test_utils/FaceCreationTestUtils.hpp>

namespace
{

const SideTestUtil::TestCaseData allSidesTestCasesUniversal =
{
  //  filename,      max     #side    sideset
  //               #procs,  entities,
  {"AA.e",        2,      11,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}}},
  {"AB.e",        2,      11,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}}},
  {"Ae.e",        2,       7,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}}},
  {"AeA.e",       3,      12,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}, {3, 0}, {3, 1}}},
  {"AeB.e",       3,      12,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 4}, {3, 5}}},

  {"AB_doubleKissing.e",  2, 10,  {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}}},

  {"Tg.e",        2,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}}},
  {"ZY.e",        2,      11,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}}},
  {"AP.e",        2,      11,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}}},
};

const SideTestUtil::TestCaseData allSidesTestCasesBlock1 =
{
  //  filename,      max     #side    sideset
  //               #procs,  entities,
  {"AA.e",        2,      11,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}}},
  {"AB.e",        2,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 4}}},
  {"Ae.e",        2,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 1}}},
  {"AeA.e",       3,      12,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}, {3, 0}, {3, 1}}},
  {"AeB.e",       3,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 1}}},

  {"AB_doubleKissing.e",  2, 6,  {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 3}}},

  {"Tg.e",        2,       4,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}}},
  {"ZY.e",        2,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 4}}},
  {"AP.e",        2,      11,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}}},
};

const SideTestUtil::TestCaseData failingAllSidesTestCases =
{
  //Appears to be a problem with pre-existing faces in the mesh.
  {"ALRA_doubleKissing.e", 2, 10,  {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}}},
  //np3 fails consistency checks due to face having different node ordering on different procs.  Maybe due to split coincidents
  {"AefA.e",      4,      12,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}, {3, 0}, {3, 1}, {4, 0}, {4, 1}}},
};

const SideTestUtil::TestCaseData allSidesCoincidentElementTestCases =
{
  {"Aef.e",       3,       7,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {3, 0}, {3, 1}}},
  {"AefA.e",      2,      12,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}, {3, 0}, {3, 1}, {4, 0}, {4, 1}}},
  {"AefB.e",      4,      12,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {3, 0}, {3, 1}, {4, 0}, {4, 1}, {4, 2}, {4, 3}, {4, 4}, {4, 5}}},
  {"ef.e",        2,       2,     {{1, 0}, {1, 1}, {2, 0}, {2, 1}}},
};

class AllSidesTester : public SideTestUtil::SideCreationTester, public ::testing::Test
{
public:
  enum WhichSelector { UNIVERSAL, BLOCK_1 };
  AllSidesTester() : SideTestUtil::SideCreationTester(MPI_COMM_WORLD) {}
protected:
  virtual void test_side_creation(stk::mesh::BulkData& bulkData,
                                  const SideTestUtil::TestCase& testCase)
  {
    stk::mesh::Part& skinnedPart = bulkData.mesh_meta_data().declare_part("interior");
    stk::mesh::create_all_sides(bulkData, get_skin_selector(bulkData.mesh_meta_data()), stk::mesh::PartVector{&skinnedPart}, false);
    expect_all_sides_connected_as_specified_in_test_case(bulkData, testCase, bulkData.mesh_meta_data().universal_part(), skinnedPart);
  }
  stk::mesh::Selector get_skin_selector(const stk::mesh::MetaData& metaData)
  {
    if (m_whichSelector == UNIVERSAL)
    {
      return metaData.universal_part();
    }
    else
    {
      return *metaData.get_part("block_1");
    }
  }
  void set_which_selector(WhichSelector whichSelector)
  {
    m_whichSelector = whichSelector;
  }
private:
  WhichSelector m_whichSelector = UNIVERSAL;
};

TEST_F(AllSidesTester, run_all_test_cases_universal_aura)
{
  run_all_test_cases(allSidesTestCasesUniversal, stk::mesh::BulkData::AUTO_AURA);
}

TEST_F(AllSidesTester, run_all_test_cases_universal_no_aura)
{
  run_all_test_cases(allSidesTestCasesUniversal, stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST_F(AllSidesTester, run_all_test_cases_block1_aura)
{
  set_which_selector(BLOCK_1);
  run_all_test_cases(allSidesTestCasesBlock1, stk::mesh::BulkData::AUTO_AURA);
}

TEST_F(AllSidesTester, run_all_test_cases_block1_no_aura)
{
  set_which_selector(BLOCK_1);
  run_all_test_cases(allSidesTestCasesBlock1, stk::mesh::BulkData::NO_AUTO_AURA);
}
// np3 fails consistency checks due to face having different node ordering on different procs.  Maybe due to split coincidents
// (failing before split coincident element support was removed)
TEST_F(AllSidesTester, DISABLED_run_all_failing_test_cases_aura)
{
  run_all_test_cases(failingAllSidesTestCases, stk::mesh::BulkData::AUTO_AURA);
}

TEST_F(AllSidesTester, run_coincident_element_test_cases_no_aura)
{
  run_all_test_cases(allSidesCoincidentElementTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST_F(AllSidesTester, np1_run_coincident_element_test_cases_no_aura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
    run_all_test_cases(allSidesCoincidentElementTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}
}
