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
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include <stk_unit_test_utils/FaceCreationTestUtils.hpp>

namespace
{

const SideTestUtil::TestCaseData interiorBlockBoundaryTestCases =
{
  /* filename, max#procs, #side,  sideset */
  {"AB.e",      2,        1,    {{1, 5}, {2, 4}}},
  {"Ae.e",      2,        1,    {{1, 5}, {2, 1}}},
  {"AeB.e",     3,        2,    {{1, 5}, {2, 0}, {2, 1}, {3, 4}}},

  {"AB_doubleKissing.e", 2, 2,  {{1, 1}, {1, 2}, {2, 0}, {2, 3}}},

  {"Tg.e",      2,        0,    {}},
  {"ZY.e",      2,        1,    {{1, 5}, {2, 4}}},

  {"2D_AB.e",     2,       1,     {{1, 0}, {2, 2}}},
  {"2D_AeB.e",    3,       2,     {{1, 0}, {2, 2}, {3, 0}, {3, 1}}},
  {"2D_AtB.e",    3,       1,     {{1, 0}, {2, 2}, {3, 0}}},
  {"2D_ArB.e",    3,       1,     {{1, 0}, {2, 2}, {3, 0}}},
  {"2D_AmB.e",    3,       1,     {{1, 0}, {2, 2}, {3, 0}}},
};

const SideTestUtil::TestCaseData failingInteriorBlockBoundaryTestCases =
{
  //np3 fails consistency checks due to face having different node ordering on different procs.  Maybe due to split coincidents
  {"AefA.e",    4,        2,    {{1, 5}, {3, 0}, {3, 1}, {4, 0}, {4, 1}, {2, 4}}},

  //np2 fails because of split (partial) coincident elements even after attempt to fix split coincidents
  {"2D_1block_full_refine.e",     2,       0,     {}},
  {"2D_2block_full_refine.e",     2,       3,     {{1, 2}, {2, 1}, {3, 1}, {4, 1}, {5, 2}, {6, 0}}},
  {"2D_1block_full_decomp.e",     2,       2,     {{3, 2}, {4, 0}, {5, 1}, {6, 1}}},
  {"2D_2block_full_decomp.e",     2,       5,     {{1, 2}, {2, 1}, {3, 1}, {3, 2}, {4, 0}, {4, 1}, {5, 1}, {5, 2}, {6, 0}, {6, 1}}},
  {"2D_1block_partial_refine.e",  2,       1,     {{1, 2}, {2, 1}}},
  {"2D_2block_partial_refine.e",  2,       1,     {{1, 2}, {2, 1}, {4, 1}}},
  {"2D_1block_partial_decomp.e",  2,       2,     {{1, 2}, {2, 1}, {3, 2}, {4, 0}}},
  {"2D_2block_partial_decomp.e",  2,       2,     {{1, 2}, {2, 1}, {3, 2}, {4, 0}, {4, 1}}},
};

const SideTestUtil::TestCaseData interiorBlockBoundaryCoincidentElementTestCases =
{
  {"Aef.e",     3,        1,    {{1, 5}, {2, 1}, {3, 1}}},
  {"AefB.e",    4,        2,    {{1, 5}, {2, 0}, {2, 1}, {3, 0}, {3, 1}, {4, 4}}},
  {"ef.e",      2,        0,    {}},
  {"eff.e",     3,        0,    {}},

  {"AP.e",      2,        0,    {}},

  {"2D_1block_full_refine.e",     1,       0,     {}},
  {"2D_2block_full_refine.e",     1,       3,     {{1, 2}, {2, 1}, {3, 1}, {4, 1}, {5, 2}, {6, 0}}},
  {"2D_1block_full_decomp.e",     1,       2,     {{3, 2}, {4, 0}, {5, 1}, {6, 1}}},
  {"2D_2block_full_decomp.e",     1,       5,     {{1, 2}, {2, 1}, {3, 1}, {3, 2}, {4, 0}, {4, 1}, {5, 1}, {5, 2}, {6, 0}, {6, 1}}},
  {"2D_1block_partial_refine.e",  1,       1,     {{1, 2}, {2, 1}}},
  {"2D_2block_partial_refine.e",  1,       1,     {{1, 2}, {2, 1}, {4, 1}}},
  {"2D_1block_partial_decomp.e",  1,       2,     {{1, 2}, {2, 1}, {3, 2}, {4, 0}}},
  {"2D_2block_partial_decomp.e",  1,       2,     {{1, 2}, {2, 1}, {3, 2}, {4, 0}, {4, 1}}},
};

const SideTestUtil::TestCaseData createInteriorBoundaryForOneBlockTestCases =
{
  /* filename, max#procs, #side,  sideset */
  {"AB.e",      2,        0,    {}},
  {"Ae.e",      2,        0,    {}},
  {"AeB.e",     3,        0,    {}},

  {"AB_doubleKissing.e", 2, 0,  {}},

  {"Tg.e",      2,        0,    {}},
  {"ZY.e",      2,        0,    {}}
};

const SideTestUtil::TestCaseData createInteriorBoundaryForTwoBlockTestCases =
{
  /* filename, max#procs, #side,  sideset */
  {"AB.e",      2,        1,    {{1, 5}, {2, 4}}},
  {"Ae.e",      2,        1,    {{1, 5}, {2, 1}}},
  {"AeB.e",     3,        1,    {{1, 5}, {2, 1}}},

  {"AB_doubleKissing.e", 2, 2,  {{1, 1}, {1, 2}, {2, 0}, {2, 3}}},

  {"Tg.e",      2,        0,    {}},
  {"ZY.e",      2,        1,    {{1, 5}, {2, 4}}}
};

const SideTestUtil::TestCaseData createInteriorBoundaryForOneBlockCoincidentElementTestCases =
{
  {"Aef.e",     3,        0,    {}},
  {"AefB.e",    4,        0,    {}},
  {"ef.e",      2,        0,    {}},
  {"eff.e",     3,        0,    {}},
};

const SideTestUtil::TestCaseData createInteriorBoundaryForTwoBlockCoincidentElementTestCases =
{
  {"Aef.e",     3,        1,    {{1, 5}, {2, 1}, {3, 1}}}, // Is this {3, 1} entry correct?
  {"AefB.e",    4,        1,    {{1, 5}, {2, 1}, {3, 1}}}, // Is this {3, 1} entry correct?
  {"ef.e",      2,        0,    {}},
  {"eff.e",     3,        0,    {}},

  {"2D_1block_full_refine.e",     1,       0,     {}},
  {"2D_2block_full_refine.e",     1,       2,     {{3, 1}, {4, 1}, {5, 2}, {6, 0}}},
  {"2D_1block_full_decomp.e",     1,       2,     {{3, 2}, {4, 0}, {5, 1}, {6, 1}}},
  {"2D_2block_full_decomp.e",     1,       1,     {{3, 2}, {4, 0}}},
  {"2D_1block_partial_refine.e",  1,       0,     {}},
  {"2D_2block_partial_refine.e",  1,       1,     {{1, 2}, {2, 1}, {4, 1}}},
  {"2D_1block_partial_decomp.e",  1,       1,     {{3, 2}, {4, 0}}},
  {"2D_2block_partial_decomp.e",  1,       0,     {}}
};

const SideTestUtil::TestCaseData createInteriorBoundaryForDegenerateElementTestCases =
{
  {"quadInteriorDegenerateConnection.g",  3,  1, {}},
  {"hexInteriorDegeneratedConnectionA.g", 3,  2, {}},
  {"hexInteriorDegeneratedConnectionB.g", 3,  1, {}}
};

class InteriorBlockBoundaryTester : public SideTestUtil::SideCreationTester
{
public:
  InteriorBlockBoundaryTester() : SideTestUtil::SideCreationTester(MPI_COMM_WORLD) {}
protected:
  virtual void test_side_creation(stk::mesh::BulkData& bulkData,
                                  const SideTestUtil::TestCase& testCase)
  {
    stk::mesh::Part &skinnedPart = bulkData.mesh_meta_data().declare_part("interior", bulkData.mesh_meta_data().side_rank());
    stk::mesh::create_interior_block_boundary_sides(bulkData, get_things_to_skin(bulkData), stk::mesh::PartVector{&skinnedPart});
    SideTestUtil::expect_interior_sides_connected_as_specified_in_test_case(bulkData, testCase, get_things_to_skin(bulkData), skinnedPart);
  }

  virtual stk::mesh::Selector get_things_to_skin(const stk::mesh::BulkData& bulkData)
  {
    return bulkData.mesh_meta_data().universal_part();
  }
};

class OneBlockInteriorBlockBoundaryTester : public InteriorBlockBoundaryTester
{
protected:
  virtual stk::mesh::Selector get_things_to_skin(const stk::mesh::BulkData& bulkData)
  {
    return *bulkData.mesh_meta_data().get_part("block_1");
  }
};

class TwoBlockInteriorBlockBoundaryTester : public InteriorBlockBoundaryTester
{
protected:
  virtual stk::mesh::Selector get_things_to_skin(const stk::mesh::BulkData& bulkData)
  {
    stk::mesh::Selector selector;
    if (bulkData.mesh_meta_data().get_part("block_1") != nullptr) {
      selector = *bulkData.mesh_meta_data().get_part("block_1");
    }
    if (bulkData.mesh_meta_data().get_part("block_2") != nullptr) {
      selector |= *bulkData.mesh_meta_data().get_part("block_2");
    }
    return selector;
  }
};

TEST(InteriorBlockBoundaryTest, run_all_test_cases_aura)
{
  InteriorBlockBoundaryTester().run_all_test_cases(interiorBlockBoundaryTestCases, stk::mesh::BulkData::AUTO_AURA);
}

TEST(InteriorBlockBoundaryTest, run_all_test_cases_no_aura)
{
  InteriorBlockBoundaryTester().run_all_test_cases(interiorBlockBoundaryTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST(InteriorBlockBoundaryTest, run_interior_block_boundary_degenerate_cases_aura)
{
  InteriorBlockBoundaryTester().run_all_test_cases(createInteriorBoundaryForDegenerateElementTestCases, stk::mesh::BulkData::AUTO_AURA);
}

TEST(InteriorBlockBoundaryTest, run_interior_block_boundary_degenerate_cases_no_aura)
{
  InteriorBlockBoundaryTester().run_all_test_cases(createInteriorBoundaryForDegenerateElementTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST(CreateInteriorBoundaryForSingleBlockTest, run_all_test_cases_aura)
{
  OneBlockInteriorBlockBoundaryTester().run_all_test_cases(createInteriorBoundaryForOneBlockTestCases, stk::mesh::BulkData::AUTO_AURA);
}

TEST(CreateInteriorBoundaryForSingleBlockTest, run_all_test_cases_no_aura)
{
  OneBlockInteriorBlockBoundaryTester().run_all_test_cases(createInteriorBoundaryForOneBlockTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST(CreateInteriorBoundaryForTwoBlockTest, run_all_test_cases_aura)
{
  TwoBlockInteriorBlockBoundaryTester().run_all_test_cases(createInteriorBoundaryForTwoBlockTestCases, stk::mesh::BulkData::AUTO_AURA);
}

TEST(CreateInteriorBoundaryForTwoBlockTest, run_all_test_cases_no_aura)
{
  TwoBlockInteriorBlockBoundaryTester().run_all_test_cases(createInteriorBoundaryForTwoBlockTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}

//np3 fails consistency checks due to face having different node ordering on different procs.  Maybe due to split coincidents
// (failing before split coincident element support was removed)
TEST(InteriorBlockBoundaryTest, DISABLED_failing_run_all_test_cases_aura)
{
  InteriorBlockBoundaryTester().run_all_test_cases(failingInteriorBlockBoundaryTestCases, stk::mesh::BulkData::AUTO_AURA);
}
TEST(InteriorBlockBoundaryTest, run_coincident_element_test_cases_no_aura)
{
  InteriorBlockBoundaryTester().run_all_test_cases(interiorBlockBoundaryCoincidentElementTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST(CreateInteriorBoundaryForSingleBlockTest, run_coincident_element_test_cases_no_aura)
{
  OneBlockInteriorBlockBoundaryTester().run_all_test_cases(createInteriorBoundaryForOneBlockCoincidentElementTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST(CreateInteriorBoundaryForTwoBlockTest, run_coincident_element_test_cases_no_aura)
{
  TwoBlockInteriorBlockBoundaryTester().run_all_test_cases(createInteriorBoundaryForTwoBlockCoincidentElementTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST(InteriorBlockBoundaryTest, np1_run_coincident_element_test_cases_no_aura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
    InteriorBlockBoundaryTester().run_all_test_cases(interiorBlockBoundaryCoincidentElementTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST(CreateInteriorBoundaryForSingleBlockTest, np1_run_coincident_element_test_cases_no_aura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
    OneBlockInteriorBlockBoundaryTester().run_all_test_cases(createInteriorBoundaryForOneBlockCoincidentElementTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST(CreateInteriorBoundaryForTwoBlockTest, np1_run_coincident_element_test_cases_no_aura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
    TwoBlockInteriorBlockBoundaryTester().run_all_test_cases(createInteriorBoundaryForTwoBlockCoincidentElementTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}
}
