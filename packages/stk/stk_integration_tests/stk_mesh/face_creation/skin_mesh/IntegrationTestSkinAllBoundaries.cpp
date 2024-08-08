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

const SideTestUtil::TestCaseData allBoundarySidesTestCases =
{
  //  filename,      max     #side    sideset
  //               #procs,  entities,
  {"AA.e",        2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"AB.e",        2,      11,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}}},
  {"Ae.e",        2,       7,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}}},
  {"Aef.e",       3,       7,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {3, 0}, {3, 1}}},
  {"AeA.e",       3,      12,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}, {3, 0}, {3, 1}}},
  {"AeB.e",       3,      12,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 4}, {3, 5}}},
  {"AefA.e",      4,      12,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}, {3, 0}, {3, 1}, {4, 0}, {4, 1}}},
  {"AefB.e",      4,      12,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {3, 0}, {3, 1}, {4, 0}, {4, 1}, {4, 2}, {4, 3}, {4, 4}, {4, 5}}},
  {"ef.e",        2,       2,     {{1, 0}, {1, 1}, {2, 0}, {2, 1}}},

  {"ALRA_doubleKissing.e", 2, 8,  {{1, 0}, {1, 3}, {1, 4}, {1, 5}, {2, 1}, {2, 2}, {2, 4}, {2, 5}}},
  {"AB_doubleKissing.e",  2, 10,  {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}}},

  {"Tg.e",        2,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}}},
  {"ZY.e",        2,      11,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}}},
  {"AP.e",        2,      11,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}}},
};

void create_all_boundary_sides(stk::mesh::BulkData& bulkData, stk::mesh::Selector stuffToSkin, const stk::mesh::PartVector& partsToAdd)
{
}

class AllBoundarySidesTester : public SideTestUtil::SideCreationTester
{
public:
  AllBoundarySidesTester() : SideTestUtil::SideCreationTester(MPI_COMM_WORLD) {}
protected:
  virtual void test_side_creation(stk::mesh::BulkData& bulkData,
                                  const SideTestUtil::TestCase& testCase)
  {
    stk::mesh::Part& skinnedPart = bulkData.mesh_meta_data().declare_part("interior");
    create_all_boundary_sides(bulkData, bulkData.mesh_meta_data().universal_part(), stk::mesh::PartVector{&skinnedPart});
    expect_all_boundary_sides_connected_as_specified_in_test_case(bulkData, testCase, bulkData.mesh_meta_data().universal_part(), skinnedPart);
  }

};

//Disabled because we haven't had a story to implement create_all_boundary_sides yet.  Trivial implementation
//of calling create_exposed then create_interior did not work.
TEST(AllBoundarySidesTest, DISABLED_run_all_test_cases_aura)
{
  AllBoundarySidesTester().run_all_test_cases(allBoundarySidesTestCases, stk::mesh::BulkData::AUTO_AURA);
}

TEST(AllBoundarySidesTest, DISABLED_run_all_test_cases_no_aura)
{
  AllBoundarySidesTester().run_all_test_cases(allBoundarySidesTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}

}
