/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <Ioss_IOFactory.h>             // for IOFactory
#include <Ioss_Region.h>                // for Region
#include <init/Ionit_Initializer.h>     // for Initializer
#include <stddef.h>                     // for size_t, nullptr
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <string>                       // for string
#include "mpi.h"                        // for MPI_COMM_WORLD
#include "stk_io/DatabasePurpose.hpp"
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include "../FaceCreationTestUtils.hpp"

namespace
{

const SideTestUtil::TestCaseData exposedBoundaryTestCases =
{
  /* filename, max#procs, #side,  sideset */
    {"A.e",         1,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}}},
    {"AA.e",        2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"AB.e",        2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ADA.e",       2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ADB.e",       2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ADDA.e",      2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ADDB.e",      2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ADReA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ADReB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
    {"ADe.e",       2,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}}},
    {"ADeDA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ADeDB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
    {"ADeLA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ADeLB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
    {"ADeRA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ADeRB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
    {"AL.e",        1,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}}},
    {"ALA.e",       2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALB.e",       2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALJ.e",       3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
    {"ALRA.e",      2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALRB.e",      2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALReA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALReB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
    {"ALe.e",       2,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}}},
    {"ALeDA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALeDB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
    {"ALeDfRA.e",   4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALeDfRB.e",   4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALeLA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALeLB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
    {"ALeRA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALeRB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
    {"ALeXfRA.e",   4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALefRA.e",    4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ARA.e",       2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ARB.e",       2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ARReA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ARReB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
    {"ARe.e",       2,      6,      {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}}},
    {"AReDA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"AReDB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
    {"AReLA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"AReLB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
    {"AReRA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"AReRB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
    {"ARefLA.e",    4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"Ae.e",        2,      6,      {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}}},
    {"AeA.e",       3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"AeB.e",       3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
    {"AeDfA.e",     4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"Aef.e",       3,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {3, 0}}},
    {"AefA.e",      4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"AefB.e",      4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {4, 0}, {4, 1}, {4, 2}, {4, 3}, {4, 5}}},
    {"Re.e",        1,      2,      {{1, 0}, {1, 1}}},
    {"ReL.e",       1,      2,      {{1, 0}, {1, 1}}},
    {"e.e",         1,      2,      {{1, 0}, {1, 1}}},
    {"eL.e",        1,      2,      {{1, 0}, {1, 1}}},
    {"ef.e",        2,      2,      {{1, 0}, {1, 1}, {2, 0}, {2, 1}}},

    {"AB_doubleKissing.e",  2,  8,  {{1, 0}, {1, 3}, {1, 4}, {1, 5}, {2, 1}, {2, 2}, {2, 4}, {2, 5}}},
    {"Tg.e",        2,      6,      {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}}},
    {"ZY.e",        2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}}
};

const SideTestUtil::TestCaseData createExposedBoundaryForOneBlockTestCases =
{
  /* filename, max#procs, #side,  sideset */
    {"AB.e",        2,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"ADB.e",       2,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"ADDB.e",      2,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"ADReB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"ADeDB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"ADeLB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"ADeRB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"ALB.e",       2,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"ALRB.e",      2,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"ALReB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"ALeDB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"ALeDfRB.e",   4,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"ALeLB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"ALeRB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"ARB.e",       2,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"ARReB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"AReDB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"AReLB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"AReRB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"AeB.e",       3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
    {"AefB.e",      4,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},

    {"ALReA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALeDA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALeDfRA.e",   4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALeLA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALeRA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALeXfRA.e",   4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ALefRA.e",    4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ARReA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"AReDA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"AReLA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"AReRA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"ARefLA.e",    4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"AeA.e",       3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"AeDfA.e",     4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
    {"AefA.e",      4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},

    {"Ae.e",        2,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},

    // Throws on 2 procs with stk::CommBuffer::unpack<T>(...){ overflow by -12 bytes. }
//    {"ef.e",        2,       2,     {{1, 0}, {1, 1}, {2, 0}, {2, 1}}},

    {"AB_doubleKissing.e",  2,  4,  {{1, 0}, {1, 3}, {1, 4}, {1, 5}}},
    {"Tg.e",        2,      4,      {{1, 0}, {1, 1}, {1, 2}, {1, 3}}},
    {"ZY.e",        2,      5,      {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}}
};

class SkinnedMesh : public SideTestUtil::SideCreationTester
{
public:
    SkinnedMesh() : SideTestUtil::SideCreationTester(MPI_COMM_WORLD) {}
protected:
    virtual void test_side_creation(stk::mesh::BulkData& bulkData,
                                    const SideTestUtil::TestCase& testCase)
    {
        stk::mesh::Part& skinnedPart = SideTestUtil::run_skin_mesh(bulkData, get_things_to_skin(bulkData));
        expect_exposed_sides_connected_as_specified_in_test_case(bulkData, testCase, skinnedPart);
    }

    void expect_exposed_sides_connected_as_specified_in_test_case(stk::mesh::BulkData& bulkData,
                                                                  const SideTestUtil::TestCase& testCase,
                                                                  stk::mesh::Part &skinnedPart)
    {
        SideTestUtil::expect_global_num_sides_in_part(bulkData, testCase, skinnedPart);
        SideTestUtil::expect_all_sides_exist_for_elem_side(bulkData, testCase.filename, testCase.sideSet);
        EXPECT_TRUE(stk::mesh::check_exposed_boundary_sides(bulkData, get_things_to_skin(bulkData), skinnedPart));
    }

    virtual stk::mesh::Selector get_things_to_skin(const stk::mesh::BulkData& bulkData)
    {
        return bulkData.mesh_meta_data().universal_part();
    }
};

class SkinSingleBlock : public SkinnedMesh
{
protected:
    virtual stk::mesh::Selector get_things_to_skin(const stk::mesh::BulkData& bulkData)
    {
        return *bulkData.mesh_meta_data().get_part("block_1");
    }
};

TEST(ExposedBlockBoundaryTest, run_all_test_cases_aura)
{
    SkinnedMesh().run_all_test_cases(exposedBoundaryTestCases, stk::mesh::BulkData::AUTO_AURA);
}

TEST(ExposedBlockBoundaryTest, run_all_test_cases_no_aura)
{
    SkinnedMesh().run_all_test_cases(exposedBoundaryTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST(CreateExposedBoundaryForSingleBlockTest, run_all_test_cases_aura)
{
    SkinSingleBlock().run_all_test_cases(createExposedBoundaryForOneBlockTestCases, stk::mesh::BulkData::AUTO_AURA);
}

TEST(CreateExposedBoundaryForSingleBlockTest, run_all_test_cases_no_aura)
{
    SkinSingleBlock().run_all_test_cases(createExposedBoundaryForOneBlockTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}

} //namespace
