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

const SideTestUtil::TestCaseData userCreatedFaceTestCases = {
    /* filename, max#procs, #side,   sideset */
    {"AL.e",      1,        1,    {{1, 5}}},

    {"Re.e",      1,        1,    {{1, 1}}},

    {"eL.e",      1,        1,    {{1, 0}}},

    {"ReL.e",     1,        2,    {{1, 1}, {1, 0}}},

    {"ALA.e",     2,        1,    {{1, 5}, {2, 4}}},
    {"ARA.e",     2,        1,    {{1, 5}, {2, 4}}},
    {"ADA.e",     2,        1,    {{1, 5}, {2, 4}}},
    {"ADB.e",     2,        1,    {{1, 5}, {2, 4}}},
    {"ADDA.e",    2,        1,    {{1, 5}, {2, 4}}},
    {"ADDB.e",    2,        1,    {{1, 5}, {2, 4}}},
    {"ALB.e",     2,        1,    {{1, 5}, {2, 4}}},
    {"ALRA.e",    2,        1,    {{1, 5}, {2, 4}}},
    {"ALRB.e",    2,        1,    {{1, 5}, {2, 4}}},
    {"ARB.e",     2,        1,    {{1, 5}, {2, 4}}},

    {"ADe.e",     2,        1,    {{1, 5}, {2, 1}}},
    {"ALe.e",     2,        1,    {{1, 5}, {2, 1}}},
    {"ARe.e",     2,        1,    {{1, 5}, {2, 1}}},

    {"ALJ.e",     3,        1,    {{1, 5}, {2, 4}, {3, 4}}},

    //A=1, e=3, A=2 (global ids)
    {"ADReA.e",   3,        1,    {{1, 5}, {3, 1}}},
    {"ALReA.e",   3,        1,    {{1, 5}, {3, 1}}},
    {"ARReA.e",   3,        1,    {{1, 5}, {3, 1}}},

    //A=1, e=2, B=3 (global ids)
    {"ADReB.e",   3,        1,    {{1, 5}, {2, 1}}},
    {"ALReB.e",   3,        1,    {{1, 5}, {2, 1}}},
    {"ARReB.e",   3,        1,    {{1, 5}, {2, 1}}},

    //A=1, e=3, A=2 (global ids)
    {"ALeRA.e",   3,        2,    {{1, 5}, {3, 0}, {3, 1}, {2, 4}}},
    {"ADeDA.e",   3,        2,    {{1, 5}, {3, 0}, {3, 1}, {2, 4}}},
    {"ADeLA.e",   3,        2,    {{1, 5}, {3, 0}, {3, 1}, {2, 4}}},
    {"ADeRA.e",   3,        2,    {{1, 5}, {3, 0}, {3, 1}, {2, 4}}},
    {"ALeDA.e",   3,        2,    {{1, 5}, {3, 0}, {3, 1}, {2, 4}}},
    {"ALeLA.e",   3,        2,    {{1, 5}, {3, 0}, {3, 1}, {2, 4}}},
    {"ALeRA.e",   3,        2,    {{1, 5}, {3, 0}, {3, 1}, {2, 4}}},
    {"AReDA.e",   3,        2,    {{1, 5}, {3, 0}, {3, 1}, {2, 4}}},
    {"AReLA.e",   3,        2,    {{1, 5}, {3, 0}, {3, 1}, {2, 4}}},
    {"AReRA.e",   3,        2,    {{1, 5}, {3, 0}, {3, 1}, {2, 4}}},

    //A=1, e=2, B=3 (global ids)
    {"ADeDB.e",   3,        2,    {{1, 5}, {2, 0}, {2, 1}, {3, 4}}},
    {"ADeLB.e",   3,        2,    {{1, 5}, {2, 0}, {2, 1}, {3, 4}}},
    {"ADeRB.e",   3,        2,    {{1, 5}, {2, 0}, {2, 1}, {3, 4}}},
    {"ALeDB.e",   3,        2,    {{1, 5}, {2, 0}, {2, 1}, {3, 4}}},
    {"ALeLB.e",   3,        2,    {{1, 5}, {2, 0}, {2, 1}, {3, 4}}},
    {"ALeRB.e",   3,        2,    {{1, 5}, {2, 0}, {2, 1}, {3, 4}}},
    {"AReDB.e",   3,        2,    {{1, 5}, {2, 0}, {2, 1}, {3, 4}}},
    {"AReLB.e",   3,        2,    {{1, 5}, {2, 0}, {2, 1}, {3, 4}}},
    {"AReRB.e",   3,        2,    {{1, 5}, {2, 0}, {2, 1}, {3, 4}}},

    //A=1, e=3, f=4, A/B=2 (global ids)
    {"ALefRA.e",  4,        2,    {{1, 5}, {3, 0}, {3, 1}, {4, 0}, {4, 1}, {2, 4}}},
    {"ARefLA.e",  4,        2,    {{1, 5}, {3, 0}, {3, 1}, {4, 0}, {4, 1}, {2, 4}}},
    {"ALeDfRA.e", 4,        2,    {{1, 5}, {3, 0}, {3, 1}, {4, 0}, {4, 1}, {2, 4}}},
    {"ALeXfRA.e", 4,        2,    {{1, 5}, {3, 0}, {3, 1}, {4, 0}, {4, 1}, {2, 4}}},
    {"ALeDfRB.e", 4,        2,    {{1, 5}, {3, 0}, {3, 1}, {4, 0}, {4, 1}, {2, 4}}},
    {"AeDfA.e",   4,        2,    {{1, 5}, {3, 0}, {3, 1}, {4, 0}, {4, 1}, {2, 4}}},

    {"ADDA_doubleKissing.e", 2, 2, {{1, 1}, {1, 2}, {2, 0}, {2, 3}}},
    {"ARRA_doubleKissing.e", 2, 2, {{1, 1}, {1, 2}, {2, 0}, {2, 3}}},
    {"ALLA_doubleKissing.e", 2, 2, {{1, 1}, {1, 2}, {2, 0}, {2, 3}}},
    {"ALRA_doubleKissing.e", 2, 2, {{1, 1}, {1, 2}, {2, 0}, {2, 3}}},

    // {"basic.e", 4, 2, {{3,3}, {4,3}, {5,1}, {6,1}}} // Ticket 13009 - Disabled. Get this working. 2D example.

    {"TDg.e", 2, 2, {{1, 1}, {2, 1}}}, // Tet adjacent to degenerate quad
    {"ZDZ.e", 2, 1, {{1, 5}, {2, 4}}}, // degenerate Hex adjacent to degenerate Hex
};

class UserCreatedSidesTester : public SideTestUtil::SideCreationTester
{
public:
    UserCreatedSidesTester() : SideTestUtil::SideCreationTester(MPI_COMM_WORLD) {}
protected:
    virtual void test_side_creation(stk::mesh::BulkData& bulkData,
                                    const SideTestUtil::TestCase& testCase)
    {
        SideTestUtil::expect_global_num_sides_correct(bulkData, testCase);
        SideTestUtil::expect_all_sides_exist_for_elem_side(bulkData, testCase.filename, testCase.sideSet);
    }
};

TEST(UserCreatedFaces, read_all_files_aura)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
        UserCreatedSidesTester().run_all_test_cases(userCreatedFaceTestCases, stk::mesh::BulkData::AUTO_AURA);
}

TEST(UserCreatedFaces, read_all_files_no_aura)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
        UserCreatedSidesTester().run_all_test_cases(userCreatedFaceTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST(UserCreatedFaces, DISABLED_parallel_read_all_files_aura)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1)
        UserCreatedSidesTester().run_all_test_cases(userCreatedFaceTestCases, stk::mesh::BulkData::AUTO_AURA);
}

TEST(UserCreatedFaces, DISABLED_parallel_read_all_files_no_aura)
{
    if(stk::parallel_machine_size(MPI_COMM_WORLD) > 1)
        UserCreatedSidesTester().run_all_test_cases(userCreatedFaceTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}

} //namespace
