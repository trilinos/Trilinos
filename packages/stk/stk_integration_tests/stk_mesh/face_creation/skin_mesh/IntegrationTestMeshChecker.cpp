/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t, nullptr
#include <string>                       // for string
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include "../FaceCreationTestUtils.hpp"
#include <mpi.h>
#include <map>
#include <string>
#include <stk_mesh/baseImpl/EquivalentEntityBlocks.hpp>
#include <stk_mesh/baseImpl/elementGraph/MeshDiagnostics.hpp>

namespace
{

struct split_element_info
{
    stk::mesh::EntityId localElementId;
    stk::mesh::EntityId remoteElementId;
    int neighboringProc;
};

struct TestCase
{
    std::string filename;
    int maxNumProcs;
    std::vector<split_element_info> expected_split_elements;
};

typedef std::vector<TestCase> TestCaseData;

const TestCaseData badDecomps =
{
  /* filename, #procs,      {local element id, remoted element id, remote proc} */
    {"Aef.e",     2,        { {3u, 2u, 1}, {2u, 3u, 0} }},
    {"ef.e",      2,        { {1u, 2u, 1}, {2u, 1u, 0} }},
    {"AefB.e",    3,        { {}, {2u, 3u, 2}, {3u, 2u, 1} }},
    {"AP.e",      2,        { {1u, 2u, 1}, {2u, 1u, 0}} }
};

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
            if(stk::parallel_machine_size(get_comm()) == testCase.maxNumProcs)
                EXPECT_THROW(test_one_case(testCase, auraOption), std::logic_error);
    }

    void test_one_case(const TestCase &testCase,
                       stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        stk::mesh::MetaData metaData;
        stk::mesh::BulkData bulkData(metaData, get_comm(), auraOption);
        SideTestUtil::read_and_decompose_mesh(testCase.filename, bulkData);

        std::map<stk::mesh::EntityId, std::pair<stk::mesh::EntityId, int> > splitCoincidentElements = stk::mesh::get_split_coincident_elements(bulkData);

        for(const auto& item : splitCoincidentElements)
        {
            stk::mesh::EntityId localElementId = item.first;
            stk::mesh::EntityId remoteElementId = item.second.first;
            int neighboringProc = item.second.second;

            EXPECT_EQ(testCase.expected_split_elements[bulkData.parallel_rank()].localElementId, localElementId);
            EXPECT_EQ(testCase.expected_split_elements[bulkData.parallel_rank()].remoteElementId, remoteElementId);
            EXPECT_EQ(testCase.expected_split_elements[bulkData.parallel_rank()].neighboringProc, neighboringProc);
        }

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


}
