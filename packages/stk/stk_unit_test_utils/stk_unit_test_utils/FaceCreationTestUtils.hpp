// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef STK_FACE_CREATION_TEST_UTILS_H
#define STK_FACE_CREATION_TEST_UTILS_H

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t, nullptr
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <string>                       // for string

#include "stk_mesh/base/MeshDiagnostics.hpp"
#include <stk_balance/fixSplitCoincidentElements.hpp>

using stk::unit_test_util::build_mesh;

namespace SideTestUtil
{

struct Side
{
    stk::mesh::EntityId elementId;
    unsigned sideOrdinal;
};

struct TestCase
{
    std::string filename;
    int maxNumProcs;
    size_t globalNumSides;
    std::vector<Side> sideSet;
};

typedef std::vector<TestCase> TestCaseData;

inline stk::mesh::Part & run_skin_mesh(stk::mesh::BulkData& bulkData, stk::mesh::Selector blocksToSkin)
{
    stk::mesh::Part &skin = bulkData.mesh_meta_data().declare_part("skin", bulkData.mesh_meta_data().side_rank());
    EXPECT_NO_FATAL_FAILURE(stk::mesh::create_exposed_block_boundary_sides(bulkData, blocksToSkin, {&skin}));
    return skin;
}

inline bool can_find_face_for_elem_side(const stk::mesh::BulkData& bulkData, stk::mesh::Entity element, unsigned sideOrdinal)
{
    stk::mesh::EntityRank sideRank = bulkData.mesh_meta_data().side_rank();
    auto iter = std::find(bulkData.begin_ordinals(element, sideRank),
                          bulkData.end_ordinals(element, sideRank),
                          sideOrdinal);
    return iter != bulkData.end_ordinals(element, sideRank);
}

inline void expect_side_exists_for_elem_side(const stk::mesh::BulkData& bulkData, const std::string &meshDescription, const Side& side)
{
    stk::mesh::Entity element = bulkData.get_entity(stk::topology::ELEM_RANK, side.elementId);
    if(bulkData.is_valid(element) && bulkData.bucket(element).owned())
    {
        EXPECT_TRUE(can_find_face_for_elem_side(bulkData, element, side.sideOrdinal)) << meshDescription << " couldn't find face for side: " << side.elementId << ", " << side.sideOrdinal << ".";
    }
}

inline void expect_all_sides_exist_for_elem_side(const stk::mesh::BulkData& bulkData,
                                          const std::string& meshDescription,
                                          const std::vector<Side> &sideSet)
{
    for(const Side &side : sideSet)
        expect_side_exists_for_elem_side(bulkData, meshDescription, side);
}


inline void read_and_decompose_mesh(const std::string &filename, stk::mesh::BulkData &bulkData)
{
    if(bulkData.parallel_rank() == 0)
        std::cout << "\t***** reading " << filename << " *****" << std::endl;
    stk::unit_test_util::read_from_serial_file_and_decompose(filename, bulkData, "cyclic");
}


inline void expect_global_num_sides_correct(const stk::mesh::BulkData& bulkData, const TestCase& testCase)
{
    std::vector<size_t> countsPerRank;
    stk::mesh::comm_mesh_counts(bulkData, countsPerRank);
    EXPECT_EQ(testCase.globalNumSides, countsPerRank[bulkData.mesh_meta_data().side_rank()]) << testCase.filename;
}

inline void expect_global_num_sides_in_part(const stk::mesh::BulkData& bulkData, size_t goldGlobalNumSides, const stk::mesh::Part &part)
{
    unsigned numLocalSides = stk::mesh::count_selected_entities(part & bulkData.mesh_meta_data().locally_owned_part(), bulkData.buckets(bulkData.mesh_meta_data().side_rank()));
    unsigned numGlobalSides = 0;
    stk::all_reduce_sum(bulkData.parallel(), &numLocalSides, &numGlobalSides, 1);
    EXPECT_EQ(goldGlobalNumSides, numGlobalSides) << "in part " << part.name();
}

inline void expect_exposed_sides_connected_as_specified_in_test_case(stk::mesh::BulkData& bulkData,
                                                                     const SideTestUtil::TestCase& testCase,
                                                                     stk::mesh::Selector skinnedThings,
                                                                     stk::mesh::Part &skinnedPart)
{
    SideTestUtil::expect_global_num_sides_in_part(bulkData, testCase.globalNumSides, skinnedPart);
    SideTestUtil::expect_all_sides_exist_for_elem_side(bulkData, testCase.filename, testCase.sideSet);
    EXPECT_TRUE(stk::mesh::check_exposed_block_boundary_sides(bulkData, skinnedThings, skinnedPart));
}

inline void expect_interior_sides_connected_as_specified_in_test_case(stk::mesh::BulkData& bulkData,
                                                                      const SideTestUtil::TestCase& testCase,
                                                                      stk::mesh::Selector skinnedThings,
                                                                      stk::mesh::Part &skinnedPart)
{
    SideTestUtil::expect_global_num_sides_in_part(bulkData, testCase.globalNumSides, skinnedPart);
    SideTestUtil::expect_all_sides_exist_for_elem_side(bulkData, testCase.filename, testCase.sideSet);
    EXPECT_TRUE(stk::mesh::check_interior_block_boundary_sides(bulkData, skinnedThings, skinnedPart));
}

inline void expect_all_sides_connected_as_specified_in_test_case(stk::mesh::BulkData& bulkData,
                                                                 const SideTestUtil::TestCase& testCase,
                                                                 stk::mesh::Selector skinnedThings,
                                                                 stk::mesh::Part &skinnedPart)
{
    SideTestUtil::expect_global_num_sides_in_part(bulkData, testCase.globalNumSides, skinnedPart);
    SideTestUtil::expect_all_sides_exist_for_elem_side(bulkData, testCase.filename, testCase.sideSet);
    EXPECT_TRUE(stk::mesh::check_all_sides(bulkData, skinnedThings, skinnedPart));
}

inline void expect_all_boundary_sides_connected_as_specified_in_test_case(stk::mesh::BulkData& bulkData,
                                                                 const SideTestUtil::TestCase& testCase,
                                                                 stk::mesh::Selector skinnedThings,
                                                                 stk::mesh::Part &skinnedPart)
{
    SideTestUtil::expect_global_num_sides_in_part(bulkData, testCase.globalNumSides, skinnedPart);
    SideTestUtil::expect_all_sides_exist_for_elem_side(bulkData, testCase.filename, testCase.sideSet);
//    EXPECT_TRUE(stk::mesh::check_all_boundary_sides(bulkData, skinnedThings, skinnedPart));
}


class SideCreationTester
{
public:
    SideCreationTester(MPI_Comm comm) : communicator(comm) {}
    virtual ~SideCreationTester() {}

    void run_all_test_cases(const SideTestUtil::TestCaseData &testCases, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        for(const SideTestUtil::TestCase& testCase : testCases)
            if(stk::parallel_machine_size(communicator) <= testCase.maxNumProcs)
                test_one_case(testCase, auraOption);
    }
protected:
    virtual void test_one_case(const SideTestUtil::TestCase &testCase,
                       stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        std::shared_ptr<stk::mesh::BulkData> bulkData = build_mesh(communicator, auraOption);
        SideTestUtil::read_and_decompose_mesh(testCase.filename, *bulkData);
        stk::balance::make_mesh_consistent_with_parallel_mesh_rule1(*bulkData);

        stk::mesh::SplitCoincidentInfo splitCoincidentElementsAfter = stk::mesh::get_split_coincident_elements(*bulkData);
        bool allOkAfterThisProc = splitCoincidentElementsAfter.size()==0;
        ASSERT_TRUE(allOkAfterThisProc);
        bool allOkEverywhereAfter = stk::is_true_on_all_procs(bulkData->parallel(), allOkAfterThisProc);
        EXPECT_TRUE(allOkEverywhereAfter);

        test_side_creation(*bulkData, testCase);
    }

    virtual void test_side_creation(stk::mesh::BulkData& bulkData,
                                    const SideTestUtil::TestCase& testCase) = 0;

    MPI_Comm communicator;
};

namespace simple_fields {

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
SideCreationTester
{
public:
    SideCreationTester(MPI_Comm comm) : communicator(comm) {}
    virtual ~SideCreationTester() {}

    void run_all_test_cases(const SideTestUtil::TestCaseData &testCases, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        for(const SideTestUtil::TestCase& testCase : testCases)
            if(stk::parallel_machine_size(communicator) <= testCase.maxNumProcs)
                test_one_case(testCase, auraOption);
    }
protected:
    virtual void test_one_case(const SideTestUtil::TestCase &testCase,
                       stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        std::shared_ptr<stk::mesh::BulkData> bulkData = build_mesh(communicator, auraOption);
        SideTestUtil::read_and_decompose_mesh(testCase.filename, *bulkData);
        stk::balance::make_mesh_consistent_with_parallel_mesh_rule1(*bulkData);

        stk::mesh::SplitCoincidentInfo splitCoincidentElementsAfter = stk::mesh::get_split_coincident_elements(*bulkData);
        bool allOkAfterThisProc = splitCoincidentElementsAfter.size()==0;
        ASSERT_TRUE(allOkAfterThisProc);
        bool allOkEverywhereAfter = stk::is_true_on_all_procs(bulkData->parallel(), allOkAfterThisProc);
        EXPECT_TRUE(allOkEverywhereAfter);

        test_side_creation(*bulkData, testCase);
    }

    virtual void test_side_creation(stk::mesh::BulkData& bulkData,
                                    const SideTestUtil::TestCase& testCase) = 0;

    MPI_Comm communicator;
};

} // namespace simple_fields

}

#endif

