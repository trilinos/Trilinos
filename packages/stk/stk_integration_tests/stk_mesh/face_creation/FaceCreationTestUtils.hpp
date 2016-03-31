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
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture

#include "penso/penso.hpp"
#include "stk_mesh/baseImpl/elementGraph/MeshDiagnostics.hpp"
#include "stk_util/parallel/ParallelReduceBool.hpp"

namespace stk {
namespace mesh {
    class SideSetEntry;
    EntityVector get_locally_owned_sides_from_sideset(BulkData &bulkData, std::vector<SideSetEntry> &skinnedSideSet);
    bool is_sideset_equivalent_to_skin(BulkData &bulkData, EntityVector &sidesetSides, const Part& skinnedPart);
}}

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

inline void expect_side_exists_for_elem_side(const stk::mesh::BulkData& bulkData, const std::string &filename, const Side& side)
{
    stk::mesh::Entity element = bulkData.get_entity(stk::topology::ELEM_RANK, side.elementId);
    if(bulkData.is_valid(element) && bulkData.bucket(element).owned())
        EXPECT_TRUE(can_find_face_for_elem_side(bulkData, element, side.sideOrdinal))
                << filename << " couldn't find face for side: " << side.elementId << ", " << side.sideOrdinal << ".";
}

inline void expect_all_sides_exist_for_elem_side(const stk::mesh::BulkData& bulkData,
                                          const std::string& filename,
                                          const std::vector<Side> &sideSet)
{
    for(const Side &side : sideSet)
        expect_side_exists_for_elem_side(bulkData, filename, side);
}


inline void read_and_decompose_mesh(const std::string &filename, stk::mesh::BulkData &bulkData)
{
    if(bulkData.parallel_rank() == 0)
        std::cerr << "\t***** reading " << filename << " *****" << std::endl;
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
            {
                test_one_case(testCase, auraOption);
            }
    }
protected:
    virtual void test_one_case(const SideTestUtil::TestCase &testCase,
                       stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        stk::mesh::MetaData metaData;
        stk::mesh::BulkData bulkData(metaData, communicator, auraOption);
        SideTestUtil::read_and_decompose_mesh(testCase.filename, bulkData);
        penso::make_mesh_consistent_with_parallel_mesh_rule1(bulkData);

        stk::mesh::SplitCoincidentInfo splitCoincidentElementsAfter = stk::mesh::get_split_coincident_elements(bulkData);
        bool allOkAfterThisProc = splitCoincidentElementsAfter.size()==0;
        ASSERT_TRUE(allOkAfterThisProc);
        bool allOkEverywhereAfter = stk::is_true_on_all_procs(bulkData.parallel(), allOkAfterThisProc);
        EXPECT_TRUE(allOkEverywhereAfter);

        test_side_creation(bulkData, testCase);
    }

    virtual void test_side_creation(stk::mesh::BulkData& bulkData,
                                    const SideTestUtil::TestCase& testCase) = 0;

    MPI_Comm communicator;
};

}
