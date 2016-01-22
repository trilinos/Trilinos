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
    if(bulkData.is_valid(element))
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

inline void expect_global_num_sides_in_part(const stk::mesh::BulkData& bulkData, const TestCase& testCase, const stk::mesh::Part &skinnedPart)
{
    unsigned numLocalSkinnedSides = stk::mesh::count_selected_entities(skinnedPart & bulkData.mesh_meta_data().locally_owned_part(), bulkData.buckets(bulkData.mesh_meta_data().side_rank()));
    unsigned numGlobalSkinnedSides = 0;
    stk::all_reduce_sum(bulkData.parallel(), &numLocalSkinnedSides, &numGlobalSkinnedSides, 1);
    EXPECT_EQ(testCase.globalNumSides, numGlobalSkinnedSides);
}

class SideCreationTester
{
public:
    SideCreationTester(MPI_Comm comm) : communicator(comm) {}
    virtual ~SideCreationTester() {}

    void run_all_test_cases(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        for(const SideTestUtil::TestCase& testCase : get_test_cases())
            if(stk::parallel_machine_size(communicator) <= testCase.maxNumProcs)
                test_one_case(testCase, auraOption);
    }

protected:
    void test_one_case(const SideTestUtil::TestCase &testCase,
                               stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        stk::mesh::MetaData metaData;
        stk::mesh::BulkData bulkData(metaData, communicator, auraOption);
        SideTestUtil::read_and_decompose_mesh(testCase.filename, bulkData);
        test_side_creation(bulkData, testCase);
    }

    virtual SideTestUtil::TestCaseData get_test_cases() = 0;

    virtual void test_side_creation(stk::mesh::BulkData& bulkData,
                                    const SideTestUtil::TestCase& testCase) = 0;
private:
    MPI_Comm communicator;
};

}
