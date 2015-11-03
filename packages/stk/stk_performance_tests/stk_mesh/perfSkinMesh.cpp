#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_mesh/base/SkinMesh.hpp>   // for skin_mesh
#include <stk_mesh/base/Comm.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/getOption.h>

class StkPerformance : public stk::unit_test_util::MeshFixture
{
protected:
    StkPerformance() :
            skinPart(get_meta().declare_part("skinPart")),
            numSkinFaces(0),
            duration(0.0)
    {}

    void test_skin_mesh()
    {
        stk::mesh::Selector thingToSkin = get_meta().universal_part();
        create_faces_with_skin_mesh(thingToSkin);
        EXPECT_GT(numSkinFaces, 0u);
    }

    void create_faces_with_skin_mesh(stk::mesh::Selector thingToSkin)
    {
        unsigned numFacesBefore = stk::mesh::count_selected_entities(skinPart, get_bulk().buckets(stk::topology::FACE_RANK));
        time_skin_mesh(thingToSkin);
        unsigned numFacesAfter = stk::mesh::count_selected_entities(skinPart, get_bulk().buckets(stk::topology::FACE_RANK));
        numSkinFaces = numFacesAfter - numFacesBefore;
    }

    void time_skin_mesh(stk::mesh::Selector thingToSkin)
    {
        double startTime = stk::wall_time();
        stk::mesh::skin_mesh(get_bulk(), thingToSkin, {&skinPart});
        duration = stk::wall_time() - startTime;
    }

    void print_stats()
    {
        double maxTime = stk::get_max_time_across_procs(duration, get_comm());
        double maxHwmInMB = stk::get_max_hwm_across_procs(get_comm()) / (1024.0 * 1024.0);
        stk::print_stats_for_performance_compare(std::cerr, maxTime, maxHwmInMB, get_num_global_faces(), get_comm());
    }

    size_t get_num_global_faces()
    {
        std::vector<size_t> meshCounts;
        stk::mesh::comm_mesh_counts(get_bulk(), meshCounts);
        return meshCounts[stk::topology::FACE_RANK];
    }

    stk::mesh::Part &skinPart;
    unsigned numSkinFaces;
    double duration;
};

TEST_F(StkPerformance, skin_mesh)
{
    const std::string meshSpec = unitTestUtils::getOption("-file", "NO_FILE_SPECIFIED");
    allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
    stk::unit_test_util::read_from_serial_file_and_decompose(meshSpec, get_bulk(), "rcb");

    test_skin_mesh();
    print_stats();
}
