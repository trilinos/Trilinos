#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_mesh/base/SkinMesh.hpp>   // for skin_mesh
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
        if(get_bulk().parallel_rank() == 0)
            std::cerr << "skin faces = " << numSkinFaces << std::endl;
        stk::parallel_print_time_for_performance_compare(get_comm(), duration, std::cerr);
    }

    stk::mesh::Part &skinPart;
    unsigned numSkinFaces;
    double duration;
};

TEST_F(StkPerformance, skin_mesh)
{
    const std::string meshSpec = unitTestUtils::getOption("-file", "NO_FILE_SPECIFIED");
    setup_mesh(meshSpec, stk::mesh::BulkData::AUTO_AURA);

    test_skin_mesh();
    print_stats();
}
