#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/GetEntities.hpp>  // for count_selected_entities
#include <stk_mesh/base/SkinMesh.hpp>   // for skin_mesh
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>

namespace
{

size_t get_num_global_faces(const stk::mesh::BulkData &bulk)
{
    std::vector<size_t> meshCounts;
    stk::mesh::comm_mesh_counts(bulk, meshCounts);
    return meshCounts[stk::topology::FACE_RANK];
}

class SkinMeshPerformance : public stk::unit_test_util::MeshFixture, public stk::unit_test_util::PerformanceTester
{
protected:
    SkinMeshPerformance() :
            stk::unit_test_util::MeshFixture(),
            stk::unit_test_util::PerformanceTester(get_comm()),
            thingToSkin(get_meta().universal_part()),
            skinPart(get_meta().declare_part("skinPart"))
    {
    }
    virtual void run_algorithm_to_time()
    {
        stk::mesh::skin_mesh(get_bulk(), thingToSkin, {&skinPart});
    }
    virtual size_t get_value_to_output_as_iteration_count()
    {
        return get_num_global_faces(get_bulk());
    }

    stk::mesh::Selector thingToSkin;
    stk::mesh::Part &skinPart;
};

std::string get_mesh_spec()
{
    return unitTestUtils::getOption("-file", "NO_FILE_SPECIFIED");
}

TEST_F(SkinMeshPerformance, read_mesh_with_auto_decomp)
{
    allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
    stk::unit_test_util::read_from_serial_file_and_decompose(get_mesh_spec(), get_bulk(), "rcb");

    run_performance_test();
}

TEST_F(SkinMeshPerformance, read_mesh)
{
    setup_mesh(get_mesh_spec(), stk::mesh::BulkData::AUTO_AURA);

    run_performance_test();
}

}

