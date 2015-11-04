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
            duration(0.0)
    {}

    void run_skin_mesh_performance_test()
    {
        stk::mesh::Selector thingToSkin = get_meta().universal_part();
        time_skin_mesh(thingToSkin);
        print_stats();
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

    std::string get_mesh_spec()
    {
        return unitTestUtils::getOption("-file", "NO_FILE_SPECIFIED");
    }

    stk::mesh::Part &skinPart;
    double duration;
};

TEST_F(StkPerformance, skin_mesh_with_auto_decomp)
{
    allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
    stk::unit_test_util::read_from_serial_file_and_decompose(get_mesh_spec(), get_bulk(), "rcb");

    run_skin_mesh_performance_test();
}

TEST_F(StkPerformance, skin_mesh)
{
    setup_mesh(get_mesh_spec(), stk::mesh::BulkData::AUTO_AURA);

    run_skin_mesh_performance_test();
}

