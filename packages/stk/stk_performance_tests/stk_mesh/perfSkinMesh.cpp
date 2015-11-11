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
            enabledTimerSet(CHILDMASK1),
            rootTimer(createRootTimer("totalTestRuntime", enabledTimerSet)),
            childTimer1("skin_mesh", CHILDMASK1, rootTimer),
            duration(0.0)
    {
    }

    ~StkPerformance()
    {
        stk::diag::deleteRootTimer(rootTimer);
    }

    void run_skin_mesh_performance_test()
    {
        stk::mesh::Selector thingToSkin = get_meta().universal_part();
        time_skin_mesh(thingToSkin);
        print_stats(std::cerr);
    }

    void time_skin_mesh(stk::mesh::Selector thingToSkin)
    {
        stk::diag::TimeBlockSynchronized timerStartSynchronizedAcrossProcessors(childTimer1, communicator);
        rootTimer.start();
        double startTime = stk::wall_time();
        stk::mesh::skin_mesh(get_bulk(), thingToSkin, {&skinPart});
        duration = stk::wall_time() - startTime;
    }

    void print_stats(std::ostream& out)
    {
        print_output_for_pass_fail_test(out);
        print_output_for_graph_generation(out);
    }

    void print_output_for_graph_generation(std::ostream& out)
    {
        bool printTimingsOnlySinceLastPrint = false;
        stk::diag::printTimersTable(out, rootTimer, stk::diag::METRICS_ALL, printTimingsOnlySinceLastPrint, get_comm());
        stk::parallel_print_time_without_output_and_hwm(get_comm(), duration, out);
    }

    void print_output_for_pass_fail_test(std::ostream& out)
    {
        double maxTime = stk::get_max_time_across_procs(duration, get_comm());
        double maxHwmInMB = stk::get_max_hwm_across_procs(get_comm()) / (1024.0 * 1024.0);
        stk::print_stats_for_performance_compare(out, maxTime, maxHwmInMB, get_num_global_faces(), get_comm());
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

    const int CHILDMASK1 = 1;
    stk::diag::TimerSet enabledTimerSet;
    stk::diag::Timer rootTimer;
    stk::diag::Timer childTimer1;

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

