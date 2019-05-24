#include <mpi.h>
#include <stk_util/environment/Env.hpp>
#include <stk_util/environment/memory_util.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/util/human_bytes.hpp>
#include <iomanip>
#include <string>

namespace stk {

void log_with_time_and_memory(MPI_Comm communicator, const std::string &message)
{
    static double startTime = stk::wall_time();
    double now = stk::wall_time();

    size_t hwm_max = 0, hwm_min = 0, hwm_avg = 0;
    stk::get_memory_high_water_mark_across_processors(communicator, hwm_max, hwm_min, hwm_avg);

    sierra::Env::outputP0() << "[time:" << std::fixed << std::setprecision(3) << std::setw(10) << now-startTime << " s, hwm:"
              << std::setfill(' ') << std::right << std::setw(8) << stk::human_bytes(hwm_avg) << "] "<< message << std::endl;
}

}
