#include "stk_util/environment/LogWithTimeAndMemory.hpp"
#include "stk_util/environment/Env.hpp"          // for outputP0
#include "stk_util/environment/WallTime.hpp"     // for wall_time
#include "stk_util/environment/memory_util.hpp"  // for get_memory_high_water_mark_across_proces...
#include "stk_util/parallel/Parallel.hpp"        // for MPI_Comm, ompi_communicator_t
#include "stk_util/util/human_bytes.hpp"         // for human_bytes
#include <cstddef>                               // for size_t
#include <iomanip>                               // for operator<<, setw, setfill, setprecision
#include <ostream>                               // for basic_ostream, operator<<, basic_ostream...
#include <string>                                // for operator<<, char_traits, string

namespace stk {

void log_with_time_and_memory(MPI_Comm communicator, const std::string &message, std::ostream& ostrm)
{
    size_t hwm_max = 0, hwm_min = 0, hwm_avg = 0;
    stk::get_memory_high_water_mark_across_processors(communicator, hwm_max, hwm_min, hwm_avg);

    static double startTime = stk::wall_time();
    double now = stk::wall_time();

    ostrm << "[time:" << std::right << std::fixed << std::setprecision(3) << std::setw(10)
          << now-startTime << " s, hwm:" << std::setfill(' ') << std::right << std::setw(8)
          << stk::human_bytes(hwm_avg) << "] "<< message << std::endl;
}

}
