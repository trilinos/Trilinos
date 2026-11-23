/*
 * Copyright(C) 1999-2025 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
// concatenates EXODUS/GENESIS output from parallel processors to a single file

#include <algorithm>
#include <array>
#include <cfloat>
#include <climits>
#include <cmath>
#include <copy_string_cpp.h>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <fmt/chrono.h>
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "format_time.h"
#include "open_file_limit.h"
#include "sys_info.h"
#include "time_stamp.h"
#include "vector_data.h"

#define USE_STD_SORT 1
#if !USE_STD_SORT
#include "pdqsort.h"
#endif
#include "hwm.h"

// Enable SMART_ASSERT even in Release mode...
#define SMART_ASSERT_DEBUG_MODE 1
#include "smart_assert.h"

#include <exodusII.h>

using StringVector = std::vector<std::string>;

#include "EP_ExodusEntity.h"
#include "EP_ExodusFile.h"
#include "EP_Internals.h"
#include "EP_ObjectType.h"
#include "EP_SystemInterface.h"
#include "EP_Variables.h"
#include "EP_Version.h"

#if EX_API_VERS_NODOT <= 467
#error "Requires exodusII version 4.68 or later"
#endif

#include "add_to_log.h"

namespace suplib_cpp {

  extern double seacas_timer();

}

using suplib_cpp::copy_string;
using suplib_cpp::Data;
using suplib_cpp::format_time;
using suplib_cpp::get_hwm_memory_info;
using suplib_cpp::seacas_timer;
using suplib_cpp::sys_info;
using suplib_cpp::time_stamp;

// The main program templated to permit float/double transfer.
template <typename T, typename INT>
int epu(Excn::SystemInterface &interFace, int start_part, int part_count, int cycle);

class mpi
{
public:
  mpi(int argc, char *argv[])
  {
#if ENABLE_PARALLEL_EPU
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);           // CHECK: ALLOW MPI_COMM_WORLD
    MPI_Comm_size(MPI_COMM_WORLD, &epu_proc_count); // CHECK: ALLOW MPI_COMM_WORLD
#else
    (void)(argc);
    (void)(argv);
#endif
  }

  ~mpi()
  {
#if ENABLE_PARALLEL_EPU
    MPI_Finalize();
#endif
  }

  int rank{0};
  int epu_proc_count{1};
};

using namespace Excn;

int main(int argc, char *argv[])
{
  mpi my_mpi(argc, argv);
  int rank = my_mpi.rank;

  try {
    time_t begin_time = std::time(nullptr);
    SystemInterface::show_version(rank);
    if (rank == 0) {
#if ENABLE_PARALLEL_EPU
      fmt::print("\tParallel Capability Enabled.\n");
#else
      fmt::print("\tParallel Capability Not Enabled.\n");
#endif
    }

    SystemInterface interFace(rank);
    bool            execute = interFace.parse_options(argc, argv);

    if (!execute) {
      return EXIT_SUCCESS;
    }

    // Debug Options: (can be or'd together)
    //   1 -- time stamp
    //   2 -- check nodal variable consistency
    //   4 -- Element Blocks
    //   8 -- Nodes
    //  16 -- Sidesets
    //  32 -- Nodesets
    //  64 -- Edge Blocks
    // 128 -- Face Blocks
    // 256 -- exodus verbose.
    // 512 -- Check consistent global field values between processors
    unsigned int debug_level = interFace.debug();

    if ((debug_level & 256) != 0U) {
      ex_opts(EX_VERBOSE | EX_DEBUG);
    }
    else {
      ex_opts(0);
    }

    int start_part      = interFace.start_part();
    int processor_count = interFace.processor_count();

    int part_count = interFace.part_count();
    if (part_count <= 1) {
      fmt::print("INFO: Only one processor or part, no concatenation needed.\n");
      return EXIT_SUCCESS;
    }

    int error = 0;
    if (my_mpi.epu_proc_count > 1) {
      interFace.subcycle(my_mpi.epu_proc_count);

      int per_proc = processor_count / my_mpi.epu_proc_count;
      int extra    = processor_count % my_mpi.epu_proc_count;

      part_count = per_proc + (rank < extra ? 1 : 0);

      if (rank < extra) {
        start_part = (per_proc + 1) * rank;
      }
      else {
        start_part = (per_proc + 1) * extra + per_proc * (rank - extra);
      }

      SMART_ASSERT(start_part + part_count <= processor_count);

      ExodusFile::initialize(interFace, start_part, part_count, rank, false);

      if (ExodusFile::io_word_size() == 4) { // Reals are floats
        if (interFace.int64()) {
          error = epu<float, int64_t>(interFace, start_part, part_count, rank);
        }
        else {
          error = epu<float, int>(interFace, start_part, part_count, rank);
        }
      }
      else { // Reals are doubles
        if (interFace.int64()) {
          error = epu<double, int64_t>(interFace, start_part, part_count, rank);
        }
        else {
          error = epu<double, int>(interFace, start_part, part_count, rank);
        }
      }

      ExodusFile::close_all();
#if ENABLE_PARALLEL_EPU
      MPI_Barrier(MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD
#endif
    }
    else {
      int max_open_file = open_file_limit() - 1; // -1 for output exodus file.

      // Only used to test the auto subcycle without requiring thousands of files...
      if (interFace.max_open_files() > 0) {
        max_open_file = interFace.max_open_files();
      }

      if (interFace.is_auto() && interFace.subcycle() < 0 && processor_count > max_open_file &&
          part_count == processor_count && interFace.cycle() == -1) {
        // Rule of thumb -- number of subcycles = cube_root(processor_count);
        // if that value > max_open_file, then use square root.
        // if that is still too large, just do no subcycles... and implement
        // a recursive subcycling capability at some point...
        int sub_cycle_count = static_cast<int>((std::pow(processor_count, 1.0 / 3) + 0.9));
        if (((processor_count + sub_cycle_count - 1) / sub_cycle_count) > max_open_file) {
          sub_cycle_count = static_cast<int>(std::sqrt(processor_count));
        }

        if (((processor_count + sub_cycle_count - 1) / sub_cycle_count) < max_open_file) {
          interFace.subcycle(sub_cycle_count);
          if (rank == 0) {
            fmt::print("\tAutomatically activating subcyle mode\n\tNumber of processors ({}) "
                       "exceeds open file limit ({}).\n"
                       "\tUsing --subcycle={}\n\n",
                       processor_count, max_open_file, sub_cycle_count);
          }
          interFace.subcycle_join(true);
        }
      }

      int cycle = interFace.cycle();
      if (interFace.subcycle() >= 0) {
        start_part = 0;
        int cycles = interFace.subcycle();
        if (cycles > 0) {
          // use the specified number of cycles...
          part_count = (processor_count + cycles - 1) / cycles;
          if (cycle >= 0) {
            start_part = cycle * part_count;
          }
        }

        // Sanity check...
        if (part_count < 1) {
          throw std::runtime_error(
              "ERROR: (EPU) The subcycle specification results in less than 1 part per "
              "cycle which is not allowed.\n");
        }
        interFace.subcycle((processor_count + part_count - 1) / part_count);

        if (start_part + part_count > processor_count) {
          part_count = processor_count - start_part;
        }
      }

      if (cycle < 0) {
        cycle = 0;
      }
      while (start_part < processor_count) {

        if (start_part + part_count > processor_count) {
          part_count = processor_count - start_part;
        }

        SMART_ASSERT(part_count > 0);
        SMART_ASSERT(start_part + part_count <= processor_count);

        ExodusFile::initialize(interFace, start_part, part_count, cycle, false);

        if (ExodusFile::io_word_size() == 4) { // Reals are floats
          if (interFace.int64()) {
            error = epu<float, int64_t>(interFace, start_part, part_count, cycle++);
          }
          else {
            error = epu<float, int>(interFace, start_part, part_count, cycle++);
          }
        }
        else { // Reals are doubles
          if (interFace.int64()) {
            error = epu<double, int64_t>(interFace, start_part, part_count, cycle++);
          }
          else {
            error = epu<double, int>(interFace, start_part, part_count, cycle++);
          }
        }

        start_part += part_count;
        ExodusFile::close_all();
        if (interFace.subcycle() < 0 || (interFace.subcycle() > 0 && interFace.cycle() >= 0)) {
          break;
        }
      }
    }

    if (interFace.subcycle() > 0 && interFace.cycle() < 0 && interFace.subcycle_join() &&
        rank == 0) {
      // Now, join the subcycled parts into a single file...
      start_part = 0;
      part_count = interFace.subcycle();
      interFace.subcycle(0);
      interFace.processor_count(part_count);
      interFace.step_min(1);
      interFace.step_max(INT_MAX);
      interFace.step_interval(1);

      ExodusFile::initialize(interFace, start_part, part_count, 0, true);
      if (ExodusFile::io_word_size() == 4) { // Reals are floats
        if (interFace.int64()) {
          error = epu<float, int64_t>(interFace, start_part, part_count, 0);
        }
        else {
          error = epu<float, int>(interFace, start_part, part_count, 0);
        }
      }
      else { // Reals are doubles
        if (interFace.int64()) {
          error = epu<double, int64_t>(interFace, start_part, part_count, 0);
        }
        else {
          error = epu<double, int>(interFace, start_part, part_count, 0);
        }
      }

      bool delete_temp = error == 0 && !interFace.keep_temporary();
      ExodusFile::handle_temporary_files(delete_temp);
    }

    if (error == 0 && interFace.remove_file_per_rank_files()) {
      ExodusFile::unlink_input_files();
    }

    time_t end_time = std::time(nullptr);
    if (rank == 0) {
      add_to_log(argv[0], static_cast<int>(end_time - begin_time));
    }
    return error;
  }
  catch (std::exception &e) {
    fmt::print(stderr, "{}\n", e.what());
    return EXIT_FAILURE;
  }
}
