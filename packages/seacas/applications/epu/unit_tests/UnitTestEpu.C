// Copyright(C) 1999-2020, 2022, 2025 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifdef SEACAS_HAVE_MPI
#include "mpi.h"
#endif
#include "gtest/gtest.h"

#include "ArgvBuilder.h"
#include "FileUtils.h"
#include "MeshFixture.h"

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
#include <fmt/ranges.h>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <exodusII.h>

#include <Ionit_Initializer.h>
#include <Ioss_CopyDatabase.h>
#include <Ioss_Enumerate.h>
#include <Ioss_MeshCopyOptions.h>
#include <Ioss_SmartAssert.h>
#include <Ioss_SubSystem.h>
#include <Ioss_Transform.h>
#include <Ioss_Utils.h>
#include <tokenize.h>

#include "format_time.h"
#include "open_file_limit.h"
#include "sys_info.h"
#include "time_stamp.h"
#include "vector_data.h"

#include <exodusII.h>

using StringVector = std::vector<std::string>;
using RegionVector = std::vector<Ioss::Region *>;

#include "EP_ExodusFile.h"
#include "EP_SystemInterface.h"

#if EX_API_VERS_NODOT <= 467
#error "Requires exodusII version 4.68 or later"
#endif

// The main program templated to permit float/double transfer.
template <typename T, typename INT>
int epu(Excn::SystemInterface &interFace, int start_part, int part_count, int cycle);

namespace {
  Ioss::MeshCopyOptions set_mesh_copy_options(Excn::SystemInterface &interFace)
  {
    Ioss::MeshCopyOptions options{};
    //  options.selected_times       = interFace.selected_times;
    //  options.selected_steps       = interFace.selected_steps;
    //  options.rel_tolerance        = interFace.rel_tolerance;
    //  options.abs_tolerance        = interFace.abs_tolerance;
    //  options.tol_floor            = interFace.tol_floor;
    options.verbose        = false;
    options.output_summary = true;
    //  options.memory_statistics    = interFace.memory_statistics;
    options.debug       = false; // interFace.debug() > 0;
    options.ints_64_bit = interFace.int64();
    //  options.delete_timesteps     = interFace.delete_timesteps;
    //  options.sort_times           = interFace.sort_times;
    //  options.shuffle_times        = interFace.shuffle_times;
    //  options.minimum_time         = interFace.minimum_time;
    //  options.maximum_time         = interFace.maximum_time;
    //  options.time_scale           = interFace.time_scale;
    //  options.time_offset          = interFace.time_offset;
    options.data_storage_type = 1;
    //  options.delay                = interFace.timestep_delay;
    //  options.reverse              = interFace.reverse;
    options.add_proc_id = interFace.add_processor_id_field();
    //  options.boundary_sideset     = interFace.boundary_sideset;
    //  options.ignore_qa_info       = interFace.ignore_qa_info;
    //  options.omitted_blocks       = !interFace.omitted_blocks.empty();
    //  options.selected_change_sets = interFace.selectedChangeSets;
    //
    //  options.omitted_sets = interFace.omitted_sets;
    Ioss::sort(options.omitted_sets);
    for (auto &name : options.omitted_sets) {
      name = Ioss::Utils::lowercase(name);
    }
    return options;
  }

  Ioss::PropertyManager set_properties(Excn::SystemInterface &interFace)
  {
    Ioss::PropertyManager properties;

    if (interFace.int64()) {
      properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
      properties.add(Ioss::Property("INTEGER_SIZE_API", 8));
    }

    if (interFace.compress_data() > 0 || interFace.szip() || interFace.quantize() ||
        interFace.zlib() || interFace.zstd() || interFace.bz2()) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
      properties.add(Ioss::Property("COMPRESSION_LEVEL", interFace.compress_data()));

      if (interFace.zlib()) {
        properties.add(Ioss::Property("COMPRESSION_METHOD", "zlib"));
      }
      else if (interFace.szip()) {
        properties.add(Ioss::Property("COMPRESSION_METHOD", "szip"));
      }
      else if (interFace.zstd()) {
        properties.add(Ioss::Property("COMPRESSION_METHOD", "zstd"));
      }
      else if (interFace.bz2()) {
        properties.add(Ioss::Property("COMPRESSION_METHOD", "bzip2"));
      }

      if (interFace.quantize()) {
        properties.add(Ioss::Property("COMPRESSION_QUANTIZE_NSD", interFace.quantize_nsd()));
      }
    }

    if (interFace.use_netcdf4()) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
    }

    if (interFace.use_netcdf5()) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf5"));
    }

    if (interFace.debug()) {
      properties.add(Ioss::Property("LOGGING", 1));
    }

    if (!interFace.map_element_ids()) {
      std::cout << "Ignore elem map\n";
      properties.add(Ioss::Property("IGNORE_ELEM_MAP", true));
    }

    if (!interFace.map_edge_ids()) {
      std::cout << "Ignore edge map\n";
      properties.add(Ioss::Property("IGNORE_EDGE_MAP", true));
    }
    if (!interFace.map_face_ids()) {
      std::cout << "Ignore face map\n";
      properties.add(Ioss::Property("IGNORE_FACE_MAP", true));
    }

    return properties;
  }

  void append_properties_from_region(Ioss::Region *inputRegion, Excn::SystemInterface &interFace,
                                     Ioss::PropertyManager &properties)
  {
    auto inputDB = inputRegion->get_database();

    // Get length of longest name on input file...
    int max_name_length = inputDB->maximum_symbol_length();
    if (max_name_length > 0) {
      properties.add(Ioss::Property("MAXIMUM_NAME_LENGTH", max_name_length));
    }

    // Get integer size being used on the input file and propagate
    // to output file...
    int int_byte_size_api = inputDB->int_byte_size_api();
    if (!properties.exists("INTEGER_SIZE_API")) {
      if (!interFace.int64()) {
        properties.add(Ioss::Property("INTEGER_SIZE_DB", 4));
      }
      properties.add(Ioss::Property("INTEGER_SIZE_API", int_byte_size_api));
    }
  }

} // namespace
namespace {

  using namespace Excn;

  class EpuTester : public utest_util::MeshFixture
  {
  public:
    EpuTester() : utest_util::MeshFixture(3)
    {
      Ioss::Init::Initializer io;
      setup_empty_mesh();
    }

    ~EpuTester() {}

  protected:
    std::string m_outputFile;
    std::string m_propertyName{"MATERIAL_PROPERTY"};
    std::string m_propertyValue{"KRYPTONITE"};

  protected:
    int run_parallel(SystemInterface &interFace)
    {
      int rank = interFace.processor_rank();

      int start_part      = interFace.start_part();
      int processor_count = interFace.processor_count();

      int part_count = interFace.part_count();
      interFace.subcycle(get_parallel_size());

      int per_proc = processor_count / get_parallel_size();
      int extra    = processor_count % get_parallel_size();

      part_count = per_proc + (rank < extra ? 1 : 0);

      if (rank < extra) {
        start_part = (per_proc + 1) * rank;
      }
      else {
        start_part = (per_proc + 1) * extra + per_proc * (rank - extra);
      }

      SMART_ASSERT(start_part + part_count <= processor_count);

      ExodusFile::initialize(interFace, start_part, part_count, rank, false);

      int error = 0;
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
      MPI_Barrier(get_comm());
#endif

      return error;
    }

    int run_serial_on_proc_0(SystemInterface &interFace)
    {
      int error           = 0;
      int start_part      = interFace.start_part();
      int processor_count = interFace.processor_count();
      int part_count      = interFace.part_count();
      int rank            = interFace.processor_rank();

      if (get_parallel_rank() == 0) {
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

      return error;
    }

    int finish_subcycle(SystemInterface &interFace)
    {
      int rank  = interFace.processor_rank();
      int error = 0;

      if (interFace.subcycle() > 0 && interFace.cycle() < 0 && interFace.subcycle_join() &&
          rank == 0) {
        // Now, join the subcycled parts into a single file...
        int start_part = 0;
        int part_count = interFace.subcycle();
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

      return error;
    }

    void setup_two_proc_epu(SystemInterface &interFace)
    {
      std::string baseName       = "twoHexEpu";
      std::string uniqueBaseName = utest_util::unique_basename(baseName);

      m_outputFile = utest_util::unique_filename(baseName, "g");

      int rank = get_parallel_rank();
      SystemInterface::show_version(rank);

      std::string meshDesc = "textmesh:"
                             "0,1,HEX_8,1,2,3,4,5,6,7,8,block_1\n"
                             "1,2,HEX_8,5,6,7,8,9,10,11,12,block_1"
                             "|coordinates:   0,0,0, 1,0,0, 1,1,0, 0,1,0, 0,0,1, 1,0,1, 1,1,1, "
                             "0,1,1, 0,0,2, 1,0,2, 1,1,2, 0,1,2";

      setup_mesh(meshDesc);

      // Add a material property to block_1 from textmesh
      auto inputRegion = get_mesh().get_region();
      auto inputDB     = get_mesh().get_database();

      add_material_property_to_element_block(inputRegion, "block_1", m_propertyName,
                                             m_propertyValue);

      utest_util::ArgvBuilder builder;
      builder.addArgument("epu_unit_test");
      builder.addArgument("-extension");
      builder.addArgument("g");
      builder.addArgument("-processor_count");
      builder.addArgument("2");
      builder.addArgument(uniqueBaseName.c_str());
      builder.addArgument("-remove_file_per_rank_files");
      builder.addArgument("-verify_valid_file");
      builder.addArgument("-debug");
      builder.addArgument("512");

      EXPECT_TRUE(interFace.parse_options(builder.argc(), builder.argv()));
      interFace.set_output_filename(m_outputFile);

      Ioss::PropertyManager properties = set_properties(interFace);
      append_properties_from_region(inputRegion, interFace, properties);

      if (inputDB->int_byte_size_api() == 8) {
        interFace.set_int64();
      }

      // Write the textmesh to file
      Ioss::MeshCopyOptions options = set_mesh_copy_options(interFace);
      write_region_to_file(inputRegion, properties, options, m_outputFile);

      // Load each individual file and test for material property
      test_property_from_file(m_outputFile, m_propertyName, m_propertyValue);

      if ((interFace.debug() & 256) != 0U) {
        ex_opts(EX_VERBOSE | EX_DEBUG);
      }
      else {
        ex_opts(0);
      }
    }
  };

  TEST_F(EpuTester, emptyTest) {}

  TEST_F(EpuTester, epuTwoHexMesh_serial)
  {
    if (get_parallel_size() != 2)
      GTEST_SKIP();

    SystemInterface interFace(get_parallel_rank());
    setup_two_proc_epu(interFace);

    bool runInParallel = false;

    int error = 0;
    if (runInParallel) {
      error = run_parallel(interFace);
    }
    else {
      error = run_serial_on_proc_0(interFace);
    }
    EXPECT_TRUE(error == 0);

    error = finish_subcycle(interFace);
    EXPECT_TRUE(error == 0);

    if (error == 0 && interFace.remove_file_per_rank_files()) {
      ExodusFile::unlink_input_files();
    }

    if (get_parallel_rank() == 0) {
      // Load epu'd file and test for material property
      test_property_from_file(Ioss::ParallelUtils::comm_self(), m_outputFile, m_propertyName,
                              m_propertyValue);

      // Delete epu'd file
      unlink(m_outputFile.c_str());
    }

#ifdef SEACAS_HAVE_MPI
    // Prevent early deletion of files
    MPI_Barrier(get_comm());
#endif

    // Delete parallel output files from textmesh
    Ioss::ParallelUtils util(get_comm());
    std::string decodedFileName = util.decode_filename(m_outputFile, (get_parallel_size() > 1));

    unlink(decodedFileName.c_str());
  }

} // namespace
