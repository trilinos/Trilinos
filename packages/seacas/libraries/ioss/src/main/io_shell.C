// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include <Ionit_Initializer.h>
#include <Ioss_CodeTypes.h>
#include <Ioss_FileInfo.h>
#include <Ioss_MeshCopyOptions.h>
#include <Ioss_MeshType.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_ScopeGuard.h>
#include <Ioss_SerializeIO.h>
#include <Ioss_SubSystem.h>
#include <Ioss_SurfaceSplit.h>
#include <Ioss_Utils.h>
#include <fmt/format.h>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <unistd.h>

#include "shell_interface.h"

// ========================================================================

namespace {
  std::string codename;
  std::string version = "5.1";

  bool mem_stats = false;

  void file_copy(IOShell::Interface &interFace, int rank);

  Ioss::PropertyManager set_properties(IOShell::Interface &interFace);
} // namespace

int main(int argc, char *argv[])
{
  int rank     = 0;
  int num_proc = 1;
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  ON_BLOCK_EXIT(MPI_Finalize);
#endif

#ifdef SEACAS_HAVE_KOKKOS
  Kokkos::ScopeGuard kokkos(argc, argv);
#endif

  IOShell::Interface interFace;
  bool               success = interFace.parse_options(argc, argv);
  if (!success) {
    exit(EXIT_FAILURE);
  }

  codename = interFace.options_.basename(argv[0]);

  Ioss::SerializeIO::setGroupFactor(interFace.serialize_io_size);
  mem_stats = interFace.memory_statistics;

  Ioss::Init::Initializer io;

  std::string in_file  = interFace.inputFile[0];
  std::string out_file = interFace.outputFile;

  if (rank == 0 && !interFace.quiet) {
    fmt::print(stderr,
               "Input:    '{}', Type: {}\n"
               "Output:   '{}', Type: {}\n\n",
               in_file, interFace.inFiletype, out_file, interFace.outFiletype);
  }

#ifdef SEACAS_HAVE_KOKKOS
  if (rank == 0)
    fmt::print(stderr, "Kokkos default execution space configuration:\n");
  Kokkos::DefaultExecutionSpace::print_configuration(std::cerr, false);
  if (rank == 0)
    fmt::print(stderr, "\n");
#endif

  double begin = Ioss::Utils::timer();
  try {
    file_copy(interFace, rank);
  }
  catch (std::exception &e) {
    if (rank == 0) {
      fmt::print(stderr, "\n{}\n\nio_shell terminated due to exception\n", e.what());
    }
    exit(EXIT_FAILURE);
  }

#ifdef SEACAS_HAVE_MPI
  Ioss::ParallelUtils parallel(MPI_COMM_WORLD);
  parallel.barrier();
#endif
  double end = Ioss::Utils::timer();

  if (rank == 0 && !interFace.quiet) {
    if (num_proc > 1) {
      fmt::print(stderr, "\n\n\tTotal Execution time = {:.5} seconds on {} processors.\n",
                 end - begin, num_proc);
    }
    else {
      fmt::print(stderr, "\n\n\tTotal Execution time = {:.5} seconds.\n", end - begin);
    }
  }
  if (mem_stats) {
    int64_t MiB = 1024 * 1024;
#ifdef SEACAS_HAVE_MPI
    int64_t min, max, avg;
    int64_t hwmin, hwmax, hwavg;
    parallel.memory_stats(min, max, avg);
    parallel.hwm_memory_stats(hwmin, hwmax, hwavg);
    if (rank == 0) {
      fmt::print(stderr, "\n\tCurrent Memory: {:n}M  {:n}M  {:n}M\n", min / MiB, max / MiB,
                 avg / MiB);
      fmt::print(stderr, "\tHigh Water Memory: {:n}M  {:n}M  {:n}M\n", hwmin / MiB, hwmax / MiB,
                 hwavg / MiB);
    }
#else
    int64_t mem = Ioss::Utils::get_memory_info();
    int64_t hwm = Ioss::Utils::get_hwm_memory_info();
    if (rank == 0) {
      fmt::print(stderr,
                 "\n\tCurrent Memory:    {:n}M\n"
                 "\tHigh Water Memory: {:n}M\n",
                 mem / MiB, hwm / MiB);
    }
#endif
  }
  if (rank == 0) {
    fmt::print(stderr, "\n{} execution successful.\n", codename);
  }
  return EXIT_SUCCESS;
}

namespace {
  void file_copy(IOShell::Interface &interFace, int rank)
  {
    Ioss::PropertyManager properties = set_properties(interFace);

    bool first = true;
    for (const auto &inpfile : interFace.inputFile) {

      //========================================================================
      // INPUT Database...
      //========================================================================
      Ioss::DatabaseIO *dbi = Ioss::IOFactory::create(
          interFace.inFiletype, inpfile, Ioss::READ_MODEL, (MPI_Comm)MPI_COMM_WORLD, properties);
      if (dbi == nullptr || !dbi->ok(true)) {
        std::exit(EXIT_FAILURE);
      }

      if (mem_stats) {
        dbi->progress("Database Creation");
      }
      if (!interFace.lower_case_variable_names) {
        dbi->set_lower_case_variable_names(false);
      }
      if (interFace.outFiletype == "cgns") {
        // CGNS stores BCs (SideSets) on the zones which
        // correspond to element blocks.  If split input sideblocks
        // by element block, then output is much easier.
        dbi->set_surface_split_type(Ioss::SPLIT_BY_ELEMENT_BLOCK);
      }
      else {
        dbi->set_surface_split_type(Ioss::int_to_surface_split(interFace.surface_split_type));
      }
      dbi->set_field_separator(interFace.fieldSuffixSeparator);

      dbi->set_field_recognition(!interFace.disable_field_recognition);

      if (interFace.ints_64_bit) {
        dbi->set_int_byte_size_api(Ioss::USE_INT64_API);
      }

      if (!interFace.groupName.empty()) {
        bool success = dbi->open_group(interFace.groupName);
        if (!success) {
          if (rank == 0) {
            fmt::print(stderr, "ERROR: Unable to open group '{}' in file '{}'\n",
                       interFace.groupName, inpfile);
          }
          return;
        }
      }

      // NOTE: 'region' owns 'db' pointer at this time...
      Ioss::Region region(dbi, "region_1");

      if (region.mesh_type() == Ioss::MeshType::HYBRID) {
        fmt::print(stderr,
                   "\nERROR: io_shell does not support '{}' meshes. Only 'Unstructured' or "
                   "'Structured' mesh is supported at this time.\n",
                   region.mesh_type_string());
        return;
      }

      // Get length of longest name on input file...
      int max_name_length = dbi->maximum_symbol_length();
      if (max_name_length > 0) {
        properties.add(Ioss::Property("MAXIMUM_NAME_LENGTH", max_name_length));
      }

      // Get integer size being used on the input file and propgate
      // to output file...
      int int_byte_size_api = dbi->int_byte_size_api();
      if (!properties.exists("INTEGER_SIZE_API")) {
        if (interFace.ints_32_bit) {
          properties.add(Ioss::Property("INTEGER_SIZE_DB", 4));
        }
        properties.add(Ioss::Property("INTEGER_SIZE_API", int_byte_size_api));
      }
      if (int_byte_size_api == 8) {
        interFace.ints_64_bit = true;
      }
      //========================================================================
      // OUTPUT Database...
      //========================================================================
      bool append = false;
      if (interFace.append_step < std::numeric_limits<int>::max()) {
        properties.add(Ioss::Property("APPEND_OUTPUT_AFTER_STEP", interFace.append_step));
        append = true;
      }

      if (interFace.append_time < std::numeric_limits<double>::max()) {
        properties.add(Ioss::Property("APPEND_OUTPUT_AFTER_TIME", interFace.append_time));
        append = true;
      }

      if (append) {
        properties.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_APPEND));
      }

      if (interFace.minimize_open_files) {
        properties.add(Ioss::Property("MINIMIZE_OPEN_FILES", "ON"));
      }

      Ioss::MeshCopyOptions options{};
      options.verbose           = !interFace.quiet;
      options.memory_statistics = interFace.memory_statistics;
      options.debug             = interFace.debug;
      options.ints_64_bit       = interFace.ints_64_bit;
      options.delete_timesteps  = interFace.delete_timesteps;
      options.minimum_time      = interFace.minimum_time;
      options.maximum_time      = interFace.maximum_time;
      options.data_storage_type = interFace.data_storage_type;
      options.delay             = interFace.timestep_delay;
      options.reverse           = interFace.reverse;
      options.add_proc_id       = interFace.add_processor_id_field;

      size_t ts_count = 0;
      if (region.property_exists("state_count") &&
          region.get_property("state_count").get_int() > 0) {
        ts_count = region.get_property("state_count").get_int();
      }

      int flush_interval = interFace.flush_interval; // Default is zero -- do not flush until end
      properties.add(Ioss::Property("FLUSH_INTERVAL", flush_interval));

      if (interFace.split_times == 0 || interFace.delete_timesteps || ts_count == 0 || append ||
          interFace.inputFile.size() > 1) {
        Ioss::DatabaseIO *dbo =
            Ioss::IOFactory::create(interFace.outFiletype, interFace.outputFile,
                                    Ioss::WRITE_RESTART, (MPI_Comm)MPI_COMM_WORLD, properties);
        if (dbo == nullptr || !dbo->ok(true)) {
          std::exit(EXIT_FAILURE);
        }

        // NOTE: 'output_region' owns 'dbo' pointer at this time
        Ioss::Region output_region(dbo, "region_2");
        // Set the qa information...
        output_region.property_add(Ioss::Property(std::string("code_name"), codename));
        output_region.property_add(Ioss::Property(std::string("code_version"), version));

        if (interFace.inputFile.size() > 1) {
          properties.add(Ioss::Property("APPEND_OUTPUT", Ioss::DB_APPEND_GROUP));

          if (!first) {
            // Putting each file into its own output group...
            // The name of the group will be the basename portion of the filename...
            Ioss::FileInfo file(inpfile);
            dbo->create_subgroup(file.tailname());
          }
          else {
            first = false;
          }
        }

        // Do normal copy...
        Ioss::Utils::copy_database(region, output_region, options);
        if (mem_stats) {
          dbo->release_memory();
        }
      }
      else {
        // We are splitting out the timesteps into separate files.
        // Each file will contain `split_times` timesteps. If
        // 'split_cyclic` is > 0, then recycle filenames
        // (A,B,C,A,B,C,...) otherwise
        // keep creating new filenames (0001, 0002, 0003, ...)

        // Get list of all times on input database...
        std::vector<double> times;
        for (size_t step = 0; step < ts_count; step++) {
          double time = region.get_state_time(step + 1);
          if (time < interFace.minimum_time) {
            continue;
          }
          if (time > interFace.maximum_time) {
            break;
          }
          times.push_back(time);
        }
        ts_count = times.size();

        int splits = (ts_count + interFace.split_times - 1) / interFace.split_times;
        int width  = std::to_string(splits).length();
        for (int split = 0; split < splits; split++) {
          int step_min = split * interFace.split_times;
          int step_max = step_min + interFace.split_times - 1;
          if (step_max >= (int)times.size()) {
            step_max = (int)times.size() - 1;
          }
          options.minimum_time = times[step_min];
          options.maximum_time = times[step_max];

          std::string filename = interFace.outputFile;
          if (interFace.split_cyclic > 0) {
            static const std::string suffix{"ABCDEFGHIJKLMNOPQRSTUVWXYZ"};
            filename += "." + suffix.substr(split % interFace.split_cyclic, 1);
          }
          else {
            filename = fmt::format("{0}_{1:0{2}}", filename, split + 1, width);
          }

          if (rank == 0 && !interFace.quiet) {
            if (step_min == step_max) {
              fmt::print(stderr, "\tWriting step {:n} to {}\n", step_min + 1, filename);
            }
            else {
              fmt::print(stderr, "\tWriting steps {:n}..{:n} to {}\n", step_min + 1, step_max + 1,
                         filename);
            }
          }

          Ioss::DatabaseIO *dbo =
              Ioss::IOFactory::create(interFace.outFiletype, filename, Ioss::WRITE_RESTART,
                                      (MPI_Comm)MPI_COMM_WORLD, properties);
          if (dbo == nullptr || !dbo->ok(true)) {
            std::exit(EXIT_FAILURE);
          }

          // NOTE: 'output_region' owns 'dbo' pointer at this time
          Ioss::Region output_region(dbo, "region_2");
          // Set the qa information...
          output_region.property_add(Ioss::Property(std::string("code_name"), codename));
          output_region.property_add(Ioss::Property(std::string("code_version"), version));

          Ioss::Utils::copy_database(region, output_region, options);
          if (mem_stats) {
            dbo->release_memory();
          }
          options.verbose = false;
        }
      }
      if (mem_stats) {
        dbi->progress("Prior to Memory Released... ");
        dbi->release_memory();
        dbi->progress("Memory Released... ");
      }
    } // loop over input files
  }

  Ioss::PropertyManager set_properties(IOShell::Interface &interFace)
  {
    Ioss::PropertyManager properties;

    if (interFace.ints_64_bit) {
      properties.add(Ioss::Property("INTEGER_SIZE_DB", 8));
      properties.add(Ioss::Property("INTEGER_SIZE_API", 8));
    }

    if (interFace.ints_32_bit) {
      properties.add(Ioss::Property("INTEGER_SIZE_DB", 4));
    }

    if (interFace.reals_32_bit) {
      properties.add(Ioss::Property("REAL_SIZE_DB", 4));
    }

    if (interFace.in_memory_read) {
      properties.add(Ioss::Property("MEMORY_READ", 1));
    }

    if (interFace.in_memory_write) {
      properties.add(Ioss::Property("MEMORY_WRITE", 1));
    }

    if (interFace.compression_level > 0 || interFace.shuffle) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
      properties.add(Ioss::Property("COMPRESSION_LEVEL", interFace.compression_level));
      properties.add(Ioss::Property("COMPRESSION_SHUFFLE", static_cast<int>(interFace.shuffle)));

      if (interFace.szip) {
        properties.add(Ioss::Property("COMPRESSION_METHOD", "szip"));
      }
      else if (interFace.zlib) {
        properties.add(Ioss::Property("COMPRESSION_METHOD", "zlib"));
      }
    }

    if (interFace.compose_output == "default") {
      if (interFace.outFiletype == "cgns") {
        properties.add(Ioss::Property("COMPOSE_RESULTS", "YES"));
        properties.add(Ioss::Property("COMPOSE_RESTART", "YES"));
      }
      else {
        properties.add(Ioss::Property("COMPOSE_RESULTS", "NO"));
        properties.add(Ioss::Property("COMPOSE_RESTART", "NO"));
      }
    }
    else if (interFace.compose_output == "external") {
      properties.add(Ioss::Property("COMPOSE_RESULTS", "NO"));
      properties.add(Ioss::Property("COMPOSE_RESTART", "NO"));
    }
    else if (interFace.compose_output != "none") {
      properties.add(Ioss::Property("COMPOSE_RESULTS", "YES"));
      properties.add(Ioss::Property("COMPOSE_RESTART", "YES"));
    }

    if (interFace.file_per_state) {
      properties.add(Ioss::Property("FILE_PER_STATE", "YES"));
    }

    if (interFace.netcdf4) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf4"));
    }

    if (interFace.netcdf5) {
      properties.add(Ioss::Property("FILE_TYPE", "netcdf5"));
    }

    if (interFace.inputFile.size() > 1) {
      properties.add(Ioss::Property("ENABLE_FILE_GROUPS", 1));
    }

    if (interFace.debug) {
      properties.add(Ioss::Property("LOGGING", 1));
    }

    if (interFace.memory_statistics) {
      properties.add(Ioss::Property("ENABLE_TRACING", 1));
    }

    if (!interFace.decomp_method.empty()) {
      properties.add(Ioss::Property("DECOMPOSITION_METHOD", interFace.decomp_method));
    }

    if (interFace.retain_empty_blocks) {
      properties.add(Ioss::Property("RETAIN_EMPTY_BLOCKS", "YES"));
    }
    return properties;
  }
} // namespace
