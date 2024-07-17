// Copyright(C) 1999-2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "Ionit_Initializer.h"
#include "Ioss_Compare.h"
#include "Ioss_CopyDatabase.h"
#include "Ioss_FileInfo.h"
#include "Ioss_MemoryUtils.h"
#include "Ioss_MeshCopyOptions.h"
#include "Ioss_MeshType.h"
#include "Ioss_ParallelUtils.h"
#include "Ioss_SerializeIO.h"
#include "Ioss_SurfaceSplit.h"
#include "Ioss_Utils.h"
#include <cstdlib>
#include <exception>
#include <fmt/core.h>
#include <fmt/format.h>
#include <limits>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <tokenize.h>
#include <vector>

#include "Ioss_DBUsage.h"
#include "Ioss_DataSize.h"
#include "Ioss_DatabaseIO.h"
#include "Ioss_GetLongOpt.h"
#include "Ioss_IOFactory.h"
#include "Ioss_Property.h"
#include "Ioss_PropertyManager.h"
#include "Ioss_Region.h"
#include "Ioss_ScopeGuard.h"
#include "Ioss_VariableType.h"
#include "shell_interface.h"

// ========================================================================

namespace {
  std::string codename;
  std::string version = "6.8 (2024/05/31)";

  bool mem_stats = false;

  void file_copy(IOShell::Interface &interFace, int rank);
  bool file_compare(IOShell::Interface &interFace, int rank);

  Ioss::PropertyManager set_properties(IOShell::Interface &interFace);
  Ioss::MeshCopyOptions set_mesh_copy_options(IOShell::Interface &interFace)
  {
    Ioss::MeshCopyOptions options{};
    options.selected_times    = interFace.selected_times;
    options.rel_tolerance     = interFace.rel_tolerance;
    options.abs_tolerance     = interFace.abs_tolerance;
    options.tol_floor         = interFace.tol_floor;
    options.verbose           = !interFace.quiet;
    options.output_summary    = true;
    options.memory_statistics = interFace.memory_statistics;
    options.debug             = interFace.debug;
    options.ints_64_bit       = interFace.ints_64_bit;
    options.delete_timesteps  = interFace.delete_timesteps;
    options.minimum_time      = interFace.minimum_time;
    options.maximum_time      = interFace.maximum_time;
    options.time_scale        = interFace.time_scale;
    options.time_offset       = interFace.time_offset;
    options.data_storage_type = interFace.data_storage_type;
    options.delay             = interFace.timestep_delay;
    options.reverse           = interFace.reverse;
    options.add_proc_id       = interFace.add_processor_id_field;
    options.boundary_sideset  = interFace.boundary_sideset;
    options.ignore_qa_info    = interFace.ignore_qa_info;
    options.omitted_blocks    = !interFace.omitted_blocks.empty();

    options.omitted_sets = interFace.omitted_sets;
    Ioss::sort(options.omitted_sets);
    for (auto &name : options.omitted_sets) {
      name = Ioss::Utils::lowercase(name);
    }
    return options;
  }

#ifdef SEACAS_HAVE_MPI
  void mpi_finalize()
  {
    MPI_Comm parentcomm;
    MPI_Comm_get_parent(&parentcomm);
    if (parentcomm != MPI_COMM_NULL) {
      int istatus = EXIT_SUCCESS;
      MPI_Send(&istatus, 1, MPI_INT, 0, 0, parentcomm);
    }
    MPI_Finalize();
  }
#endif
} // namespace

int main(int argc, char *argv[])
{
#ifdef SEACAS_HAVE_MPI
  MPI_Init(&argc, &argv);
  ON_BLOCK_EXIT(mpi_finalize);
#endif
  Ioss::ParallelUtils pu{};
  int                 rank     = pu.parallel_rank();
  int                 num_proc = pu.parallel_size();

#ifdef SEACAS_HAVE_KOKKOS
  Kokkos::ScopeGuard kokkos(argc, argv);
#endif

  IOShell::Interface interFace(version);
  bool               success = interFace.parse_options(argc, argv, rank);
  if (!success) {
    exit(EXIT_FAILURE);
  }

  codename = Ioss::GetLongOption::basename(argv[0]);

  Ioss::SerializeIO::setGroupFactor(interFace.serialize_io_size);
  mem_stats = interFace.memory_statistics;

  Ioss::Init::Initializer io;

  // See if a custom field is defined...
  if (!interFace.customField.empty()) {
    auto suffices = Ioss::tokenize(interFace.customField, ",");
    if (suffices.size() > 1) {
      Ioss::VariableType::create_named_suffix_type("UserDefined", suffices);
    }
  }
  std::string in_file  = interFace.inputFile[0];
  std::string out_file = interFace.outputFile;

  if (rank == 0 && !interFace.quiet) {
    if (interFace.compare) {
      fmt::print(stderr,
                 "Input 1:   '{}', Type: {}\n"
                 "Input 2:   '{}', Type: {}\n"
                 "\tTolerances: Absolute = {}, Relative = {}, Floor = {}\n\n",
                 in_file, interFace.inFiletype, out_file, interFace.outFiletype,
                 interFace.abs_tolerance, interFace.rel_tolerance, interFace.tol_floor);
    }
    else {
      fmt::print(stderr,
                 "Input:    '{}', Type: {}\n"
                 "Output:   '{}', Type: {}\n\n",
                 in_file, interFace.inFiletype, out_file, interFace.outFiletype);
    }
  }

#ifdef SEACAS_HAVE_KOKKOS
  if (rank == 0)
    fmt::print(stderr, "Kokkos default execution space configuration:\n");
  Kokkos::DefaultExecutionSpace().print_configuration(std::cerr, false);
  if (rank == 0)
    fmt::print(stderr, "\n");
#endif

  double begin = Ioss::Utils::timer();

  try {
    if (interFace.compare) {
      success = file_compare(interFace, rank);
    }
    else {
      file_copy(interFace, rank);
    }
  }
  catch (std::exception &e) {
    if (rank == 0) {
      fmt::print(stderr, "\n{}\n\nio_shell terminated due to exception\n", e.what());
    }
    exit(EXIT_FAILURE);
  }

  pu.barrier();
  double end = Ioss::Utils::timer();

  if (rank == 0 && !interFace.quiet) {
    if (num_proc > 1) {
      fmt::print(stderr, "\n\n\tTotal Execution Time = {:.5} seconds on {} processors.\n",
                 end - begin, num_proc);
    }
    else {
      fmt::print(stderr, "\n\n\tTotal Execution Time = {:.5} seconds.\n", end - begin);
    }
  }
  if (mem_stats) {
    int64_t MiB = 1024 * 1024;
#ifdef SEACAS_HAVE_MPI
    int64_t min, max, avg;
    int64_t hwmin, hwmax, hwavg;
    pu.memory_stats(min, max, avg);
    pu.hwm_memory_stats(hwmin, hwmax, hwavg);
    if (rank == 0) {
      fmt::print(stderr, "\n\tCurrent Memory: {}M  {}M  {}M\n", fmt::group_digits(min / MiB),
                 fmt::group_digits(max / MiB), fmt::group_digits(avg / MiB));
      fmt::print(stderr, "\tHigh Water Memory: {}M  {}M  {}M\n", fmt::group_digits(hwmin / MiB),
                 fmt::group_digits(hwmax / MiB), fmt::group_digits(hwavg / MiB));
    }
#else
    int64_t mem = Ioss::MemoryUtils::get_memory_info();
    int64_t hwm = Ioss::MemoryUtils::get_hwm_memory_info();
    if (rank == 0) {
      fmt::print(stderr,
                 "\n\tCurrent Memory:    {}M\n"
                 "\tHigh Water Memory: {}M\n",
                 fmt::group_digits(mem / MiB), fmt::group_digits(hwm / MiB));
    }
#endif
  }
  if (rank == 0) {
    fmt::print(stderr, "\n{} execution successful.\n", codename);
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
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
      Ioss::DatabaseIO *dbi =
          Ioss::IOFactory::create(interFace.inFiletype, inpfile, Ioss::READ_MODEL,
                                  Ioss::ParallelUtils::comm_world(), properties);
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
      else if (interFace.surface_split_type != Ioss::SPLIT_INVALID) {
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

      if (!interFace.omitted_blocks.empty()) {
        std::vector<std::string> inclusions{};
        dbi->set_block_omissions(interFace.omitted_blocks, inclusions);
      }

      // NOTE: 'region' owns 'db' pointer at this time...
      Ioss::Region region(dbi, "region_1");

      if (region.mesh_type() == Ioss::MeshType::HYBRID) {
        if (rank == 0) {
          fmt::print(stderr,
                     "\nERROR: io_shell does not support '{}' meshes. Only 'Unstructured' or "
                     "'Structured' mesh is supported at this time.\n",
                     region.mesh_type_string());
        }
        return;
      }

      // Get length of longest name on input file...
      int max_name_length = dbi->maximum_symbol_length();
      if (max_name_length > 0) {
        properties.add(Ioss::Property("MAXIMUM_NAME_LENGTH", max_name_length));
      }

      // Get integer size being used on the input file and propagate
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

      Ioss::MeshCopyOptions options = set_mesh_copy_options(interFace);

      size_t ts_count = region.get_optional_property("state_count", 0);

      int flush_interval = interFace.flush_interval; // Default is zero -- do not flush until end
      properties.add(Ioss::Property("FLUSH_INTERVAL", flush_interval));

      if (interFace.split_times == 0 || interFace.delete_timesteps || ts_count == 0 || append ||
          interFace.inputFile.size() > 1) {
        Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(
            interFace.outFiletype, interFace.outputFile, Ioss::WRITE_RESTART,
            Ioss::ParallelUtils::comm_world(), properties);
        if (dbo == nullptr || !dbo->ok(true)) {
          std::exit(EXIT_FAILURE);
        }

        dbo->set_field_separator(interFace.fieldSuffixSeparator);

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
        Ioss::copy_database(region, output_region, options);

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
              fmt::print(stderr, "\tWriting step {} to {}\n", fmt::group_digits(step_min + 1),
                         filename);
            }
            else {
              fmt::print(stderr, "\tWriting steps {}..{} to {}\n", fmt::group_digits(step_min + 1),
                         fmt::group_digits(step_max + 1), filename);
            }
          }

          Ioss::DatabaseIO *dbo =
              Ioss::IOFactory::create(interFace.outFiletype, filename, Ioss::WRITE_RESTART,
                                      Ioss::ParallelUtils::comm_world(), properties);
          if (dbo == nullptr || !dbo->ok(true)) {
            std::exit(EXIT_FAILURE);
          }

          dbo->set_field_separator(interFace.fieldSuffixSeparator);

          // NOTE: 'output_region' owns 'dbo' pointer at this time
          Ioss::Region output_region(dbo, "region_2");
          // Set the qa information...
          output_region.property_add(Ioss::Property(std::string("code_name"), codename));
          output_region.property_add(Ioss::Property(std::string("code_version"), version));

          Ioss::copy_database(region, output_region, options);
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

  bool file_compare(IOShell::Interface &interFace, int rank)
  {
    Ioss::PropertyManager properties = set_properties(interFace);
    const auto           &inpfile    = interFace.inputFile[0];

    //========================================================================
    // INPUT Database #1...
    //========================================================================
    Ioss::DatabaseIO *dbi1 =
        Ioss::IOFactory::create(interFace.inFiletype, inpfile, Ioss::READ_MODEL,
                                Ioss::ParallelUtils::comm_world(), properties);
    if (dbi1 == nullptr || !dbi1->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    if (mem_stats) {
      dbi1->progress("Database #1 Open");
    }
    if (!interFace.lower_case_variable_names) {
      dbi1->set_lower_case_variable_names(false);
    }
    if (interFace.outFiletype == "cgns") {
      // CGNS stores BCs (SideSets) on the zones which
      // correspond to element blocks.  If split input sideblocks
      // by element block, then output is much easier.
      dbi1->set_surface_split_type(Ioss::SPLIT_BY_ELEMENT_BLOCK);
    }
    else if (interFace.surface_split_type != Ioss::SPLIT_INVALID) {
      dbi1->set_surface_split_type(Ioss::int_to_surface_split(interFace.surface_split_type));
    }
    dbi1->set_field_separator(interFace.fieldSuffixSeparator);

    dbi1->set_field_recognition(!interFace.disable_field_recognition);

    if (interFace.ints_64_bit) {
      dbi1->set_int_byte_size_api(Ioss::USE_INT64_API);
    }

    if (!interFace.groupName.empty()) {
      bool success = dbi1->open_group(interFace.groupName);
      if (!success) {
        if (rank == 0) {
          fmt::print(stderr, "ERROR: Unable to open group '{}' in file '{}'\n", interFace.groupName,
                     inpfile);
        }
        return false;
      }
    }

    // NOTE: 'input_region1' owns 'dbi1' pointer at this time...
    Ioss::Region input_region1(dbi1, "region_1");

    if (input_region1.mesh_type() == Ioss::MeshType::HYBRID) {
      fmt::print(stderr,
                 "\nERROR: io_shell does not support '{}' meshes. Only 'Unstructured' or "
                 "'Structured' mesh is supported at this time.\n",
                 input_region1.mesh_type_string());
      return false;
    }

    // Get integer size being used on input file #1 and set it in
    // the interFace.
    int int_byte_size_api = dbi1->int_byte_size_api();
    if (int_byte_size_api == 8) {
      interFace.ints_64_bit = true;
    }

    //========================================================================
    // INPUT Database #2...
    //========================================================================
    Ioss::DatabaseIO *dbi2 =
        Ioss::IOFactory::create(interFace.outFiletype, interFace.outputFile, Ioss::READ_MODEL,
                                Ioss::ParallelUtils::comm_world(), properties);
    if (dbi2 == nullptr || !dbi2->ok(true)) {
      std::exit(EXIT_FAILURE);
    }

    if (mem_stats) {
      dbi2->progress("Database #2 Open");
    }
    if (!interFace.lower_case_variable_names) {
      dbi2->set_lower_case_variable_names(false);
    }
    if (interFace.outFiletype == "cgns") {
      // CGNS stores BCs (SideSets) on the zones which
      // correspond to element blocks.  If split input sideblocks
      // by element block, then output is much easier.
      dbi2->set_surface_split_type(Ioss::SPLIT_BY_ELEMENT_BLOCK);
    }
    else if (interFace.surface_split_type != Ioss::SPLIT_INVALID) {
      dbi2->set_surface_split_type(Ioss::int_to_surface_split(interFace.surface_split_type));
    }
    dbi2->set_field_separator(interFace.fieldSuffixSeparator);

    dbi2->set_field_recognition(!interFace.disable_field_recognition);

    if (interFace.ints_64_bit) {
      dbi2->set_int_byte_size_api(Ioss::USE_INT64_API);
    }

    if (!interFace.groupName.empty()) {
      bool success = dbi2->open_group(interFace.groupName);
      if (!success) {
        if (rank == 0) {
          fmt::print(stderr, "ERROR: Unable to open group '{}' in file '{}'\n", interFace.groupName,
                     inpfile);
        }
        return false;
      }
    }

    // NOTE: 'input_region2' owns 'dbi2' pointer at this time...
    Ioss::Region input_region2(dbi2, "region_2");

    if (input_region2.mesh_type() == Ioss::MeshType::HYBRID) {
      fmt::print(stderr,
                 "\nERROR: io_shell does not support '{}' meshes. Only 'Unstructured' or "
                 "'Structured' mesh is supported at this time.\n",
                 input_region2.mesh_type_string());
      return false;
    }

    // Get integer size being used on input file #1 and set it in
    // the interFace.
    int_byte_size_api = dbi2->int_byte_size_api();
    if (int_byte_size_api == 8) {
      interFace.ints_64_bit = true;
    }

    //========================================================================
    // COMPARE the databases...
    //========================================================================
    auto options = set_mesh_copy_options(interFace);

    bool result = Ioss::Compare::compare_database(input_region1, input_region2, options);
    if (result) {
      fmt::print(stderr, "\n\nDATABASES are EQUAL");
    }
    else {
      fmt::print(stderr, "\n\nDATABASES are NOT equal");
    }
    return result;
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

    if (interFace.delete_qa) {
      properties.add(Ioss::Property("IGNORE_QA_RECORDS", "YES"));
    }
    if (interFace.delete_info) {
      properties.add(Ioss::Property("IGNORE_INFO_RECORDS", "YES"));
    }

    if (interFace.compression_level > 0 || interFace.shuffle || interFace.szip) {
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

    if (interFace.detect_nans) {
      properties.add(Ioss::Property("NAN_DETECTION", 1));
    }

    if (interFace.memory_statistics) {
      properties.add(Ioss::Property("ENABLE_TRACING", 1));
    }

    if (interFace.outFiletype == "cgns" && interFace.inFiletype == "exodus") {
      properties.add(Ioss::Property("IGNORE_NODE_MAP", true));
      properties.add(Ioss::Property("IGNORE_ELEMENT_MAP", true));
    }
    else {
      if (interFace.ignore_node_map) {
        properties.add(Ioss::Property("IGNORE_NODE_MAP", true));
      }
      if (interFace.ignore_elem_map) {
        properties.add(Ioss::Property("IGNORE_ELEM_MAP", true));
      }
    }
    if (interFace.ignore_edge_map) {
      properties.add(Ioss::Property("IGNORE_EDGE_MAP", true));
    }
    if (interFace.ignore_face_map) {
      properties.add(Ioss::Property("IGNORE_FACE_MAP", true));
    }

    if (!interFace.decomp_method.empty()) {
      properties.add(Ioss::Property("DECOMPOSITION_METHOD", interFace.decomp_method));
      if (interFace.decomp_method == "MAP" || interFace.decomp_method == "VARIABLE") {
        properties.add(Ioss::Property("DECOMPOSITION_EXTRA", interFace.decomp_extra));
      }
      if (interFace.line_decomp) {
        properties.add(Ioss::Property("LINE_DECOMPOSITION", interFace.decomp_extra));
      }
    }

    if (interFace.retain_empty_blocks) {
      properties.add(Ioss::Property("RETAIN_EMPTY_BLOCKS", "YES"));
    }
    return properties;
  }
} // namespace
