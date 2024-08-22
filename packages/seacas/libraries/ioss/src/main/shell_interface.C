/*
 * Copyright(C) 1999-2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#include <cstdlib> // for exit, strtod, EXIT_SUCCESS, etc
#include <cstring> // for strcmp
#include <fmt/core.h>
#include <iostream> // for operator<<, basic_ostream, etc
#include <stdio.h>
#include <string> // for string, char_traits
#include <vector> // for vector

#include "Ioss_GetLongOpt.h" // for GetLongOption, etc
#include "Ioss_Sort.h"
#include "Ioss_Utils.h" // for Utils
#include "shell_interface.h"
#include "tokenize.h"

IOShell::Interface::Interface(std::string app_version) : version(std::move(app_version))
{
  enroll_options();
}

void IOShell::Interface::enroll_options()
{
  options_.usage("[options] input_file[s] output_file");

  options_.enroll("help", Ioss::GetLongOption::NoValue, "Print this summary and exit", nullptr);

  options_.enroll("in_type", Ioss::GetLongOption::MandatoryValue,
                  "Database type for input file: generated"
#if defined(SEACAS_HAVE_PAMGEN)
                  "|pamgen"
#endif
#if defined(SEACAS_HAVE_EXODUS)
                  "|exodus"
#endif
#if defined(SEACAS_HAVE_CGNS)
                  "|cgns"
#endif
                  ".\n\t\tIf not specified, guess from extension or exodus is the default.",
                  "unknown");

  options_.enroll("out_type", Ioss::GetLongOption::MandatoryValue,
                  "Database type for output file:"
#if defined(SEACAS_HAVE_EXODUS)
                  " exodus"
#if defined(SEACAS_HAVE_EXONULL)
                  " exonull"
#endif
#endif
#if defined(SEACAS_HAVE_CGNS)
                  " cgns"
#endif
#if defined(SEACAS_HAVE_FAODEL)
                  " faodel"
#endif
                  " null.\n\t\tIf not specified, guess from extension or exodus is the default.",
                  "unknown");
  options_.enroll("compare", Ioss::GetLongOption::NoValue,
                  "Compare the contents of the INPUT and OUTPUT files.", nullptr);
  options_.enroll(
      "relative", Ioss::GetLongOption::MandatoryValue,
      "Relative tolerance to use if comparing real field data. (diff > abs && diff > rel).",
      nullptr);
  options_.enroll(
      "absolute", Ioss::GetLongOption::MandatoryValue,
      "Absolute tolerance to use if comparing real field data. (diff > abs && diff > rel)",
      nullptr);
  options_.enroll("floor", Ioss::GetLongOption::MandatoryValue,
                  "Only compare values if `|a| > floor || |b| > floor`", nullptr);
  options_.enroll("ignore_qa_info", Ioss::GetLongOption::NoValue,
                  "If comparing databases, do not compare the qa and info records.", nullptr,
                  nullptr, true);

  options_.enroll("ignore_node_map", Ioss::GetLongOption::NoValue,
                  "Do not read the global node id map (if any) from the input database.", nullptr,
                  nullptr);
  options_.enroll("ignore_element_map", Ioss::GetLongOption::NoValue,
                  "Do not read the global element id map (if any) from the input database.",
                  nullptr, nullptr);
  options_.enroll("ignore_edge_map", Ioss::GetLongOption::NoValue,
                  "Do not read the global edge id map (if any) from the input database.", nullptr,
                  nullptr);
  options_.enroll("ignore_face_map", Ioss::GetLongOption::NoValue,
                  "Do not read the global face id map (if any) from the input database.", nullptr,
                  nullptr, true);

  options_.enroll("64-bit", Ioss::GetLongOption::NoValue, "Use 64-bit integers on output database",
                  nullptr);

  options_.enroll("32-bit", Ioss::GetLongOption::NoValue,
                  "Use 32-bit integers on output database."
                  " This is the default unless input database uses 64-bit integers",
                  nullptr);

  options_.enroll("float", Ioss::GetLongOption::NoValue,
                  "Use 32-bit floating point values on output database; default is 64-bits",
                  nullptr);

  options_.enroll("netcdf3", Ioss::GetLongOption::NoValue,
                  "Output database will be a classical netcdf (CDF3) file.", nullptr);

  options_.enroll("netcdf4", Ioss::GetLongOption::NoValue,
                  "Output database will be a netcdf4 "
                  "hdf5-based file instead of the "
                  "classical netcdf file format",
                  nullptr);

  options_.enroll("netcdf5", Ioss::GetLongOption::NoValue,
                  "Output database will be a netcdf5 (CDF5) "
                  "file instead of the classical netcdf file format",
                  nullptr);

  options_.enroll("shuffle", Ioss::GetLongOption::NoValue,
                  "Use a netcdf4 hdf5-based file and use hdf5s shuffle mode with compression.",
                  nullptr);

  options_.enroll("compress", Ioss::GetLongOption::MandatoryValue,
                  "Specify the hdf5 zlib compression level [0..9] or szip [even, 4..32] to be used "
                  "on the output file.",
                  nullptr);

  options_.enroll(
      "zlib", Ioss::GetLongOption::NoValue,
      "Use the Zlib / libz compression method if compression is enabled (default) [exodus only].",
      nullptr);

  options_.enroll(
      "szip", Ioss::GetLongOption::NoValue,
      "Use the SZip library if compression is enabled. Not as portable as zlib [exodus only]",
      nullptr, nullptr, true);

#if defined(SEACAS_HAVE_MPI)
  options_.enroll(
      "compose", Ioss::GetLongOption::OptionalValue,
      "If no argument, specify single-file output; if 'external', then file-per-processor.\n"
      "\t\tAll other options are ignored and just exist for backward-compatibility",
      nullptr, "true");

#if !defined(NO_ZOLTAN_SUPPORT)
  options_.enroll(
      "rcb", Ioss::GetLongOption::NoValue,
      "Use recursive coordinate bisection method to decompose the input mesh in a parallel run.",
      nullptr);
  options_.enroll(
      "rib", Ioss::GetLongOption::NoValue,
      "Use recursive inertial bisection method to decompose the input mesh in a parallel run.",
      nullptr);

  options_.enroll(
      "hsfc", Ioss::GetLongOption::NoValue,
      "Use hilbert space-filling curve method to decompose the input mesh in a parallel run.",
      nullptr);
#endif

#if !defined(NO_PARMETIS_SUPPORT)
  options_.enroll(
      "metis_sfc", Ioss::GetLongOption::NoValue,
      "Use the metis space-filling-curve method to decompose the input mesh in a parallel run.",
      nullptr);

  options_.enroll(
      "kway", Ioss::GetLongOption::NoValue,
      "Use the metis kway graph-based method to decompose the input mesh in a parallel run.",
      nullptr);

  options_.enroll("kway_geom", Ioss::GetLongOption::NoValue,
                  "Use the metis kway graph-based method with geometry speedup to decompose the "
                  "input mesh in a parallel run.",
                  nullptr);
#endif

  options_.enroll("linear", Ioss::GetLongOption::NoValue,
                  "Use the linear method to decompose the input mesh in a parallel run.\n"
                  "\t\tElements in order first n/p to proc 0, next to proc 1.",
                  nullptr);

  options_.enroll("cyclic", Ioss::GetLongOption::NoValue,
                  "Use the cyclic method to decompose the input mesh in a parallel run.\n"
                  "\t\tElements handed out to id % proc_count",
                  nullptr);

  options_.enroll("random", Ioss::GetLongOption::NoValue,
                  "Use the random method to decompose the input mesh in a parallel run.\n"
                  "\t\tElements assigned randomly to processors in a way that preserves balance\n"
                  "\t\t(do *not* use for a real run)",
                  nullptr);

  options_.enroll("map", Ioss::GetLongOption::OptionalValue,
                  "Read the decomposition data from the specified element map.\n"
                  "\t\tIf no map name is specified, then `processor_id` will be used.\n"
                  "\t\tIf the name is followed by a ',' and an integer or 'auto', then\n"
                  "\t\tthe entries in the map will be divided by the integer value or\n"
                  "\t\t(if auto) by `int((max_entry+1)/proc_count)`.",
                  nullptr);

  options_.enroll("variable", Ioss::GetLongOption::OptionalValue,
                  "Read the decomposition data from the specified element variable.\n"
                  "\t\tIf no variable name is specified, then `processor_id` will be used.\n"
                  "\t\tIf the name is followed by a ',' and an integer or 'auto', then\n"
                  "\t\tthe entries in the variable will be divided by the integer value or\n"
                  "\t\t(if auto) by `int((max_entry+1)/proc_count)`.",
                  nullptr);

  options_.enroll("line_decomp", Ioss::GetLongOption::OptionalValue,
                  "Generate the `lines` or `columns` of elements from the specified surface(s).\n"
                  "\t\tSpecify a comma-separated list of surface/sideset names from which the "
                  "lines will grow.\n"
                  "\t\tDo not split a line/column across processors.\n"
                  "\t\tOmit or enter 'ALL' for all surfaces in model.",
                  nullptr, "ALL");

  options_.enroll("external", Ioss::GetLongOption::NoValue,
                  "Files are decomposed externally into a file-per-processor in a parallel run.",
                  nullptr);

  options_.enroll("add_processor_id_field", Ioss::GetLongOption::NoValue,
                  "Add a cell-centered field whose value is the processor id of that cell",
                  nullptr);

  options_.enroll("serialize_io_size", Ioss::GetLongOption::MandatoryValue,
                  "Number of processors that can perform simultaneous IO operations in "
                  "a parallel run;\n\t\t0 to disable",
                  nullptr, nullptr, true);
#endif

  options_.enroll("extract_group", Ioss::GetLongOption::MandatoryValue,
                  "Write the data from the specified group to the output file.", nullptr);

  options_.enroll(
      "split_times", Ioss::GetLongOption::MandatoryValue,
      "If non-zero, then put <$val> timesteps in each file. Then close file and start new file.",
      nullptr);

  options_.enroll(
      "split_cyclic", Ioss::GetLongOption::MandatoryValue,
      "If non-zero, then the `split_times` timesteps will be put into <$val> files\n\t\tand "
      "then recycle filenames.",
      nullptr);

  options_.enroll("file_per_state", Ioss::GetLongOption::NoValue,
                  "put transient data for each timestep in separate file (EXPERIMENTAL)", nullptr);

  options_.enroll("minimize_open_files", Ioss::GetLongOption::NoValue,
                  "close output file after each timestep", nullptr, nullptr, true);

  options_.enroll("Maximum_Time", Ioss::GetLongOption::MandatoryValue,
                  "Maximum time on input database to transfer to output database", nullptr);

  options_.enroll("Minimum_Time", Ioss::GetLongOption::MandatoryValue,
                  "Minimum time on input database to transfer to output database", nullptr);

  options_.enroll("time_scale", Ioss::GetLongOption::MandatoryValue,
                  "The output time = input_time * time_scale + time_offset", nullptr);

  options_.enroll("time_offset", Ioss::GetLongOption::MandatoryValue,
                  "The output time = input_time * time_scale + time_offset", nullptr);

  options_.enroll("select_times", Ioss::GetLongOption::MandatoryValue,
                  "comma-separated list of times that should be transferred to output database",
                  nullptr);

  options_.enroll("append_after_time", Ioss::GetLongOption::MandatoryValue,
                  "add steps on input database after specified time on output database", nullptr);

  options_.enroll("append_after_step", Ioss::GetLongOption::MandatoryValue,
                  "add steps on input database after specified step on output database", nullptr);

  options_.enroll("flush_interval", Ioss::GetLongOption::MandatoryValue,
                  "Specify the number of steps between database flushes.\n"
                  "\t\tIf not specified, then the default database-dependent setting is used.\n"
                  "\t\tA value of 0 disables flushing.",
                  nullptr);
  options_.enroll("delete_timesteps", Ioss::GetLongOption::NoValue,
                  "Do not transfer any timesteps or transient data to the output database",
                  nullptr);

  options_.enroll("delete_qa_records", Ioss::GetLongOption::NoValue,
                  "Do not output qa records to output database.", nullptr);
  options_.enroll("delete_info_records", Ioss::GetLongOption::NoValue,
                  "Do not output info records to output database.", nullptr, nullptr, true);

  options_.enroll("field_suffix_separator", Ioss::GetLongOption::MandatoryValue,
                  "Character used to separate a field suffix from the field basename\n"
                  "\t\twhen recognizing vector, tensor fields. Enter '0' for no separator",
                  "_");

  options_.enroll("disable_field_recognition", Ioss::GetLongOption::NoValue,
                  "Do not combine fields into vector, tensor fields based on basename and suffix.\n"
                  "\t\tKeep all fields on database as scalars",
                  nullptr);

  options_.enroll("custom_field", Ioss::GetLongOption::MandatoryValue,
                  "A comma-separated list of field suffices defining a custom field that should be "
                  "recognized.\n"
                  "\t\tPrimarily used for testing",
                  nullptr);

  options_.enroll("surface_split_scheme", Ioss::GetLongOption::MandatoryValue,
                  "Method used to split sidesets into homogeneous blocks\n"
                  "\t\tOptions are: TOPOLOGY(default), BLOCK, NO_SPLIT",
                  nullptr);

  options_.enroll("native_variable_names", Ioss::GetLongOption::NoValue,
                  "Do not lowercase variable names and replace spaces with underscores.\n"
                  "\t\tVariable names are left as they appear in the input mesh file",
                  nullptr);

  options_.enroll("retain_empty_blocks", Ioss::GetLongOption::NoValue,
                  "If any empty element blocks on input file, keep them and write to output file.\n"
                  "\t\tDefault is to ignore empty blocks.",
                  nullptr);

  options_.enroll("omit_blocks", Ioss::GetLongOption::MandatoryValue,
                  "comma-separated list of element block names that should NOT be transferred to "
                  "output database\n"
                  "\t\tNote that currently any nodes connected to only empty blocks will be "
                  "retained in the output.",
                  nullptr);

  options_.enroll("omit_sets", Ioss::GetLongOption::MandatoryValue,
                  "comma-separated list of nodeset/edgeset/faceset/elemset/sideset names\n"
                  "\t\tthat should NOT be transferred to output database",
                  nullptr);

  options_.enroll("boundary_sideset", Ioss::GetLongOption::NoValue,
                  "Output a sideset for all boundary faces of the model", nullptr, nullptr, true);

  options_.enroll(
      "delay", Ioss::GetLongOption::MandatoryValue,
      "Sleep for <$val> seconds between timestep output to simulate application calculation time",
      nullptr);

#ifdef SEACAS_HAVE_KOKKOS
  options_.enroll("data_storage", Ioss::GetLongOption::MandatoryValue,
                  "Data type used internally to store field data\n"
                  "\t\tOptions are: POINTER, STD_VECTOR, KOKKOS_VIEW_1D, KOKKOS_VIEW_2D, "
                  "KOKKOS_VIEW_2D_LAYOUTRIGHT_HOSTSPACE",
                  "POINTER");
#else
  options_.enroll("data_storage", Ioss::GetLongOption::MandatoryValue,
                  "Data type used internally to store field data\n"
                  "\t\tOptions are: POINTER, STD_VECTOR",
                  "POINTER");
#endif

  options_.enroll("debug", Ioss::GetLongOption::NoValue, "turn on debugging output", nullptr);
  options_.enroll("detect_nans", Ioss::GetLongOption::NoValue, "check all real field data for NaNs",
                  nullptr);

  options_.enroll("quiet", Ioss::GetLongOption::NoValue, "minimize output", nullptr);

  options_.enroll("statistics", Ioss::GetLongOption::NoValue,
                  "output parallel io timing statistics", nullptr);

  options_.enroll("memory_statistics", Ioss::GetLongOption::NoValue,
                  "output memory usage throughout code execution", nullptr);

  options_.enroll(
      "memory_read", Ioss::GetLongOption::NoValue,
      "EXPERIMENTAL: file read into memory by netcdf library; ioss accesses memory version",
      nullptr);

  options_.enroll(
      "memory_write", Ioss::GetLongOption::NoValue,
      "EXPERIMENTAL: file written to memory, netcdf library streams to disk at file close",
      nullptr);

  options_.enroll("reverse", Ioss::GetLongOption::NoValue,
                  "define CGNS zones in reverse order. Used for testing (TEST)", nullptr);

  options_.enroll("version", Ioss::GetLongOption::NoValue, "Print version and exit", nullptr);

  options_.enroll("copyright", Ioss::GetLongOption::NoValue, "Show copyright and license data.",
                  nullptr);
}

bool IOShell::Interface::parse_options(int argc, char **argv, int my_processor)
{
  // Get options from environment variable also...
  char *options = getenv("IO_SHELL_OPTIONS");
  if (options != nullptr) {
    fmt::print(
        stderr,
        "\nThe following options were specified via the IO_SHELL_OPTIONS environment variable:\n"
        "\t{}\n\n",
        options);
    options_.parse(options, Ioss::GetLongOption::basename(*argv));
  }

  int option_index = options_.parse(argc, argv);
  if (option_index < 1) {
    return false;
  }

  if (options_.retrieve("help") != nullptr) {
    if (my_processor == 0) {
      options_.usage(std::cerr);
      fmt::print(stderr, "\n\tCan also set options via IO_SHELL_OPTIONS environment variable.\n\n");
      fmt::print(stderr,
                 "\tDocumentation: "
                 "https://sandialabs.github.io/seacas-docs/sphinx/html/index.html#io-shell\n\n");
      fmt::print(stderr, "\t->->-> Send email to gdsjaar@sandia.gov for {} support.<-<-<-\n",
                 options_.program_name());
    }
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    fmt::print(stderr, "IO_SHELL\tVersion: {}\n", version);
    exit(0);
  }

  ints_64_bit  = (options_.retrieve("64-bit") != nullptr);
  ints_32_bit  = (options_.retrieve("32-bit") != nullptr);
  reals_32_bit = (options_.retrieve("float") != nullptr);

  if (options_.retrieve("netcdf3") != nullptr) {
    netcdf3     = true;
    netcdf4     = false;
    netcdf5     = false;
    ints_32_bit = true;
  }

  if (options_.retrieve("netcdf4") != nullptr) {
    netcdf3 = false;
    netcdf4 = true;
    netcdf5 = false;
  }

  if (options_.retrieve("netcdf5") != nullptr) {
    netcdf3 = false;
    netcdf4 = false;
    netcdf5 = true;
  }

  shuffle = (options_.retrieve("shuffle") != nullptr);
  if (options_.retrieve("szip") != nullptr) {
    szip = true;
    zlib = false;
  }
  zlib = (options_.retrieve("zlib") != nullptr);

  if (szip && zlib) {
    if (my_processor == 0) {
      fmt::print(stderr, "ERROR: Only one of 'szip' or 'zlib' can be specified.\n");
    }
    return false;
  }
  compare         = (options_.retrieve("compare") != nullptr);
  ignore_qa_info  = (options_.retrieve("ignore_qa_info") != nullptr);
  ignore_node_map = (options_.retrieve("ignore_node_map") != nullptr);
  ignore_elem_map = (options_.retrieve("ignore_element_map") != nullptr);
  ignore_edge_map = (options_.retrieve("ignore_edge_map") != nullptr);
  ignore_face_map = (options_.retrieve("ignore_face_map") != nullptr);
  delete_qa       = (options_.retrieve("delete_qa_records") != nullptr);
  delete_info     = (options_.retrieve("delete_info_records") != nullptr);

  {
    const char *temp = options_.retrieve("absolute");
    if (temp != nullptr) {
      abs_tolerance = std::strtod(temp, nullptr);
    }
  }
  {
    const char *temp = options_.retrieve("relative");
    if (temp != nullptr) {
      rel_tolerance = std::strtod(temp, nullptr);
    }
  }
  {
    const char *temp = options_.retrieve("floor");
    if (temp != nullptr) {
      tol_floor = std::strtod(temp, nullptr);
    }
  }

  {
    const char *temp = options_.retrieve("compress");
    if (temp != nullptr) {
      compression_level = std::strtol(temp, nullptr, 10);

      if (zlib) {
        if (compression_level < 0 || compression_level > 9) {
          if (my_processor == 0) {
            fmt::print(stderr,
                       "ERROR: Bad compression level {}, valid value is between 0 and 9 inclusive "
                       "for gzip compression.\n",
                       compression_level);
          }
          return false;
        }
      }
      else if (szip) {
        if (compression_level % 2 != 0) {
          if (my_processor == 0) {
            fmt::print(
                stderr,
                "ERROR: Bad compression level {}. Must be an even value for szip compression.\n",
                compression_level);
          }
          return false;
        }
        if (compression_level < 4 || compression_level > 32) {
          if (my_processor == 0) {
            fmt::print(stderr,
                       "ERROR: Bad compression level {}, valid value is between 4 and 32 inclusive "
                       "for szip compression.\n",
                       compression_level);
          }
          return false;
        }
      }
    }
  }

#if defined(SEACAS_HAVE_MPI)
  add_processor_id_field = (options_.retrieve("add_processor_id_field") != nullptr);

#if !defined(NO_ZOLTAN_SUPPORT)
  if (options_.retrieve("rcb") != nullptr) {
    decomp_method = "RCB";
  }

  if (options_.retrieve("rib") != nullptr) {
    decomp_method = "RIB";
  }

  if (options_.retrieve("hsfc") != nullptr) {
    decomp_method = "HSFC";
  }
#endif

#if !defined(NO_PARMETIS_SUPPORT)
  if (options_.retrieve("metis_sfc") != nullptr) {
    decomp_method = "METIS_SFC";
  }

  if (options_.retrieve("kway") != nullptr) {
    decomp_method = "KWAY";
  }

  if (options_.retrieve("kway_geom") != nullptr) {
    decomp_method = "KWAY_GEOM";
  }
#endif

  if (options_.retrieve("linear") != nullptr) {
    decomp_method = "LINEAR";
  }

  if (options_.retrieve("map") != nullptr) {
    decomp_method = "MAP";
    decomp_extra  = options_.get_option_value("map", decomp_extra);
  }

  if (options_.retrieve("variable") != nullptr) {
    decomp_method = "VARIABLE";
    decomp_extra  = options_.get_option_value("variable", decomp_extra);
  }

  if (options_.retrieve("line_decomp") != nullptr) {
    line_decomp  = true;
    decomp_extra = options_.get_option_value("line_decomp", decomp_extra);
  }

  if (options_.retrieve("cyclic") != nullptr) {
    decomp_method = "CYCLIC";
  }

  if (options_.retrieve("random") != nullptr) {
    decomp_method = "RANDOM";
  }

  if (options_.retrieve("external") != nullptr) {
    decomp_method = "EXTERNAL";
  }

  serialize_io_size = options_.get_option_value("serialize_io_size", serialize_io_size);

#endif

  split_times  = options_.get_option_value("split_times", split_times);
  split_cyclic = options_.get_option_value("split_cyclic", split_cyclic);
  if (split_cyclic > 26) {
    split_cyclic = 26;
  }

  minimize_open_files       = (options_.retrieve("minimize_open_files") != nullptr);
  debug                     = (options_.retrieve("debug") != nullptr);
  detect_nans               = (options_.retrieve("detect_nans") != nullptr);
  file_per_state            = (options_.retrieve("file_per_state") != nullptr);
  reverse                   = (options_.retrieve("reverse") != nullptr);
  quiet                     = (options_.retrieve("quiet") != nullptr);
  statistics                = (options_.retrieve("statistics") != nullptr);
  memory_statistics         = (options_.retrieve("memory_statistics") != nullptr);
  in_memory_read            = (options_.retrieve("memory_read") != nullptr);
  in_memory_write           = (options_.retrieve("memory_write") != nullptr);
  delete_timesteps          = (options_.retrieve("delete_timesteps") != nullptr);
  lower_case_variable_names = (options_.retrieve("native_variable_names") == nullptr);
  disable_field_recognition = (options_.retrieve("disable_field_recognition") != nullptr);
  retain_empty_blocks       = (options_.retrieve("retain_empty_blocks") != nullptr);
  boundary_sideset          = (options_.retrieve("boundary_sideset") != nullptr);

  inFiletype  = options_.get_option_value("in_type", inFiletype);
  outFiletype = options_.get_option_value("out_type", outFiletype);

#if defined(SEACAS_HAVE_MPI)
  // Should be only for parallel-aware-exodus, but not sure yet how to avoid the coupling to get
  // that define here
  compose_output = options_.get_option_value("compose", compose_output);
#endif

  {
    const char *temp = options_.retrieve("custom_field");
    if (temp != nullptr) {
      customField = temp;
    }
  }

  groupName = options_.get_option_value("extract_group", groupName);

  {
    const char *temp = options_.retrieve("field_suffix_separator");
    if (temp != nullptr) {
      fieldSuffixSeparator = temp[0];
    }
  }

  {
    const char *temp = options_.retrieve("omit_blocks");
    if (temp != nullptr) {
      auto omit_str = Ioss::tokenize(std::string(temp), ",");
      for (const auto &str : omit_str) {
        omitted_blocks.push_back(str);
      }
    }
  }

  {
    const char *temp = options_.retrieve("omit_sets");
    if (temp != nullptr) {
      auto omit_str = Ioss::tokenize(std::string(temp), ",");
      for (const auto &str : omit_str) {
        omitted_sets.push_back(str);
      }
    }
  }

  {
    const char *temp = options_.retrieve("surface_split_scheme");
    if (temp != nullptr) {
      if (std::strcmp(temp, "TOPOLOGY") == 0) {
        surface_split_type = 1;
      }
      else if (std::strcmp(temp, "ELEMENT_BLOCK") == 0) {
        surface_split_type = 2;
      }
      else if (std::strcmp(temp, "BLOCK") == 0) {
        surface_split_type = 2;
      }
      else if (std::strcmp(temp, "NO_SPLIT") == 0) {
        surface_split_type = 3;
      }
    }
  }

  {
    const char *temp = options_.retrieve("data_storage");
    if (temp != nullptr) {
      data_storage_type = 0;
      if (std::strcmp(temp, "POINTER") == 0) {
        data_storage_type = 1;
      }
#ifdef SEACAS_HAVE_KOKKOS
      else if (std::strcmp(temp, "KOKKOS_VIEW_1D") == 0) {
        data_storage_type = 3;
      }
      else if (std::strcmp(temp, "KOKKOS_VIEW_2D") == 0) {
        data_storage_type = 4;
      }
      else if (std::strcmp(temp, "KOKKOS_VIEW_2D_LAYOUTRIGHT_HOSTSPACE") == 0) {
        data_storage_type = 5;
      }
#endif

      if (data_storage_type == 0) {
        if (my_processor == 0) {
          fmt::print(stderr, "ERROR: Option data_storage must be one of\n");
#ifdef SEACAS_HAVE_KOKKOS
          fmt::print(stderr, "       POINTER, KOKKOS_VIEW_1D, KOKKOS_VIEW_2D, or "
                             "KOKKOS_VIEW_2D_LAYOUTRIGHT_HOSTSPACE\n");
#else
          fmt::print(stderr, "       POINTER\n");
#endif
        }
        return false;
      }
    }
  }

  maximum_time = options_.get_option_value("Maximum_Time", maximum_time);
  minimum_time = options_.get_option_value("Minimum_Time", minimum_time);
  time_scale   = options_.get_option_value("time_scale", time_scale);
  time_offset  = options_.get_option_value("time_offset", time_offset);

  {
    const char *temp = options_.retrieve("select_times");
    if (temp != nullptr) {
      auto time_str = Ioss::tokenize(std::string(temp), ",");
      for (const auto &str : time_str) {
        auto time = std::stod(str);
        selected_times.push_back(time);
      }
      Ioss::sort(selected_times.begin(), selected_times.end());
    }
  }

  append_time    = options_.get_option_value("append_after_time", append_time);
  flush_interval = options_.get_option_value("flush_interval", flush_interval);
  timestep_delay = options_.get_option_value("delay", timestep_delay);
  append_step    = options_.get_option_value("append_after_step", append_step);

  if (options_.retrieve("copyright") != nullptr) {
    if (my_processor == 0) {
      Ioss::Utils::copyright(std::cerr, "1999-2022");
    }
    exit(EXIT_SUCCESS);
  }

  // Parse remaining options as directory paths.
  if (option_index < argc - 1) {
    while (option_index < argc - 1) {
      inputFile.emplace_back(argv[option_index++]);
    }
    outputFile = argv[option_index];
  }
  else {
    if (my_processor == 0) {
      fmt::print(stderr, "\nERROR: input and output filename not specified\n\n");
    }
    return false;
  }

  // If inFileType and/or outFileType not specified, see if can infer from file suffix type...
  if (inFiletype == "unknown") {
    inFiletype = Ioss::Utils::get_type_from_file(inputFile[0]);
  }
  if (outFiletype == "unknown") {
    outFiletype = Ioss::Utils::get_type_from_file(outputFile);
  }
  return true;
}
