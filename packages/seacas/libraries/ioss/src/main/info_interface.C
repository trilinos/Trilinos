/*
 * Copyright(C) 1999-2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#include <cstdlib> // for exit, EXIT_SUCCESS, getenv
#include <cstring>
#include <fmt/core.h>
#include <iostream> // for operator<<, basic_ostream, etc
#include <stdio.h>
#include <string> // for char_traits, string

#include "Ioss_GetLongOpt.h" // for GetLongOption, etc
#include "Ioss_Utils.h"
#include "info_interface.h"

Info::Interface::Interface() { enroll_options(); }

void Info::Interface::enroll_options()
{
  options_.usage("[options] basename");

  options_.enroll("help", Ioss::GetLongOption::NoValue, "Print this summary and exit", nullptr);

  options_.enroll("configuration", Ioss::GetLongOption::NoValue,
                  "Show configuration of IOSS library (TPL versions)", nullptr);

  options_.enroll("db_type", Ioss::GetLongOption::MandatoryValue,
                  "Database Type: generated"
#if defined(SEACAS_HAVE_PAMGEN)
                  ", pamgen"
#endif
#if defined(SEACAS_HAVE_EXODUS)
                  ", exodus"
#endif
#if defined(SEACAS_HAVE_CGNS)
                  ", cgns"
#endif
#if defined(SEACAS_HAVE_FAODEL)
                  ", faodel"
#endif
                  ".",
                  "unknown");
  options_.enroll("in_type", Ioss::GetLongOption::MandatoryValue, "(alias for db_type)", nullptr,
                  nullptr, true);

  options_.enroll("summary", Ioss::GetLongOption::NoValue,
                  "Only output counts of nodes, elements, and entities", nullptr);

  options_.enroll("check_node_status", Ioss::GetLongOption::NoValue,
                  "Check whether there are any nodes not connected to any elements", nullptr);
  options_.enroll("adjacencies", Ioss::GetLongOption::NoValue,
                  "Calculate which element blocks touch which surfaces and other element blocks",
                  nullptr);
  options_.enroll("compute_volume", Ioss::GetLongOption::NoValue,
                  "Compute the volume of all hex elements in the mesh. Outputs min/max and count",
                  nullptr);
  options_.enroll("compute_bbox", Ioss::GetLongOption::NoValue,
                  "Compute the bounding box of all element blocks in the mesh.", nullptr, nullptr,
                  true);

  options_.enroll("disable_field_recognition", Ioss::GetLongOption::NoValue,
                  "Do not combine fields into vector, tensor fields based on basename and suffix.\n"
                  "\t\tKeep all fields on database as scalars",
                  nullptr);

  options_.enroll("field_suffix_separator", Ioss::GetLongOption::MandatoryValue,
                  "Character used to separate a field suffix from the field basename\n"
                  "\t\t when recognizing vector, tensor fields. Enter '0' for no separaor",
                  "_");

  options_.enroll("custom_field", Ioss::GetLongOption::MandatoryValue,
                  "A comma-separated list of field suffices defining a custom field that should be "
                  "recognized.\n"
                  "\t\tPrimarily used for testing",
                  nullptr);

  options_.enroll("detailed_field_info", Ioss::GetLongOption::NoValue,
                  "Output very detailed information about each field", nullptr);

  options_.enroll("use_generic_names", Ioss::GetLongOption::NoValue,
                  "Use generic names (type_id) instead of names in database", nullptr);

  options_.enroll("surface_split_scheme", Ioss::GetLongOption::MandatoryValue,
                  "Method used to split sidesets into homogeneous blocks\n"
                  "\t\tOptions are: TOPOLOGY, BLOCK, NO_SPLIT",
                  nullptr, nullptr, true);

#if defined(SEACAS_HAVE_MPI)
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
                  "\t\tElements assigned randomly to processors in a way that preserves balance.\n"
                  "\t\t(do not use for a real run)",
                  nullptr);
  options_.enroll("serialize_io_size", Ioss::GetLongOption::MandatoryValue,
                  "Number of processors that can perform simultaneous IO operations in "
                  "a parallel run; 0 to disable",
                  nullptr, nullptr, true);

#endif
  options_.enroll("list_groups", Ioss::GetLongOption::NoValue,
                  "Print a list of the names of all groups in this file and then exit.", nullptr);

  options_.enroll("group_name", Ioss::GetLongOption::MandatoryValue,
                  "List information only for the specified group.", nullptr, nullptr, true);

  options_.enroll("query_timesteps_only", Ioss::GetLongOption::NoValue,
                  "Only read and output the timestep data on the file", nullptr);
  options_.enroll("64-bit", Ioss::GetLongOption::NoValue, "Use 64-bit integers", nullptr);
  options_.enroll("version", Ioss::GetLongOption::NoValue, "Print version and exit", nullptr);

  options_.enroll("copyright", Ioss::GetLongOption::NoValue, "Show copyright and license data.",
                  nullptr);
}

bool Info::Interface::parse_options(int argc, char **argv)
{
  // Get options from environment variable also...
  char *options = getenv("IO_INFO_OPTIONS");
  if (options != nullptr) {
    fmt::print(
        stderr,
        "\nThe following options were specified via the IO_INFO_OPTIONS environment variable:\n"
        "\t{}\n\n",
        options);
    options_.parse(options, Ioss::GetLongOption::basename(*argv));
  }

  int option_index = options_.parse(argc, argv);
  if (option_index < 1) {
    return false;
  }

  if (options_.retrieve("help") != nullptr) {
    options_.usage(std::cerr);
    fmt::print(stderr,
               "\n\tCan also set options via IO_INFO_OPTIONS environment variable.\n\n"
               "\tDocumentation: "
               "https://sandialabs.github.io/seacas-docs/sphinx/html/index.html#io-info\n\n"
               "\t->->-> Send email to gdsjaar@sandia.gov for {} support.<-<-<-\n",
               options_.program_name());
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    // Version is printed up front, just exit...
    exit(0);
  }

  checkNodeStatus_ = options_.retrieve("check_node_status") != nullptr;
  adjacencies_     = options_.retrieve("adjacencies") != nullptr;
  ints64Bit_       = options_.retrieve("64-bit") != nullptr;
  computeVolume_   = options_.retrieve("compute_volume") != nullptr;
  computeBBox_     = options_.retrieve("compute_bbox") != nullptr;
  listGroups_      = options_.retrieve("list_groups") != nullptr;
  useGenericNames_ = options_.retrieve("use_generic_names") != nullptr;
  summary_         = options_.retrieve("summary") != nullptr;
  showConfig_      = options_.retrieve("configuration") != nullptr;
  queryTimeOnly_   = options_.retrieve("query_timesteps_only") != nullptr;
  fieldDetails_    = options_.retrieve("detailed_field_info") != nullptr;

  filetype_  = options_.get_option_value("db_type", filetype_);
  filetype_  = options_.get_option_value("in_type", filetype_);
  groupname_ = options_.get_option_value("group_name", groupname_);

  {
    const char *temp = options_.retrieve("surface_split_scheme");
    if (temp != nullptr) {
      if (std::strcmp(temp, "TOPOLOGY") == 0) {
        surfaceSplitScheme_ = 1;
      }
      else if (std::strcmp(temp, "ELEMENT_BLOCK") == 0) {
        surfaceSplitScheme_ = 2;
      }
      else if (std::strcmp(temp, "BLOCK") == 0) {
        surfaceSplitScheme_ = 2;
      }
      else if (std::strcmp(temp, "NO_SPLIT") == 0 || std::strcmp(temp, "NOSPLIT") == 0) {
        // Backward compatibility using "NOSPLIT" after bug fix...
        surfaceSplitScheme_ = 3;
      }
    }
  }

  disableFieldRecognition_ = options_.retrieve("disable_field_recognition") != nullptr;

  {
    const char *temp = options_.retrieve("custom_field");
    if (temp != nullptr) {
      customField_ = temp;
    }
  }

  {
    const char *temp = options_.retrieve("field_suffix_separator");
    if (temp != nullptr) {
      fieldSuffixSeparator_ = temp[0];
    }
  }

#if defined(SEACAS_HAVE_MPI)
#if !defined(NO_ZOLTAN_SUPPORT)
  if (options_.retrieve("rcb") != nullptr) {
    decompMethod_ = "RCB";
  }

  if (options_.retrieve("rib") != nullptr) {
    decompMethod_ = "RIB";
  }

  if (options_.retrieve("hsfc") != nullptr) {
    decompMethod_ = "HSFC";
  }
#endif

#if !defined(NO_PARMETIS_SUPPORT)
  if (options_.retrieve("metis_sfc") != nullptr) {
    decompMethod_ = "METIS_SFC";
  }

  if (options_.retrieve("kway") != nullptr) {
    decompMethod_ = "KWAY";
  }

  if (options_.retrieve("kway_geom") != nullptr) {
    decompMethod_ = "KWAY_GEOM";
  }
#endif

  if (options_.retrieve("linear") != nullptr) {
    decompMethod_ = "LINEAR";
  }

  if (options_.retrieve("cyclic") != nullptr) {
    decompMethod_ = "CYCLIC";
  }

  if (options_.retrieve("random") != nullptr) {
    decompMethod_ = "RANDOM";
  }
#endif

  if (options_.retrieve("copyright") != nullptr) {
    Ioss::Utils::copyright(std::cerr, "1999-2022");
    exit(EXIT_SUCCESS);
  }

  // Parse remaining options as directory paths.
  if (!show_config()) {
    if (option_index < argc) {
      filename_ = argv[option_index];
    }
    else {
      fmt::print(stderr, "\nERROR: filename not specified\n\n");
      return false;
    }

    if (filetype_ == "unknown") {
      filetype_ = Ioss::Utils::get_type_from_file(filename_);
    }
  }

  return true;
}
