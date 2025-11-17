/*
 * Copyright(C) 1999-2025 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#include <cstdlib> // for exit, EXIT_SUCCESS, getenv
#include <cstring>
#include <fmt/format.h>
#include <iostream> // for operator<<, basic_ostream, etc
#include <stdio.h>
#include <string> // for char_traits, string

#include "Ioss_GetLongOpt.h" // for GetLongOption, etc
#include "Ioss_Utils.h"
#include "info_interface.h"

Info::Interface::Interface(std::string app_version) : version(std::move(app_version))
{
  enroll_options();
}

void Info::Interface::enroll_options()
{
  options_.usage("[options] basename");

  options_.enroll("help", Ioss::GetLongOption::OptType::NoValue, "Print this summary and exit",
                  nullptr);

  options_.enroll("configuration", Ioss::GetLongOption::OptType::NoValue,
                  "Show configuration of IOSS library (TPL versions)", nullptr);

  options_.enroll("db_type", Ioss::GetLongOption::OptType::MandatoryValue,
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
#if defined(SEACAS_HAVE_S3)
                  ", s3"
#endif
                  ".",
                  "unknown");
  options_.enroll("in_type", Ioss::GetLongOption::OptType::MandatoryValue, "(alias for db_type)",
                  nullptr, nullptr, true);

  options_.enroll("summary", Ioss::GetLongOption::OptType::NoValue,
                  "Only output counts of nodes, elements, and entities", nullptr);

  options_.enroll("check_node_status", Ioss::GetLongOption::OptType::NoValue,
                  "Check whether there are any nodes not connected to any elements", nullptr);
  options_.enroll("adjacencies", Ioss::GetLongOption::OptType::NoValue,
                  "Calculate which element blocks touch which surfaces and other element blocks",
                  nullptr);
  options_.enroll("compute_volume", Ioss::GetLongOption::OptType::NoValue,
                  "Compute the volume of all hex elements in the mesh. Outputs min/max and count",
                  nullptr);
  options_.enroll("compute_bbox", Ioss::GetLongOption::OptType::NoValue,
                  "Compute the bounding box of all element blocks in the mesh.", nullptr, nullptr,
                  true);

  options_.enroll("disable_field_recognition", Ioss::GetLongOption::OptType::NoValue,
                  "Do not combine fields into vector, tensor fields based on basename and suffix.\n"
                  "\t\tKeep all fields on database as scalars",
                  nullptr);

  options_.enroll("field_suffix_separator", Ioss::GetLongOption::OptType::MandatoryValue,
                  "Character used to separate a field suffix from the field basename\n"
                  "\t\t when recognizing vector, tensor fields. Enter '0' for no separaor",
                  "_");

  options_.enroll("custom_field", Ioss::GetLongOption::OptType::MandatoryValue,
                  "A comma-separated list of field suffices defining a custom field that should be "
                  "recognized.\n"
                  "\t\tPrimarily used for testing",
                  nullptr);

  options_.enroll("detailed_field_info", Ioss::GetLongOption::OptType::NoValue,
                  "Output very detailed information about each field", nullptr);

  options_.enroll("use_generic_names", Ioss::GetLongOption::OptType::NoValue,
                  "Use generic names (type_id) instead of names in database", nullptr);

  options_.enroll("surface_split_scheme", Ioss::GetLongOption::OptType::MandatoryValue,
                  "Method used to split sidesets into homogeneous blocks\n"
                  "\t\tOptions are: TOPOLOGY, BLOCK, NO_SPLIT",
                  nullptr, nullptr, true);

#if defined(SEACAS_HAVE_MPI)
#if !defined(NO_ZOLTAN_SUPPORT)
  options_.enroll(
      "rcb", Ioss::GetLongOption::OptType::NoValue,
      "Use recursive coordinate bisection method to decompose the input mesh in a parallel run.",
      nullptr);
  options_.enroll(
      "rib", Ioss::GetLongOption::OptType::NoValue,
      "Use recursive inertial bisection method to decompose the input mesh in a parallel run.",
      nullptr);

  options_.enroll(
      "hsfc", Ioss::GetLongOption::OptType::NoValue,
      "Use hilbert space-filling curve method to decompose the input mesh in a parallel run.",
      nullptr);
#endif

#if !defined(NO_PARMETIS_SUPPORT)
  options_.enroll(
      "metis_sfc", Ioss::GetLongOption::OptType::NoValue,
      "Use the metis space-filling-curve method to decompose the input mesh in a parallel run.",
      nullptr);

  options_.enroll(
      "kway", Ioss::GetLongOption::OptType::NoValue,
      "Use the metis kway graph-based method to decompose the input mesh in a parallel run.",
      nullptr);

  options_.enroll("kway_geom", Ioss::GetLongOption::OptType::NoValue,
                  "Use the metis kway graph-based method with geometry speedup to decompose the "
                  "input mesh in a parallel run.",
                  nullptr);
#endif

  options_.enroll("linear", Ioss::GetLongOption::OptType::NoValue,
                  "Use the linear method to decompose the input mesh in a parallel run.\n"
                  "\t\tElements in order first n/p to proc 0, next to proc 1.",
                  nullptr);

  options_.enroll("cyclic", Ioss::GetLongOption::OptType::NoValue,
                  "Use the cyclic method to decompose the input mesh in a parallel run.\n"
                  "\t\tElements handed out to id % proc_count",
                  nullptr);

  options_.enroll("random", Ioss::GetLongOption::OptType::NoValue,
                  "Use the random method to decompose the input mesh in a parallel run.\n"
                  "\t\tElements assigned randomly to processors in a way that preserves balance.\n"
                  "\t\t(do not use for a real run)",
                  nullptr);
  options_.enroll("serialize_io_size", Ioss::GetLongOption::OptType::MandatoryValue,
                  "Number of processors that can perform simultaneous IO operations in "
                  "a parallel run; 0 to disable",
                  nullptr, nullptr, true);

#endif

  options_.enroll("list_change_sets", Ioss::GetLongOption::OptType::NoValue,
                  "Print a list of the names of all change_sets (previuosly groups) in this file "
                  "and then exit.",
                  nullptr);
  options_.enroll("list_groups", Ioss::GetLongOption::OptType::NoValue,
                  "[deprecated] Use --list_change_sets", nullptr);

  options_.enroll("change_set_name", Ioss::GetLongOption::OptType::MandatoryValue,
                  "List information only for the specified comma-separated list of change_set(s) "
                  "or `ALL` to list for all.",
                  nullptr);
  options_.enroll("group_name", Ioss::GetLongOption::OptType::MandatoryValue,
                  "[deprecated] Use --change_set_name.", nullptr, nullptr, true);

  options_.enroll("query_timesteps_only", Ioss::GetLongOption::OptType::NoValue,
                  "Only read and output the timestep data on the file", nullptr);
  options_.enroll(
      "show_timestep_times", Ioss::GetLongOption::OptType::NoValue,
      "Show the times for all timesteps. By default only shows minimum and maximum time.", nullptr);
  options_.enroll("64-bit", Ioss::GetLongOption::OptType::NoValue, "Use 64-bit integers", nullptr);
  options_.enroll("version", Ioss::GetLongOption::OptType::NoValue, "Print version and exit",
                  nullptr);

  options_.enroll("copyright", Ioss::GetLongOption::OptType::NoValue,
                  "Show copyright and license data.", nullptr);
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
               "\t->->-> Send email to sierra-help@sandia.gov for {} support.<-<-<-\n",
               options_.program_name());
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    fmt::print(stderr, "IO_INFO\tVersion: {}\n", version);
    exit(0);
  }

  checkNodeStatus_ = options_.retrieve("check_node_status") != nullptr;
  adjacencies_     = options_.retrieve("adjacencies") != nullptr;
  ints64Bit_       = options_.retrieve("64-bit") != nullptr;
  computeVolume_   = options_.retrieve("compute_volume") != nullptr;
  computeBBox_     = options_.retrieve("compute_bbox") != nullptr;
  listChangeSets_  = options_.retrieve("list_change_sets") != nullptr ||
                    options_.retrieve("list_groups") != nullptr;
  useGenericNames_ = options_.retrieve("use_generic_names") != nullptr;
  summary_         = options_.retrieve("summary") != nullptr;
  showConfig_      = options_.retrieve("configuration") != nullptr;
  queryTimeOnly_   = options_.retrieve("query_timesteps_only") != nullptr;
  showTimes_       = options_.retrieve("show_timestep_times") != nullptr;
  fieldDetails_    = options_.retrieve("detailed_field_info") != nullptr;

  filetype_      = options_.get_option_value("db_type", filetype_);
  filetype_      = options_.get_option_value("in_type", filetype_);
  changeSetName_ = options_.get_option_value("change_set_name", changeSetName_);
  changeSetName_ = options_.get_option_value("group_name", changeSetName_);

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
