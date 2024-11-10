/*
 * Copyright(C) 1999-2024 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#include <cstdlib> // for exit, EXIT_SUCCESS, getenv
#include <fmt/core.h>
#include <iostream> // for operator<<, basic_ostream, etc
#include <stdio.h>
#include <string> // for char_traits, string

#include "Ioss_GetLongOpt.h" // for GetLongOption, etc
#include "Ioss_Sort.h"
#include "Ioss_Utils.h"
#include "skinner_interface.h"
#include "tokenize.h"

Skinner::Interface::Interface() { enroll_options(); }

void Skinner::Interface::enroll_options()
{
  options_.usage("[options] input_file[s] output_file");

  options_.enroll("help", Ioss::GetLongOption::NoValue, "Print this summary and exit", nullptr);

  options_.enroll("version", Ioss::GetLongOption::NoValue, "Print version and exit", nullptr);

  options_.enroll("64-bit", Ioss::GetLongOption::NoValue, "True if using 64-bit integers", nullptr);
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
#if defined(SEACAS_HAVE_DATAWAREHOUSE)
                  "|data_warehouse"
#endif
                  ".\n\t\tIf not specified, guess from extension or exodus is the default.",
                  "unknown");

  options_.enroll("out_type", Ioss::GetLongOption::MandatoryValue,
                  "Database type for output file:"
#if defined(SEACAS_HAVE_EXODUS)
                  " exodus"
#endif
#if defined(SEACAS_HAVE_CGNS)
                  " cgns"
#endif
                  ".\n\t\tIf not specified, guess from extension or exodus is the default.",
                  "unknown");

  options_.enroll("no_output", Ioss::GetLongOption::NoValue,
                  "Do not produce output file, just generate the faces", nullptr, nullptr, true);

  options_.enroll(
      "output_transient", Ioss::GetLongOption::NoValue,
      "Transfer nodal and element transient data from the input mesh to the output mesh.", nullptr);

  options_.enroll("Maximum_Time", Ioss::GetLongOption::MandatoryValue,
                  "Maximum time on input database to transfer to output database", nullptr);

  options_.enroll("Minimum_Time", Ioss::GetLongOption::MandatoryValue,
                  "Minimum time on input database to transfer to output database", nullptr);

  options_.enroll("select_times", Ioss::GetLongOption::MandatoryValue,
                  "comma-separated list of times that should be transferred to output database",
                  nullptr, nullptr, true);

  options_.enroll("ignore_face_hash_ids", Ioss::GetLongOption::NoValue,
                  "Don't use face ids from hash of node ids; just use 1..num_face", nullptr);

  options_.enroll("blocks", Ioss::GetLongOption::NoValue,
                  "Skin block-by-block instead of entire model boundary", nullptr, nullptr, true);

  options_.enroll("netcdf4", Ioss::GetLongOption::NoValue,
                  "Output database will be a netcdf4 "
                  "hdf5-based file instead of the "
                  "classical netcdf file format",
                  nullptr);

  options_.enroll("shuffle", Ioss::GetLongOption::NoValue,
                  "Use a netcdf4 hdf5-based file and use hdf5s shuffle mode with compression.",
                  nullptr);

  options_.enroll("compress", Ioss::GetLongOption::MandatoryValue,
                  "Specify the hdf5 compression level [0..9] to be used on the output file.",
                  nullptr, nullptr, true);

#if defined(SEACAS_HAVE_MPI)
  options_.enroll(
      "compose", Ioss::GetLongOption::OptionalValue,
      "If no argument, specify single-file output; if 'external', then file-per-processor.\n"
      "\t\tAll other options are ignored and just exist for backward-compatibility",
      nullptr, "nullptr, true");

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

  options_.enroll("linear", Ioss::GetLongOption::NoValue,
                  "Use the linear method to decompose the input mesh in a parallel run.\n"
                  "\t\telements in order first n/p to proc 0, next to proc 1.",
                  nullptr);

  options_.enroll("cyclic", Ioss::GetLongOption::NoValue,
                  "Use the cyclic method to decompose the input mesh in a parallel run.\n"
                  "\t\telements handed out to id % proc_count",
                  nullptr);

  options_.enroll(
      "random", Ioss::GetLongOption::NoValue,
      "Use the random method to decompose the input mesh in a parallel run.\n"
      "\t\telements assigned randomly to processors in a way that preserves balance (do "
      "not use for a real run)",
      nullptr);

  options_.enroll("external", Ioss::GetLongOption::NoValue,
                  "Files are decomposed externally into a file-per-processor in a parallel run.",
                  nullptr, nullptr, true);
#endif

  options_.enroll("debug", Ioss::GetLongOption::NoValue, "turn on debugging output", nullptr);

  options_.enroll("statistics", Ioss::GetLongOption::NoValue,
                  "output parallel io timing statistics", nullptr);

  options_.enroll("copyright", Ioss::GetLongOption::NoValue, "Show copyright and license data.",
                  nullptr);
}

bool Skinner::Interface::parse_options(int argc, char **argv)
{
  // Get options from environment variable also...
  char *options = getenv("IO_SKINNER_OPTIONS");
  if (options != nullptr) {
    fmt::print(stderr,
               "\nThe following options were specified via the IO_SKINNER_OPTIONS environment "
               "variable:\n"
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
               "\n\tCan also set options via IO_SKINNER_OPTIONS environment variable.\n\n"
               "\tDocumentation: "
               "https://sandialabs.github.io/seacas-docs/sphinx/html/index.html#skinner\n\n"
               "\t->->-> Send email to gdsjaar@sandia.gov for {} support.<-<-<-\n",
               options_.program_name());
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    // Version is printed up front, just exit...
    exit(0);
  }

  ints64Bit_       = options_.retrieve("64-bit") != nullptr;
  netcdf4_         = options_.retrieve("netcdf4") != nullptr;
  shuffle          = options_.retrieve("shuffle") != nullptr;
  noOutput_        = options_.retrieve("no_output") != nullptr;
  outputTransient_ = options_.retrieve("output_transient") != nullptr;
  useFaceHashIds_  = options_.retrieve("ignore_face_hash_ids") == nullptr;
  debug            = options_.retrieve("debug") != nullptr;
  statistics       = options_.retrieve("statistics") != nullptr;
  blocks_          = options_.retrieve("blocks") != nullptr;

  maximum_time = options_.get_option_value("Maximum_Time", maximum_time);
  minimum_time = options_.get_option_value("Minimum_Time", minimum_time);

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

  {
    const char *temp = options_.retrieve("compress");
    if (temp != nullptr) {
      compression_level = std::strtol(temp, nullptr, 10);
    }
  }

#if defined(SEACAS_HAVE_MPI)
  if (options_.retrieve("rcb") != nullptr) {
    decomp_method = "RCB";
  }

  if (options_.retrieve("rib") != nullptr) {
    decomp_method = "RIB";
  }

  if (options_.retrieve("hsfc") != nullptr) {
    decomp_method = "HSFC";
  }

  if (options_.retrieve("metis_sfc") != nullptr) {
    decomp_method = "METIS_SFC";
  }

  if (options_.retrieve("kway") != nullptr) {
    decomp_method = "KWAY";
  }

  if (options_.retrieve("kway_geom") != nullptr) {
    decomp_method = "KWAY_GEOM";
  }

  if (options_.retrieve("linear") != nullptr) {
    decomp_method = "LINEAR";
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

  {
    const char *temp = options_.retrieve("compose");
    if (temp != nullptr) {
      compose_output = Ioss::Utils::lowercase(temp);
    }
  }

#endif

  {
    const char *temp = options_.retrieve("in_type");
    if (temp != nullptr) {
      inFiletype_ = temp;
    }
  }

  {
    const char *temp = options_.retrieve("out_type");
    if (temp != nullptr) {
      outFiletype_ = temp;
    }
  }

  if (options_.retrieve("copyright") != nullptr) {
    Ioss::Utils::copyright(std::cerr, "1999-2023");
    exit(EXIT_SUCCESS);
  }

  // Parse remaining options as directory paths.
  if (option_index < argc) {
    inputFile_ = argv[option_index++];
  }

  if (option_index < argc && !noOutput_) {
    outputFile_ = argv[option_index];
  }

  if (inputFile_.empty() || (!noOutput_ && outputFile_.empty())) {
    fmt::print(stderr, "\nERROR: input and output filename not specified\n\n");
    return false;
  }

  // If inFileType and/or outFileType not specified, see if can infer from file suffix type...
  if (inFiletype_ == "unknown") {
    inFiletype_ = Ioss::Utils::get_type_from_file(inputFile_);
  }
  if (!noOutput_ && outFiletype_ == "unknown") {
    outFiletype_ = Ioss::Utils::get_type_from_file(outputFile_);
  }
  return true;
}
