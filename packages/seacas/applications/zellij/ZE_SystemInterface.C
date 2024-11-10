// Copyright(C) 2021, 2022, 2023, 2024 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "ZE_SystemInterface.h"
#include "ZE_Version.h" // for qainfo

#include <Ioss_GetLongOpt.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_Utils.h>
#include <copyright.h>
#include <cstdlib> // for exit, strtod, strtoul, abs, etc
#include <fmt/color.h>
#include <fmt/format.h>
#include <tokenize.h>

//! \file

SystemInterface::SystemInterface(int my_rank) : myRank_(my_rank) { enroll_options(); }

namespace {
  void parse_offset(const char *tokens, vector3d &offset, int myRank);
}

void SystemInterface::enroll_options()
{
  options_.usage("[options] -lattice <lattice_definition_file>");

  options_.enroll("lattice", Ioss::GetLongOption::MandatoryValue,
                  "Name of file to read lattice definition from. [required]", "");

  options_.enroll("output", Ioss::GetLongOption::MandatoryValue,
                  "Name of output file to create. Default is `zellij-out.e`", "zellij-out.e",
                  nullptr, true);

  options_.enroll("rcb", Ioss::GetLongOption::NoValue,
                  "Use recursive coordinate bisection method to decompose the input lattice for "
                  "parallel output.",
                  nullptr);
  options_.enroll(
      "rib", Ioss::GetLongOption::NoValue,
      "Use recursive inertial bisection method to decompose the input lattice for parallel output.",
      nullptr);

  options_.enroll("hsfc", Ioss::GetLongOption::NoValue,
                  "Use hilbert space-filling curve method to decompose the input lattice for "
                  "parallel output. [default]",
                  nullptr);

  options_.enroll("linear", Ioss::GetLongOption::NoValue,
                  "Use the linear method to decompose the input lattice for parallel output.\n"
                  "\t\tElements in order first n/p to proc 0, next to proc 1.",
                  nullptr);

  options_.enroll("cyclic", Ioss::GetLongOption::NoValue,
                  "Use the cyclic method to decompose the input lattice for parallel output.\n"
                  "\t\tElements handed out to id % proc_count",
                  nullptr);

  options_.enroll("random", Ioss::GetLongOption::NoValue,
                  "Use the random method to decompose the input lattice for parallel output.\n"
                  "\t\tElements assigned randomly to processors in a way that preserves balance\n"
                  "\t\t(do *not* use for a real run)",
                  nullptr, nullptr, true);

  options_.enroll("ranks", Ioss::GetLongOption::MandatoryValue,
                  "Number of ranks to decompose mesh/lattice across", "1");

  options_.enroll("start_rank", Ioss::GetLongOption::MandatoryValue,
                  "In partial output mode, start outputting decomposed files at this rank", "0");

  options_.enroll("rank_count", Ioss::GetLongOption::MandatoryValue,
                  "In partial output or subcycle modes, output this number of ranks", "0");

  options_.enroll("subcycle", Ioss::GetLongOption::NoValue,
                  "Process cells in groups of '-rank_count'.  Helps minimize open files,\n"
                  "\t\tbut is faster than only having a single file open.",
                  nullptr);

  options_.enroll("offset", Ioss::GetLongOption::MandatoryValue,
                  "Comma-separated x,y,z offset for coordinates of the output mesh.\n"
                  "\t\tThe output coordinates will be xyz_out = xyz * scale + offset_xyz",
                  nullptr);

  options_.enroll("scale", Ioss::GetLongOption::MandatoryValue,
                  "Scale the output mesh coordinates by the specified value", "1");

  options_.enroll(
      "minimize_open_files", Ioss::GetLongOption::OptionalValue,
      "Close files after accessing them to avoid issues with too many open files.\n"
      "\t\tIf argument is 'output' then close output, if 'unit' then close unit cells;\n"
      "\t\tif 'all' or no argument close all.\n"
      "\t\tShould not need to use this option unless you get an error message "
      "indicating this issue.",
      nullptr, "all", true);

  options_.enroll("ignore_sidesets", Ioss::GetLongOption::NoValue,
                  "Do not copy any sidesets in the unit cells to the output file.", nullptr);

  options_.enroll("generate_sidesets", Ioss::GetLongOption::MandatoryValue,
                  "Which surfaces on the output mesh should have sidesets generated,\n"
                  "\t\t Valid options are:\n"
                  "\t\t 'x' or 'i' for surface on minimum X coordinate, default name = `min_i`\n"
                  "\t\t 'y' or 'j' for surface on minimum Y coordinate, default name = `min_j`\n"
                  "\t\t 'z' or 'k' for surface on minimum Z coordinate, default name = `min_k`\n"
                  "\t\t 'X' or 'I' for surface on maximum X coordinate, default name = `max_i`\n"
                  "\t\t 'Y' or 'J' for surface on maximum Y coordinate, default name = `max_j`\n"
                  "\t\t 'Z' or 'K' for surface on maximum Z coordinate, default name = `max_k`\n"
                  "\t\t For example `xyXY` would generate sidesets on min/max X and Y surfaces.",
                  "");

  options_.enroll(
      "sideset_names", Ioss::GetLongOption::MandatoryValue,
      "Specify names for one or more of the generated sidesets.\n"
      "\t\t Form is `axis:name,axis:name,...`\n"
      "\t\t where 'axis' is one of 'ijkIJKxyzXYZ', and 'name' is the name of the sideset.\n"
      "\t\t The default names are 'min_i', 'max_i', 'min_j', 'max_j', 'min_k', 'max_k'.\n"
      "\t\t For example `x:left,X:right` would name the sideset on the min x face 'left' and the "
      "max X face 'right'.",
      "", nullptr, true);

  options_.enroll("netcdf3", Ioss::GetLongOption::NoValue,
                  "Output database will be a netcdf3 "
                  "native classical netcdf file format (32-bit only)",
                  nullptr);

  options_.enroll("netcdf4", Ioss::GetLongOption::NoValue,
                  "Output database will be a netcdf4 "
                  "hdf5-based file instead of the "
                  "classical netcdf file format (default)",
                  nullptr);

  options_.enroll("netcdf5", Ioss::GetLongOption::NoValue,
                  "Output database will be a netcdf5 (CDF5) "
                  "file instead of the classical netcdf file format",
                  nullptr, nullptr, true);

  options_.enroll("32-bit", Ioss::GetLongOption::NoValue,
                  "True if forcing the use of 32-bit integers for the output file", nullptr);

  options_.enroll("64-bit", Ioss::GetLongOption::NoValue,
                  "True if forcing the use of 64-bit integers for the output file (default)",
                  nullptr, nullptr, true);

  options_.enroll("zlib", Ioss::GetLongOption::NoValue,
                  "Use the Zlib / libz compression method if compression is enabled "
                  "(default) [exodus only].",
                  nullptr);

  options_.enroll("szip", Ioss::GetLongOption::NoValue,
                  "Use SZip compression. [exodus only, enables netcdf-4]", nullptr);

  options_.enroll("compress", Ioss::GetLongOption::MandatoryValue,
                  "Specify the hdf5 zlib compression level [0..9] or szip [even, 4..32] to be used "
                  "on the output file.",
                  nullptr, nullptr, true);

  options_.enroll("separate_cells", Ioss::GetLongOption::NoValue,
                  "Do not equivalence the nodes between adjacent unit cells.", nullptr);

  options_.enroll(
      "repeat", Ioss::GetLongOption::MandatoryValue,
      "Each lattice entry will be used the specified number of times as will\n"
      "\t\teach row in the lattice (for debugging). `-repeat 2` would double the lattice.",
      "1");

  options_.enroll(
      "skip", Ioss::GetLongOption::MandatoryValue,
      "Skip the specified number of lattice entries and rows. For example, -skip 1\n"
      "\t\twould read every other entry on the row and every other row. (for debugging)",
      "1");

  options_.enroll("help", Ioss::GetLongOption::NoValue, "Print this summary and exit", nullptr);

  options_.enroll("version", Ioss::GetLongOption::NoValue, "Print version and exit", nullptr);

  options_.enroll("debug", Ioss::GetLongOption::MandatoryValue,
                  "debug level (values are or'd)\n"
                  "\t\t   1 = Exodus Verbose mode.\n"
                  "\t\t   2 = Memory and time stamp information.\n"
                  "\t\t   4 = Verbose Unit Cell information.\n"
                  "\t\t   8 = Verbose output of Grid finalization calculations.\n"
                  "\t\t  16 = Put exodus library into verbose mode.\n"
                  "\t\t  32 = Verbose decomposition information.\n"
                  "\t\t  64 = Verbose output database summary information.\n"
                  "\t\t 128 = Verbose sideset generation information.",
                  "0");

  options_.enroll("copyright", Ioss::GetLongOption::NoValue, "Show copyright and license data.",
                  nullptr);
}

bool SystemInterface::parse_options(int argc, char **argv)
{
  int option_index = options_.parse(argc, argv);
  if (option_index < 1) {
    return false;
  }

  if (options_.retrieve("help") != nullptr) {
    if (myRank_ == 0) {
      options_.usage();
      fmt::print("\n\tCan also set options via ZELLIJ_OPTIONS environment variable.\n"
                 "\n\tDocumentation: "
                 "https://sandialabs.github.io/seacas-docs/sphinx/html/index.html#zellij\n"
                 "\n\t->->-> Send email to gdsjaar@sandia.gov for zellij support.<-<-<-\n");
    }
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    // Version is printed up front, just exit...
    exit(0);
  }

  if (options_.retrieve("rcb") != nullptr) {
    decompMethod_ = "RCB";
  }

  if (options_.retrieve("rib") != nullptr) {
    decompMethod_ = "RIB";
  }

  if (options_.retrieve("hsfc") != nullptr) {
    decompMethod_ = "HSFC";
  }

  if (options_.retrieve("linear") != nullptr) {
    decompMethod_ = "LINEAR";
  }

  if (options_.retrieve("cyclic") != nullptr) {
    decompMethod_ = "CYCLIC";
  }

  if (options_.retrieve("random") != nullptr) {
    decompMethod_ = "RANDOM";
  }

  subcycle_ = (options_.retrieve("subcycle") != nullptr);

  equivalenceNodes_ = options_.retrieve("separate_cells") == nullptr;

  {
    const char *temp = options_.retrieve("minimize_open_files");
    if (temp != nullptr) {
      auto mode = Ioss::Utils::lowercase(temp);
      if (mode == "all") {
        minimizeOpenFiles_ = Minimize::ALL;
      }
      else if (mode == "unit") {
        minimizeOpenFiles_ = Minimize::UNIT;
      }
      else if (mode == "output") {
        minimizeOpenFiles_ = Minimize::OUTPUT;
      }
      else if (mode == "none") {
        minimizeOpenFiles_ = Minimize::NONE;
      }
    }
  }

  {
    const char *temp = options_.retrieve("offset");
    if (temp != nullptr) {
      parse_offset(temp, offset_, myRank_);
    }
  }

  scaleFactor_ = options_.get_option_value("scale", scaleFactor_);
  ranks_       = options_.get_option_value("ranks", ranks_);
  startRank_   = options_.get_option_value("start_rank", startRank_);
  rankCount_   = options_.get_option_value("rank_count", rankCount_);
  debugLevel_  = options_.get_option_value("debug", debugLevel_);

  if (options_.retrieve("copyright") != nullptr) {
    if (myRank_ == 0) {
      fmt::print("{}", copyright("2021"));
    }
    exit(EXIT_SUCCESS);
  }

  // Get options from environment variable also...
  char *options = getenv("ZELLIJ_OPTIONS");
  if (options != nullptr) {
    if (myRank_ == 0) {
      fmt::print(
          "\nThe following options were specified via the ZELLIJ_OPTIONS environment variable:\n"
          "\t{}\n\n",
          options);
    }
    options_.parse(options, Ioss::GetLongOption::basename(*argv));
  }

  outputName_ = options_.get_option_value("output", outputName_);
  lattice_    = options_.get_option_value("lattice", lattice_);

  ignoreInternalSidesets_ = options_.retrieve("ignore_sidesets") != nullptr;

  sidesetSurfaces_ = options_.get_option_value("generate_sidesets", sidesetSurfaces_);
  sidesetNames_    = options_.get_option_value("sideset_names", sidesetNames_);

  // Default to 64...
  ints32bit_ = options_.retrieve("32-bit") != nullptr;
  if (options_.retrieve("64-bit") != nullptr) {
    ints32bit_ = false;
  }

  if (options_.retrieve("netcdf3") != nullptr) {
    ints32bit_  = true;
    useNetcdf4_ = false;
    useNetcdf5_ = false;
  }

  if (options_.retrieve("netcdf4") != nullptr) {
    useNetcdf4_ = true;
    useNetcdf5_ = false;
  }

  if (options_.retrieve("netcdf5") != nullptr) {
    useNetcdf4_ = false;
    useNetcdf5_ = true;
  }

  if (options_.retrieve("szip") != nullptr) {
    szip_ = true;
    zlib_ = false;
  }
  zlib_ = (options_.retrieve("zlib") != nullptr);

  if (szip_ && zlib_) {
    if (myRank_ == 0) {
      fmt::print(stderr, fmt::fg(fmt::color::red),
                 "\nERROR: Only one of 'szip' or 'zlib' can be specified.\n");
    }
  }

  compressionLevel_ = options_.get_option_value("compress", compressionLevel_);

  skip_   = options_.get_option_value("skip", skip_);
  repeat_ = options_.get_option_value("repeat", repeat_);

  // Adjust start_rank and rank_count if running in parallel...
  Ioss::ParallelUtils pu{};
  if (pu.parallel_size() > 1) {
    if (subcycle_) {
      if (myRank_ == 0) {
        fmt::print(stderr, fmt::fg(fmt::color::yellow),
                   "\nWARNING: The `subcycle` option is ignored if running in parallel.\n");
      }
      subcycle_ = false;
    }
    auto size          = pu.parallel_size();
    auto ranks_per_mpi = ranks_ / size;
    auto extra         = ranks_ % size;

    auto my_rank = pu.parallel_rank();
    if (my_rank < extra) {
      startRank_ = (ranks_per_mpi + 1) * my_rank;
    }
    else {
      startRank_ = (ranks_per_mpi + 1) * extra + ranks_per_mpi * (my_rank - extra);
    }
    rankCount_ = ranks_per_mpi + (my_rank < extra ? 1 : 0);
  }
  if ((rankCount_ == 0) || (startRank_ + rankCount_ > ranks_)) {
    rankCount_ = ranks_ - startRank_;
  }

  if (lattice().empty()) {
    if (myRank_ == 0) {
      fmt::print(stderr, fmt::fg(fmt::color::red),
                 "\nERROR: Missing specification of lattice file.\n");
      options_.usage();
    }
    return false;
  }
  return true;
}

void SystemInterface::show_version()
{
  fmt::print(fmt::fg(fmt::color::cyan),
             "Zellij\n"
             "\t(A code for tiling 1 or more template databases into a single output database.)\n"
             "\t(Version: {}) Modified: {}\n",
             qainfo[2], qainfo[1]);
}

namespace {
  void parse_offset(const char *tokens, vector3d &offset, int myRank)
  {
    // Break into tokens separated by ","
    if (tokens != nullptr) {
      std::string token_string(tokens);
      auto        var_list = Ioss::tokenize(token_string, ",");

      // At this point, var_list should contain 1,2,or 3 strings
      // corresponding to the x, y, and z coordinate offsets.
      offset[0] = offset[1] = offset[2] = 0.0;

      std::string offx = var_list[0];
      double      x    = std::stod(offx);
      offset[0]        = x;

      if (var_list.size() >= 2) {
        std::string offy = var_list[1];
        double      y    = std::stod(offy);
        offset[1]        = y;
      }

      if (var_list.size() == 3) {
        std::string offz = var_list[2];
        double      z    = std::stod(offz);
        offset[2]        = z;
      }

      if (var_list.size() > 3) {
        if (myRank == 0) {
          fmt::print(stderr, fmt::fg(fmt::color::red),
                     "ERROR: Incorrect number of offset components ({}) specified -- max is 3; "
                     "ignoring extras.\n\n",
                     var_list.size());
        }
      }
    }
  }
} // namespace
