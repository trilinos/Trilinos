/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#include "Ioss_CodeTypes.h"
#include "Ioss_FileInfo.h"
#include "Ioss_GetLongOpt.h" // for GetLongOption, etc
#include "fmt/ostream.h"
#include "info_interface.h"
#include <cstddef>  // for nullptr
#include <cstdlib>  // for exit, EXIT_SUCCESS, getenv
#include <iostream> // for operator<<, basic_ostream, etc
#include <string>   // for char_traits, string

namespace {
  std::string get_type_from_file(const std::string &filename)
  {
    Ioss::FileInfo file(filename);
    auto           extension = file.extension();
    if (extension == "e" || extension == "g" || extension == "gen" || extension == "exo") {
      return "exodus";
    }
    else if (extension == "cgns") {
      return "cgns";
    }
    else {
      // "exodus" is default...
      return "exodus";
    }
  }
} // namespace

Info::Interface::Interface() { enroll_options(); }

Info::Interface::~Interface() = default;

void Info::Interface::enroll_options()
{
  options_.usage("[options] basename");

  options_.enroll("help", Ioss::GetLongOption::NoValue, "Print this summary and exit", nullptr);

  options_.enroll("version", Ioss::GetLongOption::NoValue, "Print version and exit", nullptr);

  options_.enroll("check_node_status", Ioss::GetLongOption::NoValue,
                  "Check whether there are any nodes not connected to any elements", nullptr);
  options_.enroll("adjacencies", Ioss::GetLongOption::NoValue,
                  "Calculate which element blocks touch which surfaces and other element blocks",
                  nullptr);
  options_.enroll("64-bit", Ioss::GetLongOption::NoValue, "True if using 64-bit integers", nullptr);
  options_.enroll("compute_volume", Ioss::GetLongOption::NoValue,
                  "Compute the volume of all hex elements in the mesh. Outputs min/max and count",
                  nullptr);
  options_.enroll("compute_bbox", Ioss::GetLongOption::NoValue,
                  "Compute the bounding box of all element blocks in the mesh.", nullptr);

  options_.enroll("list_groups", Ioss::GetLongOption::NoValue,
                  "Print a list of the names of all groups in this file and then exit.", nullptr);

  options_.enroll("disable_field_recognition", Ioss::GetLongOption::NoValue,
                  "Do not combine fields into vector, tensor fields based on basename and suffix.\n"
                  "\t\tKeep all fields on database as scalars",
                  nullptr);

  options_.enroll("field_suffix_separator", Ioss::GetLongOption::MandatoryValue,
                  "Character used to separate a field suffix from the field basename\n"
                  "\t\t when recognizing vector, tensor fields. Enter '0' for no separaor",
                  "_");

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
#if defined(SEACAS_HAVE_DATAWAREHOUSE)
                  ", data_warehouse"
#endif
                  ".",
                  "unknown");
  options_.enroll("in_type", Ioss::GetLongOption::MandatoryValue, "(alias for db_type)", nullptr);

  options_.enroll("group_name", Ioss::GetLongOption::MandatoryValue,
                  "List information only for the specified group.", nullptr);

  options_.enroll("use_generic_names", Ioss::GetLongOption::NoValue,
                  "True to use generic names (type_id) instead of names in database", nullptr);

  options_.enroll("summary", Ioss::GetLongOption::NoValue,
                  "Only output counts of nodes, elements, and entities", nullptr);

  options_.enroll("configuration", Ioss::GetLongOption::NoValue,
                  "Show configuration of IOSS library (TPL versions)", nullptr);

  options_.enroll("surface_split_scheme", Ioss::GetLongOption::MandatoryValue,
                  "Method used to split sidesets into homogeneous blocks\n"
                  "\t\tOptions are: TOPOLOGY, BLOCK, NOSPLIT",
                  "TOPOLOGY");

  options_.enroll("copyright", Ioss::GetLongOption::NoValue, "Show copyright and license data.",
                  nullptr);

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
                  nullptr);
#endif
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
               "\t->->-> Send email to gdsjaar@sandia.gov for {} support.<-<-<-\n",
               options_.program_name());
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    // Version is printed up front, just exit...
    exit(0);
  }

  if (options_.retrieve("check_node_status") != nullptr) {
    checkNodeStatus_ = true;
  }

  if (options_.retrieve("adjacencies") != nullptr) {
    adjacencies_ = true;
  }

  if (options_.retrieve("64-bit") != nullptr) {
    ints64Bit_ = true;
  }

  if (options_.retrieve("compute_volume") != nullptr) {
    computeVolume_ = true;
  }

  if (options_.retrieve("compute_bbox") != nullptr) {
    computeBBox_ = true;
  }

  if (options_.retrieve("list_groups") != nullptr) {
    listGroups_ = true;
  }

  if (options_.retrieve("use_generic_names") != nullptr) {
    useGenericNames_ = true;
  }

  if (options_.retrieve("summary") != nullptr) {
    summary_ = true;
  }

  if (options_.retrieve("configuration") != nullptr) {
    showConfig_ = true;
  }

  {
    const char *temp = options_.retrieve("db_type");
    if (temp != nullptr) {
      filetype_ = temp;
    }
  }

  {
    const char *temp = options_.retrieve("in_type");
    if (temp != nullptr) {
      filetype_ = temp;
    }
  }

  {
    const char *temp = options_.retrieve("group_name");
    if (temp != nullptr) {
      groupname_ = temp;
    }
  }

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
      else if (std::strcmp(temp, "NO_SPLIT") == 0) {
        surfaceSplitScheme_ = 3;
      }
    }
  }

  if (options_.retrieve("disable_field_recognition") != nullptr) {
    disableFieldRecognition_ = true;
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
    fmt::print(stderr, "\n"
                       "Copyright(C) 1999-2017 National Technology & Engineering Solutions\n"
                       "of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with\n"
                       "NTESS, the U.S. Government retains certain rights in this software.\n\n"
                       "Redistribution and use in source and binary forms, with or without\n"
                       "modification, are permitted provided that the following conditions are\n"
                       "met:\n\n "
                       "    * Redistributions of source code must retain the above copyright\n"
                       "      notice, this list of conditions and the following disclaimer.\n\n"
                       "    * Redistributions in binary form must reproduce the above\n"
                       "      copyright notice, this list of conditions and the following\n"
                       "      disclaimer in the documentation and/or other materials provided\n"
                       "      with the distribution.\n\n"
                       "    * Neither the name of NTESS nor the names of its\n"
                       "      contributors may be used to endorse or promote products derived\n"
                       "      from this software without specific prior written permission.\n\n"
                       "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
                       "\" AS IS \" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
                       "LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR\n"
                       "A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT\n"
                       "OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,\n"
                       "SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT\n"
                       "LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,\n"
                       "DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY\n"
                       "THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n"
                       "(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE\n"
                       "OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n\n");
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
      filetype_ = get_type_from_file(filename_);
    }
  }

  return true;
}
