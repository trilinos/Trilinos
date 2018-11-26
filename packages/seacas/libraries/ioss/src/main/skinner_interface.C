/*
 * Copyright(C) 1999-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "Ioss_GetLongOpt.h" // for GetLongOption, etc
#include "Ioss_Utils.h"
#include "skinner_interface.h"
#include <cstddef>  // for nullptr
#include <cstdlib>  // for exit, EXIT_SUCCESS, getenv
#include <iostream> // for operator<<, basic_ostream, etc
#include <string>   // for char_traits, string

Skinner::Interface::Interface() { enroll_options(); }

Skinner::Interface::~Interface() = default;

void Skinner::Interface::enroll_options()
{
  options_.usage("[options] input_file[s] output_file");

  options_.enroll("help", Ioss::GetLongOption::NoValue, "Print this summary and exit", nullptr);

  options_.enroll("version", Ioss::GetLongOption::NoValue, "Print version and exit", nullptr);

  options_.enroll("64-bit", Ioss::GetLongOption::NoValue, "True if using 64-bit integers", nullptr);
  options_.enroll("in_type", Ioss::GetLongOption::MandatoryValue,
                  "Database type for input file: pamgen|generated|exodus. exodus is the default.",
                  "exodus");

  options_.enroll("out_type", Ioss::GetLongOption::MandatoryValue,
                  "Database type for output file: exodus. exodus is the default.", "exodus");

  options_.enroll("no_output", Ioss::GetLongOption::NoValue,
                  "Do not produce output file, just generate the faces", nullptr);

  options_.enroll("ignore_face_ids", Ioss::GetLongOption::NoValue,
                  "Ignore internal face ids and just use 1..num_face", nullptr);

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
                  nullptr);

  options_.enroll(
      "compose", Ioss::GetLongOption::MandatoryValue,
      "Specify the parallel-io method to be used to output a single file in a parallel run. "
      "Options are default, mpiio, mpiposix, pnetcdf",
      nullptr);

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
                  "Use the linear method to decompose the input mesh in a parallel run. "
                  "elements in order first n/p to proc 0, next to proc 1.",
                  nullptr);

  options_.enroll("cyclic", Ioss::GetLongOption::NoValue,
                  "Use the cyclic method to decompose the input mesh in a parallel run. "
                  "elements handed out to id % proc_count",
                  nullptr);

  options_.enroll("random", Ioss::GetLongOption::NoValue,
                  "Use the random method to decompose the input mesh in a parallel run."
                  "elements assigned randomly to processors in a way that preserves balance (do "
                  "not use for a real run)",
                  nullptr);

  options_.enroll("external", Ioss::GetLongOption::NoValue,
                  "Files are decomposed externally into a file-per-processor in a parallel run.",
                  nullptr);

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
    std::cerr << "\nThe following options were specified via the IO_SKINNER_OPTIONS environment "
                 "variable:\n"
              << "\t" << options << "\n\n";
    options_.parse(options, options_.basename(*argv));
  }

  int option_index = options_.parse(argc, argv);
  if (option_index < 1) {
    return false;
  }

  if (options_.retrieve("help") != nullptr) {
    options_.usage();
    std::cerr << "\n\tCan also set options via IO_SKINNER_OPTIONS environment variable.\n\n";
    std::cerr << "\n\t->->-> Send email to gdsjaar@sandia.gov for epu support.<-<-<-\n";
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version") != nullptr) {
    // Version is printed up front, just exit...
    exit(0);
  }

  if (options_.retrieve("64-bit") != nullptr) {
    ints64Bit_ = true;
  }

  if (options_.retrieve("netcdf4") != nullptr) {
    netcdf4_ = true;
  }

  if (options_.retrieve("shuffle") != nullptr) {
    shuffle = true;
  }

  if (options_.retrieve("no_output") != nullptr) {
    noOutput_ = true;
  }

  if (options_.retrieve("ignore_face_ids") != nullptr) {
    ignoreFaceIds_ = true;
  }

  {
    const char *temp = options_.retrieve("compress");
    if (temp != nullptr) {
      compression_level = std::strtol(temp, nullptr, 10);
    }
  }

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

  if (options_.retrieve("debug") != nullptr) {
    debug = true;
  }

  if (options_.retrieve("statistics") != nullptr) {
    statistics = true;
  }

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

  {
    const char *temp = options_.retrieve("compose");
    if (temp != nullptr) {
      compose_output = Ioss::Utils::lowercase(temp);
    }
  }

  if (options_.retrieve("copyright") != nullptr) {
    std::cerr << "\n"
              << "Copyright(C) 1999-2017 National Technology & Engineering Solutions\n"
              << "of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with\n"
              << "NTESS, the U.S. Government retains certain rights in this software.\n\n"
              << "Redistribution and use in source and binary forms, with or without\n"
              << "modification, are permitted provided that the following conditions are\n"
              << "met:\n\n "
              << "    * Redistributions of source code must retain the above copyright\n"
              << "      notice, this list of conditions and the following disclaimer.\n\n"
              << "    * Redistributions in binary form must reproduce the above\n"
              << "      copyright notice, this list of conditions and the following\n"
              << "      disclaimer in the documentation and/or other materials provided\n"
              << "      with the distribution.\n\n"
              << "    * Neither the name of NTESS nor the names of its\n"
              << "      contributors may be used to endorse or promote products derived\n"
              << "      from this software without specific prior written permission.\n\n"
              << "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
              << "\" AS IS \" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
              << "LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR\n"
              << "A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT\n"
              << "OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,\n"
              << "SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT\n"
              << "LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,\n"
              << "DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY\n"
              << "THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n"
              << "(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE\n"
              << "OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n\n";
    exit(EXIT_SUCCESS);
  }

  // Parse remaining options as directory paths.
  if (option_index < argc - 1) {
    inputFile_  = argv[option_index++];
    outputFile_ = argv[option_index];
  }
  else {
    std::cerr << "\nERROR: input and output filename not specified\n\n";
    return false;
  }
  return true;
}
