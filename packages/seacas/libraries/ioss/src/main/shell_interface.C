/*
 * Copyright(C) 2013 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
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
 *     * Neither the name of Sandia Corporation nor the names of its
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
 * 
 */
#include "shell_interface.h"
#include <stddef.h>                     // for NULL
#include <cstdlib>                      // for exit, strtod, EXIT_SUCCESS, etc
#include <cstring>                      // for strcmp
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <string>                       // for string, char_traits
#include <vector>                       // for vector
#include "Ioss_GetLongOpt.h"            // for GetLongOption, etc
#include "Ioss_Utils.h"                 // for Utils



#define NPOS std::string::npos

IOShell::Interface::Interface()
  :   compose_output("none"), maximum_time(0.0), minimum_time(0.0), surface_split_type(1), compression_level(0),
      shuffle(false), debug(false), statistics(false), do_transform_fields(false), ints_64_bit(false),
      reals_32_bit(false), netcdf4(false), fieldSuffixSeparator('_')
{
  enroll_options();
}

IOShell::Interface::~Interface() {}

void IOShell::Interface::enroll_options()
{
  options_.usage("[options] input_file[s] output_file");

  options_.enroll("help", Ioss::GetLongOption::NoValue,
		  "Print this summary and exit", 0);

  options_.enroll("version", Ioss::GetLongOption::NoValue,
		  "Print version and exit", NULL);

  options_.enroll("in_type", Ioss::GetLongOption::MandatoryValue,
		  "Database type for input file: pamgen|generated|exodus. exodus is the default.",
		  "exodus");

  options_.enroll("out_type", Ioss::GetLongOption::MandatoryValue,
		  "Database type for output file: exodus. exodus is the default.",
		  "exodus");

  options_.enroll("extract_group", Ioss::GetLongOption::MandatoryValue,
		  "Write the data from the specified group to the output file.\n"
		  "\t\tUse 'ALL' to extract all groups in the file to separate output files.",
		  NULL);

  options_.enroll("64-bit", Ioss::GetLongOption::NoValue,
		  "Use 64-bit integers on output database",
		  NULL);

  options_.enroll("float", Ioss::GetLongOption::NoValue,
		  "Use 32-bit floating point values on output database; default is 64-bits",
		  NULL);

  options_.enroll("netcdf4", Ioss::GetLongOption::NoValue,
		  "Output database will be a netcdf4 hdf5-based file instead of the classical netcdf file format",
		  NULL);
  
  options_.enroll("shuffle", Ioss::GetLongOption::NoValue,
		  "Use a netcdf4 hdf5-based file and use hdf5s shuffle mode with compression.",
		  NULL);
  
  options_.enroll("compress", Ioss::GetLongOption::MandatoryValue,
		  "Specify the hdf5 compression level [0..9] to be used on the output file.",
		  NULL);
  
  options_.enroll("compose", Ioss::GetLongOption::MandatoryValue,
		  "Specify the parallel-io method to be used to output a single file in a parallel run. "
		  "Options are default, mpiio, mpiposix, pnetcdf",
		  NULL);

  options_.enroll("rcb", Ioss::GetLongOption::NoValue,
		  "Use recursive coordinate bisection method to decompose the input mesh in a parallel run.",
		  NULL);
  options_.enroll("rib", Ioss::GetLongOption::NoValue,
		  "Use recursive inertial bisection method to decompose the input mesh in a parallel run.",
		  NULL);

  options_.enroll("hsfc", Ioss::GetLongOption::NoValue,
		  "Use hilbert space-filling curve method to decompose the input mesh in a parallel run.",
		  NULL);

  options_.enroll("metis_sfc", Ioss::GetLongOption::NoValue,
		  "Use the metis space-filling-curve method to decompose the input mesh in a parallel run.",
		  NULL);
  
  options_.enroll("kway", Ioss::GetLongOption::NoValue,
		  "Use the metis kway graph-based method to decompose the input mesh in a parallel run.",
		  NULL);

  options_.enroll("kway_geom", Ioss::GetLongOption::NoValue,
		  "Use the metis kway graph-based method with geometry speedup to decompose the input mesh in a parallel run.",
		  NULL);

  options_.enroll("linear", Ioss::GetLongOption::NoValue,
		  "Use the linear method to decompose the input mesh in a parallel run. "
		  "elements in order first n/p to proc 0, next to proc 1.",
		  NULL);

  options_.enroll("cyclic", Ioss::GetLongOption::NoValue,
		  "Use the cyclic method to decompose the input mesh in a parallel run. "
		  "elements handed out to id % proc_count",
		  NULL);

  options_.enroll("random", Ioss::GetLongOption::NoValue,
		  "Use the random method to decompose the input mesh in a parallel run."
		  "elements assigned randomly to processors in a way that preserves balance (do not use for a real run)",
		  NULL);

  options_.enroll("external", Ioss::GetLongOption::NoValue,
		  "Files are decomposed externally into a file-per-processor in a parallel run.",
		  NULL);

  options_.enroll("debug" , Ioss::GetLongOption::NoValue,
		  "turn on debugging output",
		  NULL);

  options_.enroll("statistics" , Ioss::GetLongOption::NoValue,
		  "output parallel io timing statistics",
		  NULL);

  options_.enroll("Maximum_Time", Ioss::GetLongOption::MandatoryValue,
		  "Maximum time on input database to transfer to output database",
		  NULL);

  options_.enroll("Minimum_Time", Ioss::GetLongOption::MandatoryValue,
		  "Minimum time on input database to transfer to output database",
		  NULL);

  options_.enroll("field_suffix_separator", Ioss::GetLongOption::MandatoryValue,
		  "Character used to separate a field suffix from the field basename\n"
		  "\t\t when recognizing vector, tensor fields. Enter '0' for no separator", "_");

  options_.enroll("surface_split_scheme", Ioss::GetLongOption::MandatoryValue,
		  "Method used to split sidesets into homogenous blocks\n"
		  "\t\tOptions are: TOPOLOGY, BLOCK, NOSPLIT",
		  "TOPOLOGY");

  options_.enroll("copyright", Ioss::GetLongOption::NoValue,
		  "Show copyright and license data.",
		  NULL);
}

bool IOShell::Interface::parse_options(int argc, char **argv)
{
#if (__SUNPRO_CC == 0x500)
  using namespace std;
#endif

  // Get options from environment variable also...
  char *options = getenv("IO_SHELL_OPTIONS");
  if (options != NULL) {
    std::cerr << "\nThe following options were specified via the IO_SHELL_OPTIONS environment variable:\n"
	      << "\t" << options << "\n\n";
    options_.parse(options, options_.basename(*argv));
  }

  int option_index = options_.parse(argc, argv);
  if ( option_index < 1 )
    return false;

  if (options_.retrieve("help")) {
    options_.usage();
    std::cerr << "\n\tCan also set options via IO_SHELL_OPTIONS environment variable.\n\n";
    std::cerr << "\n\t->->-> Send email to gdsjaar@sandia.gov for io_shell support.<-<-<-\n";
    exit(EXIT_SUCCESS);
  }

  if (options_.retrieve("version")) {
    // Version is printed up front, just exit...
    exit(0);
  }
  
  if (options_.retrieve("64-bit")) {
    ints_64_bit = true;
  }

  if (options_.retrieve("float")) {
    reals_32_bit = true;
  }

  if (options_.retrieve("netcdf4")) {
    netcdf4 = true;
  }

  if (options_.retrieve("shuffle")) {
    shuffle = true;
  }

  if (options_.retrieve("rcb")) {
    decomp_method = "RCB";
  }

  if (options_.retrieve("rib")) {
    decomp_method = "RIB";
  }

  if (options_.retrieve("hsfc")) {
    decomp_method = "HSFC";
  }

  if (options_.retrieve("metis_sfc")) {
    decomp_method = "METIS_SFC";
  }

  if (options_.retrieve("kway")) {
    decomp_method = "KWAY";
  }

  if (options_.retrieve("kway_geom")) {
    decomp_method = "KWAY_GEOM";
  }

  if (options_.retrieve("linear")) {
    decomp_method = "LINEAR";
  }

  if (options_.retrieve("cyclic")) {
    decomp_method = "CYCLIC";
  }

  if (options_.retrieve("random")) {
    decomp_method = "RANDOM";
  }

  if (options_.retrieve("external")) {
    decomp_method = "EXTERNAL";
  }

  if (options_.retrieve("debug")) {
    debug = true;
  }

  if (options_.retrieve("statistics")) {
    statistics = true;
  }

  {
    const char *temp = options_.retrieve("in_type");
    if (temp != NULL) {
      inFiletype = temp;
    }
  }

  {
    const char *temp = options_.retrieve("out_type");
    if (temp != NULL) {
      outFiletype = temp;
    }
  }

  {
    const char *temp = options_.retrieve("compose");
    if (temp != NULL) {
      compose_output = Ioss::Utils::lowercase(temp);
    }
  }

  {
    const char *temp = options_.retrieve("extract_group");
    if (temp != NULL) {
      groupName = temp;
    }
  }

  {
    const char *temp = options_.retrieve("field_suffix_separator");
    if (temp != NULL) {
      fieldSuffixSeparator = temp[0];
    }
  }

  {
    const char *temp = options_.retrieve("surface_split_scheme");
    if (temp != NULL) {
      if (std::strcmp(temp, "TOPOLOGY") == 0)
        surface_split_type = 1;
      else if (std::strcmp(temp, "ELEMENT_BLOCK") == 0)
        surface_split_type = 2;
      else if (std::strcmp(temp, "NO_SPLIT") == 0)
        surface_split_type = 3;
    }
  }

  {
    const char *temp = options_.retrieve("Maximum_Time");
    if (temp != NULL) {
      maximum_time = std::strtod(temp, NULL);
    }
  }

  {
    const char *temp = options_.retrieve("Minimum_Time");
    if (temp != NULL) {
      minimum_time = std::strtod(temp, NULL);
    }
  }

  if (options_.retrieve("copyright")) {
    std::cerr << "\n"
	      << "Copyright(C) 2013 Sandia Corporation.  Under the terms of Contract\n"
	      << "DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains\n"
	      << "certain rights in this software\n"
	      << "\n"
	      << "Redistribution and use in source and binary forms, with or without\n"
	      << "modification, are permitted provided that the following conditions are\n"
	      << "met:\n"
	      << "\n"
	      << "    * Redistributions of source code must retain the above copyright\n"
	      << "      notice, this list of conditions and the following disclaimer.\n"
	      << "\n"
	      << "    * Redistributions in binary form must reproduce the above\n"
	      << "      copyright notice, this list of conditions and the following\n"
	      << "      disclaimer in the documentation and/or other materials provided\n"
	      << "      with the distribution.\n"
	      << "\n"
	      << "    * Neither the name of Sandia Corporation nor the names of its\n"
	      << "      contributors may be used to endorse or promote products derived\n"
	      << "      from this software without specific prior written permission.\n"
	      << "\n"
	      << "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
	      << "'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
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
  if (option_index < argc-1) {
    while (option_index < argc-1) {
      inputFile.push_back(argv[option_index++]);
    }
    outputFile = argv[option_index];
  } else {
    std::cerr << "\nERROR: input and output filename not specified\n\n";
    return false;
  }
  return true;
}

void IOShell::Interface::dump(std::ostream &) const
{
}

void IOShell::Interface::show_version()
{
  std::cout << "\tExodusII Information Program\n";
}

