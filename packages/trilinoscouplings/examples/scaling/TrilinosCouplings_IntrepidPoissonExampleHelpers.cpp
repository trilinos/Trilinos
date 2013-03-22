// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

#include "TrilinosCouplings_IntrepidPoissonExampleHelpers.hpp"
#include <iostream>
#include <sstream>


namespace TrilinosCouplings {
namespace IntrepidPoissonExample {

std::string
makeMeshInput (const int nx, const int ny, const int nz)
{
  using std::endl;
  std::ostringstream os;

  TEUCHOS_TEST_FOR_EXCEPTION( nx <= 0 || ny <= 0 || nz <= 0,
    std::invalid_argument, "nx, ny, and nz must all be positive.");

  os << "mesh" << endl
     << "\trectilinear" << endl
     << "\t\tnx = " << nx << endl
     << "\t\tny = " << ny << endl
     << "\t\tnz = " << nz << endl
     << "\t\tbx = 1" << endl
     << "\t\tby = 1" << endl
     << "\t\tbz = 1" << endl
     << "\t\tgmin = 0 0 0" << endl
     << "\t\tgmax = 1 1 1" << endl
     << "\tend" << endl
     << "\tset assign" << endl
     << "\t\tsideset, ilo, 1" << endl
     << "\t\tsideset, jlo, 2" << endl
     << "\t\tsideset, klo, 3" << endl
     << "\t\tsideset, ihi, 4" << endl
     << "\t\tsideset, jhi, 5" << endl
     << "\t\tsideset, khi, 6" << endl
     << "\tend" << endl
     << "end";
  return os.str ();
}

void
setCommandLineArgumentDefaults (int& nx,
                                int& ny,
                                int& nz,
                                std::string& xmlInputParamsFile,
                                bool& verbose,
                                bool& debug)
{
  nx = 20;
  ny = 20;
  nz = 20;
  xmlInputParamsFile = "";
  verbose = false;
  debug = false;
}

void
setUpCommandLineArguments (Teuchos::CommandLineProcessor& cmdp,
                           int& nx,
                           int& ny,
                           int& nz,
                           std::string& xmlInputParamsFile,
                           bool& verbose,
                           bool& debug)
{
  cmdp.setOption ("nx", &nx, "Number of cells along the x dimension");
  cmdp.setOption ("ny", &ny, "Number of cells along the y dimension");
  cmdp.setOption ("nz", &nz, "Number of cells along the z dimension");
  cmdp.setOption ("inputParams", &xmlInputParamsFile, "XML file of input "
                  "parameters, which we read if specified and not \"\".  "
                  "If it has a \"meshInput\" parameter, we use its "
                  "std::string value as the Pamgen mesh specification.  "
                  "Otherwise, we tell Pamgen to make a cube, using "
                  "nx, ny, and nz.");
  cmdp.setOption ("verbose", "quiet", &verbose,
                  "Whether to print verbose status output.");
  cmdp.setOption ("debug", "release", &debug,
                  "Whether to print copious debugging output to stderr.");
}

void
parseCommandLineArguments (Teuchos::CommandLineProcessor& cmdp,
                           bool& printedHelp,
                           int argc,
                           char* argv[],
                           int& nx,
                           int& ny,
                           int& nz,
                           std::string& xmlInputParamsFile,
                           bool& verbose,
                           bool& debug)
{
  using Teuchos::CommandLineProcessor;

  const CommandLineProcessor::EParseCommandLineReturn parseResult =
    cmdp.parse (argc, argv);
  if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED) {
    printedHelp = true;
  }
  else {
    printedHelp = false;
    TEUCHOS_TEST_FOR_EXCEPTION(
      parseResult != CommandLineProcessor::PARSE_SUCCESSFUL,
      std::invalid_argument, "Failed to parse command-line arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      xmlInputParamsFile == "" && (nx <= 0 || ny <= 0 || nz <= 0),
      std::invalid_argument, "If no XML parameters filename is specified (via "
      "--inputParams), then the number of cells along each dimension of the "
      "mesh (--nx, --ny, and --nz) must be positive.");
  }
}

} // namespace IntrepidPoissonExample
} // namespace TrilinosCouplings
