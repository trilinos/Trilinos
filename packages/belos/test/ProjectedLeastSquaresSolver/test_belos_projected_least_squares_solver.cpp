// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <BelosProjectedLeastSquaresSolver.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_ParameterList.hpp>

int 
main (int argc, char *argv[]) 
{
  using Belos::details::ERobustness;
  using Belos::details::ProjectedLeastSquaresSolver;
  using Belos::details::robustnessEnumToString;
  using Belos::details::robustnessStringToEnum;
  using Teuchos::CommandLineProcessor;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;

  typedef double scalar_type;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  Teuchos::oblackholestream blackHole;
  // Initialize MPI using Teuchos wrappers, if Trilinos was built with
  // MPI support.  Otherwise, initialize a communicator with one
  // process.
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
  const int myRank = mpiSession.getRank();
  // Output stream only prints on MPI Proc 0.
  std::ostream& out = (myRank == 0) ? std::cout : blackHole;
  
  // Command-line arguments
  std::string robustnessLevel ("None");
  bool testBlockGivens = false;
  bool testGivensRotations = false;
  bool verbose = false;
  bool debug = false;
  int testProblemSize = 10;

  // Parse command-line arguments
  CommandLineProcessor cmdp (false,true);
  cmdp.setOption ("robustness", &robustnessLevel,
		  "Robustness level: \"None\", \"Some\", or \"Lots\".");
  cmdp.setOption ("testGivensRotations", "dontTestGivensRotations", 
		  &testGivensRotations, 
		  "Test the implementation of Givens rotations.");
  cmdp.setOption ("testBlockGivens", "dontTestBlockGivens", &testBlockGivens,
		  "Test the panel version of the Givens rotations - based "
		  "update.");
  cmdp.setOption ("verbose", "quiet", &verbose, "Print messages and results.");
  cmdp.setOption ("debug", "release", &debug, "Print copious debug output.");
  cmdp.setOption ("testProblemSize", &testProblemSize, 
		  "Number of columns in the projected least-squares test "
		  "problem.");
  if (cmdp.parse (argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    out << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
  // Verbose output stream only prints on MPI Proc 0, and only in
  // verbose mode.
  std::ostream& verboseOut = verbose ? out : blackHole;

  // Robustness level for triangular solves.
  const ERobustness robustness = robustnessStringToEnum (robustnessLevel);
  verboseOut << "-- Robustness level: " << robustnessLevel << endl;

  bool success = true; // Innocent until proven guilty.

  // Seed the pseudorandom number generator with the same seed each
  // time, to ensure repeatable results.
  STS::seedrandom (0);
  ProjectedLeastSquaresSolver<scalar_type> solver (out, robustness);

  if (testGivensRotations) {
    solver.testGivensRotations (verboseOut);
  }
  if (testProblemSize > 0) {
    const bool extraVerbose = debug;
    success = success && 
      solver.testUpdateColumn (verboseOut, testProblemSize, 
			       testBlockGivens, extraVerbose);
    success = success &&
      solver.testTriangularSolves (verboseOut, testProblemSize, 
				   robustness, extraVerbose);
  }

  if (success) {
    out << "End Result: TEST PASSED" << endl;  
    return EXIT_SUCCESS;
  } else {
    out << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
}
