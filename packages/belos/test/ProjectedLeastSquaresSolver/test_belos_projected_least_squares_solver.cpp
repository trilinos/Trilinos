//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

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
