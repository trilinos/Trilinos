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
  using Belos::details::robustnessEnumToString;
  using Belos::details::robustnessStringToEnum;
  using Teuchos::CommandLineProcessor;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  typedef double scalar_type;
  typedef Teuchos::ScalarTraits<scalar_type> STS;

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

  const ERobustness robustness = 
    robustnessStringToEnum (robustnessLevel);

  // Verbose output stream only prints on MPI Proc 0, and only in
  // verbose mode.
  std::ostream& verboseOut = verbose ? out : blackHole;

  // Robustness level for triangular solves.
  verboseOut << "-- Robustness level: " << robustnessLevel << endl;

  bool success = true; // Innocent until proven guilty.

  STS::seedrandom (0);
  Belos::details::ProjectedLeastSquaresSolver<scalar_type> solver;

  if (testGivensRotations) {
    solver.testGivensRotations (verboseOut);
  }
  if (testProblemSize > 0) {
    const bool extraVerbose = debug;
    success = success && 
      solver.testUpdateColumn (verboseOut, testProblemSize, 
			       testBlockGivens, extraVerbose);
  }

  if (testProblemSize > 0) {
    verboseOut << "Testing upper triangular solves" << endl;
    //
    // Construct an upper triangular linear system to solve.
    //
    verboseOut << "-- Generating test matrix" << endl;
    const int N = testProblemSize;
    typedef Teuchos::SerialDenseMatrix<int, scalar_type> mat_type;
    mat_type R (N, N);
    // Fill the upper triangle of R with random numbers.
    for (int j = 0; j < N; ++j) {
      for (int i = 0; i <= j; ++i) {
	R(i,j) = STS::random ();
      }
    }
    mat_type B (N, 1);
    B.random ();

    // Save a copy of the original upper triangular system.
    mat_type R_copy (Teuchos::Copy, R, N, N);
    mat_type B_copy (Teuchos::Copy, B, N, 1);

    // Solution vector.
    mat_type X (N, 1);

    // Solve RX = B.
    verboseOut << "-- Solving RX=B" << endl;
    (void) solver.solveUpperTriangularSystem (Teuchos::LEFT_SIDE, X, R, B, 
					      robustness);
    // Test the residual error.
    mat_type Resid (N, 1);
    Resid.assign (B_copy);
    Belos::details::LocalDenseMatrixOps<scalar_type> ops;
    ops.matMatMult (STS::one(), Resid, -STS::one(), R_copy, X);
    verboseOut << "---- ||R*X - B||_F = " << Resid.normFrobenius() << endl;
    verboseOut << "---- ||R||_F ||X||_F + ||B||_F = " 
	<< (R_copy.normFrobenius() * X.normFrobenius() + B_copy.normFrobenius())
	<< endl;

    // Restore R and B.
    R.assign (R_copy);
    B.assign (B_copy);

    // 
    // Set up a right-side test problem: YR = B^*.
    //
    mat_type Y (1, N);
    mat_type B_star (1, N);
    ops.conjugateTranspose (B_star, B);
    mat_type B_star_copy (1, N);
    B_star_copy.assign (B_star);
    // Solve YR = B^*.
    verboseOut << "-- Solving YR=B^*" << endl;
    (void) solver.solveUpperTriangularSystem (Teuchos::RIGHT_SIDE, Y, R, B_star, 
					      robustness);
    // Test the residual error.
    mat_type Resid2 (1, N);
    Resid2.assign (B_star_copy);
    ops.matMatMult (STS::one(), Resid2, -STS::one(), Y, R_copy);
    verboseOut << "---- ||Y*R - B^*||_F = " << Resid2.normFrobenius() << endl;
    verboseOut << "---- ||Y||_F ||R||_F + ||B^*||_F = " 
	       << (Y.normFrobenius() * R_copy.normFrobenius() + B_star_copy.normFrobenius())
	       << endl;
  }
  
  if (success) {
    out << "End Result: TEST PASSED" << endl;  
    return EXIT_SUCCESS;
  } else {
    out << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
}
