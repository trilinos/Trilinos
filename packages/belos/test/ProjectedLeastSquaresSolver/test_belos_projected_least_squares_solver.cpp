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
  using Teuchos::CommandLineProcessor;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;
  typedef double scalar_type;

  Teuchos::oblackholestream blackHole;
  // MPI is crashing on my machine for some reason, but this test
  // doesn't depend on MPI.
#if 0
  // Initialize MPI using Teuchos wrappers, if Trilinos was built with
  // MPI support.  Otherwise, initialize a communicator with one
  // process.
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  const int myRank = comm->getRank();
  //const int numProcs = comm->getSize();
#else
  const int myRank = 0;
#endif // 0
  
  // Command-line arguments
  bool testBlockGivens = false;
  bool verbose = false;
  int testProblemSize = 10;

  // Parse command-line arguments
  CommandLineProcessor cmdp (false,true);
  cmdp.setOption ("testBlockGivens", "dontTestBlockGivens", &testBlockGivens, 
		  "Print messages and results.");
  cmdp.setOption ("verbose", "quiet", &verbose, "Print messages and results.");
  cmdp.setOption ("testProblemSize", &testProblemSize, 
		  "Number of columns in the projected least-squares test problem.");
  if (cmdp.parse (argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
    if (myRank == 0) {
      std::cout << "End Result: TEST FAILED" << std::endl;
    }
    return EXIT_FAILURE;
  }
  // Output stream only prints on MPI Proc 0, and only in verbose mode.
  std::ostream& out = (myRank == 0 && verbose) ? std::cout : blackHole;

  bool success = true;
  if (testProblemSize > 0) {
    // Test the projected least-squares solver.
    Belos::details::ProjectedLeastSquaresSolver<scalar_type> solver;
    success = solver.testUpdateColumn (out, testProblemSize, 
				       testBlockGivens, verbose);
  }
  
  if (success) {
    if (myRank == 0) {
      std::cout << "End Result: TEST PASSED" << std::endl;  
    } 
    return EXIT_SUCCESS;
  } else {
    if (myRank == 0) {
      std::cout << "End Result: TEST FAILED" << std::endl;
    }
    return EXIT_FAILURE;
  }
}
