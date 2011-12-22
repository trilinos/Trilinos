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

//
// Test for BlockGCRODRSolMgr (Block Recycling GMRES, by Kirk
// Soodhalter and Michael Parks).  This just tests compilation and
// setParameters() for now.  Later, we'll test actually solving linear
// systems.
// 
#include <BelosConfigDefs.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosBlockGCRODRSolMgr.hpp>

#include <Tpetra_DefaultPlatform.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

int 
main (int argc, char *argv[]) 
{
  using Teuchos::Comm;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using std::cout;
  using std::endl;
  typedef double scalar_type;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Tpetra::Operator<scalar_type, int> OP;
  typedef Tpetra::MultiVector<scalar_type, int> MV;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &cout);
  RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  Teuchos::oblackholestream blackHole;
  const int myRank = comm->getRank();
  std::ostream& out = (myRank == 0) ? cout : blackHole;

  //
  // Get test parameters from command-line processor
  //  
  // bool verbose = false;
  //bool debug = false;
  Teuchos::CommandLineProcessor cmdp (false, true);
  //cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  //cmdp.setOption("debug","nodebug",&debug,"Run debugging checks.");

  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    out << "\nEnd Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }

  const bool success = true;

  // Test whether it's possible to instantiate the solver.
  // This is a minimal compilation test.
  Belos::BlockGCRODRSolMgr<scalar_type, MV, OP> solver;
  //
  // Test setting solver parameters.  For now, we just use an empty
  // (but non-null) parameter list, which the solver should fill in
  // with defaults.
  //
  RCP<ParameterList> params = parameterList ();
  solver.setParameters (params);

  if (success) {
    out << "\nEnd Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  } 
  else {
    out << "\nEnd Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
}


