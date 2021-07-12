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
#include <BelosTpetraTestFramework.hpp>
#include <BelosBlockGCRODRSolMgr.hpp>

#include <Tpetra_Core.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

int
main (int argc, char *argv[])
{
  using Teuchos::Comm;
  using Teuchos::FancyOStream;
  using Teuchos::getFancyOStream;
  using Teuchos::oblackholestream;
  using Teuchos::OSTab;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  using std::cout;
  using std::endl;
  //
  // Typedefs for Tpetra template arguments.
  //
  typedef double scalar_type;
  typedef long int global_ordinal_type;
  typedef int local_ordinal_type;
  //
  // Tpetra objects which are the MV and OP template parameters of the
  // Belos specialization which we are testing.
  //
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type> MV;
  typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type> OP;
  //
  // Other typedefs.
  //
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type> sparse_matrix_type;
  typedef MV::node_type node_type;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &cout);
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  RCP<oblackholestream> blackHole (new oblackholestream);
  const int myRank = comm->getRank ();

  // Output stream that prints only on Process 0.
  RCP<FancyOStream> out;
  if (myRank == 0) {
    out = Teuchos::getFancyOStream (rcpFromRef (cout));
  } else {
    out = Teuchos::getFancyOStream (blackHole);
  }

  //
  // Get test parameters from command-line processor.
  //
  // CommandLineProcessor always understands int, but may not
  // understand global_ordinal_type.  We convert to the latter below.
  int numRows = comm->getSize() * 100;
  bool tolerant = false;
  bool verbose = false;
  bool debug = false;
  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption("numRows", &numRows,
                 "Global number of rows (and columns) in the sparse matrix to generate.");
  cmdp.setOption("tolerant", "intolerant", &tolerant,
                 "Whether to parse files tolerantly.");
  cmdp.setOption("verbose", "quiet", &verbose,
                 "Print messages and results.");
  cmdp.setOption("debug", "release", &debug,
                 "Run debugging checks and print copious debugging output.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    *out << "\nEnd Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
  // Output stream for verbose output.
  RCP<FancyOStream> verbOut = verbose ? out : getFancyOStream (blackHole);

  const bool success = true;

  // Test whether it's possible to instantiate the solver.
  // This is a minimal compilation test.
  *verbOut << "Instantiating Block GCRODR solver" << endl;
  Belos::BlockGCRODRSolMgr<scalar_type, MV, OP> solver;
  //
  // Test setting solver parameters.  For now, we just use an empty
  // (but non-null) parameter list, which the solver should fill in
  // with defaults.
  //
  *verbOut << "Setting solver parameters" << endl;
  RCP<ParameterList> solverParams = parameterList ();
  solver.setParameters (solverParams);
  //
  // Create a linear system to solve.
  //
  *verbOut << "Creating linear system" << endl;
  RCP<sparse_matrix_type> A;
  RCP<MV> X_guess, X_exact, B;
  {
    Teuchos::RCP<node_type> node; // can be null; only for type deduction
    typedef Belos::Tpetra::ProblemMaker<sparse_matrix_type> factory_type;
    factory_type factory (comm, node, out, tolerant, debug);

    RCP<ParameterList> problemParams = parameterList ();
    problemParams->set ("Global number of rows",
                        static_cast<global_ordinal_type> (numRows));
    problemParams->set ("Problem type", std::string ("Nonsymmetric"));
    factory.makeProblem (A, X_guess, X_exact, B, problemParams);
  }
  // Approximate solution vector is a copy of the guess vector.
  RCP<MV> X (new MV (*X_guess));

  TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(), std::logic_error,
                             "The sparse matrix is null!");
  TEUCHOS_TEST_FOR_EXCEPTION(X_guess.is_null(), std::logic_error,
                             "The initial guess X_guess is null!");
  TEUCHOS_TEST_FOR_EXCEPTION(X_exact.is_null(), std::logic_error,
                             "The exact solution X_exact is null!");
  TEUCHOS_TEST_FOR_EXCEPTION(B.is_null(), std::logic_error,
                             "The right-hand side B is null!");
  TEUCHOS_TEST_FOR_EXCEPTION(X.is_null(), std::logic_error,
                             "The approximate solution vector X is null!");

  typedef Belos::LinearProblem<scalar_type, MV, OP> problem_type;
  RCP<problem_type> problem (new problem_type (A, X, B));
  problem->setProblem ();
  solver.setProblem (problem);

  *verbOut << "Solving linear system" << endl;
  Belos::ReturnType result = solver.solve ();

  *verbOut << "Result of solve: "
           << Belos::convertReturnTypeToString (result)
           << endl;

  // Make sure that all the processes finished.
  comm->barrier ();

  if (success) {
    *out << "\nEnd Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  }
  else {
    *out << "\nEnd Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
}


