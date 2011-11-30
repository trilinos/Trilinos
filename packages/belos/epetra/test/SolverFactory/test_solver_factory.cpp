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

/// \file test_solver_factory.cpp
/// \brief Test Belos::SolverFactory with Epetra
/// \author Mark Hoemmen
///
#include "BelosConfigDefs.hpp"
#include "BelosOutputManager.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosEpetraAdapter.hpp"

#include <Galeri_Maps.h>
#include <Galeri_CrsMatrices.h>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

#ifdef EPETRA_MPI
#  include <mpi.h>
#  include <Epetra_MpiComm.h>
#else
#  include <Epetra_SerialComm.h>
#endif
#include <Epetra_CrsMatrix.h>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using std::cout;
using std::endl;
using std::vector;

namespace {

  int
  selectVerbosity (const bool verbose, const bool debug)
  {
    // NOTE Calling this a "MsgType" (its correct type) or even an
    // "enum MsgType" confuses the compiler.
    int theType = Belos::Errors; // default (always print errors)
    if (verbose) 
      {
	// "Verbose" also means printing out Debug messages (as well
	// as everything else).
	theType = theType | 
	  Belos::Warnings | 
	  Belos::IterationDetails |
	  Belos::OrthoDetails | 
	  Belos::FinalSummary | 
	  Belos::TimingDetails |
	  Belos::StatusTestDetails | 
	  Belos::Debug;
      }
    if (debug)
      // "Debug" doesn't necessarily mean the same thing as
      // "Verbose".  We interpret "Debug" to mean printing out
      // messages marked as Debug (as well as Error messages).
      theType = theType | Belos::Debug;
    return theType;
  }

  /// The problem is defined on a 2D grid; global size is nx * nx.
  Teuchos::RCP<Epetra_Operator>
  makeMatrix (const Teuchos::RCP<Epetra_Comm>& comm,
	      Teuchos::RCP<const Epetra_Map>& domainMap,
	      Teuchos::RCP<const Epetra_Map>& rangeMap,
	      const int nx = 30)
  {
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_implicit_cast;

    ParameterList GaleriList;
    GaleriList.set("n", nx * nx);
    GaleriList.set("nx", nx);
    GaleriList.set("ny", nx);
    RCP<const Epetra_Map> rowMap = rcp (Galeri::CreateMap("Linear", *comm, GaleriList));
    // "&*rowMap" turns an RCP<Epetra_Map> into a raw pointer
    // (Epetra_Map*), which is what Galeri::CreateCrsMatrix() wants.
    RCP<Epetra_RowMatrix> A = 
      rcp (Galeri::CreateCrsMatrix("Laplace2D", &*rowMap, GaleriList));
    TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(), std::runtime_error,
			       "Galeri returned a null operator A.");
    domainMap = rowMap;
    rangeMap = rowMap;
    return rcp_implicit_cast<Epetra_Operator> (A);
  }

} // namespace (anonymous)


int 
main (int argc, char *argv[]) 
{
  using Belos::OutputManager;
  using Belos::SolverFactory;
  using Teuchos::CommandLineProcessor;
  using Teuchos::null;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_implicit_cast;
  using std::cerr;
  using std::cout;
  using std::endl;
  typedef double scalar_type;
  typedef Epetra_MultiVector MV;
  typedef Epetra_Operator OP;
  typedef Belos::SolverManager<scalar_type, MV, OP> solver_base_type;
  typedef Belos::PseudoBlockGmresSolMgr<scalar_type, MV, OP> solver_impl_type;

  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);

  RCP<Epetra_Comm> comm;
  {
#ifdef EPETRA_MPI
    RCP<Epetra_MpiComm> commSpecific (new Epetra_MpiComm (MPI_COMM_WORLD));
#else
    RCP<Epetra_SerialComm> commSpecific (new Epetra_SerialComm);
#endif // EPETRA_MPI
    comm = rcp_implicit_cast<Epetra_Comm> (commSpecific);
  }

  int numRHS = 1;
  bool verbose = false;
  bool debug = false;

  // Define command-line arguments.
  CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("numRHS", &numRHS, "Number of right-hand sides in the linear "
		  "system to solve.");
  cmdp.setOption ("verbose", "quiet", &verbose, "Print messages and results.");
  cmdp.setOption ("debug", "nodebug", &debug, "Print debugging information.");

  // Parse the command-line arguments.
  {
    const CommandLineProcessor::EParseCommandLineReturn parseResult = 
      cmdp.parse (argc,argv);
    if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED) {
      if (comm->MyPID() == 0)
	std::cout << "End Result: TEST PASSED" << endl;
      return EXIT_SUCCESS;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL, 
      std::invalid_argument, "Failed to parse command-line arguments.");
  }

  // Declare an output manager for handling local output.  Initialize,
  // using the caller's desired verbosity level.
  RCP<OutputManager<scalar_type> > outMan = 
    rcp (new OutputManager<scalar_type> (selectVerbosity (verbose, debug)));

  // Stream for debug output.  If debug output is not enabled, then
  // this stream doesn't print anything sent to it (it's a "black
  // hole" stream).
  //std::ostream& debugOut = outMan->stream (Belos::Debug);

  // Create the operator to test, with domain and range maps.
  RCP<const Epetra_Map> domainMap, rangeMap;
  RCP<Epetra_Operator> A = makeMatrix (comm, domainMap, rangeMap);
  // "Solution" input/output multivector.
  RCP<MV> X_exact = rcp (new MV (*domainMap, numRHS));
  X_exact->Seed ();
  X_exact->Random ();
  RCP<MV> X = rcp (new MV (*domainMap, numRHS));
  X->PutScalar (0.0);
  // "Right-hand side" input multivector.
  RCP<MV> B = rcp (new MV (*rangeMap, numRHS));
  A->Apply (*X_exact, *B);
  
  // Solver parameters.
  RCP<ParameterList> solverParams = parameterList ();
  //
  // Create a solver instance using a solver factory.
  //
  const std::string solverName ("Pseudoblock GMRES");
  const std::string solverType ("PseudoBlockGmresSolMgr");
  SolverFactory<scalar_type, MV, OP> factory;
  RCP<solver_base_type> solver = factory.create (solverName, solverParams);
  TEUCHOS_TEST_FOR_EXCEPTION(solver.is_null(), std::logic_error,
			     "Failed to create solver with valid name \"" 
			     << solverName << "\".");
  if (verbose && comm->MyPID() == 0) {
    cout << "Solver description: " << solver->description() << endl;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(rcp_dynamic_cast<solver_impl_type> (solver).is_null(), 
    std::logic_error, 
    "Solver does not have the correct type \"" << solverType << "\".");
  if (comm->MyPID() == 0) {
    cout << "End Result: TEST PASSED" << endl;
  }
  return EXIT_SUCCESS;
}


