/** \file  example_IntrepidPoisson_Pamgen_EpetraAztecOO.cpp
    \brief Example setup of a discretization of a Poisson equation on
           a hexahedral mesh using nodal (Hgrad) elements.  The
           system is assembled but not solved.

           This example uses the following Trilinos packages:
    \li    Pamgen to generate a Hexahedral mesh.
    \li    Sacado to form the source term from user-specified manufactured solution.
    \li    Intrepid to build the discretization matrix and right-hand side.
    \li    Tpetra to handle the global sparse matrix and dense vector.

    \verbatim

     Poisson system:

            div A grad u = f in Omega
                       u = g on Gamma

      where
             A is a material tensor (typically symmetric positive definite)
             f is a given source term

     Corresponding discrete linear system for nodal coefficients(x):

                 Kx = b

            K - HGrad stiffness matrix
            b - right hand side vector

    \endverbatim

    \author Created by P. Bochev, D. Ridzal, K. Peterson,
            D. Hensinger, C. Siefert.  Converted to Tpetra by
            I. Kalashnikova (ikalash@sandia.gov).  Modified by Mark
            Hoemmen (mhoemme@sandia.gov) and back-ported to Epetra for
            a fair comparison with Tpetra.

    \remark Use the "--help" command-line argument for usage info.

    \remark Example driver has an option to use an input file (XML
            serialization of a Teuchos::ParameterList) containing a
            Pamgen mesh description.  A version, Poisson.xml, is
            included in the same directory as this driver.  If not
            included, this program will use a default mesh
            description.

    \remark The exact solution (u) and material tensor (A) are set in
            the functions "exactSolution" and "materialTensor" and may
            be modified by the user.  We compute the source term f
            from u using Sacado automatic differentiation, so that u
            really is the exact solution.  The current implementation
            of exactSolution() has notes to guide your choice of
            solution.  For example, you might want to pick a solution
            in the finite element space, so that the discrete solution
            is exact if the solution of the linear system is exact
            (modulo rounding error when constructing the linear
            system).
*/

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include "Epetra_Comm.h"
#ifdef EPETRA_MPI
#  include "mpi.h"
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif // EPETRA_MPI

#include "TrilinosCouplings_EpetraIntrepidPoissonExample.hpp"
#include "TrilinosCouplings_IntrepidPoissonExampleHelpers.hpp"
#include "TrilinosCouplings_EpetraIntrepidPoissonExample_SolveWithAztecOO.hpp"

// ML
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"

// MueLu includes
#include "MueLu.hpp"
#include "MueLu_ParameterListInterpreter.hpp"
#include "MueLu_EpetraOperator.hpp"

int
main (int argc, char *argv[])
{
  using namespace TrilinosCouplings; // Yes, this means I'm lazy.

  using EpetraIntrepidPoissonExample::exactResidualNorm;
  using EpetraIntrepidPoissonExample::makeMatrixAndRightHandSide;
  using EpetraIntrepidPoissonExample::solveWithAztecOO;
  using IntrepidPoissonExample::makeMeshInput;
  using IntrepidPoissonExample::setCommandLineArgumentDefaults;
  using IntrepidPoissonExample::setUpCommandLineArguments;
  using IntrepidPoissonExample::parseCommandLineArguments;
  //using Tpetra::DefaultPlatform;
  //using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::rcpFromRef;
  using Teuchos::getFancyOStream;
  using Teuchos::FancyOStream;
  using std::endl;
  // Pull in typedefs from the example's namespace.
  typedef EpetraIntrepidPoissonExample::ST ST;
  //typedef EpetraIntrepidPoissonExample::Node Node;
  typedef Teuchos::ScalarTraits<ST> STS;
  typedef STS::magnitudeType MT;
  typedef Teuchos::ScalarTraits<MT> STM;
  typedef EpetraIntrepidPoissonExample::sparse_matrix_type sparse_matrix_type;
  typedef EpetraIntrepidPoissonExample::vector_type vector_type;
  typedef EpetraIntrepidPoissonExample::operator_type operator_type;

  bool success = true;
  try {

  Epetra_Object::SetTracebackMode(2);

  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
  const int myRank = mpiSession.getRank ();
  //const int numProcs = mpiSession.getNProc ();

  // Get the default communicator and Kokkos Node instance
  RCP<Epetra_Comm> comm;
#ifdef EPETRA_MPI
  comm = rcp_implicit_cast<Epetra_Comm> (rcp (new Epetra_MpiComm (MPI_COMM_WORLD)));
#else
  comm = rcp_implicit_cast<Epetra_Comm> (rcp (new Epetra_SerialComm));
#endif // EPETRA_MPI

  // Did the user specify --help at the command line to print help
  // with command-line arguments?
  bool printedHelp = false;
  // Values of command-line arguments.
  int nx, ny, nz;
  std::string xmlInputParamsFile;
  bool verbose, debug;

  // Set default values of command-line arguments.
  setCommandLineArgumentDefaults (nx, ny, nz, xmlInputParamsFile,
                                  verbose, debug);
  // Parse and validate command-line arguments.
  Teuchos::CommandLineProcessor cmdp (false, true);
  setUpCommandLineArguments (cmdp, nx, ny, nz, xmlInputParamsFile,
                             verbose, debug);
  parseCommandLineArguments (cmdp, printedHelp, argc, argv, nx, ny, nz,
                             xmlInputParamsFile, verbose, debug);
  if (printedHelp) {
    // The user specified --help at the command line to print help
    // with command-line arguments.  We printed help already, so quit
    // with a happy return code.
    return EXIT_SUCCESS;
  }

  // Both streams only print on MPI Rank 0.  "out" only prints if the
  // user specified --verbose.
  RCP<FancyOStream> out =
    getFancyOStream (rcpFromRef ((myRank == 0 && verbose) ? std::cout : blackHole));
  RCP<FancyOStream> err =
    getFancyOStream (rcpFromRef ((myRank == 0 && debug) ? std::cerr : blackHole));

#ifdef HAVE_MPI
  *out << "PARALLEL executable" << endl;
#else
  *out << "SERIAL executable" << endl;
#endif

/**********************************************************************************/
/********************************** GET XML INPUTS ********************************/
/**********************************************************************************/
  ParameterList inputList;
  if (xmlInputParamsFile != "") {
    *out << "Reading parameters from XML file \""
         << xmlInputParamsFile << "\"..." << endl;
    Teuchos::updateParametersFromXmlFile (xmlInputParamsFile, 
					  outArg (inputList));
    if (myRank == 0) {
      inputList.print (*out, 2, true, true);
      *out << endl;
    }
  }

  // Get Pamgen mesh definition string, either from the input
  // ParameterList or from our function that makes a cube and fills in
  // the number of cells along each dimension.
  std::string meshInput = inputList.get("meshInput", "");
  if (meshInput == "") {
    *out << "Generating mesh input string: nx = " << nx
         << ", ny = " << ny
         << ", nz = " << nz << endl;
    meshInput = makeMeshInput (nx, ny, nz);
  }

  // Total application run time
  {
  TEUCHOS_FUNC_TIME_MONITOR_DIFF("Total Time", total_time);

  // Construct linear system
  RCP<sparse_matrix_type> A;
  RCP<vector_type> B, X_exact, X;
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF("Total Assembly", total_assembly);
    makeMatrixAndRightHandSide (A, B, X_exact, X, comm, meshInput,
				out, err, verbose, debug);
  }

  const std::vector<MT> norms = exactResidualNorm (A, B, X_exact);
  // X_exact is the exact solution of the PDE, projected onto the
  // discrete mesh.  It may not necessarily equal the exact solution
  // of the linear system.
  *out << "||B - A*X_exact||_2 = " << norms[0] << endl
       << "||B||_2 = " << norms[1] << endl
       << "||A||_F = " << norms[2] << endl;

  // Setup preconditioner
  std::string prec_type = inputList.get("Preconditioner", "None");
  RCP<operator_type> M;
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF("Total Preconditioner Setup", total_prec);

    if (prec_type == "ML") {
      ParameterList mlParams;
      if (inputList.isSublist("ML"))
	mlParams = inputList.sublist("ML");
      else {
	ML_Epetra::SetDefaults("SA", mlParams);
	mlParams.set("ML output", 0);
      }
      M = rcp(new ML_Epetra::MultiLevelPreconditioner(*A, mlParams));
    }

    else if (prec_type == "MueLu") {
      // Turns a Epetra_CrsMatrix into a MueLu::Matrix
      RCP<Xpetra::CrsMatrix<ST> > mueluA_ = 
	rcp(new Xpetra::EpetraCrsMatrix(A));
      RCP<Xpetra::Matrix <ST> > mueluA  = 
	rcp(new Xpetra::CrsMatrixWrap<ST>(mueluA_));
      
      // Multigrid Hierarchy
      ParameterList mueluParams;
      if (inputList.isSublist("MueLu"))
	mueluParams = inputList.sublist("MueLu");
      MueLu::ParameterListInterpreter<ST> mueLuFactory(mueluParams);
      RCP<MueLu::Hierarchy<ST> > H = 
	mueLuFactory.CreateHierarchy();
      H->setVerbLevel(Teuchos::VERB_HIGH);
      H->GetLevel(0)->Set("A", mueluA);
      
      // Multigrid setup phase
      H->Setup();

      // Wrap MueLu Hierarchy as a Tpetra::Operator
      M = rcp(new MueLu::EpetraOperator(H));
    }
  }

  bool converged = false;
  int numItersPerformed = 0;
  const MT tol = inputList.get("Convergence Tolerance",
			       STM::squareroot (STM::eps ()));
  const int maxNumIters = inputList.get("Maximum Iterations", 200);
  {
    TEUCHOS_FUNC_TIME_MONITOR_DIFF("Total Solve", total_solve);
    solveWithAztecOO (converged, numItersPerformed, tol, maxNumIters,
		      X, A, B, M);
  }

  // Compute ||X-X_exact||_2
  MT norm_x, norm_error;
  X_exact->Norm2(&norm_x);
  X_exact->Update(-1.0, *X, 1.0);
  X_exact->Norm2(&norm_error);
  *out << endl
       << "||X-X_exact||_2 / ||X_exact||_2 = " << norm_error / norm_x 
       << endl;

  } // total time block

  // Summarize timings
  // RCP<ParameterList> reportParams = parameterList ("TimeMonitor::report");
  // reportParams->set ("Report format", std::string ("YAML"));
  // reportParams->set ("writeGlobalStats", true);
  // Teuchos::TimeMonitor::report (*out, reportParams);
  Teuchos::TimeMonitor::summarize(std::cout);

  } //try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  if (success)
    return EXIT_SUCCESS;
  return EXIT_FAILURE;
}

