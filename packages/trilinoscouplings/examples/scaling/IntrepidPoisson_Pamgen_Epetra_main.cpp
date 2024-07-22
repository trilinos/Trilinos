// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file  IntrepidPoisson_Pamgen_Epetra_main.cpp
    \brief Example: Discretize Poisson's equation with Dirichlet
           boundary conditions on a hexahedral mesh using nodal
           (Hgrad) elements.  The system is assembled into Epetra data
           structures, and optionally solved.

           This example uses the following Trilinos packages:
    \li    Pamgen to generate a Hexahedral mesh.
    \li    Sacado to form the source term from user-specified manufactured solution.
    \li    Intrepid to build the discretization matrix and right-hand side.
    \li    Epetra to handle the global sparse matrix and dense vector.

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
            Hoemmen (mhoemme@sandia.gov) and a bunch of other people,
            and back-ported to Epetra for a fair comparison with
            Tpetra.

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
#  include "Teuchos_DefaultMpiComm.hpp"
#else
#  include "Epetra_SerialComm.h"
#  include "Teuchos_DefaultSerialComm.hpp"
#endif // EPETRA_MPI

#include "TrilinosCouplings_config.h"
#include "TrilinosCouplings_EpetraIntrepidPoissonExample.hpp"
#include "TrilinosCouplings_IntrepidPoissonExampleHelpers.hpp"

#ifdef HAVE_TRILINOSCOUPLINGS_ML
#  include "ml_include.h"
#  include "ml_MultiLevelPreconditioner.h"
#endif // HAVE_TRILINOSCOUPLINGS_ML

#ifdef HAVE_TRILINOSCOUPLINGS_MUELU
#  include "MueLu.hpp"
#  include "MueLu_ParameterListInterpreter.hpp"
#  include "MueLu_EpetraOperator.hpp"
#endif // HAVE_TRILINOSCOUPLINGS_MUELU

int
main (int argc, char *argv[])
{
  using namespace TrilinosCouplings; // Yes, this means I'm lazy.

  using EpetraIntrepidPoissonExample::exactResidualNorm;
  using EpetraIntrepidPoissonExample::makeMatrixAndRightHandSide;
  using EpetraIntrepidPoissonExample::solveWithBelos;
  using IntrepidPoissonExample::makeMeshInput;
  using IntrepidPoissonExample::parseCommandLineArguments;
  using IntrepidPoissonExample::setCommandLineArgumentDefaults;
  using IntrepidPoissonExample::setMaterialTensorOffDiagonalValue;
  using IntrepidPoissonExample::setUpCommandLineArguments;
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

    // Get the default communicator
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
    int maxNumItersFromCmdLine = -1; // -1 means "read from XML file"
    double tolFromCmdLine = -1.0; // -1 means "read from XML file"
    std::string solverName = "GMRES";
    ST materialTensorOffDiagonalValue = 0.0;

    // Set default values of command-line arguments.
    setCommandLineArgumentDefaults (nx, ny, nz, xmlInputParamsFile,
                                    solverName, verbose, debug);
    // Parse and validate command-line arguments.
    Teuchos::CommandLineProcessor cmdp (false, true);
    setUpCommandLineArguments (cmdp, nx, ny, nz, xmlInputParamsFile,
                               solverName, tolFromCmdLine,
                               maxNumItersFromCmdLine,
                               verbose, debug);
    cmdp.setOption ("materialTensorOffDiagonalValue",
                    &materialTensorOffDiagonalValue, "Off-diagonal value in "
                    "the material tensor.  This controls the iteration count.  "
                    "Be careful with this if you use CG, since you can easily "
                    "make the matrix indefinite.");

    parseCommandLineArguments (cmdp, printedHelp, argc, argv, nx, ny, nz,
                               xmlInputParamsFile, solverName, verbose, debug);
    if (printedHelp) {
      // The user specified --help at the command line to print help
      // with command-line arguments.  We printed help already, so quit
      // with a happy return code.
      return EXIT_SUCCESS;
    }

    setMaterialTensorOffDiagonalValue (materialTensorOffDiagonalValue);

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
      std::string prec_type = inputList.get ("Preconditioner", "None");
      RCP<operator_type> M;
      {
        TEUCHOS_FUNC_TIME_MONITOR_DIFF("Total Preconditioner Setup", total_prec);

        if (prec_type == "ML") {
#ifdef HAVE_TRILINOSCOUPLINGS_ML
          ParameterList mlParams;
          if (inputList.isSublist("ML"))
            mlParams = inputList.sublist("ML");
          else {
            ML_Epetra::SetDefaults("SA", mlParams);
            mlParams.set("ML output", 0);
          }
          M = rcp(new ML_Epetra::MultiLevelPreconditioner(*A, mlParams));
#else // NOT HAVE_TRILINOSCOUPLINGS_ML
          TEUCHOS_TEST_FOR_EXCEPTION(
            prec_type == "ML", std::runtime_error, "Epetra scaling example: "
            "In order to precondition with ML, you must have built Trilinos "
            "with the ML package enabled.");
#endif // HAVE_TRILINOSCOUPLINGS_ML
        }
        else if (prec_type == "MueLu") {
#ifdef HAVE_TRILINOSCOUPLINGS_MUELU
          // Turns a Epetra_CrsMatrix into a MueLu::Matrix
          RCP<Xpetra::CrsMatrix<ST, int, int, Xpetra::EpetraNode> > mueluA_ =
            rcp(new Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>(A));
          RCP<Xpetra::Matrix <ST, int, int, Xpetra::EpetraNode> > mueluA  =
            rcp(new Xpetra::CrsMatrixWrap<ST, int, int, Xpetra::EpetraNode>(mueluA_));

          // Multigrid Hierarchy
          ParameterList mueluParams;
          if (inputList.isSublist("MueLu"))
            mueluParams = inputList.sublist("MueLu");
          MueLu::ParameterListInterpreter<ST, int, int, Xpetra::EpetraNode> mueLuFactory(mueluParams);
          RCP<MueLu::Hierarchy<ST, int, int, Xpetra::EpetraNode> > H =
            mueLuFactory.CreateHierarchy();
          H->setVerbLevel(Teuchos::VERB_HIGH);
          H->GetLevel(0)->Set("A", mueluA);

          // Multigrid setup phase
          H->Setup();

          // Wrap MueLu Hierarchy as a Tpetra::Operator
          M = rcp(new MueLu::EpetraOperator(H));
#else // NOT HAVE_TRILINOSCOUPLINGS_MUELU
          TEUCHOS_TEST_FOR_EXCEPTION(
            prec_type == "MueLu", std::runtime_error, "Epetra scaling example: "
            "In order to precondition with MueLu, you must have built Trilinos "
            "with the MueLu package enabled.");
#endif // HAVE_TRILINOSCOUPLINGS_MUELU
        }
      } // setup preconditioner

      // Get the convergence tolerance for each linear solve.
      // If the user provided a nonnegative value at the command
      // line, it overrides any value in the input ParameterList.
      MT tol = STM::squareroot (STM::eps ()); // default value
      if (tolFromCmdLine < STM::zero ()) {
        tol = inputList.get ("Convergence Tolerance", tol);
      } else {
        tol = tolFromCmdLine;
      }

      // Get the maximum number of iterations for each linear solve.
      // If the user provided a value other than -1 at the command
      // line, it overrides any value in the input ParameterList.
      int maxNumIters = 200; // default value
      if (maxNumItersFromCmdLine == -1) {
        maxNumIters = inputList.get ("Maximum Iterations", maxNumIters);
      } else {
        maxNumIters = maxNumItersFromCmdLine;
      }

      // Get the number of "time steps."  We imitate a time-dependent
      // PDE by doing this many linear solves.
      const int num_steps = inputList.get ("Number of Time Steps", 1);

      // Do the linear solve(s).
      bool converged = false;
      int numItersPerformed = 0;
      {
        TEUCHOS_FUNC_TIME_MONITOR_DIFF("Total Solve", total_solve);
        solveWithBelos (converged, numItersPerformed, solverName, tol,
                        maxNumIters, num_steps, X, A, B, Teuchos::null, M);
      }

      // Compute ||X-X_exact||_2
      MT norm_x, norm_error;
      X_exact->Norm2 (&norm_x);
      X_exact->Update (-1.0, *X, 1.0);
      X_exact->Norm2 (&norm_error);
      *out << endl
           << "||X - X_exact||_2 / ||X_exact||_2 = " << norm_error / norm_x
           << endl;
    } // total time block

    // Summarize timings
    {
      // Make a Teuchos::Comm corresponding to the Epetra_Comm.
      RCP<const Teuchos::Comm<int> > teuchosComm;
#ifdef EPETRA_MPI
      teuchosComm = rcp (new Teuchos::MpiComm<int> (MPI_COMM_WORLD));
#else
      teuchosComm = rcp (new Teuchos::SerialComm<int> ());
#endif // EPETRA_MPI

      // Use the Teuchos::Comm to report timings.
      Teuchos::TimeMonitor::report (teuchosComm.ptr (), std::cout);
    }
  } // try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  if (success) {
    return EXIT_SUCCESS;
  } else {
    return EXIT_FAILURE;
  }
}

