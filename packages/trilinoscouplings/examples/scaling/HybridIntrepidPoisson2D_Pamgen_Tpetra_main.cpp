// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file  HybridIntrepidPoisson2D_Pamgen_Tpetra_main.cpp
    \brief Example: Discretize Poisson's equation with Dirichlet
           boundary conditions on a quadrilateral mesh using nodal
           (Hgrad) elements.  The system is assembled into Tpetra data
           structures, and optionally solved.

           This example uses the following Trilinos packages:
    \li    Pamgen to generate a Quadrilateral mesh.
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
            Hoemmen (mhoemme@sandia.gov) and a bunch of other people.

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

#include "Tpetra_Core.hpp"

#include "TrilinosCouplings_config.h"
#include "TrilinosCouplings_TpetraIntrepidHybridPoisson2DExample.hpp"
#include "TrilinosCouplings_IntrepidPoissonExampleHelpers.hpp"

#ifdef HAVE_TRILINOSCOUPLINGS_MUELU
#  include "MueLu_CreateTpetraPreconditioner.hpp"
#endif // HAVE_TRILINOSCOUPLINGS_MUELU

#include <MatrixMarket_Tpetra.hpp>

int
main (int argc, char *argv[])
{
  using namespace TrilinosCouplings; // Yes, this means I'm lazy.

  using TpetraIntrepidPoissonExample::exactResidualNorm;
  using TpetraIntrepidPoissonExample::makeMatrixAndRightHandSide;
  using TpetraIntrepidPoissonExample::solveWithBelos;
  using TpetraIntrepidPoissonExample::solveWithBelosGPU;
  using IntrepidPoissonExample::makeMeshInput;
  using IntrepidPoissonExample::parseCommandLineArguments;
  using IntrepidPoissonExample::setCommandLineArgumentDefaults;
  using IntrepidPoissonExample::setMaterialTensorOffDiagonalValue;
  using IntrepidPoissonExample::setUpCommandLineArguments;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::getFancyOStream;
  using Teuchos::FancyOStream;
  using std::endl;
  // Pull in typedefs from the example's namespace.
  typedef TpetraIntrepidPoissonExample::ST ST;
#ifdef HAVE_TRILINOSCOUPLINGS_MUELU
  typedef TpetraIntrepidPoissonExample::LO LO;
  typedef TpetraIntrepidPoissonExample::GO GO;
#else
  typedef TpetraIntrepidPoissonExample::sparse_matrix_type::local_ordinal_type LO;
  typedef TpetraIntrepidPoissonExample::sparse_matrix_type::global_ordinal_type GO;
#endif // HAVE_TRILINOSCOUPLINGS_MUELU
  typedef TpetraIntrepidPoissonExample::Node Node;
  typedef Teuchos::ScalarTraits<ST> STS;
  typedef STS::magnitudeType MT;
  typedef Teuchos::ScalarTraits<MT> STM;
  typedef TpetraIntrepidPoissonExample::sparse_matrix_type sparse_matrix_type;
  typedef TpetraIntrepidPoissonExample::vector_type vector_type;
  typedef TpetraIntrepidPoissonExample::operator_type operator_type;

  bool success = true;
  try {
    Teuchos::oblackholestream blackHole;
    Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
    const int myRank = mpiSession.getRank ();
    //const int numProcs = mpiSession.getNProc ();

    // Get the default communicator and Kokkos Node instance
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();

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

    // Additional command-line arguments for GPU experimentation.
    bool gpu = false;
    cmdp.setOption ("gpu", "no-gpu", &gpu,
                    "Run example using GPU node (if supported)");
    int ranks_per_node = 1;
    cmdp.setOption ("ranks_per_node", &ranks_per_node,
                    "Number of MPI ranks per node");
    int gpu_ranks_per_node = 1;
    cmdp.setOption ("gpu_ranks_per_node", &gpu_ranks_per_node,
                    "Number of MPI ranks per node for GPUs");
    int device_offset = 0;
    cmdp.setOption ("device_offset", &device_offset,
                    "Offset for attaching MPI ranks to CUDA devices");

    // Additional command-line arguments for dumping the generated
    // matrix or its row Map to output files.
    //
    // FIXME (mfh 09 Apr 2014) Need to port these command-line
    // arguments to the Epetra version.

    // If matrixFilename is nonempty, dump the matrix to that file
    // in MatrixMarket format.
    std::string matrixFilename;
    cmdp.setOption ("matrixFilename", &matrixFilename, "If nonempty, dump the "
                    "generated matrix to that file in MatrixMarket format.");

    // If rowMapFilename is nonempty, dump the matrix's row Map to
    // that file in MatrixMarket format.
    std::string rowMapFilename;
    cmdp.setOption ("rowMapFilename", &rowMapFilename, "If nonempty, dump the "
                    "generated matrix's row Map to that file in a format that "
                    "Tpetra::MatrixMarket::Reader can read.");
    // Option to exit after building A and b (and dumping stuff to
    // files, if requested).
    bool exitAfterAssembly = false;
    cmdp.setOption ("exitAfterAssembly", "dontExitAfterAssembly",
                    &exitAfterAssembly, "If true, exit after building the "
                    "sparse matrix and dense right-hand side vector.  If either"
                    " --matrixFilename or --rowMapFilename are nonempty strings"
                    ", dump the matrix resp. row Map to their respective files "
                    "before exiting.");

    // Option to exit after building the preconditioner
    bool exitAfterPrecond = false;
    cmdp.setOption ("exitAfterPrecond", "dontExitAfterPrecond",
                    &exitAfterPrecond, "If true, exit after building the "
                    "preconditioner.");


     // If matrixFilename is nonempty, dump the matrix to that file
    // in MatrixMarket format.
    int numMueluRebuilds=0;
    cmdp.setOption ("rebuild", &numMueluRebuilds, "Number of times to rebuild the MueLu hierarchy.");

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

      RCP<sparse_matrix_type> A;
      RCP<vector_type> B, X_exact, X;
      Teuchos::Array<Teuchos::Array<ST> > coordsArray(3);
      Teuchos::Array<LO> lNodesPerDim(3);
      {
        TEUCHOS_FUNC_TIME_MONITOR_DIFF("Total Assembly", total_assembly);
        makeMatrixAndRightHandSide (A, B, X_exact, X, coordsArray, lNodesPerDim,
                                    comm, meshInput, out, err, verbose, debug);
      }

      // Optionally dump the matrix and/or its row Map to files.
      {
        typedef Tpetra::MatrixMarket::Writer<sparse_matrix_type> writer_type;
        if (matrixFilename != "") {
          writer_type::writeSparseFile (matrixFilename, A);
        }
        if (rowMapFilename != "") {
          writer_type::writeMapFile (rowMapFilename, * (A->getRowMap ()));
        }
      }

      if (exitAfterAssembly) {
        // Users might still be interested in assembly time.
        Teuchos::TimeMonitor::report (comm.ptr (), std::cout);
        return EXIT_SUCCESS;
      }

      const std::vector<MT> norms = exactResidualNorm (A, B, X_exact);
      // X_exact is the exact solution of the PDE, projected onto the
      // discrete mesh.  It may not necessarily equal the exact solution
      // of the linear system.
      *out << "||B - A*X_exact||_2 = " << norms[0] << endl
           << "||B||_2 = " << norms[1] << endl
           << "||A||_F = " << norms[2] << endl;

      // Convert coordsArray to multivector
      int dim = 2;
      Teuchos::Array<Teuchos::ArrayView<const double>> coordsArrayView(dim);
      for(int i=0; i<dim; ++i) { coordsArrayView[i] = coordsArray[i](); }
      RCP< Tpetra::MultiVector<double, LO, GO, Node> > coordinates = rcp (new Tpetra::MultiVector<double,LO,GO,Node>(A->getRowMap(),coordsArrayView,dim));

      // Set region type
      std::string regionType;
      if(myRank % 2 == 0){
        regionType = "uncoupled";
      } else {
        regionType = "structured";
      }

      // Setup preconditioner
      std::string prec_type = inputList.get ("Preconditioner", "None");
      RCP<operator_type> M;
      RCP<operator_type> opA(A);
      {
        TEUCHOS_FUNC_TIME_MONITOR_DIFF("Total Preconditioner Setup", total_prec);

        if (prec_type == "MueLu") {
#ifdef HAVE_TRILINOSCOUPLINGS_MUELU
	  for(int i=0; i<numMueluRebuilds+1; i++) {
	    if (inputList.isSublist("MueLu")) {
	      ParameterList mueluParams = inputList.sublist("MueLu");
              const std::string userName = "user data";
              Teuchos::ParameterList& userParamList = mueluParams.sublist(userName);
              userParamList.set<Teuchos::Array<LO> >("Array<LO> lNodesPerDim", lNodesPerDim);
              userParamList.set< RCP<Tpetra::MultiVector<double,LO,GO,Node>> >("Coordinates", coordinates);
              userParamList.set< std::string >("string aggregationRegionType", regionType);
	      M = MueLu::CreateTpetraPreconditioner<ST,LO,GO,Node>(opA, mueluParams);
	    } else {
	      M = MueLu::CreateTpetraPreconditioner<ST,LO,GO,Node>(opA);
	    }
	  }
#else // NOT HAVE_TRILINOSCOUPLINGS_MUELU
          TEUCHOS_TEST_FOR_EXCEPTION(
            prec_type == "MueLu", std::runtime_error, "Tpetra scaling example: "
            "In order to precondition with MueLu, you must have built Trilinos "
            "with the MueLu package enabled.");
#endif // HAVE_TRILINOSCOUPLINGS_MUELU
        }
      } // setup preconditioner


      // exit if we only wanted the preconditioner
      if(exitAfterPrecond) {
        Teuchos::TimeMonitor::report (comm.ptr (), std::cout);
        return EXIT_SUCCESS;
      }


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
      if (gpu) {
        TEUCHOS_FUNC_TIME_MONITOR_DIFF("Total GPU Solve", total_solve);
        solveWithBelosGPU (converged, numItersPerformed, tol, maxNumIters,
                           num_steps, ranks_per_node, gpu_ranks_per_node,
                           device_offset, prec_type, X, A, B, Teuchos::null, M);
      }
      else {
        TEUCHOS_FUNC_TIME_MONITOR_DIFF("Total Solve", total_solve);
        solveWithBelos (converged, numItersPerformed, solverName, tol,
                        maxNumIters, num_steps, X, A, B, Teuchos::null, M);
      }

      // Compute ||X-X_exact||_2
      const MT norm_x = X_exact->norm2 ();
      X_exact->update (-1.0, *X, 1.0);
      const MT norm_error = X_exact->norm2 ();
      *out << endl
           << "||X - X_exact||_2 / ||X_exact||_2 = " << norm_error / norm_x
           << endl;
    } // total time block

    // Summarize timings
    Teuchos::TimeMonitor::report (comm.ptr (), std::cout);
  } // try
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  if (success) {
    return EXIT_SUCCESS;
  } else {
    return EXIT_FAILURE;
  }
}
