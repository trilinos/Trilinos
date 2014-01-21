/** \file shylu_iqr_driver.cpp

    \brief Factors and solves a sparse matrix using ShyLU and IQR, if enabled.

    \author Radu Popescu <i.radu.popescu@gmail.com>

	* The example needs the following Ifpack variables enabled when compiling Trilinos:
	  -DIfpack_ENABLE_DYNAMIC_FACTORY:BOOL=ON (allows the use of the new factory which can register
	  new preconditioners within the client application)
	  -DIfpack_ENABLE_PARALLEL_SUBDOMAIN_SOLVERS:BOOL=ON (allows using MPI parallel subdomain
	  solvers for Ifpack AAS)
	* The example doesn't fail if these variables are disabled, but it doesn't do anything
	* interesting, it just solves the problem with Ifpack AAS with serial Amesos subdomain solvers.

	* Highlights:
	  * There is a builder function for Ifpack_ShyLU declared at the beginning of the file.
	  * Ifpack_ShyLU is registered with the Ifpack_DynamicFactory class at the beginning of the
		main() function, using the builder function that was declared earlier.
	  * An alternate XML file is provided in case the CMake variables mentioned earlier are not
	    enabled, to prevent the test from failing.
	  * A global Teuchos parameter list is read from an XML file, with sublists for the matrix
	    partitioner, the preconditioner and the linear solver
	  * The ML list contains a sublist for Ifpack (as usual), which in turn contains a sublist for
	  	all ShyLU parameters (this behaviour is enabled with the Ifpack_DynamicFactory switch, to
	  	avoid parameter list polution)
	  * Through Ifpack_DynamicFactory, ML is be able to build the Ifpack_AdditiveSchwarz<Ifpack_ShyLU>
	    preconditioner.

	* Key parameters:
	  * Isorropia parameters:
	    * we set the partitioning method to HIER_GRAPH, to ensure that the parts which make up each
	      AAS subdomain are connected
	    * we set TOPOLOGY to the number of processor per AAS subdomain
	  * ML parameters:
		* smoother: ifpack type - string - should be set to ShyLU
		* smoother: ifpack overlap - int - 0 - using multiple processors per AAS subdomain forces
		  this. Overlap is not supported.
	  * Ifpack parameters:
		* subdomain: number-of-processors - int - number of processors per AAS subdomain (must be a
		  divisor of the total number of MPI processes used).
	  * ShyLU parameters:
		* Schur Approximation Method - usually A22AndBlockDiagonals; can be set to IQR or G. IQR
		  means that we use IQR to solve the Schur complement system inexactly (Krylov subspace
		  reuse method), G means that we just approximate the inverse of the Schur complement with
		  an AAS of the G subblock of the subdomain matrix (which is faster, but leads to a looser
		  coupling at the subdomain level, for some problems it leads to fewer outer GMRES iterations
		  than IQR).
		* IQR Initial Prec Type - string - default: Amesos (this is the actual string given to an
		  Ifpack factory; it means AAS with serial Amesos on subdomains) - the preconditioner used
		  for the GMRES solver within IQR. This is also used for the approximation of the Schur
		  complement inverse when Schur Approximation Method is set to G.
		* IQR Initial Prec Amesos Type - string - default: Amesos_Klu - which Amesos solver to use
		  for the IQR preconditioner
	  * Amesos_Klu is given as a default in multiple places of the XML file. Should be substituted
	    with faster alternatives: Amesos_Pardiso, Amesos_Umfpack etc.

	\remark Usage:
    \code mpirun -n 2np ShyLU_iqr_driver.exe

*/

#include <string>

#include <Epetra_ConfigDefs.h>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_CrsMatrix.h>

#include <EpetraExt_CrsMatrixIn.h>

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosBlockGmresSolMgr.hpp>

#include <Ifpack_ConfigDefs.h>
#if defined(HAVE_IFPACK_DYNAMIC_FACTORY) && defined(HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS)
#include <Ifpack_DynamicFactory.h>
#include <Ifpack_AdditiveSchwarz.h>
#include <Ifpack_ShyLU.h>
#else
#include <Ifpack.h>
#endif

#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>

#include <ml_MultiLevelPreconditioner.h>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_Time.hpp>

namespace
{
typedef Epetra_MultiVector  MV;
typedef Epetra_Operator     OP;

#if defined(HAVE_IFPACK_DYNAMIC_FACTORY) && defined(HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS)
// Declare Ifpack_ShyLU builder function for later registration
Ifpack_Preconditioner* buildShyLU(Epetra_RowMatrix* Matrix,
								  int Overlap,
								  bool /*serial*/,
								  bool /*overrideSerialDefault*/)
{
	return new Ifpack_AdditiveSchwarz<Ifpack_ShyLU>(Matrix, Overlap);
}
#endif
}

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;

int main(int argc, char** argv)
{
	// Initialize MPI environment
#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif

    int nProcs = Comm.NumProc();
    int myPID = Comm.MyPID();
    bool verbose = (myPID == 0);

    if (verbose) {
    	std::cout << "--- Running ShyLU-IQR test with " << nProcs << " processes." << std::endl;
    }

    // Register Ifpack_ShyLU with the Ifpack factory.
#if defined(HAVE_IFPACK_DYNAMIC_FACTORY) && defined(HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS)
    {
    	int err = Ifpack_DynamicFactory::RegisterPreconditioner("ShyLU", buildShyLU);
    	if (! err) {
    		if (verbose) {
    			std::cout << "--- Ifpack_ShyLU was registered with Ifpack_DynamicFactory" << std::endl;
    		}
    	} else {
    		if (verbose) {
    			std::cout << "!!! Error registering preconditioner ShyLU with"
    					  << " Ifpack_DynamicFactory. Exiting."
    					  << std::endl;
    		}
    		MPI_Finalize();
    		return -1;
    	}
    }
#endif

	// Read the file containing the test parameters
#if defined(HAVE_IFPACK_DYNAMIC_FACTORY) && defined(HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS)
    std::string parameterFileName = "shylu_iqr_parameters.xml";
#else
    std::string parameterFileName = "shylu_no_iqr_parameters.xml";
#endif

    // Accept an alternate XML parameter file name, read in from the command line
    {
		Teuchos::CommandLineProcessor cl(false, false);
		cl.setDocString ("");
		cl.setOption ("parameter-file",
					  &parameterFileName,
					  "The name of the XML test parameter file");
		Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn = cl.parse (argc, argv);
		if ( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED ) {
			MPI_Finalize();
			return 0;
		}
		if( parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL   ) {
			MPI_Finalize();
			return -2;
		}
    }
    if (verbose) {
    	std::cout << "--- Using parameter file: " << parameterFileName << std::endl;
    }

    // Read the XML parameter file. The resulting Teuchos parameter list contains sublists for
    // the partitioner, the preconditioner, linear solver etc.
    RCP<ParameterList> globalParams	= rcp(new ParameterList());
    globalParams = Teuchos::getParametersFromXmlFile(parameterFileName);

	// Read the matrix
    std::string matrixFileName
    		= globalParams->get<std::string>("matrix file name", "wathenSmall.mtx");

    if (verbose) {
    	std::cout << "--- Using matrix file: " << matrixFileName << std::endl;
    }

    Epetra_CrsMatrix *A;
    int err = EpetraExt::MatrixMarketFileToCrsMatrix(matrixFileName.c_str(), Comm, A);
    if (err) {
    	if (verbose) {
    		std::cout << "!!! Matrix file could not be read in, info = "<< err << std::endl;
    	}
    	MPI_Finalize();
    	return -3;
    }

    // Partition the matrix
    ParameterList isorropiaParams = globalParams->sublist("Isorropia parameters");
    RCP<Isorropia::Epetra::Partitioner> partitioner
    		= rcp(new Isorropia::Epetra::Partitioner(A, isorropiaParams, false), true);
    partitioner->partition();
    RCP<Isorropia::Epetra::Redistributor> rd
    		= rcp(new Isorropia::Epetra::Redistributor(partitioner));

    // Generate a RHS and LHS
    int numRows = A->NumGlobalRows();
    RCP<Epetra_Map> vectorMap = rcp(new Epetra_Map(numRows, 0, Comm));
    Epetra_MultiVector *RHS = new Epetra_MultiVector(*vectorMap, 1, false);
    Epetra_MultiVector *LHS = new Epetra_MultiVector(*vectorMap, 1, true);
    RHS->Random();

    // Redistribute matrix and vectors
    RCP<Epetra_CrsMatrix> rcpA;
    RCP<Epetra_MultiVector> rcpRHS;
    RCP<Epetra_MultiVector> rcpLHS;

    {
		Epetra_CrsMatrix *newA;
		Epetra_MultiVector *newLHS, *newRHS;
		rd->redistribute(*A, newA);
		delete A;
		A = newA;
		rcpA = rcp(A, true);

		rd->redistribute(*RHS, newRHS);
		delete RHS;
		RHS = newRHS;
		rcpRHS = rcp(RHS, true);

		rd->redistribute(*LHS, newLHS);
		delete LHS;
		LHS = newLHS;
		rcpLHS = rcp(LHS, true);
    }

	// Build ML preconditioner
    if (! globalParams->isSublist("ML parameters")) {
    	if (verbose) {
    		std::cout << "!!! ML parameter list not found. Exiting." << std::endl;
    	}
    	MPI_Finalize();
    	return -4;
    }
    ParameterList mlParameters = globalParams->sublist("ML parameters");

    Teuchos::Time setupPrecTimer("preconditioner setup timer", false);
    setupPrecTimer.start();
    RCP<ML_Epetra::MultiLevelPreconditioner> MLprec
    		= rcp(new ML_Epetra::MultiLevelPreconditioner(*A, mlParameters, true), true);
    setupPrecTimer.stop();

	// Build linear solver
    if (! globalParams->isSublist("Belos parameters")) {
    	if (verbose) {
    		std::cout << "!!! Belos parameter list not found. Exiting." << std::endl;
    	}
    	MPI_Finalize();
    	return -5;
    }
    ParameterList belosParams = globalParams->sublist("Belos parameters");
    int belosBlockSize = belosParams.get<int>("Block Size", 1);
    int belosMaxRestarts = belosParams.get<int>("Maximum Restarts", 0);
    int belosMaxIterations = belosParams.get<int>("Maximum Iterations", numRows);
    double belosTolerance = belosParams.get<double>("Convergence Tolerance", 1e-10);

    RCP<Belos::EpetraPrecOp> belosPrec = rcp(new Belos::EpetraPrecOp(MLprec), false);

    // Construct a preconditioned linear problem
    RCP<Belos::LinearProblem<double, MV, OP> > problem
    		= rcp( new Belos::LinearProblem<double, MV, OP>( rcpA, rcpLHS, rcpRHS ) );
    problem->setRightPrec( belosPrec );

    if (! problem->setProblem()) {
		if (verbose) {
		  std::cout << "!!! Belos::LinearProblem failed to set up correctly. Exiting." << std::endl;
		}
		MPI_Finalize();
		return -6;
    }

    // Create an iterative solver manager
    RCP< Belos::SolverManager<double, MV, OP> > solver
    = rcp( new Belos::BlockGmresSolMgr<double, MV, OP>(problem, rcp(&belosParams,false)));

	// Solve linear system
	if (verbose) {
		std::cout << std::endl;
		std::cout << "--- Dimension of matrix: " << numRows << std::endl;
		std::cout << "--- Block size used by solver: " << belosBlockSize << std::endl;
		std::cout << "--- Number of restarts allowed: " << belosMaxRestarts << std::endl;
		std::cout << "--- Max number of Gmres iterations per restart cycle: "
				  << belosMaxIterations << std::endl;
		std::cout << "--- Relative residual tolerance: " << belosTolerance << std::endl;
	}

    Teuchos::Time linearSolverTimer("linear solver timer");
    linearSolverTimer.start();
	Belos::ReturnType ret = solver->solve();
	linearSolverTimer.stop();
	if (ret == Belos::Unconverged) {
		if (verbose) {
			std::cout << "!!! Linear solver did not converge to prescribed precision. Test failed." << std::endl;
		}
		MPI_Finalize();
		return -7;
	}

	// Print time measurements
	int numIters = solver->getNumIters();
	double timeSetupPrec = setupPrecTimer.totalElapsedTime();
	double timeLinearSolver = linearSolverTimer.totalElapsedTime();
	if (verbose) {
		std::cout << "--- Preconditioner setup time: " << timeSetupPrec << std::endl;
		std::cout << "--- Number of iterations performed: " << numIters << std::endl;
		std::cout << "--- Time to GMRES convergence: " << timeLinearSolver << std::endl;
		std::cout << "--- Average time per GMRES iteration: "
				  << timeLinearSolver / numIters << std::endl;
		std::cout << "--- Total time to solution (prec + GMRES): "
				  << timeSetupPrec + timeLinearSolver << std::endl;
	}

	// Release the ML preconditioner, destroying MPI subcommunicators before the call to
	// MPI_Finalize() in order to avoid MPI errors (this is related to ParMETIS, I think).
	MLprec.reset();

    MPI_Finalize();
	return 0;
}
