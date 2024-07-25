// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

#include <Epetra_Vector.h>
#include <EpetraExt_CrsMatrixIn.h>
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_BlockMapIn.h"

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

int InitMatValues( const Epetra_CrsMatrix& newA, Epetra_CrsMatrix* A );
int InitMVValues( const Epetra_MultiVector& newb, Epetra_MultiVector* b );

int main(int argc, char** argv)
{
    // Initialize MPI environment
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif

    bool success = true;
    std::string pass = "End Result: TEST PASSED";
    std::string fail = "End Result: TEST FAILED";

    int nProcs = Comm.NumProc();
    int myPID = Comm.MyPID();
    bool verbose = (myPID == 0);

    if (verbose) {
        std::cout << "--- Running ShyLU-IQR test with " << nProcs <<
         " processes." << std::endl;
    }

    // Register Ifpack_ShyLU with the Ifpack factory.
#if defined(HAVE_IFPACK_DYNAMIC_FACTORY) && defined(HAVE_IFPACK_PARALLEL_SUBDOMAIN_SOLVERS)
    {
        int err = Ifpack_DynamicFactory::RegisterPreconditioner("ShyLU",
                        buildShyLU);
        if (! err) {
            if (verbose) {
                std::cout << "--- Ifpack_ShyLU was registered with"
                <<" Ifpack_DynamicFactory" << std::endl;
            }
        } else {
            if (verbose) {
                std::cout << "!!! Error registering preconditioner ShyLU with"
                          << " Ifpack_DynamicFactory. Exiting." << std::endl;
            }
            std::cout << fail << std::endl;
#ifdef HAVE_MPI
            MPI_Finalize();
#endif
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

    // Accept an alternate XML parameter file name, read from the command line
    {
        Teuchos::CommandLineProcessor cl(false, false);
        cl.setDocString ("");
        cl.setOption ("parameter-file",
                      &parameterFileName,
                      "The name of the XML test parameter file");
        Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn =
                                                     cl.parse (argc, argv);
        if ( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
          std::cout << fail << std::endl;
#ifdef HAVE_MPI
            MPI_Finalize();
#endif
            return 0;
        }
        if( parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL   ) {
          std::cout << fail << std::endl;
#ifdef HAVE_MPI
            MPI_Finalize();
#endif
            return -2;
        }
    }
    if (verbose) {
        std::cout << "--- Using parameter file: " << parameterFileName <<
         std::endl;
    }

    // Read the XML parameter file. The resulting Teuchos parameter list
    // contains sublists for
    // the partitioner, the preconditioner, linear solver etc.
    RCP<ParameterList> globalParams    = rcp(new ParameterList());
    globalParams = Teuchos::getParametersFromXmlFile(parameterFileName);

    // Read the matrix
    std::string matrixFileName
        = globalParams->get<std::string>("matrix file name", "wathenSmall.mtx");
    std::string rhsFileName = globalParams->get<std::string>("rhs_file", "");
    std::string mapFileName = globalParams->get<std::string>("map_file", "");

    int maxFiles = globalParams->get<int>("Maximum number of files to read in",
                             1);
    int startFile = globalParams->get<int>("Number of initial file", 1);
    int file_number = startFile;

    char file_name[200];
    if (mapFileName != "" && maxFiles > 1) {
        mapFileName += "%d.mm";
        sprintf( file_name, mapFileName.c_str(), file_number );
    }
    else {
        strcpy( file_name, mapFileName.c_str());
    }
    if (verbose) {
        std::cout << "--- Using map file: " << file_name << std::endl;
    }

    bool mapAvail = false;
    Epetra_BlockMap *vecMap2 = NULL;
    if (mapFileName != "")
    {
        mapAvail = true;
        int err = EpetraExt::MatrixMarketFileToBlockMap(file_name, Comm, vecMap2);
        if (err) {
            if (verbose) {
                std::cout << "!!! Matrix file could not be read in, info = "<<
                     err << std::endl;
            }
            std::cout << fail << std::endl;
#ifdef HAVE_MPI
            MPI_Finalize();
#endif
            return -3;
        }
    }
    Epetra_Map* vecMap = static_cast<Epetra_Map *>( vecMap2 );


    if (maxFiles > 1) {
        matrixFileName += "%d.mm";
        sprintf( file_name, matrixFileName.c_str(), file_number );
    }
    else {
        strcpy( file_name, matrixFileName.c_str());
    }

    if (verbose) {
        std::cout << "--- Using matrix file: " << file_name << std::endl;
    }

    Epetra_CrsMatrix *A;
    int err;
    if (mapAvail) {
        err = EpetraExt::MatrixMarketFileToCrsMatrix(file_name, *vecMap,  A);
    }
    else {
        err = EpetraExt::MatrixMarketFileToCrsMatrix(file_name, Comm, A);
    }
    if (err) {
        if (verbose) {
            std::cout << "!!! Matrix file could not be read in, info = "<<
             err << std::endl;
        }
        std::cout << fail << std::endl;
#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        return -3;
    }


    if (rhsFileName != "" && maxFiles > 1)
    {
        rhsFileName += "%d.mm";
        sprintf( file_name, rhsFileName.c_str(), file_number );
    }
    else
    {
        strcpy( file_name, rhsFileName.c_str());
    }
    if (verbose) {
        std::cout << "--- Using rhs file: " << file_name << std::endl;
    }

    // Partition the matrix
    ParameterList isorropiaParams = globalParams->sublist
                                    ("Isorropia parameters");
    RCP<Isorropia::Epetra::Partitioner> partitioner
            = rcp(new Isorropia::Epetra::Partitioner(A, isorropiaParams, false),
                    true);
    partitioner->partition();
    RCP<Isorropia::Epetra::Redistributor> rd
            = rcp(new Isorropia::Epetra::Redistributor(partitioner));

    int numRows = A->NumGlobalRows();
    RCP<Epetra_Map> vectorMap = rcp(new Epetra_Map(numRows, 0, Comm));
    Epetra_MultiVector *RHS;
    bool allOneRHS = false;

    if (rhsFileName != "")
    {
        if (mapAvail) {
        err = EpetraExt::MatrixMarketFileToMultiVector(file_name, *vecMap, RHS);
        }
        else {
        err = EpetraExt::MatrixMarketFileToMultiVector(file_name, *vectorMap, RHS);
        }
    }
    else
    {
        // Generate a RHS and LHS
        allOneRHS = true;
        if (mapAvail) {
            RHS = new Epetra_MultiVector(*vecMap, 1, false);
        }
        else {
            RHS = new Epetra_MultiVector(*vectorMap, 1, false);
        }
        RHS->Random();
    }

    Epetra_MultiVector *LHS;
    if (mapAvail) {
        LHS = new Epetra_MultiVector(*vecMap, 1, true);
    }
    else {
        LHS = new Epetra_MultiVector(*vectorMap, 1, true);
    }
    LHS->PutScalar(0.0);

    // Redistribute matrix and vectors
    RCP<Epetra_CrsMatrix> rcpA;
    RCP<Epetra_MultiVector> rcpRHS;
    RCP<Epetra_MultiVector> rcpLHS;
    Epetra_CrsMatrix *newA;
    Epetra_MultiVector *newLHS, *newRHS;

    if (mapAvail) {
        rcpA = rcp(A, true);
        rcpRHS = rcp(RHS, true);
        rcpLHS = rcp(LHS, true);
    }
    else {
        rd->redistribute(*A, newA);
        delete A;
        rcpA = rcp(newA, true);

        rd->redistribute(*RHS, newRHS);
        delete RHS;
        rcpRHS = rcp(newRHS, true);

        rd->redistribute(*LHS, newLHS);
        delete LHS;
        rcpLHS = rcp(newLHS, true);
    }

    if (! globalParams->isSublist("ML parameters")) {
        if (verbose) {
            std::cout << "!!! ML parameter list not found. Exiting."
                << std::endl;
        }
        std::cout << fail << std::endl;
#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        return -4;
    }
    ParameterList mlParameters = globalParams->sublist("ML parameters");

    if (! globalParams->isSublist("Belos parameters")) {
        if (verbose) {
            std::cout << "!!! Belos parameter list not found. Exiting." <<
                     std::endl;
        }
        std::cout << fail << std::endl;
#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        return -5;
    }
    ParameterList belosParams = globalParams->sublist("Belos parameters");
    int belosBlockSize = belosParams.get<int>("Block Size", 1);
    int belosMaxRestarts = belosParams.get<int>("Maximum Restarts", 0);
    int belosMaxIterations = belosParams.get<int>("Maximum Iterations",
                                 numRows);
    double belosTolerance = belosParams.get<double>("Convergence Tolerance",
                             1e-10);
    if (verbose) {
        std::cout << std::endl;
        std::cout << "--- Dimension of matrix: " << numRows << std::endl;
        std::cout << "--- Block size used by solver: " << belosBlockSize <<
                         std::endl;
        std::cout << "--- Number of restarts allowed: " << belosMaxRestarts
                             << std::endl;
        std::cout << "--- Max number of Gmres iterations per restart cycle: "
                  << belosMaxIterations << std::endl;
        std::cout << "--- Relative residual tolerance: " << belosTolerance
                     << std::endl;
    }

    Epetra_CrsMatrix *iterA = 0;
    Epetra_CrsMatrix *redistA = 0;
    Epetra_MultiVector *iterb1 = 0;
    RCP<ML_Epetra::MultiLevelPreconditioner> MLprec;
    Teuchos::Time setupPrecTimer("preconditioner setup timer", false);
    Teuchos::Time linearSolverTimer("linear solver timer");
    RCP< Belos::SolverManager<double, MV, OP> > solver;
    while(file_number < maxFiles+startFile)
    {
        /*if (file_number == startFile)
        {*/
            // Build ML preconditioner
            if (mapAvail)
            {
                setupPrecTimer.start();
                MLprec = rcp(new ML_Epetra::MultiLevelPreconditioner(*A,
                             mlParameters, false), true);
                MLprec->ComputePreconditioner();
                setupPrecTimer.stop();
            }
            else
            {
                setupPrecTimer.start();
                MLprec = rcp(new ML_Epetra::MultiLevelPreconditioner(*newA,
                             mlParameters, false), true);
                MLprec->ComputePreconditioner();
                setupPrecTimer.stop();
            }
        /*}
        else
        {
            setupPrecTimer.start();
            MLprec->ReComputePreconditioner();
            setupPrecTimer.stop();
        }*/

    // Build linear solver

    RCP<Belos::EpetraPrecOp> belosPrec = rcp(new Belos::EpetraPrecOp(MLprec),
                                         false);

    // Construct a preconditioned linear problem
    RCP<Belos::LinearProblem<double, MV, OP> > problem
    = rcp( new Belos::LinearProblem<double, MV, OP>( rcpA, rcpLHS, rcpRHS ) );
    problem->setRightPrec( belosPrec );

    if (! problem->setProblem()) {
        if (verbose) {
          std::cout << "!!! Belos::LinearProblem failed to set up correctly."
          << " Exiting." << std::endl;
        }
        std::cout << fail << std::endl;
#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        return -6;
    }

    // Create an iterative solver manager
    solver = rcp( new Belos::BlockGmresSolMgr<double, MV, OP>(problem,
                     rcp(&belosParams,false)));
    // Solve linear system
    linearSolverTimer.start();
    Belos::ReturnType ret = solver->solve();
    linearSolverTimer.stop();
    if (ret == Belos::Unconverged) {
        if (verbose) {
            std::cout << "!!! Linear solver did not converge to prescribed"
            <<" precision. Test failed." << std::endl;
        }
        std::cout << fail << std::endl;
#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        return -7;
    }
    // Print time measurements
    int numIters = solver->getNumIters();
    double timeSetupPrec = setupPrecTimer.totalElapsedTime();
    double timeLinearSolver = linearSolverTimer.totalElapsedTime();
    if (verbose) {
        std::cout << "--- Preconditioner setup time: " << timeSetupPrec <<
         std::endl;
        std::cout << "--- Number of iterations performed: " << numIters <<
         std::endl;
        std::cout << "--- Time to GMRES convergence: " << timeLinearSolver <<
         std::endl;
        std::cout << "--- Average time per GMRES iteration: "
                  << timeLinearSolver / numIters << std::endl;
        std::cout << "--- Total time to solution (prec + GMRES): "
                  << timeSetupPrec + timeLinearSolver << std::endl;
    }

    file_number++;
    if (file_number >= maxFiles+startFile)
    {
      if (redistA != NULL) { delete redistA; redistA = NULL; }
      if (iterb1 != NULL){delete iterb1; iterb1 = NULL;}
      break;
    }
    else
    {
        sprintf(file_name, matrixFileName.c_str(), file_number);

        if (redistA != NULL) { delete redistA; redistA = NULL; }
        // Load the new matrix
        if (mapAvail) {
            err = EpetraExt::MatrixMarketFileToCrsMatrix
                            (file_name, *vecMap, redistA);
        }
        else {
            err = EpetraExt::MatrixMarketFileToCrsMatrix(file_name, Comm,
                        iterA);
        }
        if (err != 0) {
          if (myPID == 0) std::cout << "Could not open file: "<< file_name << std::endl;
          success = false;
        }
        else
        {
            if (mapAvail)
            {
                InitMatValues(*redistA, A);
                if (redistA != NULL) { delete redistA; redistA = NULL; }
            }
            else {
                rd->redistribute(*iterA, redistA);
                if (iterA != NULL){delete iterA; iterA = NULL;}
                InitMatValues(*redistA, newA);
            }
        }

        // Load the new rhs
        if (!allOneRHS)
        {
            sprintf(file_name, rhsFileName.c_str(), file_number);

            if (iterb1 != NULL){delete iterb1; iterb1 = NULL;}
            if (mapAvail) {
                err = EpetraExt::MatrixMarketFileToMultiVector(file_name,
                                *vecMap, iterb1);
            }
            else {
                err = EpetraExt::MatrixMarketFileToMultiVector(file_name,
                                     *vectorMap, RHS);
            }
            if (err != 0) {
                if (myPID==0)
                    std::cout << "Could not open file: "<< file_name << std::endl;
                success = false;
            }
            else {
                if (mapAvail) {
                    InitMVValues( *iterb1, RHS );
                    if (iterb1 != NULL){ delete iterb1; iterb1 = NULL; }
                }
                else {
                    rd->redistribute(*RHS, iterb1);
                    if (RHS != NULL){ delete RHS; RHS = NULL; }
                    InitMVValues( *iterb1, newRHS );
                    // Should we delete iterb1
                }
            }
        }
    }
    }

    // Release the ML preconditioner, destroying MPI subcommunicators before the call to
    // MPI_Finalize() in order to avoid MPI errors (this is related to ParMETIS, I think).



    if(success)
      std::cout << pass << std::endl;
    else
      std::cout << fail << std::endl;



    MLprec.reset();
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
}


int InitMatValues( const Epetra_CrsMatrix& newA, Epetra_CrsMatrix* A )
{
  int numMyRows = newA.NumMyRows();
  int maxNum = newA.MaxNumEntries();
  int numIn;
  int *idx = 0;
  double *vals = 0;

  idx = new int[maxNum];
  vals = new double[maxNum];

  // For each row get the values and indices, and replace the values in A.
  for (int i=0; i<numMyRows; ++i) {

    // Get the values and indices from newA.
    EPETRA_CHK_ERR( newA.ExtractMyRowCopy(i, maxNum, numIn, vals, idx) );

    // Replace the values in A
    EPETRA_CHK_ERR( A->ReplaceMyValues(i, numIn, vals, idx) );

  }

  // Clean up.
  delete [] idx;
  delete [] vals;

  return 0;
}

int InitMVValues( const Epetra_MultiVector& newb, Epetra_MultiVector* b )
{
  int length = newb.MyLength();
  int numVecs = newb.NumVectors();
  const Epetra_Vector *tempnewvec;
  Epetra_Vector *tempvec = 0;

  for (int i=0; i<numVecs; ++i) {
    tempnewvec = newb(i);
    tempvec = (*b)(i);
    for (int j=0; j<length; ++j)
      (*tempvec)[j] = (*tempnewvec)[j];
  }

  return 0;
}
