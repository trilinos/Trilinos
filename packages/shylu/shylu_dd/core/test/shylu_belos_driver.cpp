// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file shylu_driver.cpp

    \brief Factors and solves a sparse matrix using LU factorization.

    \author Siva Rajamanickam

    \remark Usage:
    \code mpirun -n np shylu_driver.exe

*/

#include <assert.h>
#include <iostream>
#include <sstream>

#include "Isorropia_config.h" // Just for HAVE_MPI

// Epetra includes
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#endif // HAVE_MPI
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LinearProblem.h"

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_MultiVectorIn.h"

// Amesos includes
#include "Amesos.h"
#include "Amesos_BaseSolver.h"

// AztecOO includes
#include "AztecOO.h"

// ML includes
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

#include "Ifpack_ConfigDefs.h"
#include "shylu.h"
#include "shylu_util.h"
#include "Ifpack_ShyLU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_ILU.h"

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

#include "EpetraExt_readEpetraLinearSystem.h"

//TODO Combine Belos driver and AztecOO driver into one driver.

int InitMatValues( const Epetra_CrsMatrix& newA, Epetra_CrsMatrix* A );
int InitMVValues( const Epetra_MultiVector& newb, Epetra_MultiVector* b );

using namespace std;

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, 0);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif
    typedef double                            ST;
    typedef Teuchos::ScalarTraits<ST>        SCT;
    typedef SCT::magnitudeType                MT;
    typedef Epetra_MultiVector                MV;
    typedef Epetra_Operator                   OP;
    typedef Belos::MultiVecTraits<ST,MV>     MVT;
    typedef Belos::OperatorTraits<ST,MV,OP>  OPT;
    using Teuchos::RCP;
    using Teuchos::rcp;


    bool success = true;
    string pass = "End Result: TEST PASSED";
    string fail = "End Result: TEST FAILED";



    bool verbose = false, proc_verbose = true;
    bool leftprec = false;      // left preconditioning or right.
    int frequency = -1;        // frequency of status test output.
    int blocksize = 1;         // blocksize
    int numrhs = 1;            // number of right-hand sides to solve for
    int maxrestarts = 15;      // maximum number of restarts allowed
    int maxsubspace = 25;      // maximum number of blocks the solver can use
                               // for the subspace
    char file_name[100];

    int nProcs, myPID ;
    Teuchos::RCP <Teuchos::ParameterList> pLUList ;        // ParaLU parameters
    Teuchos::ParameterList isoList ;        // Isorropia parameters
    Teuchos::ParameterList shyLUList ;    // ShyLU parameters
    string ipFileName = "ShyLU.xml";       // TODO : Accept as i/p

#ifdef HAVE_MPI
    nProcs = mpiSession.getNProc();
    myPID = Comm.MyPID();
#else
    nProcs = 1;
    myPID = 0;
#endif

    if (myPID == 0)
    {
        cout <<"Parallel execution: nProcs="<< nProcs << endl;
    }

    // =================== Read input xml file =============================
    pLUList = Teuchos::getParametersFromXmlFile(ipFileName);
    isoList = pLUList->sublist("Isorropia Input");
    shyLUList = pLUList->sublist("ShyLU Input");
    shyLUList.set("Outer Solver Library", "Belos");
    // Get matrix market file name
    string MMFileName = Teuchos::getParameter<string>(*pLUList, "mm_file");
    string prec_type = Teuchos::getParameter<string>(*pLUList, "preconditioner");
    int maxiters = Teuchos::getParameter<int>(*pLUList, "Outer Solver MaxIters");
    MT tol = Teuchos::getParameter<double>(*pLUList, "Outer Solver Tolerance");
    string rhsFileName = pLUList->get<string>("rhs_file", "");


    int maxFiles = pLUList->get<int>("Maximum number of files to read in", 1);
    int startFile = pLUList->get<int>("Number of initial file", 1);
    int file_number = startFile;

    if (myPID == 0)
    {
        cout << "Input :" << endl;
        cout << "ParaLU params " << endl;
        pLUList->print(std::cout, 2, true, true);
        cout << "Matrix market file name: " << MMFileName << endl;
    }

    if (maxFiles > 1)
    {
        MMFileName += "%d.mm";
        sprintf( file_name, MMFileName.c_str(), file_number );
    }
    else
    {
        strcpy( file_name, MMFileName.c_str());
    }

    // ==================== Read input Matrix ==============================
    Epetra_CrsMatrix *A;
    Epetra_MultiVector *b1;

    int err = EpetraExt::MatrixMarketFileToCrsMatrix(file_name, Comm, A);
    if (err != 0 && myPID == 0)
      {
        cout << "Matrix file could not be read in!!!, info = "<< err << endl;
        success = false;
      }

    int n = A->NumGlobalRows();

    // ==================== Read input rhs  ==============================
    if (rhsFileName != "" && maxFiles > 1)
    {
        rhsFileName += "%d.mm";
        sprintf( file_name, rhsFileName.c_str(), file_number );
    }
    else
    {
        strcpy( file_name, rhsFileName.c_str());
    }

    Epetra_Map vecMap(n, 0, Comm);
    bool allOneRHS = false;
    if (rhsFileName != "")
    {
        err = EpetraExt::MatrixMarketFileToMultiVector(file_name, vecMap, b1);
    }
    else
    {
        b1 = new Epetra_MultiVector(vecMap, 1, false);
        b1->Random();
        allOneRHS = true;
    }

    Epetra_MultiVector x(vecMap, 1);

    // Partition the matrix with hypergraph partitioning and redisstribute
    Isorropia::Epetra::Partitioner *partitioner = new
                            Isorropia::Epetra::Partitioner(A, isoList, false);

    Teuchos::Time ptime("Partition time");
    ptime.start();
    partitioner->partition();
    ptime.stop();
    if (myPID == 0)
    {
        cout << "Time to partition   : " << ptime.totalElapsedTime() << endl << endl;
    }

    Teuchos::Time rtime("Redistribute time");
    Isorropia::Epetra::Redistributor rd(partitioner);

    Epetra_CrsMatrix *newA;
    Epetra_MultiVector *newX, *newB;
    rd.redistribute(*A, newA);
    delete A;
    A = newA;
    RCP<Epetra_CrsMatrix> rcpA(A, false);

    rtime.start();
    rd.redistribute(x, newX);
    rd.redistribute(*b1, newB);
    rtime.stop();

    delete b1;
    RCP<Epetra_MultiVector> rcpx (newX, false);
    RCP<Epetra_MultiVector> rcpb (newB, false);
    //OPT::Apply(*rcpA, *rcpx, *rcpb );
    if (myPID == 0)
    {
        cout << "Time to redistribute: " << rtime.totalElapsedTime() << endl << endl;
    }

    Epetra_CrsMatrix *iterA = 0;
    Epetra_CrsMatrix *redistA = 0;
    Epetra_MultiVector *iterb1 = 0;
    Ifpack_Preconditioner *prec = nullptr;
    ML_Epetra::MultiLevelPreconditioner *MLprec = nullptr;
    while(file_number < maxFiles+startFile)
    {

        if (prec_type.compare("ShyLU") == 0)
        {
            if (file_number == startFile)
            {
                Teuchos::Time itime("Initialize time");
                itime.start();
                prec = new Ifpack_ShyLU(A);
#ifdef HAVE_IFPACK_DYNAMIC_FACTORY
                Teuchos::ParameterList shyluParameters;
                shyluParameters.set<Teuchos::ParameterList>("ShyLU list", shyLUList);
                prec->SetParameters(shyluParameters);
#else
                prec->SetParameters(shyLUList);
#endif
                prec->Initialize();
                itime.stop();
                if (myPID == 0)
                {
                    cout << "Time to initialiize : " << itime.totalElapsedTime() << endl << endl;
                }
            }
            Teuchos::Time ftime("Compute    time");
            ftime.start();
            prec->Compute();
            ftime.stop();
            if (myPID == 0)
            {
                cout << "Time to compute     : " << ftime.totalElapsedTime() << endl << endl;
            }
            //cout << " Going to set it in solver" << endl ;
            //solver.SetPrecOperator(prec);
            //cout << " Done setting the solver" << endl ;
        }
        else if (prec_type.compare("ILU") == 0)
        {
            prec = new Ifpack_ILU(A);
            prec->Initialize();
            prec->Compute();
            //solver.SetPrecOperator(prec);
        }
        else if (prec_type.compare("ILUT") == 0)
        {
            prec = new Ifpack_ILUT(A);
            prec->Initialize();
            prec->Compute();
            //solver.SetPrecOperator(prec);
        }
        else if (prec_type.compare("ML") == 0)
        {
            Teuchos::ParameterList mlList; // TODO : Take it from i/p
            MLprec = new ML_Epetra::MultiLevelPreconditioner(*A, mlList, true);
            //solver.SetPrecOperator(MLprec);
        }

        RCP<Ifpack_Preconditioner> rcpPrec(prec, false);
        RCP<Belos::EpetraPrecOp> belosPrec = rcp(new Belos::EpetraPrecOp(rcpPrec));

        const int NumGlobalElements = rcpb->GlobalLength();
        Teuchos::ParameterList belosList;
         //belosList.set( "Flexible Gmres", true );
        belosList.set( "Num Blocks", maxsubspace );// Maximum number of blocks in Krylov factorization
        belosList.set( "Block Size", blocksize );  // Blocksize to be used by iterative solver
        belosList.set( "Maximum Iterations", maxiters ); // Maximum number of iterations allowed
        belosList.set( "Maximum Restarts", maxrestarts );// Maximum number of restarts allowed
        belosList.set( "Convergence Tolerance", tol );   // Relative convergence tolerance requested
        if (numrhs > 1) {
        belosList.set( "Show Maximum Residual Norm Only", true );  // Show only the maximum residual norm
        }
        if (verbose) {
        belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
               Belos::TimingDetails + Belos::StatusTestDetails );
        if (frequency > 0)
          belosList.set( "Output Frequency", frequency );
        }
        else
        belosList.set( "Verbosity", Belos::Errors + Belos::Warnings );
        //
        // *******Construct a preconditioned linear problem********
        //

        rcpx->PutScalar(0.0);
        RCP<Belos::LinearProblem<double,MV,OP> > problem
        = rcp( new Belos::LinearProblem<double,MV,OP>( rcpA, rcpx, rcpb ) );
        if (leftprec) {
            problem->setLeftPrec( belosPrec );
        }
        else {
            problem->setRightPrec( belosPrec );
        }
        bool set = problem->setProblem();
        if (set == false) {
            if (proc_verbose)
            {
                cout << endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << endl;
            }
            cout << fail << endl;
            success = false;
            return -1;
        }

        // Create an iterative solver manager.
        RCP< Belos::SolverManager<double,MV,OP> > solver
        = rcp( new Belos::BlockGmresSolMgr<double,MV,OP>(problem,
                rcp(&belosList,false)));

        //
        // *******************************************************************
        // *************Start the block Gmres iteration*************************
        // *******************************************************************
        //
        if (proc_verbose && myPID == 0)
        {
            cout << std::endl << std::endl;
            cout << "Dimension of matrix: " << NumGlobalElements << endl;
            cout << "Number of right-hand sides: " << numrhs << endl;
            cout << "Block size used by solver: " << blocksize << endl;
            cout << "Number of restarts allowed: " << maxrestarts << endl;
            cout << "Max number of Gmres iterations per restart cycle: " <<
                        maxiters << endl;
            cout << "Relative residual tolerance: " << tol << endl;
            cout << endl;
        }

        if(tol > 1e-5)
          {
            success = false;
          }



        //
        // Perform solve
        //
        Teuchos::Time stime("Solve      time");
        stime.start();
        solver->solve ();
        stime.stop();
        if (myPID == 0)
        {
            cout << "Time to solve       : " << stime.totalElapsedTime() << endl << endl;
        }

        //
        // Get the number of iterations for this solve.
        //
        int numIters = solver->getNumIters();
        if (proc_verbose && myPID == 0)
        {
            cout << "Number of iterations performed for this solve: " <<
                     numIters << endl;
        }
        //
        // Compute actual residuals.
        //
        //bool badRes = false; // unused
        std::vector<double> actual_resids( numrhs );
        std::vector<double> rhs_norm( numrhs );
        Epetra_MultiVector resid((*rcpA).RowMap(), numrhs);
        OPT::Apply( *rcpA, *rcpx, resid );
        MVT::MvAddMv( -1.0, resid, 1.0, *rcpb, resid );
        MVT::MvNorm( resid, actual_resids );
        MVT::MvNorm( *rcpb, rhs_norm );
        if (proc_verbose && myPID == 0)
        {
            cout<< "------ Actual Residuals (normalized) -------"<<endl;
            for ( int i=0; i<numrhs; i++)
            {
                double actRes = actual_resids[i]/rhs_norm[i];
                std::cout<<"Problem "<<i<<" : \t"<< actRes;
                if (actRes > tol) {
                  //badRes = true; // unused
                  std::cout<<" (NOT CONVERGED)"<< std::endl;
                  success = false;
                } else {
                  std::cout<<" (CONVERGED)"<< std::endl;
                }
            }
        }

        file_number++;
        if (file_number >= maxFiles+startFile)
        {
          break;
        }
        else
        {
            sprintf(file_name, MMFileName.c_str(), file_number);

            if (redistA != NULL) delete redistA;
            // Load the new matrix
            err = EpetraExt::MatrixMarketFileToCrsMatrix(file_name,
                            Comm, iterA);
            if (err != 0)
            {
                if (myPID == 0)
                  {
                    cout << "Could not open file: "<< file_name << endl;

                  }
                success = false;
            }
            else
            {
                rd.redistribute(*iterA, redistA);
                delete iterA;
                InitMatValues(*redistA, A);
            }

            // Load the new rhs
            if (!allOneRHS)
            {
                sprintf(file_name, rhsFileName.c_str(), file_number);

                if (iterb1 != NULL) delete iterb1;
                err = EpetraExt::MatrixMarketFileToMultiVector(file_name,
                        vecMap, b1);
                if (err != 0)
                {
                    if (myPID==0)
                      {
                        cout << "Could not open file: "<< file_name << endl;
                        success = false;
                      }
                }
                else
                {
                    rd.redistribute(*b1, iterb1);
                    delete b1;
                    InitMVValues( *iterb1, newB );
                }
            }
        }
    }
    if (myPID == 0)
    {
        if(success)
          {
            cout << pass << endl;
          }
        else
          {
            cout << fail << endl;
          }
    }

    if (redistA != NULL) delete redistA;
    if (iterb1 != NULL) delete iterb1;


    if (prec_type.compare("ML") == 0)
    {
        delete MLprec;
    }
    else
    {
        delete prec;
    }
    delete newX;
    delete newB;
    delete A;
    delete partitioner;

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
