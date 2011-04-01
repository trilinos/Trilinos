/** \file hyperlu_driver.cpp

    \brief Factors and solves a sparse matrix using LU factorization.

    \author Siva Rajamanickam

    \remark Usage:
    \code mpirun -n np hyperlu_driver.exe

*/

#include <assert.h>
#include <iostream>
#include <sstream>

#include "Isorropia_config.h" // Just for HAVE_MPI

// Epetra includes
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
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

// Amesos includes
#include "Amesos.h"
#include "Amesos_BaseSolver.h"

// AztecOO includes
#include "AztecOO.h"

// ML includes
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

#include "hyperlu.h"
#include "hyperlu_util.h"
#include "Ifpack_HyperLU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_ILU.h"

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

#include "EpetraExt_readEpetraLinearSystem.h"

//TODO Combine Belos driver and AztecOO driver into one driver.


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


    bool verbose = false, proc_verbose = true;
    bool leftprec = true;      // left preconditioning or right.
    int frequency = -1;        // frequency of status test output.
    int blocksize = 1;         // blocksize
    int numrhs = 1;            // number of right-hand sides to solve for
    int maxrestarts = 15;      // maximum number of restarts allowed 
    int maxsubspace = 25;      // maximum number of blocks the solver can use 
                               // for the subspace


    int nProcs, myPID ;
    Teuchos::ParameterList pLUList ;        // ParaLU parameters
    Teuchos::ParameterList isoList ;        // Isorropia parameters
    Teuchos::ParameterList hyperLUList ;    // HyperLU parameters
    string ipFileName = "HyperLU.xml";       // TODO : Accept as i/p

    nProcs = mpiSession.getNProc();
    myPID = Comm.MyPID();

    if (myPID == 0)
    {
        cout <<"Parallel execution: nProcs="<< nProcs << endl;
    }

    // =================== Read input xml file =============================
    Teuchos::updateParametersFromXmlFile(ipFileName, &pLUList);
    isoList = pLUList.sublist("Isorropia Input");
    hyperLUList = pLUList.sublist("HyperLU Input");
    hyperLUList.set("Outer Solver Library", "Belos");
    // Get matrix market file name
    string MMFileName = Teuchos::getParameter<string>(pLUList, "mm_file");
    string prec_type = Teuchos::getParameter<string>(pLUList, "preconditioner");
    int maxiters = Teuchos::getParameter<int>(pLUList, "Outer Solver MaxIters");
    MT tol = Teuchos::getParameter<double>(pLUList, "Outer Solver Tolerance");

    if (myPID == 0)
    {
        cout << "Input :" << endl;
        cout << "ParaLU params " << endl;
        pLUList.print(std::cout, 2, true, true);
        cout << "Matrix market file name: " << MMFileName << endl;
    }

    // ==================== Read input Matrix ==============================
    Epetra_CrsMatrix *A;

    int err = EpetraExt::MatrixMarketFileToCrsMatrix(MMFileName.c_str(), Comm, 
                                    A);
    //cout <<"Done reading the matrix"<< endl;
    int n = A->NumGlobalRows();
    //cout <<"n="<< n << endl;

    // Create input vectors
    Epetra_Map vecMap(n, 0, Comm);
    Epetra_MultiVector x(vecMap, 1);
    Epetra_MultiVector b(vecMap, 1, false);

    b.PutScalar(1.0); // TODO : Accept it as input

    // Partition the matrix with hypergraph partitioning and redisstribute
    Isorropia::Epetra::Partitioner *partitioner = new
                            Isorropia::Epetra::Partitioner(A, isoList, false);
    partitioner->partition();
    Isorropia::Epetra::Redistributor rd(partitioner);

    Epetra_CrsMatrix *newA;
    Epetra_MultiVector *newX, *newB; 
    rd.redistribute(*A, newA);
    delete A;
    A = newA;
    RCP<Epetra_CrsMatrix> rcpA(A, false);

    rd.redistribute(x, newX);
    rd.redistribute(b, newB);
    RCP<Epetra_MultiVector> rcpx (newX, false);
    RCP<Epetra_MultiVector> rcpb (newB, false);
    //OPT::Apply(*rcpA, *rcpx, *rcpb );


    Ifpack_Preconditioner *prec;
    ML_Epetra::MultiLevelPreconditioner *MLprec;
    if (prec_type.compare("HyperLU") == 0)
    {
        prec = new Ifpack_HyperLU(A);
        prec->SetParameters(hyperLUList);
        prec->Initialize();
        prec->Compute();
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
      std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return -1;
    }

    // Create an iterative solver manager.
    RCP< Belos::SolverManager<double,MV,OP> > solver
    = rcp( new Belos::BlockGmresSolMgr<double,MV,OP>(problem, rcp(&belosList,false)));

    //
    // *******************************************************************
    // *************Start the block Gmres iteration*************************
    // *******************************************************************
    //
    if (proc_verbose) {
    std::cout << std::endl << std::endl;
    std::cout << "Dimension of matrix: " << NumGlobalElements << std::endl;
    std::cout << "Number of right-hand sides: " << numrhs << std::endl;
    std::cout << "Block size used by solver: " << blocksize << std::endl;
    std::cout << "Number of restarts allowed: " << maxrestarts << std::endl;
    std::cout << "Max number of Gmres iterations per restart cycle: " << maxiters << std::endl;
    std::cout << "Relative residual tolerance: " << tol << std::endl;
    std::cout << std::endl;
    }
    //
    // Perform solve
    //
    Belos::ReturnType ret = solver->solve();
    //
    // Get the number of iterations for this solve.
    //
    int numIters = solver->getNumIters();
    if (proc_verbose)
    std::cout << "Number of iterations performed for this solve: " << numIters << std::endl;
    //
    // Compute actual residuals.
    //
    bool badRes = false;
                                         
    std::vector<double> actual_resids( numrhs );
    std::vector<double> rhs_norm( numrhs );
    Epetra_MultiVector resid((*rcpA).RowMap(), numrhs);
    OPT::Apply( *rcpA, *rcpx, resid );
    MVT::MvAddMv( -1.0, resid, 1.0, *rcpb, resid );
    MVT::MvNorm( resid, actual_resids );
    MVT::MvNorm( *rcpb, rhs_norm );
    if (proc_verbose) {
    std::cout<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
    for ( int i=0; i<numrhs; i++) {
      double actRes = actual_resids[i]/rhs_norm[i];
      std::cout<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
      if (actRes > tol) badRes = true;
    }
    }


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
}
