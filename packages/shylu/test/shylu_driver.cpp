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

#include "shylu.h"
#include "shylu_util.h"
#include "Ifpack_ShyLU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_ILU.h"

using namespace std;

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, 0);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif
    int nProcs, myPID ;
    Teuchos::ParameterList pLUList ;        // ParaLU parameters
    Teuchos::ParameterList isoList ;        // Isorropia parameters
    Teuchos::ParameterList shyLUList ;    // shyLU parameters
    Teuchos::ParameterList ifpackList ;    // shyLU parameters
    string ipFileName = "ShyLU.xml";       // TODO : Accept as i/p

    nProcs = mpiSession.getNProc();
    myPID = Comm.MyPID();

    if (myPID == 0)
    {
        cout <<"Parallel execution: nProcs="<< nProcs << endl;
    }

    // =================== Read input xml file =============================
    Teuchos::updateParametersFromXmlFile(ipFileName, &pLUList);
    isoList = pLUList.sublist("Isorropia Input");
    shyLUList = pLUList.sublist("ShyLU Input");
    shyLUList.set("Outer Solver Library", "AztecOO");
    // Get matrix market file name
    string MMFileName = Teuchos::getParameter<string>(pLUList, "mm_file");
    string prec_type = Teuchos::getParameter<string>(pLUList, "preconditioner");
    int maxiters = Teuchos::getParameter<int>(pLUList, "Outer Solver MaxIters");
    double tol = Teuchos::getParameter<double>(pLUList, "Outer Solver Tolerance");
    string rhsFileName = pLUList.get<string>("rhs_file", "");

    if (myPID == 0)
    {
        cout << "Input :" << endl;
        cout << "ParaLU params " << endl;
        pLUList.print(std::cout, 2, true, true);
        cout << "Matrix market file name: " << MMFileName << endl;
    }

    // ==================== Read input Matrix ==============================
    Epetra_CrsMatrix *A;
    Epetra_MultiVector *b1;

    int err = EpetraExt::MatrixMarketFileToCrsMatrix(MMFileName.c_str(), Comm,
                                                        A);
    //EpetraExt::MatlabFileToCrsMatrix(MMFileName.c_str(), Comm, A);
    //assert(err != 0);
    //cout <<"Done reading the matrix"<< endl;
    int n = A->NumGlobalRows();
    //cout <<"n="<< n << endl;

    // Create input vectors
    Epetra_Map vecMap(n, 0, Comm);
    if (rhsFileName != "")
    {
        err = EpetraExt::MatrixMarketFileToMultiVector(rhsFileName.c_str(),
                                         vecMap, b1);
    }
    else
    {
        b1 = new Epetra_MultiVector(vecMap, 1, false);
        b1->PutScalar(1.0);
    }

    Epetra_MultiVector x(vecMap, 1);
    //cout << "Created the vectors" << endl;

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

    rd.redistribute(x, newX);
    rd.redistribute(*b1, newB);

    Epetra_LinearProblem problem(A, newX, newB);

    AztecOO solver(problem);

    ifpackList ;
    Ifpack_Preconditioner *prec;
    ML_Epetra::MultiLevelPreconditioner *MLprec;
    if (prec_type.compare("ShyLU") == 0)
    {
        prec = new Ifpack_ShyLU(A);
        prec->SetParameters(shyLUList);
        prec->Initialize();
        prec->Compute();
        //(dynamic_cast<Ifpack_ShyLU *>(prec))->JustTryIt();
        //cout << " Going to set it in solver" << endl ;
        solver.SetPrecOperator(prec);
        //cout << " Done setting the solver" << endl ;
    }
    else if (prec_type.compare("ILU") == 0)
    {
        ifpackList.set( "fact: level-of-fill", 1 );
        prec = new Ifpack_ILU(A);
        prec->SetParameters(ifpackList);
        prec->Initialize();
        prec->Compute();
        solver.SetPrecOperator(prec);
    }
    else if (prec_type.compare("ILUT") == 0)
    {
        ifpackList.set( "fact: ilut level-of-fill", 2 );
        ifpackList.set( "fact: drop tolerance", 1e-8);
        prec = new Ifpack_ILUT(A);
        prec->SetParameters(ifpackList);
        prec->Initialize();
        prec->Compute();
        solver.SetPrecOperator(prec);
    }
    else if (prec_type.compare("ML") == 0)
    {
        Teuchos::ParameterList mlList; // TODO : Take it from i/p
        MLprec = new ML_Epetra::MultiLevelPreconditioner(*A, mlList, true);
        solver.SetPrecOperator(MLprec);
    }

    solver.SetAztecOption(AZ_solver, AZ_gmres);
    solver.SetMatrixName(333);
    //solver.SetAztecOption(AZ_output, 1);
    //solver.SetAztecOption(AZ_conv, AZ_Anorm);
    //cout << "Going to iterate for the global problem" << endl;

    solver.Iterate(maxiters, tol);

    // compute ||Ax - b||
    double Norm;
    Epetra_MultiVector Ax(vecMap, 1);

    Epetra_MultiVector *newAx; 
    rd.redistribute(Ax, newAx);
    A->Multiply(false, *newX, *newAx);
    newAx->Update(1.0, *newB, -1.0);
    newAx->Norm2(&Norm);
    double ANorm = A->NormOne();

    cout << "|Ax-b |/|A| = " << Norm/ANorm << endl;

    delete newAx;
    if (prec_type.compare("ML") == 0)
    {
        delete MLprec;
    }
    else
    {
        delete prec;
    }

    delete b1;
    delete newX;
    delete newB;
    delete A;
    delete partitioner;
}
