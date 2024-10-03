// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file amesos_driver.cpp

    \brief Factors and solves a sparse matrix using LU factorization.

    \author Siva Rajamanickam

    \remark Usage:
    \code mpirun -n np amesos_driver.exe

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
#include "Teuchos_Time.hpp"

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
    // Get matrix market file name
    string MMFileName = Teuchos::getParameter<string>(pLUList, "mm_file");
    string prec_type = Teuchos::getParameter<string>(pLUList, "preconditioner");

    if (myPID == 0)
    {
        cout << "Input :" << endl;
        cout << "ParaLU params " << endl;
        pLUList.print(std::cout, 2, true, true);
        cout << "Matrix market file name: " << MMFileName << endl;
    }

    // ==================== Read input Matrix ==============================
    Epetra_CrsMatrix *A;

    int err = EpetraExt::MatrixMarketFileToCrsMatrix(MMFileName.c_str(), Comm, A);
    //EpetraExt::MatlabFileToCrsMatrix(MMFileName.c_str(), Comm, A);
    //assert(err != 0);
    cout <<"Done reading the matrix"<< endl;
    int n = A->NumGlobalRows();
    cout <<"n="<< n << endl;

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

    rd.redistribute(x, newX);
    rd.redistribute(b, newB);

    Epetra_LinearProblem problem(A, newX, newB);

    Amesos Factory;
    char* SolverType = "Amesos_Klu";
    bool IsAvailable = Factory.Query(SolverType);

    Epetra_LinearProblem *LP = new Epetra_LinearProblem();
    LP->SetOperator(A);
    LP->SetLHS(newX);
    LP->SetRHS(newB);
    Amesos_BaseSolver *Solver = Factory.Create(SolverType, *LP);


    Solver->SymbolicFactorization();
  Teuchos::Time ftime("setup time");
      ftime.start();
    Solver->NumericFactorization();
    cout << "Numeric Factorization" << endl;
    Solver->Solve();
    cout << "Solve done" << endl;

    ftime.stop();
    cout << "Time to setup" << ftime.totalElapsedTime() << endl;

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
    delete newX;
    delete newB;
    delete A;
    delete partitioner;
}
