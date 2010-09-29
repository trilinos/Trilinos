
/** \file hyperlu.cpp

    \brief Factors and solves a sparse matrix using LU factorization.

    \author Siva Rajamanickam

    This file is NOT yet complete. This is a place holder for the method
    that iterates only on the Schur complement. Lot of TODOs.


*/

#include <assert.h>
#include <iostream>
#include <sstream>

/* This define will make S block diagonal, To make this happen even some
   zero columns/rows of C and R will be stored.
   */
#define BLOCK_DIAGONAL_Si

#include "TrilinosCouplings_config.h" // Just for HAVE_MPI

// Epetra includes
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Time.h" 
#include "Epetra_CrsMatrix.h" 
#include "Epetra_Map.h" 
#include "Epetra_MultiVector.h" 
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h" 
#include "Epetra_Export.h" 

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_CrsMatrixIn.h"

// Amesos includes
#include "Amesos.h"
#include "Amesos_BaseSolver.h"

#include "hyperlu.h"
#include "hyperlu_util.h"

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
    string ipFileName = "ParaLU.xml";       // TODO : Accept as i/p
    int sym = 1; // TODO: Need rectangular factorization for unsymmetric case

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

    if (myPID == 0)
    {
        cout << "Input :" << endl;
        cout << "ParaLU params " << endl;
        pLUList.print(std::cout, 2, true, true);
        cout << "Matrix market file name: " << MMFileName << endl;
    }



    // ==================== Read input Matrix ==============================
    Epetra_CrsMatrix *tempA;
    EpetraExt::MatrixMarketFileToCrsMatrix(MMFileName.c_str(), Comm, tempA);

    Epetra_CrsMatrix *A = balanceAndRedistribute(tempA, isoList);
    delete tempA;

    Epetra_MultiVector *localS;
    Epetra_MultiVector *CMV;
    Epetra_LinearProblem *LP;
    Amesos_BaseSolver *Solver;
    int Dnr, Snr;
    int *DRowElems, *SRowElems, *piv;

    HyperLU_factor(A, sym, localS, LP, Solver, Dnr, DRowElems, Snr, SRowElems,
                    piv, CMV);

    delete localS;
    delete LP;
    delete Solver;
    delete[] DRowElems;
    delete[] SRowElems;
    delete[] piv;
    delete CMV;
}
