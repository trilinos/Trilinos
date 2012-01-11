/** \file test_orderer.cpp

    \brief Tests the Orderer interface in Isorropia

    \author Siva Rajamanickam

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

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"

// EpetraExt includes
#include "EpetraExt_CrsMatrixIn.h"

#include "Isorropia_EpetraOrderer.hpp"


using namespace std;

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, 0);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    using Teuchos::RCP;
    using Teuchos::rcp;

    char file_name[100];

    int nProcs, myPID ;
    Teuchos::ParameterList isoList ;        // Isorropia parameters
    Teuchos::ParameterList driverList ;        // Driver parameters
    string ipFileName = "orderer.xml";       // TODO : Accept as i/p

    nProcs = mpiSession.getNProc();
    myPID = Comm.MyPID();

    if (myPID == 0)
    {
        cout <<"Parallel execution: nProcs="<< nProcs << endl;
    }

    // =================== Read input xml file =============================
    Teuchos::updateParametersFromXmlFile(ipFileName, &driverList);
    isoList = driverList.sublist("Isorropia Input");

    // Get matrix market file name
    string MMFileName = Teuchos::getParameter<string>(driverList, "mm_file");


    if (myPID == 0)
    {
        cout << "Input :" << endl;
        cout << "Driver params " << endl;
        driverList.print(std::cout, 2, true, true);
        cout << "Matrix market file name: " << MMFileName << endl;
    }

    sprintf(file_name, MMFileName.c_str());

    // ==================== Read input Matrix ==============================
    Epetra_CrsMatrix *A;

    int err = EpetraExt::MatrixMarketFileToCrsMatrix(file_name, Comm, A);
    if (err != 0 && myPID == 0)
        cout << "Matrix file could not be read in!!!, info = "<< err << endl;

    // Order the matrix
    Isorropia::Epetra::Orderer *orderer = new
                            Isorropia::Epetra::Orderer(A, isoList, false);
    orderer->order();

    const int *perm;
    int perm_size;
    orderer->extractPermutationView(perm_size, perm);

    cout << "Permutation in PID=" << myPID << " is: ";
    for (int i = 0; i < perm_size; i++)
        cout << perm[i] << " ";
    cout << endl;

    delete A;
    delete orderer;

#endif
    return 0;
}
