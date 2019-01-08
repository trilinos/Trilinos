//
//  main.cpp
//  
//
//  Created by Friederike Röver on 18.12.18.
//

#include "main.hpp"
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_SerialDenseMatrix.h"

using namespace std;
// function declaration
void compute_loc_matrix( double *x_triangle, double *y_triangle,
                        Epetra_SerialDenseMatrix &Ke );
int find( const int list[], const int length, const int index);
// main driver
int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm Comm;
#endif
    if (Comm.NumProc() != 16) {
#ifdef HAVE_MPI
        MPI_Finalize();
#endif
        return(0);
    }
    
    Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm = Teuchos::rcp(new Teuchos::MpiComm<int> (MPI_COMM_WORLD));

    int NumProcs = TeuchosComm->getSize(); //Number of Procs->16 required;
    int nInRow = 9; // Numer of nodes in a row;
    
    
    int NumMyElements = 0;         // NODES assigned to this processor
    int NumMyExternalElements = 0; // nodes used by this proc, but not hosted
    int NumMyTotalElements = 0;
    int FE_NumMyElements = 0;      // TRIANGLES assigned to this processor
    int * MyGlobalElements = 0;    // nodes assigned to this processor
    Epetra_IntSerialDenseMatrix T; // store the grid connectivity
    
    int nline = (nInRow-1)/sqrt(NumProcs);
    
    
    
}
