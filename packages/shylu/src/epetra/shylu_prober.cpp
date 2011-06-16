/* Prober example */

#include <Isorropia_config.h>

#if defined(HAVE_MPI) && defined(HAVE_EPETRA)
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_CrsMatrix.h" 
#include "Epetra_Import.h" 
#include "Epetra_BlockMap.h" 

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"

// Isorropia includes
#include <Isorropia_EpetraProber.hpp>

int main(int argc, char *argv[])
{
#if defined(HAVE_MPI) && defined(HAVE_EPETRA)
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, 0);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);

    int rank, MyPID, NumProc ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MyPID = Comm.MyPID();
    NumProc = Comm.NumProc();

    int NumMyEquations = 5;
    int NumGlobalEquations = NumMyEquations * NumProc;

    Epetra_Map Map(NumGlobalEquations, NumMyEquations, 0, Comm);
    int NumMyElements = Map.NumMyElements();

    std::vector<int> MyGlobalElements(NumMyElements);
    Map.MyGlobalElements(&MyGlobalElements[0]);

    std::vector<int> NumNz(NumMyElements);

    int i, ierr;
    for (i=0; i<NumMyElements; i++)
    {
        if (MyGlobalElements[i]==0 || 
                MyGlobalElements[i] == NumGlobalEquations-1)
            NumNz[i] = 2;
        else
            NumNz[i] = 3;
    }

    // Create a tridiagonal Epetra_Matrix
    Epetra_CrsMatrix A(Copy, Map, &NumNz[0]);

    // Add  rows one-at-a-time
    // Off diagonal Values will always be -1
    std::vector<double> Values(2);
    Values[0] = -1.0; Values[1] = -1.0;
    std::vector<int> Indices(2);
    double two = 2.0;
    int NumEntries;
      
    for (i=0; i<NumMyElements; i++)
    {
        if (MyGlobalElements[i]==0)
        {
            Indices[0] = 1;
            NumEntries = 1;
        }
        else if (MyGlobalElements[i] == NumGlobalEquations-1)
        {
            Indices[0] = NumGlobalEquations-2;
            NumEntries = 1;
        }
        else
        {
            Indices[0] = MyGlobalElements[i]-1;
            Indices[1] = MyGlobalElements[i]+1;
            NumEntries = 2;
        }
        ierr = A.InsertGlobalValues(MyGlobalElements[i], NumEntries,
                        &Values[0], &Indices[0]);
        assert(ierr==0);
        // Put in the diagonal entry
        ierr = A.InsertGlobalValues(MyGlobalElements[i], 1, &two,
                        &MyGlobalElements[i]);
        assert(ierr==0);
    }
       
    // Finish up matrix construction
    ierr = A.FillComplete();
    assert(ierr==0);

    cout << A << endl;

    // Construct the graph for a diagonal matrix
    for (i=0; i<NumMyElements; i++)
    {
        if (MyGlobalElements[i]==0 || 
                MyGlobalElements[i] == NumGlobalEquations-1)
            NumNz[i] = 1;
        else
            NumNz[i] = 2;
    }
    //Epetra_CrsGraph G1 (Copy, Map, 1) ;
    Epetra_CrsGraph G1 (Copy, Map, &NumNz[0]) ;

    // Insert the indices
    for (i=0; i<NumMyElements; i++)
    {
        if (MyGlobalElements[i]==0)
        {
            Indices[0] = 1;
            NumEntries = 1;
        }
        else if (MyGlobalElements[i] == NumGlobalEquations-1)
        {
            Indices[0] = NumGlobalEquations-2;
            NumEntries = 1;
        }
        else
        {
            Indices[0] = MyGlobalElements[i]-1;
            Indices[1] = MyGlobalElements[i]+1;
            NumEntries = 2;
        }
        ierr = G1.InsertGlobalIndices(MyGlobalElements[i], NumEntries,
                        &Indices[0]);
        assert(ierr==0);
        // Put in the diagonal entry
        //G1.InsertGlobalIndices(MyGlobalElements[i], 1, &MyGlobalElements[i]);
        //assert(ierr==0);
    }

    G1.FillComplete();
    Teuchos::RCP<const Epetra_CrsGraph> RCPG1 = Teuchos::rcpFromRef(G1);

    //Set up a prober
    Teuchos::ParameterList pList;
    Isorropia::Epetra::Prober prober(RCPG1, pList, false);

    cout << "Created prober" << endl;
    cout << G1 << endl;

    prober.color();

    cout << "Importer = " << (G1.Importer())->TargetMap().MinMyGID() ;
    cout << "Done Coloring" << endl;
    Teuchos::RCP<Epetra_CrsMatrix> D = prober.probe(A);
    cout << "Done Probing" << endl;

    cout << *D << endl;

#else
    printf("Need MPI and Epetra !!!! \n")
    return 1;
#endif
}
