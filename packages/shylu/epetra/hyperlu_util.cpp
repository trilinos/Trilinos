/** \file hyperlu_util.cpp

    \brief Utilities for HyperLU

    \author Siva Rajamanickam

*/

#include <assert.h>

#include "TrilinosCouplings_config.h" // Just for HAVE_MPI

#include "Epetra_MpiComm.h"
#include "Epetra_CrsMatrix.h" 

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"

//Teuchos includes
#include "Teuchos_XMLParameterListHelpers.hpp"

//Isorropia includes
#include "Isorropia_Epetra.hpp"
#include "Isorropia_EpetraRedistributor.hpp"
#include "Isorropia_EpetraPartitioner.hpp"

#include "hyperlu_util.h"


using namespace std;


// Currently takes onle MpiComm
Epetra_CrsMatrix *balanceAndRedistribute(Epetra_CrsMatrix *A, 
                Teuchos::ParameterList isoList)
{
    int myPID = A->Comm().MyPID();

    // Debug [
    Epetra_Map ARowMap = A->RowMap();
    int nrows = ARowMap.NumMyElements();
    int *rows = ARowMap.MyGlobalElements();
    // ]

    // ==================== Symbolic factorization =========================
    // 1. Partition and redistribute [
    Isorropia::Epetra::Partitioner *partitioner = new 
                            Isorropia::Epetra::Partitioner(A, isoList, false);
    partitioner->partition();

    Isorropia::Epetra::Redistributor rd(partitioner);
    Epetra_CrsMatrix *newA;
    rd.redistribute(*A, newA);
    // ]
    EpetraExt::RowMatrixToMatlabFile("A.mat", *newA);

    delete partitioner;
    return newA;
}

/* TODO : Do this only for Debug ? */
void checkMaps(Epetra_CrsMatrix *A)
{
    // Get column map
    Epetra_Map AColMap = A->ColMap();
    int ncols = AColMap.NumMyElements();
    int *cols = AColMap.MyGlobalElements();

    // Get domain map
    Epetra_Map ADomainMap =  A->DomainMap();
    int nelems = ADomainMap.NumMyElements();
    int *dom_cols = ADomainMap.MyGlobalElements();

    // Get range map
    Epetra_Map ARangeMap =  A->RangeMap();
    int npts = ARangeMap.NumMyElements();
    int *ran_cols = ARangeMap.MyGlobalElements();

    // Get row map
    Epetra_Map ARowMap = A->RowMap();
    int nrows = ARowMap.NumMyElements();
    int *rows = ARowMap.MyGlobalElements();

    cout <<"In PID ="<< A->Comm().MyPID() <<" #cols="<< ncols << " #rows="<< 
        nrows <<" #domain elems="<< nelems <<" #range elems="<< npts << endl;
    // See if domain map == range map == row map
    for (int i = 0; i < nelems ; i++)
    {
        // Will this always be the case ? We will find out if assertion fails !
        assert(dom_cols[i] == ran_cols[i]);
        assert(rows[i] == ran_cols[i]);
    }
}

void findLocalColumns(Epetra_CrsMatrix *A, int *gvals)
{

    int n = A->NumGlobalRows();
    // Get column map
    Epetra_Map AColMap = A->ColMap();
    int ncols = AColMap.NumMyElements();
    int *cols = AColMap.MyGlobalElements();

    // 2. Find column permutation [
    // Find all columns in this proc
    int *vals = new int[n];       // vector of size n, not ncols
    for (int i = 0; i < n ; i++) 
    {
        vals[i] = 0;
        gvals[i] = 0;
    }

    // Null columns in A are not part of any proc
    for (int i = 0; i < ncols ; i++)
    {
        vals[cols[i]] = 1;        // Set to 1 for locally owned columns
    }

    // Bottleneck?: Compute the column permutation
    A->Comm().SumAll(vals, gvals, n);

    delete vals;

}

void findBlockElems(int nrows, int *rows, int *gvals, int Lnr, int *LeftElems, 
        int Rnr, int *RightElems, string s1, string s2)
{
 
    int gid;
    int rcnt = 0, lcnt = 0;
    // Assemble ids in two arrays
    ostringstream ssmsg1;
    ostringstream ssmsg2;

    ssmsg1 << s1;
    ssmsg2 << s2;
    for (int i = 0; i < nrows; i++)
    {
        gid = rows[i];
        assert (gvals[gid] >= 1);
        if (gvals[gid] == 1)
        {
            LeftElems[lcnt++] = gid;
            ssmsg1 << gid << " ";
        }
        else
        {
            RightElems[rcnt++] = gid; 
            ssmsg2 << gid << " ";
        }
    }
    //cout << ssmsg1.str() << endl;
    //cout << ssmsg2.str() << endl; // TODO: Enable it only in debug mode
    ssmsg1.clear(); ssmsg1.str("");
    ssmsg2.clear(); ssmsg2.str("");

    assert(lcnt == Lnr);
    assert(rcnt == Rnr);
}
