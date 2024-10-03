// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file shylu_util.cpp

    \brief Utilities for ShyLU

    \author Siva Rajamanickam

*/

#include <assert.h>
#include <fstream>

#include "shylu_util.h"

#ifdef HAVE_SHYLU_DDCORE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"

//Teuchos includes
#include "Teuchos_XMLParameterListHelpers.hpp"

//Isorropia includes
#include "Isorropia_Epetra.hpp"
#include "Isorropia_EpetraRedistributor.hpp"
#include "Isorropia_EpetraPartitioner.hpp"





// Currently takes onle MpiComm
Epetra_CrsMatrix *balanceAndRedistribute(Epetra_CrsMatrix *A,
                Teuchos::ParameterList isoList)
{
    // int myPID = A->Comm().MyPID(); // unused

    // Debug [
    Epetra_Map ARowMap = A->RowMap();
    // int nrows = ARowMap.NumMyElements(); // unused
    // int *rows = ARowMap.MyGlobalElements(); // unused
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
    // int ncols = AColMap.NumMyElements(); // unused
    // int *cols = AColMap.MyGlobalElements(); // unused

    // Get domain map
    Epetra_Map ADomainMap =  A->DomainMap();

    int nelems = ADomainMap.NumMyElements();
#ifndef NDEBUG
    // mfh 25 May 2015: Only used in an assert() below.
    // assert() is defined to nothing in a release build.
    int *dom_cols = ADomainMap.MyGlobalElements();
#endif // NDEBUG

    // Get range map
    Epetra_Map ARangeMap =  A->RangeMap();
    // int npts = ARangeMap.NumMyElements(); // unused
#ifndef NDEBUG
    // mfh 25 May 2015: Only used in an assert() below.
    // assert() is defined to nothing in a release build.
    int *ran_cols = ARangeMap.MyGlobalElements();
#endif // NDEBUG

    // Get row map
    Epetra_Map ARowMap = A->RowMap();
    // int nrows = ARowMap.NumMyElements(); // unused
#ifndef NDEBUG
    // mfh 25 May 2015: Only used in an assert() below.
    // assert() is defined to nothing in a release build.
    int *rows = ARowMap.MyGlobalElements();
#endif // NDEBUG

    //cout <<"In PID ="<< A->Comm().MyPID() <<" #cols="<< ncols << " #rows="<<
        //nrows <<" #domain elems="<< nelems <<" #range elems="<< npts << endl;
    // See if domain map == range map == row map
    for (int i = 0; i < nelems ; i++)
    {
        // Will this always be the case ? We will find out if assertion fails !
        assert(dom_cols[i] == ran_cols[i]);
        assert(rows[i] == ran_cols[i]);
    }
}

// TODO: SNumGlobalCols never used
void findLocalColumns(Epetra_CrsMatrix *A, int *gvals, int &SNumGlobalCols)
{

    int n = A->NumGlobalRows();
    // Get column map
    Epetra_Map AColMap = A->ColMap();
    int ncols = AColMap.NumMyElements();
    int *cols = AColMap.MyGlobalElements(); // TODO : Indexing using GID !

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

    SNumGlobalCols = 0;
    for (int i = 0; i < n ; i++)
    {
        //cout << gvals[i] ;
        if (gvals[i] > 1)
            SNumGlobalCols++;
    }
    //cout << endl;
    //cout << "Snum Global cols=" << SNumGlobalCols << endl;

    delete[] vals;
    return;
}

// This function uses a very simple tie-breaking heuristic to find a
// "narrow" separator from a wide separator. The vertices in the proc with
// smaller procID will become part of the separator
// This is not a true narrow separator, which needs a vertex cover algorithm.
// This is like a medium separator !
// TODO : This assumes symmetry I guess, Check
void findNarrowSeparator(Epetra_CrsMatrix *A, int *gvals)
{
    int nentries;
    double *values;
    int *indices;
    int n = A->NumGlobalRows();

    int myPID = A->Comm().MyPID();

    // Get row map
    Epetra_Map rMap = A->RowMap();
    Epetra_Map cMap = A->ColMap();
    int *rows = rMap.MyGlobalElements();
    int relems = rMap.NumMyElements();

    int *vals = new int[n];       // vector of size n, not ncols
    int *allGIDs = new int[n];       // vector of size n, not ncols
    for (int i = 0; i < n ; i++) // initialize to zero
    {
        vals[i] = 0;
    }

    // Rows are uniquely owned, so this will work
    for (int i = 0; i < relems ; i++)
    {
        vals[rows[i]] = myPID;        // I own relems[i]
    }


    // **************** Collective communication **************
    // This will not scale well for very large number of nodes
    // But on the node this should be fine
    A->Comm().SumAll(vals, allGIDs, n);

    // At this point all procs know who owns what rows
    for (int i = 0; i < n ; i++) // initialize to zero
        vals[i] = 0;

    int gid, cgid;
    for (int i = 0; i < relems; i++)
    {
        gid = rows[i];
        //cout << "PID=" << myPID << " " << "rowid=" << gid ;
        if (gvals[gid] != 1)
        {
            //cout << " in the sep ";
            bool movetoBlockDiagonal = false;
            // mfh 25 May 2015: This call used to assign its (int)
            // return value to 'err'.  I got rid of this, because
            // 'err' was unused.  This resulted in a "set but unused
            // variable" warning.
            (void) A->ExtractMyRowView(i, nentries, values, indices);
            //cout << " with nentries= "<< nentries;

            assert(nentries != 0);
            for (int j = 0; j < nentries; j++)
            {
                cgid = cMap.GID(indices[j]);
                assert(cgid != -1);
                if (gvals[cgid] == 1 || allGIDs[cgid] == myPID)
                    continue; // simplify the rest

                /*if (numProcs == 2)
                {*/
                    if (allGIDs[cgid] < myPID)
                    {
                        // row cgid is owned by a proc with smaller PID
                        movetoBlockDiagonal = true;
                        //cout << "\t mving to diag because of column" << cgid;
                    }
                    else
                    {
                        // There is at least one edge from this vertex to a
                        // vertex in a proc with PID > myPID, cannot move
                        // to diagonal. This is too restrictive, but
                        // important for correctness, until we can use a
                        // vertex cover algorithm.
                        movetoBlockDiagonal = false;
                        break;
                        //cout << "\tNo problem with cgid=" << cgid << "in sep";
                    }
                /*}
                else
                {
                    if (myPID == 0 && allGIDs[cgid] == numProcs-1)
                    {
                        // row cgid is owned by a proc with smaller PID
                        movetoBlockDiagonal = true;
                        cout << "\t I AM HERE mving to diag because of column" << cgid;
                    }
                    else if (myPID == numProcs-1 && allGIDs[cgid] == 0)
                    {
                        cout << "\t I AM HERE to continue " << cgid;
                        continue;
                    }
                    else if (allGIDs[cgid] < myPID)
                    {
                        // row cgid is owned by a proc with smaller PID
                        movetoBlockDiagonal = true;
                        cout << "\t mving to diag because of column" << cgid;
                    }
                    else
                    {
                        //cout << "\tNo problem with cgid=" << cgid << "in sep";
                    }
                }*/
            }
            if (movetoBlockDiagonal)
            {
                //cout << "Moving to Diagonal";
                vals[gid] = 1;
                gvals[gid] = 1; // The smaller PIDs have to know about this
                                // change. Send the change using gvals.
            }
        }
        else
        {
            // do nothing, in the diagonal block already
            //cout << "In the diagonal block";
        }
        //cout << endl;
    }

    // Reuse allGIDs to propagate the result of moving to diagonal
    for (int i = 0; i < n ; i++) // initialize to zero
        allGIDs[i] = 0;

    A->Comm().SumAll(vals, allGIDs, n);
    for (int i = 0; i < n ; i++)
    {
        if (allGIDs[i] == 1)
        {
            // Some interface columns will have gvals[1] after this
            // as the separator is narrow now.
            gvals[i] = 1; // GIDs as indices assumption
        }
    }

    delete[] vals;
    delete[] allGIDs;

}

void findBlockElems(Epetra_CrsMatrix *A, int nrows, int *rows, int *gvals,
        int Lnr, int *LeftElems,
        int Rnr, int *RightElems, std::string s1, std::string s2, bool cols)
{

    int gid;
    int rcnt = 0; int lcnt = 0;
    // Assemble ids in two arrays
    std::ostringstream ssmsg1;
    std::ostringstream ssmsg2;

#ifdef DUMP_MATRICES
    std::ostringstream fnamestr;
    fnamestr << s1 << ".mat";
    std::string fname = fnamestr.str();
    ofstream outD(fname.c_str());

    std::ostringstream fnamestrR;
    fnamestrR << s2 << ".mat";
    std::string fnameR = fnamestrR.str();
    ofstream outR(fnameR.c_str());
#endif

    ssmsg1 << s1;
    ssmsg2 << s2;
    for (int i = 0; i < nrows; i++)
    {
        gid = rows[i];
        assert (gvals[gid] >= 1);
        // If the row is local & row/column is not shared then row/column
        // belongs to D (this is not true for R, which can have more columns
        // than D)
        //if (A->LRID(gid) != -1 && gvals[gid] == 1)
        if (gvals[gid] == 1)
        {
            if (cols && A->LRID(gid) == -1) continue;
            assert(lcnt < Lnr);
            LeftElems[lcnt++] = gid;
            ssmsg1 << gid << " ";
#ifdef DUMP_MATRICES
            outD << gid << endl;
#endif
        }
        else
        {
            assert(rcnt < Rnr);
            RightElems[rcnt++] = gid;
            ssmsg2 << gid << " ";
#ifdef DUMP_MATRICES
            outR << gid << endl;
#endif
        }
    }

#ifdef DUMP_MATRICES
    outD.close();
    outR.close();
#endif

#ifdef DEBUG
    cout << ssmsg1.str() << endl;
    cout << ssmsg2.str() << endl;
#endif
    ssmsg1.clear(); ssmsg1.str("");
    ssmsg2.clear(); ssmsg2.str("");

    assert(lcnt == Lnr);
    assert(rcnt == Rnr);
}
