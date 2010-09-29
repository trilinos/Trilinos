/** \file parallel_LU.cpp

    \brief Factors sparse matrix using LU factorization.

    \author Siva Rajamanickam

    \remark Usage:
    \code ./parallel_LU.exe \endcode

    This version extracts the non zero rows/columns of R/C, hence Si is smaller
    the complete diagonal block. This version also
    extracts the entire As matrix in each proc. This is not a huge problem.
    When one needs vec = S * v = As * v - Si * v(corres rows) and minus updates
    the corresponding rows of vec correctly. This also needs storing both As and
    Si. The preconditioner for this method is not yet writted TODO, Currently
    keep BLOCK_DIAGONAL_Si always defined.

    When BLOCK_DIAGONAL_Si is defined:
    This version extracts all rows/columns of R/C, 
    This version also extracts the block diagonal of As matrix in each proc. 
    We store As - Si.

*/

#include <assert.h>
#include <iostream>
#include <sstream>

/* This define will make S block diagonal, To make this happen even some
   zero columns/rows of C and R will be stored.
   */
#define BLOCK_DIAGONAL_Si

// To dump all the matrices into files.
//#define DUMP_MATRICES

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


int HyperLU_factor(Epetra_CrsMatrix *A, int sym, Epetra_MultiVector *&localS,
        Epetra_LinearProblem *&LP, Amesos_BaseSolver *&Solver, int &Dnr, 
        int *&DRowElems, int &Snr, int *&SRowElems, int *&piv,
        Epetra_MultiVector *&CMV)
{
    int myPID = A->Comm().MyPID();
    int n = A->NumGlobalRows();

    checkMaps(A);

    // Get column map
    Epetra_Map AColMap = A->ColMap();
    int ncols = AColMap.NumMyElements();
    int *cols = AColMap.MyGlobalElements();

    // Get row map
    Epetra_Map ARowMap = A->RowMap();
    int nrows = ARowMap.NumMyElements();
    int *rows = ARowMap.MyGlobalElements();

    // Find all columns in this proc
    int *gvals = new int[n];       // vector of size n, not ncols !
    // gvals[local cols] = 1, gvals[shared cols] > 1.
    findLocalColumns(A, gvals);

    // 3. Assemble diagonal block and the border in convenient form [
    /* In each processor, we have (in a permuted form)
     *  | D_i    C_i   |
     *  | R_i    S_i   |
     * D_i - diagonal block, C_i - Column Separator, R_i - Row separator
     * S_i - Schur complement part of A
     * Assemble all four blocks in local matrices. */

     ostringstream ssmsg1;
     ssmsg1 << "PID =" << myPID << " ";
     string msg = ssmsg1.str();
     ssmsg1.clear(); ssmsg1.str("");

    // Find #cols in each block
    int Dnc = 0;        // #cols in diagonal block
    int Snc = 0;        // #cols in the col. separator
    for (int i = 0; i < ncols ; i++)
    {
        if (gvals[cols[i]] == 1)
            Dnc++;
        else
            Snc++;
    }

    // Find #rows in each block 
    Dnr = Dnc;          // #rows in square diagonal block
    Snr = nrows - Dnr;  // #rows in the row separator

    cout << msg << " #columns in diagonal blk ="<< Dnc << endl;
    cout << msg << " #rows in diagonal blk ="<< Dnr << endl;
#ifdef BLOCK_DIAGONAL_Si
    cout << msg << " #columns in S ="<< Snr << endl;
#else
    cout << msg << " #columns in S ="<< Snc << endl;
#endif
    cout << msg << " #rows in S ="<< Snr << endl;

    // Create a row map for the D and S blocks [
    DRowElems = new int[Dnr];
    SRowElems = new int[Snr];
    int gid;
    // Assemble row ids in two arrays (for D and R blocks)
    if (sym)
    {
        findBlockElems(nrows, rows, gvals, Dnr, DRowElems, Snr, SRowElems, 
                    "D Rows ", "S Rows") ;
    }
    else
    {
        // SRowElems are not known until factorization, TODO
        assert(0 == 1);
    }

    // Create the local row map 
    Epetra_SerialComm LComm;        // Use Serial Comm for the local blocks.
    Epetra_Map LocalDRowMap(-1, Dnr, DRowElems, 0, LComm);
    Epetra_Map LocalSRowMap(-1, Snr, SRowElems, 0, LComm);
    // ]

    // Create a column map for the D and S blocks [
    int *DColElems = new int[Dnc]; // Elems in column map of D 
    int *SColElems = new int[Snc]; // Elems in column map of C
    // Assemble column ids in two arrays (for D and C blocks)
    findBlockElems(ncols, cols, gvals, Dnc, DColElems, Snc, SColElems, 
                    "D Cols ", "S Cols") ;

    // Create the local column map 
    Epetra_Map LocalDColMap(-1, Dnc, DColElems, 0, LComm);
    Epetra_Map LocalSColMap(-1, Snc, SColElems, 0, LComm);
    for (int i = 0; i < Snr; i++)
    {
        // Epetra guarentees columns corresponding to local rows will be first
        // in the column map.
        assert(SRowElems[i] == SColElems[i]);
    }
    // ]

    // Compute nnz per row.
    int *Ai;
    double *Ax;
    int *DNumEntriesPerRow = new int[Dnr];
    int *SNumEntriesPerRow = new int[Snr];
    int Dmaxnnz=0, Cmaxnnz=0, Rmaxnnz=0, Smaxnnz=0; // Max nnz in any row

    int *CColMax = new int[n]; // Possible Elems in column map of C
    int *RRowMax = new int[n]; // Possible Elems in row map of R
    int Cnc = 0; // nonzero columns in C
    int Rnr = 0; // nonzero rows in R

    int dcol, ccol, rcol, scol;
    int *CColElems, *RRowElems;

    if (sym)
    {
        for (int i = 0; i < n; i++)
        {
            CColMax[i] = -1; RRowMax[i] = -1;
        }

        int NumEntries;
        int dcnt, scnt;
        dcol = 0 ; ccol = 0; rcol = 0; scol = 0;
        // Three things to worry about here, Num entries in each row of S,
        // size of C and R and the corresponding cols/rows.
        for (int i = 0; i < nrows; i++) 
        {
            dcnt = 0; scnt = 0;
            // Need to pass local id i to this function 
            int err = A->ExtractMyRowView(i, NumEntries, Ax, Ai);
            if (err != 0)
            { 
                cout << msg << "\t" << i << "\t" << rows[i] << endl 
                << "Trouble " << " " << err << endl; 
            }

            gid = rows[i];
            if (gvals[gid] == 1)
            { // D or C row
                int gcid;
                for (int j = 0 ; j < NumEntries ; j++)
                { // O(nnz) ! Careful what you do inside
                    gcid = A->GCID(Ai[j]); 
                    if (gvals[gcid] == 1) 
                    {
                        dcnt++;
                    }
                    else
                    { 
                        // Give a mapping from global column ids to vector ids 
                        // in C
                        if (CColMax[gcid] ==  -1)
                            CColMax[gcid] = Cnc++;
                    }
                }
                // Assign 0 and increase dcol even if it is just C row. 
                // There should not be a null row in D.
                assert(dcnt != 0);
                DNumEntriesPerRow[dcol++] = dcnt;
            }
            else
            { // R or S row
#ifdef BLOCK_DIAGONAL_Si
                if (RRowMax[gid] == -1)
                    RRowMax[gid] = Rnr++;
#endif
                int gcid;
                //cout << "proc/row " << myPID << "/" << gid;
                //cout << " has Cols = " ;
                for (int j = 0 ; j < NumEntries ; j++)
                { // O(nnz) ! Careful what you do inside
                    gcid = A->GCID(Ai[j]); 
                    if (gvals[gcid] == 1)
                    {
#ifndef BLOCK_DIAGONAL_Si
                        if (RRowMax[gid] == -1)
                            RRowMax[gid] = Rnr++;
#endif
                    }
                    else
                    {
#ifdef BLOCK_DIAGONAL_Si
                        if (A->LRID(gcid) != -1)
                        {
#endif
                            //cout << gcid << " " ;
                            scnt++;
#ifdef BLOCK_DIAGONAL_Si
                        }
#endif
                    }
                }
                //cout << " count of " << scnt << endl;
                // Assign 0 and increase scol even if it is just R row. 
                // There should not be a null row in S.
                assert(scnt != 0);
                SNumEntriesPerRow[scol++] = scnt;
            }
            Dmaxnnz = max(Dmaxnnz, dcnt);
            Smaxnnz = max(Smaxnnz, scnt); 
        }
        assert( dcol == Dnr);
        assert( scol == Snr);

        RRowElems = new int[Rnr]; // Possible Elems in row map of R
        int temp = 0 ;

#ifdef BLOCK_DIAGONAL_Si
        assert(Snr == Rnr);
        CColElems = new int[Rnr]; // Possible Elems in column map of C
        for (int i = 0; i < Snr; i++)
        {
            if (RRowMax[SRowElems[i]] != -1)
            {
                CColElems[temp++] = SRowElems[i]; // keep the gid for maps
            }
        }
        assert(temp == Snr);
#else
        CColElems = new int[Cnc]; // Possible Elems in column map of C
        for (int i = 0; i < Snc; i++)
        {
            if (CColMax[SColElems[i]] != -1)
            {
                CColElems[temp++] = SColElems[i]; // keep the gid for maps
            }
        }
#endif

        temp = 0;
        for (int i = 0; i < Snr; i++)
        {
            if (RRowMax[SRowElems[i]] != -1)
            {
                RRowElems[temp++] = SRowElems[i]; // keep the gid for maps
            }
        }
        cout << msg << " #cols in C=" << Cnc << "#rows in R=" << Rnr << endl;
    }
    else
    {
        assert (0 == 1);
        /*assert(Dnr == nrows);
        // TODO : Need to compute the entries for D and C in one array of size
        // nrow*/
    }


    int Sbd_nc; // #cols in block diagonal of S.
#ifdef BLOCK_DIAGONAL_Si
    Sbd_nc = Snr;
#else
    Sbd_nc = Cnc;
#endif
    Epetra_Map LocalCColMap(-1, Sbd_nc, CColElems, 0, LComm);

#if 0
        // Create the local column map 
    Epetra_Map LocalRColMap(-1, rcol, RColElems, 0, LComm);
#endif

    //Create the local matrices
    Epetra_CrsMatrix D(Copy, LocalDRowMap, LocalDColMap, DNumEntriesPerRow, 
                            true);
#ifdef BLOCK_DIAGONAL_Si
    Epetra_CrsMatrix S(Copy, LocalSRowMap, LocalSRowMap, SNumEntriesPerRow, 
                            true);
#else
    Epetra_CrsMatrix S(Copy, LocalSRowMap, LocalSColMap, SNumEntriesPerRow, 
                            true);
#endif
    //Epetra_CrsMatrix C(Copy, LocalDRowMap, LocalSColMap, CNumEntriesPerRow, 
                            //true);
    //Epetra_CrsMatrix R(Copy, LocalSRowMap, LocalDColMap, RNumEntriesPerRow, 
                            //true);

    CMV = new Epetra_MultiVector(LocalDRowMap, Sbd_nc);

    // Assuming D's column map is a super set of R's columns
    // This is R' and not for the solve later to work right.
    Epetra_MultiVector *RMV = new Epetra_MultiVector(LocalDColMap, Rnr); 

    // TODO : Use S or SMV, not both
    Epetra_MultiVector *SMV = new Epetra_MultiVector(LocalSRowMap, Sbd_nc); 

    int *LeftIndex = new int[max(Dmaxnnz, Rmaxnnz)];
    double *LeftValues = new double[max(Dmaxnnz, Rmaxnnz)];
    int *RightIndex = new int[max(Cmaxnnz, Smaxnnz)];
    double *RightValues = new double[max(Cmaxnnz, Smaxnnz)];

    int err;
    int lcnt, rcnt;
    int debug_row = 0;
    for (int i = 0; i < nrows ; i++)
    {
        int NumEntries;
        err = A->ExtractMyRowView(i, NumEntries, Ax, Ai);

        lcnt = 0, rcnt = 0;
        // Place the entry in the correct sub matrix, Works only for sym
        gid = rows[i];
        for (int j = 0 ; j < NumEntries ; j++)
        { // O(nnz) ! Careful what you do inside
            // Row permutation does not matter here 
            if (gvals[A->GCID(Ai[j])] == 1)
            {
                assert(lcnt < max(Dmaxnnz, Rmaxnnz));
                LeftIndex[lcnt] = A->GCID(Ai[j]);
                LeftValues[lcnt++] = Ax[j];
            }
            else
            {
#ifdef BLOCK_DIAGONAL_Si
                if (A->LRID(A->GCID(Ai[j])) != -1)
                {
#endif
                if (rcnt >= max(Cmaxnnz, Smaxnnz)) cout << "rcnt =" << rcnt << "Smaxnnz =" << Smaxnnz << endl;
                assert(rcnt < max(Cmaxnnz, Smaxnnz));
                RightIndex[rcnt] = A->GCID(Ai[j]);
                RightValues[rcnt++] = Ax[j];
#ifdef BLOCK_DIAGONAL_Si
                }
#endif
            }
        }

        if (gvals[gid] == 1)
        { // D or C row
            D.InsertGlobalValues(gid, lcnt, LeftValues, LeftIndex);
            //C.InsertGlobalValues(gid, rcnt, RightValues, RightIndex);
            // This is ugly code for cache TODO
            for (int k = 0; k < rcnt ; k++)
            {
#ifdef BLOCK_DIAGONAL_Si
                CMV->ReplaceGlobalValue(gid, RRowMax[RightIndex[k]], 
                                        RightValues[k]);
#else
                CMV->ReplaceGlobalValue(gid, CColMax[RightIndex[k]], 
                                        RightValues[k]);
#endif
            }
        }
        else
        { // R or S row
            //R.InsertGlobalValues(gid, lcnt, LeftValues, LeftIndex);
            assert(rcnt == SNumEntriesPerRow[debug_row]);
            debug_row++;
            for (int k = 0; k < lcnt ; k++)
            {
                RMV->ReplaceGlobalValue(LeftIndex[k], RRowMax[gid], 
                    LeftValues[k]);
            }
            err = S.InsertGlobalValues(gid, rcnt, RightValues, RightIndex);

            // This is ugly code for cache TODO
            for (int k = 0; k < rcnt ; k++)
            {
#ifdef BLOCK_DIAGONAL_Si
                SMV->ReplaceGlobalValue(gid, RRowMax[RightIndex[k]],
                                        RightValues[k]);
#else
                SMV->ReplaceGlobalValue(gid, CColMax[RightIndex[k]],
                                        RightValues[k]);
#endif
            }
            assert(err == 0);
        }
    }
    D.FillComplete();
#ifdef BLOCK_DIAGONAL_Si
    S.FillComplete();
#else
    S.FillComplete(LocalSColMap, LocalSRowMap);
#endif
    //C.FillComplete(LocalCDomainMap, LocalDRowMap);
    //R.FillComplete(LocalDDomainMap, LocalSRowMap);
    // A is no longer needed
    delete[] LeftIndex;
    delete[] LeftValues;
    delete[] RightIndex;
    delete[] RightValues;

    cout << msg << "S rows=" << S.NumGlobalRows() << " S cols=" << 
        S.NumGlobalCols() << "#cols in column map="<< 
        S.ColMap().NumMyElements() << endl;

    // ]

    // ======================= Numeric factorization =========================
    //Amesos_BaseSolver *Solver;
    Amesos Factory;
    char* SolverType = "Amesos_Klu";
    bool IsAvailable = Factory.Query(SolverType);

    cout << msg <<" #columns in CMV ="<< Cnc << endl;
    cout << msg <<" #rows in CMV ="<< Dnr << endl;

    Epetra_MultiVector *X = new Epetra_MultiVector(LocalDRowMap, Sbd_nc);
    Epetra_MultiVector *residual = new Epetra_MultiVector(LocalDRowMap, Sbd_nc);
    double *resid = new double[Sbd_nc];

    LP = new Epetra_LinearProblem();
    LP->SetOperator(&D);
    LP->SetLHS(X);
    LP->SetRHS(CMV);
    Solver = Factory.Create(SolverType, *LP);

    Solver->SymbolicFactorization();
    Solver->NumericFactorization();
    Solver->Solve();

    D.Multiply(false, *X, *residual);
    residual->Update(-1.0, *CMV, 1.0);
    residual->Norm2(resid);
    double max_resid = 0;
    for (int i = 0 ; i < Sbd_nc; i++)
        max_resid = max(max_resid, abs(resid[i]));

    cout << "Maximum residual = " << max_resid << endl;

#ifdef DUMP_MATRICES
    if (myPID == 1)
    {
        EpetraExt::MultiVectorToMatlabFile("RHS.mat", *CMV);
        EpetraExt::MultiVectorToMatlabFile("LHS.mat", *X);
        EpetraExt::RowMatrixToMatlabFile("D.mat", D);
    }
#endif

    localS = new Epetra_MultiVector(LocalCColMap, Sbd_nc);
    localS->Multiply('T', 'N', 1.0, *RMV, *X, 1.0);
    localS->Update(1.0, *SMV, -1.0);

#ifdef DUMP_MATRICES
    if (myPID == 0)
    {
        EpetraExt::MultiVectorToMatlabFile("S0.mat", *localS);
    }
#endif

    Teuchos::LAPACK<int, double> lapack;
    int info;
    piv = new int[Sbd_nc];
    lapack.GETRF(Sbd_nc, Sbd_nc, localS->Values(), Sbd_nc, piv, &info);

    // ======================= Solve =========================================

    // Clean up
    delete X;
    //delete B;
    //delete CMV;
    delete RMV;
    delete SMV;
    delete[] DNumEntriesPerRow;
    delete[] SNumEntriesPerRow;
    delete[] DColElems;
    delete[] SColElems;
    delete[] CColMax;
    delete[] RRowMax;
    delete[] CColElems;
    delete[] RRowElems;
    //delete[] DRowElems;
    //delete[] SRowElems;

    delete[] gvals;
    //delete A;
    return 0;
}
