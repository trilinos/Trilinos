
//@HEADER
// ************************************************************************
// 
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

/** \file shylu_factor.cpp

    \brief Factors sparse matrix using LU factorization.

    \author Siva Rajamanickam

    This version extracts the non zero rows/columns of R/C, hence Si is smaller
    the complete diagonal block. This version also
    extracts the entire As matrix in each proc. This is not a huge problem.
    When one needs vec = S * v = As * v - Si * v(corres rows) and minus updates
    the corresponding rows of vec correctly. This also needs storing both As and
    Si. The preconditioner for this method is not yet written.


*/

/* This define will make S block diagonal, To make this happen even some
   zero columns/rows of C and R will be stored.
   */
#define BLOCK_DIAGONAL_Si

#include "shylu.h"
#include "shylu_util.h"
#include <EpetraExt_Reindex_LinearProblem2.h>

int create_matrices
(
    Epetra_CrsMatrix *A,    // i/p: A matrix
    shylu_symbolic *ssym,   // symbolic structure
    shylu_data *data,       // numeric structure, TODO: Required ?
    shylu_config *config    // i/p: library configuration
)
{
    int Dnr = data->Dnr;
    int Snr = data->Snr;
    int Dnc = data->Dnc;
    int *DRowElems = data->DRowElems;
    int *SRowElems = data->SRowElems;
    int *DColElems = data->DColElems;
    int *gvals = data->gvals;

    double Sdiagfactor = config->Sdiagfactor;
    int sym = config->sym;

    /* --------------- Create the maps for the DBBD form ------------------- */
    // Create the local and distributed row map
    Epetra_MpiComm LComm(MPI_COMM_SELF);

	// Use Serial Comm for the local blocks.
    Epetra_Map LocalDRowMap(-1, Dnr, DRowElems, 0, LComm);
    Epetra_Map DRowMap(-1, Dnr, DRowElems, 0, A->Comm());
    Epetra_Map SRowMap(-1, Snr, SRowElems, 0, A->Comm());
    Epetra_Map LocalSRowMap(-1, Snr, SRowElems, 0, LComm);

    // Create the local column map
    Epetra_Map LocalDColMap(-1, Dnc, DColElems, 0, LComm);

    /*----------- Compute nnz per row for all five matrices ----------------- */
    // Compute nnz per row.
    int *Ai;
    double *Ax;
    int *DNumEntriesPerRow = new int[Dnr];
    int *GNumEntriesPerRow = new int[Snr];
    int *CNumEntriesPerRow = new int[Dnr];
    int *RNumEntriesPerRow = new int[Snr];
    int *SBarNumEntriesPerRow = new int[Snr];
    int Dmaxnnz=0, Cmaxnnz=0, Rmaxnnz=0, Smaxnnz=0; // Max nnz in any row

    // Find the required no of diagonals
    /*int Sdiag = (int) SNumGlobalCols * Sdiagfactor;
    //cout << "No of diagonals in Sbar =" << Sdiag << endl;
    Sdiag = MIN(Sdiag, SNumGlobalCols-1);*/
    int Sdiag = (int) Snr * Sdiagfactor;
    Sdiag = MIN(Sdiag, Snr-1);
    Sdiag = MAX(Sdiag, 0);
    //cout << "No of diagonals in Sbar =" << Sdiag << endl;
    //assert (Sdiag <= SNumGlobalCols-1);
    if (Snr != 0) assert (Sdiag <= Snr-1);

    int dcol, ccol, rcol, scol;
    int gcid;
    int gid;

    int nrows = A->RowMap().NumMyElements();
    int *rows = A->RowMap().MyGlobalElements();
    if (sym)
    {

        int NumEntries;
        int dcnt, scnt, ccnt, rcnt;
        dcol = 0 ; ccol = 0; rcol = 0; scol = 0;
        // Three things to worry about here, Num entries in each row of S,
        // size of C and R and the corresponding cols/rows.
        for (int i = 0; i < nrows; i++) 
        {
            dcnt = 0; scnt = 0; ccnt = 0; rcnt = 0;
            // Need to pass local id i to this function 
            int err = A->ExtractMyRowView(i, NumEntries, Ax, Ai);
            if (err != 0)
            { 
                assert (err == 0);
                //config->dm.error("create_matrices: extract_my_row failed");
            }

            gid = rows[i];
            if (gvals[gid] == 1)
            { // D or C row
                for (int j = 0 ; j < NumEntries ; j++)
                { // O(nnz) ! Careful what you do inside
                    gcid = A->GCID(Ai[j]);
                    // Only cols corresponding to local rows belong to D.
                    if (A->LRID(gcid) != -1 && gvals[gcid] == 1)
                    {
                        dcnt++;
                    }
                    else
                    { 
                        ccnt++;
                    }
                }
                // There should not be a null row in D.
                assert(dcnt != 0);
                DNumEntriesPerRow[dcol++] = dcnt;
                CNumEntriesPerRow[ccol++] = ccnt;
            }
            else
            { // R or S row
                //cout << "proc/row " << myPID << "/" << gid;
                //cout << " has Cols = " ;
                for (int j = 0 ; j < NumEntries ; j++)
                { // O(nnz) ! Careful what you do inside
                    gcid = A->GCID(Ai[j]); 
                    //if (A->LRID(gcid) != -1 && gvals[gcid] == 1) // TBD : Need to change here
                    // No need to check for local rows as above.
                    // All cols with gvals[cols] == 1 belong to R.
                    if (gvals[gcid] == 1)
                    {
                        rcnt++;
                    }
                    else
                    {
                        scnt++;
                    }
                }
                // There should not be a null row in S.
                assert(scnt != 0);
                assert(scol < Snr);
                assert(rcol < Snr);
                GNumEntriesPerRow[scol++] = scnt;
                RNumEntriesPerRow[rcol++] = rcnt;
            }
            Dmaxnnz = max(Dmaxnnz, dcnt);
            Smaxnnz = max(Smaxnnz, scnt); 
            Rmaxnnz = max(Rmaxnnz, rcnt); 
            Cmaxnnz = max(Cmaxnnz, ccnt); 
        }
        assert( dcol == Dnr);
        assert( scol == Snr);
        assert( ccol == Dnr);
        assert( rcol == Snr);
        int sbarcol = 0;
        for (int i = 0; i < nrows; i++) 
        {
            gid = rows[i];
            if (gvals[gid] != 1)
            { // S row
                // Add the number of required diagonals in approximate Schur 
                // complement , +1 below for the main diagonal
                assert(sbarcol < Snr);
                SBarNumEntriesPerRow[sbarcol] = GNumEntriesPerRow[sbarcol] 
                                                + Sdiag*2 + 1;
                sbarcol++;
                //SBarNumEntriesPerRow[i] = GNumEntriesPerRow[i] ;

            }
        }
        assert( sbarcol == Snr);
        Smaxnnz = Smaxnnz + Sdiag * 2 + 1;
    }
    else
    {
        assert (0 == 1);
        /*assert(Dnr == nrows);
        // TODO : Need to compute the entries for D and C in one array of size
        // nrow*/
    }


    //Create the local matrices
    ssym->D = Teuchos::RCP<Epetra_CrsMatrix> (new Epetra_CrsMatrix(Copy,
                        LocalDRowMap, LocalDColMap, DNumEntriesPerRow, true));
    //config->dm.print(5, "Created D matrix");

    // Leave the column map out, Let Epetra do the work in the rest of the cases
    ssym->G = Teuchos::RCP<Epetra_CrsMatrix> (new Epetra_CrsMatrix(Copy,
                                SRowMap, GNumEntriesPerRow, true));
    //config->dm.print(5, "Created S matrix");
    if (config->sep_type != 1)
    {

        ssym->C = Teuchos::RCP<Epetra_CrsMatrix> (new Epetra_CrsMatrix(Copy,
                                DRowMap, CNumEntriesPerRow, true));
        //config->dm.print(5, "Created C matrix");

        ssym->R = Teuchos::RCP<Epetra_CrsMatrix> (new Epetra_CrsMatrix(Copy,
                                SRowMap, RNumEntriesPerRow, true));
        //config->dm.print(5, "Created R matrix");
    }
    else
    {
        ssym->C = Teuchos::RCP<Epetra_CrsMatrix> (new Epetra_CrsMatrix(Copy,
                                    LocalDRowMap, CNumEntriesPerRow, true));
        //config->dm.print(5, "Created C matrix");

        ssym->R = Teuchos::RCP<Epetra_CrsMatrix> (new Epetra_CrsMatrix(Copy,
                                    LocalSRowMap, RNumEntriesPerRow, true));
        //config->dm.print(5, "Created R matrix");
    }

    if (config->schurApproxMethod == 1)
    {
        // TODO: Check for the local case.
        ssym->Sg = Teuchos::RCP<Epetra_CrsGraph> (new Epetra_CrsGraph(Copy,
                                SRowMap, SBarNumEntriesPerRow, false));
        //config->dm.print(5, "Created Sg graph");
    }

    data->lmax = max(Dmaxnnz, Rmaxnnz);
    data->rmax = max(Cmaxnnz, Smaxnnz);

    delete[] DNumEntriesPerRow;
    delete[] GNumEntriesPerRow;
    delete[] CNumEntriesPerRow;
    delete[] RNumEntriesPerRow;
    delete[] SBarNumEntriesPerRow;
    return 0;
}

int extract_matrices
(
    Epetra_CrsMatrix *A,    // i/p: A matrix
    shylu_symbolic *ssym,   // symbolic structure
    shylu_data *data,       // numeric structure, TODO: Required ?
    shylu_config *config,   // i/p: library configuration
    bool insertValues       // true implies values will be inserted and fill
                            // complete will be called. false implies values
                            // will be replaced.
)
{
    Teuchos::RCP<Epetra_CrsMatrix> D = ssym->D;
    Teuchos::RCP<Epetra_CrsMatrix> C = ssym->C;
    Teuchos::RCP<Epetra_CrsMatrix> R = ssym->R;
    Teuchos::RCP<Epetra_CrsMatrix> G = ssym->G;
    Teuchos::RCP<Epetra_CrsGraph> Sg = ssym->Sg;
    int *DColElems = data->DColElems;
    int *gvals = data->gvals;
    double Sdiagfactor = config->Sdiagfactor;

    int *LeftIndex = new int[data->lmax];
    double *LeftValues = new double[data->lmax];
    int *RightIndex = new int[data->rmax];
    double *RightValues = new double[data->rmax];
    int err;
    int lcnt, rcnt ;
    int gcid;
    int gid;
    int *Ai;
    double *Ax;

    int nrows = A->RowMap().NumMyElements();
    int *rows = A->RowMap().MyGlobalElements();

    for (int i = 0; i < nrows ; i++)
    {
        int NumEntries;
        err = A->ExtractMyRowView(i, NumEntries, Ax, Ai);

        lcnt = 0; rcnt = 0;
        // Place the entry in the correct sub matrix, Works only for sym
        gid = rows[i];
        int lcid;
        for (int j = 0 ; j < NumEntries ; j++)
        { // O(nnz) ! Careful what you do inside
            // Row permutation does not matter here 
            gcid = A->GCID(Ai[j]);
            assert(gcid != -1);
            //Either in D or R 
            if ((gvals[gid] != 1 && gvals[gcid] == 1)
               || (gvals[gid] == 1 && A->LRID(gcid) != -1 && gvals[gcid] == 1))
            {
                assert(lcnt < data->lmax);
                if (insertValues)
                    LeftIndex[lcnt] = gcid;
                else
                {
                    //local column id
                    lcid = (gvals[gid] == 1 ? D->LCID(gcid) : R->LCID(gcid));
                    assert(lcid != -1);
                    LeftIndex[lcnt] = lcid;
                }
                LeftValues[lcnt++] = Ax[j];
            }
            else
            {
                assert(rcnt < data->rmax);
                if (insertValues)
                    RightIndex[rcnt] = gcid;
                else
                {
                    //local column id
                    lcid = (gvals[gid] == 1 ? C->LCID(gcid) : G->LCID(gcid));
                    assert(lcid != -1);
                    RightIndex[rcnt] = lcid;
                }
                RightValues[rcnt++] = Ax[j];
            }
        }

        if (gvals[gid] == 1)
        { // D or C row
            if (insertValues)
            {
                err = D->InsertGlobalValues(gid, lcnt, LeftValues, LeftIndex);
                assert(err == 0);
                err = C->InsertGlobalValues(gid, rcnt, RightValues, RightIndex);
                assert(err == 0);
            }
            else
            {
                err = D->ReplaceMyValues(D->LRID(gid), lcnt, LeftValues,
                                    LeftIndex);
                assert(err == 0);
                err = C->ReplaceMyValues(C->LRID(gid), rcnt, RightValues,
                                    RightIndex);
                assert(err == 0);
            }
        }
        else
        { // R or S row
            //assert(lcnt > 0); // TODO: Enable this once using narrow sep.
            if (insertValues)
            {
                assert(rcnt > 0);
                err = R->InsertGlobalValues(gid, lcnt, LeftValues, LeftIndex);
                assert(err == 0);
                err = G->InsertGlobalValues(gid, rcnt, RightValues, RightIndex);
                assert(err == 0);
                if (config->schurApproxMethod == 1)
                {
                    err = Sg->InsertGlobalIndices(gid, rcnt, RightIndex);
                    assert(err == 0);
                }
            }
            else
            {
                assert(rcnt > 0);
                err = R->ReplaceMyValues(R->LRID(gid), lcnt, LeftValues,
                                    LeftIndex);
                assert(err == 0);
                err = G->ReplaceMyValues(G->LRID(gid), rcnt, RightValues,
                                    RightIndex);
                assert(err == 0);
            }
        }
    }

    if (insertValues)
    {
        /* ------------- Create the maps for the DBBD form ------------------ */
        Epetra_Map *DRowMap, *SRowMap, *DColMap;
        Epetra_SerialComm LComm;
        if (config->sep_type != 1)
        {
            DRowMap = new Epetra_Map(-1, data->Dnr, data->DRowElems, 0,
                             A->Comm());
            SRowMap = new Epetra_Map(-1, data->Snr, data->SRowElems, 0,
                             A->Comm());
            DColMap = new Epetra_Map(-1, data->Dnc, DColElems, 0,
                             A->Comm());
        }
        else
        {
            DRowMap = new Epetra_Map(-1, data->Dnr, data->DRowElems, 0, LComm);
            SRowMap = new Epetra_Map(-1, data->Snr, data->SRowElems, 0, LComm);
            DColMap = new Epetra_Map(-1, data->Dnc, DColElems, 0, LComm);
        }

        D->FillComplete();
        //config->dm.print(5, "Done D fillcomplete");

        G->FillComplete();
        //config->dm.print(5, "Done G fillcomplete");

        C->FillComplete(*SRowMap, *DRowMap); //TODO:Won't work if permutation is
                                                // unsymmetric SRowMap
        //config->dm.print(5, "Done C fillcomplete");

        R->FillComplete(*DColMap, *SRowMap);
        //config->dm.print(5, "Done R fillcomplete");

        int Sdiag = (int) data->Snr * Sdiagfactor;
        Sdiag = MIN(Sdiag, data->Snr-1);
        Sdiag = MAX(Sdiag, 0);

        // Add the diagonals to Sg
        for (int i = 0; config->schurApproxMethod == 1 && i < nrows ; i++)
        {
            gid = rows[i];
            if (gvals[gid] == 1) continue; // not a row in S
            if (data->Snr == 0) assert(0 == 1);

            rcnt = 0;
            //TODO Will be trouble if SNumGlobalCols != Snc
            //assert(SNumGlobalCols == Snc);
            //for (int j = MAX(i-Sdiag,0) ; j<MIN(SNumGlobalCols, i+Sdiag); j++)
            for (int j = MAX(i-Sdiag, 0) ; j < MIN(data->Snr, i+Sdiag); j++)
            {
                // find the adjacent columns from the row map of S
                //assert (j >= 0 && j < Snr);
                RightIndex[rcnt++] = data->SRowElems[j];
            }
            err = Sg->InsertGlobalIndices(gid, rcnt, RightIndex);
            assert(err == 0);
            // Always insert the diagonals, if it is added twice that is fine.
            err = Sg->InsertGlobalIndices(gid, 1, &gid);
            assert(err == 0);
        }

        if (config->schurApproxMethod == 1)
            Sg->FillComplete();

        delete DRowMap;
        delete SRowMap;
        delete DColMap;
    }

#if 0
    if (insertValues)
    {
#ifdef TIMING_OUTPUT
    Teuchos::Time ttime("transpose time");
    ttime.start();
#endif
        bool MakeDataContiguous = true;
        ssym->transposer = Teuchos::RCP<EpetraExt::RowMatrix_Transpose>(new EpetraExt::RowMatrix_Transpose(MakeDataContiguous));
        ssym->DT = Teuchos::rcp( dynamic_cast<Epetra_CrsMatrix *>(&(*ssym->transposer)(*D)), false);

#ifdef TIMING_OUTPUT
    ttime.stop();
    cout << "Transpose Time" << ttime.totalElapsedTime() << endl;
    ttime.reset();
#endif

    }
    else
    {
        ssym->transposer->fwd();
        //ssym->ReIdx_LP->fwd(); // TODO: Needed ?
    }
#endif

    // A is no longer needed
    delete[] LeftIndex;
    delete[] LeftValues;
    delete[] RightIndex;
    delete[] RightValues;

    //cout << msg << "S rows=" << S.NumGlobalRows() << " S cols=" <<
        //S.NumGlobalCols() << "#cols in column map="<<
        //S.ColMap().NumMyElements() << endl;
    //cout << msg << "C rows=" << Cptr->NumGlobalRows() << " C cols=" <<
        //Cptr->NumGlobalCols() << "#cols in column map="<<
        //Cptr->ColMap().NumMyElements() << endl;
    //cout << msg << "D rows=" << D.NumGlobalRows() << " D cols=" <<
        //D.NumGlobalCols() << "#cols in column map="<<
        //D.ColMap().NumMyElements() << endl;
    //cout << msg << "R rows=" << Rptr->NumGlobalRows() << " R cols=" <<
        //Rptr->NumGlobalCols() << "#cols in column map="<<
        //Rptr->ColMap().NumMyElements() << endl;
    // ]

    return 0;

}

/* Find the DBBD form */
int shylu_symbolic_factor
(
    Epetra_CrsMatrix *A,    // i/p: A matrix
    shylu_symbolic *ssym,   // symbolic structure
    shylu_data *data,       // numeric structure, TODO: Required ?
    shylu_config *config    // i/p: library configuration
)
{
#ifdef TIMING_OUTPUT
    Teuchos::Time symtime("symbolic time");
    symtime.start();
#endif
    int myPID = A->Comm().MyPID();
    int n = A->NumGlobalRows();

    int Dnr;
    int Snr;
    int *DRowElems;
    int *SRowElems;
    int sym = config->sym;

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
    int SNumGlobalCols;
    findLocalColumns(A, gvals, SNumGlobalCols);

    // See if you can shrink the separator by assigning more rows/columns to
    // the block diagonals
    // TODO: This is because of a bug in coloring remove the if once that is
    // fixed
    //if (config->schurApproxMethod == 2)
    if (config->sep_type == 2)
        findNarrowSeparator(A, gvals);

    // 3. Assemble diagonal block and the border in convenient form [
    /* In each processor, we have (in a permuted form)
     *  | D_i    C_i   |
     *  | R_i    S_i   |
     * D_i - diagonal block, C_i - Column Separator, R_i - Row separator
     * S_i - A22 block corresponding to Schur complement part of A
     * Assemble all four blocks in local matrices. */

     ostringstream ssmsg1;
     ssmsg1 << "PID =" << myPID << " ";
     string msg = ssmsg1.str();
     ssmsg1.clear(); ssmsg1.str("");

    // Find #cols in each block
    int Dnc = 0;        // #cols in diagonal block
    int Snc = 0;        // #cols in the col. separator
    /* Looping on cols will work only for wide separator
     * as for narrow sep there will be some sep cols with gvals[col] ==1
     * */
    /*for (int i = 0; i < ncols ; i++)
    {
        if (gvals[cols[i]] == 1)
            Dnc++;
        else
            Snc++;
    }
    // Find #rows in each block 
    Dnr = Dnc;          // #rows in square diagonal block
    Snr = nrows - Dnr;  // #rows in the row separator*/

    // Find #rows in each block
    Dnr = 0;
    Snr = 0;
    for (int i = 0; i < nrows ; i++)
    {
        if (gvals[rows[i]] == 1)
            Dnr++;
        else
            Snr++;
    }
    Dnc = Dnr;
    // TODO: Snc is no longer useful, should remove it
    for (int i = 0; i < ncols ; i++)
    {
        if (gvals[cols[i]] != 1)
            Snc++;
    }

    assert(Snc >= 0);

    // TODO : The above assignment may not be correct in the unsymetric case

    ////config->dm.print(2, msg + " Mycols=");
    cout << msg << " Mycols="<< ncols << "Myrows ="<< nrows << endl;
    cout << msg << " #rows and #cols in diagonal blk ="<< Dnr << endl;
    cout << msg << " #columns in S ="<< Snc << endl;
    cout << msg << " #rows in S ="<< Snr << endl;

    ostringstream pidstr;
    pidstr <<  myPID ;
    // Create a row map for the D and S blocks [
    DRowElems = new int[Dnr];
    SRowElems = new int[Snr];
    int gid;
    // Assemble row ids in two arrays (for D and R blocks)
    if (sym)
    {
        findBlockElems(A, nrows, rows, gvals, Dnr, DRowElems, Snr, SRowElems,
                    "D"+pidstr.str()+"Rows", "S"+pidstr.str()+"Rows", false) ;
    }
    else
    {
        // SRowElems are not known until factorization, TODO
        assert(0 == 1);
    }

    data->Dnr = Dnr;
    data->Snr = Snr;
    data->Dnc = Dnc;
    data->DRowElems = DRowElems;
    data->SRowElems = SRowElems;

    // Create a column map for the D and S blocks [
    int *DColElems = new int[Dnc]; // Elems in column map of D 
    int *SColElems = new int[Snc]; // Elems in column map of C TODO: Unused
    // Assemble column ids in two arrays (for D and C blocks)
    findBlockElems(A, ncols, cols, gvals, Dnc, DColElems, Snc, SColElems,
                    "D"+pidstr.str()+"Cols", "S"+pidstr.str()+"Cols", true) ;

    data->DColElems = DColElems;
    data->gvals = gvals;

    for (int i = 0; i < Snr; i++)
    {
        // Epetra guarentees columns corresponding to local rows will be first
        // in the column map.
        assert(SRowElems[i] == SColElems[i]);
    }
    // ]

    /*--Create the Epetra Matrices with the maps (does not insert values) --- */
    create_matrices(A, ssym, data, config);

    /*--Extract the Epetra Matrices and call fillComplete --- */
    extract_matrices(A, ssym, data, config, true);

    delete[] SColElems;

    Amesos Factory;
    const char* SolverType = config->diagonalBlockSolver.c_str();
    bool IsAvailable = Factory.Query(SolverType);
    assert(IsAvailable == true);

    Teuchos::RCP<Epetra_LinearProblem> LP = Teuchos::RCP<Epetra_LinearProblem> 
                                        (new Epetra_LinearProblem());
    LP->SetOperator((ssym->D).getRawPtr());
    //LP->SetOperator((ssym->DT).getRawPtr()); // for transpose

    // Create temp vectors
    ssym->Dlhs = Teuchos::RCP<Epetra_MultiVector>
                    (new Epetra_MultiVector(ssym->D->RowMap(), 16));
    ssym->Drhs = Teuchos::RCP<Epetra_MultiVector>
                    (new Epetra_MultiVector(ssym->D->RowMap(), 16));
    ssym->Gvec = Teuchos::RCP<Epetra_MultiVector>
                    (new Epetra_MultiVector(ssym->G->RowMap(), 16));

    LP->SetRHS(ssym->Drhs.getRawPtr());
    LP->SetLHS(ssym->Dlhs.getRawPtr());

    ssym->ReIdx_LP = Teuchos::RCP<
                    EpetraExt::ViewTransform<Epetra_LinearProblem> >
                    (new EpetraExt::LinearProblem_Reindex2(0));
    ssym->LP = Teuchos::RCP<Epetra_LinearProblem>(&((*(ssym->ReIdx_LP))(*LP)),
                                        false);

    Teuchos::RCP<Amesos_BaseSolver> Solver = Teuchos::RCP<Amesos_BaseSolver>
                                    (Factory.Create(SolverType, *(ssym->LP)));
    //config->dm.print(5, "Created the diagonal solver");

#ifdef TIMING_OUTPUT
    Teuchos::Time ftime("setup time");
    ftime.start();
#endif
    //Solver->SetUseTranspose(true); // for transpose
    Teuchos::ParameterList aList;
    aList.set("TrustMe", true);
    Solver->SetParameters(aList);
    Solver->SymbolicFactorization();
    //config->dm.print(3, "Symbolic Factorization done");

#ifdef TIMING_OUTPUT
    ftime.stop();
    cout << "Symbolic Factorization Time" << ftime.totalElapsedTime() << endl;
    ftime.reset();
#endif

    ssym->OrigLP = LP;
    //ssym->LP = LP;
    ssym->Solver = Solver;

    if (config->schurApproxMethod == 1)
    {
        Teuchos::ParameterList pList;
        Teuchos::RCP<Isorropia::Epetra::Prober> prober = 
                         Teuchos::RCP<Isorropia::Epetra::Prober> (new
                          Isorropia::Epetra::Prober((ssym->Sg).getRawPtr(),
                                                     pList, false));
        //config->dm.print(3, "Doing Coloring");
#ifdef TIMING_OUTPUT
        ftime.start();
#endif
        prober->color();
#ifdef TIMING_OUTPUT
        ftime.stop();
        cout << "Time to color" << ftime.totalElapsedTime() << endl;
        ftime.reset();
        ftime.start();
#endif
        ssym->prober = prober;
    }
#ifdef TIMING_OUTPUT
    symtime.stop();
    cout << "Symbolic Time" << symtime.totalElapsedTime() << endl;
    symtime.reset();
#endif
}

int shylu_factor(Epetra_CrsMatrix *A, shylu_symbolic *ssym, shylu_data *data,
             shylu_config *config)
{
#ifdef TIMING_OUTPUT
    Teuchos::Time fact_time("factor time");
    fact_time.start();
#endif

    Teuchos::RCP<Epetra_LinearProblem> LP = ssym->LP;
    Teuchos::RCP<Amesos_BaseSolver> Solver = ssym->Solver;
    Teuchos::RCP<Isorropia::Epetra::Prober> prober = ssym->prober;

    /*--Extract the Epetra Matrices into already existing matrices --- */
    extract_matrices(A, ssym, data, config, false);

    // ======================= Numeric factorization =========================
#ifdef TIMING_OUTPUT
    Teuchos::Time ftime("setup time");
    ftime.start();
#endif

    //config->dm.print(3, "In Numeric Factorization");
    Solver->NumericFactorization();
    //config->dm.print(3, "Numeric Factorization done");

#ifdef SHYLU_DEBUG
    Solver->PrintStatus();
#endif

#ifdef TIMING_OUTPUT
    ftime.stop();
    cout << "Time to factor" << ftime.totalElapsedTime() << endl;
    ftime.reset();
#endif

    // Create the local and distributed row map
    Epetra_MpiComm LComm(MPI_COMM_SELF);
    Epetra_Map LocalDRowMap(-1, data->Dnr, data->DRowElems, 0, LComm);
    Teuchos::RCP<Epetra_CrsMatrix> Sbar;

    data->schur_op = Teuchos::RCP<ShyLU_Probing_Operator> (new
             ShyLU_Probing_Operator(ssym, (ssym->G).getRawPtr(),
             (ssym->R).getRawPtr(), (ssym->LP).getRawPtr(),
             (ssym->Solver).getRawPtr(), (ssym->C).getRawPtr(), &LocalDRowMap,
             1));

    if (config->schurApproxMethod == 1)
    {
        int nvectors = prober->getNumOrthogonalVectors();
        //Set up the probing operator
        // TODO: Change to RCPs. Call Set vectors on schur_op and remove
        // probeop
        ShyLU_Probing_Operator probeop(ssym, (ssym->G).getRawPtr(),
         (ssym->R).getRawPtr(), (ssym->LP).getRawPtr(),
         (ssym->Solver).getRawPtr(), (ssym->C).getRawPtr(), &LocalDRowMap,
         nvectors);

        //cout << "Doing probing" << endl;
        Sbar = prober->probe(probeop);
        //cout << "SIZE of SBAR = " << (*Sbar).NumGlobalRows() << endl;
#ifdef TIMING_OUTPUT
        ftime.stop();
        cout << "Time to probe" << ftime.totalElapsedTime() << endl;
        ftime.reset();
#endif
    }
    else if (config->schurApproxMethod == 2)
    {
        // Compute and drop the entries in the Schur complement
        // Ignore the structure of the Schur complement
        //cout << "Computing the Approx Schur complement" << endl;
#ifdef TIMING_OUTPUT
        ftime.start();
#endif

        // TODO: Pass the operator rather than all the matrices
        if (config->sep_type == 2)
            Sbar = computeApproxSchur(config, ssym, (ssym->G).getRawPtr(),
             (ssym->R).getRawPtr(), (ssym->LP).getRawPtr(),
             (ssym->Solver).getRawPtr(), (ssym->C).getRawPtr(),
             &LocalDRowMap);
        else
            Sbar = computeApproxWideSchur(config, ssym, (ssym->G).getRawPtr(),
             (ssym->R).getRawPtr(), (ssym->LP).getRawPtr(),
             (ssym->Solver).getRawPtr(), (ssym->C).getRawPtr(),
             &LocalDRowMap);

        //cout << *Sbar << endl;
#ifdef TIMING_OUTPUT
        ftime.stop();
        cout << "Time to Compute Approx Schur Complement" << ftime.totalElapsedTime() << endl;
        ftime.reset();
#endif
        //cout << "Computed Approx Schur complement" << endl;
    }
    else if (config->schurApproxMethod == 3)
    {
        Sbar = computeSchur_GuidedProbing(config, ssym, data, &LocalDRowMap);

        // TODO: The above call does not work with narrow separator. We need
        // distribute probing here, save the Sbar from the iteration in
        // data in a different variable other than Sbar (as Sbar gets
        // overwritten after every iteration, but we need the original Sbar
        // graph for the probing in later iterations) Set the prober, when
        // new dropped is computed in previous iteration, otherwise just do
        // probing, same idea from above function this time it is distributed.
        // This was slower when I implemented it for wide sep (and I deleted
        // it !).
    }


    data->Sbar = Sbar;

    if (config->schurSolver == "Amesos")
    {
        Teuchos::ParameterList aList;
        //aList.set("Reindex", true);

        data->Sbarlhs = Teuchos::RCP<Epetra_MultiVector>
                        (new Epetra_MultiVector(data->Sbar->RowMap(), 16));
        data->Sbarrhs = Teuchos::RCP<Epetra_MultiVector>
                        (new Epetra_MultiVector(data->Sbar->RowMap(), 16));

        Teuchos::RCP<Epetra_LinearProblem> LP2 =
        		Teuchos::rcp(new Epetra_LinearProblem);

        LP2->SetOperator(data->Sbar.get());
        LP2->SetRHS(data->Sbarrhs.getRawPtr());
        LP2->SetLHS(data->Sbarlhs.getRawPtr());

        data->ReIdx_LP2 = Teuchos::RCP
        		<EpetraExt::ViewTransform<Epetra_LinearProblem> >
				(new EpetraExt::LinearProblem_Reindex2(0));
        data->LP2 =
        		Teuchos::RCP<Epetra_LinearProblem>
        		(&((*(data->ReIdx_LP2))(*LP2)), false);

        data->OrigLP2 = LP2;

        Amesos Factory;
        bool IsAvailable = Factory.Query(config->schurAmesosSolver);
        assert(IsAvailable == true);
        data->dsolver = Factory.Create(
        		config->schurAmesosSolver, *(data->LP2));

        data->dsolver->SetParameters(aList);
        //cout << "Created the direct Schur  Solver" << endl;

        data->dsolver->SymbolicFactorization();
        //cout << "Symbolic Factorization done for schur complement" << endl;

        //cout << "In Numeric Factorization of Schur complement" << endl;
        data->dsolver->NumericFactorization();
        //cout << "Numeric Factorization done for schur complement" << endl;

    }
    else if (config->schurSolver == "AztecOO-Exact")
    {
        data->schur_prec = Teuchos::RCP<AmesosSchurOperator> (new
                                     AmesosSchurOperator(Sbar.getRawPtr()));
        data->schur_prec->Initialize();
        data->schur_prec->Compute();

        int err = data->innersolver->SetUserOperator
                    (data->schur_op.get());
        assert (err == 0);

        err = data->innersolver->SetPrecOperator
                    (data->schur_prec.get());
        assert (err == 0);

        data->innersolver->SetAztecOption(AZ_solver, AZ_gmres);
        data->innersolver->SetMatrixName(999);
    }
    else if (config->schurSolver == "AztecOO-Inexact")
    {
        // Doing an inexact Schur complement
        if (config->libName == "Belos")
        {
            int err = data->innersolver->SetUserMatrix
                        (data->Sbar.get());
            assert (err == 0);

            data->innersolver->SetAztecOption(AZ_solver, AZ_gmres);
            data->innersolver->SetAztecOption(AZ_precond, AZ_dom_decomp);
            data->innersolver->SetAztecOption(AZ_keep_info, 1);
#ifndef DEBUG
            data->innersolver->SetAztecOption(AZ_output, AZ_none);
#endif
            //solver->SetAztecOption(AZ_overlap, 3);
            //solver->SetAztecOption(AZ_subdomain_solve, AZ_ilu);
            data->innersolver->SetMatrixName(999);

            double condest;
            err = data->innersolver->ConstructPreconditioner(condest);
            assert (err == 0);
            //cout << "Condition number of inner Sbar" << condest << endl;
        }
        else
        {
            // The outer solver is Aztec
            // I suspect there is a bug in AztecOO. Doing what we do in the if
            // case
            // here will cause an error when we use the solver in ApplyInverse
            // The error will not happen when we call the dummy JustTryIt()
            // below
            data->innersolver = NULL;
        }
    }
    else
    {
        assert (0 == 1);
    }

    //cout << " Out of factor" << endl ;
#ifdef TIMING_OUTPUT
    fact_time.stop();
    cout << "Factor Time" << fact_time.totalElapsedTime() << endl;
    fact_time.reset();
#endif
    return 0;
}
