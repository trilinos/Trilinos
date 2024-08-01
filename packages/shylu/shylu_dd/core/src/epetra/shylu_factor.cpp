// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
#ifdef HAVE_SHYLU_DDCORE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <EpetraExt_Reindex_LinearProblem2.h>

#include "Ifpack_config.h"

#ifdef HAVE_IFPACK_DYNAMIC_FACTORY
#include "Ifpack_DynamicFactory.h"
#else
#include "Ifpack.h"
#endif

#include "shylu_internal_gmres.h"
#include "shylu_internal_gmres_tools.h"

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
#ifdef HAVE_SHYLU_DDCORE_MPI
    Epetra_MpiComm LComm(MPI_COMM_SELF);
#else
    Epetra_SerialComm LComm;
#endif // HAVE_SHYLU_DDCORE_MPI

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
    //std::cout << "No of diagonals in Sbar =" << Sdiag << std::endl;
    Sdiag = SHYLU_CORE_MIN(Sdiag, SNumGlobalCols-1);*/
    int Sdiag = (int) Snr * Sdiagfactor;
    Sdiag = SHYLU_CORE_MIN(Sdiag, Snr-1);
    Sdiag = SHYLU_CORE_MAX(Sdiag, 0);
    //std::cout << "No of diagonals in Sbar =" << Sdiag << std::endl;
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
                //std::cout << "proc/row " << myPID << "/" << gid;
                //std::cout << " has Cols = " ;
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
            Dmaxnnz = std::max(Dmaxnnz, dcnt);
            Smaxnnz = std::max(Smaxnnz, scnt);
            Rmaxnnz = std::max(Rmaxnnz, rcnt);
            Cmaxnnz = std::max(Cmaxnnz, ccnt);
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

    data->lmax = std::max(Dmaxnnz, Rmaxnnz);
    data->rmax = std::max(Cmaxnnz, Smaxnnz);

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
    int err = 0;
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
        Sdiag = SHYLU_CORE_MIN(Sdiag, data->Snr-1);
        Sdiag = SHYLU_CORE_MAX(Sdiag, 0);

        // Add the diagonals to Sg
        for (int i = 0; config->schurApproxMethod == 1 && i < nrows ; i++)
        {
            gid = rows[i];
            if (gvals[gid] == 1) continue; // not a row in S
            if (data->Snr == 0) assert(0 == 1);

            rcnt = 0;
            //TODO Will be trouble if SNumGlobalCols != Snc
            //assert(SNumGlobalCols == Snc);
            //for (int j = SHYLU_CORE_MAX(i-Sdiag,0) ; j<SHYLU_CORE_MIN(SNumGlobalCols, i+Sdiag); j++)
            for (int j = SHYLU_CORE_MAX(i-Sdiag, 0) ; j < SHYLU_CORE_MIN(data->Snr, i+Sdiag); j++)
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

    // A is no longer needed
    delete[] LeftIndex;
    delete[] LeftValues;
    delete[] RightIndex;
    delete[] RightValues;

    //std::cout << msg << "S rows=" << S.NumGlobalRows() << " S cols=" <<
        //S.NumGlobalCols() << "#cols in column map="<<
        //S.ColMap().NumMyElements() << std::endl;
    //std::cout << msg << "C rows=" << Cptr->NumGlobalRows() << " C cols=" <<
        //Cptr->NumGlobalCols() << "#cols in column map="<<
        //Cptr->ColMap().NumMyElements() << std::endl;
    //std::cout << msg << "D rows=" << D.NumGlobalRows() << " D cols=" <<
        //D.NumGlobalCols() << "#cols in column map="<<
        //D.ColMap().NumMyElements() << std::endl;
    //std::cout << msg << "R rows=" << Rptr->NumGlobalRows() << " R cols=" <<
        //Rptr->NumGlobalCols() << "#cols in column map="<<
        //Rptr->ColMap().NumMyElements() << std::endl;
    // ]

    return err;
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

     std::ostringstream ssmsg1;
     ssmsg1 << "PID =" << myPID << " ";
     std::string msg = ssmsg1.str();
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
    std::cout << msg << " Mycols="<< ncols << "Myrows ="<< nrows << std::endl;
    std::cout << msg << " #rows and #cols in diagonal blk ="<< Dnr << std::endl;
    std::cout << msg << " #columns in S ="<< Snc << std::endl;
    std::cout << msg << " #rows in S ="<< Snr << std::endl;

    std::ostringstream pidstr;
    pidstr <<  myPID ;
    // Create a row map for the D and S blocks [
    DRowElems = new int[Dnr];
    SRowElems = new int[Snr];
    // int gid; // unused
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

    // Are these even needed, plan to remove them ?
    data->temp3 = Teuchos::RCP<Epetra_MultiVector>
                    (new Epetra_MultiVector(View, *(ssym->Drhs), 0, 1));
    data->locallhs = Teuchos::RCP<Epetra_MultiVector>
                    (new Epetra_MultiVector(View, *(ssym->Dlhs), 0, 1));

    LP->SetRHS(ssym->Drhs.getRawPtr());
    LP->SetLHS(ssym->Dlhs.getRawPtr());

    ssym->ReIdx_LP = Teuchos::RCP<
                    EpetraExt::ViewTransform<Epetra_LinearProblem> >
                    (new EpetraExt::LinearProblem_Reindex2(0));
    ssym->LP = Teuchos::RCP<Epetra_LinearProblem>(&((*(ssym->ReIdx_LP))(*LP)),
                                        false);

    Teuchos::RCP<Amesos_BaseSolver> Solver;
    Teuchos::RCP<Ifpack_Preconditioner> ifpackSolver;
    std::size_t found = (config->diagonalBlockSolver).find("Amesos");
    if (found == 0)
    {
        Amesos Factory;
        const char* SolverType = config->diagonalBlockSolver.c_str();
        // mfh 25 May 2015: assert() gets defined to nothing in a
        // release build.  Thus, the return value IsAvailable won't
        // get used, which causes a build warning.
#ifdef NDEBUG
        (void) Factory.Query(SolverType);
#else
        bool IsAvailable = Factory.Query(SolverType);
        assert(IsAvailable == true);
#endif // NDEBUG
        config->amesosForDiagonal = true;
        Solver = Teuchos::RCP<Amesos_BaseSolver>
                    (Factory.Create(SolverType, *(ssym->LP)));
    }
    else
    {
        config->amesosForDiagonal = false;
#ifdef HAVE_IFPACK_DYNAMIC_FACTORY
    Ifpack_DynamicFactory IfpackFactory;
#else
    Ifpack IfpackFactory;
#endif
        ifpackSolver = Teuchos::rcp<Ifpack_Preconditioner>(IfpackFactory.Create
         (config->diagonalBlockSolver, (ssym->D).getRawPtr(), 0, false));
    }
    //config->dm.print(5, "Created the diagonal solver");

#ifdef TIMING_OUTPUT
    Teuchos::Time ftime("setup time");
    ftime.start();
#endif
    //Solver->SetUseTranspose(true); // for transpose
    if (config->amesosForDiagonal)
    {
        Teuchos::ParameterList aList;
        aList.set("TrustMe", true);
        Solver->SetParameters(aList);
        Solver->SymbolicFactorization();
    }
    else
    {
        ifpackSolver->Initialize();
    }

    //config->dm.print(3, "Symbolic Factorization done");

#ifdef TIMING_OUTPUT
    ftime.stop();
    std::cout << " Shylu_Symbolic_Factor(" << myPID << ") :: Symbolic Factorization Time : " << ftime.totalElapsedTime() << std::endl;
    ftime.reset();
#endif

    ssym->OrigLP = LP;
    //ssym->LP = LP;
    ssym->Solver = Solver;
    ssym->ifSolver = ifpackSolver;

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
        std::cout << " Shylu_Symbolic_Factor(" << myPID << ") :: Time to color " << ftime.totalElapsedTime() << std::endl;
        ftime.reset();
        ftime.start();
#endif
        ssym->prober = prober;
    }

    // Set the maps, importers and multivectors to be used in the solve once.
#ifdef HAVE_SHYLU_DDCORE_MPI
    Epetra_MpiComm LComm(MPI_COMM_SELF);
#else
    Epetra_SerialComm LComm;
#endif // HAVE_SHYLU_DDCORE_MPI
    data->LDRowMap = Teuchos::rcp(new Epetra_Map(-1, data->Dnr,
                                 data->DRowElems, 0, LComm));
    data->LGRowMap = Teuchos::rcp(new Epetra_Map(-1, data->Snr,
                                 data->SRowElems, 0, LComm));
    data->GMap = Teuchos::rcp(new Epetra_Map(-1, data->Snr,
                             data->SRowElems, 0, A->Comm()));

    // Assuming X and A will have the same rowmap. Should it be domain map ?
    data->BdImporter = Teuchos::rcp(new Epetra_Import(*(data->LDRowMap),
                     A->RowMap()));
    data->DistImporter = Teuchos::rcp(new Epetra_Import(*(data->GMap),
                         *(data->LGRowMap)));;
    // Assuming X and A will have the same rowmap. Should it be domain map ?
    data->BsImporter = Teuchos::rcp(new Epetra_Import(*(data->GMap),
                     A->RowMap()));
    data->XsImporter = Teuchos::rcp(new Epetra_Import(*(data->LGRowMap),
                         *(data->GMap)));
    // Assuming Y and A will have the same rowmap. Should it be range map ?
    data->XdExporter = Teuchos::rcp(new Epetra_Export(*(data->LDRowMap),
                 A->RowMap()));
    data->XsExporter = Teuchos::rcp(new Epetra_Export(*(data->LGRowMap),
             A->RowMap()));

    // Create multivectors for solve, TODO: Can we do with fewer
    data->localrhs = Teuchos::rcp(new Epetra_MultiVector(*(data->LDRowMap), 1));
    data->temp1 = Teuchos::rcp(new Epetra_MultiVector(*(data->LGRowMap), 1));
    data->temp2 = Teuchos::rcp(new Epetra_MultiVector(*(data->GMap), 1));
    data->Bs = Teuchos::rcp(new Epetra_MultiVector(*(data->GMap), 1));
    data->Xs = Teuchos::rcp(new Epetra_MultiVector(*(data->GMap), 1));
    data->LocalXs = Teuchos::rcp(new Epetra_MultiVector(*(data->LGRowMap), 1));
    data->Xs->PutScalar(0.0);

    //data->importExportTime = Teuchos::rcp(new Teuchos::Time("import export time"));
    //data->innerIterTime = Teuchos::rcp(new Teuchos::Time("innertIter time"));
    //data->fwdTime = Teuchos::rcp(new Teuchos::Time("reindex fwd time"));
    //data->amesosSchurTime = Teuchos::rcp(new Teuchos::Time("amesos schur time"));
#ifdef TIMING_OUTPUT
    symtime.stop();
    std::cout << " Shylu_Symbolic_Factor(" << myPID << ") :: Total Time " << symtime.totalElapsedTime() << std::endl;
    symtime.reset();
#endif

    // FIXME (mfh 25 May 2015) Uh, shouldn't this function return an int???
    return 0; // note this test does not currently return any failed error codes
}

int shylu_factor(Epetra_CrsMatrix *A, shylu_symbolic *ssym, shylu_data *data,
             shylu_config *config)
{
#ifdef TIMING_OUTPUT
    int myPID = A->Comm().MyPID();
    Teuchos::Time fact_time("factor time");
    fact_time.start();
#endif

    int err = 0;
    Teuchos::RCP<Epetra_LinearProblem> LP = ssym->LP;
    Teuchos::RCP<Amesos_BaseSolver> Solver = ssym->Solver;
    Teuchos::RCP<Ifpack_Preconditioner> ifpackSolver = ssym->ifSolver;
    Teuchos::RCP<Isorropia::Epetra::Prober> prober = ssym->prober;

    /*--Extract the Epetra Matrices into already existing matrices --- */
    extract_matrices(A, ssym, data, config, false);

    // ======================= Numeric factorization =========================
#ifdef TIMING_OUTPUT
    Teuchos::Time ftime("setup time");
    ftime.start();
#endif

    //config->dm.print(3, "In Numeric Factorization");
    if (config->amesosForDiagonal)
        Solver->NumericFactorization();
    else
        ifpackSolver->Compute();

    //config->dm.print(3, "Numeric Factorization done");

#ifdef SHYLU_DEBUG
    if (config->amesosForDiagonal)
        Solver->PrintStatus();
#endif

#ifdef TIMING_OUTPUT
    ftime.stop();
    std::cout << " Shylu_Factor(" << myPID << ") :: Time to factor " << ftime.totalElapsedTime() << std::endl;
    ftime.reset();
#endif

    // Create the local and distributed row map
#ifdef HAVE_SHYLU_DDCORE_MPI
    Epetra_MpiComm LComm(MPI_COMM_SELF);
#else
    Epetra_SerialComm LComm;
#endif // HAVE_SHYLU_DDCORE_MPI
    Epetra_Map LocalDRowMap(-1, data->Dnr, data->DRowElems, 0, LComm);
    Teuchos::RCP<Epetra_CrsMatrix> Sbar;

   data->schur_op = Teuchos::RCP<ShyLU_Probing_Operator> (new
             ShyLU_Probing_Operator(config, ssym, (ssym->G).getRawPtr(),
             (ssym->R).getRawPtr(), (ssym->LP).getRawPtr(),
             (ssym->Solver).getRawPtr(), (ssym->ifSolver).getRawPtr(),
             (ssym->C).getRawPtr(), &LocalDRowMap, 1));

    if (config->schurApproxMethod == 1)
    {
        int nvectors = prober->getNumOrthogonalVectors();
        //Set up the probing operator
        // TODO: Change to RCPs. Call Set vectors on schur_op and remove
        // probeop
        ShyLU_Probing_Operator probeop(config, ssym, (ssym->G).getRawPtr(),
         (ssym->R).getRawPtr(), (ssym->LP).getRawPtr(),
         (ssym->Solver).getRawPtr(), (ssym->ifSolver).getRawPtr(),
         (ssym->C).getRawPtr(), &LocalDRowMap, nvectors);

        //std::cout << "Doing probing" << std::endl;
        Sbar = prober->probe(probeop);
        //std::cout << "SIZE of SBAR = " << (*Sbar).NumGlobalRows() << std::endl;
#ifdef TIMING_OUTPUT
        ftime.stop();
        std::cout << " Shylu_Factor(" << myPID << ") :: Time to probe " << ftime.totalElapsedTime() << std::endl;
        ftime.reset();
#endif
    }
    else if (config->schurApproxMethod == 2)
    {
        // Compute and drop the entries in the Schur complement
        // Ignore the structure of the Schur complement
        //std::cout << "Computing the Approx Schur complement" << std::endl;
#ifdef TIMING_OUTPUT
        ftime.start();
#endif

        // TODO: Pass the operator rather than all the matrices
        if (config->sep_type == 2)
            Sbar = computeApproxSchur(config, ssym, (ssym->G).getRawPtr(),
             (ssym->R).getRawPtr(), (ssym->LP).getRawPtr(),
             (ssym->Solver).getRawPtr(), (ssym->ifSolver).getRawPtr(),
             (ssym->C).getRawPtr(), &LocalDRowMap);
        else
            Sbar = computeApproxWideSchur(config, ssym, (ssym->G).getRawPtr(),
             (ssym->R).getRawPtr(), (ssym->LP).getRawPtr(),
             (ssym->Solver).getRawPtr(), (ssym->ifSolver).getRawPtr(),
             (ssym->C).getRawPtr(),
             &LocalDRowMap);

        //std::cout << *Sbar << std::endl;
#ifdef TIMING_OUTPUT
        ftime.stop();
        std::cout << " Shylu_Factor(" << myPID << ") :: Time to Compute Approx Schur Complement " << ftime.totalElapsedTime() << std::endl;
        ftime.reset();
#endif
        //std::cout << "Computed Approx Schur complement" << std::endl;
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
        // mfh 25 May 2015: assert() gets defined to nothing in a
        // release build.  This causes a build warning that
        // IsAvailable is set but not used.  That's why I protect its
        // declaration with the NDEBUG macro.
#ifdef NDEBUG
        (void) Factory.Query(config->schurAmesosSolver);
#else
        bool IsAvailable = Factory.Query(config->schurAmesosSolver);
        assert(IsAvailable == true);
#endif // NDEBUG
        data->dsolver = Factory.Create(
                        config->schurAmesosSolver, *(data->LP2));

        data->dsolver->SetParameters(aList);
        //std::cout << "Created the direct Schur  Solver" << std::endl;

        data->dsolver->SymbolicFactorization();
        //std::cout << "Symbolic Factorization done for schur complement" << std::endl;

        //std::cout << "In Numeric Factorization of Schur complement" << std::endl;
        data->dsolver->NumericFactorization();
        //std::cout << "Numeric Factorization done for schur complement" << std::endl;

    }
    else if (config->schurSolver == "AztecOO-Exact")
    {
        if (config->libName == "AztecOO") {
                data->innersolver = new AztecOO();
        }
        std::string schurPrec = config->schurPreconditioner;

        err = data->innersolver->SetUserOperator
                    (data->schur_op.get());
        assert (err == 0);

#ifdef HAVE_IFPACK_DYNAMIC_FACTORY
    Ifpack_DynamicFactory IfpackFactory;
#else
    Ifpack IfpackFactory;
#endif
        data->schur_prec = Teuchos::rcp<Ifpack_Preconditioner>
                                (IfpackFactory.Create(schurPrec,
                                         Sbar.getRawPtr(), config->overlap, false));

        data->schur_prec->Initialize();
        data->schur_prec->Compute();

        err = data->innersolver->SetPrecOperator
                    (data->schur_prec.get());
        assert (err == 0);

        data->innersolver->SetAztecOption(AZ_solver, AZ_gmres);

        if (config->silent_subiter) {
                data->innersolver->SetAztecOption(AZ_output, AZ_none);
        }

        data->innersolver->SetMatrixName(999);

    }
    else if (config->schurSolver == "AztecOO-Inexact")
    {
        // Doing an inexact Schur complement
        if (config->libName == "Belos")
        {
            err = data->innersolver->SetUserMatrix
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
            //std::cout << "Condition number of inner Sbar" << condest << std::endl;
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
    else if ((config->schurSolver == "IQR") || (config->schurSolver == "G"))
    {
        // DO NOTHING
    } else
    {
        assert (0 == 1);
    }

    //std::cout << " Out of factor" << std::endl ;
#ifdef TIMING_OUTPUT
    fact_time.stop();
    std::cout << " Shylu_Factor(" << myPID << ") :: Total Time" << fact_time.totalElapsedTime() << std::endl;
    fact_time.reset();
#endif
    return err;
}
