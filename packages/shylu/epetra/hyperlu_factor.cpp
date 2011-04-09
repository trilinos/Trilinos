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

/* This define will make S block diagonal, To make this happen even some
   zero columns/rows of C and R will be stored.
   */
#define BLOCK_DIAGONAL_Si

#include "hyperlu.h"
#include "hyperlu_util.h"
#include "hyperlu_probing_operator.h"


int HyperLU_factor(Epetra_CrsMatrix *A, hyperlu_data *data, hyperlu_config
                *config)
{
    int myPID = A->Comm().MyPID();
    int n = A->NumGlobalRows();

    Epetra_LinearProblem *LP;
    Amesos_BaseSolver *Solver;
    Epetra_CrsMatrix *Cptr;
    int Dnr;
    int Snr;
    int *DRowElems;
    int *SRowElems;
    Teuchos::RCP<Epetra_CrsMatrix> Sbar;
    int sym = config->sym;
    double Sdiagfactor = config->Sdiagfactor;

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
    if (config->schurApproxMethod == 2)
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

    // TODO : The above assignment may not be correct in the unsymetric case

    cout << msg << " #columns in diagonal blk ="<< Dnc << endl;
    cout << msg << " #rows in diagonal blk ="<< Dnr << endl;
    cout << msg << " #columns in S ="<< Snc << endl;
    cout << msg << " #rows in S ="<< Snr << endl;

    // Create a row map for the D and S blocks [
    DRowElems = new int[Dnr];
    SRowElems = new int[Snr];
    int gid;
    // Assemble row ids in two arrays (for D and R blocks)
    if (sym)
    {
        findBlockElems(nrows, rows, gvals, Dnr, DRowElems, Snr, SRowElems, 
                    msg + "D Rows ", msg + "S Rows") ;
    }
    else
    {
        // SRowElems are not known until factorization, TODO
        assert(0 == 1);
    }

    data->Dnr = Dnr;
    data->Snr = Snr;
    data->DRowElems = DRowElems;
    data->SRowElems = SRowElems;

    // Create the local row map 
    Epetra_SerialComm *LComm = new Epetra_SerialComm; // TODO: Fix this
    // Somehow this Comm is not available in Amesos_Pardiso's solve
    //
	// Use Serial Comm for the local blocks.
    Epetra_Map *LocalDRowMap = new Epetra_Map(-1, Dnr, DRowElems, 0, *LComm);
    Epetra_Map DRowMap(-1, Dnr, DRowElems, 0, A->Comm());
    //Epetra_Map LocalSRowMap(-1, Snr, SRowElems, 0, *LComm);
    Epetra_Map SRowMap(-1, Snr, SRowElems, 0, A->Comm());
    // ]

    // Create a column map for the D and S blocks [
    int *DColElems = new int[Dnc]; // Elems in column map of D 
    int *SColElems = new int[Snc]; // Elems in column map of C
    // Assemble column ids in two arrays (for D and C blocks)
    findBlockElems(ncols, cols, gvals, Dnc, DColElems, Snc, SColElems, 
                    "D Cols ", "S Cols") ;

    // Create the local column map 
    Epetra_Map *LocalDColMap = new Epetra_Map(-1, Dnc, DColElems, 0, *LComm);
    Epetra_Map DColMap(-1, Dnc, DColElems, 0, A->Comm());
    //Epetra_Map LocalSColMap(-1, Snc, SColElems, 0, *LComm);
    //Epetra_Map SColMap(SNumGlobalCols, Snc, SColElems, 0, A->Comm());
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
    cout << "No of diagonals in Sbar =" << Sdiag << endl;
    //assert (Sdiag <= SNumGlobalCols-1);
    if (Snr != 0) assert (Sdiag <= Snr-1);

    int dcol, ccol, rcol, scol;
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
                int gcid;
                //cout << "proc/row " << myPID << "/" << gid;
                //cout << " has Cols = " ;
                for (int j = 0 ; j < NumEntries ; j++)
                { // O(nnz) ! Careful what you do inside
                    gcid = A->GCID(Ai[j]); 
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
                SNumEntriesPerRow[scol++] = scnt;
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
                SBarNumEntriesPerRow[sbarcol] = SNumEntriesPerRow[sbarcol] 
                                                + Sdiag*2 + 1;
                sbarcol++;
                //SBarNumEntriesPerRow[i] = SNumEntriesPerRow[i] ;

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


    //cout << " I am here" << endl;
    //Create the local matrices
    Epetra_CrsMatrix *D = new Epetra_CrsMatrix(Copy, *LocalDRowMap,
				*LocalDColMap, DNumEntriesPerRow, true);
    //cout << " Created D matrix" << endl;
    //MPI_Barrier(MPI_COMM_WORLD);
    //Epetra_CrsMatrix S(Copy, SRowMap, SColMap, SNumEntriesPerRow, true);
    Epetra_CrsMatrix S(Copy, SRowMap, SNumEntriesPerRow, true);
    //cout << " Created S matrix" << endl;
    //Cptr = new Epetra_CrsMatrix(Copy, DRowMap, SColMap, 
                                        //CNumEntriesPerRow, true);
    Cptr = new Epetra_CrsMatrix(Copy, DRowMap, CNumEntriesPerRow, true);
    //cout << " Created C matrix" << endl;
    Epetra_CrsMatrix R(Copy, SRowMap, DColMap, RNumEntriesPerRow, true);
    //cout << " Created all the matrices" << endl;
    //Epetra_CrsGraph Sg(Copy, SRowMap, SColMap, SNumEntriesPerRow, true);
    // Leave the column map out, Let Epetra do the work.
    Epetra_CrsGraph Sg(Copy, SRowMap, SBarNumEntriesPerRow, false);
    // for all diagonals
    //Epetra_CrsGraph Sg(Copy, SRowMap, SRowMap, 0, false);
    //cout << " Created Sbar graph" << endl;

    int *LeftIndex = new int[max(Dmaxnnz, Rmaxnnz)];
    double *LeftValues = new double[max(Dmaxnnz, Rmaxnnz)];
    int *RightIndex = new int[max(Cmaxnnz, Smaxnnz)];
    double *RightValues = new double[max(Cmaxnnz, Smaxnnz)];

    //cout << "Inserting values in matrices" << endl;
    int err;
    int lcnt, rcnt ;
    //int debug_row = 0;
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
                if (rcnt >= max(Cmaxnnz, Smaxnnz)) 
                    cout << "rcnt =" << rcnt << "Smaxnnz =" << Smaxnnz << endl;
                assert(rcnt < max(Cmaxnnz, Smaxnnz));
                RightIndex[rcnt] = A->GCID(Ai[j]);
                RightValues[rcnt++] = Ax[j];
            }
        }

        if (gvals[gid] == 1)
        { // D or C row
            err = D->InsertGlobalValues(gid, lcnt, LeftValues, LeftIndex);
            assert(err == 0);
            err = Cptr->InsertGlobalValues(gid, rcnt, RightValues, RightIndex);
            assert(err == 0);
        }
        else
        { // R or S row
            //assert(lcnt > 0); // TODO: Enable this once using narrow sep.
            assert(rcnt > 0);
            err = R.InsertGlobalValues(gid, lcnt, LeftValues, LeftIndex);
            assert(err == 0);
            err = S.InsertGlobalValues(gid, rcnt, RightValues, RightIndex);
            assert(err == 0);
            if (config->schurApproxMethod == 1)
            {
                err = Sg.InsertGlobalIndices(gid, rcnt, RightIndex);
                assert(err == 0);
            }
        }
    }
    //cout << "Done Inserting values in matrices" << endl;
    D->FillComplete();
    //cout << "Done D fill complete" << endl;
    S.FillComplete();
    //cout << "Done S fill complete" << endl;
    //Cptr->FillComplete(SColMap, DRowMap);
    Cptr->FillComplete(SRowMap, DRowMap); //TODO: Won't work if permutation is
                                            // unsymmetric SRowMap
    //cout << "Done C fill complete" << endl;
    R.FillComplete(DColMap, SRowMap);
    //cout << "Done R fill complete" << endl;

    // Add the diagonals to Sg
    for (int i = 0; config->schurApproxMethod == 1 && i < nrows ; i++)
    {
        gid = rows[i];
        if (gvals[gid] == 1) continue; // not a row in S
        if (Snr == 0) assert(0 == 1);

        rcnt = 0;
        //TODO Will be trouble if SNumGlobalCols != Snc
        //assert(SNumGlobalCols == Snc);
        //for (int j = MAX(i-Sdiag, 0) ; j < MIN(SNumGlobalCols, i+Sdiag); j++)
        for (int j = MAX(i-Sdiag, 0) ; j < MIN(Snr, i+Sdiag); j++)
        {
            // find the adjacent columns from the row map of S
            //assert (j >= 0 && j < Snr);
            RightIndex[rcnt++] = SRowElems[j];
        }
        err = Sg.InsertGlobalIndices(gid, rcnt, RightIndex);
        assert(err == 0);
        // Always insert the diagonals, if it is added twice that is fine.
        err = Sg.InsertGlobalIndices(gid, 1, &gid);
        assert(err == 0);
    }
    /*// This is for all diagonals, Should remove later
    int totalElems = SRowMap.NumGlobalElements();
    cout << "totalElems" << totalElems << endl;

    int *allSGID = new int[totalElems];   // vector of size Schur complement !
    getAllSGIDs(Snr, SRowElems, allSGID);

    // OUT OF DATE, 2 is used for something else
    for (int i = 0; config->schurApproxMethod == 2 && i < nrows ; i++)
    {
        rcnt = 0;
        cout << msg << "i=" << i ;
        //for (int j = MAX(i-Sdiag, 0) ; j < MIN(SNumGlobalCols, i+Sdiag); j++)
        for (int j = MAX(i-Sdiag, 0) ; j < MIN(totalElems, i+Sdiag); j++)
        {
            // find the adjacent columns from the row map of S
            //assert (j >= 0 && j < Snr);
            RightIndex[rcnt++] = j;
            cout << "j=" << j ;
        }
        err = Sg.InsertMyIndices(i, rcnt, RightIndex);
        cout << "Error = " << err << endl;
        assert(err == 0);
        // Always insert the diagonals, if it is added twice that is fine.
        err = Sg.InsertMyIndices(i, 1, &i);
        assert(err == 0);
    }*/

    if (config->schurApproxMethod == 1)
        Sg.FillComplete();
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
    //cout << msg << "D rows=" << D->NumGlobalRows() << " D cols=" <<
        //D->NumGlobalCols() << "#cols in column map="<<
        //D->ColMap().NumMyElements() << endl;
    //cout << msg << "R rows=" << R.NumGlobalRows() << " R cols=" <<
        //R.NumGlobalCols() << "#cols in column map="<<
        //R.ColMap().NumMyElements() << endl;
    // ]

    // ======================= Numeric factorization =========================
    //Amesos_BaseSolver *Solver;
    Amesos Factory;
    char* SolverType = "Amesos_Klu";
    bool IsAvailable = Factory.Query(SolverType);
    assert(IsAvailable == true);

    // TODO : All these three vectors and the solve can be removed
    /*Epetra_MultiVector *X = new Epetra_MultiVector(LocalDRowMap, 1);
    Epetra_MultiVector *residual = new Epetra_MultiVector(LocalDRowMap, 1);
    Epetra_MultiVector *B = new Epetra_MultiVector(LocalDRowMap, 1);
    B->PutScalar(1.0);

    double *resid = new double[Snr];*/

    LP = new Epetra_LinearProblem();
    LP->SetOperator(D);
    /*LP->SetLHS(X);
    LP->SetRHS(B);*/
    Solver = Factory.Create(SolverType, *LP);
    cout << "Created the diagonal Solver" << endl;

    Solver->SymbolicFactorization();
    cout << "Symbolic Factorization done" << endl;
    Teuchos::Time ftime("setup time");
    ftime.start();

    cout << "In Numeric Factorization" << endl;
    Solver->NumericFactorization();
    cout << "Numeric Factorization done" << Solver << endl;

    data->SComm = LComm;
    data->LP = LP;
    data->Solver = Solver;
    data->D = D;
    data->Cptr = Cptr;
    data->LDRowMap = LocalDRowMap;
    data->LDColMap = LocalDColMap;

    /*Solver->Solve();
    cout << "Solve done" << endl;

    ftime.stop();
    cout << "Time to factor" << ftime.totalElapsedTime() << endl;
    ftime.reset();

    D->Multiply(false, *X, *residual);
    residual->Update(-1.0, *B, 1.0);
    residual->Norm2(resid);
    double max_resid = 0;
    for (int i = 0 ; i < Snr; i++)
        max_resid = max(max_resid, abs(resid[i]));

    cout << "Maximum residual = " << max_resid << endl;*/


#ifdef DUMP_MATRICES
    if (myPID == 1)
    {
        EpetraExt::MultiVectorToMatlabFile("LHS.mat", *X);
        EpetraExt::RowMatrixToMatlabFile("D.mat", D);
    }
#endif

    if (config->schurApproxMethod == 1)
    {
        //Set up the probing operator
        HyperLU_Probing_Operator probeop(&S, &R, LP, Solver, Cptr,
                                        LocalDRowMap);

        Teuchos::ParameterList pList;
        Teuchos::RCP<const Epetra_CrsGraph> rSg = Teuchos::rcpFromRef(Sg);
        Isorropia::Epetra::Prober prober(rSg, pList, false);
        cout << Sg << endl;
        cout << "Doing coloring" << endl;
        ftime.start();
        prober.color();
        ftime.stop();
        cout << "Time to color" << ftime.totalElapsedTime() << endl;
        ftime.reset();
        cout << "Doing probing" << endl;
        ftime.start();
        Sbar = prober.probe(probeop);
        cout << "SIZE of SBAR = " << (*Sbar).NumGlobalRows() << endl;
        ftime.stop();
        cout << "Time to probe" << ftime.totalElapsedTime() << endl;
        ftime.reset();
    }
    else if (config->schurApproxMethod == 2)
    {
        // Compute and drop the entries in the Schur complement
        // Ignore the structure of the Schur complement
        cout << "Computing the Approx Schur complement" << endl;
        Sbar = computeApproxSchur(config, &S, &R, LP, Solver, Cptr,
                                        LocalDRowMap);
        cout << "Approx Schur complement is done" << endl;
    }

    data->Sbar  = Sbar;

    // ======================= Solve =========================================

    // Clean up
    /*delete residual;
    delete[] resid;
    delete X;*/
    //delete B;
    delete[] DNumEntriesPerRow;
    delete[] SNumEntriesPerRow;
    delete[] CNumEntriesPerRow;
    delete[] RNumEntriesPerRow;
    delete[] SBarNumEntriesPerRow;
    delete[] DColElems;
    delete[] SColElems;
    //delete[] CColMax;
    //delete[] RRowMax;
    //delete[] CColElems;
    //delete[] RRowElems;
    //delete[] DRowElems;
    //delete[] SRowElems;

    delete[] gvals;
    //delete A;
    //cout << " Out of factor" << endl ;
    return 0;
}
