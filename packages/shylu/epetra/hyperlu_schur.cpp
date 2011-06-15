
#include "hyperlu.h"
#include "hyperlu_util.h"
//#include "EpetraExt_Transpose_RowMatrix.h"
#include "Epetra_Vector.h"
#include "hyperlu_probing_operator.h"

/* Apply an identity matrix to the Schur complement operator. Drop the entries
   entries using a relative threshold. Assemble the result in a Crs Matrix
   which will be our approximate Schur complement.
   */
Teuchos::RCP<Epetra_CrsMatrix> computeApproxSchur(hyperlu_config *config,
    Epetra_CrsMatrix *G, Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver, Epetra_CrsMatrix *C,
    Epetra_Map *localDRowMap)
{
    double relative_thres = config->relative_threshold;
    int nvectors = 16;

    HyperLU_Probing_Operator probeop(G, R, LP, solver, C, localDRowMap,
                                    nvectors);

    // Get row map
    Epetra_Map rMap = G->RowMap();
    int *rows = rMap.MyGlobalElements();
    int totalElems = rMap.NumGlobalElements();
    int localElems = rMap.NumMyElements();
    //cout << " totalElems in Schur Complement" << totalElems << endl;
    //cout << myPID << " localElems" << localElems << endl;

    // **************** Two collectives here *********************
#ifdef TIMING_OUTPUT
    Teuchos::Time ftime("setup time");
    ftime.start();
#endif
    int prefixSum;
    G->Comm().ScanSum(&localElems, &prefixSum, 1);
    //cout << " prefixSum" << prefixSum << endl;
    // Start the index in prefixSum-localElems
    int *mySGID = new int[totalElems];   // vector of size Schur complement !
    int *allSGID = new int[totalElems];   // vector of size Schur complement !
    int i, j;
    for (i = 0, j = 0; i < totalElems ; i++)
    {
        if (i >= prefixSum - localElems && i < prefixSum)
        {
            mySGID[i] = rows[j];
            j++;
        }
        else
        {
            mySGID[i] = 0;
        }
        allSGID[i] = 0;
    }

    C->Comm().SumAll(mySGID, allSGID, totalElems);

#ifdef TIMING_OUTPUT
    ftime.stop();
    cout << "Time to Compute RowIDS" << ftime.totalElapsedTime() << endl;
    ftime.reset();
#endif
    // Now everyone knows the GIDs in the Schur complement

    //cout << rMap << endl;
    j = 0;
    Teuchos::RCP<Epetra_CrsMatrix> Sbar = Teuchos::rcp(new Epetra_CrsMatrix(
                                            Copy, rMap, localElems));
    int nentries;
    double *values = new double[localElems]; // Need to adjust this for more
    int *indices = new int[localElems];      // than one vector
    double *vecvalues;
    int dropped = 0;
    double *maxvalue = new double[nvectors];
#ifdef TIMING_OUTPUT
    ftime.start();
#endif
#ifdef TIMING_OUTPUT
    Teuchos::Time app_time("Apply time");
#endif
    int findex = totalElems / nvectors ;
    for (i = 0 ; i < findex*nvectors ; i+=nvectors)
    {
        Epetra_MultiVector probevec(rMap, nvectors);
        Epetra_MultiVector Scol(rMap, nvectors);

        probevec.PutScalar(0.0);
        int cindex;
        for (int k = 0; k < nvectors; k++)
        {
            cindex = k+i;
            if (cindex >= prefixSum - localElems && cindex < prefixSum)
            {
                probevec.ReplaceGlobalValue(allSGID[cindex], k, 1.0);
            }
        }

#ifdef TIMING_OUTPUT
        app_time.start();
#endif
        probeop.Apply(probevec, Scol);
#ifdef TIMING_OUTPUT
        app_time.stop();
#endif
        //cout << Scol << endl;

        Scol.MaxValue(maxvalue);
        for (int k = 0; k < nvectors; k++) //TODO:Need to switch these loops
        {
            cindex = k+i;
            vecvalues = Scol[k];
            //cout << "MAX" << maxvalue << endl;
            for (j = 0 ; j < localElems ; j++)
            {
                nentries = 0; // inserting one entry in each row for now
                if (allSGID[cindex] == rows[j]) // diagonal entry
                {
                    values[nentries] = vecvalues[j];
                    indices[nentries] = allSGID[cindex];
                    nentries++;
                    Sbar->InsertGlobalValues(rows[j], nentries, values, indices);
                }
                else if (abs(vecvalues[j]/maxvalue[k]) > relative_thres)
                {
                    values[nentries] = vecvalues[j];
                    indices[nentries] = allSGID[cindex];
                    nentries++;
                    Sbar->InsertGlobalValues(rows[j], nentries, values, indices);
                }
                else
                {
                    if (vecvalues[j] != 0.0) dropped++;
                }
            }
        }
    }

    probeop.ResetTempVectors(1);

    for ( ; i < totalElems ; i++)
    {
        Epetra_MultiVector probevec(rMap, 1); // TODO: Try doing more than one
        Epetra_MultiVector Scol(rMap, 1);     // vector at a time

        probevec.PutScalar(0.0);
        if (i >= prefixSum - localElems && i < prefixSum)
        {
            probevec.ReplaceGlobalValue(allSGID[i], 0, 1.0);
        }

#ifdef TIMING_OUTPUT
        app_time.start();
#endif
        probeop.Apply(probevec, Scol);
#ifdef TIMING_OUTPUT
        app_time.stop();
#endif
        //cout << Scol << endl;

        vecvalues = Scol[0];
        Scol.MaxValue(maxvalue);
        //cout << "MAX" << maxvalue << endl;
        for (j = 0 ; j < localElems ; j++)
        {
            nentries = 0; // inserting one entry in each row for now
            if (allSGID[i] == rows[j]) // diagonal entry
            {
                values[nentries] = vecvalues[j];
                indices[nentries] = allSGID[i];
                nentries++;
                Sbar->InsertGlobalValues(rows[j], nentries, values, indices);
            }
            else if (abs(vecvalues[j]/maxvalue[0]) > relative_thres)
            {
                values[nentries] = vecvalues[j];
                indices[nentries] = allSGID[i];
                nentries++;
                Sbar->InsertGlobalValues(rows[j], nentries, values, indices);
            }
            else
            {
                if (vecvalues[j] != 0.0) dropped++;
            }
        }
    }
#ifdef TIMING_OUTPUT
    ftime.stop();
    cout << "Time in finding and dropping entries" << ftime.totalElapsedTime() << endl;
    ftime.reset();
#endif
#ifdef TIMING_OUTPUT
    cout << "Time in Apply of probing" << app_time.totalElapsedTime() << endl;
#endif
    Sbar->FillComplete();
    cout << "#dropped entries" << dropped << endl;
    delete[] allSGID;
    delete[] mySGID;
    delete[] values;
    delete[] indices;
    delete[] maxvalue;

    return Sbar;
}

/* Computes the approximate Schur complement for the wide separator */
Teuchos::RCP<Epetra_CrsMatrix> computeApproxWideSchur(hyperlu_config *config,
    Epetra_CrsMatrix *G, Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver, Epetra_CrsMatrix *C,
    Epetra_Map *localDRowMap)
{
    int i;
    double relative_thres = config->relative_threshold;

    // Need to create local G (block diagonal portion) , R, C

    // Get row map of G
    Epetra_Map CrMap = C->RowMap();
    int *c_rows = CrMap.MyGlobalElements();
    int *c_cols = (C->ColMap()).MyGlobalElements();
    //int c_totalElems = CrMap.NumGlobalElements();
    int c_localElems = CrMap.NumMyElements();
    int c_localcolElems = (C->ColMap()).NumMyElements();

    Epetra_Map GrMap = G->RowMap();
    int *g_rows = GrMap.MyGlobalElements();
    //int g_totalElems = GrMap.NumGlobalElements();
    int g_localElems = GrMap.NumMyElements();

    Epetra_Map RrMap = R->RowMap();
    int *r_rows = RrMap.MyGlobalElements();
    int *r_cols = (R->ColMap()).MyGlobalElements();
    //int r_totalElems = RrMap.NumGlobalElements();
    int r_localElems = RrMap.NumMyElements();
    int r_localcolElems = (R->ColMap()).NumMyElements();

    Epetra_SerialComm LComm;
    Epetra_Map C_localRMap (-1, c_localElems, c_rows, 0, LComm);
    Epetra_Map C_localCMap (-1, c_localcolElems, c_cols, 0, LComm);
    Epetra_Map G_localRMap (-1, g_localElems, g_rows, 0, LComm);
    Epetra_Map R_localRMap (-1, r_localElems, r_rows, 0, LComm);
    Epetra_Map R_localCMap (-1, r_localcolElems, r_cols, 0, LComm);

    cout << "#local rows" << g_localElems << "#non zero local cols" << c_localcolElems << endl;

    int nentries1, gid;
    // maxentries is the maximum of all three possible matrices as the arrays
    // are reused between the three
    int maxentries = max(C->MaxNumEntries(), R->MaxNumEntries());
    maxentries = max(maxentries, G->MaxNumEntries());

    double *values1 = new double[maxentries];
    double *values2 = new double[maxentries];
    double *values3 = new double[maxentries];
    int *indices1 = new int[maxentries];
    int *indices2 = new int[maxentries];
    int *indices3 = new int[maxentries];

    //cout << *C << endl;
    //cout << "Creating local matrices" << endl;
    int err;
    Epetra_CrsMatrix localC(Copy, C_localRMap, C->MaxNumEntries(), false);
    for (i = 0; i < c_localElems ; i++)
    {
        gid = c_rows[i];
        err = C->ExtractGlobalRowCopy(gid, maxentries, nentries1, values1, indices1);
        assert (err == 0);
        //if (nentries1 > 0) // TODO: Later
        //{
        err = localC.InsertGlobalValues(gid, nentries1, values1, indices1);
        assert (err == 0);
        //}
    }
    localC.FillComplete(C_localCMap, C_localRMap);
    //cout << "Created local C matrix" << endl;

    Epetra_CrsMatrix localR(Copy, R_localRMap, R->MaxNumEntries(), false);
    for (i = 0; i < r_localElems ; i++)
    {
        gid = r_rows[i];
        R->ExtractGlobalRowCopy(gid, maxentries, nentries1, values1, indices1);
        localR.InsertGlobalValues(gid, nentries1, values1, indices1);
    }
    localR.FillComplete(R_localCMap, R_localRMap);
    //cout << "Created local R matrix" << endl;

    // Sbar - Approximate Schur complement
    Teuchos::RCP<Epetra_CrsMatrix> Sbar = Teuchos::rcp(new Epetra_CrsMatrix(
                                            Copy, GrMap, g_localElems));

    // Include only the diagonal elements of G in localG
    Epetra_CrsMatrix localG(Copy, G_localRMap, G->MaxNumEntries(), false);
    int cnt, scnt;
    for (i = 0; i < g_localElems ; i++)
    {
        gid = g_rows[i];
        G->ExtractGlobalRowCopy(gid, maxentries, nentries1, values1, indices1);

        cnt = 0;
        scnt = 0;
        for (int j = 0 ; j < nentries1 ; j++)
        {
            if (G->LRID(indices1[j]) != -1)
            {
                values2[cnt] = values1[j];
                indices2[cnt++] = indices1[j];
            }
            else
            {
                // Add it to Sbar immediately
                values3[scnt] = values1[j];
                indices3[scnt++] = indices1[j];
            }
        }

        localG.InsertGlobalValues(gid, cnt, values2, indices2);
        Sbar->InsertGlobalValues(gid, scnt, values3, indices3);
    }
    localG.FillComplete();
    cout << "Created local G matrix" << endl;

    int nvectors = 16;
    HyperLU_Probing_Operator probeop(&localG, &localR, LP, solver, &localC,
                                        localDRowMap, nvectors);

    //cout << " totalElems in Schur Complement" << totalElems << endl;
    //cout << myPID << " localElems" << localElems << endl;

    // **************** Two collectives here *********************
#ifdef TIMING_OUTPUT
    Teuchos::Time ftime("setup time");
    ftime.start();
#endif
#ifdef TIMING_OUTPUT
    Teuchos::Time app_time("setup time");
#endif

    int nentries;
    // size > maxentries as there could be fill
    // TODO: Currently the size of the two arrays can be one, Even if we switch
    // the loop below the size of the array required is nvectors. Fix it
    double *values = new double[g_localElems];
    int *indices = new int[g_localElems];
    double *vecvalues;
    int dropped = 0;
    double *maxvalue = new double[nvectors];
#ifdef TIMING_OUTPUT
    ftime.start();
#endif
    int findex = g_localElems / nvectors ;

    int cindex;
    int mypid = C->Comm().MyPID();
    Epetra_MultiVector probevec(G_localRMap, nvectors);
    Epetra_MultiVector Scol(G_localRMap, nvectors);
    for (i = 0 ; i < findex*nvectors ; i+=nvectors)
    {
        probevec.PutScalar(0.0);
        for (int k = 0; k < nvectors; k++)
        {
            cindex = k+i;
            // TODO: Can do better than this, just need to go to the column map
            // of C, there might be null columns in C
            probevec.ReplaceGlobalValue(g_rows[cindex], k, 1.0);
            //if (mypid == 0)
            //cout << "Changing row to 1.0 " << g_rows[cindex] << endl;
        }

        //if (mypid == 0)
        //cout << probevec << endl;

#ifdef TIMING_OUTPUT
        app_time.start();
#endif
        probeop.Apply(probevec, Scol);
#ifdef TIMING_OUTPUT
        app_time.stop();
#endif
        //if (mypid == 0)
        //cout << Scol << endl;

        Scol.MaxValue(maxvalue);
        for (int k = 0; k < nvectors; k++) //TODO:Need to switch these loops
        {
            cindex = k+i;
            vecvalues = Scol[k];
            //cout << "MAX" << maxvalue << endl;
            for (int j = 0 ; j < g_localElems ; j++)
            {
                nentries = 0; // inserting one entry in each row for now
                if (g_rows[cindex] == g_rows[j]) // diagonal entry
                {
                    values[nentries] = vecvalues[j];
                    indices[nentries] = g_rows[cindex];
                    nentries++;
                    Sbar->InsertGlobalValues(g_rows[j], nentries, values, indices);
                }
                else if (abs(vecvalues[j]/maxvalue[k]) > relative_thres)
                {
                    values[nentries] = vecvalues[j];
                    indices[nentries] = g_rows[cindex];
                    nentries++;
                    Sbar->InsertGlobalValues(g_rows[j], nentries, values, indices);
                }
                else
                {
                    if (vecvalues[j] != 0.0) dropped++;
                }
            }
        }
    }

    probeop.ResetTempVectors(1);

    for ( ; i < g_localElems ; i++)
    {
        // TODO: Can move the next two decalarations outside the loop
        Epetra_MultiVector probevec(G_localRMap, 1);
        Epetra_MultiVector Scol(G_localRMap, 1);

        probevec.PutScalar(0.0);
        // TODO: Can do better than this, just need to go to the column map
        // of C, there might be null columns in C
        probevec.ReplaceGlobalValue(g_rows[i], 0, 1.0);

#ifdef TIMING_OUTPUT
        app_time.start();
#endif
        probeop.Apply(probevec, Scol);
#ifdef TIMING_OUTPUT
        app_time.stop();
#endif
        //cout << Scol << endl;

        vecvalues = Scol[0];
        Scol.MaxValue(maxvalue);
        //cout << "MAX" << maxvalue << endl;
        for (int j = 0 ; j < g_localElems ; j++)
        {
            nentries = 0; // inserting one entry in each row for now
            if (g_rows[i] == g_rows[j]) // diagonal entry
            {
                values[nentries] = vecvalues[j];
                indices[nentries] = g_rows[i];
                nentries++;
                Sbar->InsertGlobalValues(g_rows[j], nentries, values, indices);
            }
            else if (abs(vecvalues[j]/maxvalue[0]) > relative_thres)
            {
                values[nentries] = vecvalues[j];
                indices[nentries] = g_rows[i];
                nentries++;
                Sbar->InsertGlobalValues(g_rows[j], nentries, values, indices);
            }
            else
            {
                if (vecvalues[j] != 0.0) dropped++;
            }
        }
    }

    // Still need to add the off diagonal entries of G
    // Should we do dropping here.

#ifdef TIMING_OUTPUT
    ftime.stop();
    cout << "Time in finding and dropping entries" << ftime.totalElapsedTime() << endl;
    ftime.reset();
#endif
#ifdef TIMING_OUTPUT
    cout << "Time in Apply of probing" << app_time.totalElapsedTime() << endl;
#endif
    probeop.PrintTimingInfo();
    Sbar->FillComplete();
    cout << "#dropped entries" << dropped << endl;
    delete[] values;
    delete[] indices;
    delete[] values1;
    delete[] indices1;
    delete[] values2;
    delete[] indices2;
    delete[] values3;
    delete[] indices3;
    delete[] maxvalue;

    return Sbar;
}
