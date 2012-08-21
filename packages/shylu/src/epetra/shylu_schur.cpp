
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


#include "shylu.h"
#include "shylu_util.h"
//#include "EpetraExt_Transpose_RowMatrix.h"
#include "Epetra_Vector.h"
#include "shylu_probing_operator.h"
#include "shylu_local_schur_operator.h"
#include <EpetraExt_Reindex_CrsMatrix.h>

/* Apply an identity matrix to the Schur complement operator. Drop the entries
   entries using a relative threshold. Assemble the result in a Crs Matrix
   which will be our approximate Schur complement.
   */
Teuchos::RCP<Epetra_CrsMatrix> computeApproxSchur(shylu_config *config,
    shylu_symbolic *sym,
    Epetra_CrsMatrix *G, Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver, Epetra_CrsMatrix *C,
    Epetra_Map *localDRowMap)
{
    double relative_thres = config->relative_threshold;
    int nvectors = 16;

    ShyLU_Probing_Operator probeop(sym, G, R, LP, solver, C, localDRowMap,
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
Teuchos::RCP<Epetra_CrsMatrix> computeApproxWideSchur(shylu_config *config,
    shylu_symbolic *ssym,   // symbolic structure
    Epetra_CrsMatrix *G, Epetra_CrsMatrix *R,
    Epetra_LinearProblem *LP, Amesos_BaseSolver *solver, Epetra_CrsMatrix *C,
    Epetra_Map *localDRowMap)
{
    int i;
    double relative_thres = config->relative_threshold;

    // Need to create local G (block diagonal portion) , R, C

    // Get row map of G
    //Epetra_Map CrMap = C->RowMap();
    //int *c_rows = CrMap.MyGlobalElements();
    //int *c_cols = (C->ColMap()).MyGlobalElements();
    //int c_totalElems = CrMap.NumGlobalElements();
    //int c_localElems = CrMap.NumMyElements();
    //int c_localcolElems = (C->ColMap()).NumMyElements();

    Epetra_Map GrMap = G->RowMap();
    int *g_rows = GrMap.MyGlobalElements();
    //int g_totalElems = GrMap.NumGlobalElements();
    int g_localElems = GrMap.NumMyElements();

    //Epetra_Map RrMap = R->RowMap();
    //int *r_rows = RrMap.MyGlobalElements();
    //int *r_cols = (R->ColMap()).MyGlobalElements();
    //int r_totalElems = RrMap.NumGlobalElements();
    //int r_localElems = RrMap.NumMyElements();
    //int r_localcolElems = (R->ColMap()).NumMyElements();

    Epetra_SerialComm LComm;
    Epetra_Map G_localRMap (-1, g_localElems, g_rows, 0, LComm);

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

    // Sbar - Approximate Schur complement
    Teuchos::RCP<Epetra_CrsMatrix> Sbar = Teuchos::rcp(new Epetra_CrsMatrix(
                                            Copy, GrMap, g_localElems));

    // Include only the block diagonal elements of G in localG
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
    //cout << "Created local G matrix" << endl;

    int nvectors = 16;
    /*ShyLU_Probing_Operator probeop(&localG, &localR, LP, solver, &localC,
                                        localDRowMap, nvectors);*/
    ShyLU_Local_Schur_Operator probeop(ssym, &localG, R, LP, solver, C,
                                        localDRowMap, nvectors);

#ifdef DUMP_MATRICES
    //ostringstream fnamestr;
    //fnamestr << "localC" << C->Comm().MyPID() << ".mat";
    //string Cfname = fnamestr.str();
    //EpetraExt::RowMatrixToMatlabFile(Cfname.c_str(), localC);

    //Epetra_Map defMapg(-1, g_localElems, 0, localG.Comm());
    //EpetraExt::ViewTransform<Epetra_CrsMatrix> * ReIdx_MatTransg =
                        //new EpetraExt::CrsMatrix_Reindex( defMapg );
    //Epetra_CrsMatrix t2G = (*ReIdx_MatTransg)( localG );
    //ReIdx_MatTransg->fwd();
    //EpetraExt::RowMatrixToMatlabFile("localG.mat", t2G);
#endif

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
    double *values = new double[nvectors];
    int *indices = new int[nvectors];
    double *vecvalues;
    int dropped = 0;
    double *maxvalue = new double[nvectors];
#ifdef TIMING_OUTPUT
    ftime.start();
#endif
    int findex = g_localElems / nvectors ;

    int cindex;
    int mypid = C->Comm().MyPID();
    Epetra_MultiVector probevec (G_localRMap, nvectors);
    Epetra_MultiVector Scol (G_localRMap, nvectors);
    probevec.PutScalar(0.0);
    for (i = 0 ; i < findex*nvectors ; i+=nvectors)
    {
        // Set the probevec to find block columns of S.
        for (int k = 0; k < nvectors; k++)
        {
            cindex = k+i;
            // TODO: Can do better than this, just need to go to the column map
            // of C, there might be null columns in C
            probevec.ReplaceGlobalValue(g_rows[cindex], k, 1.0);
            //if (mypid == 0)
            //cout << "Changing row to 1.0 " << g_rows[cindex] << endl;
        }

#ifdef TIMING_OUTPUT
        app_time.start();
#endif
        probeop.Apply(probevec, Scol);
#ifdef TIMING_OUTPUT
        app_time.stop();
#endif

        // Reset the probevec to all zeros.
        for (int k = 0; k < nvectors; k++)
        {
            cindex = k+i;
            probevec.ReplaceGlobalValue(g_rows[cindex], k, 0.0);
        }

        Scol.MaxValue(maxvalue);
        nentries = 0;
        for (int j = 0 ; j < g_localElems ; j++)
        {
            for (int k = 0; k < nvectors; k++)
            {
                cindex = k+i;
                vecvalues = Scol[k];
                if ((g_rows[cindex] == g_rows[j])  ||
                (abs(vecvalues[j]/maxvalue[k]) > relative_thres))
                // diagonal entry or large entry.
                {
                    values[nentries] = vecvalues[j];
                    indices[nentries++] = g_rows[cindex];
                }
#ifdef SHYLU_DEBUG
                else
                {
                    if (vecvalues[j] != 0.0)
                    {
                        dropped++;
                    }
                }
#endif
            }
            Sbar->InsertGlobalValues(g_rows[j], nentries, values,
                        indices);
            nentries = 0;
        }
    }

    if (i < g_localElems)
    {
        nvectors = g_localElems - i;
        probeop.ResetTempVectors(nvectors);
        Epetra_MultiVector probevec1 (G_localRMap, nvectors);
        Epetra_MultiVector Scol1 (G_localRMap, nvectors);

        probevec1.PutScalar(0.0);
        for (int k = 0; k < nvectors; k++)
        {
            cindex = k+i;
            // TODO: Can do better than this, just need to go to the column map
            // of C, there might be null columns in C
            probevec1.ReplaceGlobalValue(g_rows[cindex], k, 1.0);
        }

#ifdef TIMING_OUTPUT
        app_time.start();
#endif
        probeop.Apply(probevec1, Scol1);
#ifdef TIMING_OUTPUT
        app_time.stop();
#endif
        Scol1.MaxValue(maxvalue);
        nentries = 0;
        for (int j = 0 ; j < g_localElems ; j++)
        {
            //cout << "MAX" << maxvalue << endl;
            for (int k = 0; k < nvectors; k++)
            {
                cindex = k+i;
                vecvalues = Scol1[k];
                //nentries = 0; // inserting one entry in each row for now
                if ((g_rows[cindex] == g_rows[j])  ||
                (abs(vecvalues[j]/maxvalue[k]) > relative_thres))
                // diagonal entry or large entry.
                {
                    values[nentries] = vecvalues[j];
                    indices[nentries++] = g_rows[cindex];
                }
#ifdef SHYLU_DEBUG
                else
                {
                    if (vecvalues[j] != 0.0)
                    {
                        dropped++;
                    }
                }
#endif
            }
            Sbar->InsertGlobalValues(g_rows[j], nentries, values,
                        indices);
            nentries = 0;
        }
    }

#ifdef TIMING_OUTPUT
    ftime.stop();
    cout << "Time in finding and dropping entries" << ftime.totalElapsedTime()
                     << endl;
    ftime.reset();
    cout << "Time in Apply of probing" << app_time.totalElapsedTime() << endl;
    probeop.PrintTimingInfo();
#endif
    Sbar->FillComplete();

#ifdef DUMP_MATRICES
    Epetra_Map defMap2(-1, g_localElems, 0, C->Comm());
    EpetraExt::ViewTransform<Epetra_CrsMatrix> * ReIdx_MatTrans2 =
                        new EpetraExt::CrsMatrix_Reindex( defMap2 );
    Epetra_CrsMatrix t2S = (*ReIdx_MatTrans2)( *Sbar );
    ReIdx_MatTrans2->fwd();
    EpetraExt::RowMatrixToMatlabFile("Schur.mat", t2S);
#endif

#ifdef SHYLU_DEBUG
    cout << "#dropped entries" << dropped << endl;
#endif
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

/* Computes the approximate Schur complement for the wide separator
   using guided probing*/
Teuchos::RCP<Epetra_CrsMatrix> computeSchur_GuidedProbing
(
    shylu_config *config,
    shylu_symbolic *ssym,   // symbolic structure
    shylu_data *data,       // numeric structure
    Epetra_Map *localDRowMap
)
{
    int i;
    double relative_thres = config->relative_threshold;

    Epetra_CrsMatrix *G = ssym->G.getRawPtr();
    Epetra_CrsMatrix *R = ssym->R.getRawPtr();
    Epetra_LinearProblem *LP = ssym->LP.getRawPtr();
    Amesos_BaseSolver *solver = ssym->Solver.getRawPtr();
    Epetra_CrsMatrix *C = ssym->C.getRawPtr();

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

    //cout << "#local rows" << g_localElems << "#non zero local cols" << c_localcolElems << endl;

#ifdef DEBUG
    cout << "DEBUG MODE" << endl;
    int nrows = C->RowMap().NumMyElements();
    assert(nrows == localDRowMap->NumGlobalElements());

    int gids[nrows], gids1[nrows];
    C_localRMap.MyGlobalElements(gids);
    localDRowMap->MyGlobalElements(gids1);
    cout << "Comparing R's domain map with D's row map" << endl;

    for (int i = 0; i < nrows; i++)
    {
       assert(gids[i] == gids1[i]);
    }
#endif

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

    //cout << "Creating local matrices" << endl;
    int err;
    Epetra_CrsMatrix localC(Copy, C_localRMap, C->MaxNumEntries(), false);
    for (i = 0; i < c_localElems ; i++)
    {
        gid = c_rows[i];
        err = C->ExtractGlobalRowCopy(gid, maxentries, nentries1, values1,
                                        indices1);
        assert (err == 0);
        //if (nentries1 > 0) // TODO: Later
        //{
        err = localC.InsertGlobalValues(gid, nentries1, values1, indices1);
        assert (err == 0);
        //}
    }
    localC.FillComplete(G_localRMap, C_localRMap);

    //cout << "Created local C matrix" << endl;

    Epetra_CrsMatrix localR(Copy, R_localRMap, R->MaxNumEntries(), false);
    for (i = 0; i < r_localElems ; i++)
    {
        gid = r_rows[i];
        R->ExtractGlobalRowCopy(gid, maxentries, nentries1, values1, indices1);
        localR.InsertGlobalValues(gid, nentries1, values1, indices1);
    }
    localR.FillComplete(*localDRowMap, R_localRMap);
    //cout << "Created local R matrix" << endl;

    // Sbar - Approximate Schur complement
    Teuchos::RCP<Epetra_CrsMatrix> Sbar = Teuchos::rcp(new Epetra_CrsMatrix(
                                            Copy, GrMap, g_localElems));

    // Include only the block diagonal elements of G in localG
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
    ShyLU_Probing_Operator probeop(ssym, &localG, &localR, LP, solver, &localC,
                                        localDRowMap, nvectors);

#ifdef DUMP_MATRICES
    //ostringstream fnamestr;
    //fnamestr << "localC" << C->Comm().MyPID() << ".mat";
    //string Cfname = fnamestr.str();
    //EpetraExt::RowMatrixToMatlabFile(Cfname.c_str(), localC);

    //Epetra_Map defMapg(-1, g_localElems, 0, localG.Comm());
    //EpetraExt::ViewTransform<Epetra_CrsMatrix> * ReIdx_MatTransg =
                        //new EpetraExt::CrsMatrix_Reindex( defMapg );
    //Epetra_CrsMatrix t2G = (*ReIdx_MatTransg)( localG );
    //ReIdx_MatTransg->fwd();
    //EpetraExt::RowMatrixToMatlabFile("localG.mat", t2G);
#endif

    //cout << " totalElems in Schur Complement" << totalElems << endl;
    //cout << myPID << " localElems" << localElems << endl;

    // **************** Two collectives here *********************
#ifdef TIMING_OUTPUT
    Teuchos::Time ftime("setup time");
#endif
#ifdef TIMING_OUTPUT
    Teuchos::Time app_time("setup time");
#endif

    Teuchos::RCP<Epetra_CrsGraph> lSGraph = Teuchos::RCP<Epetra_CrsGraph> (
                    new Epetra_CrsGraph(Copy, G_localRMap, maxentries));

    if (data->num_compute % config->reset_iter == 0)
    {
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
            probevec.PutScalar(0.0); // TODO: Move it out
            for (int k = 0; k < nvectors; k++)
            {
                cindex = k+i;
                // TODO: Can do better than this, just need to go to the column map
                // of C, there might be null columns in C
                // Not much of use for Shasta 2x2 .. Later.
                probevec.ReplaceGlobalValue(g_rows[cindex], k, 1.0);
                //if (mypid == 0)
                //cout << "Changing row to 1.0 " << g_rows[cindex] << endl;
            }

#ifdef TIMING_OUTPUT
            app_time.start();
#endif
            probeop.Apply(probevec, Scol);
#ifdef TIMING_OUTPUT
            app_time.stop();
#endif

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
                        err = Sbar->InsertGlobalValues(g_rows[j], nentries, values,
                                                        indices);
                        assert(err >= 0);
                        err = lSGraph->InsertGlobalIndices(g_rows[j], nentries,
                                                        indices);
                        assert(err >= 0);
                    }
                    else if (abs(vecvalues[j]/maxvalue[k]) > relative_thres)
                    {
                        values[nentries] = vecvalues[j];
                        indices[nentries] = g_rows[cindex];
                        nentries++;
                        err = Sbar->InsertGlobalValues(g_rows[j], nentries, values,
                                                        indices);
                        assert(err >= 0);
                        err = lSGraph->InsertGlobalIndices(g_rows[j], nentries,
                                                        indices);
                        assert(err >= 0);
                    }
                    else
                    {
                        if (vecvalues[j] != 0.0)
                        {
                            dropped++;
                            //cout << "vecvalues[j]" << vecvalues[j] <<
                                    // " max" << maxvalue[k] << endl;
                        }
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
                    err = Sbar->InsertGlobalValues(g_rows[j], nentries, values, indices);
                    assert(err >= 0);
                    err = lSGraph->InsertGlobalIndices(g_rows[j], nentries, indices);
                    assert(err >= 0);
                }
                else if (abs(vecvalues[j]/maxvalue[0]) > relative_thres)
                {
                    values[nentries] = vecvalues[j];
                    indices[nentries] = g_rows[i];
                    nentries++;
                    err = Sbar->InsertGlobalValues(g_rows[j], nentries, values, indices);
                    assert(err >= 0);
                    err = lSGraph->InsertGlobalIndices(g_rows[j], nentries, indices);
                    assert(err >= 0);
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
        probeop.PrintTimingInfo();
        Sbar->FillComplete();
        lSGraph->FillComplete();

        data->localSbargraph = lSGraph;

#ifdef DUMP_MATRICES
        Epetra_Map defMap2(-1, g_localElems, 0, C->Comm());
        EpetraExt::ViewTransform<Epetra_CrsMatrix> * ReIdx_MatTrans2 =
                            new EpetraExt::CrsMatrix_Reindex( defMap2 );
        Epetra_CrsMatrix t2S = (*ReIdx_MatTrans2)( *Sbar );
        ReIdx_MatTrans2->fwd();
        EpetraExt::RowMatrixToMatlabFile("Schur.mat", t2S);
#endif

        cout << "#dropped entries" << dropped << endl;
        delete[] values;
        delete[] indices;
        delete[] maxvalue;
    }
    else
    {
        if (((data->num_compute-1) % config->reset_iter) == 0)
        {
            // We recomputed the Schur complement with dropping for the last
            // compute. Reset the prober with the new orthogonal vectors for
            // the Sbar from the previous iteration.
            Teuchos::ParameterList pList;
            Teuchos::RCP<Isorropia::Epetra::Prober> gprober =
                         Teuchos::RCP<Isorropia::Epetra::Prober> (new
                          Isorropia::Epetra::Prober(
                            data->localSbargraph.getRawPtr(), pList, false));
            gprober->color();
            data->guided_prober = gprober;

        }
        // Use the prober to probe the probeop for the sparsity pattern
        // add that to Sbar and call Fill complete
        int nvectors = data->guided_prober->getNumOrthogonalVectors();
        cout << "Number of Orthogonal Vectors for guided probing" << nvectors
                << endl;

        probeop.ResetTempVectors(nvectors);
        Teuchos::RCP<Epetra_CrsMatrix> blockdiag_Sbar =
                                 data->guided_prober->probe(probeop);
        int maxentries = blockdiag_Sbar->GlobalMaxNumEntries();
        int *indices = new int[maxentries];
        double *values = new double[maxentries];

        int numentries;
        for (int i = 0; i < blockdiag_Sbar->NumGlobalRows() ; i++)
        {
            int gid = blockdiag_Sbar->GRID(i);
            blockdiag_Sbar->ExtractGlobalRowCopy(gid, maxentries, numentries,
                                            values, indices);
            Sbar->InsertGlobalValues(gid, numentries, values, indices);
        }

        Sbar->FillComplete();
        delete[] indices;
        delete[] values;
    }

    delete[] values1;
    delete[] indices1;
    delete[] values2;
    delete[] indices2;
    delete[] values3;
    delete[] indices3;
    return Sbar;
}
