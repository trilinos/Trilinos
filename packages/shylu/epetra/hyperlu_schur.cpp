
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

    HyperLU_Probing_Operator probeop(G, R, LP, solver, C, localDRowMap);

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
    int nvectors = 16;
    double *maxvalue = new double[nvectors];
#ifdef TIMING_OUTPUT
    ftime.start();
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

        probeop.Apply(probevec, Scol);
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

    for ( ; i < totalElems ; i++)
    {
        Epetra_MultiVector probevec(rMap, 1); // TODO: Try doing more than one
        Epetra_MultiVector Scol(rMap, 1);     // vector at a time

        probevec.PutScalar(0.0);
        if (i >= prefixSum - localElems && i < prefixSum)
        {
            probevec.ReplaceGlobalValue(allSGID[i], 0, 1.0);
        }

        probeop.Apply(probevec, Scol);
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
    Sbar->FillComplete();
    cout << "#dropped entries" << dropped << endl;
    delete[] allSGID;
    delete[] mySGID;
    delete[] values;
    delete[] indices;
    delete[] maxvalue;

    return Sbar;
}
