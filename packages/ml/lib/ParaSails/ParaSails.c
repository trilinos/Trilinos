/*BHEADER**********************************************************************
 * (c) 1999   The Regents of the University of California
 *
 * See the file COPYRIGHT_and_DISCLAIMER for a complete copyright
 * notice, contact person, and disclaimer.
 *
 * $Revision$
 *********************************************************************EHEADER*/
/******************************************************************************
 *
 * ParaSails - Parallel sparse approximate inverse least squares.
 *
 *****************************************************************************/

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "Common.h"
#include "Matrix.h"
#include "Numbering.h"
#include "RowPatt.h"
#include "StoredRows.h"
#include "PrunedRows.h"
#include "OrderStat.h"
#include "LoadBal.h"
#include "ParaSails.h"

#define ROW_REQ_TAG        222
#define ROW_REPI_TAG       223
#define ROW_REPV_TAG       224

#ifdef ESSL
#include <essl.h>
#else
void hypre_F90_NAME(dpotrf)(char *, int *, double *, int *, int *);
void hypre_F90_NAME(dpotrs)(char *, int *, int *, double *, int *, double *,
  int *, int *);
void hypre_F90_NAME(dgels)(char *, int *, int *, int *, double *, int *,
  double *, int *, double *, int *, int *);
#endif

/******************************************************************************
 *
 * ParaSails private functions
 *
 *****************************************************************************/

/*--------------------------------------------------------------------------
 * FindNumReplies - Find the number of replies that this processor should
 * expect.  The input "replies_list" is an array that indicates what
 * processors were sent a message from the local processor.  An Allreduce
 * operation determines the total number of messages sent to the local
 * processor.
 *--------------------------------------------------------------------------*/

int FindNumReplies(MPI_Comm comm, int *replies_list)
{
    int num_replies;
    int npes, mype;
    int *replies_list2;

    MPI_Comm_rank(comm, &mype);
    MPI_Comm_size(comm, &npes);

    replies_list2 = (int *) malloc(npes * sizeof(int));

    MPI_Allreduce(replies_list, replies_list2, npes, MPI_INT, MPI_SUM, comm);
    num_replies = replies_list2[mype];

    free(replies_list2);

    return num_replies;
}

/*--------------------------------------------------------------------------
 * SendRequests - Given a list of indices "reqind" of length "reqlen",
 * send a sublist to the appropriate processors, thereby requesting
 * the rows (for example) corresponding to these indices.  The number
 * of requests made is returned in "num_requests".
 *
 * comm   - MPI communicator (input)
 * mat    - matrix used to map row and column numbers to processors (input)
 * reqlen - length of request list (input)
 * reqind - list of indices (input)
 * num_requests - number of requests made (output)
 * replies_list - if non-null, on input this should be a buffer initialized
 *          to zero of size the number of nonzero entries.  On output this
 *          buffer contains a 1 in position i if a request was made to
 *          processor i.  This array can be used to count (using
 *          MPI_AllReduce) the number of requests made to the current
 *          processor when the communication pattern is nonsymmetric.
 *--------------------------------------------------------------------------*/

static void SendRequests(MPI_Comm comm, Matrix *mat, int reqlen, int *reqind,
  int *num_requests, int *replies_list)
{
    MPI_Request request;
    int i, j, this_pe;

    shell_sort(reqlen, reqind);

    *num_requests = 0;

    for (i=0; i<reqlen; i=j) /* j is set below */
    {
        /* The processor that owns the row with index reqind[i] */
        this_pe = MatrixRowPe(mat, reqind[i]);

        /* Figure out other rows we need from this_pe */
        for (j=i+1; j<reqlen; j++)
        {
            /* if row is on different pe */
            if (reqind[j] < mat->beg_rows[this_pe] ||
                reqind[j] > mat->end_rows[this_pe])
                   break;
        }

        /* Request rows in reqind[i..j-1] */
        MPI_Isend(&reqind[i], j-i, MPI_INT, this_pe, ROW_REQ_TAG,
            comm, &request);
        MPI_Request_free(&request);
        (*num_requests)++;

        if (replies_list != NULL)
            replies_list[this_pe] = 1;
    }
}

/*--------------------------------------------------------------------------
 * ReceiveRequest - Receive a request sent with SendRequests by another
 * processor.  This function should be placed inside a loop which is
 * executed once for every request that this processor expects to receive.
 * This is the number of requests this processor made in SendRequests
 * in the symmetric case.
 *
 * comm   - MPI communicator (input)
 * source - number of the processor that sent the message (output)
 * buffer - buffer provided by the user.  On output, it contains a
 *          list of indices.  Buffer will be reallocated if too small
 *          (input/output)
 * buflen - size of the buffer (input).  Size will be updated if buffer
 *          is too small (input/output)
 * count  - number of indices in the output buffer (output)
 *--------------------------------------------------------------------------*/

static void ReceiveRequest(MPI_Comm comm, int *source, int **buffer,
  int *buflen, int *count)
{
    MPI_Status status;

    MPI_Probe(MPI_ANY_SOURCE, ROW_REQ_TAG, comm, &status);
    *source = status.MPI_SOURCE;
    MPI_Get_count(&status, MPI_INT, count);

    if (*count > *buflen)
    {
        free(*buffer);
        *buflen = *count;
        *buffer = (int *) malloc(*buflen * sizeof(int));
    }

    MPI_Recv(*buffer, *count, MPI_INT, *source, ROW_REQ_TAG, comm, &status);
}

/*--------------------------------------------------------------------------
 * SendReplyPrunedRows - Send a reply of pruned rows for each request
 * received by this processor using ReceiveRequest.
 *
 * comm    - MPI communicator (input)
 * dest    - pe to send to (input)
 * buffer  - list of indices (input)
 * count   - number of indices in buffer (input)
 * pruned_rows - the pruned_rows object where the pruned rows reside (input)
 * mem     - pointer to memory object used for reply buffers (input)
 * request - request handle of send (output)
 *
 * The function allocates space for each send buffer using "mem", and the
 * caller must free this space when all the sends have completed..
 *
 * The reply has the following structure for the integer data in indbuf:
 * num_rows, index_1, ..., index_n, len_1, row_1_indices, len_2, indices, ...
 *--------------------------------------------------------------------------*/

static void SendReplyPrunedRows(MPI_Comm comm, Numbering *numb,
  int dest, int *buffer, int count,
  PrunedRows *pruned_rows, Mem *mem, MPI_Request *request)
{
    int sendbacksize, j;
    int len, *ind, *indbuf, *indbufp;
    int temp;

    /* Determine the size of the integer message we need to send back */
    sendbacksize = count+1; /* length of header part */
    for (j=0; j<count; j++)
    {
        NumberingGlobalToLocal(numb, 1, &buffer[j], &temp);
        PrunedRowsGet(pruned_rows, temp, &len, &ind);
        sendbacksize += (len+1);  /* add one for the row length */
    }

    /* Reply buffer - will be freed by caller */
    indbuf = (int *) MemAlloc(mem, sendbacksize * sizeof(int));

    /* Pointer used to construct reply message */
    indbufp = indbuf;

    /* Construct integer reply message in local buffer, with this format:
       number of rows to send, row numbers, indices of each row */

    *indbufp++ = count; /* number of rows to send */

    for (j=0; j<count; j++)
        *indbufp++ = buffer[j]; /* row numbers */

    for (j=0; j<count; j++)
    {
        NumberingGlobalToLocal(numb, 1, &buffer[j], &temp);
        PrunedRowsGet(pruned_rows, temp, &len, &ind);

        *indbufp++ = len;
        /* memcpy(indbufp, ind, sizeof(int)*len); */
        NumberingLocalToGlobal(numb, len, ind, indbufp);
        indbufp += len;
    }

    MPI_Isend(indbuf, indbufp-indbuf, MPI_INT, dest, ROW_REPI_TAG,
        comm, request);
}

/*--------------------------------------------------------------------------
 * ReceiveReplyPrunedRows - Receive a reply sent by SendReplyPrunedRows
 *
 * comm    - MPI communicator (input)
 * pruned_rows - the pruned_rows object where the rows should be stored
 * patt    - each pruned row is merged into patt before returning (input).
 *           Only the external indices of the pattern is merged
 * mat     - Matrix argument used for determining the external indices
 *--------------------------------------------------------------------------*/

static void ReceiveReplyPrunedRows(MPI_Comm comm, Numbering *numb,
  PrunedRows *pruned_rows, RowPatt *patt)
{
    MPI_Status status;
    int source, count;
    int len, *ind, num_rows, *row_nums, j;

    /* Don't know the size of reply, so use probe and get count */
    MPI_Probe(MPI_ANY_SOURCE, ROW_REPI_TAG, comm, &status);
    source = status.MPI_SOURCE;
    MPI_Get_count(&status, MPI_INT, &count);

    /* Allocate space in stored rows data structure */
    ind = PrunedRowsAlloc(pruned_rows, count);
    MPI_Recv(ind, count, MPI_INT, source, ROW_REPI_TAG, comm, &status);

    /* Parse the message */
    num_rows = *ind++; /* number of rows */
    row_nums = ind;    /* row numbers */
    ind += num_rows;

    /* Convert global row numbers to local row numbers */
    NumberingGlobalToLocal(numb, num_rows, row_nums, row_nums);

    /* Set the pointers to the individual rows */
    for (j=0; j<num_rows; j++)
    {
        len = *ind++;
        NumberingGlobalToLocal(numb, len, ind, ind);
        PrunedRowsPut(pruned_rows, row_nums[j], len, ind);
        RowPattMergeExt(patt, len, ind, numb->num_loc);
        ind += len;
    }
}

/*--------------------------------------------------------------------------
 * SendReplyStoredRows - Send a reply of stored rows for each request
 * received by this processor using ReceiveRequest.
 *
 * comm    - MPI communicator (input)
 * dest    - pe to send to (input)
 * buffer  - list of indices (input)
 * count   - number of indices in buffer (input)
 * stored_rows - the stored_rows object where the rows reside (input)
 * mem     - pointer to memory object used for reply buffers (input)
 * request - request handle of send (output)
 *
 * The function allocates space for each send buffer using "mem", and the
 * caller must free this space when all the sends have completed..
 *
 * The reply has the following structure for the integer data in indbuf:
 * num_rows, index_1, ..., index_n, len_1, row_1_indices, len_2, indices, ...
 *
 * The reply has the following structure for the value data:
 * row_1_values, row_2_values, ...
 *--------------------------------------------------------------------------*/

static void SendReplyStoredRows(MPI_Comm comm, Numbering *numb,
  int dest, int *buffer, int count,
  StoredRows *stored_rows, Mem *mem, MPI_Request *request)
{
    int sendbacksize, j;
    int len, *ind, *indbuf, *indbufp;
    double *val, *valbuf, *valbufp;
    int temp;

    /* Determine the size of the integer message we need to send back */
    sendbacksize = count+1; /* length of header part */
    for (j=0; j<count; j++)
    {
        NumberingGlobalToLocal(numb, 1, &buffer[j], &temp);
        StoredRowsGet(stored_rows, temp, &len, &ind, &val);
        sendbacksize += (len+1);  /* add one for the row length */
    }

    /* Reply buffers - will be freed by caller */
    indbuf = (int *)    MemAlloc(mem, sendbacksize * sizeof(int));
    valbuf = (double *) MemAlloc(mem, sendbacksize * sizeof(double));

    /* Pointers used to construct reply messages */
    indbufp = indbuf;
    valbufp = valbuf;

    /* Construct integer reply message in local buffer, with this format:
       number of rows to send, row numbers, len of row, indices each row,
       len of next row, indices of row, etc. */

    *indbufp++ = count; /* number of rows to send */

    for (j=0; j<count; j++)
        *indbufp++ = buffer[j]; /* row numbers */

    for (j=0; j<count; j++)
    {
        NumberingGlobalToLocal(numb, 1, &buffer[j], &temp);
        StoredRowsGet(stored_rows, temp, &len, &ind, &val);

        *indbufp++ = len;
        /* memcpy(indbufp, ind, sizeof(int)*len); */
        NumberingLocalToGlobal(numb, len, ind, indbufp);
        memcpy(valbufp, val, sizeof(double)*len);
        indbufp += len;
        valbufp += len;
    }

    MPI_Isend(indbuf, indbufp-indbuf, MPI_INT, dest, ROW_REPI_TAG,
        comm, request);

    MPI_Request_free(request);

    MPI_Isend(valbuf, valbufp-valbuf, MPI_DOUBLE, dest, ROW_REPV_TAG,
        comm, request);
}

/*--------------------------------------------------------------------------
 * ReceiveReplyStoredRows - Receive a reply sent by SendReplyStoredRows
 *
 * comm    - MPI communicator (input)
 * numb    - Numbering object (input)
 * stored_rows - the stored_rows object where the rows should be stored
 *--------------------------------------------------------------------------*/

static void ReceiveReplyStoredRows(MPI_Comm comm, Numbering *numb,
  StoredRows *stored_rows)
{
    MPI_Status status;
    int source, count;
    int len, *ind, num_rows, *row_nums, j;
    double *val;

    /* Don't know the size of reply, so use probe and get count */
    MPI_Probe(MPI_ANY_SOURCE, ROW_REPI_TAG, comm, &status);
    source = status.MPI_SOURCE;
    MPI_Get_count(&status, MPI_INT, &count);

    /* Allocate space in stored rows data structure */
    ind = StoredRowsAllocInd(stored_rows, count);
    MPI_Recv(ind, count, MPI_INT, source, ROW_REPI_TAG, comm, &status);
    val = StoredRowsAllocVal(stored_rows, count);
    MPI_Recv(val, count, MPI_DOUBLE, source, ROW_REPV_TAG, comm, &status);

    /* Parse the message */
    num_rows = *ind++; /* number of rows */
    row_nums = ind;    /* row numbers */
    ind += num_rows;

    /* Convert global row numbers to local row numbers */
    NumberingGlobalToLocal(numb, num_rows, row_nums, row_nums);

    /* Set the pointers to the individual rows */
    for (j=0; j<num_rows; j++)
    {
        len = *ind++;
        NumberingGlobalToLocal(numb, len, ind, ind);
        StoredRowsPut(stored_rows, row_nums[j], len, ind, val);
        ind += len;
        val += len;
    }
}

/*--------------------------------------------------------------------------
 * ExchangePrunedRows
 *--------------------------------------------------------------------------*/

static void ExchangePrunedRows(MPI_Comm comm, Matrix *M, Numbering *numb,
  PrunedRows *pruned_rows, int num_levels)
{
    RowPatt *patt;
    int row, len, *ind;

    int num_requests;
    int source;

    int bufferlen;
    int *buffer;

    int level;

    int i;
    int count;
    MPI_Request *requests;
    MPI_Status *statuses;
    int npes;
    int num_replies, *replies_list;

    Mem *mem;

    MPI_Comm_size(comm, &npes);
    requests = (MPI_Request *) malloc(npes * sizeof(MPI_Request));
    statuses = (MPI_Status *) malloc(npes * sizeof(MPI_Status));

    /* Merged pattern of pruned rows on this processor */

    patt = RowPattCreate(PARASAILS_MAXLEN);

    for (row=0; row<=M->end_row - M->beg_row; row++)
    {
        PrunedRowsGet(pruned_rows, row, &len, &ind);
        RowPattMergeExt(patt, len, ind, numb->num_loc);
    }

    /* Loop to construct pattern of pruned rows on this processor */

    bufferlen = 10; /* size will grow if get a long msg */
    buffer = (int *) malloc(bufferlen * sizeof(int));

    for (level=1; level<=num_levels; level++)
    {
        mem = (Mem *) MemCreate();

        /* Get list of indices that were just merged */
        RowPattPrevLevel(patt, &len, &ind);

        /* Convert local row numbers to global row numbers */
        NumberingLocalToGlobal(numb, len, ind, ind);

        replies_list = (int *) calloc(npes, sizeof(int));

        SendRequests(comm, M, len, ind, &num_requests, replies_list);

        num_replies = FindNumReplies(comm, replies_list);
        free(replies_list);

        for (i=0; i<num_replies; i++)
        {
            /* Receive count indices stored in buffer */
            ReceiveRequest(comm, &source, &buffer, &bufferlen, &count);

            SendReplyPrunedRows(comm, numb, source, buffer, count,
                pruned_rows, mem, &requests[i]);
        }

        for (i=0; i<num_requests; i++)
        {
            /* Will also merge the pattern of received rows into "patt" */
            ReceiveReplyPrunedRows(comm, numb, pruned_rows, patt);
        }

        MPI_Waitall(num_replies, requests, statuses);
        MemDestroy(mem);
    }

    RowPattDestroy(patt);
    free(buffer);
    free(requests);
    free(statuses);
}

/*--------------------------------------------------------------------------
 * ExchangePrunedRowsExt
 *--------------------------------------------------------------------------*/

static void ExchangePrunedRowsExt(MPI_Comm comm, Matrix *M, Numbering *numb,
  PrunedRows *pruned_rows_global, PrunedRows *pruned_rows_local, int num_levels)
{
    RowPatt *patt;
    int row, len, *ind;

    int num_requests;
    int source;

    int bufferlen;
    int *buffer;

    int level;

    int i;
    int count;
    MPI_Request *requests;
    MPI_Status *statuses;
    int npes;
    int num_replies, *replies_list;

    Mem *mem;

    MPI_Comm_size(comm, &npes);
    requests = (MPI_Request *) malloc(npes * sizeof(MPI_Request));
    statuses = (MPI_Status *) malloc(npes * sizeof(MPI_Status));

    /* Merged pattern of pruned rows on this processor */

    patt = RowPattCreate(PARASAILS_MAXLEN);

    for (row=0; row<=M->end_row - M->beg_row; row++)
    {
        PrunedRowsGet(pruned_rows_global, row, &len, &ind);
        RowPattMergeExt(patt, len, ind, numb->num_loc);
    }

    /* Loop to construct pattern of pruned rows on this processor */

    bufferlen = 10; /* size will grow if get a long msg */
    buffer = (int *) malloc(bufferlen * sizeof(int));

    for (level=0; level<=num_levels; level++)  /* MUST DO THIS AT LEAST ONCE */
    {
        mem = (Mem *) MemCreate();

        /* Get list of indices that were just merged */
        RowPattPrevLevel(patt, &len, &ind);

        /* Convert local row numbers to global row numbers */
        NumberingLocalToGlobal(numb, len, ind, ind);

        replies_list = (int *) calloc(npes, sizeof(int));

        SendRequests(comm, M, len, ind, &num_requests, replies_list);

        num_replies = FindNumReplies(comm, replies_list);
        free(replies_list);

        for (i=0; i<num_replies; i++)
        {
            /* Receive count indices stored in buffer */
            ReceiveRequest(comm, &source, &buffer, &bufferlen, &count);

            SendReplyPrunedRows(comm, numb, source, buffer, count,
                pruned_rows_local, mem, &requests[i]);
        }

        for (i=0; i<num_requests; i++)
        {
            /* Will also merge the pattern of received rows into "patt" */
            ReceiveReplyPrunedRows(comm, numb, pruned_rows_local, patt);
        }

        MPI_Waitall(num_replies, requests, statuses);
        MemDestroy(mem);
    }

    RowPattDestroy(patt);
    free(buffer);
    free(requests);
    free(statuses);
}

/*--------------------------------------------------------------------------
 * ExchangePrunedRowsExt2 - part 2 of the algorithm
 *--------------------------------------------------------------------------*/

static void ExchangePrunedRowsExt2(MPI_Comm comm, Matrix *M, Numbering *numb,
  PrunedRows *pruned_rows_global, PrunedRows *pruned_rows_local, int num_levels)
{
    RowPatt *patt;
    int row, len, *ind;

    int num_requests;
    int source;

    int bufferlen;
    int *buffer;

    int level;

    int i;
    int count;
    MPI_Request *requests;
    MPI_Status *statuses;
    int npes;
    int num_replies, *replies_list;

    Mem *mem;

    MPI_Comm_size(comm, &npes);
    requests = (MPI_Request *) malloc(npes * sizeof(MPI_Request));
    statuses = (MPI_Status *) malloc(npes * sizeof(MPI_Status));

    /* Merged pattern of pruned rows on this processor */

    patt = RowPattCreate(PARASAILS_MAXLEN);

    for (row=0; row<=M->end_row - M->beg_row; row++)
    {
        PrunedRowsGet(pruned_rows_local, row, &len, &ind);
        RowPattMergeExt(patt, len, ind, numb->num_loc);
    }

    /* Compute powers with local matrix - no communication is needed */

    for (level=1; level<=num_levels; level++)
    {
        int lenprev, *indprev;

        /* Get the indices that were just added */
        RowPattPrevLevel(patt, &lenprev, &indprev);

        for (i=0; i<lenprev; i++)
        {
            PrunedRowsGet(pruned_rows_local, indprev[i], &len, &ind);
            RowPattMergeExt(patt, len, ind, numb->num_loc);
        }
    }

    /* Now get rows from pruned_rows_global */

    bufferlen = 10; /* size will grow if get a long msg */
    buffer = (int *) malloc(bufferlen * sizeof(int));

    /* DO THIS ONCE */
    {
        mem = (Mem *) MemCreate();

	/* Get list of indices - these are all nonlocal indices */
        RowPattGet(patt, &len, &ind);

        /* Convert local row numbers to global row numbers */
        NumberingLocalToGlobal(numb, len, ind, ind);

        replies_list = (int *) calloc(npes, sizeof(int));

        SendRequests(comm, M, len, ind, &num_requests, replies_list);

        num_replies = FindNumReplies(comm, replies_list);
        free(replies_list);

        for (i=0; i<num_replies; i++)
        {
            /* Receive count indices stored in buffer */
            ReceiveRequest(comm, &source, &buffer, &bufferlen, &count);

            SendReplyPrunedRows(comm, numb, source, buffer, count,
                pruned_rows_global, mem, &requests[i]);
        }

        for (i=0; i<num_requests; i++)
        {
            /* Will also merge the pattern of received rows into "patt" */
            ReceiveReplyPrunedRows(comm, numb, pruned_rows_global, patt);
        }

        MPI_Waitall(num_replies, requests, statuses);
        MemDestroy(mem);
    }

    RowPattDestroy(patt);
    free(buffer);
    free(requests);
    free(statuses);
}

/*--------------------------------------------------------------------------
 * ExchangeStoredRows
 *--------------------------------------------------------------------------*/

static void ExchangeStoredRows(MPI_Comm comm, Matrix *A, Matrix *M,
  Numbering *numb, StoredRows *stored_rows, LoadBal *load_bal)
{
    RowPatt *patt;
    int row, len, *ind;
    double *val;

    int num_requests;
    int source;

    int bufferlen;
    int *buffer;

    int i;
    int count;
    MPI_Request *requests;
    MPI_Status *statuses;
    int npes;
    int num_replies, *replies_list;

    Mem *mem = (Mem *) MemCreate();

    MPI_Comm_size(comm, &npes);

    /* Merge the patterns of all the rows of M on this processor */
    /* The merged pattern is not already known, since M is triangular */

    patt = RowPattCreate(PARASAILS_MAXLEN);

    /* for (row=load_bal->beg_row; row<=M->end_row; row++) */
    /* We need the additional rows if we need to Rescale */
    /* i.e., if filter is nonzero and we are in symmetric case */

    for (row=M->beg_row; row<=M->end_row; row++)
    {
        MatrixGetRow(M, row - M->beg_row, &len, &ind, &val);
        RowPattMergeExt(patt, len, ind, numb->num_loc);
    }

    /* Merge patterns for load balancing recipient rows */

    for (i=0; i<load_bal->num_taken; i++)
    {
      for (row=0; row <= load_bal->recip_data[i].mat->end_row -
                         load_bal->recip_data[i].mat->beg_row; row++)
      {
        MatrixGetRow(load_bal->recip_data[i].mat, row, &len, &ind, &val);
        RowPattMergeExt(patt, len, ind, numb->num_loc);
      }
    }

    RowPattGet(patt, &len, &ind);

    /* Convert local row numbers to global row numbers */
    NumberingLocalToGlobal(numb, len, ind, ind);

    replies_list = (int *) calloc(npes, sizeof(int));

    SendRequests(comm, A, len, ind, &num_requests, replies_list);

    num_replies = FindNumReplies(comm, replies_list);
    free(replies_list);

    requests = (MPI_Request *) malloc(num_replies * sizeof(MPI_Request));
    statuses = (MPI_Status *) malloc(num_replies * sizeof(MPI_Status));

    bufferlen = 10; /* size will grow if get a long msg */
    buffer = (int *) malloc(bufferlen * sizeof(int));

    for (i=0; i<num_replies; i++)
    {
        /* Receive count indices stored in buffer */
        ReceiveRequest(comm, &source, &buffer, &bufferlen, &count);

        SendReplyStoredRows(comm, numb, source, buffer, count,
            stored_rows, mem, &requests[i]);
    }

    for (i=0; i<num_requests; i++)
    {
        ReceiveReplyStoredRows(comm, numb, stored_rows);
    }

    MPI_Waitall(num_replies, requests, statuses);

    /* Free all send buffers */
    MemDestroy(mem);

    RowPattDestroy(patt);
    free(buffer);
    free(requests);
    free(statuses);
}

/*--------------------------------------------------------------------------
 * ConstructPatternForEachRow
 *
 * pruned_rows - pruned rows, used for constructing row patterns (input)
 * num_levels  - number of levels in pattern (input)
 * M           - matrix where the row patterns will be stored (input/output).
 *               This is the approximate inverse with lower triangular pattern
 *--------------------------------------------------------------------------*/

static void ConstructPatternForEachRow(int symmetric, PrunedRows *pruned_rows,
  int num_levels, Numbering *numb, Matrix *M, double *costp)
{
    int row, len, *ind, level, lenprev, *indprev;
    int i, j;
    RowPatt *row_patt;
    int nnz = 0;
    int npes;

    MPI_Comm_size(M->comm, &npes);
    *costp = 0.0;

    row_patt = RowPattCreate(PARASAILS_MAXLEN);

    for (row=0; row<=M->end_row - M->beg_row; row++)
    {
        /* Get initial pattern for row */
        PrunedRowsGet(pruned_rows, row, &len, &ind);
        RowPattMerge(row_patt, len, ind);

        /* Loop */
        for (level=1; level<=num_levels; level++)
        {
            /* Get the indices that were just added */
            RowPattPrevLevel(row_patt, &lenprev, &indprev);

            for (i=0; i<lenprev; i++)
            {
                PrunedRowsGet(pruned_rows, indprev[i], &len, &ind);
                RowPattMerge(row_patt, len, ind);
            }
        }

        RowPattGet(row_patt, &len, &ind);

        /* do reset here, because now we mess with ind array */
        RowPattReset(row_patt);

        if (symmetric)
        {
            /* Store the lower triangular part of row pattern into the matrix */
            j = 0;
            for (i=0; i<len; i++)
            {
                if (numb->local_to_global[ind[i]] <= numb->local_to_global[row])
                    ind[j++] = ind[i];
            }
            len = j;
        }

        /* Store structure of row in matrix M */
        /* Following statement allocates space but does not store values */
        MatrixSetRow(M, row+M->beg_row, len, ind, NULL);

        nnz += len;
        (*costp) += (double) len*len*len;
    }

#if 0
    {
    int mype;
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    printf("%d: nnz: %10d  ********* cost %7.1e\n", mype, nnz, *costp);
    fflush(NULL);
    }
#endif

    RowPattDestroy(row_patt);
}

/*--------------------------------------------------------------------------
 * ConstructPatternForEachRowExt - extended version
 *
 * pruned_rows - pruned rows, used for constructing row patterns (input)
 * num_levels  - number of levels in pattern (input)
 * M           - matrix where the row patterns will be stored (input/output).
 *               This is the approximate inverse with lower triangular pattern
 *--------------------------------------------------------------------------*/

static void ConstructPatternForEachRowExt(int symmetric, 
  PrunedRows *pruned_rows_global, PrunedRows *pruned_rows_local, 
  int num_levels, Numbering *numb, Matrix *M, double *costp)
{
    int row, len, *ind, level, lenprev, *indprev;
    int i, j;
    RowPatt *row_patt;
    RowPatt *row_patt2;
    int nnz = 0;
    int npes;

    MPI_Comm_size(M->comm, &npes);
    *costp = 0.0;

    row_patt = RowPattCreate(PARASAILS_MAXLEN);
    row_patt2 = RowPattCreate(PARASAILS_MAXLEN);

    for (row=0; row<=M->end_row - M->beg_row; row++)
    {
        /* Get initial pattern for row */
        PrunedRowsGet(pruned_rows_global, row, &len, &ind);
        RowPattMerge(row_patt, len, ind);

        /* Loop */
        for (level=0; level<=num_levels; level++) /* at least once */
        {
            /* Get the indices that were just added */
            RowPattPrevLevel(row_patt, &lenprev, &indprev);

            for (i=0; i<lenprev; i++)
            {
                PrunedRowsGet(pruned_rows_local, indprev[i], &len, &ind);
                RowPattMerge(row_patt, len, ind);
            }
        }

        /***********************
	 * Now do the transpose 
	 ***********************/

        /* Get initial pattern for row */
        PrunedRowsGet(pruned_rows_local, row, &len, &ind);
        RowPattMerge(row_patt2, len, ind);

        /* Loop */
        for (level=1; level<=num_levels; level++)
        {
            /* Get the indices that were just added */
            RowPattPrevLevel(row_patt2, &lenprev, &indprev);

            for (i=0; i<lenprev; i++)
            {
                PrunedRowsGet(pruned_rows_local, indprev[i], &len, &ind);
                RowPattMerge(row_patt2, len, ind);
            }
        }

	/* One more merge, with pruned_rows_global */
        RowPattGet(row_patt2, &lenprev, &indprev);
        for (i=0; i<lenprev; i++)
        {
            PrunedRowsGet(pruned_rows_global, indprev[i], &len, &ind);
            RowPattMerge(row_patt2, len, ind);
        }


        /****************************
	 * Merge the two row patterns
	 ****************************/

        RowPattGet(row_patt2, &len, &ind);
        RowPattMerge(row_patt, len, ind);

        /****************************
	 * Done computing pattern!
	 ****************************/

	/* get the indices in the pattern */
        RowPattGet(row_patt, &len, &ind);

        /* do reset here, because now we mess with ind array */
        RowPattReset(row_patt);
        RowPattReset(row_patt2);

        if (symmetric)
        {
            /* Store the lower triangular part of row pattern into the matrix */
            j = 0;
            for (i=0; i<len; i++)
            {
                if (numb->local_to_global[ind[i]] <= numb->local_to_global[row])
                    ind[j++] = ind[i];
            }
            len = j;
        }

        /* Store structure of row in matrix M */
        /* Following statement allocates space but does not store values */
        MatrixSetRow(M, row+M->beg_row, len, ind, NULL);

        nnz += len;
        (*costp) += (double) len*len*len;
    }

#if 0
    {
    int mype;
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    printf("%d: nnz: %10d  ********* cost %7.1e\n", mype, nnz, *costp);
    fflush(NULL);
    }
#endif

    RowPattDestroy(row_patt);
    RowPattDestroy(row_patt2);
}

/*--------------------------------------------------------------------------
 * ComputeValuesSym
 *--------------------------------------------------------------------------*/

static void ComputeValuesSym(StoredRows *stored_rows, Matrix *mat,
  int local_beg_row, Numbering *numb, int symmetric)
{
    int *marker;
    int row, maxlen, len, *ind;
    double *val;

    double *ahat, *ahatp;
    int i, j, len2, *ind2, loc;
    double *val2, temp;
    double time0, time1, timet = 0.0, timea = 0.0;

    double ahatcost = 0.0;

#ifndef ESSL
    char uplo = 'L';
    int one = 1;
    int info;
#endif

    /* Allocate and initialize full length marker array */
    marker = (int *) malloc(numb->num_ind * sizeof(int));
    for (i=0; i<numb->num_ind; i++)
        marker[i] = -1;

    /* Determine the length of the longest row of M on this processor */
    /* This determines the maximum storage required for the ahat matrix */
    maxlen = 0;
    for (row=local_beg_row; row<=mat->end_row; row++)
    {
        MatrixGetRow(mat, row - mat->beg_row, &len, &ind, &val);
        maxlen = (len > maxlen ? len : maxlen);
    }

#ifdef ESSL
    ahat = (double *) malloc(maxlen*(maxlen+1)/2 * sizeof(double));
#else
    ahat = (double *) malloc(maxlen*maxlen * sizeof(double));
#endif

    /* Compute values for row "row" of approximate inverse */
    for (row=local_beg_row; row<=mat->end_row; row++)
    {
        /* Retrieve local indices */
        MatrixGetRow(mat, row - mat->beg_row, &len, &ind, &val);

        /* Fill marker array in locations of local indices */
        for (i=0; i<len; i++)
            marker[ind[i]] = i;

        /* Initialize ahat to zero */
#ifdef ESSL
        bzero((char *) ahat, len*(len+1)/2 * sizeof(double));
#else
        bzero((char *) ahat, len*len * sizeof(double));
#endif

        time0 = MPI_Wtime();

        /* Form ahat matrix, entries correspond to indices in "ind" only */
        ahatp = ahat;
        for (i=0; i<len; i++)
        {
            StoredRowsGet(stored_rows, ind[i], &len2, &ind2, &val2);
            assert(len2 > 0);

#ifdef ESSL
            for (j=0; j<len2; j++)
            {
                loc = marker[ind2[j]];

                if (loc != -1) /* redundant */
                    if (loc >= i)
                        ahatp[loc - i] = val2[j];
            }

            ahatp += (len-i);
#else
            for (j=0; j<len2; j++)
            {
                loc = marker[ind2[j]];

                if (loc != -1)
                    ahatp[loc] = val2[j];
            }

            ahatp += len;
#endif
        }

        if (symmetric == 2)
        {
#ifdef ESSL
            printf("Symmetric precon for nonsym problem not yet available\n");
            printf("for ESSL version.  Please contact the author.\n");
            PARASAILS_EXIT;
#else
            int k, kk;
            k = 0;
            for (i=0; i<len; i++)
            {
                for (j=0; j<len; j++)
                {
                    kk = j*len + i;
                    ahat[k] = (ahat[k] + ahat[kk]) / 2.0;
                    k++;
                }
            }
#endif
        }

        time1 = MPI_Wtime();
        timea += (time1-time0);
        ahatcost += (double) (len*len2);

        /* Set the right-hand side */
        bzero((char *) val, len*sizeof(double));
        NumberingGlobalToLocal(numb, 1, &row, &loc);
        loc = marker[loc];
        assert(loc != -1);
        val[loc] = 1.0;

        /* Reset marker array */
        for (i=0; i<len; i++)
            marker[ind[i]] = -1;

        time0 = MPI_Wtime();

#ifdef ESSL
        dppf(ahat, len, 1);
        dpps(ahat, len, val, 1);
#else
        /* Solve local linear system - factor phase */
        hypre_F90_NAME(dpotrf)(&uplo, &len, ahat, &len, &info);
        if (info != 0)
        {
            printf("Matrix may not be symmetric positive definite.\n");
            printf("ParaSails: row %d, dpotrf returned %d.\n", row, info);
            printf("ParaSails: len %d, ahat: %f %f %f %f\n", len,
                ahat[0], ahat[1], ahat[2], ahat[3]);
            PARASAILS_EXIT;
        }

        /* Solve local linear system - solve phase */
        hypre_F90_NAME(dpotrs)(&uplo, &len, &one, ahat, &len, val, &len, &info);
        if (info != 0)
        {
            printf("ParaSails: row %d, dpotrs returned %d.\n", row, info);
            printf("ParaSails: len %d, ahat: %f %f %f %f\n", len,
                ahat[0], ahat[1], ahat[2], ahat[3]);
            PARASAILS_EXIT;
        }
#endif
        time1 = MPI_Wtime();
        timet += (time1-time0);

        /* Scale the result */
        temp = 1.0 / sqrt(ABS(val[loc]));
        for (i=0; i<len; i++)
            val[i] = val[i] * temp;
    }

    free(marker);
    free(ahat);

#if 0
    {
    int mype;
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    printf("%d: Time for ahat: %f, for local solves: %f\n", mype, timea, timet);
    printf("%d: ahatcost: %7.1e, numrows: %d, maxlen: %d\n",
        mype, ahatcost, mat->end_row-local_beg_row+1, maxlen);
    fflush(NULL);
    }
#endif
}

/*--------------------------------------------------------------------------
 * ComputeValuesNonsym
 *--------------------------------------------------------------------------*/

static void ComputeValuesNonsym(StoredRows *stored_rows, Matrix *mat,
  int local_beg_row, Numbering *numb)
{
    int *marker;
    double *ahat, *ahatp, *bhat;
    double *work;
    int ahat_size = 10000, bhat_size = 1000, work_size = 2000*64;

    int row, len, *ind;
    double *val;

    int i, j, len2, *ind2, loc;
    double *val2;
    double time0, time1, timet = 0.0, timea = 0.0;

    int npat;
    int pattsize = 1000;
    int *patt = (int *) malloc(pattsize*sizeof(int));

    int info;

#ifndef ESSL
    char trans = 'N';
    int one = 1;
#endif

    /* Allocate and initialize marker array */
    /* Since numb already knows about the indices of the external rows that
       will be needed, numb_ind is the maximum size of the marker array */
    marker = (int *) malloc(numb->num_ind * sizeof(int));
    for (i=0; i<numb->num_ind; i++)
        marker[i] = -1;

    bhat = (double *) malloc(bhat_size * sizeof(double));
    ahat = (double *) malloc(ahat_size * sizeof(double));
    work = (double *) calloc(work_size,  sizeof(double));

    /* Compute values for row "row" of approximate inverse */
    for (row=local_beg_row; row<=mat->end_row; row++)
    {
        time0 = MPI_Wtime();

        /* Retrieve local indices */
        MatrixGetRow(mat, row - mat->beg_row, &len, &ind, &val);

        npat = 0;

        /* Put the diagonal entry into the marker array */
        NumberingGlobalToLocal(numb, 1, &row, &loc);
        marker[loc] = npat;
        patt[npat++] = loc;

        /* Fill marker array */
        for (i=0; i<len; i++)
        {
            StoredRowsGet(stored_rows, ind[i], &len2, &ind2, &val2);
            assert(len2 > 0);

            for (j=0; j<len2; j++)
            {
                loc = marker[ind2[j]];

                if (loc == -1)
                {
                    marker[ind2[j]] = npat;
                    if (npat >= pattsize)
                    {
                        pattsize = npat*2;
                        patt = (int *) realloc(patt, pattsize*sizeof(int));
                    }
                    patt[npat++] = ind2[j];
                }
            }
        }

        if (len*npat > ahat_size)
        {
            free(ahat);
            ahat_size = len*npat;
            ahat = (double *) malloc(ahat_size * sizeof(double));
        }

        /* Initialize ahat to zero */
        bzero((char *) ahat, len*npat * sizeof(double));

        /* Form ahat matrix, entries correspond to indices in "ind" only */
        ahatp = ahat;
        for (i=0; i<len; i++)
        {
            StoredRowsGet(stored_rows, ind[i], &len2, &ind2, &val2);

            for (j=0; j<len2; j++)
            {
                loc = marker[ind2[j]];
                ahatp[loc] = val2[j];
            }
            ahatp += npat;
        }

        time1 = MPI_Wtime();
        timea += (time1-time0);

        /* Reallocate bhat if necessary */
        if (npat > bhat_size)
        {
            free(bhat);
            bhat_size = npat;
            bhat = (double *) malloc(bhat_size * sizeof(double));
        }

        /* Set the right-hand side, bhat */
        bzero((char *) bhat, npat*sizeof(double));
        NumberingGlobalToLocal(numb, 1, &row, &loc);
        loc = marker[loc];
        assert(loc != -1);
        bhat[loc] = 1.0;

        /* Reset marker array */
        for (i=0; i<npat; i++)
            marker[patt[i]] = -1;

        time0 = MPI_Wtime();

#ifdef ESSL
        /* rhs in bhat, and put solution in val */
        dgells(0, ahat, npat, bhat, npat, val, len, NULL, 1.e-12, npat, len, 1,
            &info, work, work_size);
#else
        /* rhs in bhat, and put solution in bhat */
        hypre_F90_NAME(dgels)(&trans, &npat, &len, &one, ahat, &npat,
            bhat, &npat, work, &work_size, &info);

        if (info != 0)
        {
            printf("ParaSails: row %d, dgels returned %d.\n", row, info);
            printf("ParaSails: len %d, ahat: %f %f %f %f\n", len,
                ahat[0], ahat[1], ahat[2], ahat[3]);
            PARASAILS_EXIT;
        }

        /* Copy result into row */
        for (j=0; j<len; j++)
            val[j] = bhat[j];
#endif
        time1 = MPI_Wtime();
        timet += (time1-time0);
    }

    free(patt);
    free(marker);
    free(bhat);
    free(ahat);
    free(work);

#if 0
    {
    int mype;
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    printf("%d: Time for ahat: %f, for local solves: %f\n", mype, timea, timet);
    fflush(NULL);
    }
#endif
}

/*--------------------------------------------------------------------------
 * SelectThresh - select a threshold for the preconditioner pattern.
 * The threshold attempts to be chosen such that approximately (1-param) of
 * all the matrix elements is larger than this threshold.  This is accomplished
 * by finding the element in each row that is smaller than (1-param) of the
 * elements in that row, and averaging these elements over all rows.  The
 * threshold is selected on the diagonally scaled matrix.
 *--------------------------------------------------------------------------*/

static double SelectThresh(MPI_Comm comm, Matrix *A, DiagScale *diag_scale,
  double param)
{
    int row, len, *ind, i, npes;
    double *val;
    double localsum = 0.0, sum;
    double temp;

    /* Buffer for storing the values in each row when computing the
       i-th smallest element - buffer will grow if necessary */
    double *buffer;
    int buflen = 10;
    buffer = (double *) malloc(buflen * sizeof(double));

    for (row=0; row<=A->end_row - A->beg_row; row++)
    {
        MatrixGetRow(A, row, &len, &ind, &val);

        if (len > buflen)
        {
            free(buffer);
            buflen = len;
            buffer = (double *) malloc(buflen * sizeof(double));
        }

        /* Copy the scaled absolute values into a work buffer */
        temp = DiagScaleGet(diag_scale, row);
        for (i=0; i<len; i++)
        {
            buffer[i] = temp*ABS(val[i])*DiagScaleGet(diag_scale, ind[i]);
            if (ind[i] == row)
                buffer[i] = 0.0; /* diagonal is not same scale as off-diag */
        }

        /* Compute which element to select */
        i = (int) (len * param) + 1;

        /* Select the i-th smallest element */
        localsum += randomized_select(buffer, 0, len-1, i);
    }

    /* Find the average across all processors */
    MPI_Allreduce(&localsum, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Comm_size(comm, &npes);

    free(buffer);
    return sum / (A->end_rows[npes-1] - A->beg_rows[0] + 1);
}

/*--------------------------------------------------------------------------
 * SelectFilter - Similar to SelectThresh, but on the preconditioner.
 * Assumes matrix is in local indexing.
 *--------------------------------------------------------------------------*/

static double SelectFilter(MPI_Comm comm, Matrix *M, DiagScale *diag_scale,
  double param, int symmetric)
{
    int row, len, *ind, i, npes;
    double *val;
    double localsum = 0.0, sum;
    double temp = 1.0;

    /* Buffer for storing the values in each row when computing the
       i-th smallest element - buffer will grow if necessary */
    double *buffer;
    int buflen = 10;
    buffer = (double *) malloc(buflen * sizeof(double));

    for (row=0; row<=M->end_row - M->beg_row; row++)
    {
        MatrixGetRow(M, row, &len, &ind, &val);

        if (len > buflen)
        {
            free(buffer);
            buflen = len;
            buffer = (double *) malloc(buflen * sizeof(double));
        }

        if (symmetric == 0)
            temp = 1. / DiagScaleGet(diag_scale, row);

        /* Copy the scaled absolute values into a work buffer */
        for (i=0; i<len; i++)
        {
            buffer[i] = temp * ABS(val[i]) / DiagScaleGet(diag_scale, ind[i]);
            if (ind[i] == row)
                buffer[i] = 0.0;
        }

        /* Compute which element to select */
        i = (int) (len * param) + 1;

        /* Select the i-th smallest element */
        localsum += randomized_select(buffer, 0, len-1, i);
    }

    /* Find the average across all processors */
    MPI_Allreduce(&localsum, &sum, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Comm_size(comm, &npes);

    free(buffer);
    return sum / (M->end_rows[npes-1] - M->beg_rows[0] + 1);
}

/*--------------------------------------------------------------------------
 * FilterValues - Filter the values in a preconditioner matrix.
 * M - original matrix, in local ordering
 * F - new matrix, that has been created already
 * Also, return the cost estimate, in case SetupValues is called again
 * with load balancing - the old cost estimate would be incorrect.
 *--------------------------------------------------------------------------*/

static void FilterValues(Matrix *M, Matrix *F, DiagScale *diag_scale,
  double filter, int symmetric, double *newcostp)
{
    int i, j;
    int row, len, *ind;
    double *val, temp = 1.0;
    double cost = 0.0;

    for (row=0; row<=M->end_row - M->beg_row; row++)
    {
        MatrixGetRow(M, row, &len, &ind, &val);

        j = 0;
        for (i=0; i<len; i++)
        {
            if (symmetric == 0)
                temp = 1. / DiagScaleGet(diag_scale, row);

            if (temp * ABS(val[i]) / DiagScaleGet(diag_scale, ind[i]) >= filter
              || row == ind[i])
            {
                val[j] = val[i];
                ind[j] = ind[i];
                j++;
            }
        }

        MatrixSetRow(F, row+F->beg_row, j, ind, val);

        cost += (double) j*j*j;
    }

    *newcostp = cost;
}

/*--------------------------------------------------------------------------
 * Rescale - Rescaling to be used after filtering, in symmetric case.
 *--------------------------------------------------------------------------*/

static void Rescale(Matrix *M, StoredRows *stored_rows, int num_ind)
{
    int len, *ind, len2, *ind2;
    double *val, *val2, *w;
    int row, j, i;
    double accum, prod;

    /* Allocate full-length workspace */
    w = (double *) calloc(num_ind, sizeof(double));

    /* Loop over rows */
    for (row=0; row<=M->end_row - M->beg_row; row++)
    {
        MatrixGetRow(M, row, &len, &ind, &val);

        accum = 0.0;

        /* Loop over nonzeros in row */
        for (j=0; j<len; j++)
        {
            /* Get the row of A corresponding to current nonzero */
            StoredRowsGet(stored_rows, ind[j], &len2, &ind2, &val2);

            /* Scatter nonzeros of A */
            for (i=0; i<len2; i++)
            {
                assert(ind2[i] < num_ind);
                w[ind2[i]] = val2[i];
            }

            /* Form inner product of current row with this row */
            prod = 0.0;
            for (i=0; i<len; i++)
            {
                assert(ind[i] < num_ind);
                prod += val[i] * w[ind[i]];
            }

            accum += val[j] * prod;

            /* Reset workspace */
            for (i=0; i<len2; i++)
                w[ind2[i]] = 0.0;
        }

        /* Scale the row */
        accum = 1./sqrt(accum);
        for (j=0; j<len; j++)
            val[j] *= accum;
    }
}

/******************************************************************************
 *
 * ParaSails public functions
 *
 * After creating a ParaSails object, the preconditioner requires two set up
 * steps:  one for the pattern, and one for the numerical values.  Once the
 * pattern has been set up, the numerical values can be set up for different
 * matrices, i.e., ParaSailsSetupValues can be called again with a different
 * matrix, and used in another iterative solve.
 *
 *****************************************************************************/

/*--------------------------------------------------------------------------
 * ParaSailsCreate - Allocate, initialize, and return a pointer to a
 * ParaSails preconditioner data structure.
 *--------------------------------------------------------------------------*/

ParaSails *ParaSailsCreate(MPI_Comm comm, int beg_row, int end_row, int sym)
{
    ParaSails *ps = (ParaSails *) malloc(sizeof(ParaSails));
    int npes;

    ps->symmetric          = sym;
    ps->thresh             = 0.1;
    ps->num_levels         = 1;
    ps->filter             = 0.0;
    ps->loadbal_beta       = 0.0;
    ps->cost               = 0.0;
    ps->setup_pattern_time = 0.0;
    ps->setup_values_time  = 0.0;
    ps->numb               = NULL;
    ps->M                  = NULL;
    ps->comm               = comm;
    ps->beg_row            = beg_row;
    ps->end_row            = end_row;

    MPI_Comm_size(comm, &npes);

    ps->beg_rows = (int *) malloc(npes * sizeof(int));
    ps->end_rows = (int *) malloc(npes * sizeof(int));

    MPI_Allgather(&beg_row, 1, MPI_INT, ps->beg_rows, 1, MPI_INT, comm);
    MPI_Allgather(&end_row, 1, MPI_INT, ps->end_rows, 1, MPI_INT, comm);

    return ps;
}

/*--------------------------------------------------------------------------
 * ParaSailsDestroy - Deallocate a ParaSails data structure.
 *--------------------------------------------------------------------------*/

void ParaSailsDestroy(ParaSails *ps)
{
    if (ps == NULL)
        return;

    if (ps->numb)
        NumberingDestroy(ps->numb);

    if (ps->M)
        MatrixDestroy(ps->M);

    free(ps->beg_rows);
    free(ps->end_rows);

    free(ps);
}

/*--------------------------------------------------------------------------
 * ParaSailsSetupPattern - Set up a pattern for the ParaSails preconditioner.
 *--------------------------------------------------------------------------*/

void ParaSailsSetupPattern(ParaSails *ps, Matrix *A,
  double thresh, int num_levels)
{
    DiagScale  *diag_scale;
    PrunedRows *pruned_rows;
    double time0, time1;

    time0 = MPI_Wtime();

    ps->thresh     = thresh;
    ps->num_levels = num_levels;

    if (ps->numb) NumberingDestroy(ps->numb);
    ps->numb = NumberingCreateCopy(A->numb);

    if (ps->M) MatrixDestroy(ps->M);
    ps->M = MatrixCreate(ps->comm, ps->beg_row, ps->end_row);

    diag_scale = DiagScaleCreate(A, A->numb);

    if (ps->thresh < 0.0)
        ps->thresh = SelectThresh(ps->comm, A, diag_scale, -ps->thresh);

    pruned_rows = PrunedRowsCreate(A, PARASAILS_NROWS, diag_scale, ps->thresh);

    ExchangePrunedRows(ps->comm, A, ps->numb, pruned_rows, ps->num_levels);

    ConstructPatternForEachRow(ps->symmetric, pruned_rows, ps->num_levels,
        ps->numb, ps->M, &ps->cost);

    DiagScaleDestroy(diag_scale);
    PrunedRowsDestroy(pruned_rows);

    time1 = MPI_Wtime();
    ps->setup_pattern_time = time1 - time0;
}

/*--------------------------------------------------------------------------
 * ParaSailsSetupPatternExt - Set up a pattern for the ParaSails preconditioner.
 * Extended version.
 *--------------------------------------------------------------------------*/

void ParaSailsSetupPatternExt(ParaSails *ps, Matrix *A,
  double thresh_global, double thresh_local, int num_levels)
{
    DiagScale  *diag_scale;
    PrunedRows *pruned_rows_global;
    PrunedRows *pruned_rows_local;
    double time0, time1;

    time0 = MPI_Wtime();

    ps->thresh     = thresh_global*1000000.+thresh_local; /* dummy */
    ps->num_levels = num_levels;

    if (ps->numb) NumberingDestroy(ps->numb);
    ps->numb = NumberingCreateCopy(A->numb);

    if (ps->M) MatrixDestroy(ps->M);
    ps->M = MatrixCreate(ps->comm, ps->beg_row, ps->end_row);

    diag_scale = DiagScaleCreate(A, A->numb);

    if (ps->thresh < 0.0)
        ps->thresh = SelectThresh(ps->comm, A, diag_scale, -ps->thresh);

    pruned_rows_global = PrunedRowsCreate(A, PARASAILS_NROWS, diag_scale, 
         thresh_global);
    pruned_rows_local = PrunedRowsCreate(A, PARASAILS_NROWS, diag_scale, 
         thresh_local);

    ExchangePrunedRowsExt(ps->comm, A, ps->numb, 
        pruned_rows_global, pruned_rows_local, ps->num_levels);

    ExchangePrunedRowsExt2(ps->comm, A, ps->numb, 
        pruned_rows_global, pruned_rows_local, ps->num_levels);

    ConstructPatternForEachRowExt(ps->symmetric, pruned_rows_global,
	pruned_rows_local, ps->num_levels, ps->numb, ps->M, &ps->cost);

    DiagScaleDestroy(diag_scale);
    PrunedRowsDestroy(pruned_rows_global);
    PrunedRowsDestroy(pruned_rows_local);

    time1 = MPI_Wtime();
    ps->setup_pattern_time = time1 - time0;
}

/*--------------------------------------------------------------------------
 * ParaSailsSetupValues - Compute the numerical values of the ParaSails
 * preconditioner, for the pattern set up using ParaSailsSetupPattern.
 * This function may be called repeatedly with different input matrices
 * "A", for which a preconditioner is constructed.
 *--------------------------------------------------------------------------*/

void ParaSailsSetupValues(ParaSails *ps, Matrix *A, double filter)
{
    LoadBal    *load_bal;
    StoredRows *stored_rows;
    int row, len, *ind;
    double *val;
    int i;
    double time0, time1;

    time0 = MPI_Wtime();

    /*
     * If the preconditioner matrix has its own numbering object, then we
     * assume it is in its own local numbering, and we change the numbering
     * in the matrix to the ParaSails numbering.
     */

    if (ps->M->numb != NULL)
    {
        for (row=0; row<=ps->M->end_row - ps->M->beg_row; row++)
        {
           MatrixGetRow(ps->M, row, &len, &ind, &val);
           NumberingLocalToGlobal(ps->M->numb, len, ind, ind);
           NumberingGlobalToLocal(ps->numb,    len, ind, ind);
        }
    }

    load_bal = LoadBalDonate(ps->comm, ps->M, ps->numb, ps->cost,
        ps->loadbal_beta);

    stored_rows = StoredRowsCreate(A, PARASAILS_NROWS);

    ExchangeStoredRows(ps->comm, A, ps->M, ps->numb, stored_rows, load_bal);

    if (ps->symmetric)
    {
        ComputeValuesSym(stored_rows, ps->M, load_bal->beg_row, ps->numb,
            ps->symmetric);

        for (i=0; i<load_bal->num_taken; i++)
        {
            ComputeValuesSym(stored_rows,
                load_bal->recip_data[i].mat,
                load_bal->recip_data[i].mat->beg_row, ps->numb,
                ps->symmetric);
        }
    }
    else
    {
        ComputeValuesNonsym(stored_rows, ps->M, load_bal->beg_row, ps->numb);

        for (i=0; i<load_bal->num_taken; i++)
        {
            ComputeValuesNonsym(stored_rows,
                load_bal->recip_data[i].mat,
                load_bal->recip_data[i].mat->beg_row, ps->numb);
        }
    }

    time1 = MPI_Wtime();
    ps->setup_values_time = time1 - time0;

    LoadBalReturn(load_bal, ps->comm, ps->M);

    /* Filtering */

    ps->filter = filter;

    if (ps->filter != 0.0)
    {
        DiagScale *diag_scale = DiagScaleCreate(A, ps->numb);
        Matrix    *filtered_matrix = MatrixCreate(ps->comm,
                                         ps->beg_row, ps->end_row);

        if (ps->filter < 0.0)
            ps->filter = SelectFilter(ps->comm, ps->M, diag_scale, -ps->filter,
                ps->symmetric);

        FilterValues(ps->M, filtered_matrix, diag_scale, ps->filter,
            ps->symmetric, &ps->cost);

        DiagScaleDestroy(diag_scale);
        MatrixDestroy(ps->M);
        ps->M = filtered_matrix;

        /* Rescale if factored preconditioner */
        if (ps->symmetric != 0)
            Rescale(ps->M, stored_rows, ps->numb->num_ind);
    }

    /*
     * If the preconditioner matrix has its own numbering object, then we
     * change the numbering in the matrix to this numbering.  If not, then
     * we put the preconditioner matrix in global numbering, and call
     * MatrixComplete (to create numbering object, convert the indices,
     * and create the matvec info).
     */

    if (ps->M->numb != NULL)
    {
        /* Convert to own numbering system */
        for (row=0; row<=ps->M->end_row - ps->M->beg_row; row++)
        {
            MatrixGetRow(ps->M, row, &len, &ind, &val);
            NumberingLocalToGlobal(ps->numb,    len, ind, ind);
            NumberingGlobalToLocal(ps->M->numb, len, ind, ind);
        }
    }
    else
    {
        /* Convert to global numbering system and call MatrixComplete */
        for (row=0; row<=ps->M->end_row - ps->M->beg_row; row++)
        {
            MatrixGetRow(ps->M, row, &len, &ind, &val);
            NumberingLocalToGlobal(ps->numb, len, ind, ind);
        }

        MatrixComplete(ps->M);
    }

    StoredRowsDestroy(stored_rows);
}

/*--------------------------------------------------------------------------
 * ParaSailsApply - Apply the ParaSails preconditioner
 *
 * ps - input ParaSails object
 * u  - input array of doubles
 * v  - output array of doubles
 *
 * Although this computation can be done in place, it typically will not
 * be used this way, since the caller usually needs to preserve the input
 * vector.
 *--------------------------------------------------------------------------*/

void ParaSailsApply(ParaSails *ps, double *u, double *v)
{
    if (ps->symmetric)
    {
        MatrixMatvec(ps->M, u, v);      /* need to preserve u */
        MatrixMatvecTrans(ps->M, v, v); /* do the second mult in place */
    }
    else
    {
        MatrixMatvec(ps->M, u, v);
    }
}

/*--------------------------------------------------------------------------
 * ParaSailsStatsPattern - Print some statistics about ParaSailsSetupPattern.
 * Returns a cost, which can be used to preempt ParaSailsSetupValues if the
 * cost is too high.
 *--------------------------------------------------------------------------*/

double ParaSailsStatsPattern(ParaSails *ps, Matrix *A)
{
    int mype, npes;
    int n, nnzm, nnza;
    MPI_Comm comm = ps->comm;
    double max_pattern_time, max_cost, ave_cost;

    MPI_Comm_rank(comm, &mype);
    MPI_Comm_size(comm, &npes);

    nnzm = MatrixNnz(ps->M);
    nnza = MatrixNnz(A);
    if (ps->symmetric)
    {
        n = ps->end_rows[npes-1] - ps->beg_rows[0] + 1;
	nnza = (nnza - n) / 2 + n;
    }

    MPI_Allreduce(&ps->setup_pattern_time, &max_pattern_time, 
	1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&ps->cost, &max_cost, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&ps->cost, &ave_cost, 1, MPI_DOUBLE, MPI_SUM, comm);
    ave_cost = ave_cost / (double) npes;

    if (mype)
	return ave_cost;

    if (ps->symmetric == 0)
        max_cost *= 8.0;  /* nonsymmetric method is harder */

    printf("** ParaSails Setup Pattern Statistics ***********\n");
    printf("symmetric             : %d\n", ps->symmetric);
    printf("thresh                : %f\n", ps->thresh);
    printf("num_levels            : %d\n", ps->num_levels);
    printf("Max cost (average)    : %7.1e (%7.1e)\n", max_cost, ave_cost);
    printf("Nnz (ratio)           : %d (%5.2f)\n", nnzm, nnzm/(double)nnza);
    printf("Max setup pattern time: %8.1f\n", max_pattern_time);
    printf("*************************************************\n");
    fflush(NULL);

    return ave_cost;
}

/*--------------------------------------------------------------------------
 * ParaSailsStatsValues - Print some statistics about ParaSailsSetupValues.
 *--------------------------------------------------------------------------*/

void ParaSailsStatsValues(ParaSails *ps, Matrix *A)
{
    int mype, npes;
    int n, nnzm, nnza;
    MPI_Comm comm = ps->comm;
    double max_values_time;
    double temp, *setup_times = NULL;
    int i;

    MPI_Comm_rank(comm, &mype);
    MPI_Comm_size(comm, &npes);

    nnzm = MatrixNnz(ps->M);
    nnza = MatrixNnz(A);
    if (ps->symmetric)
    {
        n = ps->end_rows[npes-1] - ps->beg_rows[0] + 1;
        nnza = (nnza - n) / 2 + n;
    }

    MPI_Allreduce(&ps->setup_values_time, &max_values_time, 
	1, MPI_DOUBLE, MPI_MAX, comm);

    if (!mype)
        setup_times = (double *) malloc(npes * sizeof(double));

    temp = ps->setup_pattern_time + ps->setup_values_time;
    MPI_Gather(&temp, 1, MPI_DOUBLE, setup_times, 1, MPI_DOUBLE, 0, comm);

    if (mype)
        return;

    printf("** ParaSails Setup Values Statistics ************\n");
    printf("filter                : %f\n", ps->filter);
    printf("loadbal               : %f\n", ps->loadbal_beta);
    printf("Final Nnz (ratio)     : %d (%5.2f)\n", nnzm, nnzm/(double)nnza);
    printf("Max setup values time : %8.1f\n", max_values_time);
    printf("*************************************************\n");
    printf("Setup (pattern and values) times:\n");

    temp = 0.0;
    for (i=0; i<npes; i++)
    {
        printf("%3d: %8.1f\n", i, setup_times[i]);
        temp += setup_times[i];
    }
    printf("ave: %8.1f\n", temp / (double) npes);
    printf("*************************************************\n");

    free(setup_times);

    fflush(NULL);
}
