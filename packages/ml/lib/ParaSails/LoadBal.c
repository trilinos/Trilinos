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
 * LoadBal - Load balancing module for ParaSails.
 *
 *****************************************************************************/

#include <assert.h>
#include <stdlib.h>
#include "Common.h"
#include "Matrix.h"
#include "Numbering.h"
#include "LoadBal.h"

/*--------------------------------------------------------------------------
 * LoadBalInit - determine the amount of work to be donated and received by 
 * each processor, given the amount of work that each processor has 
 * ("local_cost").  The number of processors that this processor will donate
 * to is "num_given" and the number of processors from which this processor
 * will receive is "num_taken".  Additional donor information is stored in
 * "donor_data_pe" and "donor_data_cost".
 *
 * local_cost - amount of work that this processor has
 * beta - target load balance factor
 *--------------------------------------------------------------------------*/

void LoadBalInit(MPI_Comm comm, double local_cost, double beta, 
  int *num_given, int *donor_data_pe, double *donor_data_cost,
  int *num_taken)
{
    int mype, npes;
    double *cost, average, upper, move, accept;
    int i, jj, j;

    *num_given = 0;
    *num_taken = 0;

    if (beta == 0.0)
	return;

    MPI_Comm_rank(comm, &mype);
    MPI_Comm_size(comm, &npes);

    cost = (double *) malloc(npes * sizeof(double));

    MPI_Allgather(&local_cost, 1, MPI_DOUBLE, cost, 1, MPI_DOUBLE, comm);

    /* Compute the average cost */
    average = 0.0;
    for (i=0; i<npes; i++)
        average += cost[i];
    average = average / npes;

    /* Maximum cost allowed by load balancer */
    upper = average / beta;

    for (i=0; i<npes; i++)
    {
        if (cost[i] > upper)
        {
            move = cost[i] - upper;

            /* for j=[i+1:n 1:i-1] */
            for (jj=i+1; jj<=i+npes; jj++)
            {
		j = jj % npes;
		if (j == i)
		    continue;

                if (cost[j] < average)
                {
                    accept = upper - cost[j];

                    /* If we are sender, record it */
                    if (mype == i)
                    {
                        donor_data_pe[*num_given] = j;
                        donor_data_cost[*num_given] = MIN(move, accept);
                        (*num_given)++;
                    }

                    /* If we are receiver, record it */
                    if (mype == j)
                    {
                        (*num_taken)++;
                    }

                    if (move <= accept)
                    {
                        cost[i] = cost[i] - move;
                        cost[j] = cost[j] + move;
#ifdef PARASAILS_DEBUG
			if (mype == 0)
                            printf("moved from %d to %d (%7.1e)\n", i,j,move);
#endif
                        /*nummoves = nummoves + 1;*/
                        break;
                    }
                    else
                    {
                        cost[i] = cost[i] - accept;
                        cost[j] = cost[j] + accept;
#ifdef PARASAILS_DEBUG
			if (mype == 0)
                            printf("moved from %d to %d (%7.1e)\n", i,j,accept);
#endif
                        /*nummoves = nummoves + 1;*/
                        move = cost[i] - upper;
                    }
                }
            }
        }
    }

    free(cost);
}

/*--------------------------------------------------------------------------
 * LoadBalDonorSend - send the indices of the donated rows.
 * The message structure is: beg_row, end_row, len1, indices1, len2, ....
 * Caller must free the allocated buffers.
 *--------------------------------------------------------------------------*/

void LoadBalDonorSend(MPI_Comm comm, Matrix *mat, Numbering *numb,
  int num_given, const int *donor_data_pe, const double *donor_data_cost, 
  DonorData *donor_data, int *local_beg_row, MPI_Request *request)
{
    int send_beg_row, send_end_row;
    int i, row;
    double accum;
    int buflen;
    int *bufferp;
    int len, *ind;
    double *val;

    send_end_row = mat->beg_row - 1; /* imaginary end of previous block */

    for (i=0; i<num_given; i++)
    {
	send_beg_row = send_end_row + 1;
        send_end_row = send_beg_row - 1;

        /* Portion out rows that add up to the workload to be sent out */
	/* and determine the size of the buffer needed */

        accum = 0.0; /* amount of work portioned out so far */
        buflen = 2;  /* front of buffer will contain beg_row, end_row */

        do
        {
            send_end_row++;
            assert(send_end_row <= mat->end_row);
            MatrixGetRow(mat, send_end_row - mat->beg_row, &len, &ind, &val);
            accum += (double) len*len*len;
            buflen += (len+1); /* additional one for row length */
        }
        while (accum < donor_data_cost[i]);

        /* Create entry in donor_data structure */

        donor_data[i].pe      = donor_data_pe[i];
        donor_data[i].beg_row = send_beg_row;
        donor_data[i].end_row = send_end_row;
        donor_data[i].buffer  = (int *) malloc((buflen) * sizeof(int));

	/* Construct send buffer */

         bufferp   = donor_data[i].buffer;
        *bufferp++ = send_beg_row;
        *bufferp++ = send_end_row;

        for (row=send_beg_row; row<=send_end_row; row++)
        {
            MatrixGetRow(mat, row - mat->beg_row, &len, &ind, &val);
            *bufferp++ = len;
            /* memcpy(bufferp, ind, len*sizeof(int)); */ /* copy into buffer */
	    NumberingLocalToGlobal(numb, len, ind, bufferp);
            bufferp += len;
        }

        MPI_Isend(donor_data[i].buffer, buflen, MPI_INT, donor_data[i].pe,
            LOADBAL_REQ_TAG, comm, &request[i]);
    }

    *local_beg_row = send_end_row + 1;
}

/*--------------------------------------------------------------------------
 * LoadBalRecipRecv - receive the indices of the donated rows.
 * The message structure is: beg_row, end_row, len1, indices1, len2, ....
 *--------------------------------------------------------------------------*/

void LoadBalRecipRecv(MPI_Comm comm, Numbering *numb,
  int num_taken, RecipData *recip_data)
{
    int i, row;
    int count;
    MPI_Status status;
    int *buffer, *bufferp;
    int beg_row, end_row;
    int len;

    for (i=0; i<num_taken; i++)
    {
        MPI_Probe(MPI_ANY_SOURCE, LOADBAL_REQ_TAG, comm, &status);
        recip_data[i].pe = status.MPI_SOURCE;
        MPI_Get_count(&status, MPI_INT, &count);

        buffer = (int *) malloc(count * sizeof(int));
        MPI_Recv(buffer, count, MPI_INT, recip_data[i].pe, LOADBAL_REQ_TAG, 
           comm, &status);

	bufferp =  buffer;
        beg_row = *bufferp++;
        end_row = *bufferp++;

        recip_data[i].mat = MatrixCreateLocal(beg_row, end_row);

	/* Set the indices of the local matrix containing donated rows */

        for (row=beg_row; row<=end_row; row++)
        {
            len = *bufferp++;
	    NumberingGlobalToLocal(numb, len, bufferp, bufferp);
            MatrixSetRow(recip_data[i].mat, row, len, bufferp, NULL);
            bufferp += len;
        }

	free(buffer);
    }
}

/*--------------------------------------------------------------------------
 * LoadBalRecipSend - send back the computed values of the donated rows.
 * Traverse all the donated local matrices.
 * Assume indices are in the same order.
 * Caller must free the allocated buffers.
 *--------------------------------------------------------------------------*/

void LoadBalRecipSend(MPI_Comm comm, int num_taken, 
  RecipData *recip_data, MPI_Request *request)
{
    int i, row, buflen;
    double *bufferp;
    Matrix *mat;
    int len, *ind;
    double *val;

    for (i=0; i<num_taken; i++)
    {
        mat = recip_data[i].mat;

        /* Find size of output buffer */
	buflen = 0;
        for (row=0; row<=mat->end_row - mat->beg_row; row++)
        {
            MatrixGetRow(mat, row, &len, &ind, &val);
	    buflen += len;
	}

	recip_data[i].buffer = (double *) malloc(buflen * sizeof(double));

	/* Construct send buffer */

	bufferp = recip_data[i].buffer;
        for (row=0; row<=mat->end_row - mat->beg_row; row++)
        {
            MatrixGetRow(mat, row, &len, &ind, &val);
            memcpy(bufferp, val, len*sizeof(double)); /* copy into buffer */
            bufferp += len;
        }

        MPI_Isend(recip_data[i].buffer, buflen, MPI_DOUBLE, recip_data[i].pe,
            LOADBAL_REP_TAG, comm, &request[i]);

        MatrixDestroy(mat);
    }
}

/*--------------------------------------------------------------------------
 * LoadBalDonorRecv - receive the computed values of the donated rows.
 * Traverse all the donated local matrices.
 * Assume indices are in the same order.
 *--------------------------------------------------------------------------*/

void LoadBalDonorRecv(MPI_Comm comm, Matrix *mat, 
  int num_given, DonorData *donor_data)
{
    int i, j, row;
    int source, count;
    MPI_Status status;
    double *buffer, *bufferp;
    int len, *ind;
    double *val;

    for (i=0; i<num_given; i++)
    {
        MPI_Probe(MPI_ANY_SOURCE, LOADBAL_REP_TAG, comm, &status);
        source = status.MPI_SOURCE;
        MPI_Get_count(&status, MPI_DOUBLE, &count);

        buffer = (double *) malloc(count * sizeof(double));
        MPI_Recv(buffer, count, MPI_DOUBLE, source, LOADBAL_REP_TAG, 
           comm, &status);

	/* search for which entry in donor_data this message corresponds to */
	for (j=0; j<num_given; j++)
	{
	    if (donor_data[j].pe == source)
		break;
	}
	assert(j < num_given);

        /* Parse the message and put row values into local matrix */
	bufferp = buffer;
        for (row=donor_data[j].beg_row; row<=donor_data[j].end_row; row++)
        {
            MatrixGetRow(mat, row - mat->beg_row, &len, &ind, &val);
	    memcpy(val, bufferp, len*sizeof(double)); /* copy into matrix */
            bufferp += len;
        }

	free(buffer);
    }
}

/*--------------------------------------------------------------------------
 * LoadBalDonate
 *--------------------------------------------------------------------------*/

LoadBal *LoadBalDonate(MPI_Comm comm, Matrix *mat, Numbering *numb,
  double local_cost, double beta)
{
    LoadBal *p;
    int i, npes;
    int    *donor_data_pe;
    double *donor_data_cost;
    MPI_Request *requests;
    MPI_Status  *statuses;

    p = (LoadBal *) malloc(sizeof(LoadBal));

    MPI_Comm_size(comm, &npes);

    donor_data_pe   = (int *)    malloc(npes * sizeof(int));
    donor_data_cost = (double *) malloc(npes * sizeof(double));

    LoadBalInit(comm, local_cost, beta, &p->num_given, 
        donor_data_pe, donor_data_cost, &p->num_taken);

    p->donor_data = (DonorData *) malloc(p->num_given * sizeof(DonorData));
    p->recip_data = (RecipData *) malloc(p->num_taken * sizeof(RecipData));

    requests = (MPI_Request *) malloc(p->num_given * sizeof(MPI_Request));
    statuses = (MPI_Status  *) malloc(p->num_given * sizeof(MPI_Status));

    LoadBalDonorSend(comm, mat, numb, p->num_given,
        donor_data_pe, donor_data_cost, p->donor_data, &p->beg_row, requests);

    free(donor_data_pe);
    free(donor_data_cost);

    LoadBalRecipRecv(comm, numb, p->num_taken, p->recip_data);

    MPI_Waitall(p->num_given, requests, statuses);

    free(requests);
    free(statuses);

    /* Free the send buffers which were allocated by LoadBalDonorSend */
    for (i=0; i<p->num_given; i++)
	free(p->donor_data[i].buffer);

    return p;
}

/*--------------------------------------------------------------------------
 * LoadBalReturn
 *--------------------------------------------------------------------------*/

void LoadBalReturn(LoadBal *p, MPI_Comm comm, Matrix *mat)
{
    int i;

    MPI_Request *requests;
    MPI_Status  *statuses;

    requests = (MPI_Request *) malloc(p->num_taken * sizeof(MPI_Request));
    statuses = (MPI_Status  *) malloc(p->num_taken * sizeof(MPI_Status));

    LoadBalRecipSend(comm, p->num_taken, p->recip_data, requests);

    LoadBalDonorRecv(comm, mat, p->num_given, p->donor_data);

    MPI_Waitall(p->num_taken, requests, statuses);

    free(requests);
    free(statuses);

    /* Free the send buffers which were allocated by LoadBalRecipSend */
    for (i=0; i<p->num_taken; i++)
	free(p->recip_data[i].buffer);

    free(p->donor_data);
    free(p->recip_data);

    free(p);
}

