/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


#include <stdio.h>
#include <memory.h>
#include <mpi.h>
#include "comm.h"
#include "zoltan_mem.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#if 0
static void show_int_buffers(int me, int procs, char *buf, int *bcounts, int *boffsets)
{
  int i,p;
  int *ibuf = (int *)buf;

  printf("%d  count offset (values)\n",me);

  for (p=0; p < procs; p++){
    printf("%d %d (",bcounts[p]/4, boffsets[p]/4);
    for (i=0; i < bcounts[p]/4; i++){
      printf("%d ",*ibuf);
      ibuf++;
    }
    printf(")\n\n");
    fflush(stdout);
  }
}
#endif

/* Given the communication object, perform a communication operation as
   efficiently as possible. */

/* Several cases are supported:

   I. All objects are of same size (sizes == NULL)
     A. Data is already packed by destination proc (indices_to == NULL)
        Can just send straight out of user space.
     B. Data not packed by destination proc (indices_to != NULL)
	Need to allocate send buff and copy it in before sending

   II. Objects of variable sizes (sizes != NULL)
       Same two cases above, but offsets are more complicated.
       Need stars_to_ptr/_from_ptr to address memory for (A), and
       indices_to_ptr/_from_ptr to address memory for (B).

   III. On receive, we need to do opposite of the send.
	Can receive directly into user space if  indices_from == NULL
	Otherwise need to receive in buffer and copy.
*/

/*****************************************************************************/

int       Zoltan_Comm_Do(
ZOLTAN_COMM_OBJ * plan,		/* communication data structure */
int tag,			/* message tag for communicating */
char *send_data,		/* array of data I currently own */
int nbytes,			/* multiplier for sizes */
char *recv_data)		/* array of data I'll own after comm */
{
  int status = ZOLTAN_OK;

  if (!plan->maxed_recvs){
    status = Zoltan_Comm_Do_Post (plan, tag, send_data, nbytes, recv_data);
    if (status == ZOLTAN_OK)
       status = Zoltan_Comm_Do_Wait (plan, tag, send_data, nbytes, recv_data);
  }
  else{
    status = Zoltan_Comm_Do_AlltoAll(plan, send_data, nbytes, recv_data);
  }

  return status;
}

/*****************************************************************************/

int       Zoltan_Comm_Do_Post(
ZOLTAN_COMM_OBJ * plan,		/* communication data structure */
int tag,			/* message tag for communicating */
char *send_data,		/* array of data I currently own */
int nbytes,			/* multiplier for sizes */
char *recv_data)		/* array of data I'll own after comm */
{
    char     *send_buff;	/* space to buffer outgoing data */
    int       my_proc;		/* processor ID */
    size_t    self_recv_address = 0;/* where in recv_data self info starts */
    int       self_num=0;       /* where in send list my_proc appears */
    size_t    offset;		/* offset into array I'm copying into */
    int       self_index = 0;	/* send offset for data I'm keeping */
    int       out_of_mem;	/* am I out of memory? */
    int       nblocks;		/* number of procs who need my data */
    int       proc_index;	/* loop counter over procs to send to */
    int       i, j, k, jj;	/* loop counters */

    static char *yo = "Zoltan_Comm_Do_Post";


    /* Check input parameters */
    if (!plan) {
        MPI_Comm_rank(MPI_COMM_WORLD, &my_proc);
	ZOLTAN_COMM_ERROR("Communication plan = NULL", yo, my_proc);
	return ZOLTAN_FATAL;
    }

    /* If not point to point, currently we do synchroneous communications */
    if (plan->maxed_recvs){
      int status;
      status = Zoltan_Comm_Do_AlltoAll(plan, send_data, nbytes, recv_data);
      return (status);
    }

    MPI_Comm_rank(plan->comm, &my_proc);

    if ((plan->nsends + plan->self_msg) && !send_data) {
        size_t sum = 0;
        if (plan->sizes_to)   /* Not an error if all sizes_to == 0 */
            for (i = 0; i < (plan->nsends + plan->self_msg); i++) 
                sum += plan->sizes_to[i];
        if (!plan->sizes_to || (plan->sizes_to && sum)) {
            ZOLTAN_COMM_ERROR("nsends not zero, but send_data = NULL", 
                              yo, my_proc);
            return ZOLTAN_FATAL;
        }
    }
    if ((plan->nrecvs + plan->self_msg) && !recv_data) {
        size_t sum = 0;
        if (plan->sizes_from)   /* Not an error if all sizes_from == 0 */
            for (i = 0; i < (plan->nrecvs + plan->self_msg); i++) 
                sum += plan->sizes_from[i];
        if (!plan->sizes_from || (plan->sizes_from && sum)) {
            ZOLTAN_COMM_ERROR("nrecvs not zero, but recv_data = NULL", 
                              yo, my_proc);
            return ZOLTAN_FATAL;
        }
    }
    if (nbytes < 0) {
	ZOLTAN_COMM_ERROR("Scale factor nbytes is negative", yo, my_proc);
	return ZOLTAN_FATAL;
    }


    /* Post irecvs */

    out_of_mem = 0;

    if (plan->indices_from == NULL) {
	/* Data can go directly into user space. */
	plan->recv_buff = recv_data;
    }
    else {			/* Need to buffer receive to reorder */
        size_t rsize = (size_t) (plan->total_recv_size) * (size_t) nbytes;
	plan->recv_buff = (char *) ZOLTAN_MALLOC(rsize);
	if (plan->recv_buff == NULL && rsize != 0)
	    out_of_mem = 1;
    }

    if (!out_of_mem) {
	if (plan->sizes == NULL) {	/* All data the same size */
	    k = 0;
	    for (i = 0; i < plan->nrecvs + plan->self_msg; i++) {
		if (plan->procs_from[i] != my_proc) {
		    MPI_Irecv((void *)
                              &plan->recv_buff[(size_t)(plan->starts_from[i]) * (size_t)nbytes],
			      plan->lengths_from[i] * nbytes,
			      (MPI_Datatype) MPI_BYTE, plan->procs_from[i], tag,
			      plan->comm, &plan->request[k]);
		    k++;
		}
		else {
		    self_recv_address = (size_t)(plan->starts_from[i]) * (size_t)nbytes;
		}
	    }
	}

	else {			/* Data of varying sizes */
	    k = 0;
	    for (i = 0; i < plan->nrecvs + plan->self_msg; i++) {
		if (plan->procs_from[i] != my_proc) {
                    if (plan->sizes_from[i]) {
		        MPI_Irecv((void *)
                            &plan->recv_buff[(size_t)(plan->starts_from_ptr[i]) 
                                              * (size_t)nbytes],
		                  plan->sizes_from[i] * nbytes,
			          (MPI_Datatype) MPI_BYTE, plan->procs_from[i], 
			          tag, plan->comm, &plan->request[k]);
                    }
                    else
                        plan->request[k] = MPI_REQUEST_NULL;
	            k++;
		}
		else {
		    self_recv_address = (size_t)(plan->starts_from_ptr[i]) * (size_t)nbytes;
		}
	    }
	}
    }


    /* Do remaining allocation to check for any mem problems. */
    if (plan->indices_to != NULL) {	/* can't sent straight from input */
        size_t ssize = (size_t)(plan->max_send_size) * (size_t)nbytes;
	send_buff = (char *) ZOLTAN_MALLOC(ssize);
	if (send_buff == 0 && ssize != 0)
	    out_of_mem = 1;
    }
    else
	send_buff = NULL;

    /* Barrier to ensure irecvs are posted before doing any sends. */
    /* Simultaneously see if anyone out of memory */

    MPI_Allreduce(&out_of_mem, &j, 1, MPI_INT, MPI_SUM, plan->comm);

    if (j > 0) {		/* Some proc is out of memory -> Punt */
	ZOLTAN_FREE(&send_buff);
	if (plan->indices_from != NULL)
	    ZOLTAN_FREE(&plan->recv_buff);
	return (ZOLTAN_MEMERR);
    }

    /* Send out data */

    /* Scan through procs_to list to start w/ higher numbered procs */
    /* This should balance message traffic. */

    nblocks = plan->nsends + plan->self_msg;
    proc_index = 0;
    while (proc_index < nblocks && plan->procs_to[proc_index] < my_proc)
	proc_index++;
    if (proc_index == nblocks)
	proc_index = 0;

    if (plan->sizes == NULL) {	/* Data all of same size */
	if (plan->indices_to == NULL) {	/* data already blocked by processor. */
	    for (i = proc_index, j = 0; j < nblocks; j++) {
		if (plan->procs_to[i] != my_proc) {
		    MPI_Rsend((void *) &send_data[(size_t)(plan->starts_to[i]) * (size_t)nbytes],
			      plan->lengths_to[i] * nbytes,
			      (MPI_Datatype) MPI_BYTE, plan->procs_to[i], tag,
			      plan->comm);
		}
		else
		    self_num = i;
		if (++i == nblocks)
		    i = 0;
	    }

	    if (plan->self_msg) {	/* Copy data to self. */
		/* I use array+offset instead of &(array[offset]) because of
		   a bug with PGI v9 */
		/* I use memmove because I'm not sure that the pointer are not
		   overlapped. */
		memmove(
                  plan->recv_buff+self_recv_address,
                  send_data+(size_t)(plan->starts_to[self_num])*(size_t)nbytes,
                  (size_t) (plan->lengths_to[self_num]) * (size_t) nbytes);
	    }
	}

	else {			/* Not blocked by processor.  Need to buffer. */
	    for (i = proc_index, jj = 0; jj < nblocks; jj++) {
		if (plan->procs_to[i] != my_proc) {
		    /* Need to pack message first. */
		    offset = 0;
		    j = plan->starts_to[i];
		    for (k = 0; k < plan->lengths_to[i]; k++) {
			memcpy(&send_buff[offset],
			       &send_data[(size_t)(plan->indices_to[j++]) * (size_t)nbytes], nbytes);
			offset += nbytes;
		    }
		    MPI_Rsend((void *) send_buff, plan->lengths_to[i] * nbytes,
		      (MPI_Datatype) MPI_BYTE, plan->procs_to[i], tag, plan->comm);
		}
		else {
		    self_num = i;
		    self_index = plan->starts_to[i];
		}
		if (++i == nblocks)
		    i = 0;
	    }
	    if (plan->self_msg) {	/* Copy data to self. */
		for (k = 0; k < plan->lengths_to[self_num]; k++) {
		    memcpy(&plan->recv_buff[self_recv_address],
		      &send_data[(size_t)(plan->indices_to[self_index++]) * (size_t)nbytes], nbytes);
		    self_recv_address += nbytes;
		}
	    }

	    ZOLTAN_FREE(&send_buff);
	}
    }
    else {			/* Data of differing sizes */
	if (plan->indices_to == NULL) {	/* data already blocked by processor. */
	    for (i = proc_index, j = 0; j < nblocks; j++) {

		if (plan->procs_to[i] != my_proc) {
                    if (plan->sizes_to[i]) {
		        MPI_Rsend((void *)
                                  &send_data[(size_t)(plan->starts_to_ptr[i]) * (size_t)nbytes],
			          plan->sizes_to[i] * nbytes,
			          (MPI_Datatype) MPI_BYTE, plan->procs_to[i],
			          tag, plan->comm);
                    }
		}
		else
		    self_num = i;
		if (++i == nblocks)
		    i = 0;
	    }

	    if (plan->self_msg) {	/* Copy data to self. */
                if (plan->sizes_to[self_num]) {
                    char* lrecv = &plan->recv_buff[self_recv_address];
                    char* lsend = &send_data[(size_t)(plan->starts_to_ptr[self_num]) * (size_t)nbytes];
                    int sindex = plan->sizes_to[self_num], idx;
                    for (idx=0; idx<nbytes; idx++) {
                        memcpy(lrecv, lsend, sindex);
                        lrecv += sindex;
                        lsend += sindex;
                    }
                }
	    }
	}

	else {			/* Not blocked by processor.  Need to buffer. */
	    for (i = proc_index, jj = 0; jj < nblocks; jj++) {
		if (plan->procs_to[i] != my_proc) {
		    /* Need to pack message first. */
		    offset = 0;
		    j = plan->starts_to[i];
		    for (k = 0; k < plan->lengths_to[i]; k++) {
                        if (plan->sizes[plan->indices_to[j]]) {
			    memcpy(&send_buff[offset],
			       &send_data[(size_t)(plan->indices_to_ptr[j]) * (size_t)nbytes],
			       (size_t)(plan->sizes[plan->indices_to[j]]) * (size_t)nbytes);
			    offset += (size_t)(plan->sizes[plan->indices_to[j]]) * (size_t)nbytes;
                        }
			j++;
		    }
                    if (plan->sizes_to[i]) {
		        MPI_Rsend((void *) send_buff, 
                                  plan->sizes_to[i] * nbytes,
		                  (MPI_Datatype) MPI_BYTE, plan->procs_to[i],
                                  tag, plan->comm);
                    }
		}
		else
		    self_num = i;
		if (++i == nblocks)
		    i = 0;
	    }
	    if (plan->self_msg) {	/* Copy data to self. */
                if (plan->sizes_to[self_num]) {
		    j = plan->starts_to[self_num];
		    for (k = 0; k < plan->lengths_to[self_num]; k++) {
		        int kk = plan->indices_to_ptr[j];
                        char* lrecv = &plan->recv_buff[self_recv_address];
                        size_t send_idx = (size_t)kk * (size_t)nbytes;
                        char* lsend = &send_data[send_idx];
                        int sindex = plan->sizes[plan->indices_to[j]], idx;
                        for (idx=0; idx<nbytes; idx++) {
                            memcpy(lrecv, lsend, sindex);
                            lrecv += sindex;
                            lsend += sindex;
                        }
		        self_recv_address += (size_t)(plan->sizes[plan->indices_to[j]])
                                           * (size_t) nbytes;
		        j++;
		    }
		}
	    }

	    ZOLTAN_FREE(&send_buff);
	}
    }
    return (ZOLTAN_OK);
}


/*****************************************************************************/

int       Zoltan_Comm_Do_Wait(
ZOLTAN_COMM_OBJ * plan,		/* communication data structure */
int tag,			/* message tag for communicating */
char *send_data,		/* array of data I currently own */
int nbytes,			/* multiplier for sizes */
char *recv_data)		/* array of data I'll own after comm */
{
    MPI_Status status;		/* return from Waitany */
    int       my_proc;		/* processor ID */
    int       self_num;		/* where in send list my_proc appears */
    int       i, j, k, jj;	/* loop counters */

    /* If not point to point, currently we do synchroneous communications */
    if (plan->maxed_recvs){
      /* Do nothing */
      return (ZOLTAN_OK);
    }

    MPI_Comm_rank(plan->comm, &my_proc);

    /* Wait for messages to arrive & unpack them if necessary. */
    /* Note: since request is in plan, could wait in later routine. */

    if (plan->indices_from == NULL) {	/* No copying required */
        if (plan->nrecvs > 0) {
	    MPI_Waitall(plan->nrecvs, plan->request, plan->status);
	}
    }

    else {			 	/* Need to copy into recv_data. */
	if (plan->self_msg) {		/* Unpack own data before waiting */
	    for (self_num = 0; self_num < plan->nrecvs + plan->self_msg; self_num++) 
		if (plan->procs_from[self_num] == my_proc) break;
	    k = plan->starts_from[self_num];
            if (!plan->sizes_from || plan->sizes_from[self_num]) {
	        for (j = plan->lengths_from[self_num]; j; j--) {
		    memcpy(&recv_data[(size_t)(plan->indices_from[k]) * (size_t)nbytes],
		        &plan->recv_buff[(size_t)k * (size_t)nbytes], nbytes);
		    k++;
	        }
	    }
	}
	else
	    self_num = plan->nrecvs;

	for (jj = 0; jj < plan->nrecvs; jj++) {

	    MPI_Waitany(plan->nrecvs, plan->request, &i, &status);

            if (i == MPI_UNDEFINED) break;  /* No more receives */

	    if (i >= self_num) i++;

	    k = plan->starts_from[i];
	    for (j = plan->lengths_from[i]; j; j--) {
		memcpy(&recv_data[(size_t)(plan->indices_from[k]) * (size_t)nbytes],
		    &plan->recv_buff[(size_t)k * (size_t)nbytes], nbytes);
		k++;
	    }
	}

	ZOLTAN_FREE(&plan->recv_buff);
    }

    return (ZOLTAN_OK);
}

/*****************************************************************************/

/* Do_Post would require posting more receives than allowed on this platform.
*  We use MPI_AlltoAllv instead, which is probably implemented such that each
*  process does one receive at a time.
*/

int       Zoltan_Comm_Do_AlltoAll(
ZOLTAN_COMM_OBJ * plan,		/* communication data structure */
char *send_data,		/* array of data I currently own */
int nbytes,			/* multiplier for sizes */
char *recv_data)		/* array of data I'll own after comm */
{
  static char *yo = "Zoltan_Comm_Do_AlltoAll";
  char *outbuf=NULL, *inbuf=NULL, *buf=NULL;
  int *outbufCounts=NULL, *outbufOffsets=NULL; 
  int *inbufCounts=NULL, *inbufOffsets=NULL;
  int nprocs, me, rc, i, j, k, p, sorted;
  int nSendMsgs, nSendItems, nRecvMsgs, nRecvItems;
  int length, offset, itemSize, outbufLen;

  int sm = (plan->self_msg > 0) ? 1 : 0;

  nSendMsgs = plan->nsends + sm;
  nRecvMsgs = plan->nrecvs + sm;

  for (i=0, nSendItems=0; i <nSendMsgs; i++){
    nSendItems += plan->lengths_to[i];
  }
  for (i=0, nRecvItems=0; i <nRecvMsgs; i++){
    nRecvItems += plan->lengths_from[i];
  }

  MPI_Comm_size(plan->comm, &nprocs);
  MPI_Comm_rank(plan->comm, &me);

  outbufCounts = (int *) ZOLTAN_CALLOC(nprocs , sizeof(int));
  outbufOffsets = (int *) ZOLTAN_CALLOC(nprocs , sizeof(int));
  inbufCounts = (int *) ZOLTAN_CALLOC(nprocs , sizeof(int));
  inbufOffsets = (int *) ZOLTAN_CALLOC(nprocs , sizeof(int));

  if (!outbufCounts || !outbufOffsets || !inbufCounts || !inbufOffsets){
    ZOLTAN_COMM_ERROR("memory error", yo, me);
  }

  /* The *_to fields of the plan refer to the items in the send_data buffer,
   * and how to pull out the correct items for each receiver.  The
   * *_from fields of the plan refer to the recv_data buffer.  Items 
   * arrive in process rank order, and these fields tell us where to
   * put them in the recv_data buffer.
   */

  /* CREATE SEND BUFFER */

  sorted = 0;
  if (plan->indices_to == NULL){
    sorted = 1;
    for (i=1; i< nSendMsgs; i++){
      if (plan->starts_to[i] < plan->starts_to[i-1]){
        sorted = 0;
        break;
      }
    }
  }

  if (plan->sizes_to){
    /*
     * Each message contains items for a process, and each item may be
     * a different size.
     */

    for (i=0, outbufLen=0; i < nSendMsgs; i++){
      outbufLen += plan->sizes_to[i];
    }

    if (plan->indices_to){
      /*
       * items are not grouped by message
       */

      buf = outbuf = (char *)ZOLTAN_MALLOC((size_t)outbufLen * (size_t)nbytes);
      if (outbufLen && nbytes && !outbuf){
        ZOLTAN_COMM_ERROR("memory error", yo, me);
      }

      for (p=0, i=0, k=0; p < nprocs; p++){

        length = 0;

        if (i < nSendMsgs){
          if (plan->procs_to[i] == p){   /* procs_to is sorted */
  
            for (j=0; j < plan->lengths_to[i]; j++,k++){
              itemSize = plan->sizes[plan->indices_to[k]] * nbytes;
              offset = plan->indices_to_ptr[k] * nbytes;
  
              memcpy(buf, send_data + offset, itemSize);
  
              buf += itemSize;
              length += itemSize;
            }
            i++;
          }
        }
  
        outbufCounts[p] = length;
        if (p){
          outbufOffsets[p] = outbufOffsets[p-1] + outbufCounts[p-1];
        }
      }
    }
    else{
      /*
       * items are stored contiguously for each message
       */

      if (!sorted || (plan->nvals > nSendItems) ){

        buf = outbuf = (char *)ZOLTAN_MALLOC((size_t)outbufLen * (size_t)nbytes);
        if (outbufLen && nbytes && !outbuf){
          ZOLTAN_COMM_ERROR("memory error", yo, me);
        }
      }
      else{
        /* All items in send_data are being sent, and they are sorted
         * in process rank order.
         */
        outbuf = send_data;
      }

      for (p=0, i=0; p < nprocs; p++){

        length = 0;

        if (i < nSendMsgs){
          if (plan->procs_to[i] == p){   /* procs_to is sorted */
            length = plan->sizes_to[i] * nbytes;
            offset = plan->starts_to_ptr[i] * nbytes;
  
            if ((!sorted || (plan->nvals > nSendItems)) && length){
              memcpy(buf, send_data + offset, length);
              buf += length;
            }
            i++;
          }
        }
  
        outbufCounts[p] = length;
        if (p){
          outbufOffsets[p] = outbufOffsets[p-1] + outbufCounts[p-1];
        }
      }
    }
  }
  else if (plan->indices_to){
    /*
     * item sizes are constant, however the items belonging in a given
     * message may not be contiguous in send_data
     */

    buf = outbuf = (char *)ZOLTAN_MALLOC((size_t)nSendItems * (size_t)nbytes);
    if (nSendMsgs && nbytes && !outbuf){
      ZOLTAN_COMM_ERROR("memory error", yo, me);
    }

    for (p=0, i=0, k=0; p < nprocs; p++){

      length = 0;
     
      if (i < nSendMsgs){
        if (plan->procs_to[i] == p){   /* procs_to is sorted */
          for (j=0; j < plan->lengths_to[i]; j++,k++){
            offset = plan->indices_to[k] * nbytes;
            memcpy(buf, send_data + offset, nbytes);
            buf += nbytes;
          }
          length = plan->lengths_to[i] * nbytes;
          i++;
        }
      }

      outbufCounts[p] = length;
      if (p){
        outbufOffsets[p] = outbufOffsets[p-1] + outbufCounts[p-1];
      }
    }
  }
  else{                          

    /* item sizes are constant, and items belonging to a
     * given message are always stored contiguously in send_data
     */

    if (!sorted || (plan->nvals > nSendItems)){
      buf = outbuf = (char *)ZOLTAN_MALLOC((size_t)nSendItems * (size_t)nbytes);
      if (nSendItems && nbytes && !outbuf){
        ZOLTAN_COMM_ERROR("memory error", yo, me);
      }
    }
    else{
      /* send_data is sorted by process, and we don't skip
       * any of the data in the buffer, so we can use send_data 
       * in the alltoall call
       */
      outbuf = send_data;
    }

    for (p=0,i=0; p < nprocs; p++){

      length = 0;
     
      if (i < nSendMsgs){
        if (plan->procs_to[i] == p){    /* procs_to is sorted */
          offset = plan->starts_to[i] * nbytes;
          length = plan->lengths_to[i] * nbytes;
  
          if ((!sorted || (plan->nvals > nSendItems)) && length){
            memcpy(buf, send_data + offset, length);
            buf += length;
          }
          i++;
        }
      }

      outbufCounts[p] = length;
      if (p){
        outbufOffsets[p] = outbufOffsets[p-1] + outbufCounts[p-1];
      }
    }
  }

  /* CREATE RECEIVE BUFFER */

  sorted = 0;
  if (plan->indices_from == NULL){
    sorted = 1;
    for (i=1; i< nRecvMsgs; i++){
      if (plan->starts_from[i] < plan->starts_from[i-1]){
        sorted = 0;
        break;
      }
    }
  }

  if (sorted){
    /* Caller already expects received data to be ordered by
     * the sending process rank.
     */
    inbuf = recv_data;
  }
  else{
    inbuf = (char *)ZOLTAN_MALLOC((size_t)(plan->total_recv_size) * (size_t)nbytes);
    if (plan->total_recv_size && nbytes && !inbuf){
      ZOLTAN_COMM_ERROR("memory error", yo, me);
    }
  }

  for (p=0, i=0; p < nprocs; p++){
    length = 0;

    if (i < nRecvMsgs){
      if (plan->procs_from[i] == p){
  
        if (plan->sizes == NULL){
          length = plan->lengths_from[i] * nbytes;
        }
        else{
          length = plan->sizes_from[i] * nbytes;
        }  
        i++;
      }
    }

    inbufCounts[p] = length;
    if (p){
      inbufOffsets[p] = inbufOffsets[p-1] + inbufCounts[p-1];
    }
  }

  /* EXCHANGE DATA */

  rc = MPI_Alltoallv(outbuf, outbufCounts, outbufOffsets, MPI_BYTE,
                     inbuf, inbufCounts, inbufOffsets, MPI_BYTE,
                     plan->comm);

  if (outbuf != send_data){
    ZOLTAN_FREE(&outbuf);
  }
  ZOLTAN_FREE(&outbufCounts);
  ZOLTAN_FREE(&outbufOffsets);


  /* WRITE RECEIVED DATA INTO USER'S BUFFER WHERE IT'S EXPECTED */

  if (!sorted){

    buf = inbuf;

    if (plan->sizes == NULL){

      /* each item in each message is nbytes long */

      if (plan->indices_from == NULL){
        for (i=0; i < nRecvMsgs; i++){
          offset = plan->starts_from[i] * nbytes;
          length = plan->lengths_from[i] * nbytes;
          memcpy(recv_data + offset, buf, length);
          buf += length;
        }
      }
      else{
        for (i=0,k=0; i < nRecvMsgs; i++){

          for (j=0; j < plan->lengths_from[i]; j++,k++){
            offset = plan->indices_from[k] * nbytes;
            memcpy(recv_data + offset, buf, nbytes);
            buf += nbytes;
          }
        }
      }
    }
    else{  /* (sizes!=NULL) && (indices_from!=NULL) not allowed by Zoltan_Comm_Resize */

      /* items can be different sizes */

      for (i=0; i < nRecvMsgs; i++){
        offset = plan->starts_from_ptr[i] * nbytes;
        length = plan->sizes_from[i] * nbytes;
        memcpy(recv_data + offset, buf, length);
        buf += length;
      }
    }

    ZOLTAN_FREE(&inbuf);
  }

  ZOLTAN_FREE(&inbufCounts);
  ZOLTAN_FREE(&inbufOffsets);

  return ZOLTAN_OK;
}
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
