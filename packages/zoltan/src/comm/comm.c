/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

/* irregular communiction functions for C and F77 - Apr 1998

   Steve Plimpton, MS 1111, Dept 9221, Sandia National Labs  (505) 845-7873
   sjplimp@cs.sandia.gov
*/

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "comm_const.h"
#include "all_allo_const.h"

#define MAX(A,B) ((A) > (B)) ? (A) : (B)
/* KDDKDD */
static void debug_sort2_int_int(int n, int ra[], int rb[]);

/* ------------------------------------------------------------------- */
/* perform irregular communication using a pre-computed plan */

void LB_Comm_Do(struct Comm_Obj *plan,          /* plan from create_comm */
	     char *sendbuf,                  /* list of datums to send */
	     int nsize,                      /* size in bytes of each datum */
	     char *recvbuf)                  /* space for datums to recv */

{
  char *yo = "LB_Comm_Do";
  int i,isend,irecv,index,send_offset,recv_offset;
  char *buf;

/* post all receives */

  recv_offset = 0;
  for (irecv = 0; irecv < plan->nrecv; irecv++) {
    MPI_Irecv(&recvbuf[recv_offset],nsize*plan->lengths_from[irecv],MPI_BYTE,
	      plan->procs_from[irecv],0,plan->comm,&plan->request[irecv]);
    recv_offset += nsize * plan->lengths_from[irecv];
  }

/* malloc buf for largest send */

  buf = (char *) LB_MALLOC(nsize*plan->nsendmax*sizeof(char));
  if ((nsize*plan->nsendmax > 0) && !buf) {
      fprintf(stderr, "Error in %s: Insufficient memory\n", yo);
      exit(-1);
  }

/* send each message, packing buf with needed datums */

  index = 0;
  for (isend = 0; isend < plan->nsend; isend++) {
    send_offset = 0;
    for (i = 0; i < plan->lengths_to[isend]; i++) {
      memcpy(&buf[send_offset],
	     &sendbuf[nsize*plan->indices_to[index++]],nsize);
      send_offset += nsize;
    }
    MPI_Send(buf,nsize*plan->lengths_to[isend],MPI_BYTE,
	     plan->procs_to[isend],0,plan->comm);
  }       

/* copy self data directly from sendbuf to recvbuf */

  if (plan->nself)
    for (i = 0; i < plan->lengths_to[plan->nsend]; i++) {
      memcpy(&recvbuf[recv_offset],
	     &sendbuf[nsize*plan->indices_to[index++]],nsize);
      recv_offset += nsize;
    }

/* free temporary send buffer */

  LB_FREE(&buf);

/* wait on all incoming messages */

  if (plan->nrecv) MPI_Waitall(plan->nrecv,plan->request,plan->status);
}


/* ------------------------------------------------------------------- */
/* create an irregular communication plan */

struct Comm_Obj *LB_Comm_Create(
  int n,                       /* # of datums */
  int *proclist,               /* which proc each datum is to be sent to */
  MPI_Comm original_comm,      /* communicator for all participating procs */
  int *nreturn)                /* # of datums I will recv including self */

{
  char *yo = "LB_Comm_Create";
  MPI_Comm new_comm;
  struct Comm_Obj *plan;
  int me,nprocs;
  int i,iproc,isend,nsend,nrecv,nself,nsendmax;
  int *procs_to,*lengths_to,*indices_to;
  int *procs_from,*lengths_from;
  MPI_Request *request;
  MPI_Status *status;
  int *list,*counts;

/* create new MPI communicator for plan */

  MPI_Comm_dup(original_comm,&new_comm);
  MPI_Comm_rank(new_comm,&me);
  MPI_Comm_size(new_comm,&nprocs);

/* allocate plan and work vectors */

  plan = (struct Comm_Obj *) LB_MALLOC(sizeof(struct Comm_Obj));
  if (!plan) {
      fprintf(stderr, "[%d] Error in %s: Insufficient memory (plan)\n",
              me, yo);
      exit(-1);
  }
  list = (int *) LB_MALLOC(nprocs*sizeof(int));
  if (!list) {
      fprintf(stderr, "[%d] Error in %s: Insufficient memory (list)\n",
              me, yo);
      exit(-1);
  }
  counts = (int *) LB_MALLOC(nprocs*sizeof(int));
  if (!counts) {
      fprintf(stderr, "[%d] Error in %s: Insufficient memory (counts)\n",
              me, yo);
      exit(-1);
  }

/* nrecv = # of messages I receive, not including self
   nself = 0 if no data for self, 1 if there is */

  for (i = 0; i < nprocs; i++) {
    list[i] = 0;
    counts[i] = 1;
  }

  for (i = 0; i < n; i++) list[proclist[i]] = 1;

  MPI_Reduce_scatter(list,&nrecv,counts,MPI_INT,MPI_SUM,new_comm);

  nself = 0;
  if (list[me]) nself = 1;
  if (nself) nrecv--;

/* storage for recv info, not including self */

  procs_from = (int *) LB_MALLOC(nrecv*sizeof(int));
  if ((nrecv > 0) && !procs_from) {
      fprintf(stderr, "[%d] Error in %s: Insufficient memory (procs_from)\n",
              me, yo);
      exit(-1);
  }
  lengths_from = (int *) LB_MALLOC(nrecv*sizeof(int));
  if ((nrecv > 0) && !lengths_from) {
      fprintf(stderr, "[%d] Error in %s: Insufficient memory (lengths_from)\n",
              me, yo);
      exit(-1);
  }
  request = (MPI_Request *) LB_MALLOC(nrecv*sizeof(MPI_Request));
  if ((nrecv > 0) && !request) {
      fprintf(stderr, "[%d] Error in %s: Insufficient memory (request)\n",
              me, yo);
      exit(-1);
  }
  status = (MPI_Status *) LB_MALLOC(nrecv*sizeof(MPI_Status));
  if ((nrecv > 0) && !status) {
      fprintf(stderr, "[%d] Error in %s: Insufficient memory (status)\n",
              me, yo);
      exit(-1);
  }

/* nsend = # of messages I send, not including self */

  for (i = 0; i < nprocs; i++) list[i] = 0;

  for (i = 0; i < n; i++) list[proclist[i]]++;

  nsend = 0;
  for (i = 0; i < nprocs; i++) if (list[i] > 0) nsend++;
  if (nself) nsend--;

/* storage for send info, including self */

  procs_to = (int *) LB_MALLOC((nsend+nself)*sizeof(int));
  if ((nsend+nself > 0) && !procs_to) {
      fprintf(stderr, "[%d] Error in %s: Insufficient memory (procs_to)\n",
              me, yo);
      exit(-1);
  }
  lengths_to = (int *) LB_MALLOC((nsend+nself)*sizeof(int));
  if ((nsend+nself > 0) && !lengths_to) {
      fprintf(stderr, "[%d] Error in %s: Insufficient memory (lengths_to)\n",
              me, yo);
      exit(-1);
  }
  indices_to = (int *) LB_MALLOC(n*sizeof(int));
  if ((n > 0) && !indices_to) {
      fprintf(stderr, "[%d] Error in %s: Insufficient memory (indices_to)\n",
              me, yo);
      exit(-1);
  }


/* set send info in procs_to and lengths_to, including self
   each proc begins with iproc > me, and continues until iproc = me
   store pointer to sends in list for later use in setting indices_to */

  iproc = me;
  isend = 0;
  for (i = 0; i < nprocs; i++) {
    iproc++;
    if (iproc == nprocs) iproc = 0;
    if (list[iproc] > 0) {
      procs_to[isend] = iproc;
      lengths_to[isend] = list[iproc];
      list[iproc] = isend;
      isend++;
    }
  }

/* tell receivers how many I'm sending, not including self
   nsendmax = largest # of datums I send */

  nsendmax = 0;
  for (i = 0; i < nsend; i++) {
    MPI_Send(&lengths_to[i],1,MPI_INT,procs_to[i],0,new_comm);
    nsendmax = MAX(nsendmax,lengths_to[i]);
  }

/* receive sizes and sources of incoming data, not including self
   nreturn = total # of datums I recv, including self */

  *nreturn = 0;
  for (i = 0; i < nrecv; i++) {
    MPI_Recv(&lengths_from[i],1,MPI_INT,MPI_ANY_SOURCE,0,new_comm,status);
    procs_from[i] = status->MPI_SOURCE;
    *nreturn += lengths_from[i];
  }

  if (nself) *nreturn += lengths_to[nsend];

/* barrier to insure all my MPI_ANY_SOURCE messages are received
   else some procs could proceed to LB_Comm_Do and start sending to me */

  MPI_Barrier(new_comm);

/* KDDKDD  sort the procs_from array in ascending order.
           also reorder lengths_from.  
           do this ONLY to remove non-determinism; non-determinism
           may cause different results due to different order of
           migration. 

  debug_sort2_int_int(nrecv, procs_from-1, lengths_from-1);
*/

/* setup indices_to, including self
   counts = current offset into indices_to for each proc I send to */

  counts[0] = 0;
  for (i = 1; i < nsend+nself; i++) counts[i] = counts[i-1] + lengths_to[i-1];

  for (i = 0; i < n; i++) {
    isend = list[proclist[i]];
    indices_to[counts[isend]++] = i;
  }

/* free work vectors */

  LB_FREE(&counts);
  LB_FREE(&list);
    
/* initialize plan and return it */

  plan->nsend = nsend;
  plan->nrecv = nrecv;
  plan->nself = nself;
  plan->nsendmax = nsendmax;

  plan->procs_to = procs_to;
  plan->lengths_to = lengths_to;
  plan->indices_to = indices_to;
  plan->procs_from = procs_from;
  plan->lengths_from = lengths_from;
  plan->request = request;
  plan->status = status;

  plan->comm = new_comm;

  return plan;
}


/* ------------------------------------------------------------------- */
/* free all memory associated with irregular communication plan */

void LB_Comm_Destroy(struct Comm_Obj **plan)

{
/* free MPI communicator */

  MPI_Comm_free(&((*plan)->comm));

/* free internal arrays */

  LB_FREE(&((*plan)->procs_to));
  LB_FREE(&((*plan)->procs_from));
  LB_FREE(&((*plan)->lengths_to));
  LB_FREE(&((*plan)->lengths_from));
  LB_FREE(&((*plan)->indices_to));
  LB_FREE(&((*plan)->request));
  LB_FREE(&((*plan)->status));

/* free plan itself */

  LB_FREE(plan);
}


/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/* F77 wrappers on 3 comm functions */
/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */

/* F77 syntax:

     double precision plan
     integer n,nreturn,nsize
     integer proclist(*)
     dimension sendbuf(*),recvbuf(*)

     call LB_Comm_Create(n,proclist,MPI_COMM_WORLD,nreturn,plan)
     ...
     call LB_Comm_Do(plan,sendbuf,nsize,recvbuf)
     ...
     call LB_Comm_Destroy(plan)
*/


/* ------------------------------------------------------------------- */
/* F77 wrapper on LB_Comm_Do */

void LB_Comm_Do_(struct Comm_Obj **plan, char *sendbuf, int *nsize, char *recvbuf)

{
  LB_Comm_Do(*plan,sendbuf,*nsize,recvbuf);
}

/* ------------------------------------------------------------------- */
/* F77 wrapper on LB_Comm_Create */

void LB_Comm_Create_(int *n, int *proclist, MPI_Comm *original_comm,
		  int *nreturn, struct Comm_Obj **plan)

{
  *plan = LB_Comm_Create(*n,proclist,*original_comm,nreturn);
}

/* ------------------------------------------------------------------- */
/* F77 wrapper on LB_Comm_Destroy */

void LB_Comm_Destroy_(struct Comm_Obj **plan)

{
  LB_Comm_Destroy(plan);
}

/* KDDKDDKDDKDDKDDKDDKDDKDDKDDKDDKDDKDDKDDKDDKDDKDDKDDKDDKDDKDDKDDKDDKDD */

/* sorting routine to sort messages to be received by their processor
 * number.
 * this is used only to test whether non-deterministic order of messages
 * is influencing results.
 * otherwise, this sort is extraneous and unnecessary.
 */


static void debug_sort2_int_int(int n, int ra[], int rb[])
{
/*
 *       Numerical Recipies in C source code
 *       modified to have first argument an integer array (JS)
 *
 *       Sorts the array ra[1,..,n] in ascending numerical order using heapsort
 *       algorithm, while making the corresponding rearrangement of the
 *       array rb[1,..,n].
 */

  int l,j,ir,i;
  int rra;
  int rrb;

  /*
   *  No need to sort if one or fewer items.
   */
  if (n <= 1) return;

  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra=ra[--l];
      rrb=rb[l];
    } else {
      rra=ra[ir];
      rrb=rb[ir];
      ra[ir]=ra[1];
      rb[ir]=rb[1];
      if (--ir == 1) {
        ra[1]=rra;
        rb[1]=rrb;
        return;
      }
    }
    i=l;
    j=l << 1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
        ra[i]=ra[j];
        rb[i]=rb[j];
        j += (i=j);
      }
      else j=ir+1;
    }
    ra[i]=rra;
    rb[i]=rrb;
  }
}
