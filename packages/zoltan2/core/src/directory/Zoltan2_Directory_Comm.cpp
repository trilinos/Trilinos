/*
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan2 Directory for Load-balancing, Partitioning, Ordering and Coloring
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
 * 3. Neither the name of the Corporation nor the names of theremove_local
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

#include "Zoltan2_Directory_Comm.hpp"
#include <stdexcept>
#include <memory>

namespace Zoltan2 {

void Zoltan2_Directory_Plan::getInvertedValues(Zoltan2_Directory_Plan * from) {
  total_recv_size = 0;
  for (int i = 0; i < from->nsends + from->self_msg; i++) {
    total_recv_size += from->lengths_to[i];
  }

  max_send_size = 0;
  for (int i = 0; i < from->nrecvs; i++) {
	  if (from->lengths_from[i] > max_send_size) {
	    max_send_size = from->lengths_from[i];
    }
  }

  nvals         = from->nvals_recv;
  nvals_recv    = from->nvals;
  lengths_to    = from->lengths_from;
  procs_to      = from->procs_from;
  indices_to    = from->indices_from;

  starts_to     = from->starts_from;
  lengths_from  = from->lengths_to;
  procs_from    = from->procs_to;
  indices_from  = from->indices_to;

  starts_from   = from->starts_to;
  nrecvs        = from->nsends;
  nsends        = from->nrecvs;
  self_msg      = from->self_msg;
  comm          = from->comm;
}


void Zoltan2_Directory_Plan::print(const std::string& headerMessage) const {

  #define PRINT_VECTOR(v)                        \
  if(v != Teuchos::null) {                       \
    std::cout << "  " << #v << " ";              \
    for(Teuchos::ArrayRCP<int>::size_type n = 0; n < v.size(); ++n) { \
      std::cout << v[n] << " ";                  \
    }                                            \
    std::cout << std::endl;                      \
  }

  #define PRINT_VAL(val) std::cout << "  " << #val << ": " << val << std::endl;

  for(int proc = 0; proc < comm->getSize(); ++proc) {
    comm->barrier();
    if(proc == comm->getRank()) {

      std::cout << "Rank " << proc << " " << headerMessage << std::endl;

      PRINT_VECTOR(procs_to)
      PRINT_VECTOR(procs_from)
      PRINT_VECTOR(lengths_to)
      PRINT_VECTOR(lengths_from)
      PRINT_VECTOR(starts_to)
      PRINT_VECTOR(starts_from)
      PRINT_VECTOR(indices_to)
      PRINT_VECTOR(indices_from)
      PRINT_VECTOR(sizes)
      PRINT_VECTOR(sizes_to)
      PRINT_VECTOR(sizes_from)
      PRINT_VECTOR(starts_to_ptr)
      PRINT_VECTOR(starts_from_ptr)
      PRINT_VECTOR(indices_to_ptr)
      PRINT_VECTOR(indices_from_ptr)

      PRINT_VAL(nvals)
      PRINT_VAL(nvals_recv)
      PRINT_VAL(nrecvs)
      PRINT_VAL(nsends)
      PRINT_VAL(self_msg)
      PRINT_VAL(max_send_size)
      PRINT_VAL(total_recv_size)
      PRINT_VAL(maxed_recvs)
    }
  }
  comm->barrier();
}

Zoltan2_Directory_Comm::Zoltan2_Directory_Comm(
  int       nvals,	                    	/* number of values I currently own */
  const Teuchos::ArrayRCP<int>  &assign,  /* processor assignment for all my values */
  Teuchos::RCP<const Teuchos::Comm<int> > comm,               /* communicator */
  int       tag) :  	                               /* message tag I can use */
  comm_(comm),
  plan_forward(NULL)
{
  if (comm == Teuchos::null){
    throw std::logic_error("Invalid communicator: MPI_COMM_NULL.");
  }

  int my_proc = comm->getRank();		/* my processor tag in communicator */
  int nprocs = comm->getSize();		/* number of  processors in communicator */

  /* First check to see if items are grouped by processor with no gaps. */
  /* If so, indices_to should be NULL (= identity) */

  /* Make data structures that will allow me to traverse arrays quickly. */
  /* Begin by determining number of objects I'm sending to each processor. */
  Teuchos::ArrayRCP<int> starts(new int[nprocs + 1], 0, nprocs + 1, true);
  for(int n = 0; n < starts.size(); ++n) {
    starts[n] = 0;
  }

  /* Note: Negative assign value means ignore item. */
  /* Non-trailing negatives mean data not packed so need send_buf. */
  /* Could (but don't) allow negatives between processor blocks w/o buf. */
  int nactive = 0;	/* number of values to remap */

  int no_send_buff = 1;	/* is data nicely grouped by processor? */

  int prev_proc = nprocs;	/* processor on previous loop pass */

  for (int i = 0; i < nvals; i++) {
    int proc = assign[i];
    if (no_send_buff && proc != prev_proc) { /* Checks if blocked by proc */
      if (proc >= 0 && (starts[proc] || prev_proc < 0)) {
        no_send_buff = 0;
      }
      else {
        prev_proc = proc;
      }
    }
    if (proc >= 0) {
      ++starts[proc];
      ++nactive;
    }
  }

  int self_msg = (starts[my_proc] != 0);	/* do I have data for myself? */

  Teuchos::ArrayRCP<int> lengths_to;   /* lengths I'll send to */
  Teuchos::ArrayRCP<int> procs_to;     /* processors I'll send to */
  Teuchos::ArrayRCP<int> starts_to;	   /* where in list my sends begin */
  Teuchos::ArrayRCP<int> indices_to;	 /* local_id values I'll be sending */

  int max_send_size = 0;	       /* size of longest message I send */
  int nsends = 0;		             /* # procs I'll send to (including self) */
  int nrecvs = 0;		             /* # procs I'll recv from (including self) */

  if (no_send_buff) {
    /* Grouped by processor.  Array indices_to can be NULL (= identity) */
    nsends = 0;
    for (int i = 0; i < nprocs; i++) {
      if (starts[i] != 0) ++nsends;
    }

    lengths_to.resize(nsends);
    starts_to.resize(nsends);
    procs_to.resize(nsends);

    int index = 0;		/* index into list of objects */
    /* Note that procs_to is in the order the data was passed in. */
    for (int i = 0; i < nsends; i++) {
      starts_to[i] = index;
      int proc = assign[index];
      procs_to[i] = proc;
      index += starts[proc];
    }

    /* Now sort the outgoing procs. */
    /* This keeps recvs deterministic if I ever invert communication */
    /* It also allows for better balance of traffic in comm_do */
    sort_ints(procs_to, starts_to);

    max_send_size = 0;
    for (int i = 0; i < nsends; i++) {
      int proc = procs_to[i];
      lengths_to[i] = starts[proc];
      if (proc != my_proc && lengths_to[i] > max_send_size) {
        max_send_size = lengths_to[i];
      }
    }
  }
  else {	/* Not grouped by processor.  More complex data structures. */
    /* Sum starts values to be offsets into indices_to array. */
    nsends = (starts[0] != 0);
    for (int i = 1; i < nprocs; i++) {
      if (starts[i] != 0)
        ++nsends;
      starts[i] += starts[i - 1];
    }

    for (int i = nprocs - 1; i; i--) {
      starts[i] = starts[i - 1];
    }

    starts[0] = 0;

    indices_to = (nactive > 0) ?
      Teuchos::arcp(new int[nactive], 0, nactive, true) : Teuchos::null;

    for (int i = 0; i < nvals; i++) {
      int proc = assign[i];
      if (proc >= 0) {
        indices_to[starts[proc]] = i;
        ++starts[proc];
      }
    }

    /* Indices_to array now has the data in clumps for each processor. */
    /* Now reconstruct starts array to index into indices_to. */
    for (int i = nprocs - 1; i; i--) {
      starts[i] = starts[i - 1];
    }
    starts[0] = 0;
    starts[nprocs] = nactive;

    /* Construct lengths_to, starts_to and procs_to arrays. */
    /* Note: If indices_to is needed, procs are in increasing order */
    lengths_to.resize(nsends);
    starts_to.resize(nsends);
    procs_to.resize(nsends);

    int j = 0;
    max_send_size = 0;
    for (int i = 0; i < nprocs; i++) {
      if (starts[i + 1] != starts[i]) {
        starts_to[j] = starts[i];
        lengths_to[j] = starts[i + 1] - starts[i];
        if (i != my_proc && lengths_to[j] > max_send_size) {
          max_send_size = lengths_to[j];
        }
        procs_to[j] = i;
        j++;
      }
    }
  }

  /* Now change nsends to count only non-self messages */
  nsends -= self_msg;

  /* Determine how many messages & what length I'll receive. */
  Teuchos::ArrayRCP<int> lengths_from; /* lengths of messages I'll receive */
  Teuchos::ArrayRCP<int> procs_from; /* processors I'll receive from */
  int out_of_mem = 0; // TODO refactor this bit

  int comm_flag = invert_map(lengths_to, procs_to, nsends, self_msg,
    lengths_from, procs_from, &nrecvs, my_proc, nprocs,
    out_of_mem, tag, comm);

  /* pointers for where to put recv data */
  Teuchos::ArrayRCP<int> starts_from(
    new int[nrecvs + self_msg], 0, nrecvs + self_msg, true);
  int j = 0;
  for (int i = 0; i < nrecvs + self_msg; i++) {
    starts_from[i] = j;
    j += lengths_from[i];
  }

  if (comm_flag != 0) {
    throw std::logic_error("Failed to construct Zoltan2_Directory_Comm");
  }

  int total_recv_size = 0;  /* total size of messages I recv */
  for (int i = 0; i < nrecvs + self_msg; i++) {
    total_recv_size += lengths_from[i];
  }

  plan_forward = new Zoltan2_Directory_Plan;
  plan_forward->lengths_to = lengths_to;
  plan_forward->starts_to = starts_to;
  plan_forward->procs_to = procs_to;
  plan_forward->indices_to = indices_to;
  plan_forward->lengths_from = lengths_from;
  plan_forward->starts_from = starts_from;
  plan_forward->procs_from = procs_from;
  plan_forward->nvals = nvals;
  plan_forward->nvals_recv = total_recv_size;
  plan_forward->nrecvs = nrecvs;
  plan_forward->nsends = nsends;
  plan_forward->self_msg = self_msg;
  plan_forward->max_send_size = max_send_size;
  plan_forward->total_recv_size = total_recv_size;
  plan_forward->maxed_recvs = 0;
  plan_forward->comm = comm;

  if (MPI_RECV_LIMIT > 0) {

    throw std::logic_error("UNTESTED COMM 2"); // needs unit testing

    /* If we have a limit to the number of posted receives we are allowed,
    ** and our plan has exceeded that, then switch to an MPI_Alltoallv so
    ** that we will have fewer receives posted when we do the communication.
    */
    int global_nrecvs;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &nrecvs, &global_nrecvs);
    if (global_nrecvs > MPI_RECV_LIMIT){
      plan_forward->maxed_recvs = 1;
    }
  }

  if (plan_forward->maxed_recvs == 0) {
    // See notes in header for MPI_Request
    plan_forward->request.resize(plan_forward->nrecvs);
    plan_forward->status.resize(plan_forward->nrecvs);
    //    plan_forward->request = Teuchos::arcp(new MPI_Request[plan_forward->nrecvs], 0, plan_forward->nrecvs, true);
    //    plan_forward->status = Teuchos::arcp(new MPI_Status[plan_forward->nrecvs], 0, plan_forward->nrecvs, true);
  }

  nrec = total_recv_size;
}

Zoltan2_Directory_Comm::~Zoltan2_Directory_Comm()
{
  delete plan_forward;
}

int Zoltan2_Directory_Comm::invert_map(
  const Teuchos::ArrayRCP<int> &lengths_to,  /* number of items I'm sending */
  const Teuchos::ArrayRCP<int> &procs_to,    /* procs I send to */
  int       nsends,                    /* number of messages I'll send */
  int       self_msg,                  /* do I copy data to myself? */
  Teuchos::ArrayRCP<int> &lengths_from,      /* number of items I'm receiving */
  Teuchos::ArrayRCP<int> &procs_from,        /* procs I recv lengths from */
  int      *pnrecvs,                   /* number of messages I receive */
  int       my_proc,                   /* my processor number */
  int       nprocs,                    /* total number of processors */
  int       /* out_of_mem */,          /* tell everyone I'm out of memory? */
  int       tag,                       /* message tag I can use */
  Teuchos::RCP<const Teuchos::Comm<int> > comm) /* communicator */
{
  Teuchos::ArrayRCP<int> msg_count(new int[nprocs], 0, nprocs, true);
  Teuchos::ArrayRCP<int> counts(new int[nprocs], 0, nprocs, true);
  for(int i = 0; i < nprocs; ++i) {
    msg_count[i] = 0;
    counts[i] = 1;
  }

  for (int i = 0; i < nsends + self_msg; i++) {
    msg_count[procs_to[i]] = 1;
  }

  /*
  *  KDDKDD:  Replaced MPI_Reduce_scatter with MPI_Reduce and MPI_Scatter
  *  KDDKDD:  to avoid reported problems with MPICH 1.5.2.1.
  *  KDDKDD:  Some sort of MPI_TYPE_INDEXED error.
  *  KDDKDD:  Bug fix suggested by Clark Dohrmann and Rob Hoekstra.
  *  KDDKDD:  July 20, 2004

    MPI_Reduce_scatter((void *) msg_count, (void *) &nrecvs, counts, MPI_INT,
      MPI_SUM, comm);
  */
  Teuchos::reduceAll<int>(*comm, Teuchos::REDUCE_SUM, nprocs,
    msg_count.getRawPtr(), counts.getRawPtr());

  int nrecvs = 0;		/* number of messages I'll receive */

  // TODO: Teuchos::scatterImpl(...)
  MPI_Scatter(&(counts[0]), 1, MPI_INT, &nrecvs, 1, MPI_INT, 0, getRawComm());

  int max_nrecvs = 0;
  if (my_proc == 0) {
    for (int i=0; i < nprocs; i++) {
      if (counts[i] > max_nrecvs) {
        max_nrecvs = counts[i];
      }
    }
  }

  // TODO: Teuchos::MpiComm<Ordinal>::broadcast(...)

  MPI_Bcast(&max_nrecvs, 1, MPI_INT, 0, getRawComm());

  if(nrecvs > 0) {
    lengths_from.resize(nrecvs);   /* number of items I'm receiving */
    procs_from.resize(nrecvs);    /* processors I'll receive from  */
  }

  if (MPI_RECV_LIMIT == 0 || max_nrecvs <= MPI_RECV_LIMIT) {
    // See notes in header for MPI_Request
    std::vector<MPI_Request> req(nrecvs);
    //Teuchos::ArrayRCP<MPI_Request> req(new MPI_Request[nrecvs], 0, nrecvs, true);

    /* Note: I'm counting on having a unique tag or some of my incoming */
    /* messages might get confused with others.                         */

    // TODO: Teuchos::ireceiveImpl(...)
    for (int i=0; i < nrecvs; i++) {
      MPI_Irecv(&(lengths_from[0]) + i, 1, MPI_INT, MPI_ANY_SOURCE,
        tag, getRawComm(), &(req[i]));
    }

    // TODO: Teuchos::sendImpl(...)
    for (int i=0; i < nsends+self_msg; i++) {
      // MPI_Send takes non-const for buf (1st param)
      // Apparently in future versions this will become true const
      // Then the const_cast could be removed
      MPI_Send(const_cast<int*>(&lengths_to[i]), 1, MPI_INT, procs_to[i], tag,
        getRawComm());
    }

    for (int i=0; i < nrecvs; i++) {
      MPI_Status status;
      MPI_Wait(&(req[i]), &status);
      procs_from[i] = status.MPI_SOURCE;
    }

  }
  else { /* some large HPC machines have a limit on number of posted receives */
    throw std::logic_error("UNTESTED COMM 3"); // needs unit testing

    Teuchos::ArrayRCP<int> sendbuf(new int[nprocs], 0, nprocs, true);
    Teuchos::ArrayRCP<int> recvbuf(new int[nprocs], 0, nprocs, true);

    for (int i=0; i < nsends + self_msg; i++) {
      sendbuf[procs_to[i]] = lengths_to[i];
    }

    // TODO: Teuchos
    MPI_Alltoall(&(sendbuf[0]), 1,  MPI_INT, &(recvbuf[0]), 1, MPI_INT,
      getRawComm());

    for (int i=0, j=0; i < nprocs; i++) {
      if (recvbuf[i] > 0){
        lengths_from[j] = recvbuf[i];
        procs_from[j] = i;
        if (++j == nrecvs) {
          break;
        }
      }
    }
  }

  /* Sort recv lists to keep execution deterministic (e.g. for debugging) */
  sort_ints(procs_from, lengths_from);

  *pnrecvs = nrecvs - self_msg;    /* Only return number of true messages */

  return 0;
}

int Zoltan2_Directory_Comm::sort_ints(
  Teuchos::ArrayRCP<int> &vals_sort,     /* values to be sorted */
  Teuchos::ArrayRCP<int> &vals_other)    /* other array to be reordered w/ sort */
{
  // TODO: Check - perhaps we can skip all of these for efficiency
  if (vals_sort == Teuchos::null || vals_sort.size() == 0) {
    return 1;
  }
  if (vals_other == Teuchos::null || vals_other.size() == 0) {
    return 1;
  }
  if (vals_sort == Teuchos::null || vals_sort.size() == 1) {
    return 0;           /* fastest way to sort 1 item is to return */
  }

  /* find largest value (sort sometimes used for non processor lists) */
  int already_sorted = 1;  /* flag indicating whether vals_sort is
                              already sorted; can exit early and skip
                              memory allocations if it is.  */
  int top = vals_sort[0]; /* largest integer to sort, smallest is assumed 0 */
  for (Teuchos::ArrayRCP<int>::size_type i = 1; i < vals_sort.size(); i++) {
    if (vals_sort[i-1] > vals_sort[i]) {
      already_sorted = 0;
    }
    if (top < vals_sort[i]) {
      top = vals_sort[i];
    }
  }

  if (already_sorted) {
    return 0;
  }

  Teuchos::ArrayRCP<int> store(new int[top+2], 0, top+2, true);
  for(int n = 0; n < store.size(); ++n) {
    store[n] = 0;
  }

  Teuchos::ArrayRCP<int> copy_sort(new int[vals_sort.size()], 0, vals_sort.size(), true);
  for(Teuchos::ArrayRCP<int>::size_type n = 0; n < copy_sort.size(); ++n) {
    copy_sort[n] = vals_sort[n]; // TODO - use deepCopy method?
  }

  Teuchos::ArrayRCP<int> copy_other(new int[vals_other.size()], 0, vals_other.size(), true);
  for(Teuchos::ArrayRCP<int>::size_type n = 0; n < copy_other.size(); ++n) {
    copy_other[n] = vals_other[n]; // TODO - use deepCopy method?
  }

  // TODO: May want to modernize this ptr handling - however I didn't want
  // to introduce inefficiencies so for now have kept the original structure
  int *p = &(store[1]);
  for (Teuchos::ArrayRCP<int>::size_type i = 0; i < vals_sort.size(); i++) {
    p[copy_sort[i]]++;                /* count number of occurances */
  }

  for (int i = 1; i < top+1; i++) {
    p[i] += p[i-1];                   /* compute partial sums */
  }
                                   /* assert: p[top] = nvals */
  p = &(store[0]);                    /* effectively shifts down by one */
  for (Teuchos::ArrayRCP<int>::size_type i = 0; i < vals_sort.size(); i++) {
    vals_sort[p[copy_sort[i]]] = copy_sort[i];
    vals_other[p[copy_sort[i]]] = copy_other[i];
    ++p[copy_sort[i]];
  }

  return 0;
}

int Zoltan2_Directory_Comm::do_forward(
  int tag,			                         /* message tag for communicating */
  const Teuchos::ArrayRCP<char> &send_data,		 /* array of data I currently own */
  int nbytes,                            /* msg size */
  Teuchos::ArrayRCP<char> &recv_data)		 /* array of data I'll own after comm */
{
  int status = 0;

  if (!plan_forward->maxed_recvs) {
    status = do_post (plan_forward, tag, send_data, nbytes, recv_data);
    if (status == 0) {
      status = do_wait (plan_forward, tag, send_data, nbytes, recv_data);
    }
  }
  else {
    status = do_all_to_all(plan_forward, send_data, nbytes, recv_data);
  }

  return status;
}

int Zoltan2_Directory_Comm::do_post(
  Zoltan2_Directory_Plan *plan,          /* communication data structure  */
  int tag,			                         /* message tag for communicating */
  const Teuchos::ArrayRCP<char> &send_data,		 /* array of data I currently own */
  int nbytes,                            /* msg size */
  Teuchos::ArrayRCP<char> &recv_data)     /* array of data I'll own after comm */
{
  /* Check input parameters */
  if (!plan) {
    throw std::logic_error("Communication plan = NULL");
  }

  /* If not point to point, currently we do synchroneous communications */
  if (plan->maxed_recvs) {
    throw std::logic_error("UNTESTED COMM 4"); // needs unit testing
    return do_all_to_all(plan, send_data, nbytes, recv_data);
  }

  int my_proc = plan->comm->getRank();		/* processor ID */
  if ((plan->nsends + plan->self_msg) && !send_data.size()) {
    throw std::logic_error("UNTESTED COMM 5"); // needs unit testing
    size_t sum = 0;
    if (plan->sizes_to.size()) {   /* Not an error if all sizes_to == 0 */
      for (int i = 0; i < (plan->nsends + plan->self_msg); i++) {
        sum += plan->sizes_to[i];
      }
    }
    if (!plan->sizes_to.size() || (plan->sizes_to.size() && sum)) {
      throw std::logic_error("nsends not zero, but send_data = NULL");
    }
  }
  if ((plan->nrecvs + plan->self_msg) && recv_data == Teuchos::null) {
    throw std::logic_error("UNTESTED COMM 6"); // needs unit testing
    size_t sum = 0;
    if (plan->sizes_from != Teuchos::null)   /* Not an error if all sizes_from == 0 */
      for (int i = 0; i < (plan->nrecvs + plan->self_msg); i++)
        sum += plan->sizes_from[i];
    if (!plan->sizes_from.size() || (plan->sizes_from.size() && sum)) {
      throw std::logic_error("nrecvs not zero, but recv_data = NULL");
    }
  }

  /* Post irecvs */
  if (plan->indices_from == Teuchos::null) {
    /* Data can go directly into user space. */
    plan->recv_buff = recv_data;
  }
  else {			/* Need to buffer receive to reorder */
    size_t rsize = (size_t) (plan->total_recv_size) * (size_t) nbytes;
    plan->recv_buff = Teuchos::arcp(new char[rsize], 0, rsize, true);
  }

  size_t self_recv_address = 0;  /* where in recv_data self info starts */
  if (!plan->using_sizes) {	/* All data the same size */
    int k = 0;
    for (int i = 0; i < plan->nrecvs + plan->self_msg; i++) {
      if (plan->procs_from[i] != my_proc) {
        // TODO: Teuchos::ireceiveImpl(...)
        MPI_Irecv((void *)
          &plan->getRecvBuff().getRawPtr()[(size_t)(plan->starts_from[i])*(size_t)nbytes],
          plan->lengths_from[i] * nbytes,
          (MPI_Datatype) MPI_BYTE, plan->procs_from[i], tag,
          getRawComm(), &plan->request[k]);
        k++;
      }
      else {
        self_recv_address = (size_t)(plan->starts_from[i]) * (size_t)nbytes;
      }
    }
  }
  else {			/* Data of varying sizes */
    int k = 0;
    for (int i = 0; i < plan->nrecvs + plan->self_msg; i++) {
      if (plan->procs_from[i] != my_proc) {
        if (plan->sizes_from[i]) {
          // TODO: Teuchos::ireceiveImpl(...)
          MPI_Irecv((void *)
            &plan->getRecvBuff().getRawPtr()[(size_t)(plan->starts_from_ptr[i])
              * (size_t)nbytes],
            plan->sizes_from[i] * nbytes,
            (MPI_Datatype) MPI_BYTE, plan->procs_from[i],
            tag, getRawComm(), &plan->request[k]);
        }
        else
          plan->request[k] = MPI_REQUEST_NULL;
        k++;
      }
      else {
        self_recv_address =
          (size_t)(plan->starts_from_ptr[i]) * (size_t)nbytes;
      }
    }
  }

  Teuchos::ArrayRCP<char> send_buff;
  if(plan->indices_to != Teuchos::null) {
    send_buff = Teuchos::arcp(new char[plan->max_send_size * nbytes], 0, plan->max_send_size * nbytes, true);
  }

  /* Barrier to ensure irecvs are posted before doing any sends. */
  /* Simultaneously see if anyone out of memory */
  int out_of_mem = 0;
  // WARNING - do not delete this without proper barrier added as replacmeent.
  // I'm refactoring memory handling so probably we won't use out_of_mem
  // in the new version but we must still preserve a barrier here or get
  // intermittent failures.
  // I'll keep the memory reduce for now since we may end up with a memory
  // handling anyways.
  int global_out_of_mem;
  Teuchos::reduceAll(*plan->comm, Teuchos::REDUCE_SUM, 1, &out_of_mem,
    &global_out_of_mem);

  /* Send out data */

  /* Scan through procs_to list to start w/ higher numbered procs */
  /* This should balance message traffic. */

  int nblocks = plan->nsends + plan->self_msg; /* # procs needing my data */
  int proc_index = 0;	/* loop counter over procs to send to */
  while (proc_index < nblocks && plan->procs_to[proc_index] < my_proc) {
    proc_index++;
  }
  if (proc_index == nblocks) {
    proc_index = 0;
  }

  if (!plan->using_sizes) {	/* Data all of same size */
    if (plan->indices_to == Teuchos::null) {	/* data already blocked by processor. */
      int self_num = 0;       /* where in send list my_proc appears */
      for (int i = proc_index, j = 0; j < nblocks; j++) {
        if (plan->procs_to[i] != my_proc) {
          // TODO: Teuchos::readySend(...)
          MPI_Rsend(
            (void *) &send_data[(size_t)(plan->starts_to[i])*(size_t)nbytes],
            plan->lengths_to[i] * nbytes, (MPI_Datatype) MPI_BYTE,
            plan->procs_to[i], tag, getRawComm());
        }
        else {
          self_num = i;
        }
        if (++i == nblocks) {
          i = 0;
        }
      }

      if (plan->self_msg) {	/* Copy data to self. */
        /* I use array+offset instead of &(array[offset]) because of
         a bug with PGI v9 */
        /* I use memmove because I'm not sure that the pointer are not
         overlapped. */

        // TODO: Refactor to make C++ and remove getRawPtr, etc
        memmove(
          plan->getRecvBuff().getRawPtr()+self_recv_address,
          send_data.getRawPtr()+(size_t)(plan->starts_to[self_num])*(size_t)nbytes,
          (size_t) (plan->lengths_to[self_num]) * (size_t) nbytes);
      }
    }
    else { /* Not blocked by processor.  Need to buffer. */
      int self_index = 0;	    /* send offset for data I'm keeping */
      int self_num = 0;       /* where in send list my_proc appears */
      for (int i = proc_index, jj = 0; jj < nblocks; jj++) {
        if (plan->procs_to[i] != my_proc) {
          /* Need to pack message first. */
          size_t offset = 0; /* offset into array I'm copying into */
          int j = plan->starts_to[i];
          for (int k = 0; k < plan->lengths_to[i]; k++) {
            memcpy(&send_buff[offset],
            &send_data[(size_t)(plan->indices_to[j++]) * (size_t)nbytes], nbytes);
            offset += nbytes;
          }
          // TODO: Teuchos::readySend(...)
          MPI_Rsend((void *) &(send_buff[0]), plan->lengths_to[i] * nbytes,
                    (MPI_Datatype) MPI_BYTE, plan->procs_to[i], tag,
                    getRawComm());
        }
        else {
          self_num = i;
          self_index = plan->starts_to[i];
        }
        if (++i == nblocks)
          i = 0;
      }

      if (plan->self_msg) {	/* Copy data to self. */
        for (int k = 0; k < plan->lengths_to[self_num]; k++) {
          memcpy(&(plan->getRecvBuff())[self_recv_address],
            &send_data[
              (size_t)(plan->indices_to[self_index++]) * (size_t)nbytes],
              nbytes);
          self_recv_address += nbytes;
        }
      }
    }
  }
  else {                       /* Data of differing sizes */
    if (plan->indices_to == Teuchos::null) {        /* data already blocked by processor. */
      int self_num = 0;       /* where in send list my_proc appears */
      for (int i = proc_index, j = 0; j < nblocks; j++) {
        if (plan->procs_to[i] != my_proc) {
          if (plan->sizes_to[i]) {
            // TODO: Teuchos::readySend(...)
            MPI_Rsend((void *) &send_data[(size_t)(
                      plan->starts_to_ptr[i]) * (size_t)nbytes],
                      plan->sizes_to[i] * nbytes,
                      (MPI_Datatype) MPI_BYTE, plan->procs_to[i],
                      tag, getRawComm());
          }
        }
        else
          self_num = i;

        if (++i == nblocks)
          i = 0;
      }
      if (plan->self_msg) {    /* Copy data to self. */
        if (plan->sizes_to[self_num]) {
          char* lrecv = &plan->getRecvBuff().getRawPtr()[self_recv_address];
          const char* lsend =
          &send_data.getRawPtr()[(size_t)(plan->starts_to_ptr[self_num]) * (size_t)nbytes];
          int sindex = plan->sizes_to[self_num], idx;
          for (idx=0; idx<nbytes; idx++) {
            memcpy(lrecv, lsend, sindex);
            lrecv += sindex;
            lsend += sindex;
          }
        }
      }
    }
    else {                     /* Not blocked by processor.  Need to buffer. */
      int self_num = 0;       /* where in send list my_proc appears */
      for (int i = proc_index, jj = 0; jj < nblocks; jj++) {
        if (plan->procs_to[i] != my_proc) {
          size_t offset = 0; /* offset into array I'm copying into */
          int j = plan->starts_to[i];
          for (int k = 0; k < plan->lengths_to[i]; k++) {
            if (plan->sizes[plan->indices_to[j]]) {
              memcpy(&send_buff[offset],
                &send_data[(size_t)(plan->indices_to_ptr[j]) * (size_t)nbytes],
                (size_t)(plan->sizes[plan->indices_to[j]]) * (size_t)nbytes);
              offset +=
                (size_t)(plan->sizes[plan->indices_to[j]]) * (size_t)nbytes;
            }
            j++;
          }
          if (plan->sizes_to[i]) {
            // TODO: Teuchos::readySend(...)
            MPI_Rsend((void *) &(send_buff[0]), plan->sizes_to[i] * nbytes,
              (MPI_Datatype) MPI_BYTE, plan->procs_to[i], tag, getRawComm());
          }
        }
        else
          self_num = i;

        if (++i == nblocks)
          i = 0;
      }
      if (plan->self_msg) {    /* Copy data to self. */
        if (plan->sizes_to[self_num]) {
          int j = plan->starts_to[self_num];
          for (int k = 0; k < plan->lengths_to[self_num]; k++) {
            int kk = plan->indices_to_ptr[j];
            char* lrecv = &(plan->getRecvBuff())[self_recv_address];
            size_t send_idx = (size_t)kk * (size_t)nbytes;
            const char* lsend = &send_data[send_idx];
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
    }
  }

  return 0;
}

int Zoltan2_Directory_Comm::do_wait(
  Zoltan2_Directory_Plan *plan,          /* communication data structure */
  int /* tag */,			                         /* message tag for communicating */
  const Teuchos::ArrayRCP<char> &/* send_data */, /* array of data I currently own */
  int nbytes,                            /* msg size */
  Teuchos::ArrayRCP<char> &recv_data)		 /* array of data I'll own after comm */
{
  /* If not point to point, currently we do synchroneous communications */
  if (plan->maxed_recvs){
    /* Do nothing */
    return 0;
  }

  int my_proc = plan->comm->getRank();		/* processor ID */

  /* Wait for messages to arrive & unpack them if necessary. */
  /* Note: since request is in plan, could wait in later routine. */
  if (plan->indices_from == Teuchos::null) {	/* No copying required */
    if (plan->nrecvs > 0) {
      // TODO: Teuchos::waitAllImpl(...)
      MPI_Waitall(plan->nrecvs, &plan->request[0], &plan->status[0]);
      //MPI_Waitall(plan->nrecvs, plan->request.getRawPtr(), plan->status.getRawPtr());
    }
  }
  else {			 	/* Need to copy into recv_data. */
    int self_num;		/* where in send list my_proc appears */
    size_t offsetDst = 0;
    if (plan->self_msg) {		/* Unpack own data before waiting */
      for (self_num = 0; self_num < plan->nrecvs + plan->self_msg;
        self_num++) {
        if (plan->procs_from[self_num] == my_proc) {
          break;
        }
      }

      if(plan->sizes_from.size()) {
        // NEW METHOD for variable sized data
        // This will NOT put the data in order but during the directory read
        // of the buffer which follows, this data gets sorted anyways since each
        // element has an index to tell where it belongs. I am not sure yet if
        // we must sort here, or if it's fine to allow the sort to happen through
        // the index value. This needs furthe discussion.
        memcpy(&recv_data[offsetDst * (size_t)nbytes],
          &(plan->getRecvBuff())[plan->starts_from_ptr[self_num] * (size_t)nbytes],
          plan->sizes_from[self_num] * (size_t)nbytes);
        offsetDst += plan->sizes_from[self_num];
      }
      else {
        int k = plan->starts_from[self_num];
        for (int j = plan->lengths_from[self_num]; j; j--) {
          memcpy(&recv_data[(size_t)(plan->indices_from[k]) * (size_t)nbytes],
            &(plan->getRecvBuff())[(size_t)k * (size_t)nbytes], nbytes);
          k++;
        }
      }
    }
    else {
      self_num = plan->nrecvs;
    }

    for (int jj = 0; jj < plan->nrecvs; jj++) {
      MPI_Status status;		/* return from Waitany */
      int index;
      MPI_Waitany(plan->nrecvs, &plan->request[0], &index, &status);
      //MPI_Waitany(plan->nrecvs, plan->request.getRawPtr(), &index, &status);

      if (index == MPI_UNDEFINED) {
        break;  /* No more receives */
      }

      if (index >= self_num) {
        index++;
      }

      if(plan->sizes_from.size()) {
        // NEW METHOD for variable sized data
        // This will NOT put the data in order but during the directory read
        // of the buffer which follows, this data gets sorted anyways since each
        // element has an index to tell where it belongs. I am not sure yet if
        // we must sort here, or if it's fine to allow the sort to happen through
        // the index value. This needs furthe discussion.
        memcpy(&recv_data[offsetDst * (size_t)nbytes],
          &plan->getRecvBuff().getRawPtr()[plan->starts_from_ptr[index] * (size_t)nbytes],
          plan->sizes_from[index] * (size_t)nbytes);
        offsetDst += plan->sizes_from[index];
      }
      else {
        int k = plan->starts_from[index];
        for (int j = plan->lengths_from[index]; j; j--) {
          memcpy(&recv_data.getRawPtr()[(size_t)(plan->indices_from[k]) * (size_t)nbytes],
                 &plan->getRecvBuff().getRawPtr()[(size_t)k * (size_t)nbytes], nbytes);
          k++;
        }
      }
    }
  }

  return 0;
}

/* Do_Post would require posting more receives than allowed on this platform.
*  We use MPI_AlltoAllv instead, which is probably implemented such that each
*  process does one receive at a time.
*/

int Zoltan2_Directory_Comm::do_all_to_all(
  Zoltan2_Directory_Plan *plan,              /* communication data structure */
  const Teuchos::ArrayRCP<char> &send_data,		     /* array of data I currently own */
  int nbytes,                                /* msg size */
  Teuchos::ArrayRCP<char> &recv_data)		 /* array of data I'll own after comm */
{
  throw std::logic_error("UNTESTED COMM 10"); // needs unit testing

  int sm = (plan->self_msg > 0) ? 1 : 0;

  int nSendMsgs = plan->nsends + sm;
  int nRecvMsgs = plan->nrecvs + sm;

  int nSendItems = 0;
  for (int i=0; i <nSendMsgs; i++) {
    nSendItems += plan->lengths_to[i];
  }
  int nRecvItems = 0;
  for (int i=0; i <nRecvMsgs; i++) {
    nRecvItems += plan->lengths_from[i];
  }

  int nprocs = plan->comm->getSize();

  Teuchos::ArrayRCP<int> outbufCounts(new int[nprocs], 0, nprocs, true);
  Teuchos::ArrayRCP<int> outbufOffsets(new int[nprocs], 0, nprocs, true);
  Teuchos::ArrayRCP<int> inbufCounts(new int[nprocs], 0, nprocs, true);
  Teuchos::ArrayRCP<int> inbufOffsets(new int[nprocs], 0, nprocs, true);

  /* The *_to fields of the plan refer to the items in the send_data buffer,
   * and how to pull out the correct items for each receiver.  The
   * *_from fields of the plan refer to the recv_data buffer.  Items
   * arrive in process rank order, and these fields tell us where to
   * put them in the recv_data buffer.
   */

  /* CREATE SEND BUFFER */

  int sorted = 0;
  if (plan->indices_to == Teuchos::null){
    sorted = 1;
    for (int i=1; i< nSendMsgs; i++){
      if (plan->starts_to[i] < plan->starts_to[i-1]){
        sorted = 0;
        break;
      }
    }
  }

  Teuchos::ArrayRCP<char> outbuf;
  Teuchos::ArrayRCP<char> inbuf;
  Teuchos::ArrayRCP<char> buf;

  if (plan->sizes_to.size()){
    /*
     * Each message contains items for a process, and each item may be
     * a different size.
     */

    int outbufLen = 0;
    for (int i = 0; i < nSendMsgs; i++){
      outbufLen += plan->sizes_to[i];
    }

    if (plan->indices_to != Teuchos::null) {
      /*
       * items are not grouped by message
       */
      buf.resize(outbufLen*nbytes);
      outbuf.resize(outbufLen*nbytes);
      char * pBufPtr = &(outbuf[0]);
      int i = 0;
      int k = 0;
      for (int p = 0; p < nprocs; p++) {

        int length = 0;

        if (i < nSendMsgs){
          if (plan->procs_to[i] == p){   /* procs_to is sorted */

            for (int j=0; j < plan->lengths_to[i]; j++,k++){
              int itemSize = plan->sizes[plan->indices_to[k]] * nbytes;
              int offset = plan->indices_to_ptr[k] * nbytes;

              memcpy(pBufPtr, &(send_data[0]) + offset, itemSize);

              pBufPtr += itemSize;
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
        buf.resize(outbufLen*nbytes);
        outbuf.resize(outbufLen*nbytes);
      }
      else{
        /* All items in send_data are being sent, and they are sorted
         * in process rank order.
         */
        // TODO: Optimize - original just set the ptr...
        for(int n = 0; n < outbufLen*nbytes; ++n) {
          outbuf[n] = send_data[n];
        }
      }

      char * pBufPtr = &(outbuf[0]);

      int i = 0;
      for (int p = 0; p < nprocs; p++) {

        int length = 0;

        if (i < nSendMsgs){
          if (plan->procs_to[i] == p){   /* procs_to is sorted */
            length = plan->sizes_to[i] * nbytes;
            int offset = plan->starts_to_ptr[i] * nbytes;

            if ((!sorted || (plan->nvals > nSendItems)) && length){
              memcpy(pBufPtr, &(send_data[0]) + offset, length);
              pBufPtr += length;
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
  else if (plan->indices_to != Teuchos::null) {
    /*
     * item sizes are constant, however the items belonging in a given
     * message may not be contiguous in send_data
     */
    buf.resize(nSendItems*nbytes);
    outbuf.resize(nSendItems*nbytes);
    char * pBufPtr = &(outbuf[0]);
    int i = 0;
    int k = 0;
    for (int p = 0; p < nprocs; p++){

      int length = 0;

      if (i < nSendMsgs){
        if (plan->procs_to[i] == p){   /* procs_to is sorted */
          for (int j=0; j < plan->lengths_to[i]; j++,k++) {
            int offset = plan->indices_to[k] * nbytes;
            memcpy(pBufPtr, &(send_data[0]) + offset, nbytes);
            pBufPtr += nbytes;
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
      buf.resize(nSendItems*nbytes);
      outbuf.resize(nSendItems*nbytes);
    }
    else{
      /* send_data is sorted by process, and we don't skip
       * any of the data in the buffer, so we can use send_data
       * in the alltoall call
       */
      // TODO: Optimize - original just set ptr
      outbuf = send_data;
    }

    char * pBufPtr = &(outbuf[0]);

    int i = 0;
    for (int p=0; p < nprocs; p++) {

      int length = 0;

      if (i < nSendMsgs){
        if (plan->procs_to[i] == p){    /* procs_to is sorted */
          int offset = plan->starts_to[i] * nbytes;
          length = plan->lengths_to[i] * nbytes;

          if ((!sorted || (plan->nvals > nSendItems)) && length){
            memcpy(pBufPtr, &(send_data[0]) + offset, length);
            pBufPtr += length;
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
  int i;
  if (plan->indices_from == Teuchos::null) {
    sorted = 1;
    for (i=1; i< nRecvMsgs; i++) {
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

    // TODO: Optimize - original just set ptr
    outbuf = send_data;
    inbuf = recv_data;
  }
  else {
    inbuf.resize(plan->total_recv_size * nbytes);
  }

  for (int p = 0; p < nprocs; p++) {
    int length = 0;

    if (i < nRecvMsgs){
      if (plan->procs_from[i] == p){

        if (!plan->using_sizes){
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
  MPI_Alltoallv(&(outbuf[0]), &(outbufCounts[0]), &(outbufOffsets[0]), MPI_BYTE,
   &(inbuf[0]), &(inbufCounts[0]), &(inbufOffsets[0]), MPI_BYTE, getRawComm());

  // TODO: Restorecan we just copy optimizations - original just set ptr
  //if (outbuf != send_data){
  //  ZOLTAN_FREE(outbuf);
  //}

  /* WRITE RECEIVED DATA INTO USER'S BUFFER WHERE IT'S EXPECTED */

  if (!sorted){

    char * pBufPtr = &(inbuf[0]);

    if (!plan->using_sizes){

      /* each item in each message is nbytes long */

      if (plan->indices_from == Teuchos::null) {
        for (i=0; i < nRecvMsgs; i++){
          int offset = plan->starts_from[i] * nbytes;
          int length = plan->lengths_from[i] * nbytes;
          memcpy(&(recv_data[0]) + offset, pBufPtr, length);
          pBufPtr += length;
        }
      }
      else{
        int k = 0;
        for (i=0; i < nRecvMsgs; i++) {

          for (int j=0; j < plan->lengths_from[i]; j++,k++){
            int offset = plan->indices_from[k] * nbytes;
            memcpy(&(recv_data[0]) + offset, pBufPtr, nbytes);
            pBufPtr += nbytes;
          }
        }
      }
    }
    else{  /* (sizes!=NULL) && (indices_from!=NULL) not allowed by Zoltan_Comm_Resize */

      /* items can be different sizes */

      for (i=0; i < nRecvMsgs; i++){
        int offset = plan->starts_from_ptr[i] * nbytes;
        int length = plan->sizes_from[i] * nbytes;
        memcpy(&(recv_data[0]) + offset, pBufPtr, length);
        pBufPtr += length;
      }
    }
  }

  return 0;
}

int Zoltan2_Directory_Comm::do_reverse(
  int tag,			                             /* message tag for communicating */
  const Teuchos::ArrayRCP<char> &send_data,	       /* array of data I currently own */
  int nbytes,                                /* msg size */
  const Teuchos::ArrayRCP<int> &sizes,
  Teuchos::ArrayRCP<char> &recv_data)	 /* array of data I'll own after reverse comm */
{
  /* create plan->plan_reverse
   */
  int status = create_reverse_plan(tag, sizes);

  // NEW METHOD
  // Set up recv_data with the proper size
  // This information is only available after the above create_reverse_plan is
  // called so we can setup the return data with proper buffer size now.
  // However should we do this here?
  size_t new_size = plan_forward->plan_reverse->total_recv_size*nbytes;
  if(new_size > 0) {
    // have to be careful with new[0] which will for example turn a
    // arrayRCP.getRawPtr() from a valid NULL read to a debug assert fail.
    recv_data = Teuchos::arcp(new char[new_size], 0, new_size, true);
  }

  if (status == 0) {

    if (plan_forward->plan_reverse->maxed_recvs) {

      throw std::logic_error("UNTESTED COMM 11"); // needs unit testing

      /* use MPI_Alltoallv to implement plan->plan_reverse, because comm_do_post
       * would post more receives that allowed on this machine
       */

      status = do_all_to_all(plan_forward->plan_reverse, send_data,
        nbytes, recv_data);
    }
    else {
      /* use post/wait which is faster when each sends to few
       */
      status = do_post(plan_forward->plan_reverse, tag, send_data,
        nbytes, recv_data);

      if (status == 0) {
        status = do_wait (plan_forward->plan_reverse, tag, send_data,
          nbytes, recv_data);
      }
    }
  }

  free_reverse_plan(plan_forward);

  return status;
}

void Zoltan2_Directory_Comm::free_reverse_plan(Zoltan2_Directory_Plan *plan)
{
  if(!plan) {
    throw std::logic_error("Plan is NULL!");
  }
  delete plan->plan_reverse;
  plan->plan_reverse = NULL;
}

int Zoltan2_Directory_Comm::create_reverse_plan(
  int tag,
  const Teuchos::ArrayRCP<int> &sizes)/* variable size of objects (if not size 0) */
{
  /* Check input parameters */
  if (!plan_forward){
    throw std::logic_error("memory error");
  }

  /* Let Zoltan_Comm_Do check the remaining parameters. */
  plan_forward->plan_reverse = new Zoltan2_Directory_Plan;
  plan_forward->plan_reverse->getInvertedValues(plan_forward);

  if (MPI_RECV_LIMIT > 0){
    /* If we have a limit to the number of posted receives we are allowed,
    ** and our plan has exceeded that, then switch to an MPI_Alltoallv so
    ** that we will have fewer receives posted when we do the communication.
    */
    int global_nsends;
    Teuchos::reduceAll<int>(*plan_forward->comm, Teuchos::REDUCE_SUM, 1,
      &plan_forward->nsends, &global_nsends);
    if (global_nsends > MPI_RECV_LIMIT){
      plan_forward->plan_reverse->maxed_recvs = 1;
    }
  }

  if (plan_forward->plan_reverse->maxed_recvs == 0) {
    // See notes in header for MPI_Request
    //   plan_forward->plan_reverse->request = Teuchos::arcp(new MPI_Request[plan_forward->plan_reverse->nrecvs], 0, plan_forward->plan_reverse->nrecvs, true);
    //   plan_forward->plan_reverse->status = Teuchos::arcp(new MPI_Status[plan_forward->plan_reverse->nrecvs], 0, plan_forward->plan_reverse->nrecvs, true);
    plan_forward->plan_reverse->request.resize(plan_forward->plan_reverse->nrecvs);
    plan_forward->plan_reverse->status.resize(plan_forward->plan_reverse->nrecvs);
  }

  int sum_recv_sizes;
  int comm_flag = resize( plan_forward->plan_reverse,
    sizes, tag, &sum_recv_sizes);

  if (comm_flag != 0) {
    return(comm_flag);
  }

  if (sum_recv_sizes != plan_forward->plan_reverse->total_recv_size){
     /* Sanity check */
     return 1;
  }

  return 0;
}

int Zoltan2_Directory_Comm::resize(
  const Teuchos::ArrayRCP<int> &sizes,	 /* size of each item I'm sending */
  int       tag,			             /* message tag I can use */
  int      *sum_recv_sizes)        /* sum of the sizes of the items I'll receive */
{
  return resize(plan_forward, sizes, tag, sum_recv_sizes);
}

int Zoltan2_Directory_Comm::resize(
  Zoltan2_Directory_Plan *plan,    /* communication plan object */
  const Teuchos::ArrayRCP<int> &sizes,	 /* size of each item I'm sending */
  int       tag,			             /* message tag I can use */
  int      *sum_recv_sizes)        /* sum of the sizes of the items I'll receive */
{
  /* If sizes vary, then I need to compute and communicate message lengths */
  /* First check if sizes array is NULL on all procs. */
  int my_proc = plan->comm->getRank();		/* my processor ID */
  int has_sizes = (sizes.size() != 0);
  int var_sizes;        /* items have variable sizes? */

  // I think we'll need to do this raw, not with Teuchos
  MPI_Allreduce(&has_sizes, &var_sizes, 1, MPI_INT, MPI_LOR, getRawComm());

  if (var_sizes && plan->indices_from != Teuchos::null) {
    // NEW METHOD
    // Allow this to run now - the below implementation is working but perhaps
    // not done correctly for other usage cases I have not considered yet.
    // throw std::logic_error("Non-blocked, variable-sized recvs not supported");
  }

  int nsends = plan->nsends; /* number of msgs I'll send */
  int nrecvs = plan->nrecvs; /* number of msgs I'll recv */
  int self_msg = plan->self_msg;

  Teuchos::ArrayRCP<int> sizes_to;
  Teuchos::ArrayRCP<int> sizes_from;
  Teuchos::ArrayRCP<int> starts_to_ptr;
  Teuchos::ArrayRCP<int> starts_from_ptr;
  Teuchos::ArrayRCP<int> indices_to_ptr;
  Teuchos::ArrayRCP<int> indices_from_ptr;

  if (!var_sizes) { /* Easy case.  Size = length */
    plan->total_recv_size = 0;
    for (int i = 0; i < nrecvs + self_msg; i++) {
      plan->total_recv_size += plan->lengths_from[i];
    }

    plan->max_send_size = 0;
    for (int i = 0; i < nsends + self_msg; i++) {
      if (plan->procs_to[i] != my_proc &&
          plan->lengths_to[i] > plan->max_send_size) {
        plan->max_send_size = plan->lengths_to[i];
      }
    }
  }
  else {		/* Need to actually compute message sizes */

    // TODO Investigate purpose of the +1 in the old code. Is that used?

    //  OLD CODE
    //  plan->sizes.resize(plan->nvals + 1);
    //  for (int i = 0; i < plan->nvals; i++) {
    //    plan->sizes[i] = sizes[i];
    //  }

    // NEW CODE
    plan->sizes = sizes; // can we just copy?
    plan->using_sizes = true;

    if(nsends + self_msg > 0) {
      sizes_to = Teuchos::arcp(
        new int[nsends + self_msg], 0, nsends + self_msg, true);
      for(int n = 0; n < sizes_to.size(); ++n) {
        sizes_to[n] = 0;
      }
    }
    if(nrecvs + self_msg > 0) {
      sizes_from = Teuchos::arcp(
        new int[nrecvs + self_msg], 0, nrecvs + self_msg, true);
    }

    /* Several cases:
     1. indices_to == NULL
     => starts_to != NULL, need to allocate, set starts_to_ptr
     2. indices_to != NULL (=> starts_to == NULL)
     need to allocate, set indices_to_ptr
     3,4. mirror cases for _from
     */
    if(nsends + self_msg > 0) {
      starts_to_ptr = Teuchos::arcp(
        new int[nsends + self_msg], 0, nsends + self_msg, true);
    }
    if (plan->indices_to == Teuchos::null) {
      /* Simpler case; sends already blocked by processor */
      Teuchos::ArrayRCP<int> index;
      Teuchos::ArrayRCP<int> sort_val;
      if(nsends + self_msg > 0) {
        index = Teuchos::arcp(
          new int[nsends + self_msg], 0, nsends + self_msg, true);
        sort_val = Teuchos::arcp(
          new int[nsends + self_msg], 0, nsends + self_msg, true);
      }
      for (int i = 0; i < nsends + self_msg; i++) {
        int j = plan->starts_to[i];

        for (int k = 0; k < plan->lengths_to[i]; k++) {
          sizes_to[i] += sizes[j++];
        }
        if (sizes_to[i] > plan->max_send_size &&
            plan->procs_to[i] != my_proc)
          plan->max_send_size = sizes_to[i];
      }
      for (int i = 0; i < nsends + self_msg; i++) {
        sort_val[i] = plan->starts_to[i];
        index[i] = i;
      }
      sort_ints(sort_val, index);
      int sum = 0;
      for (int i = 0; i < nsends + self_msg; i++) {
        starts_to_ptr[index[i]] = sum;
        sum += sizes_to[index[i]];
      }
    }
    else {		/* Harder case, sends not blocked */
      Teuchos::ArrayRCP<int> offset;
      if(plan->nvals > 0) {
        offset = Teuchos::arcp(new int[plan->nvals], 0, plan->nvals, true);
      }
      indices_to_ptr.resize(plan->nvals);

      /* Compute address for every item in send array */
      int sum = 0;
      for (int i = 0; i < plan->nvals; i++) {
        offset[i] = sum;
        sum += sizes[i];
      }

      sum = 0;
      plan->max_send_size = 0;
      for (int i = 0; i < nsends + self_msg; i++) {
        starts_to_ptr[i] = sum;
        int j = plan->starts_to[i];
        for (int k = 0; k < plan->lengths_to[i]; k++) {
          indices_to_ptr[j] = offset[plan->indices_to[j]];
          sizes_to[i] += sizes[plan->indices_to[j++]];
        }
        if (sizes_to[i] > plan->max_send_size &&
            plan->procs_to[i] != my_proc)
          plan->max_send_size = sizes_to[i];
        sum += sizes_to[i];
      }
    }

    /* Note: This routine only gets message sizes, not object sizes. */
    /*	Anything requiring item sizes requires more code */
    /*      Should such functionality reside here? */
    exchange_sizes(sizes_to, plan->procs_to, nsends, self_msg,
      sizes_from, plan->procs_from, nrecvs,
      &plan->total_recv_size, my_proc, tag, plan->comm);

    if(nrecvs + self_msg > 0) {
      starts_from_ptr = Teuchos::arcp(
        new int[nrecvs + self_msg], 0, nrecvs + self_msg, true);
    }

    if (plan->indices_from == Teuchos::null) {
      /* Simpler case; recvs already blocked by processor */
      Teuchos::ArrayRCP<int> index;
      Teuchos::ArrayRCP<int> sort_val;
      if(nrecvs + self_msg > 0) {
        index = Teuchos::arcp(
          new int[nrecvs + self_msg], 0, nrecvs + self_msg, true);
        sort_val = Teuchos::arcp<int>(
            new int[nrecvs + self_msg], 0, nrecvs + self_msg, true);
      }

      for (int i = 0; i < nrecvs + self_msg; i++) {
        sort_val[i] = plan->starts_from[i];
        index[i] = i;
      }
      sort_ints(sort_val, index);

      int sum = 0;
      for (int i = 0; i < nrecvs + self_msg; i++) {
        starts_from_ptr[index[i]] = sum;
        sum += sizes_from[index[i]];
      }
    }
    else {
      // OLD COMMENT left here for reference but to be deleted
           /* Harder case, recvs not blocked */
           /* Not currently supported */
           /* Can't do w/o individual item sizes */
           /* Currently checked for at top of file */

      // NEW METHOD
      // Note this is currently just a duplicate of above block which is working
      // since I've set up do_wait to just copy the entire block. However I am
      // not sure yet how to organize this implementation and suspect this may
      // not be correct anyways - though it seems to work. After do_wait copies
      // the data (out of order) the directory will handle resorting it since
      // each element has index to indicate where it goes in the array. In the
      // original do_wait implementation it seemed it was setup to place the
      // elemenets in order, even though they would be sorted later, so that
      // part I need to discuss further
      Teuchos::ArrayRCP<int> index;
      Teuchos::ArrayRCP<int> sort_val;
      if(nrecvs + self_msg > 0) {
        index = Teuchos::arcp(
          new int[nrecvs + self_msg], 0, nrecvs + self_msg, true);
        sort_val = Teuchos::arcp(
          new int[nrecvs + self_msg], 0, nrecvs + self_msg, true);
      }

      for (int i = 0; i < nrecvs + self_msg; i++) {
        sort_val[i] = plan->starts_from[i];
        index[i] = i;
      }
      sort_ints(sort_val, index);

      int sum = 0;
      for (int i = 0; i < nrecvs + self_msg; i++) {
        starts_from_ptr[index[i]] = sum;
        sum += sizes_from[index[i]];
      }
    }
  }
  plan->sizes_to = sizes_to;
  plan->sizes_from = sizes_from;
  plan->starts_to_ptr = starts_to_ptr;
  plan->starts_from_ptr = starts_from_ptr;
  plan->indices_to_ptr = indices_to_ptr;
  plan->indices_from_ptr = indices_from_ptr;

  if (sum_recv_sizes) {
    *sum_recv_sizes = plan->total_recv_size;
  }

  return 0;
}

int Zoltan2_Directory_Comm::exchange_sizes(
  const Teuchos::ArrayRCP<int> &sizes_to,		/* value I need to exchange (size of true msg) */
  const Teuchos::ArrayRCP<int> &procs_to,		/* procs I send to */
  int       nsends,		                /* number of messages I'll send */
  int       self_msg,		              /* do I copy data to myself? */
  Teuchos::ArrayRCP<int> &sizes_from,	    	/* (returned) size of all my receives */
  const Teuchos::ArrayRCP<int> &procs_from, /* procs I recv from */
  int       nrecvs,		                /* number of messages I receive */
  int      *total_recv_size,	        /* (returned) sum of all incoming sizes */
  int       my_proc,		              /* my processor number */
  int       tag,			                /* message tag I can use */
  Teuchos::RCP<const Teuchos::Comm<int> > /* comm */) {		/* communicator */

  /* If sizes vary, then I need to communicate messaaage lengths */
  int self_index_to = -1;  /* location of self in procs_to */
  for (int i = 0; i < nsends + self_msg; i++) {
    if (procs_to[i] != my_proc) {
      // TODO: Teuchos::send(...)
      MPI_Send((void *) &sizes_to[i], 1, MPI_INT, procs_to[i], tag, getRawComm());
    }
    else {
      self_index_to = i;
    }
  }

  *total_recv_size = 0;

  for (int i = 0; i < nrecvs + self_msg; i++) {
    if (procs_from[i] != my_proc) {
      MPI_Status status;		/* status of commuication operation */
      // TODO: Teuchos::receive(...)
      MPI_Recv((void *) &(sizes_from[i]), 1, MPI_INT, procs_from[i],
        tag, getRawComm(), &status);
    }
    else {
      sizes_from[i] = sizes_to[self_index_to];
    }
    *total_recv_size += sizes_from[i];
  }
  return 0;
}

} // end namespace Zoltan2
