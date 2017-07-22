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

#ifndef ZOLTAN2_DIRECTORY_COMM_H_
#define ZOLTAN2_DIRECTORY_COMM_H_

#include <Teuchos_CommHelpers.hpp>
#include <vector>
#include <mpi.h>

#ifdef HAVE_MPI

namespace Zoltan2 {

class Zoltan2_Directory_Plan {	/* data for mapping between decompositions */
  public:
    Zoltan2_Directory_Plan() :
      indices_to_allocated(false), indices_from_allocated(false),
      sizes_allocated(false),
      maxed_recvs(0), recv_buff(NULL), bOwnRecvBuff(false) {
    }
    ~Zoltan2_Directory_Plan() {
      if(bOwnRecvBuff) {
        delete recv_buff;
      }
    }

    void getInvertedValues(Zoltan2_Directory_Plan * from);

    void print(const std::string& headerMessage) const;

    std::vector<int> procs_to;      /* processors I'll send to */
    std::vector<int> procs_from;    /* processors I'll receive from*/
    std::vector<int> lengths_to;    /* # items I send in my messages */
    std::vector<int> lengths_from;  /* # items I recv in my messages */

    /* Following arrays used if send/recv data is packed contiguously */
    std::vector<int> starts_to;	    /* where in item lists each send starts */
    std::vector<int> starts_from;	  /* where in item lists each recv starts */

    /* Following arrays used is send/recv data not packed contiguously */
    std::vector<int> indices_to;    /* indices of items I send in my msgs */
    bool indices_to_allocated; // will need to rethink - new vector (loses NULL vs 0)

				/* ordered consistent with lengths_to */
    std::vector<int> indices_from;  /* indices for where to put arriving data */
    bool indices_from_allocated; // will need to rethink - new vector (loses NULL vs 0)

				/* ordered consistent with lengths_from */

    /* Above information is sufficient if items are all of the same size */
    /* If item sizes are variable, then need following additional arrays */
    std::vector<int> sizes;      /* size of each item to send (if items vary) */
    bool sizes_allocated; // will need to rethink - new vector (loses NULL vs 0)
				/* Note only on sending processor: */
				/* assuming recv proc can figure it out */

    std::vector<int> sizes_to;    /* size of each msg to send (if items vary) */
    std::vector<int> sizes_from;  /* size of each msg to recv (if items vary) */

    /* Following used if send/recv data is packed contiguously & items vary */
    std::vector<int> starts_to_ptr;	  /* where in dense array sends starts */
    std::vector<int> starts_from_ptr;	/* where in dense each recv starts */

    /* Following used is send/recv data not packed contiguously & items vary */
    std::vector<int> indices_to_ptr; /* where to find items I send in my msgs */
				/* ordered consistent with lengths_to */
    std::vector<int> indices_from_ptr; /* where to find items I recv */
				/* ordered consistent with lengths_from */

    /* Note: ALL above arrays include data for self-msg */

    int       nvals;		       /* number of values I own to start */
    int       nvals_recv;	     /* number of values I own after remapping */
    int       nrecvs;		       /* number of msgs I'll recv (w/o self_msg) */
    int       nsends;		       /* number of msgs I'll send (w/o self_msg) */
    int       nindices_to;
    int       nindices_from;
    int       self_msg;		     /* do I have data for myself? */
    int       max_send_size;	 /* size of longest message I send (w/o self) */
    int       total_recv_size; /* total amount of data I'll recv (w/ self) */
    int       maxed_recvs;     /* use MPI_Alltoallv if too many receives */
    Teuchos::RCP<const Teuchos::Comm<int> > comm; /* communicator */
    std::vector<MPI_Request> request;      /* MPI requests for posted recvs */
    std::vector<MPI_Status> status;		     /* MPI status for those recvs */

    Zoltan2_Directory_Plan* plan_reverse;    /* to support POST & WAIT */

    std::vector<char> *recv_buff;            /* To support POST & WAIT */
    char * getRecvBuff() const { return &((*recv_buff)[0]); }
    size_t getRecvBuffSize() const { return recv_buff->size(); }
    bool bOwnRecvBuff;
};

class Zoltan2_Directory_Comm {
  public:
    Zoltan2_Directory_Comm(
      int       nvals,		                /* number of values I currently own */
      const std::vector<int>  &assign, /* processor assignment for all values */
      Teuchos::RCP<const Teuchos::Comm<int> > comm,	          /* communicator */
      int       tag);			                           /* message tag I can use */

    ~Zoltan2_Directory_Comm();

    int do_forward(
      int tag,		      	                    /* message tag for communicating */
      const std::vector<char> &send_data,	    /* array of data I currently own */
      int nbytes,                             /* msg size */
      std::vector<char> &recv_data);          /* array of data to receive */

    int do_reverse(
      int tag,			                         /* message tag for communicating */
      const std::vector<char> &send_data,    /* array of data I currently own */
      int nbytes,                            /* msg size */
      const std::vector<int> &sizes,
      std::vector<char> &recv_data);     /* array of data owned after reverse */

    int getNRec() const { return nrec; } /* accessor for nrec */

    int get_plan_forward_recv_size() const {
      return plan_forward->total_recv_size;
    }

    int resize(const std::vector<int> &sizes, int tag,
      int *sum_recv_sizes);

  private:
    int resize(Zoltan2_Directory_Plan *plan,
      const std::vector<int> &sizes, int tag, int *sum_recv_sizes);

    int do_post(Zoltan2_Directory_Plan *plan, int tag,
      const std::vector<char> &send_data,
      int nbytes,                            /* msg size */
      std::vector<char> &recv_data);

    int do_wait(Zoltan2_Directory_Plan *plan, int tag,
      const std::vector<char> &send_data,
      int nbytes,                            /* msg size */
      std::vector<char> &recv_data);

    int do_all_to_all(Zoltan2_Directory_Plan *plan,
      const std::vector<char> &send_data,
      int nbytes,                            /* msg size */
      std::vector<char> &recv_data);

    int sort_ints(std::vector<int> &vals_sort, std::vector<int> &vals_other);

    int invert_map(const std::vector<int> &lengths_to,
      const std::vector<int> &procs_to, int nsends, int self_msg,
      std::vector<int> &lengths_from, std::vector<int> &procs_from,
      int *pnrecvs,	int my_proc,int nprocs,	int out_of_mem, int tag,
      Teuchos::RCP<const Teuchos::Comm<int> > comm);

    int exchange_sizes(const std::vector<int> &sizes_to,
      const std::vector<int> &procs_to, int nsends,
      int self_msg,	std::vector<int> *sizes_from,
      const std::vector<int> &procs_from,
      int nrecvs, int *total_recv_size,	int my_proc, int tag,
      Teuchos::RCP<const Teuchos::Comm<int> > comm);

    void free_reverse_plan(Zoltan2_Directory_Plan *plan);

    int create_reverse_plan(int tag, const std::vector<int> &sizes);

    Zoltan2_Directory_Plan * plan_forward; // for efficient MPI communication
    int nrec;
};

// -----------------------------------------------------------------------------
// TODO: Decide how to handle this code - copied from zoltan - some may be relic
    /* Red Storm MPI permits a maximum of 2048 receives.  We set our
     * limit of posted receives to 2000, leaving some for the application.
     */
    #ifndef MPI_RECV_LIMIT
    /* Decided for Trilinos v10/Zoltan v3.2 would almost always use */
    /* MPI_Alltoall communication instead of point-to-point.        */
    /* August 2009 */
    /* #define MPI_RECV_LIMIT 4 */

    /* Decided for zoltan_gid_64 branch to always used posted receives because
     * Alltoall requires that offsets be 32-bit integers.  October 2010
     */
    #define MPI_RECV_LIMIT 0
    /* #define MPI_RECV_LIMIT 2000 */
    #endif
// -----------------------------------------------------------------------------

}; // end namespace Zoltan2

#endif // #ifdef HAVE_MPI

#endif