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
#include <Teuchos_ArrayRCP.hpp>

namespace Zoltan2 {

class Zoltan2_Directory_Plan {	/* data for mapping between decompositions */
  public:
    Zoltan2_Directory_Plan() : using_sizes(false), maxed_recvs(0) {
    }

    ~Zoltan2_Directory_Plan() {
    }

    void getInvertedValues(Zoltan2_Directory_Plan * from);

    void print(const std::string& headerMessage) const;

    Teuchos::ArrayRCP<int> procs_to;      /* processors I'll send to */
    Teuchos::ArrayRCP<int> procs_from;    /* processors I'll receive from*/
    Teuchos::ArrayRCP<int> lengths_to;    /* # items I send in my messages */
    Teuchos::ArrayRCP<int> lengths_from;  /* # items I recv in my messages */

    /* Following arrays used if send/recv data is packed contiguously */
    Teuchos::ArrayRCP<int> starts_to;	    /* where in item lists each send starts */
    Teuchos::ArrayRCP<int> starts_from;	  /* where in item lists each recv starts */

    /* Following arrays used is send/recv data not packed contiguously */
    Teuchos::ArrayRCP<int> indices_to;    /* indices of items I send in my msgs */

    /* ordered consistent with lengths_to */
    Teuchos::ArrayRCP<int> indices_from;  /* indices for where to put arriving data */

				/* ordered consistent with lengths_from */

    /* Above information is sufficient if items are all of the same size */
    /* If item sizes are variable, then need following additional arrays */
    Teuchos::ArrayRCP<int> sizes;      /* size of each item to send (if items vary) */
    bool using_sizes; /* may refactor this so it's out - tracks whether we are in size mode even if 0 size */

    Teuchos::ArrayRCP<int> sizes_to;    /* size of each msg to send (if items vary) */
    Teuchos::ArrayRCP<int> sizes_from;  /* size of each msg to recv (if items vary) */

    /* Following used if send/recv data is packed contiguously & items vary */
    Teuchos::ArrayRCP<int> starts_to_ptr;	  /* where in dense array sends starts */
    Teuchos::ArrayRCP<int> starts_from_ptr;	/* where in dense each recv starts */

    /* Following used is send/recv data not packed contiguously & items vary */
    Teuchos::ArrayRCP<int> indices_to_ptr; /* where to find items I send in my msgs */
				/* ordered consistent with lengths_to */
    Teuchos::ArrayRCP<int> indices_from_ptr; /* where to find items I recv */
				/* ordered consistent with lengths_from */

    /* Note: ALL above arrays include data for self-msg */

    int       nvals;		       /* number of values I own to start */
    int       nvals_recv;	     /* number of values I own after remapping */
    int       nrecvs;		       /* number of msgs I'll recv (w/o self_msg) */
    int       nsends;		       /* number of msgs I'll send (w/o self_msg) */
    int       self_msg;		     /* do I have data for myself? */
    int       max_send_size;	 /* size of longest message I send (w/o self) */
    int       total_recv_size; /* total amount of data I'll recv (w/ self) */
    int       maxed_recvs;     /* use MPI_Alltoallv if too many receives */
    Teuchos::RCP<const Teuchos::Comm<int> > comm; /* communicator */

    // making this ArrayRCP is causing issues with name() demangling for gcc
    // sems build, but not clang... will need to work on this further
    // back to std::vector for the moment
    // there is probably an option in the sems building turning on
    // something that creates this conflict
    std::vector<MPI_Request> request;      /* MPI requests for posted recvs */
    std::vector<MPI_Status> status;		     /* MPI status for those recvs */

    Zoltan2_Directory_Plan* plan_reverse;    /* to support POST & WAIT */

    Teuchos::ArrayRCP<char> recv_buff;       /* To support POST & WAIT */
    Teuchos::ArrayRCP<char> getRecvBuff() const { return recv_buff; }
};

class Zoltan2_Directory_Comm {
  public:
    Zoltan2_Directory_Comm(
      int       nvals,		                /* number of values I currently own */
      const Teuchos::ArrayRCP<int>  &assign, /* processor assignment for all values */
      Teuchos::RCP<const Teuchos::Comm<int> > comm,	          /* communicator */
      int       tag);			                           /* message tag I can use */

    ~Zoltan2_Directory_Comm();

    int do_forward(
      int tag,		      	                    /* message tag for communicating */
      const Teuchos::ArrayRCP<char> &send_data,	    /* array of data I currently own */
      int nbytes,                             /* msg size */
      Teuchos::ArrayRCP<char> &recv_data);     /* array of data to receive */

    int do_reverse(
      int tag,			                         /* message tag for communicating */
      const Teuchos::ArrayRCP<char> &send_data,    /* array of data I currently own */
      int nbytes,                            /* msg size */
      const Teuchos::ArrayRCP<int> &sizes,
      Teuchos::ArrayRCP<char> &recv_data);     /* array of data owned after reverse */

    int getNRec() const { return nrec; } /* accessor for nrec */

    int get_plan_forward_recv_size() const {
      return plan_forward->total_recv_size;
    }

    int resize(const Teuchos::ArrayRCP<int> &sizes, int tag,
      int *sum_recv_sizes);

  private:
    int resize(Zoltan2_Directory_Plan *plan,
      const Teuchos::ArrayRCP<int> &sizes, int tag, int *sum_recv_sizes);

    int do_post(Zoltan2_Directory_Plan *plan, int tag,
      const Teuchos::ArrayRCP<char> &send_data,
      int nbytes,                            /* msg size */
      Teuchos::ArrayRCP<char> &recv_data);

    int do_wait(Zoltan2_Directory_Plan *plan, int tag,
      const Teuchos::ArrayRCP<char> &send_data,
      int nbytes,                            /* msg size */
      Teuchos::ArrayRCP<char> &recv_data);

    int do_all_to_all(Zoltan2_Directory_Plan *plan,
      const Teuchos::ArrayRCP<char> &send_data,
      int nbytes,                            /* msg size */
      Teuchos::ArrayRCP<char> &recv_data);

    int sort_ints(Teuchos::ArrayRCP<int> &vals_sort, Teuchos::ArrayRCP<int> &vals_other);

    int invert_map(const Teuchos::ArrayRCP<int> &lengths_to,
      const Teuchos::ArrayRCP<int> &procs_to, int nsends, int self_msg,
      Teuchos::ArrayRCP<int> &lengths_from, Teuchos::ArrayRCP<int> &procs_from,
      int *pnrecvs,	int my_proc,int nprocs,	int out_of_mem, int tag,
      Teuchos::RCP<const Teuchos::Comm<int> > comm);

    int exchange_sizes(const Teuchos::ArrayRCP<int> &sizes_to,
      const Teuchos::ArrayRCP<int> &procs_to, int nsends,
      int self_msg,	Teuchos::ArrayRCP<int> &sizes_from,
      const Teuchos::ArrayRCP<int> &procs_from,
      int nrecvs, int *total_recv_size,	int my_proc, int tag,
      Teuchos::RCP<const Teuchos::Comm<int> > comm);

    void free_reverse_plan(Zoltan2_Directory_Plan *plan);

    int create_reverse_plan(int tag, const Teuchos::ArrayRCP<int> &sizes);

    MPI_Comm getRawComm() {
    #ifdef HAVE_MPI
      return Teuchos::getRawMpiComm(*comm_);
    #else
      return MPI_COMM_WORLD;
    #endif
    }

    Teuchos::RCP<const Teuchos::Comm<int> > comm_;
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

} // end namespace Zoltan2

#endif
