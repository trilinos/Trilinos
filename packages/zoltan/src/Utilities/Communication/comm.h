// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef __COMM_H
#define __COMM_H

#include <mpi.h>
#include "zoltan_types.h"
#include "zoltan_util.h"
#include "zoltan_comm.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* Data structures and macros for the Zoltan Communication Package.  */
/* This file should be included only by communication package files. */
/* Communication package users should include zoltan_comm.h.          */

#ifndef TRUE
#define FALSE (0)
#define TRUE  (1)
#endif /* !TRUE */

#define ZOLTAN_COMM_ERROR(a,b,c) \
  ZOLTAN_PRINT_ERROR((c),(b),(a));

#define ZOLTAN_COMM_WARNING(a,b,c) \
  ZOLTAN_PRINT_WARNING((c),(b),(a));





/* Data structures for communication object. */

struct Zoltan_Comm_Obj {	/* data for mapping between decompositions */
    int      *procs_to;         /* processors I'll send to */
    int      *procs_from;       /* processors I'll receive from*/
    int      *lengths_to;       /* # items I send in my messages */
    int      *lengths_from;     /* # items I recv in my messages */

    /* Following arrays used if send/recv data is packed contiguously */
    int      *starts_to;	/* where in item lists each send starts */
    int      *starts_from;	/* where in item lists each recv starts */

    /* Following arrays used is send/recv data not packed contiguously */
    int      *indices_to;       /* indices of items I send in my msgs */
				/* ordered consistent with lengths_to */
    int      *indices_from;     /* indices for where to put arriving data */
				/* ordered consistent with lengths_from */

    /* Above information is sufficient if items are all of the same size */
    /* If item sizes are variable, then need following additional arrays */
    int      *sizes;            /* size of each item to send (if items vary) */
				/* Note only on sending processor: */
				/* assuming recv proc can figure it out */

    int      *sizes_to;         /* size of each msg to send (if items vary) */
    int      *sizes_from;       /* size of each msg to recv (if items vary) */

    /* Following used if send/recv data is packed contiguously & items vary */
    int      *starts_to_ptr;	/* where in dense array sends starts */
    int      *starts_from_ptr;	/* where in dense each recv starts */

    /* Following used is send/recv data not packed contiguously & items vary */
    int      *indices_to_ptr;   /* where to find items I send in my msgs */
				/* ordered consistent with lengths_to */
    int      *indices_from_ptr; /* where to find items I recv */
				/* ordered consistent with lengths_from */

    /* Note: ALL above arrays include data for self-msg */

    int       nvals;		/* number of values I own to start */
    int       nvals_recv;	/* number of values I own after remapping */
    int       nrecvs;		/* number of msgs I'll recv (w/o self_msg) */
    int       nsends;		/* number of msgs I'll send (w/o self_msg) */
    int       nindices_to;
    int       nindices_from;
    int       self_msg;		/* do I have data for myself? */
    int       max_send_size;	/* size of longest message I send (w/o self) */
    int       total_recv_size;	/* total amount of data I'll recv (w/ self) */
    int       maxed_recvs;      /* do I have to many receives to post all
                                 * at once? if so use MPI_Alltoallv */
    MPI_Comm  comm;		/* communicator for operation */
    MPI_Request *request;       /* MPI requests for posted recvs */
    MPI_Status *status;		/* MPI status for those recvs */
    
    ZOLTAN_COMM_OBJ* plan_reverse;   /* to support POST & WAIT */
    char*     recv_buff;  /* To support POST & WAIT */    
};

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


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
