/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __COMM_H
#define __COMM_H

#include <mpi.h>

#ifndef TRUE
#define FALSE (0)
#define TRUE  (1)
#endif /* !TRUE */

/*
 * Error codes for Comm library
 *   COMM_OK     - no errors
 *   COMM_WARN   - some warning occurred in COMM library;
 *                 application should be able to continue running
 *   COMM_FATAL  - a fatal error occurred
 *   COMM_MEMERR - memory allocation failed; it might be possible to
 *                 use a different, more memory-friendly, algorithm
 */
#define COMM_OK     0
#define COMM_WARN   11
#define COMM_FATAL  -11
#define COMM_MEMERR -12


#define COMM_ERROR(a,b,c) \
  fprintf(stderr, "[%d] COMM ERROR in %s (line %d): %s", \
  (c), (b), __LINE__, (a));

#define COMM_WARNING(a,b,c) \
  fprintf(stderr, "[%d] COMM WARNING in %s (line %d): %s", \
  (c), (b), __LINE__, (a));





/* Data structures for communication object. */

struct Comm_Obj {		/* data for mapping between decompositions */
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
    int       self_msg;		/* do I have data for myself? */
    int       max_send_size;	/* size of longest message I send (w/o self) */
    int       total_recv_size;	/* total amount of data I'll recv (w/ self) */
    MPI_Comm  comm;		/* communicator for operation */
    MPI_Request *request;       /* MPI requests for posted recvs */
    MPI_Status *status;		/* MPI status for those recvs */
};

typedef struct Comm_Obj COMM_OBJ;

/* function prototypes */

extern int LB_Comm_Create(COMM_OBJ **, int, int *, MPI_Comm, int, int *);

extern int LB_Comm_Destroy(COMM_OBJ **);

extern int LB_Comm_Invert_Map(int *, int *, int, int, int **, int **, int *,
    int, int, int, int, MPI_Comm);

extern int LB_Comm_Sort_Ints(int *, int *, int);

extern int LB_Comm_Exchange_Sizes(int *, int, int, int *, int, int *, int *,
    int *, int, int, MPI_Comm);

extern int LB_Comm_Resize(COMM_OBJ *, int *, int);

extern int LB_Comm_Do(COMM_OBJ *, int, char *, int, char *);

extern int LB_Comm_Do_Reverse(COMM_OBJ *, int, char *, int, int *, char *);

#endif
