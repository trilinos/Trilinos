/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * Zoltan is distributed under the GNU Lesser General Public License 2.1.    * 
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

/* Data structures for communication object. */

struct Comm_Obj {		/* data for mapping between decompositions */
    int      *lengths_to;       /* lengths of messages I'll send */
    int      *procs_to;         /* processors I'll send to (-1 ends) */
    int      *indices_to;       /* offsets for where to read sending data */
    int      *lengths_from;     /* lengths of messages I'll receive */
    int      *procs_from;       /* processors I'll receive from (-1 ends) */
    int      *indices_from;     /* offsets for where to put arriving data */
    int       nrecvs;		/* number of msgs I'll recv (w/o self_msg) */
    int       nsends;		/* number of msgs I'll send (w/o self_msg) */
    int       self_msg;		/* do I have data for myself? */
    int       max_send_length;  /* size of longest message I send */
    int       total_recv_length;/* total amount of data I'll recv */
    MPI_Comm  comm;		/* communicator for operation */
    MPI_Request *request;       /* MPI requests for posted recvs */
    MPI_Status *status;		/* MPI status for those recvs */
};

typedef struct Comm_Obj COMM_OBJ;

/* function prototypes */

extern int LB_Comm_Create(COMM_OBJ **, int, int *, MPI_Comm, int, int, int *);
extern int LB_Comm_Do(COMM_OBJ *, int, char *, int, char *);
extern int LB_Comm_Do_Reverse(COMM_OBJ *, int, char *, int, char *);
extern int LB_Comm_Destroy(COMM_OBJ **);
extern int LB_Invert_Map(int *, int *, int, int, int **, int **, int *,
    int, int, int, int, MPI_Comm);

#endif
