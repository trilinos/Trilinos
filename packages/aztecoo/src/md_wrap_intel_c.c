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
 * $Name$
 *====================================================================*/
#ifndef lint
static char *cvs_wrapint_id =
  "$Id$";
#endif

/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/

#include <stdio.h>
#include <signal.h>
#include "az_aztec_defs.h"
#ifdef AZTEC_MPI
#include <mpi.h>
#else
#define MPI_Request int
#endif

#define nCUBE 1
#define INTEL 2
#define SUN   3
#define DELTA 4
#define MACHINE INTEL
#define CUBESIZ 4096 /* a power of two >= number of processors */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void get_parallel_info(int *proc, int *nprocs, int *dim)

{
  *proc   = mynode();
  *nprocs = numnodes();
  *dim    = 0;

} /* get_parallel_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_read(char *buf, int bytes, int *source, int *type, int *flag)

{
  /*
    extern int Me, Numnodes;
    */

  long int info[8];
  long     lflag;

  crecvx(*type, buf, (long) bytes, *source, -1, info);

  if ((*source == -1) || (*type == -1)) {
     *type   = info[0];
     *source = info[2];
  }

  return ((int) info[1]);

} /* md_read */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_write(char *buf, int bytes, int dest, int type, int *flag)

{
  /*
    extern int Me;
     */
  int  Me;
  long lflag;

  Me    = mynode();
  lflag = type * CUBESIZ + Me;

  /* dest=0xffff broadcasts to the current process on all of the nodes. */

  if (dest == 0xffff)
    csend(lflag, buf,(long) bytes, -1L, 0L); /* to everyone */
  else
    csend(type, buf, (long) bytes, (long) dest, 0L);

  return 0;

} /* md_write */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_wrap_iread(void *buf, int bytes, int *source, int *type,
                  MPI_Request *request)


/*******************************************************************************

  Machine dependent wrapped message-reading communication routine for the
  Intel.  This routine is a simple no-op but is used in order to provide
  compatibility with the MPI communication routine order.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.

  source:          Source processor number.

  type:            Message type

*******************************************************************************/

{

  return 0;

} /* md_wrap_iread */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_wrap_write(void *buf, int bytes, int dest, int type, int *flag)

/*******************************************************************************

  Machine dependent wrapped message-sending communication routine for the
  Intel.  This routine is exactly the same as md_write.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.

  dest:            Destination processor number.

  type:            Message type

  flag:

*******************************************************************************/

{

  int  Me;
  long lflag;

  Me    = mynode();
  lflag = type * CUBESIZ + Me;

  /* dest=0xffff broadcasts to the current process on all of the nodes. */

  if (dest == 0xffff)
    csend(lflag, buf,(long) bytes, -1L, 0L); /* to everyone */
  else
    csend(type, buf, (long) bytes, (long) dest, 0L);

  return 0;

} /* md_wrap_write */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_wrap_wait(void *buf, int bytes, int *source, int *type, int *flag,
                 MPI_Request *request)

/*******************************************************************************

  Machine dependent wrapped message-wait communication routine for the Intel.
  This routine is identical to md_read but is put here in order to be compatible
  with the order required to do MPI communication.

  Author:          Scott A. Hutchinson, SNL, 9221
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  buf:             Beginning address of data to be sent.

  bytes:           Length of message in bytes.
  dest:            Destination processor number.

  type:            Message type

  flag:

*******************************************************************************/

{

  long int info[8];

  crecvx(*type, buf, (long) bytes, *source, -1, info);

  if ((*source == -1) || (*type == -1)) {
     *type   = info[0];
     *source = info[2];
  }

  return ((int) info[1]);

} /* md_wrap_wait */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*
 * This section is used to kill the program on all nodes when an "exit" is
 * called from the program.  The "exit" call must be redefined in a header file
 * (e.g., rf_salsa.h) to mean "osf_exit" when being compiled for OSF on the
 * Intel.  Make sure -DPARA is used when compiling.
 */


int osf_exit(int ignore, char *fname, int lno)

{
  fprintf(stderr, "global exit called from %s line number %d\n",
          fname, lno);
  fprintf(stderr, "Killing all nodes due to error.\n");
  kill(0, SIGKILL);
}

int md_wrap_iwrite(void *buf, int bytes, int dest, int type, int *flag,
                  int *request)
{
int ret_info;

ret_info = md_wrap_write(buf, bytes, dest, type, flag);
return(ret_info);

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_wrap_request_free(MPI_Request *request)

/*******************************************************************************

  Machine dependent wrapped request object deletion routine. 
  (Trivial function except for MPI version).

  Author:          Michael A. Heroux, SNL, 9214
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  request:           Pointer to an existing request object that will be freed.

*******************************************************************************/
{

  int err = 0;
  return err;

} /* md_wrap_request_free */

