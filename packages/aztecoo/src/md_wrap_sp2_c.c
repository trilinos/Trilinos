/*@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

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

#include <stdio.h>
#if defined (MPL)
#include <sys/utsname.h>
#include <mpproto.h>
#elif defined (MPI)
#include <sys/utsname.h>
#include <mpi.h>
#endif

void get_parallel_info(int *proc, int *nprocs, int *dim)

{

#if defined (MPL) || defined (MPI)
  int i;
  struct utsname sp2_host;
#endif

#if defined (MPL)
  mpc_environ(nprocs, proc);
#elif defined (MPI)
  MPI_Comm_size(MPI_COMM_WORLD, nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, proc);
#endif

#if defined (MPL) || defined (MPI)
  uname(&sp2_host);
  (void) fprintf(stderr, "Task %d is on node %s\n", *proc, sp2_host.nodename);
#else
  *proc   = 0;
  *nprocs = 1;
#endif
  *dim = 0;

} /* get_parallel_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void md_whoami(int *node, int *proc, int *host, int *dim)

{

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_read(char *buf, int bytes, int *source, int *type, int *flag)

{

#if defined (MPL)
  int err, nbytes, buffer;

  if (bytes == 0) {
    err = mpc_brecv(&buffer, 1, source, type, &nbytes);
  } else {
    err = mpc_brecv(buf, bytes, source, type, &nbytes);
  }

  if (err != 0) (void) fprintf(stderr, "mpc_brecv error = %d\n", mperrno);

#elif defined (MPI)
  int err, buffer;
  MPI_Status status;

  if (bytes == 0) {
    err = MPI_Recv(&buffer, 1, MPI_BYTE, *source, *type, MPI_COMM_WORLD,
                   &status);
  } else {
    err = MPI_Recv(buf, bytes, MPI_BYTE, *source, *type, MPI_COMM_WORLD,
                   &status);
  }
  if (err!=0) (void) fprintf(stderr, "MPI_Recv error = %d\n", err);
#endif

  return bytes;

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int md_write(char *buf, int bytes, int dest, int type, int *flag)

{
#if defined (MPL)
  int err, buffer;

  if (bytes == 0) {
      err = mpc_bsend(&buffer, 1, dest, type);
  } else {
      err = mpc_bsend(buf, bytes, dest, type);
  }

  if (err!=0) (void) fprintf(stderr, "mpc_bsend error = %d\n", mperrno);

#elif defined (MPI)
  int err, buffer;

  if (bytes == 0) {
      err = MPI_Send(&buffer, 1, MPI_BYTE, dest, type, MPI_COMM_WORLD);
  } else {
      err = MPI_Send(buf, bytes, MPI_BYTE, dest, type, MPI_COMM_WORLD);
  }

  if (err != 0) (void) fprintf(stderr, "MPI_Send error = %d\n", err);
#endif

  return 0;

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

