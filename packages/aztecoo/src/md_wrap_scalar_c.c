/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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

#include "az_aztec_defs.h"
#ifdef AZTEC_MPI
#include <mpi.h>
#else
#define MPI_Request int
#endif

#include <stdlib.h>
#include <stdio.h>

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
extern void get_parallel_info(int *proc, int *nprocs, int *dim);
extern int md_read(char *buf, int bytes, int *source, int *type, int *flag);
extern int md_write(char *buf, int bytes, int dest, int type, int *flag);
extern int md_wrap_iread(void *buf, int bytes, int *source, int *type,
                  MPI_Request *request);
extern int md_wrap_write(void *buf, int bytes, int dest, int type, int *flag);
extern int md_wrap_wait(void *buf, int bytes, int *source, int *type, int *flag,
                 MPI_Request *request);
extern int md_wrap_iwrite(void *buf, int bytes, int dest, int type, int *flag,
                  int *request);

void get_parallel_info(int *proc, int *nprocs, int *dim)

{
  *proc   = 0;
  *nprocs = 1;
  *dim    = 0;

} /* get_parallel_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*ARGSUSED*/

int md_read(char *buf, int bytes, int *source, int *type, int *flag)

{
  return bytes;

} /* md_read */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*ARGSUSED*/

int md_write(char *buf, int bytes, int dest, int type, int *flag)

{
  return 0;

} /* md_write */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*ARGSUSED*/

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
/*ARGSUSED*/

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

  return 0;

} /* md_wrap_write */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*ARGSUSED*/

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

  return bytes;

} /* md_wrap_wait */

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

