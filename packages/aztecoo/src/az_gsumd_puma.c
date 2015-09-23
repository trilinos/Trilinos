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

#include <stdlib.h>     /* these appear outside of ifdef */
#include <stdio.h>      /* to avoid a lint warning       */

#ifdef PUMA_GSUMD       /* Only compiled in on request. */

#include <sunmos.h>

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

double AZ_gsum_double(double val, int proc_config[])

/*******************************************************************************

  Sum a double value across all processors. Use a low level call and pack value
  into header for efficiency. NOTE: must be compiled WITHOUT  -Mnodepchk  and
  -Msafeptr.

  Author:          Bruce Hendrickson, SNL, 1422
  =======

  Return code:     double, result of global sum.
  ============

  Parameter list:
  ===============

  val:             Value to be summed from each processor.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int     type;                 /* type of next message */
  int     data;                 /* integer message being sent */
  int     partner;              /* processor I exchange with */
  int     mask;                 /* bit pattern identifying partner */
  int     hbit;                 /* largest nonzero bit in nprocs */
  int     nprocs_small;         /* largest power of 2 <= nprocs */
  int    *vptr;                 /* pointer into parts of val */
  int    *vptrp;                /* pointer into parts of val */
  int    *v2ptr;                /* pointer into parts of val2 */
  int    *v2ptrp;               /* pointer into parts of val2 */
  double  val2;                 /* arriving value to add */
  int     cflag;                /* dummy argument for compatability */
  void    exit();
  int     _nsend(), _nrecv();
  int     node, nprocs;

  /**************************** execution begins ******************************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  vptr   = (int *) &val;
  vptrp  = vptr + 1;
  v2ptr  = (int *) &val2;
  v2ptrp = v2ptr + 1;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    data = *vptr;
    type = *vptrp;
    if (_nsend((char *) NULL, 0, partner, type, NULL, data)) {
      (void)fprintf(stderr, "AZ_gsum_double: ERROR on node %d\n", node);
      (void)fprintf(stderr, "nwrite failed, message type = %d\n", type);
      exit(-1);
    }
  }
  else if (node + nprocs_small < nprocs) {
    type = -1;
    cflag = 0;
    if (_nrecv((char *) NULL, &cflag, &partner, &type, NULL, &data) != 0) {
      (void)fprintf(stderr, "AZ_gsum_double: ERROR on node %d\n", node);
      (void)fprintf(stderr, "nwrite failed, message type = %d\n", type);
      exit(-1);
    }

    *(v2ptr)  = data;
    *(v2ptrp) = type;
    val      += *((double *) v2ptr);
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;
      data    = *vptr;
      type    = *(vptrp);
      if (_nsend((char *) NULL, 0, partner, type, NULL, data)) {
        (void)fprintf(stderr, "AZ_gsum_double: ERROR on node %d\n", node);
        (void)fprintf(stderr, "nwrite failed, message type = %d\n", type);
        exit(-1);
      }
      type  = -1;
      cflag = 0;
      if (_nrecv((char *)NULL, &cflag, &partner, &type, NULL, &data) != 0) {
        (void) fprintf(stderr, "AZ_gsum_double: ERROR on node %d\n",
                       node);
        (void) fprintf(stderr, "nwrite failed, message type = %d\n",
                       type);
        exit(-1);
      }

      *(v2ptr)  = data;
      *(v2ptrp) = type;
      val      += *((double *) v2ptr);
    }
  }
  else type += hbit;

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    type  = -1;
    cflag = 0;
    if (_nrecv((char *) NULL, &cflag, &partner, &type, NULL, &data) != 0) {
      (void) fprintf(stderr, "AZ_gsum_double: ERROR on node %d\n", node);
      (void) fprintf(stderr, "nwrite failed, message type = %d\n", type);
      exit(-1);
    }

    *(vptr)  = data;
    *(vptrp) = type;
  }
  else if (node+nprocs_small < nprocs) {
    data = *vptr;
    type = *(vptrp);
    if (_nsend((char *) NULL, 0, partner, type, NULL, data)) {
      (void) fprintf(stderr, "AZ_gsum_double: ERROR on node %d\n", node);
      (void) fprintf(stderr, "nwrite failed, message type = %d\n", type);
      exit(-1);
    }
  }

  return val;
}

#endif /* ifdef PUMA_GSUMD */
