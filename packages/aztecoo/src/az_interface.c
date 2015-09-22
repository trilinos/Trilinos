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

#include <stdio.h>
#include <stdlib.h>
#include "az_aztec.h"

extern double machine_dependent_second(void);
extern void get_parallel_info(int *proc, int *nprocs, int *dim);
#ifdef AZTEC_MPI
extern void parallel_info(int *proc, int *nprocs, int *dim, MPI_Comm comm);
#else
extern void parallel_info(int *proc, int *nprocs, int *dim);
#endif

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_second(void)

{

  return machine_dependent_second();

} /* AZ_second */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_processor_info(int proc_config[])

{
  static int first_time = 1;

  AZ__MPI_comm_space_ok();
  get_parallel_info(&(proc_config[AZ_node]), &(proc_config[AZ_N_procs]),
                    &(proc_config[AZ_dim]));
  proc_config[AZ_Comm_MPI] = 1;
  if ((proc_config[AZ_node] == 0) && (first_time==1))
     printf("Warning: AZ_processor_info() command should be replaced with AZ_set_proc_config()\n");
  first_time = 0;

} /* AZ_processor_info */


void AZ_set_proc_config(int proc_config[], MPI_AZComm comm)

{


  AZ__MPI_comm_space_ok();
#ifdef AZTEC_MPI
  parallel_info(&(proc_config[AZ_node]), &(proc_config[AZ_N_procs]),
                &(proc_config[AZ_dim]), comm);
#else
  get_parallel_info(&(proc_config[AZ_node]), &(proc_config[AZ_N_procs]),
                    &(proc_config[AZ_dim]));
  proc_config[AZ_Comm_MPI] = 1;
#endif
  AZ_set_comm(proc_config, comm);

} /* get_parallel_info */
MPI_AZComm *AZ_get_comm(int proc_config[])
{
   if (proc_config[AZ_Comm_Set] != AZ_Done_by_User) {
      printf("Error(AZ_get_comm):Communicator not set. Use AZ_set_comm()\n   ");
      printf("              (e.g. AZ_set_comm(proc_config,MPI_COMM_WORLD)).\n");
      exit(1);
   }
#ifndef AZTEC_MPI
   printf("Warning(AZ_get_comm):Not using MPI? Returning bogus communicator\n");
#endif
   return( (MPI_AZComm *) &(proc_config[AZ_Comm_MPI]) );
}

void AZ_set_comm(int proc_config[], MPI_AZComm comm)
{
#ifdef AZTEC_MPI
  char *ptr1, *ptr2;
  int  i;
#endif

  AZ__MPI_comm_space_ok();
#ifdef AZTEC_MPI
  ptr1 = (char *) &comm;
  ptr2 = (char *) &(proc_config[AZ_Comm_MPI]);
  for (i = 0; i < sizeof(MPI_AZComm); i++) ptr2[i] = ptr1[i];

  proc_config[AZ_MPI_Tag] = AZ_MSG_TYPE;
#else
  proc_config[AZ_Comm_MPI] = 1;
#endif
  proc_config[AZ_Comm_Set] = AZ_Done_by_User;
}
