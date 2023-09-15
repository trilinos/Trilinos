/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
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

/* Program to test MPI_Comm_split, etc., on thunderbird */

#include <iostream>
#include <mpi.h>

#include "get_heap_usage.h"

#define NUM_ITER 10

/////////////////////////////////////////////////////////////////////////////
void test_function()
{
MPI_Comm local_comm, tmp_comm;
int myproc, nprocs;               // MPI info wrt comm form zoltan_get_global_comm().
int local_myproc, local_nprocs;   // MPI info wrt local_comm.
int set, procmid;
int commcnt = 0;
size_t oldheap, newheap;
static int itercnt = 0;

  MPI_Comm_size(zoltan_get_global_comm(), &nprocs);
  MPI_Comm_rank(zoltan_get_global_comm(), &myproc);

  //  Duplicate zoltan_get_global_comm() to local communicator.
  oldheap = get_heap_usage();
  std::cout << "KDD " << myproc 
            << " ITER " << itercnt
            << " BEFORE Comm_dup:  " << oldheap << std::endl;
  MPI_Comm_dup(zoltan_get_global_comm(),&local_comm);
  newheap = get_heap_usage();
  std::cout << "KDD " << myproc 
            << " ITER " << itercnt
            << " AFTER  Comm_dup:  " << newheap 
            << " Used: " << newheap - oldheap << std::endl;
  commcnt++;

  //  Set up loop for split.
  local_nprocs = nprocs;
  local_myproc = myproc;
  while (local_nprocs > 1) { 

    std::cout << "KDD " << myproc << "(" << local_myproc << ")" 
            << " In main loop: local_nprocs = " << local_nprocs << std::endl;

    //  Split communicator in half.
    procmid = local_nprocs / 2;
    if (local_myproc < procmid) set = 0; /* set = LOWERHALF */
    else set = 1;                        /* set = UPPERHALF; */

    oldheap = get_heap_usage();
    std::cout << "KDD " << myproc << "(" << local_myproc << ")" 
              << " BEFORE Comm_split: " << oldheap << std::endl;
    MPI_Comm_split(local_comm,set,local_myproc,&tmp_comm);
    commcnt++;

    //  Free old local_comm; keep new one.
    newheap = get_heap_usage();
    std::cout << "KDD " << myproc << "(" << local_myproc << ")" 
              << " AFTER  Comm_split: " << newheap 
              << " Used: " << newheap - oldheap << std::endl;
    oldheap = get_heap_usage();
    std::cout << "KDD " << myproc << "(" << local_myproc << ")" 
              << " BEFORE local Comm_free: " << oldheap << std::endl;
    MPI_Comm_free(&local_comm);
    commcnt--;
    newheap = get_heap_usage();
    std::cout << "KDD " << myproc << "(" << local_myproc << ")" 
              << " AFTER  local Comm_free: " << newheap
              << " Freed: " << oldheap - newheap << std::endl;
    local_comm = tmp_comm;

    //  Update MPI info wrt new local_comm.
    MPI_Comm_rank(local_comm, &local_myproc);
    MPI_Comm_size(local_comm, &local_nprocs);
  }

  // Free local_comm.
  oldheap = get_heap_usage();
  std::cout << "KDD " << myproc 
            << " ITER " << itercnt
            << " BEFORE final Comm_free:  " << oldheap
            << std::endl;
  MPI_Comm_free(&local_comm);
  commcnt--;
  newheap = get_heap_usage();
  std::cout << "KDD " << myproc 
            << " ITER " << itercnt
            << " AFTER  final Comm_free:  " << newheap
            << " Freed: " << oldheap - newheap 
            << " commcnt = " << commcnt << std::endl;

  itercnt++;
}

/////////////////////////////////////////////////////////////////////////////
main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  size_t initheap = get_heap_usage();
  for (int i = 0; i < NUM_ITER; i++) test_function();
  size_t finalheap = get_heap_usage();

  int myproc;
  MPI_Comm_rank(zoltan_get_global_comm(), &myproc);
  std::cout << "KDDEND " << myproc 
            << " Total leaked " << finalheap - initheap
            << " Avg per iteration " << (finalheap - initheap) / NUM_ITER
            << std::endl;

  MPI_Finalize();

  return(0);  
}

