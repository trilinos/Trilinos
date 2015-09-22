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

/////////////////////////////////////////////////////////////////
// Program to test MPI_Comm_split memory usage on thunderbird. //
/////////////////////////////////////////////////////////////////

#include <iostream>
#include <mpi.h>

#include "get_heap_usage.h"

#define NUM_ITER 1000

size_t total_leak = 0;

/////////////////////////////////////////////////////////////////////////////
void test_function()
{
MPI_Comm local_comm;
int myproc, nprocs;               // MPI info wrt MPI_COMM_WORLD.
int set;
size_t oldheap, newheap;
size_t used, freed;
static int itercnt = 0;
int ierr;

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  //  Split MPI_COMM_WORLD to half-sized local communicator.
  if (myproc < nprocs/2) set = 0;      // set = LOWERHALF;
  else set = 1;                        // set = UPPERHALF;

  oldheap = get_heap_usage();
  std::cout << "KDD " << myproc 
            << " ITER " << itercnt
            << " BEFORE Comm_split:  " << oldheap << std::endl;
  ierr = MPI_Comm_split(MPI_COMM_WORLD, set, myproc, &local_comm);
  newheap = get_heap_usage();
  used = newheap - oldheap;
  if (ierr != MPI_SUCCESS) std::cout << " ERROR SPLIT " << ierr << std::endl;
  std::cout << "KDD " << myproc 
            << " ITER " << itercnt
            << " AFTER  Comm_split:  " << newheap 
            << " Used: " << used << std::endl;

  // Free local_comm.
  oldheap = get_heap_usage();
  std::cout << "KDD " << myproc 
            << " ITER " << itercnt
            << " BEFORE final Comm_free:  " << oldheap
            << std::endl;
  ierr = MPI_Comm_free(&local_comm);
  newheap = get_heap_usage();
  freed = oldheap - newheap;
  if (ierr != MPI_SUCCESS) std::cout << " ERROR FREE " << ierr << std::endl;
  std::cout << "KDD " << myproc 
            << " ITER " << itercnt
            << " AFTER  final Comm_free:  " << newheap
            << " Freed: " << freed
            << " Leaked: " << used - freed << std::endl;

  if (itercnt) total_leak += (used - freed);
  itercnt++;
}

/////////////////////////////////////////////////////////////////////////////
main(int argc, char *argv[])
{
  size_t beginning = get_heap_usage();
  MPI_Init(&argc, &argv);

  //  The first MPI_Comm_split always leaks memory; get it over with and
  //  do not include it in the overall stats.
  size_t firstiterbef = get_heap_usage();
  test_function();
  size_t firstiteraft = get_heap_usage();

  size_t initheap = get_heap_usage();
  for (int i = 0; i < NUM_ITER; i++) test_function();
  size_t finalheap = get_heap_usage();

  int myproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  int localmax, globalmax;
  localmax = total_leak;
  MPI_Allreduce(&localmax, &globalmax, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);

  MPI_Finalize();
  size_t ending = get_heap_usage();

  std::cout << "KDDEND " << myproc
            << " First MPI_Comm_split leaked " << firstiteraft - firstiterbef
            << std::endl;
  std::cout << "KDDEND " << myproc
            << " Subsequent MPI_Comm_split leaked (total) "
            << finalheap - initheap << " = " << total_leak << std::endl;
  std::cout << "KDDEND " << myproc
            << " Avg per Subsequent MPI_Comm_split "
            << (finalheap - initheap) / NUM_ITER
            << std::endl;
  std::cout << "KDDEND " << myproc
            << " Max per Subsequent MPI_Comm_split "
            << globalmax
            << std::endl;
  std::cout << "KDDEND " << myproc
            << " Total Leak "
            << (ending - beginning) 
            << std::endl;



  return(0);  
}

