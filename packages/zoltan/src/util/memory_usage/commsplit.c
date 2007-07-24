
/////////////////////////////////////////////////////////////////
// Program to test MPI_Comm_split memory usage on thunderbird. //
/////////////////////////////////////////////////////////////////

#include <iostream>
#include <mpi.h>

#include "get_heap_usage.h"

#define NUM_ITER 10

/////////////////////////////////////////////////////////////////////////////
void test_function()
{
MPI_Comm local_comm;
int myproc, nprocs;               // MPI info wrt MPI_COMM_WORLD.
int set;
size_t oldheap, newheap;
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
  if (ierr != MPI_SUCCESS) std::cout << " ERROR SPLIT " << ierr << std::endl;
  std::cout << "KDD " << myproc 
            << " ITER " << itercnt
            << " AFTER  Comm_split:  " << newheap 
            << " Used: " << newheap - oldheap << std::endl;

  // Free local_comm.
  oldheap = get_heap_usage();
  std::cout << "KDD " << myproc 
            << " ITER " << itercnt
            << " BEFORE final Comm_free:  " << oldheap
            << std::endl;
  ierr = MPI_Comm_free(&local_comm);
  newheap = get_heap_usage();
  if (ierr != MPI_SUCCESS) std::cout << " ERROR FREE " << ierr << std::endl;
  std::cout << "KDD " << myproc 
            << " ITER " << itercnt
            << " AFTER  final Comm_free:  " << newheap
            << " Freed: " << oldheap - newheap  << std::endl;

  itercnt++;
}

/////////////////////////////////////////////////////////////////////////////
main(int argc, char *argv[])
{
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

  MPI_Finalize();
  std::cout << "KDDEND " << myproc
            << " First MPI_Comm_split leaked " << firstiteraft - firstiterbef
            << std::endl;
  std::cout << "KDDEND " << myproc
            << " Subsequent MPI_Comm_split leaked (total) "
            << finalheap - initheap << std::endl;
  std::cout << "KDDEND " << myproc
            << " Avg per Subsequent MPI_Comm_split "
            << (finalheap - initheap) / NUM_ITER
            << std::endl;

  return(0);  
}

