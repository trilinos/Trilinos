
///////////////////////////////////////////////////////////////
// Program to test MPI_Comm_dup memory usage on thunderbird. //
///////////////////////////////////////////////////////////////

#include <iostream>
#include <mpi.h>
#include "get_heap_usage.h"

#define NUM_ITER 10

/////////////////////////////////////////////////////////////////////////////
void test_function()
{
MPI_Comm local_comm, tmp_comm;
int myproc, nprocs;               // MPI info wrt MPI_COMM_WORLD.
size_t oldheap, newheap;
static int itercnt = 0;

  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  //  Duplicate MPI_COMM_WORLD to local communicator.
  oldheap = get_heap_usage();
  std::cout << "KDD " << myproc 
            << " ITER " << itercnt
            << " BEFORE Comm_dup:  " << oldheap << std::endl;
  MPI_Comm_dup(MPI_COMM_WORLD,&local_comm);
  newheap = get_heap_usage();
  std::cout << "KDD " << myproc 
            << " ITER " << itercnt
            << " AFTER  Comm_dup:  " << newheap 
            << " Used: " << newheap - oldheap << std::endl;

  // Free local_comm.
  oldheap = get_heap_usage();
  std::cout << "KDD " << myproc 
            << " ITER " << itercnt
            << " BEFORE final Comm_free:  " << oldheap
            << std::endl;
  MPI_Comm_free(&local_comm);
  newheap = get_heap_usage();
  std::cout << "KDD " << myproc 
            << " ITER " << itercnt
            << " AFTER  final Comm_free:  " << newheap
            << " Freed: " << oldheap - newheap << std::endl;

  itercnt++;
}

/////////////////////////////////////////////////////////////////////////////
main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  //  The first MPI_Comm_dup always leaks memory; get it over with and
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
            << " First MPI_Comm_dup leaked " << firstiteraft - firstiterbef
            << std::endl;
  std::cout << "KDDEND " << myproc 
            << " Subsequent MPI_Comm_dups leaked (total) " 
            << finalheap - initheap << std::endl;
  std::cout << "KDDEND " << myproc 
            << " Avg per Subsequent MPI_Comm_dup " 
            << (finalheap - initheap) / NUM_ITER
            << std::endl;

  return(0);  
}

