
#include "stk_util/parallel/Parallel.hpp" 
#include <vector> 

//
//  Take a list vector of T's on each processor.  Sum it into a single list that will be placed on all processors.
//  The list contents and order will be guaranteed identical on every processor and formed by concatenation of the list
//  in processor order.
//
//  Return Code:  An MPI error code, MPI_SUCESS if correct
//
template <typename T> inline int parallelConcat(ParallelMachine comm, const std::vector<T> & localVec, std::vector<T> & globalVec ) {
  globalVec.clear();
  
  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  unsigned int sizeT     = sizeof(T);
  unsigned int localSize = sizeT * localVec.size();

  //
  //  Determine the total number of bytes being sent by each other processor.
  //
  vector<Int> messageSizes(p_size);
  int mpiResult = MPI_SUCCESS ;  
  mpiResult = MPI_AllGather(&localSize, 1, MPI_INT, &messageSizes[0], 1, MPI_INT, 0, comm);
  if(mpiResult != MPI_SUCCESS) {
    // Unknown failure, pass error code up the chain
    return mpiResult;
  }

  //
  //  Compute the offsets into the resultant array
  //
  vector<Int> offsets(p_size+1);
  offsets[0] = 0;
  for(Int iproc=1; iproc<p_size+1; ++iproc) {
    offsets[iproc] = offsets[iproc-1] + messagesSizes[iproc-1];
  }
  unsigned int totalSize = (offsets[p_size])/sizeT;
  globalVec.resize(totalSize);

  //
  //  Do the all gather to copy the actual array data and propogate to all processors
  //
  mpiResult = MPI_AllGatherv(&localVec[0], &localSize, MPI_CHAR, &globalVec[0], &messageSizes[0], &offsets[0], MPI_CHAR, comm);
  return mpiResult;
}
