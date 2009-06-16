#ifndef KOKKOS_CUDANODE_CUH_
#define KOKKOS_CUDANODE_CUH_

#include <cuda.h>
#include <Kokkos_CUDA_util_inline_runtime.h>
#include <Kokkos_CUDA_util_sharedmem.cuh>

// must define this before including any kernels
#define KERNEL_PREFIX __device__ __host__

#include <Kokkos_CUDANode.hpp>

// The code for the reduction was taken from Mark Harris's reduction demo in the the CUDA SDK.
// This is not the most efficient version of the kernel; for simplicity, we used the slightly less efficient reduce5 kernel.
// It has been modified to run for arrays whose length is not a power of two (making it even less efficient, no doubt).
// At some point in the future, the kernel needs to be updated to a more efficient version.
// We also use the sharedmem.cuh functionality provided by that example and the cutil macro collection.


#ifdef CUDANODE_INCLUDE_EXECUTE1D
template <class WDP>
__global__ void
Tkern1D(int length, WDP wd, int h)
{
  unsigned int b = h*(blockIdx.x*blockDim.x + threadIdx.x);
  unsigned int t = b+h;
  if (t < length) 
    wd(b,t);
  else 
    wd(b,length);
}

template <class WDP>
void CUDANode::execute1D(int length, WDP wd) {
  if (length == 0) return;
  unsigned int h = length / (numThreads_ * numBlocks_);
  if ((length % (numThreads_*numBlocks_)) != 0) {
    h = h + 1;
  }
  Tkern1D<WDP> <<< numBlocks_, numThreads_ >>>(length,wd,h);
}
#endif // execute1D

#ifdef CUDANODE_INCLUDE_REDUCE1D
template <class WDP, int FirstLevel, int blockSize>
__global__ void
Treduce1D(int length, WDP wd, void *d_blkpart)
{
  typedef typename WDP::ReductionType T;
  T * out = (T*)d_blkpart;
  SharedMemory<T> smem;
  T *sdata = smem.getPointer();
  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*(blockSize*2) + threadIdx.x; // each thread gets two elements
  // each thread does the first pair
  if (FirstLevel == 1) {
    // we have either zero, two or one elements to process
    T gi = wd.identity(), gip = wd.identity();
    if (i < length) {
      gi = wd.generate(i);
      if (i+blockSize < length) {
        gip = wd.generate(i+blockSize);
      }
    }
    sdata[tid] = wd.reduce( gi, gip );
  }
  else {
    sdata[tid] = wd.reduce( out[i], out[i+blockSize] );
  }
  __syncthreads();
  // finish the fan-in with a decreasing number of threads
  if (blockSize >= 512) { if (tid < 256) { sdata[tid] = wd.reduce(sdata[tid], sdata[tid + 256]); } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] = wd.reduce(sdata[tid], sdata[tid + 128]); } __syncthreads(); }
  if (blockSize >= 128) { if (tid <  64) { sdata[tid] = wd.reduce(sdata[tid], sdata[tid +  64]); } __syncthreads(); }
  if (tid < 32)
  {
    if (blockSize >= 64) { sdata[tid] = wd.reduce( sdata[tid], sdata[tid + 32] ); __syncthreads(); }
    if (blockSize >= 32) { sdata[tid] = wd.reduce( sdata[tid], sdata[tid + 16] ); __syncthreads(); }
    if (blockSize >= 16) { sdata[tid] = wd.reduce( sdata[tid], sdata[tid +  8] ); __syncthreads(); }
    if (blockSize >=  8) { sdata[tid] = wd.reduce( sdata[tid], sdata[tid +  4] ); __syncthreads(); }
    if (blockSize >=  4) { sdata[tid] = wd.reduce( sdata[tid], sdata[tid +  2] ); __syncthreads(); }
    if (blockSize >=  2) { sdata[tid] = wd.reduce( sdata[tid], sdata[tid +  1] ); __syncthreads(); }
  }
  // write result for this block to global mem
  if (tid == 0) out[blockIdx.x] = sdata[0];
}

template <class WDP, int FirstLevel>
void CUDANode::call_reduce(int length, WDP wd, int threads, int blocks, void *d_blkpart)
{
  typedef typename WDP::ReductionType T;
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);
  int smemSize = threads * sizeof(T);
  switch (threads)
  {
    case 512:
      Treduce1D<WDP, FirstLevel, 512> <<< dimGrid, dimBlock, smemSize >>>(length, wd, d_blkpart); break;
    case 256:
      Treduce1D<WDP, FirstLevel, 256> <<< dimGrid, dimBlock, smemSize >>>(length, wd, d_blkpart); break;
    case 128:
      Treduce1D<WDP, FirstLevel, 128> <<< dimGrid, dimBlock, smemSize >>>(length, wd, d_blkpart); break;
    case 64:
      Treduce1D<WDP, FirstLevel,  64> <<< dimGrid, dimBlock, smemSize >>>(length, wd, d_blkpart); break;
    case 32:
      Treduce1D<WDP, FirstLevel,  32> <<< dimGrid, dimBlock, smemSize >>>(length, wd, d_blkpart); break;
    case 16:
      Treduce1D<WDP, FirstLevel,  16> <<< dimGrid, dimBlock, smemSize >>>(length, wd, d_blkpart); break;
    case  8:
      Treduce1D<WDP, FirstLevel,   8> <<< dimGrid, dimBlock, smemSize >>>(length, wd, d_blkpart); break;
    case  4:
      Treduce1D<WDP, FirstLevel,   4> <<< dimGrid, dimBlock, smemSize >>>(length, wd, d_blkpart); break;
    case  2:
      Treduce1D<WDP, FirstLevel,   2> <<< dimGrid, dimBlock, smemSize >>>(length, wd, d_blkpart); break;
    case  1:
      Treduce1D<WDP, FirstLevel,   1> <<< dimGrid, dimBlock, smemSize >>>(length, wd, d_blkpart); break;
  }
  cutilCheckMsg("Kernel execution failed");
}

template <class WDP>
void CUDANode::reduce1D(int length, WDP &wd) 
{
  if (length < 2) {
    if (length == 0) {
      wd.result = wd.identity();
    }
    else {
      wd.result = wd.generate(0);
    }
    return;
  }
  const int FIRST_LEVEL = 1;
  const int NOT_FIRST_LEVEL = 0;
  // we will do multiple levels of reduction
  // level 1: reduce array of size length using lclNumBlocks blocks and numThreads_ threads
  //          level 1 is implemented by a single call to Treduce1D
  //          a fan-in in every block among all threads reduces the entire array to lclNumBlocks entries
  //        INPUT:
  //          the array to be reduced is implicitly stored in the WDP object: generate(i) produces the i-th element
  //          the length is not necessarily a power-of-two, so the kernel routine will make allowances for this
  //          in the first level.
  //       OUTPUT:
  //          the result of level 1 is a global device-memory array of size lclNumBlocks, where the i-th element
  //          was produced by the i-th block. this array is called d_blkpart, and is allocated and deallocated 
  //          here.
  //       We will implicitly expand the first level so that the number of entries processed is the first power-of-two
  //       as large as the reduction length, i.e., min n such that length <= 2^n.
  // levels 2 through L: 
  //          next, the entries (of quantity lclNumBlocks) in d_blkpart must be reduced.
  //          we will use further calls to Treduce1D, where the input data comes now from the d_blkpart array
  //          (this distinction is indicated by a template argument to the Treduce1D kernel.)
  //          because there are a maximum number of 512 threads/block (per the CUDA SDK), the number of partial 
  //          block results after Level 1 reduction is effectively dictated (for this kernel) by the length of the 
  //          input array. because this may still be greater than the maximum allowable number of threads per blocks, 
  //          we may need multiple reductions from here.
  // CUDANode guarantees that numThreads_ is a power-of-two; this is the maximum number of threads that we will use. 
  // 
  typedef typename WDP::ReductionType T;
  // determine the length for the first level; "round up" length to the next power of two.
  int po2len = 1;
  while (po2len < length) po2len <<= 1; 
  int nT = (numThreads_*2 >= po2len) ? po2len/2 : numThreads_;
  // nB * 2 * nT =must be= po2len, or the work won't get done
  int nB = po2len / (2*nT); // this is a power-of-two by construction
  T *d_blkpart = 0;
  cutilSafeCallNoSync( cudaMalloc((void**) &d_blkpart, nB*sizeof(T)) );    // consider saving this for later use
  // do the first level of reduction, condensing the data 
  call_reduce<WDP,FIRST_LEVEL>(length,wd,nT,nB,(void*)d_blkpart);
  // now reduce the data in d_blkpart
  po2len = nB;
  while (po2len > 1) {
    nT = (numThreads_*2 >= po2len) ? po2len/2 : numThreads_;
    nB = po2len / (2*nT);
    // call the reduction; this is no longer the first level; work out of and into d_blkpart
    call_reduce<WDP,NOT_FIRST_LEVEL>(po2len,wd,nT,nB,(void*)d_blkpart);
    po2len = po2len / (2*nT);
  }
  // the answer is in d_blkpart[0]; get it
  cutilSafeCallNoSync( cudaMemcpy( &wd.result, d_blkpart, sizeof(T), cudaMemcpyDeviceToHost) );
  // free d_blkpart
  cutilSafeCallNoSync( cudaFree(d_blkpart) );
}
#endif // reduce1D

#endif
