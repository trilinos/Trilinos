// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This test exercises MPI on CUDA platforms, without any interaction 
// with Trilinos.  If this test fails, there is likely something wrong
// with the MPI installation on the CUDA platform.
// See Trilinos github issue #3405 and related test cudaAwareMpi.cpp
// in this directory.

#include <iostream>
#include <cstdio>
#include <mpi.h>

// This is not a Tpetra-based test, so I regret the Tpetra header file.
// But since Tpetra allows Tpetra_ASSUME_MPI_IS_GPU_AWARE to be set 
// either as a Cmake option or as an environment variable, it is easiest
// to use Tpetra's mechanism to extract the value.
#include "Tpetra_Details_Behavior.hpp"

#if defined(__NVCC__)
#include <cuda.h>
#define TEST_MALLOC(ptr, size) cudaMalloc(ptr,size)
#define TEST_MEMCPY(dst, src, size, flag) cudaMemcpy(dst,src,size,cuda##flag)
#define TEST_FREE(ptr) cudaFree(ptr)

#elif defined(__HIP__)
#include <hip/hip_runtime.h>
#define TEST_MALLOC(ptr, size) hipMalloc(ptr,size)
#define TEST_MEMCPY(dst, src, size, flag) hipMemcpy(dst,src,size,hip##flag)
#define TEST_FREE(ptr) hipFree(ptr)

#else

#error "This test can only be build when either CUDA or HIP is enabled"

#endif

int main(int narg, char **arg)
{
  MPI_Init(&narg, &arg);
  int myRank;    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  int numProcs;  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  int lclSuccess = 1; 
  int gblSuccess = 0; 
  int nfailures = 0;

  if (myRank == 0)
    std::cout << "Testing GPU-awareness of the MPI implementation" 
              << std::endl;

  if (!Tpetra::Details::Behavior::assumeMpiIsGPUAware()) {
    // MPI is not CUDA aware, so there is nothing to test here.
    // Declare success and go home
    if (myRank == 0) 
      std::cout << "PASS:  Not using GPU-aware MPI; no need to run tests."
                << std::endl;
    MPI_Finalize();
    return 0;
  }

  if (myRank == 0) {
    std::cout << "GPU-aware MPI is detected; beginning to run tests." 
              << std::endl;
    if (numProcs < 2)
      std::cout << "This test is more meaningful if run with at least 2 MPI "
                << "processes.  You ran with only 1 MPI process." << std::endl;
  }

  // Create Views for sending and receiving data.
  int length = 1;

  int *sendBuf, *recvBuf;
  TEST_MALLOC(&sendBuf, length * sizeof(int));
  TEST_MALLOC(&recvBuf, length * sizeof(int));

  int *sendBuf_h = (int*)malloc(length * sizeof(int));
  int *recvBuf_h = (int*)malloc(length * sizeof(int));

  // NEIGHBOR TEST
  if (numProcs > 1) {
    std::cout << myRank << " NEIGHBOR TEST" << std::endl;

    const int msgTag = 43;
    const int correctValue = 1;
    const int flagValue = -1;

    //We exercise only the first two processes. 
        
    if (myRank == 0) { // sending process

      for (int i = 0; i < length; i++) sendBuf_h[i] = correctValue;
      TEST_MEMCPY(sendBuf, sendBuf_h, length*sizeof(int), MemcpyHostToDevice);

      const int tgtRank = 1;

      std::cout << "     " << myRank << " Send " << std::endl;
      MPI_Request sendReq;
      MPI_Isend(sendBuf, length, MPI_INT, tgtRank, msgTag, 
                MPI_COMM_WORLD, &sendReq);

      std::cout << "     " << myRank << " Wait " << std::endl;
      MPI_Status stat;
      MPI_Wait(&sendReq, &stat);
    }
    else if (myRank == 1) { // receiving process

      for (int i = 0; i < length; i++) recvBuf_h[i] = flagValue;
      TEST_MEMCPY(recvBuf, recvBuf_h, length*sizeof(int), MemcpyHostToDevice);

      const int srcRank = 0;

      std::cout << "     " << myRank << " Recv " << std::endl;
      MPI_Request recvReq;
      MPI_Irecv(recvBuf, length, MPI_INT, srcRank, msgTag, 
                MPI_COMM_WORLD, &recvReq);

      std::cout << "     " << myRank << " Wait " << std::endl;
      MPI_Status stat;
      MPI_Wait(&recvReq, &stat);

      TEST_MEMCPY(recvBuf_h, recvBuf, length*sizeof(int), MemcpyDeviceToHost);
      if (recvBuf_h[0] == correctValue)  {
        lclSuccess = 1;
        std::cout << "     " << myRank << " Check Result:  OK " << std::endl;
      }
      else {
        lclSuccess = 0;
        std::cout << "     " << myRank << " Check Result:  BAD " << std::endl;
      }
    }

    gblSuccess = 0; 
    MPI_Allreduce(&lclSuccess, &gblSuccess, 1, MPI_INT, MPI_MIN,MPI_COMM_WORLD);
    if (gblSuccess != 1) {
      if (myRank == 0)
        std::cout << "Neighbor test failed on some process!" << std::endl;
      nfailures++;
    }
    else 
      if (myRank == 0)
        std::cout << "Neighbor test succeeded on all processes!" << std::endl;
  }

  
  // Self-message test
  {
    std::cout << myRank << " SELF-MSG TEST" << std::endl;

    const int msgTag = 42; // just some nontrivial tag for MPI messages

    // Fill send buffer with some unique positive value.
    const int correctValue = myRank + 1;
    for (int i = 0; i < length; i++) sendBuf_h[i] = correctValue;
    TEST_MEMCPY(sendBuf, sendBuf_h, length*sizeof(int), MemcpyHostToDevice);

    // Fill receive buffer with a flag value, always negative.
    const int flagValue = -(myRank + 1);
    for (int i = 0; i < length; i++) recvBuf_h[i] = flagValue;
    TEST_MEMCPY(recvBuf, recvBuf_h, length*sizeof(int), MemcpyHostToDevice);

    const int srcRank = myRank; // self-messages
    const int tgtRank = myRank;

    std::cout << "     " << myRank << " Send/Recv " << std::endl;
    MPI_Request recvReq, sendReq;
    MPI_Irecv(recvBuf, 1, MPI_INT, srcRank, msgTag, 
              MPI_COMM_WORLD, &recvReq);
    MPI_Isend(sendBuf, 1, MPI_INT, tgtRank, msgTag, 
              MPI_COMM_WORLD, &sendReq);

    std::cout << "     " << myRank << " Wait " << std::endl;
    MPI_Status stat;
    MPI_Wait(&sendReq, &stat);
    MPI_Wait(&recvReq, &stat);

    TEST_MEMCPY(recvBuf_h, recvBuf, length*sizeof(int), MemcpyDeviceToHost);
    if (recvBuf_h[0] == correctValue)  {
      lclSuccess = 1;
      std::cout << "     " << myRank << " Check Result:  OK " << std::endl;
    }
    else {
      lclSuccess = 0;
      std::cout << "     " << myRank << " Check Result:  BAD " << std::endl;
    }

    // Make sure that everybody finished and got the right answer.
    gblSuccess = 0; 
    MPI_Allreduce(&lclSuccess, &gblSuccess, 1, MPI_INT, MPI_MIN,MPI_COMM_WORLD);
   
    if (gblSuccess != 1) {
      if (myRank == 0)
        std::cout << "Self-message test failed on some process; !" << std::endl;
      nfailures++;
    }
    else
      if (myRank == 0)
        std::cout << "Self-message test succeeded on all processes!" 
                  << std::endl;
  }

  if (myRank == 0) {
    if (nfailures > 0) 
      std::cout << "FAIL:  nfailures = " << nfailures << std::endl;
    else
      std::cout << "PASS" << std::endl;
  }

  TEST_FREE(sendBuf);
  TEST_FREE(recvBuf);

  free(sendBuf_h);
  free(recvBuf_h);

  MPI_Finalize();
  return 0;
}
