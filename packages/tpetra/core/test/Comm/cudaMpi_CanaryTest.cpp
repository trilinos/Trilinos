/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER
*/

// This test exercises MPI on CUDA platforms, without any interaction 
// with Trilinos.  If this test fails, there is likely something wrong
// with the MPI installation on the CUDA platform.
// See Trilinos github issue #3405 and related test cudaAwareMpi.cpp
// in this directory.

#include <iostream>
#include <cstdio>
#include <mpi.h>

int main(int narg, char **arg)
{
  MPI_Init(&narg, &arg);

  int lclSuccess = 1; 
  int gblSuccess = 0; 

  std::cout << "Testing CUDA-awareness of the MPI implementation" << std::endl;

  int myRank;    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  int numProcs;  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if (numProcs < 2) {
    std::cout << "This test is more meaningful if run with at least 2 MPI "
              << "processes.  You ran with only 1 MPI process." << std::endl;
  }

  // Create Views for sending and receiving data.
  int length = 1;

  int *sendBuf, *recvBuf;
  cudaMalloc(&sendBuf, length * sizeof(int));
  cudaMalloc(&recvBuf, length * sizeof(int));

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
      cudaMemcpy(sendBuf, sendBuf_h, length*sizeof(int), cudaMemcpyHostToDevice);

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
      cudaMemcpy(recvBuf, recvBuf_h, length*sizeof(int), cudaMemcpyHostToDevice);

      const int srcRank = 0;

      std::cout << "     " << myRank << " Recv " << std::endl;
      MPI_Request recvReq;
      MPI_Irecv(recvBuf, length, MPI_INT, srcRank, msgTag, 
                MPI_COMM_WORLD, &recvReq);

      std::cout << "     " << myRank << " Wait " << std::endl;
      MPI_Status stat;
      MPI_Wait(&recvReq, &stat);

      cudaMemcpy(recvBuf_h, recvBuf, length*sizeof(int), cudaMemcpyDeviceToHost);
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
    MPI_Allreduce(&lclSuccess, &gblSuccess, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (gblSuccess != 1)
      std::cout << "Neighbor test failed on some process!" << std::endl;
    else 
      std::cout << "Neighbor test succeeded on all processes!" << std::endl;
  }

  
  // Self-message test
  {
    std::cout << myRank << " SELF-MSG TEST" << std::endl;

    const int msgTag = 42; // just some nontrivial tag for MPI messages

    // Fill send buffer with some unique positive value.
    const int correctValue = myRank + 1;
    for (int i = 0; i < length; i++) sendBuf_h[i] = correctValue;
    cudaMemcpy(sendBuf, sendBuf_h, length*sizeof(int), cudaMemcpyHostToDevice);

    // Fill receive buffer with a flag value, always negative.
    const int flagValue = -(myRank + 1);
    for (int i = 0; i < length; i++) recvBuf_h[i] = flagValue;
    cudaMemcpy(recvBuf, recvBuf_h, length*sizeof(int), cudaMemcpyHostToDevice);

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

    cudaMemcpy(recvBuf_h, recvBuf, length*sizeof(int), cudaMemcpyDeviceToHost);
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
    MPI_Allreduce(&lclSuccess, &gblSuccess, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   
    if (gblSuccess != 1)
      std::cout << "Self-message test failed on some process; !" << std::endl;
    else
      std::cout << "Self-message test succeeded on all processes!" << std::endl;
  }

  cudaFree(sendBuf);
  cudaFree(recvBuf);

  free(sendBuf_h);
  free(recvBuf_h);

  MPI_Finalize();
  return 0;
}
