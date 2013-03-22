// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

// TODO: doxygen comments
//     make this a real unit test that gives helpful information if it fails
//     and uses different template values

#include <Zoltan2_Environment.hpp>   
#include <Zoltan2_AlltoAll.hpp>   
#include <Zoltan2_TestHelpers.hpp>   

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>

#include <Teuchos_RCP.hpp>   
#include <Teuchos_ArrayRCP.hpp>   
#include <Teuchos_Comm.hpp>   
#include <Teuchos_DefaultComm.hpp>   

using namespace std;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = 
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();

  Teuchos::RCP<const Zoltan2::Environment> envPtr = 
    Teuchos::rcp(new Zoltan2::Environment);

  int errcode = 0;

  if (!errcode){

    // test of Zoltan2::AlltoAllv (different sized messages) using a Scalar type

    int myMsgSizeBase=rank*nprocs + 1;
    Array<int> sendCount(nprocs, 0);
    Array<int> recvCount(nprocs, 0);
    long totalOut = 0;

    for (int p=0; p < nprocs; p++){
      sendCount[p] = myMsgSizeBase + p;
      totalOut += sendCount[p];
    }
  
    Array<int> sendBuf(totalOut, 0);
  
    int *out = &(sendBuf[0]);
    for (int p=0; p < nprocs; p++){
      for (int i=0; i < sendCount[p]; i++){
        *out++ = p+rank;
      }
    }
  
    Teuchos::ArrayRCP<int> recvBuf;
    
    Zoltan2::AlltoAllv<int>(*comm, *envPtr,
                  sendBuf(),
                  sendCount(),
                  recvBuf,
                  recvCount());

    int *inBuf = recvBuf.get();

    for (int p=0; p < nprocs; p++){
      for (int i=0; i < recvCount[p]; i++){
        if (*inBuf++ != rank+p){
          errcode = 4;
          break;
        }
      }
    }
  }

  TEST_FAIL_AND_EXIT(*comm, errcode==0, "int", errcode);

  if (!errcode){

    // test of Zoltan2::AlltoAllv using strings - which can not
    //    be serialized by Teuchos.
    //  Rank p sends p messages to each process.

    int nstrings = nprocs * rank;
    string *sendStrings = NULL;

    if (nstrings > 0)
      sendStrings = new string [nstrings];

    ostringstream myMessage;
    myMessage << "message from process " << rank;

    for (int i=0; i < nstrings; i++)
      sendStrings[i] = myMessage.str();

    int *counts = new int [nprocs];
    for (int i=0; i < nprocs ; i++)
      counts[i] = rank;

    Teuchos::ArrayView<const string> sendBuf(sendStrings, nstrings);
    Teuchos::ArrayView<const int> sendCount(counts, nprocs);
    Teuchos::Array<int> recvCounts(nprocs, 0);

    Teuchos::ArrayRCP<string> recvBuf;

    Zoltan2::AlltoAllv<string>(*comm, *envPtr,
        sendBuf,    
        sendCount,   
        recvBuf,
        recvCounts());

    delete [] sendStrings;
    delete [] counts;
  
    int next = 0;
    for (int i=0; i < nprocs; i++){
      if (recvCounts[i] != i){
        errcode = 5;
        break;
      }
      ostringstream msg;
      msg << "message from process " << i;
      for (int j=0; j < recvCounts[i]; j++){
        if (recvBuf[next++] != msg.str()){
          errcode = 6;
          break;
        }
      }
    }
  }

  TEST_FAIL_AND_EXIT(*comm, errcode==0, "strings", errcode);

  if (rank == 0){
    if (errcode)
      std::cout << "FAIL" << std::endl;
    else
      std::cout << "PASS" << std::endl;
  }

  return errcode;
}
