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

/*! \file Zoltan2_AlltoAll.cpp
    \brief AlltoAll communication methods that don't require templates.
*/

#include <Zoltan2_AlltoAll.hpp>
#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Environment.hpp>

#include <vector>
#include <climits>

namespace Zoltan2
{

/*! \brief Each process sends a value to every process, an all-to-all.
 *  \param  comm   The communicator for the process group involved
 *  \param  env    The environment, required for error messages
 *  \param  sendCount The number to send to process p is in sendCount[p].
 *  \param  recvCount On return, The number received from process p 
 *                     will be in recvCount[p].
 */

void AlltoAllCount(const Comm<int> &comm, const Environment &env,
 const ArrayView<const int> &sendCount, ArrayRCP<int> &recvCount)
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  RCP<const int> *messages = new RCP<const int> [nprocs];
  for (int p=0; p < nprocs; p++)
    messages[p] = rcp(sendCount.getRawPtr()+p, false);

  ArrayRCP<RCP<const int> > messageArray(messages, 0, nprocs, true);

  int *counts = new int [nprocs];
  recvCount = arcp(counts, 0, nprocs, true);

  counts[rank] = sendCount[rank];

#ifdef HAVE_ZOLTAN2_MPI

  // I was getting hangs in Teuchos::waitAll, so I do
  // blocking receives below.

  for (int p=1; p < nprocs; p++){
    int recvFrom = (rank + nprocs - p) % nprocs;
    int sendTo = (rank + p) % nprocs;

    try{  // non blocking send
      Teuchos::isend<int, int>(comm, messageArray[sendTo], sendTo);
    }
    Z2_THROW_OUTSIDE_ERROR(env);

    try{  // blocking receive for message just sent to me
      Teuchos::receive<int, int>(comm, recvFrom, counts + recvFrom);
    }
    Z2_THROW_OUTSIDE_ERROR(env);
  }

  comm.barrier();
#endif
}

}                   // namespace Z2
