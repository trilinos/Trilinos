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

/*! \file Zoltan2_AlltoAll.hpp
    \brief AlltoAll communication methods
*/

#ifndef _ZOLTAN2_ALLTOALL_HPP_
#define _ZOLTAN2_ALLTOALL_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Environment.hpp>

#include <vector>
#include <climits>

namespace Zoltan2
{

extern void AlltoAllCount(const Comm<int> &comm, const Environment &env,
 const ArrayView<const int> &sendCount, ArrayRCP<int> &recvCount);

/*! \brief AlltoAllv sends/receives data to/from all processes.
 *
 *  \param  comm   The communicator for the process group involved
 *  \param  env    The environment, required for error messages
 *  \param  sendBuf  The data to be sent, in destination process rank order
 *  \param  sendCount The number of Ts to send to process p is in sendCount[p].
 *  \param  recvBuf  On return, recvBuf has been allocated and contains
 *                   the packets sent to this process by others.
 *  \param  recvCount On return, The number of Ts received from process p 
 *                     will be in recvCount[p].
 *  \param countsAreUniform set to true if all messages sizes are the same.
 *
 * The data type T must be a type for which either
 * Zoltan2::IdentifierTraits are defined or 
 * Teuchos::DirectSerializationTraits is defined.
 *
 * AlltoAll uses only point-to-point messages.  This is to avoid the MPI 
 * limitation of integer offsets and counters in collective operations.
 * In other words, sendBuf.size() can exceed a value that fits into
 * 32 bits. However each point to point message size must fit in an int.
 *
 * It also avoids non-scalable MPI data structures that are associated
 * with collective operations.
 *
 * In addition it can be used for Zoltan2 global ID data types that are
 * not serializable by Teuchos point-to-point messages.
 */

template <typename T>
static void AlltoAllv(const Comm<int> &comm,
              const Environment &env,  
              const ArrayView<const T> &sendBuf,
              const ArrayView<const int> &sendCount,
              ArrayRCP<T> &recvBuf,      // output, allocated here
              ArrayRCP<int> &recvCount,   // output, allocated here
              bool countsAreUniform=false)
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  if (countsAreUniform){
    int *counts = new int [nprocs];
    for (int i=0; i < nprocs; i++)
      counts[i] = sendCount[0];
    recvCount = arcp(counts, 0, nprocs, true);
  }
  else{
    try{
      Zoltan2::AlltoAllCount(comm, env, sendCount, recvCount);
    }
    Z2_FORWARD_EXCEPTIONS;
  }

  size_t *offsetIn = new size_t [nprocs+1];
  size_t *offsetOut = new size_t [nprocs+1];

  ArrayRCP<size_t> offArray1(offsetIn, 0, nprocs+1, true);
  ArrayRCP<size_t> offArray2(offsetOut, 0, nprocs+1, true);
  
  offsetIn[0] = offsetOut[0] = 0;

  int maxMsg=0;
  bool offProc = false;

  for (int i=0; i < nprocs; i++){
    offsetIn[i+1] = offsetIn[i] + recvCount[i];
    offsetOut[i+1] = offsetOut[i] + sendCount[i];
    if (recvCount[i] > maxMsg)
      maxMsg = recvCount[i];
    if (sendCount[i] > maxMsg)
      maxMsg = sendCount[i];

    if (!offProc && (i != rank) && (recvCount[i] > 0 || sendCount[i] > 0))
      offProc = true;
  }

  env.globalInputAssertion(__FILE__, __LINE__,
    "message size exceeds MPI limit (sizes, offsets, counts are ints) ",
    maxMsg*sizeof(T) <= INT_MAX, BASIC_ASSERTION, rcp(&comm, false));

  size_t totalIn = offsetIn[nprocs];

  T *rptr = NULL;

  if (totalIn)
    rptr = new T [totalIn]; 
  
  env.globalMemoryAssertion(__FILE__, __LINE__, totalIn, !totalIn||rptr, 
    rcp(&comm, false));

  recvBuf = Teuchos::arcp<T>(rptr, 0, totalIn, true);

  const T *sptr = sendBuf.getRawPtr();

  // Copy self messages

  if (recvCount[rank] > 0)
    memcpy(rptr + offsetIn[rank], sptr + offsetOut[rank], 
      recvCount[rank]*sizeof(T));

  if (nprocs < 2)
    return;

#ifdef HAVE_ZOLTAN2_MPI

  // I was getting hangs in Teuchos::waitAll, so I do
  // blocking receives below.

  if (offProc){
    Array<ArrayRCP<const T> > sendArray(nprocs);
    for (int p=0; p < nprocs; p++){
      if (p != rank && sendCount[p] > 0)
        sendArray[p] = arcp(sptr + offsetOut[p], 0, sendCount[p], false);
    }
  
    for (int p=1; p < nprocs; p++){
      int recvFrom = (rank + nprocs - p) % nprocs;
      int sendTo = (rank + p) % nprocs;
  
      if (sendCount[sendTo] > 0){
        try{  // non blocking send
          Teuchos::isend<int, T>(comm, sendArray[sendTo], sendTo);
        }
        Z2_THROW_OUTSIDE_ERROR(env);
      }
  
      if (recvCount[recvFrom] > 0){
        try{  // blocking receive for message just sent to me
          Teuchos::receive<int, T>(comm, recvFrom, recvCount[recvFrom],
             rptr + offsetIn[recvFrom]);
        }
        Z2_THROW_OUTSIDE_ERROR(env);
      }
    }
  }

  comm.barrier();

#endif
}

/////////////////
/* \brief Specialization for std::string.
  
    For string of char. Number of chars in a string limited to SCHAR_MAX.
    Send as chars: 1 char for length of string, then chars in string,
     1 char for length of next string, and so on.
    \todo error checking
 */
template <>
void AlltoAllv(const Comm<int> &comm,
              const Environment &env,  
              const ArrayView<const string> &sendBuf,
              const ArrayView<const int> &sendCount,
              ArrayRCP<string> &recvBuf,  // output, allocated here
              ArrayRCP<int> &recvCount,  // output, allocated here
              bool countsAreUniform)
{
  int nprocs = comm.getSize();
  int *newCount = new int [nprocs];
  memset(newCount, 0, sizeof(int) * nprocs);
  ArrayView<const int> newSendCount(newCount, nprocs);


  size_t numStrings = sendBuf.size();
  size_t numChars = 0;
  bool fail=false;

  for (int p=0, i=0; !fail && p < nprocs; p++){
    for (int c=0; !fail && c < sendCount[p]; c++, i++){
      size_t nchars = sendBuf[i].size();
      if (nchars > SCHAR_MAX)
        fail = true;
      else
        newCount[p] += nchars;
    }
    newCount[p] += sendCount[p];
    numChars += newCount[p];
  }

  if (fail)
    throw std::runtime_error("id string length exceeds SCHAR_MAX");

  char *sbuf = NULL;
  if (numChars > 0)
    sbuf = new char [numChars];
  char *sbufptr = sbuf;

  ArrayView<const char> newSendBuf(sbuf, numChars);

  for (size_t i=0; i < numStrings; i++){
    size_t nchars = sendBuf[i].size();
    *sbufptr++ = static_cast<char>(nchars);
    for (size_t j=0; j < nchars; j++)
      *sbufptr++ = sendBuf[i][j];
  }

  ArrayRCP<char> newRecvBuf;
  ArrayRCP<int> newRecvCount;

  AlltoAllv<char>(comm, env, newSendBuf, newSendCount, 
    newRecvBuf, newRecvCount, countsAreUniform);

  delete [] sbuf;
  delete [] newCount;

  char *inBuf = newRecvBuf.getRawPtr();

  int numNewStrings = 0;
  char *buf = inBuf;
  char *endChar = inBuf + newRecvBuf.size();
  while (buf < endChar){
    int slen = static_cast<int>(*buf++);
    buf += slen;
    numNewStrings++;
  }

  // Counts to return
  int *numStringsRecv = new int [nprocs];
  memset(numStringsRecv, 0, sizeof(int) * nprocs);

  // Data to return
  string *newStrings = new string [numNewStrings];

  buf = inBuf;
  int next = 0;

  for (int p=0; p < nprocs; p++){
    int nchars = newRecvCount[p];
    endChar = buf + nchars;
    while (buf < endChar){
      int slen = *buf++;
      string nextString;
      for (int i=0; i < slen; i++)
        nextString.push_back(*buf++);
      newStrings[next++] = nextString;
      numStringsRecv[p]++;
    }
  }

  recvBuf = arcp<string>(newStrings, 0, numNewStrings, true);
  recvCount = arcp<int>(numStringsRecv, 0, nprocs, true);
}

#ifdef HAVE_ZOLTAN2_LONG_LONG

/* \brief Specialization for unsigned long long 
 */
template <>
void AlltoAllv(const Comm<int> &comm,
              const Environment &env,  
              const ArrayView<const unsigned long long> &sendBuf,
              const ArrayView<const int> &sendCount,
              ArrayRCP<unsigned long long> &recvBuf,  // output, allocated here
              ArrayRCP<int> &recvCount,  // output, allocated here
              bool countsAreUniform)
{
  const long long *sbuf = 
    reinterpret_cast<const long long *>(sendBuf.getRawPtr());
  ArrayView<const long long> newSendBuf(sbuf, sendBuf.size());
  ArrayRCP<long long> newRecvBuf;

  AlltoAllv<long long>(comm, env, newSendBuf, sendCount, 
    newRecvBuf, recvCount, countsAreUniform);

  recvBuf = arcp_reinterpret_cast<unsigned long long>(newRecvBuf);
}
#endif

/* \brief Specialization for unsigned short 
 */
template <>
void AlltoAllv(const Comm<int> &comm,
              const Environment &env,  
              const ArrayView<const unsigned short> &sendBuf,
              const ArrayView<const int> &sendCount,
              ArrayRCP<unsigned short> &recvBuf,  // output, allocated here
              ArrayRCP<int> &recvCount,  // output, allocated here
              bool countsAreUniform)
{
  const short *sbuf = reinterpret_cast<const short *>(sendBuf.getRawPtr());
  ArrayView<const short> newSendBuf(sbuf, sendBuf.size());
  ArrayRCP<short> newRecvBuf;

  AlltoAllv<short>(comm, env, newSendBuf, sendCount, 
    newRecvBuf, recvCount, countsAreUniform);

  recvBuf = arcp_reinterpret_cast<unsigned short>(newRecvBuf);
}

/* \brief For data type unsigned char (no Teuchos::DirectSerializationTraits)
 */
template <>
void AlltoAllv(const Comm<int> &comm,
              const Environment &env,  
              const ArrayView<const unsigned char> &sendBuf,
              const ArrayView<const int> &sendCount,
              ArrayRCP<unsigned char> &recvBuf,      // output, allocated here
              ArrayRCP<int> &recvCount,  // output, allocated here
              bool countsAreUniform)
{
  const char *sbuf = reinterpret_cast<const char *>(sendBuf.getRawPtr());
  ArrayView<const char> newSendBuf(sbuf, sendBuf.size());
  ArrayRCP<char> newRecvBuf;

  AlltoAllv<char>(comm, env, newSendBuf, sendCount, 
    newRecvBuf, recvCount, countsAreUniform);

  recvBuf = arcp_reinterpret_cast<unsigned char>(newRecvBuf);
}
/////////////////
}                   // namespace Z2
#endif
