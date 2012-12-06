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


}                   // namespace Z2
