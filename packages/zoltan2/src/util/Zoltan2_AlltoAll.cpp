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
    \brief AlltoAll communication methods that don't require templates, along
           with specializations.
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
 *
 *  Note:  If Teuchos::Comm adds an AlltoAll method, we should use it
 *  instead of this function.  TODO
 */

void AlltoAllCount(
  const Comm<int> &comm,      // Communicator to use for AlltoAll
  const Environment &env,     // Needed only for error handling
  const ArrayView<const int> &sendCount,   // Input:  # of data items to
                                           //         send to each process
  const ArrayView<int> &recvCount          // Output: # of data itmes to
                                           //         receive from each process
)
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  recvCount[rank] = sendCount[rank];

  if (nprocs > 1) {

    // Post receives
    RCP<CommRequest<int> > *requests = new RCP<CommRequest<int> > [nprocs];
    for (int cnt = 0, i = 0; i < nprocs; i++) {
      if (i != rank) {
        try {
          requests[cnt++] = Teuchos::ireceive<int,int>(comm,
                                                     rcp(&(recvCount[i]),false),
                                                     i);
        }
        Z2_THROW_OUTSIDE_ERROR(env);
      }
    }

    Teuchos::barrier<int>(comm);

    // Send data; can use readySend since receives are posted.
    for (int i = 0; i < nprocs; i++) {
      if (i != rank) {
        try {
          Teuchos::readySend<int,int>(comm, sendCount[i], i);
        }
        Z2_THROW_OUTSIDE_ERROR(env);
      }
    }

    // Wait for messages to return.
    try {
      Teuchos::waitAll<int>(comm, arrayView(requests, nprocs-1));
    }
    Z2_THROW_OUTSIDE_ERROR(env);

    delete [] requests;
  }
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
               ArrayRCP<string> &recvBuf,
               const ArrayView<int> &recvCount
)
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
  Array<int> newRecvCount(nprocs, 0);

  AlltoAllv<char>(comm, env, newSendBuf, newSendCount,
                  newRecvBuf, newRecvCount());

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
      recvCount[p]++;
    }
  }

  recvBuf = arcp<string>(newStrings, 0, numNewStrings, true);
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
              const ArrayView<int> &recvCount   // output
)
{
  const long long *sbuf =
    reinterpret_cast<const long long *>(sendBuf.getRawPtr());
  ArrayView<const long long> newSendBuf(sbuf, sendBuf.size());
  ArrayRCP<long long> newRecvBuf;

  AlltoAllv<long long>(comm, env, newSendBuf, sendCount,
    newRecvBuf, recvCount);

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
              const ArrayView<int> &recvCount   // output
)
{
  const short *sbuf = reinterpret_cast<const short *>(sendBuf.getRawPtr());
  ArrayView<const short> newSendBuf(sbuf, sendBuf.size());
  ArrayRCP<short> newRecvBuf;

  AlltoAllv<short>(comm, env, newSendBuf, sendCount,
    newRecvBuf, recvCount);

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
              const ArrayView<int> &recvCount   // output
)
{
  const char *sbuf = reinterpret_cast<const char *>(sendBuf.getRawPtr());
  ArrayView<const char> newSendBuf(sbuf, sendBuf.size());
  ArrayRCP<char> newRecvBuf;

  AlltoAllv<char>(comm, env, newSendBuf, sendCount,
    newRecvBuf, recvCount);

  recvBuf = arcp_reinterpret_cast<unsigned char>(newRecvBuf);
}

}                   // namespace Z2
