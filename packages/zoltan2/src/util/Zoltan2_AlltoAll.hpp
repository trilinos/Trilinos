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
  const ArrayView<const int> &sendCount, const ArrayView<int> &recvCount);

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
void AlltoAllv(const Comm<int> &comm,
              const Environment &env,
              const ArrayView<const T> &sendBuf,
              const ArrayView<const int> &sendCount,
              ArrayRCP<T> &recvBuf,      // output, allocated here
              const ArrayView<int> &recvCount   // output
)
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  try{
    Zoltan2::AlltoAllCount(comm, env, sendCount, recvCount);
  }
  Z2_FORWARD_EXCEPTIONS;

  // Allocate the receive buffer.
  size_t totalrecv = 0;
  int maxMsg = 0;
  int nrecvranks = 0;
  for(int i = 0; i < nprocs; i++) {
    if (recvCount[i] > 0) {
      totalrecv += recvCount[i];
      nrecvranks++;
      if (recvCount[i] > maxMsg) maxMsg = recvCount[i];
    }
  }


  T *rbuf = new T[totalrecv];

  if (nprocs > 1) {

    RCP<CommRequest<int> > *requests = new RCP<CommRequest<int> > [nrecvranks];

    // Error checking for memory and message size.
    int OK[2] = {1,1};
                // OK[0] -- true/false indicating whether each message size
                //          fits in an int (for MPI).
                // OK[1] -- true/false indicating whether memory allocs are OK
    int gOK[2]; // For global reduce of OK.

    if (size_t(maxMsg) * sizeof(T) > INT_MAX && nprocs > 1) OK[0] = false;
    if (totalrecv && !rbuf) OK[1] = 0;
    if (!requests) OK[1] = 0;

    // Post receives

    size_t offset = 0;
    size_t myrecvoffset = 0;
    size_t mysendoffset = 0;

    if (OK[0] && OK[1]) {
      int rcnt = 0;
      for (int i = 0; i < nprocs; i++) {
        if (i != rank && recvCount[i]) {
          try {
            requests[rcnt++] = Teuchos::ireceive<int,T>(comm,
                             Teuchos::arcp(&rbuf[offset],0,recvCount[i],false),
                             i);
          }
          Z2_THROW_OUTSIDE_ERROR(env);
        }
        else if (i == rank) {
          myrecvoffset = offset;
        }
        offset += recvCount[i];
      }
    }

    // Use barrier for error checking
    Teuchos::reduceAll<int>(comm, Teuchos::REDUCE_MIN, 2, OK, gOK);
    if (!gOK[0] || !gOK[1]) {
      delete [] rbuf;
      delete [] requests;
      if (!gOK[0])
        throw std::runtime_error("Max single message length exceeded");
      else
        throw std::bad_alloc();
    }

    // Send data; can use readySend since receives are posted.
    offset = 0;
    for (int i = 0; i < nprocs; i++) {
      if (i != rank && sendCount[i]) {
        try {
          Teuchos::readySend<int,T>(comm,
                            Teuchos::arrayView(&sendBuf[offset],sendCount[i]),
                            i);
        }
        Z2_THROW_OUTSIDE_ERROR(env);
      }
      else if (i == rank) {
        mysendoffset = offset;
      }
      offset += sendCount[i];
    }

    // Copy local data
    for (int j = 0; j < sendCount[rank]; j++)
      rbuf[myrecvoffset++] = sendBuf[mysendoffset++];

    // Wait for messages to return.
    try {
      Teuchos::waitAll<int>(comm, Teuchos::arrayView(requests, nrecvranks));
    }
    Z2_THROW_OUTSIDE_ERROR(env);

    delete [] requests;
  }
  else { // nprocs == 1; no communication needed

    if (totalrecv && !rbuf)
      throw std::bad_alloc();

    for (int j = 0; j < sendCount[0]; j++)
      rbuf[j] = sendBuf[j];
  }

  if (totalrecv)
    recvBuf = ArrayRCP<T>(rbuf, 0, totalrecv, true);
  else
    recvBuf = Teuchos::null;
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
              const ArrayView<int> &recvCount);

#ifdef HAVE_ZOLTAN2_LONG_LONG

/* \brief Specialization for unsigned long long
 */
template <>
void AlltoAllv(const Comm<int> &comm,
              const Environment &env,
              const ArrayView<const unsigned long long> &sendBuf,
              const ArrayView<const int> &sendCount,
              ArrayRCP<unsigned long long> &recvBuf,
              const ArrayView<int> &recvCount);
#endif

/* \brief Specialization for unsigned short
 */
template <>
void AlltoAllv(const Comm<int> &comm,
              const Environment &env,
              const ArrayView<const unsigned short> &sendBuf,
              const ArrayView<const int> &sendCount,
              ArrayRCP<unsigned short> &recvBuf,
              const ArrayView<int> &recvCount);

/* \brief For data type unsigned char (no Teuchos::DirectSerializationTraits)
 */
template <>
void AlltoAllv(const Comm<int> &comm,
              const Environment &env,
              const ArrayView<const unsigned char> &sendBuf,
              const ArrayView<const int> &sendCount,
              ArrayRCP<unsigned char> &recvBuf,
              const ArrayView<int> &recvCount);

}                   // namespace Z2
#endif
