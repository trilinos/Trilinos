// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
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

/*! \brief AlltoAll sends/receives a fixed number of objects to/from all processes.
 *
 *  \param  comm   The communicator for the process group involved
 *  \param  env    The environment, required for error messages
 *  \param  sendBuf  The data to be sent, in destination process rank order
 *  \param  count    The number of Ts to sendBuf to send to each process.
 *                   This must be the same on all processes.
 *  \param  recvBuf  On return, recvBuf has been allocated and contains
 *                   the packets sent to this process by others.
 *
 * The data type T must be a type for which SerializationTraits are defined
 * in Teuchos.
 *
 * AlltoAll uses only point-to-point messages.  This is to avoid the MPI 
 * limitation of integer offsets and counters in collective operations.
 * In other words, LNO can be a 64-bit integer.  However each point to
 * point message size must fit in an int.
 *
 * It also avoids non-scalable MPI data structures that are associated
 * with collective operations.
 *
 * So this is slow, but will not encounter MPI resource limits for very 
 * large applications, some of which have already been encountered by
 * Zoltan users.
 */

template <typename T, typename LNO>
void AlltoAll(const Comm<int> &comm,
              const Environment &env,
              const ArrayView<const T> &sendBuf,
              LNO count,
              ArrayRCP<T> &recvBuf)         // output - allocated here
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  if (count == 0) return;   // count is the same on all procs

  LNO n = nprocs * count;
  T *rptr = new T [n]; 
  env.globalMemoryAssertion(__FILE__, __LINE__, n, rptr, rcp(&comm, false));
  recvBuf = Teuchos::arcp<T>(rptr, 0, n);

  const T *sptr = sendBuf.getRawPtr();

  // Do self messages

  for (LNO i=0, offset = rank*count; i < count; i++, offset++){
    rptr[offset] = sptr[offset];
  }

#ifdef HAVE_ZOLTAN2_MPI
  // Perform nprocs-1 point-to-point sends and receives.

  size_t packetSize = count * sizeof(T);
  env.globalInputAssertion(__FILE__, __LINE__,
      "message size exceeds MPI limit (sizes, offsets, counts are ints) ",
      packetSize <= INT_MAX, BASIC_ASSERTION, rcp(&comm, false));

  Array<ArrayRCP<const T> > sendArray(nprocs);
  for (int p=0; p < nprocs; p++){
    sendArray[p] = arcp(sptr + p*count, 0, count, false);
  }
  
  for (int p=1; p < nprocs; p++){
    int recvFrom = (rank + nprocs - p) % nprocs;
    int sendTo = (rank + p) % nprocs;

    try{  // Non blocking send
      Teuchos::isend<int, T>(comm, sendArray[sendTo], sendTo);
    }
    Z2_THROW_OUTSIDE_ERROR(env);

    try{  // blocking receive for msg just sent to me
      Teuchos::receive<int, T>(comm, recvFrom, count, rptr + recvFrom*count);
    }
    Z2_THROW_OUTSIDE_ERROR(env);
  }

  comm.barrier();
#endif
}

/*! \brief AlltoAllv sends/receives a variable number of objects to/from all processes.
 *
 *  \param  comm   The communicator for the process group involved
 *  \param  env    The environment, required for error messages
 *  \param  sendBuf  The data to be sent, in destination process rank order
 *  \param  sendCount The number of Ts to send to process p is in sendCount[p].
 *  \param  recvBuf  On return, recvBuf has been allocated and contains
 *                   the packets sent to this process by others.
 *  \param  recvCount On return, The number of Ts received from process p 
 *                     will be in recvCount[p].
 */

template <typename T, typename LNO>
void AlltoAllv(const Comm<int> &comm,
              const Environment &env,  
              const ArrayView<const T> &sendBuf,
              const ArrayView<const LNO> &sendCount,
              ArrayRCP<T> &recvBuf,      // output, allocated here
              ArrayRCP<LNO> &recvCount)  // output, allocated here
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  try{
    AlltoAll<LNO, LNO>(comm, env, sendCount, 1, recvCount);
  }
  Z2_FORWARD_EXCEPTIONS;

  size_t *offsetIn = new size_t [nprocs+1];
  size_t *offsetOut = new size_t [nprocs+1];
  
  offsetIn[0] = offsetOut[0] = 0;

  size_t maxMsg=0;

  for (int i=0; i < nprocs; i++){
    offsetIn[i+1] = offsetIn[i] + recvCount[i];
    offsetOut[i+1] = offsetOut[i] + sendCount[i];
    if (recvCount[i] > maxMsg)
      maxMsg = recvCount[i];
    if (sendCount[i] > maxMsg)
      maxMsg = recvCount[i];
  }

  env.globalInputAssertion(__FILE__, __LINE__,
    "message size exceeds MPI limit (sizes, offsets, counts are ints) ",
    maxMsg <= INT_MAX, BASIC_ASSERTION, rcp(&comm, false));

  size_t totalIn = offsetIn[nprocs];

  T *rptr = NULL;
  if (totalIn){
    rptr = new T [totalIn]; 
  }
  env.globalMemoryAssertion(__FILE__, __LINE__, totalIn, !totalIn||rptr, 
    rcp(&comm, false));

  recvBuf = Teuchos::arcp<T>(rptr, 0, totalIn, true);

  const T *sptr = sendBuf.getRawPtr();

  // Copy self messages

  memcpy(rptr + offsetIn[rank], sptr + offsetOut[rank], 
    sizeof(T) * recvCount[rank]);

  if (nprocs < 2)
    return;

#ifdef HAVE_ZOLTAN2_MPI

  Array<ArrayRCP<const T> > sendArray(nprocs);
  for (int p=0; p < nprocs; p++){
    if (sendCount[p] > 0)
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

  comm.barrier();
#endif
}

/*! \brief Serialization for std::vector<T>
 *
 * Teuchos::SerializationTraits exist for types that can be
 * sent in a Teuchos message, such as std::pair<T1, T2>.  It
 * does not exist for std::vector<T>, and it seems the serialization
 * interface does not make this possible.
 * 
 * These four methods are what the SerializationTraits interface
 * might look like if it could support type std::vector<T>.
 *
 * Buffer layout for std::vector<T> of size N, variable length vectors -
 *      LNO numberOfVectors      
 *      LNO offsetToStartsOfVectorElements[N] 
 *      T first element of first vector
 *        ...
 *      T last element of last vector
 *
 * Buffer layout for std::vector<T> of size N, identical length vectors -
 *      LNO numberOfVectors      
 *      T first element of first vector
 *        ...
 *      T last element of last vector
 *
 * Important: number of bytes returned is always a multiple of sizeof(T)
 */

template <typename T, typename LNO>
  LNO fromObjectsToIndirectBytes(const LNO count, 
    std::vector<T> const v[], 
    LNO vLen=0)   // set vLen to vector length if all are the same length
{
  LNO nelements=0, nbytes=0, preamble=0;

  if (vLen == 0){
    for (LNO i=0; i < count; i++)
      nelements += v[i].size();
    preamble = sizeof(LNO) * (1 + count);
  }
  else{
    nelements = vLen * count;
    preamble = sizeof(LNO);
  }

  nbytes = preamble + (preamble % sizeof(T));  // T alignment

  nbytes += nelements * sizeof(T);

  return nbytes;
}

template <typename T, typename LNO>
  void serialize(
    const LNO count, const std::vector<T> v[], const LNO bytes, char buf[], 
      LNO vLen=0)
{
  LNO preamble = sizeof(LNO);

  if (vLen == 0){
    preamble *= (1 + count);
  }

  LNO nbytes = preamble + (preamble % sizeof(T));  // T alignment

  LNO offset = nbytes / sizeof(T);

  LNO *info = reinterpret_cast<LNO *>(buf);
  T* elements = reinterpret_cast<T *>(buf) + offset;

  *info++ = count;

  if (vLen == 0){
    for (LNO i=0; i < count; i++){
      int nelements = v[i].size();
      *info++ = offset;
  
      for (LNO j=0; j < nelements; j++)
        *elements++ = v[i][j];
  
      offset += nelements;
    }
  }
  else{
    for (LNO i=0; i < count; i++){
      for (LNO j=0; j < vLen; j++){
        *elements++ = v[i][j];
      }
    }
  }
}

template <typename T, typename LNO>
LNO fromIndirectBytesToObjectCount(const LNO bytes, char buf[]) 
{
  LNO *count = reinterpret_cast<LNO *>(buf);
  return count[0];
}

template <typename T, typename LNO>
  void deserialize(const LNO bytes, const char buf[], 
     const LNO count, std::vector<T> v[], LNO vLen=0)
{
  LNO preamble = sizeof(LNO);

  if (vLen == 0){
    preamble *= (1 + count);
  }

  LNO nbytes = preamble + (preamble % sizeof(T));  // T alignment
  LNO offset = nbytes / sizeof(T);

  const T* elements = reinterpret_cast<const T *>(buf) + offset;

  if (vLen > 0){
    for (LNO i=0; i < count; i++){
      v[i].resize(vLen);
      v[i].clear();
      for (LNO j=0; j < vLen; j++){
        v[i].push_back(*elements++);
      }
    }
  }
  else{
    const LNO *info = reinterpret_cast<const LNO *>(buf) + 1;
    LNO lastOffset = LNO(bytes/sizeof(T));
  
    for (LNO i=0; i < count; i++){
  
      LNO length = ((i == count-1) ? lastOffset : info[i+1]) - info[i];
  
      v[i].resize(length);
      v[i].clear();
      for (LNO j=0; j < length; j++){
        v[i].push_back(*elements++);
      }
    }
  }

}

/*! \brief AlltoAllv sends/receives a std::vector<T> to/from all processes.
 *
 *  \param  comm   The communicator for the process group involved.
 *  \param  env    The environment, required for error messages.
 *  \param  sendBuf  The data to be sent, in destination process rank order.
 *  \param  sendCount sendCount[p] is the count of vectors to be sent 
                       to process p.
 *  \param  recvBuf  On return, recvBuf has been allocated and contains
 *                   the vectors sent to this process by others.
 *  \param  recvCount On return, The number of vectors received from process p 
 *                     will be in recvCount[p].
 *  \param  vLen     If all vectors are the same length, set vLen to the
 *                  length of the vectors to make AlltoAllv more efficient.
 *
 * The vectors need not be the same length.
 *
 * This was written to better understand Teuchos communication, but it
 * may be useful as well. 
 */

template <typename T, typename LNO>
void AlltoAllv(const Comm<int>     &comm,
  const Environment &env,
  const ArrayView<const std::vector<T> > &sendBuf,
  const ArrayView<const LNO>             &sendCount,
  ArrayRCP<std::vector<T> >        &recvBuf,
  ArrayRCP<LNO>                    &recvCount,
  LNO            vLen=0)      // set if all vectors are the same length
{
  int nprocs = comm.getSize();
  size_t totalSendSize = 0;
  LNO offset = 0;
  using Teuchos::is_null;

  LNO *sendSize = new LNO [nprocs];
  env.globalMemoryAssertion(__FILE__, __LINE__, nprocs, sendSize,
    rcp(&comm, false));

  for (int p=0; p < nprocs; p++){
    if (sendCount[p] > 0){
      sendSize[p] = 
        fromObjectsToIndirectBytes<T, LNO>(sendCount[p], 
          sendBuf.getRawPtr() + offset, vLen);

      offset += sendCount[p];
      totalSendSize += sendSize[p];
    }
    else{
      sendSize[p] = 0;
    }
  }

  size_t bufSize = totalSendSize/sizeof(T);
  
  T *buf = NULL;
  if (bufSize)
    buf = new T [bufSize];
  
  env.globalMemoryAssertion(__FILE__, __LINE__, bufSize, !bufSize || buf,
    rcp(&comm, false));

  const std::vector<T> *vptr = sendBuf.getRawPtr();

  char *charBuf = reinterpret_cast<char *>(buf);

  for (int p=0; p < nprocs; p++){
    if (sendCount[p] > 0){
      serialize<T, LNO>(sendCount[p], vptr, sendSize[p], charBuf, vLen);
      vptr += sendCount[p];
      charBuf += sendSize[p];
      sendSize[p] /= sizeof(T);
    }
  }

  ArrayRCP<T> recvT;
  ArrayRCP<LNO> recvSize;
  ArrayView<const T> bufView(buf, bufSize);
  ArrayView<const LNO> sendSizeView(sendSize, nprocs);

  try{
    AlltoAllv<T, LNO>(comm, env, bufView, sendSizeView, recvT, recvSize);
  }
  Z2_FORWARD_EXCEPTIONS;

  delete [] sendSize;
  if (bufSize)
    delete [] buf;

  LNO *vectorCount = new LNO [nprocs];
  env.globalMemoryAssertion(__FILE__, __LINE__, nprocs, vectorCount,
    rcp(&comm, false));

  LNO totalCount = 0;

  charBuf = reinterpret_cast<char *>(recvT.get());

  for (int p=0; p < nprocs; p++){
    if (recvSize[p] > 0){
      LNO bytecount = recvSize[p] * sizeof(T);
      vectorCount[p] = 
        fromIndirectBytesToObjectCount<T, LNO>(bytecount, charBuf);

      charBuf += bytecount;
      totalCount += vectorCount[p];
    }
    else{
      vectorCount[p] = 0;
    }
  }

  std::vector<T> *inVectors = NULL;
  if (totalCount)
    inVectors = new std::vector<T> [totalCount];

  env.globalMemoryAssertion(__FILE__, __LINE__, nprocs, 
    !totalCount || inVectors, rcp(&comm, false));

  charBuf = reinterpret_cast<char *>(recvT.get());
  std::vector<T> *inv = inVectors;

  for (int p=0; p < nprocs; p++){
    if (recvSize[p] > 0){
      LNO bytecount = recvSize[p] * sizeof(T);
      deserialize<T, LNO>(bytecount, charBuf, vectorCount[p], inv, vLen);

      charBuf += bytecount;
      inv += vectorCount[p];
    }
  }

  recvBuf = Teuchos::arcp(inVectors, 0, totalCount);
  recvCount = Teuchos::arcp(vectorCount, 0, nprocs);
}

}                   // namespace Z2
#endif
