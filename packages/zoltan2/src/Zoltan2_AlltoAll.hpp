// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_ALLTOALL_HPP_
#define _ZOLTAN2_ALLTOALL_HPP_

/*! \file Zoltan2_AlltoAll.hpp
*/

#include <vector>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_RCP.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_Exceptions.hpp>
#include <Zoltan2_Memory.hpp>

//
// TODO: doxygen comments and error handling and timing
//

namespace Z2
{

/*! \brief AlltoAll sends/receives a fixed number of objects to/from all processes.
 *
 * The data type T of the objects must be a type for which 
 * Teuchos::SerializationTraits are defined.  This is most likely every 
 * fundamental data type plus std::pair<T1,T2>. It does not
 * include std::vector<T2>.
 */

template <typename T, typename LNO>
void AlltoAll(const Teuchos::Comm<int> &comm,
              Zoltan2::Environment &env,
              const Teuchos::ArrayView<T> &sendBuf,  // input
              LNO count,                             // input
              Teuchos::ArrayRCP<T> &recvBuf)         // output - allocated here
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  if (count == 0) return;   // count is the same on all procs

  Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest> > req(nprocs-1);

  // Create a T-aligned receive buffer.

  T *ptr = NULL;
  Z2_SYNC_MEMORY_ALLOC(comm, env, T, ptr, nprocs * count);

  Teuchos::ArrayRCP<T> inBuf(ptr, 0, nprocs * count, true);

  // Do self messages

  for (LNO i=0, offset = rank*count; i < count; i++, offset++){
    inBuf.get()[offset] = sendBuf.getRawPtr()[offset];
  }

  // Post receives

  Teuchos::RCP<Teuchos::CommRequest> r;

  for (int p=0; p < nprocs; p++){
    if (p != rank){
      Teuchos::ArrayRCP<T> recvBufPtr(inBuf.get() + p*count, 0, count, false);
      try{
        r  = Teuchos::ireceive<int, T>(comm, recvBufPtr, p);
      }
      catch (const std::exception &e){
        Z2_THROW_OUTSIDE_ERROR(env, e);
      }
    
      req.push_back(r);
    }
  }

  // Wait until all are posted

  Teuchos::barrier(comm);

  // Do ready sends.

  for (int p=0; p < nprocs; p++){
    if (p != rank)
      try {
        Teuchos::readySend<int, T>(comm, sendBuf.view(p*count, count), p);
      } 
      catch (const std::exception &e)
        Z2_THROW_OUTSIDE_ERROR(env, e);
  }

  if (req.size() > 0){
    try {
      Teuchos::waitAll<int>(comm, req);
    }
    catch (const std::exception &e)
      Z2_THROW_OUTSIDE_ERROR(env, e);
  }

  recvBuf = inBuf;
}

/*! \brief AlltoAllv sends/receives a variable number of objects to/from all processes.
 *
 * The data type T of the objects must be a type for which 
 * Teuchos::SerializationTraits are defined.  This is most likely every 
 * fundamental data type plus std::pair<T1,T2>. It does not
 * include std::vector<T2>.
 */

template <typename T, typename LNO>
void AlltoAllv(const Teuchos::Comm<int> &comm,
              Zoltan2::Environment &env,  
              const Teuchos::ArrayView<T> &sendBuf,      // input
              const Teuchos::ArrayView<LNO> &sendCount,  // input
              Teuchos::ArrayRCP<T> &recvBuf,      // output, allocated here
              Teuchos::ArrayRCP<LNO> &recvCount)  // output, allocated here
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  try{
    AlltoAll<LNO, LNO>(comm, env, sendCount, 1, recvCount);
  }
  catch (const std::exception &e)
    Z2_THROW_ZOLTAN2_ERROR(env, e);

  size_t totalIn=0, offsetIn=0, offsetOut=0;

  for (int i=0; i < nprocs; i++){
    totalIn += recvCount[i];
    if (i < rank){
      offsetIn += recvCount[i];
      offsetOut += sendCount[i];
    }
  }

  T *ptr = NULL;
  Z2_SYNC_MEMORY_ALLOC(comm, env, T, ptr, totalIn);

  Teuchos::ArrayRCP<T> inBuf(ptr, 0, totalIn, true);

  T *in = inBuf.get() + offsetIn;           // Copy self messages
  T *out = sendBuf.getRawPtr() + offsetOut;

  for (LNO i=0; i < recvCount[rank]; i++){
    in[i] = out[i];
  }

  // Post receives

  Teuchos::RCP<Teuchos::CommRequest> r;
  Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest> > req(nprocs-1);

  offsetIn = 0;

  for (int p=0; p < nprocs; p++){
    if (p != rank && recvCount[p] > 0){
      Teuchos::ArrayRCP<T> recvBufPtr(inBuf.get() + offsetIn, 0, recvCount[p], false);

      try{
        r  = Teuchos::ireceive<int, T>(comm, recvBufPtr, p);
      }
      catch (const std::exception &e){
        Z2_THROW_OUTSIDE_ERROR(env, e);
      }
    
      req.push_back(r);
    }
    offsetIn += recvCount[p];
  }

  // Wait till all are posted

  Teuchos::barrier(comm);

  // Do all ready sends

  offsetOut = 0;

  for (int p=0; p < nprocs; p++){
    if (p != rank && sendCount[p] > 0){
      try{
        Teuchos::readySend<int, T>(comm, sendBuf.view(offsetOut, sendCount[p]), p);
      }
      catch(const std::exception &e)
        Z2_THROW_OUTSIDE_ERROR(env, e);
    }
    offsetOut += sendCount[p];
  }

  if (req.size() > 0){
    try{
      Teuchos::waitAll<int>(comm, req);
    }
    catch(const std::exception &e)
      Z2_THROW_OUTSIDE_ERROR(env, e);
  }

  recvBuf = inBuf;
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
 * The vectors need not be the same length. The data type T must be a type 
 * for which Teuchos::SerializationTraits are defined.  
 */

template <typename T, typename LNO>
void AlltoAllv(const Teuchos::Comm<int>     &comm,
  Zoltan2::Environment &env,
  const Teuchos::ArrayView<std::vector<T> > &sendBuf,
  const Teuchos::ArrayView<LNO>             &sendCount,
  Teuchos::ArrayRCP<std::vector<T> >        &recvBuf,
  Teuchos::ArrayRCP<LNO>                    &recvCount,
  LNO            vLen=0)      // set if all vectors are the same length
{
  int nprocs = comm.getSize();
  size_t totalSendSize = 0;
  LNO offset = 0;

  LNO *sendSize = NULL;
  Z2_SYNC_MEMORY_ALLOC(comm, env, LNO, sendSize, nprocs);

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

  T *buf = NULL;
  Z2_SYNC_MEMORY_ALLOC(comm, env, T, buf, totalSendSize/sizeof(T));

  std::vector<T> *vptr = sendBuf.getRawPtr();

  char *charBuf = reinterpret_cast<char *>(buf);

  for (int p=0; p < nprocs; p++){
    if (sendCount[p] > 0){
      serialize<T, LNO>(sendCount[p], vptr, sendSize[p], charBuf, vLen);
      vptr += sendCount[p];
      charBuf += sendSize[p];
      sendSize[p] /= sizeof(T);
    }
  }

  Teuchos::ArrayView<T> sendTView(buf, totalSendSize/sizeof(T));
  Teuchos::ArrayView<LNO> sendSizeView(sendSize, nprocs);
  Teuchos::ArrayRCP<T> recvT;
  Teuchos::ArrayRCP<LNO> recvSize;

  try{
    AlltoAllv<T, LNO>(comm, env, sendTView, sendSizeView, recvT, recvSize);
  }
  catch (const std::exception &e)
    Z2_THROW_ZOLTAN2_ERROR(env, e);

  if (buf)
    delete [] buf;

  delete [] sendSize;

  LNO *vectorCount = NULL;
  Z2_SYNC_MEMORY_ALLOC(comm, env, LNO, vectorCount, nprocs);

  LNO totalCount = 0;

  buf = recvT.get();
  charBuf = reinterpret_cast<char *>(buf);

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
  Z2_SYNC_MEMORY_ALLOC(comm, env, std::vector<T>, inVectors, totalCount);

  buf = recvT.get();
  charBuf = reinterpret_cast<char *>(buf);
  vptr = inVectors;

  for (int p=0; p < nprocs; p++){
    if (recvSize[p] > 0){
      LNO bytecount = recvSize[p] * sizeof(T);
      deserialize<T, LNO>(bytecount, charBuf, vectorCount[p], vptr, vLen);

      charBuf += bytecount;
      vptr += vectorCount[p];
    }
  }

  recvBuf = Teuchos::ArrayRCP<std::vector<T> >(inVectors, 0, totalCount, true);
  recvCount = Teuchos::ArrayRCP<LNO>(vectorCount, 0, nprocs, true);
}

}                   // namespace Z2
#endif
