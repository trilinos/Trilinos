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
#include <Teuchos_ArrayRCP.hpp>

#ifdef SERIALIZATION_SUPPORTS_VECTORS
//
// There is no instantiation for std::vector<> in SerializationTraits.  
//   Because this represents a fundamental global ID type in the original Zoltan 
//   (an array of arbitrary size of unsigned ints) we should support it.   
//   Message passing of std::vector<> will only occur during preprocessing 
//   and only in the event that the caller's global ID type is std::vector<T>.
//   Callers will know that this is not efficient.
//
// Layout for std::vector<T> -
//      int numberOfVectors                   /* LNO may be too small, use int */
//      int offsetToStartOfVectorElements[N]  /* offset (in Ts) from start of buffer */
//      T *vectorElements                     /* start on a sizeof(T)-aligned boundary */
//
// The SerializationTraits class, to support types that can be different sizes, needs methods :
//
//   fromObjectsToIndirectBytes
//   fromIndirectBytesToObjectCount
//
// to use instead of:
//
//   fromCountToIndirectBytes
//   fromIndirectBytesToCount
//
// The fromObjects method has the objects, not just a count of the objects.
// The ToObjectCount method has the buffer of chars, not just the size of the buffer. 
// There may be a better generalization of DirectSerialization_Traits.
//

#include<Teuchos_SerializationTraits.hpp>

template<typename Ordinal, typename T>
class SerializationTraits<Ordinal,std::vector<T> >
{
  // static Ordinal fromCountToIndirectBytes(const Ordinal count) 
  static Ordinal fromObjectsToIndirectBytes(const Ordinal count, const std::vector<T> buffer[])
  {
    Ordinal nbytes = 0;

    int preamble = sizeof(int) * (1 + count);
    int extra = preamble % sizeof(T);        // T alignment
    int numT = 0;

    for (int i=0; i < count; i++)
      numT += buffer[i].size();

    return preamble + extra + numT * sizeof(T);
  }

  static void serialize(
    const Ordinal count, const std::vector<T> buffer[], const Ordinal bytes, char charBuffer[])
  {
    int preamble = sizeof(int) * (1 + count);
    int extra = preamble % sizeof(T); 

    int *info = reinterpret_cast<int *>(charBuffer);
    int offset = (preamble + extra) / sizeof(T);
    T* elements = reinterpret_cast<T *>(charBuffer) + offset;

    *info++ = count;

    for (int i=0; i < count; i++){
      int nelements = buffer[i].size();
      *info++ = offset;

      for (j=0; j < nelements; j++)
        *elements++ = buffer[i][j];

      offset += nelements;
    }
  }

  static Ordinal fromIndirectBytesToObjectCount(const Ordinal bytes, char charBuffer[]) 
  {
    int *buf = reinterpret_cast<int *>charBuffer;
    return buf[0];
  }

  static void deserialize(
    const Ordinal bytes, const char charBuffer[], const Ordinal count, std::vector<T> buffer[])
  {
    int preamble = sizeof(int) * (1 + count);
    int extra = preamble % sizeof(T); 
    int offset = (preamble + extra) / sizeof(T);

    int *info = reinterpret_cast<int *>(charBuffer);
    T* elements = reinterpret_cast<T *>(charBuffer) + offset;

    info++;    // skip count

    for (i=0; i < count; i++){
      int offset = info[i];

      if (i == count-1)   
        nextOffset = bytes/sizeof(T);
      else
        nextOffset = info[i+1];

      int length = nextOffset - offset;

      buffer[i].resize(length);
      for (int j=0; j < length; j++){
        buffer[i].push_back(*elements++);
      }
    }
  }
};

}

#endif

namespace Z2
{

/*! \function Z2::AlltoAll
 *
 *  Interprocess communication for non-conforming data types.
 *
 *  Zoltan2 uses Teuchos for interprocess communication, which
 *  requires Scalars and GlobalOrdinals that support traits
 *  defined in Teuchos.  In general, the application global IDs
 *  may not be of Scalar or GlobalOrdinal type.  If they are not
 *  then they are converted to conforming types by Zoltan2.  But 
 *  this requires interprocess communcation during Problem pre- 
 *  and post-processing.
 *
 *  For this event, we define an AlltoAll that takes an arbitrary
 *  data type.  Once global IDs are translated to Zoltan2 internal
 *  global numbers, communication is done exclusively with Teuchos.
 */

// TODO error checking

template <typename T, typename GNO, typename LNO>
void AlltoAll(const Teuchos::Comm<int> &comm,
              const Teuchos::ArrayRCP<T> &sendBuf,  // input
              LNO count,                      // input
              Teuchos::ArrayRCP<T> &recvBuf)  // output - allocated in AlltoAll
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  if (count == 0) return;

  Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest> > req(nprocs-1);

  Teuchos::ArrayRCP<T> inBuf = Teuchos::arcp<T>(nprocs * count);

  for (LNO i=0, offset = rank*count; i < count; i++, offset++){
    inBuf.get()[offset] = sendBuf.get()[offset];
  }

  for (int p=0; p < nprocs; p++){
    if (p != rank){

      Teuchos::ArrayRCP<T> recvBufPtr(inBuf.get() + p*count, 0, count, false);

      req.push_back(Teuchos::ireceive<int, T>(comm, recvBufPtr, p));
    }
  }

  Teuchos::barrier(comm);

  for (int p=0; p < nprocs; p++){
    if (p != rank)
      Teuchos::readySend<int, T>(comm, Teuchos::ArrayView<T>(sendBuf.get() + p*count, count), p);
  }

  if (req.size() > 0){
    Teuchos::waitAll<int>(comm, req);
  }

  recvBuf = inBuf;
}

/*! \function Z2::AlltoAllv
 *
 *  Interprocess communication for non-conforming data types.
 *
 *  Zoltan2 uses Teuchos for interprocess communication, which
 *  requires Scalars and GlobalOrdinals that support traits
 *  defined in Teuchos.  In general, the application global IDs
 *  may not be of Scalar or GlobalOrdinal type.  If they are not
 *  then they are converted to conforming types by Zoltan2.  But 
 *  this requires interprocess communcation during Problem pre- 
 *  and post-processing.
 *
 *  For this event, we define an AlltoAllv that takes an arbitrary
 *  data type.  Once global IDs are translated to Zoltan2 internal
 *  global numbers, communication is done exclusively with Teuchos.
 */

// TODO error checking

template <typename T, typename GNO, typename LNO>
void AlltoAllv(const Teuchos::Comm<int> &comm,
              const Teuchos::ArrayRCP<T> &sendBuf,      // input
              const Teuchos::ArrayRCP<GNO> &sendCount,  // input
              Teuchos::ArrayRCP<T> &recvBuf,      // output, allocated in AlltoAllv
              Teuchos::ArrayRCP<GNO> &recvCount)  // output, allocated in AlltoAllv
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  AlltoAll<GNO, GNO, LNO>(comm, sendCount, 1, recvCount);

  GNO totalIn=0, offsetIn=0, offsetOut=0;

  for (int i=0; i < nprocs; i++){
    totalIn += recvCount[i];
    if (i < rank){
      offsetIn += recvCount[i];
      offsetOut += sendCount[i];
    }
  }

  Teuchos::ArrayRCP<T> inBuf = Teuchos::arcp<T>(totalIn);

  T *in = inBuf.get() + offsetIn;
  T *out = sendBuf.get() + offsetOut;

  for (LNO i=0; i < recvCount[rank]; i++){
    in[i] = out[i];
  }

  Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest> > req(nprocs-1);

  offsetIn = 0;

  for (int p=0; p < nprocs; p++){
    if (p != rank && recvCount[p] > 0){
      Teuchos::ArrayRCP<T> recvBufPtr(inBuf.get() + offsetIn, 0, recvCount[p], false);
      req.push_back(Teuchos::ireceive<int, T>(comm, recvBufPtr, p));
    }
    offsetIn += recvCount[p];
  }

  Teuchos::barrier(comm);

  offsetOut = 0;

  for (int p=0; p < nprocs; p++){
    if (p != rank && sendCount[p] > 0){
      Teuchos::ArrayView<T> sendBufView(sendBuf.get() + offsetOut, sendCount[p]);
      Teuchos::readySend<int, T>(comm, sendBufView, p);
    }
    offsetOut += sendCount[p];
  }

  Teuchos::waitAll<int>(comm, req);

  recvBuf = inBuf;
}


}                   // namespace Z2
#endif
