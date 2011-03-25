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
              const Teuchos::ArrayView<T> &sendBuf,  // input
              LNO count,                      // input
              Teuchos::ArrayRCP<T> &recvBuf)  // output - allocated here
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  if (count == 0) return;

  Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest> > req(nprocs-1);

  Teuchos::ArrayRCP<T> inBuf = Teuchos::arcp<T>(nprocs * count);

  for (LNO i=0, offset = rank*count; i < count; i++, offset++){
    inBuf.get()[offset] = sendBuf.getRawPtr()[offset];
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
      Teuchos::readySend<int, T>(comm, sendBuf.view(p*count, count), p);
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
              const Teuchos::ArrayView<T> &sendBuf,      // input
              const Teuchos::ArrayView<GNO> &sendCount,  // input
              Teuchos::ArrayRCP<T> &recvBuf,      // output, allocated here
              Teuchos::ArrayRCP<GNO> &recvCount)  // output, allocated here
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
  T *out = sendBuf.getRawPtr() + offsetOut;

  for (LNO i=0; i < recvCount[rank]; i++){
    in[i] = out[i];
  }

  Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest> > req(nprocs-1);

  offsetIn = 0;

  for (int p=0; p < nprocs; p++){
    if (p != rank && recvCount[p] > 0){
      Teuchos::ArrayRCP<T> recvBufPtr(inBuf.get() + offsetIn, 
        0, recvCount[p], false);
      req.push_back(Teuchos::ireceive<int, T>(comm, recvBufPtr, p));
    }
    offsetIn += recvCount[p];
  }

  Teuchos::barrier(comm);

  offsetOut = 0;

  for (int p=0; p < nprocs; p++){
    if (p != rank && sendCount[p] > 0){
      Teuchos::readySend<int, T>(comm, 
        sendBuf.view(offsetOut, sendCount[p]), p);
    }
    offsetOut += sendCount[p];
  }

  Teuchos::waitAll<int>(comm, req);

  recvBuf = inBuf;
}

// A version of AlltoAllv for sending std::vector<T>.
//
// There is no instantiation for std::vector<> in Teuchos::SerializationTraits. 
// These four methods are what it might look like.
//
// Layout for std::vector<T> of size N -
//      int numberOfVectors      
//      int offsetToStartOfVectorElements[N] 
//      T *vectorElements 
//

template <typename T>
  size_t fromObjectsToIndirectBytes(const LNO count, std::vector<T> const v[])
{
  size_t nbytes = 0;

  nbytes += sizeof(int) * (1 + count);   // preamble
  nbytes += preamble % sizeof(T);        // T alignment

  for (int i=0; i < count; i++)
    nbytes += v[i].size();

  return nbytes;
}

template <typename T>
  void serialize(
    const LNO count, const std::vector<T> v[], const size_t bytes, char buf[])
{
  int preamble = sizeof(int) * (1 + count);   // TODO make the ints LNOs
  int extra = preamble % sizeof(T); 

  int *info = reinterpret_cast<int *>(buf);
  int offset = (preamble + extra) / sizeof(T);
  T* elements = reinterpret_cast<T *>(buf) + offset;

  *info++ = count;

  for (int i=0; i < count; i++){
    int nelements = v[i].size();
    *info++ = offset;

    for (j=0; j < nelements; j++)
      *elements++ = v[i][j];

    offset += nelements;
  }
}

template <typename T>
int fromIndirectBytesToObjectCount(const size_t bytes, char buf[]) 
{
  int *count = reinterpret_cast<int *>buf;
  return buf[0];
}

template <typename T>
  void deserialize(const size_t bytes, const char buf[], 
     const LNO count, std::vector<T> v[])
{
  int preamble = sizeof(int) * (1 + count);
  int extra = preamble % sizeof(T); 
  int offset = (preamble + extra) / sizeof(T);

  int *info = reinterpret_cast<int *>(buf);
  T* elements = reinterpret_cast<T *>(buf) + offset;

  info++;    // skip count

  for (i=0; i < count; i++){
    int offset = info[i];

    if (i == count-1)   
      nextOffset = bytes/sizeof(T);
    else
      nextOffset = info[i+1];

    int length = nextOffset - offset;

    v[i].clear();
    v[i].resize(length);
    for (int j=0; j < length; j++){
      v[i].push_back(*elements++);
    }
  }
}

template <typename T, typename GNO, typename LNO>
void AlltoAllv(const Teuchos::Comm<int> &comm,
  const Teuchos::ArrayView<std::vector<T> > &sendBuf,
  const Teuchos::ArrayView<LNO>             &sendCount,
  Teuchos::ArrayRCP<std::vector<T> >        &recvBuf,
  Teuchos::ArrayRCP<LNO>                    &recvCount)
{
  int nprocs = comm.getSize();

  GNO *sendSize = new GNO [nprocs];
  GNO totalSendSize = 0;

  std::vector<T> *vptr = sendBuf.getRawPtr();

  for (int p=0; p < nprocs; p++){
    if (sendCount[p] > 0){
      sendSize[p] = 
        static_cast<GNO>(fromObjectsToIndirectBytes(sendCount[p], vptr));

      vptr += sendCount[p];
      totalSendSize += sendSize[p];
    }
    else{
      sendSize[p] = 0;
    }
  }

  char *sendChars = *buf = new char [totalSendSize];
  vptr = sendBuf.getRawPtr();

  for (int p=0; p < nprocs; p++){
    if (sendCount[p] > 0){
      serialize<T>(sendCount[p], vptr, sendSize[p], buf)
      vptr += sendCount[p];
      buf += sendSize[p];
    }
  }

  Teuchos::ArrayView<char> sendCharsView(sendChars, totalSendSize);
  Teuchos::ArrayView<GNO> sendSizeView(sendSize, nprocs);
  Teuchos::ArrayRCP<char> recvChars;
  Teuchos::ArrayRCP<GNO> recvSize;

  AlltoAllv<char, GNO, LNO>(comm, sendCharsView, sendSizeView,
                            recvChars, recvSize);

  delete [] sendChars;
  delete [] sendSize;

  LNO *inCount = new LNO [nprocs];
  buf = recvChars.get();
  LNO totalInCount = 0;

  for (int p=0; p < nprocs; p++){
    if (recvSize[p] > 0){
      inCount[p] = 
        fromIndirectBytesToObjectCount(recvSize[p], buf);

      buf += recvSize[p];
      totalInCount += inCount[p];
    }
    else{
      inCount[p] = 0;
    }
  }

  std::vector<T> *inVectors = new std::vector<T> [totalInCount];

  buf = recvChars.get();
  vptr = inVectors;

  for (int p=0; p < nprocs; p++){
    if (recvSize[p] > 0){
      deserialize(recvSize[p], buf, inCount[p], vptr);

      buf += recvSize[p];
      vptr += inCount[p];
    }
  }

  recvBuf = Teuchos::ArrayRCP<std::vector<T> >(inVectors, totalInCount);
  recvCount = Teuchos::ArrayRCP<LNO>(inCount, nprocs);
}

}                   // namespace Z2
#endif
