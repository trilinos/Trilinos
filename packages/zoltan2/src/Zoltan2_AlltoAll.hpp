// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_ALLTOALL_HPP_
#define _ZOLTAN2_ALLTOALL_HPP_

/*! \file Zoltan2_Hash.hpp
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
void AlltoAll(Teuchos::Comm<int> &comm,
              Teuchos::ArrayRCP<T> &sendBuf,
              LNO count,
              Teuchos::ArrayRCP<T> &recvBuf) 
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  if (count == 0) return;

  T *inbuf = new T [count * nprocs];
  Teuchos::ArrayRCP<T> inbufPtr(inbuf, 0, size_type(count * nprocs));

  T *in = inbuf;
  T *out = *sendBuf * (count*rank);

  std::vector<Teuchos::RCP<Teuchos::CommRequest> > req(nprocs-1);

  for (int p=0, next=0; p < nprocs; p++){

    if (p != rank){
      req[next++] = Teuchos::ireceive(comm, in, p);
    }
    else{
      for (LNO i=0; i < count; i++)
        in[i] = out[i];
    }

    in += count;
  }

  Teuchos::barrier(comm);

  for (int p=0, out=*sendBuf; p < nprocs; p++){

    if (p != rank)
      Teuchos::readySend(comm, out, p);

    out += count;
  }

  Teuchos::ArrayView<Teuchos::RCP<Teuchos::CommRequest> > avReq(req);

  Teuchos::waitAll(comm, avReq);

  recvBuf = inbufPtr;
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
void AlltoAllv(Teuchos::Comm<int> &comm,
              Teuchos::ArrayRCP<T> &sendBuf,
              Teuchos::ArrayRCP<GNO> &sendCount,
              Teuchos::ArrayRCP<T> &recvBuf, 
              Teuchos::ArrayRCP<GNO> &recvCount)
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  std::vector<Teuchos::RCP<Teuchos::CommRequest> > req(nprocs-1);

  AlltoAll(comm, sendCount, 1, recvCount);

  GNO totalIn=0, offset=0;

  for (int i=0; i < nprocs; i++){
    totalIn += recvCount[i];
    if (i < rank)
      offset += recvCount[i];
  }

  T *inbuf = new T [totalIn];
  Teuchos::ArrayRCP<T> inbufPtr(inbuf, 0, size_type(totalIn));

  T *in = inbuf;
  T *out = *sendBuf + offset;

  for (int p=0; p < nprocs; p++){

    if (p != rank && recvCount[p] > 0){
      req.push_back(Teuchos::ireceive(comm, in, p));
    }
    else if (p == rank && recvCount[p] > 0){
      for (LNO i=0; i < recvCount[rank]; i++){
        in[i] = out[i];
      }
    }

    in += recvCount[p];
  }

  Teuchos::barrier(comm);

  for (int p=0, out=*sendBuf; p < nprocs; p++){

    if (p != rank && sendCount[p] > 0)
      Teuchos::readySend(comm, out, p);

    out += sendCount[p];
  }

  Teuchos::ArrayView<Teuchos::RCP<Teuchos::CommRequest> > avReq(req);
  Teuchos::waitAll(comm, avReq);
  recvBuf = inbufPtr;
}


}                   // namespace Z2
#endif
