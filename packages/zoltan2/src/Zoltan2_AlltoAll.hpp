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
              const Teuchos::ArrayRCP<T> &sendBuf,  // input
              LNO count,                      // input
              Teuchos::ArrayRCP<T> &recvBuf)  // output - allocated in AlltoAll
{
  int nprocs = comm.getSize();
  int rank = comm.getRank();

  if (count == 0) return;

  T *inbuf = new T [count * nprocs];
  T *in = inbuf;
  T *out = sendBuf.get() + (count*rank);

  // using ArrayView of this (as waitAll documentation tells you to) causes waitAll runtime error
  // std::vector<Teuchos::RCP<Teuchos::CommRequest> > req(nprocs-1);

  Teuchos::Array<Teuchos::RCP<Teuchos::CommRequest> > req(nprocs-1);

  for (int p=0; p < nprocs; p++){

    if (p != rank){
      req.push_back(Teuchos::ireceive(comm, Teuchos::ArrayRCP<T>(in, 0, count, true), p));
    }
    else{
      for (LNO i=0; i < count; i++)
        in[i] = out[i];
    }

    in += count;
  }

  Teuchos::barrier(comm);

  out = sendBuf.get();

  for (int p=0; p < nprocs; p++){

    if (p != rank){
      Teuchos::ArrayView<T> view(out, count);
      Teuchos::readySend<int,T>(comm, view, p);
      //Teuchos::readySend(comm, view, p);             This doesn't compile
    }

    out += count;
  }

  if (req.size() > 0){

    // Using ArrayView of requests causes waitAll runtime error
    // Teuchos::ArrayView<Teuchos::RCP<Teuchos::CommRequest> > avReq = Teuchos::arrayViewFromVector(req);  
    // Teuchos::waitAll(comm, avReq);

    Teuchos::waitAll(comm, req);
  }

  recvBuf = Teuchos::ArrayRCP<T>(inbuf, 0, count * nprocs, true);
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

  std::vector<Teuchos::RCP<Teuchos::CommRequest> > req(nprocs-1);

  AlltoAll(comm, sendCount, 1, recvCount);

  GNO totalIn=0, offset=0;

  for (int i=0; i < nprocs; i++){
    totalIn += recvCount[i];
    if (i < rank)
      offset += recvCount[i];
  }

  T *inbuf = new T [totalIn];

  T *in = inbuf;
  T *out = sendBuf.get() + offset;

  for (int p=0; p < nprocs; p++){

    if (p != rank && recvCount[p] > 0){
      req.push_back(Teuchos::ireceive(comm, Teuchos::ArrayRCP<T>(in, 0, recvCount[p], true), p));
    }
    else if (p == rank && recvCount[p] > 0){
      for (LNO i=0; i < recvCount[rank]; i++){
        in[i] = out[i];
      }
    }

    in += recvCount[p];
  }

  Teuchos::barrier(comm);

  for (int p=0, out=sendBuf.get(); p < nprocs; p++){

    if (p != rank && sendCount[p] > 0)
      Teuchos::readySend(comm, Teuchos::ArrayView<T>(out, sendCount[p]), p);

    out += sendCount[p];
  }

  Teuchos::ArrayView<Teuchos::RCP<Teuchos::CommRequest> > avReq(req);
  Teuchos::waitAll(comm, avReq);

  recvBuf = Teuchos::ArrayRCP<T>(in, 0, totalIn, true);
}


}                   // namespace Z2
#endif
