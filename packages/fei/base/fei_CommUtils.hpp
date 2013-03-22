
/*--------------------------------------------------------------------*/
/*    Copyright 2007 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_CommUtils_hpp_
#define _fei_CommUtils_hpp_

#include <fei_macros.hpp>
#include <fei_mpi.h>
#include <fei_mpiTraits.hpp>
#include <fei_chk_mpi.hpp>
#include <fei_iostream.hpp>
#include <fei_CommMap.hpp>
#include <fei_TemplateUtils.hpp>
#include <snl_fei_RaggedTable.hpp>

#include <vector>
#include <set>
#include <map>

#include <fei_ErrMacros.hpp>
#undef fei_file
#define fei_file "fei_CommUtils.hpp"

namespace fei {

/** Return the MPI rank of the local processor.
  If the macro FEI_SER is defined, returns 0; otherwise calls MPI_Comm_rank.
*/
int localProc(MPI_Comm comm);

/** Return the number of processors (number of MPI ranks).
  If the macro FEI_SER is defined, returns 1; otherwise calls MPI_Comm_size.
*/
int numProcs(MPI_Comm comm);

void Barrier(MPI_Comm comm);

/** Scenario: The local processor has a list of processors to which data
    will be sent, but doesn't know which processors data will be received
    from. This method produces that list of processors to be received from.
    This is a collective method.
*/
int mirrorProcs(MPI_Comm comm, std::vector<int>& toProcs, std::vector<int>& fromProcs);

typedef snl_fei::RaggedTable<std::map<int,std::set<int>*>,std::set<int> > comm_map;

/** Given a comm-pattern, create and initialize its mirror...
 */
int mirrorCommPattern(MPI_Comm comm, comm_map* inPattern, comm_map*& outPattern);

/** Perform a "global OR" to reduce the local-bool to a global-bool.
   i.e., if localBool is true on any processor, then on exit globalBool will
   be true on all processors.
*/
int Allreduce(MPI_Comm comm, bool localBool, bool& globalBool);

/** Given a list of processors to send to, a scalar to send to each,
    and a list of processors to receive from, perform the exchange.

    @param sendProcs Input. List of processors to send to.
    @param sendData Input. List of data, same length as 'sendProcs', to be
    sent. (One item to be sent to each send proc.)
       Note: sendData is treated as const, but not declared const because
       it is passed to the MPI_Send routine which takes (non-const) void*.
    @param recvProcs Input. List of processors to receive from.
       Note: if you don't know which procs will be recv'd from, see the
         'mirrorProcs' function above.
    @param recvData Output. On exit, contains one item received from each
    recv proc.
    @return error-code 0 if successful
*/
int exchangeIntData(MPI_Comm comm,
                    const std::vector<int>& sendProcs,
                    std::vector<int>& sendData,
                    const std::vector<int>& recvProcs,
                    std::vector<int>& recvData);
 
//----------------------------------------------------------------------------
/** Perform an MPI_Allreduce with op = MPI_MAX
	@param local Input.
	@param global Output.
	@return MPI error-code
*/
template<class T>
int GlobalMax(MPI_Comm comm, std::vector<T>& local, std::vector<T>& global)
{
#ifdef FEI_SER
  global = local;
#else

  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  try {
    global.resize(local.size());
  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what()<<FEI_ENDL;
    return(-1);
  }

  CHK_MPI( MPI_Allreduce(&(local[0]), &(global[0]),
			 local.size(), mpi_dtype, MPI_MAX, comm) );
#endif

  return(0);
}

//----------------------------------------------------------------------------
/** Perform an MPI_Allreduce with op = MPI_MAX
	@param local Input.
	@param global Output.
	@return MPI error-code
*/
template<class T>
int GlobalMax(MPI_Comm comm, T local, T& global)
{
#ifdef FEI_SER
  global = local;
#else
  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  CHK_MPI( MPI_Allreduce(&local, &global, 1, mpi_dtype, MPI_MAX, comm) );
#endif
  return(0);
}

//----------------------------------------------------------------------------
/** Perform an MPI_Allreduce with op = MPI_MIN
	@param local Input.
	@param global Output.
	@return MPI error-code
*/
template<class T>
int GlobalMin(MPI_Comm comm, std::vector<T>& local, std::vector<T>& global)
{
#ifdef FEI_SER
  global = local;
#else

  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  try {
    global.resize(local.size());
  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what()<<FEI_ENDL;
    return(-1);
  }

  CHK_MPI( MPI_Allreduce(&(local[0]), &(global[0]),
			 local.size(), mpi_dtype, MPI_MIN, comm) );
#endif

  return(0);
}

//----------------------------------------------------------------------------
/** Perform an MPI_Allreduce with op = MPI_MIN
	@param local Input.
	@param global Output.
	@return MPI error-code
*/
template<class T>
int GlobalMin(MPI_Comm comm, T local, T& global)
{
#ifdef FEI_SER
  global = local;
#else
  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  CHK_MPI( MPI_Allreduce(&local, &global, 1, mpi_dtype, MPI_MIN, comm) );
#endif
  return(0);
}

//----------------------------------------------------------------------------
/** Perform an MPI_Allreduce with op = MPI_SUM
	@param local Input.
	@param global Output.
	@return MPI error-code
*/
template<class T>
int GlobalSum(MPI_Comm comm, std::vector<T>& local, std::vector<T>& global)
{
#ifdef FEI_SER
  global = local;
#else
  global.resize(local.size());

  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  CHK_MPI( MPI_Allreduce(&(local[0]), &(global[0]),
                      local.size(), mpi_dtype, MPI_SUM, comm) );
#endif
  return(0);
}

//----------------------------------------------------------------------------
/** Single-scalar version of the GlobalSum function. */
template<class T>
int GlobalSum(MPI_Comm comm, T local, T& global)
{
#ifdef FEI_SER
  global = local;
#else
  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  CHK_MPI( MPI_Allreduce(&local, &global, 1, mpi_dtype, MPI_SUM, comm) );
#endif
  return(0);
}


/** Allgatherv function that takes std::vectors.
*/
template<class T>
int Allgatherv(MPI_Comm comm,
               std::vector<T>& sendbuf,
               std::vector<int>& recvLengths,
               std::vector<T>& recvbuf)
{
#ifdef FEI_SER
  //If we're in serial mode, just copy sendbuf to recvbuf and return.

  recvbuf = sendbuf;
  recvLengths.resize(1);
  recvLengths[0] = sendbuf.size();
#else
  int numProcs = 1;
  MPI_Comm_size(comm, &numProcs);

  try {

  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  std::vector<int> tmpInt(numProcs, 0);

  int len = sendbuf.size();
  int* tmpBuf = &tmpInt[0];

  recvLengths.resize(numProcs);
  int* recvLenPtr = &recvLengths[0];

  CHK_MPI( MPI_Allgather(&len, 1, MPI_INT, recvLenPtr, 1, MPI_INT, comm) );

  int displ = 0;
  for(int i=0; i<numProcs; i++) {
    tmpBuf[i] = displ;
    displ += recvLenPtr[i];
  }

  if (displ == 0) {
    recvbuf.resize(0);
    return(0);
  }

  recvbuf.resize(displ);

  T* sendbufPtr = sendbuf.size()>0 ? &sendbuf[0] : NULL;
  
  CHK_MPI( MPI_Allgatherv(sendbufPtr, len, mpi_dtype,
			&recvbuf[0], &recvLengths[0], tmpBuf,
			mpi_dtype, comm) );

  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what() << FEI_ENDL;
    return(-1);
  }
#endif

  return(0);
}

//------------------------------------------------------------------------
template<class T>
int Bcast(MPI_Comm comm, std::vector<T>& sendbuf, int sourceProc)
{
#ifndef FEI_SER
  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  CHK_MPI(MPI_Bcast(&sendbuf[0], sendbuf.size(), mpi_dtype,
                    sourceProc, comm) );
#endif
  return(0);
}

//------------------------------------------------------------------------
/** Perform communication send data from one CommMap and receive into another CommMap.
  \param comm Input. MPI communicator.
  \param sendCommMap Input. Maps destination processors to vectors of data to be sent
                     to them.
  \param recvCommMap Output. On exit, will map source processors to vectors of data
                     that was received from them.
 */
template<typename T>
int exchangeCommMapData(MPI_Comm comm,
                        const typename CommMap<T>::Type& sendCommMap,
                        typename CommMap<T>::Type& recvCommMap,
                        bool recvProcsKnownOnEntry = false,
                        bool recvLengthsKnownOnEntry = false)
{
  if (!recvProcsKnownOnEntry) {
    recvCommMap.clear();
  }

#ifndef FEI_SER
  int tag = 11120;
  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  std::vector<int> sendProcs;
  fei::copyKeysToVector(sendCommMap, sendProcs);
  std::vector<int> recvProcs;

  if (recvProcsKnownOnEntry) {
    fei::copyKeysToVector(recvCommMap, recvProcs);
  }
  else {
    mirrorProcs(comm, sendProcs, recvProcs);
    for(size_t i=0; i<recvProcs.size(); ++i) {
      addItemsToCommMap<T>(recvProcs[i], 0, NULL, recvCommMap);
    }
  }

  if (!recvLengthsKnownOnEntry) {
    std::vector<int> tmpIntData(sendProcs.size());
    std::vector<int> recvLengths(recvProcs.size());
    
    typename fei::CommMap<T>::Type::const_iterator
      s_iter = sendCommMap.begin(), s_end = sendCommMap.end();

    for(size_t i=0; s_iter != s_end; ++s_iter, ++i) {
      tmpIntData[i] = s_iter->second.size();
    }

    if ( exchangeIntData(comm, sendProcs, tmpIntData, recvProcs, recvLengths) != 0) {
      return(-1);
    }
    for(size_t i=0; i<recvProcs.size(); ++i) {
      std::vector<T>& rdata = recvCommMap[recvProcs[i]];
      rdata.resize(recvLengths[i]);
    }
  }

  //launch Irecv's for recv-data:
  std::vector<MPI_Request> mpiReqs;
  mpiReqs.resize(recvProcs.size());

  typename fei::CommMap<T>::Type::iterator
    r_iter = recvCommMap.begin(), r_end = recvCommMap.end();

  size_t req_offset = 0;
  for(; r_iter != r_end; ++r_iter) {
    int rproc = r_iter->first;
    std::vector<T>& recv_vec = r_iter->second;
    int len = recv_vec.size();
    T* recv_buf = len > 0 ? &recv_vec[0] : NULL;

    CHK_MPI( MPI_Irecv(recv_buf, len, mpi_dtype, rproc,
                       tag, comm, &mpiReqs[req_offset++]) );
  }

  //send the send-data:

  typename fei::CommMap<T>::Type::const_iterator
    s_iter = sendCommMap.begin(), s_end = sendCommMap.end();

  for(; s_iter != s_end; ++s_iter) {
    int sproc = s_iter->first;
    const std::vector<T>& send_vec = s_iter->second;
    int len = send_vec.size();
    T* send_buf = len>0 ? const_cast<T*>(&send_vec[0]) : NULL;

    CHK_MPI( MPI_Send(send_buf, len, mpi_dtype, sproc, tag, comm) );
  }

  //complete the Irecvs:
  for(size_t i=0; i<mpiReqs.size(); ++i) {
    int index;
    MPI_Status status;
    CHK_MPI( MPI_Waitany(mpiReqs.size(), &mpiReqs[0], &index, &status) );
  }

#endif
  return(0);
}


//------------------------------------------------------------------------
template<class T>
int exchangeData(MPI_Comm comm,
                 std::vector<int>& sendProcs,
                 std::vector<std::vector<T> >& sendData,
                 std::vector<int>& recvProcs,
                 bool recvDataLengthsKnownOnEntry,
                 std::vector<std::vector<T> >& recvData)
{
  if (sendProcs.size() == 0 && recvProcs.size() == 0) return(0);
  if (sendProcs.size() != sendData.size()) return(-1);
#ifndef FEI_SER
  std::vector<MPI_Request> mpiReqs;
  mpiReqs.resize(recvProcs.size());

  int tag = 11119;
  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  if (!recvDataLengthsKnownOnEntry) {
    std::vector<int> tmpIntData(sendData.size());
    std::vector<int> recvLengths(recvProcs.size());
    for(unsigned i=0; i<sendData.size(); ++i) {
      tmpIntData[i] = sendData[i].size();
    }

    if ( exchangeIntData(comm, sendProcs, tmpIntData, recvProcs, recvLengths) != 0) {
      return(-1);
    }
    for(unsigned i=0; i<recvProcs.size(); ++i) {
      recvData[i].resize(recvLengths[i]);
    }
  }

  //launch Irecv's for recvData:

  size_t numRecvProcs = recvProcs.size();
  int req_offset = 0;
  int localProc = fei::localProc(comm);
  for(size_t i=0; i<recvProcs.size(); ++i) {
    if (recvProcs[i] == localProc) {--numRecvProcs; continue; }

    int len = recvData[i].size();
    std::vector<T>& recv_vec = recvData[i];
    T* recv_buf = len > 0 ? &recv_vec[0] : NULL;

    CHK_MPI( MPI_Irecv(recv_buf, len, mpi_dtype, recvProcs[i],
                       tag, comm, &mpiReqs[req_offset++]) );
  }

  //send the sendData:

  for(size_t i=0; i<sendProcs.size(); ++i) {
    if (sendProcs[i] == localProc) continue;

    std::vector<T>& send_buf = sendData[i];
    CHK_MPI( MPI_Send(&send_buf[0], sendData[i].size(), mpi_dtype,
                      sendProcs[i], tag, comm) );
  }

  //complete the Irecvs:
  for(size_t i=0; i<numRecvProcs; ++i) {
    if (recvProcs[i] == localProc) continue;
    int index;
    MPI_Status status;
    CHK_MPI( MPI_Waitany(numRecvProcs, &mpiReqs[0], &index, &status) );
  }

#endif
  return(0);
}

//------------------------------------------------------------------------
template<class T>
int exchangeData(MPI_Comm comm,
                 std::vector<int>& sendProcs,
                 std::vector<std::vector<T>*>& sendData,
                 std::vector<int>& recvProcs,
                 bool recvLengthsKnownOnEntry,
                 std::vector<std::vector<T>*>& recvData)
{
  if (sendProcs.size() == 0 && recvProcs.size() == 0) return(0);
  if (sendProcs.size() != sendData.size()) return(-1);
#ifndef FEI_SER
  int tag = 11115;
  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();
  std::vector<MPI_Request> mpiReqs;

  try {
  mpiReqs.resize(recvProcs.size());

  if (!recvLengthsKnownOnEntry) {
    std::vector<int> tmpIntData;
    tmpIntData.resize(sendData.size());
    std::vector<int> recvLens(sendData.size());
    for(unsigned i=0; i<sendData.size(); ++i) {
      tmpIntData[i] = (int)sendData[i]->size();
    }

    if (exchangeIntData(comm, sendProcs, tmpIntData, recvProcs, recvLens) != 0) {
      return(-1);
    }

    for(unsigned i=0; i<recvLens.size(); ++i) {
      recvData[i]->resize(recvLens[i]);
    }
  }
  }
  catch(std::runtime_error& exc) {
    fei::console_out() << exc.what() << FEI_ENDL;
    return(-1);
  }

  //launch Irecv's for recvData:

  size_t numRecvProcs = recvProcs.size();
  int req_offset = 0;
  int localProc = fei::localProc(comm);
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    if (recvProcs[i] == localProc) {--numRecvProcs; continue;}

    size_t len = recvData[i]->size();
    std::vector<T>& rbuf = *recvData[i];

    CHK_MPI( MPI_Irecv(&rbuf[0], (int)len, mpi_dtype,
                       recvProcs[i], tag, comm, &mpiReqs[req_offset++]) );
  }

  //send the sendData:

  for(unsigned i=0; i<sendProcs.size(); ++i) {
    if (sendProcs[i] == localProc) continue;

    std::vector<T>& sbuf = *sendData[i];
    CHK_MPI( MPI_Send(&sbuf[0], (int)sbuf.size(), mpi_dtype,
                      sendProcs[i], tag, comm) );
  }

  //complete the Irecv's:
  for(unsigned i=0; i<numRecvProcs; ++i) {
    if (recvProcs[i] == localProc) continue;
    int index;
    MPI_Status status;
    CHK_MPI( MPI_Waitany((int)numRecvProcs, &mpiReqs[0], &index, &status) );
  }

#endif
  return(0);
}


//------------------------------------------------------------------------
/** MessageHandler is an interface representing a general sparse data
    exchange among processors.

    Scenario: An object knows which processors will be sending data to the
    local processor and which processors will be receiving data from the
    local processor. Furthermore, the object can produce outgoing data on
    demand, and deal with incoming data when it arrives.

    Given the above conditions (embodied by an implementation of this
    MessageHandler interface), then the function exchange(...) can
    perform the inter-processor communication required to complete the data
    exchange.

    Note: If the object knows only which processors will be sent to but not
    which processors will be received from, or vice-versa, then the method
    mirrorProcs(...) can be used to find the missing information.
*/
template<class T>
class MessageHandler {
public:
  virtual ~MessageHandler(){}

  /** Get the list of processors to which messages are to be sent in a
      sparse data exchange.
  */
  virtual std::vector<int>& getSendProcs() = 0;

  /** Get the list of processors from which messages are to be received in
	a sparse data exchange.
  */
  virtual std::vector<int>& getRecvProcs() = 0;

  /** Get the length of a message that is to be send to the specified
	destination processor.
  */
  virtual int getSendMessageLength(int destProc, int& messageLength) = 0;

  /** Prepare (pack) a message that is to be sent to the specified
	destination processor.
  */
  virtual int getSendMessage(int destProc, std::vector<T>& message) = 0;

  /** Process a message that has been received from the specified source
	processor.
  */
  virtual int processRecvMessage(int srcProc, std::vector<T>& message) = 0;
};//class MessageHandler


//------------------------------------------------------------------------
template<class T>
int exchange(MPI_Comm comm, MessageHandler<T>* msgHandler)
{
#ifdef FEI_SER
  (void)msgHandler;
#else
  int numProcs = fei::numProcs(comm);
  if (numProcs < 2) return(0);

  std::vector<int>& sendProcs = msgHandler->getSendProcs();
  int numSendProcs = sendProcs.size();
  std::vector<int>& recvProcs = msgHandler->getRecvProcs();
  int numRecvProcs = recvProcs.size();
  int i;

  if (numSendProcs < 1 && numRecvProcs < 1) {
    return(0);
  }

  std::vector<int> sendMsgLengths(numSendProcs);

  for(i=0; i<numSendProcs; ++i) {
    CHK_ERR( msgHandler->getSendMessageLength(sendProcs[i], sendMsgLengths[i]) );
  }

  std::vector<std::vector<T> > recvMsgs(numRecvProcs);

  std::vector<std::vector<T> > sendMsgs(numSendProcs);
  for(i=0; i<numSendProcs; ++i) {
    CHK_ERR( msgHandler->getSendMessage(sendProcs[i], sendMsgs[i]) );
  }

  CHK_ERR( exchangeData(comm, sendProcs, sendMsgs,
                        recvProcs, false, recvMsgs) );

  for(i=0; i<numRecvProcs; ++i) {
    std::vector<T>& recvMsg = recvMsgs[i];
    CHK_ERR( msgHandler->processRecvMessage(recvProcs[i], recvMsg ) );
  }
#endif

  return(0);
}


} //namespace fei

#endif // _fei_CommUtils_hpp_

