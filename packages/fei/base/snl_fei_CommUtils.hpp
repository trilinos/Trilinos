/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_CommUtils_hpp_
#define _snl_fei_CommUtils_hpp_

#include "fei_fwd.hpp"
#include "fei_mpi.h"
#include "snl_fei_Utils.hpp"
#include "fei_CommUtilsBase.hpp"
#include "fei_TemplateUtils.hpp"

namespace snl_fei {

  /** MessageHandler is an interface representing a general sparse data
      exchange among processors.

      Scenario: An object knows which processors will be sending data to the
      local processor and which processors will be receiving data from the
      local processor. Furthermore, the object can produce outgoing data on
      demand, and deal with incoming data when it arrives.

      Given the above conditions (embodied by an implementation of this
      MessageHandler interface), then the method CommUtils::exchange() can
      perform the inter-processor communication required to complete the data
      exchange.

      Note: If the object knows only which processors will be sent to but not
      which processors will be received from, or vice-versa, then the method
      CommUtils::mirrorProcs() can be used to find the missing information.
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
  };

  /** Container for basic inter-processor communication utilities. */
  template<class T>
  class CommUtils : public fei::CommUtilsBase {
  public:
    /** Constructor */
    CommUtils(MPI_Comm comm);

    /** Destructor */
    virtual ~CommUtils();

   /** Call this function when each processor has data in sendbuf that needs to
       be sent to every other processor.
       On exit, recvbuf is a packed list that contains the sendbufs from all
       processors. The lengths of each processor's sendbuf are first gathered
       to all procs using MPI_Allgather, then recvbuf is allocated, and finally
       MPI_Allgatherv is called to gather the data.
       Uses the MPI_Comm that this CommUtils instance was constructed with.
       All processors in that comm group must call this function before any
       can exit.
       @param sendbuf Input. Local data to be sent to all processors.
       @param recvLengths Output. Array, length numProcs. recvLengths[i] is the
       length of the buf recv'd from proc i.
       @param recvbuf Output. Array of data from all processors.
       @return MPI error-code
    */
    int Allgatherv(std::vector<T>& sendbuf,
                   std::vector<int>& recvLengths,
                   std::vector<T>& recvbuf);

    /** Perform an MPI_Bcast using the communicator that this object was
	constructed with.
	NOTE: sendbuf must be the same length on all processors, on entry.
    */
    int Bcast(std::vector<T>& sendbuf, int sourceProc);

    /** Perform an MPI_Allreduce with op = MPI_MAX
	@param local Input.
	@param global Output.
	@return MPI error-code
    */
    int GlobalMax(std::vector<T>& local, std::vector<T>& global);

    /** Perform an MPI_Allreduce with op = MPI_MAX
	@param local Input.
	@param global Output.
	@return MPI error-code
    */
    int GlobalMax(T local, T& global);

    /** Perform an MPI_Allreduce with op = MPI_SUM
	@param local Input.
	@param global Output.
	@return MPI error-code
    */
    int GlobalSum(std::vector<T>& local, std::vector<T>& global);

    /** Single-scalar version of the GlobalSum function. */
    int GlobalSum(T local, T& global);

    /** Perform a "global OR" to reduce the local-bool to a global-bool.
	i.e., if localBool is true on any processor, then on exit globalBool will
	be true on all processors.
     */
    int Allreduce(bool localBool, bool& globalBool);

    /** Given a list of processors to send to, and corresponding list of data
	to send to each send-proc, and a list of processors to receive from,
	perform the exchange.

	@param sendProcs Input. List of processors to send to.
	@param sendData Input. Lists of data to be sent, one list per send-proc.
	@param recvProcs Input. List of processors to receive from.
	@param recvLengths Output. List, of same length as recvProcs, containing
	the length of the message received from each recv-proc.
        @param recvLengthsKnownOnEntry Input. If true, skip calculation of
        lengths.
	@param recvData Output. Packed list containing all data received.
	@return error-code 0 if successful
    */
    int exchangeData(std::vector<int>& sendProcs,
		     std::vector<std::vector<T>*>& sendData,
		     std::vector<int>& recvProcs,
		     std::vector<int>& recvLengths,
		     bool recvLengthsKnownOnEntry,
		     std::vector<T>& recvData);

    /** Given a list of processors to send to, and corresponding list of data to
        send to each send-proc, and a list of processors to receive from,
        perform the exchange.

        @param sendProcs Input. List of processors to send to.
        @param sendData Input. Lists of data to be sent, one list per send-proc.
        @param recvProcs Input. List of processors to receive from.
        @param recvLengthsKnownOnEntry Input. If true, skip calculation of
        lengths.
        @param recvData Input/Output. array of arrays, one array per recv-proc.On
        entry, the "outer" array should be of length recvProcs.length(). On exit,
        each of the "inner" arrays will be sized to the correct message length
        and will contain the data received.
        @return error-code 0 if successful
    */
    int exchangeData(std::vector<int>& sendProcs,
                     std::vector<std::vector<T>*>& sendData,
                     std::vector<int>& recvProcs,
                     bool recvLengthsKnownOnEntry,
                     std::vector<std::vector<T>*>& recvData);

    /** Perform sparse data exchange. */
    int exchange(MessageHandler<T>* msgHandler);

  private:
    bool inErrorState_;
    std::vector<T> tmpData_;
  };

} //namespace snl_fei

#include <snl_fei_Utils.hpp>

#include <fei_chk_mpi.hpp>
#include <fei_mpiTraits.hpp>

#undef fei_file
#define fei_file "snl_fei_CommUtils.h"

#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
template<class T>
snl_fei::CommUtils<T>::CommUtils(MPI_Comm comm)
  : fei::CommUtilsBase(comm),
    inErrorState_(false),
    tmpData_(0)
{
  tmpData_.resize(numProcs_);
}

//----------------------------------------------------------------------------
template<class T>
snl_fei::CommUtils<T>::~CommUtils()
{
}

//----------------------------------------------------------------------------
template<class T>
int snl_fei::CommUtils<T>::Allgatherv(std::vector<T>& sendbuf,
                                      std::vector<int>& recvLengths,
                                      std::vector<T>& recvbuf)
{
  if (inErrorState_) return(-1);

  return( fei::Allgatherv<T>(comm_, sendbuf, recvLengths, recvbuf) );
}

//----------------------------------------------------------------------------
template<class T>
int snl_fei::CommUtils<T>::Allreduce(bool localBool, bool& globalBool)
{
#ifndef FEI_SER
  int localInt = localBool ? 1 : 0;
  int globalInt = 0;

  CHK_MPI( MPI_Allreduce(&localInt, &globalInt, 1, MPI_INT, MPI_MAX, comm_) );

  globalBool = globalInt==1 ? true : false;
#else
  globalBool = localBool;
#endif

  return(0);
}

//----------------------------------------------------------------------------
template<class T>
int snl_fei::CommUtils<T>::Bcast(std::vector<T>& sendbuf, int sourceProc)
{
#ifndef FEI_SER
  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  CHK_MPI(MPI_Bcast(&sendbuf[0], sendbuf.size(), mpi_dtype,
                    sourceProc, comm_) );
#endif
  return(0);
}

//----------------------------------------------------------------------------
template<class T>
int snl_fei::CommUtils<T>::GlobalMax(std::vector<T>& local, std::vector<T>& global)
{
#ifdef FEI_SER
  global = local;
#else

  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  try {
  global.resize(local.size());
  }
  catch(std::runtime_error& exc) {
    FEI_CERR << exc.what()<<FEI_ENDL;
    return(-1);
  }

  CHK_MPI( MPI_Allreduce(&(local[0]), &(global[0]),
			 local.size(), mpi_dtype, MPI_MAX, comm_) );
#endif

  return(0);
}

//----------------------------------------------------------------------------
template<class T>
int snl_fei::CommUtils<T>::GlobalMax(T local, T& global)
{
#ifdef FEI_SER
  global = local;
#else
  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  CHK_MPI( MPI_Allreduce(&local, &global, 1, mpi_dtype, MPI_MAX, comm_) );
#endif
  return(0);
}

//----------------------------------------------------------------------------
template<class T>
int snl_fei::CommUtils<T>::GlobalSum(std::vector<T>& local, std::vector<T>& global)
{
#ifdef FEI_SER
  global = local;
#else
  global.resize(local.size());

  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  CHK_MPI( MPI_Allreduce(&(local[0]), &(global[0]),
                      local.size(), mpi_dtype, MPI_SUM, comm_) );
#endif
  return(0);
}

//----------------------------------------------------------------------------
template<class T>
int snl_fei::CommUtils<T>::GlobalSum(T local, T& global)
{
#ifdef FEI_SER
  global = local;
#else
  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  CHK_MPI( MPI_Allreduce(&local, &global, 1, mpi_dtype, MPI_SUM, comm_) );
#endif
  return(0);
}

//----------------------------------------------------------------------------
template<class T>
int snl_fei::CommUtils<T>::exchangeData(std::vector<int>& sendProcs,
                                        std::vector<std::vector<T>*>& sendData,
                                        std::vector<int>& recvProcs,
                                        std::vector<int>& recvLengths,
                                        bool recvLengthsKnownOnEntry,
                                        std::vector<T>& recvData)
{
  commCore_->increment_tag();

  if (sendProcs.size() == 0 && recvProcs.size() == 0) return(0);
  if (sendProcs.size() != sendData.size()) return(-1);
#ifndef FEI_SER
  std::vector<MPI_Request>& mpiReqs = commCore_->mpiReqs_;
  mpiReqs.resize(recvProcs.size());
  recvLengths.resize(recvProcs.size());

  int tag = commCore_->get_tag();
  std::vector<int>& tmpIntData = commCore_->tmpIntData_;
  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();

  if (!recvLengthsKnownOnEntry) {
    tmpIntData.resize(sendData.size());
    for(unsigned i=0; i<sendData.size(); ++i) {
      tmpIntData[i] = sendData[i]->size();
    }

    CHK_ERR( exchangeIntData(sendProcs, tmpIntData, recvProcs, recvLengths) );
    int totalRecvLength = 0;
    for(unsigned i=0; i<recvLengths.size(); ++i) {
      totalRecvLength += recvLengths[i];
    }

    recvData.resize(totalRecvLength);
  }

  //launch Irecv's for recvData:

  unsigned numRecvProcs = recvProcs.size();
  int recv_offset = 0;
  int req_offset = 0;
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    if (recvProcs[i] == localProc_) {--numRecvProcs; continue; }

    int len = recvLengths[i];
    T* recvBuf = len>0 ? &(recvData[recv_offset]) : NULL;

    CHK_MPI( MPI_Irecv(recvBuf, len, mpi_dtype, recvProcs[i],
                       tag, comm_, &mpiReqs[req_offset++]) );

    recv_offset += len;
  }

  //send the sendData:

  for(unsigned i=0; i<sendProcs.size(); ++i) {
    if (sendProcs[i] == localProc_) continue;

    CHK_MPI( MPI_Send(&(*(sendData[i]))[0], sendData[i]->size(), mpi_dtype,
                      sendProcs[i], tag, comm_) );
  }

  //complete the Irecvs:
  for(unsigned i=0; i<numRecvProcs; ++i) {
    if (recvProcs[i] == localProc_) continue;
    int index;
    MPI_Status status;
    CHK_MPI( MPI_Waitany(numRecvProcs, &mpiReqs[0], &index, &status) );
  }

#endif
  return(0);
}

//----------------------------------------------------------------------------
template<class T>
int snl_fei::CommUtils<T>::exchangeData(std::vector<int>& sendProcs,
                                        std::vector<std::vector<T>*>& sendData,
                                        std::vector<int>& recvProcs,
                                        bool recvLengthsKnownOnEntry,
                                        std::vector<std::vector<T>*>& recvData)
{
  commCore_->increment_tag();

  if (sendProcs.size() == 0 && recvProcs.size() == 0) return(0);
  if (sendProcs.size() != sendData.size()) return(-1);
#ifndef FEI_SER
  int tag = commCore_->get_tag();
  MPI_Datatype mpi_dtype = fei::mpiTraits<T>::mpi_type();
  std::vector<MPI_Request>& mpiReqs = commCore_->mpiReqs_;

  try {
  mpiReqs.resize(recvProcs.size());

  if (!recvLengthsKnownOnEntry) {
    std::vector<int>& tmpIntData = commCore_->tmpIntData_;
    tmpIntData.resize(sendData.size());
    std::vector<int> recvLens(sendData.size());
    for(unsigned i=0; i<sendData.size(); ++i) {
      tmpIntData[i] = (int)sendData[i]->size();
    }

    CHK_ERR( exchangeIntData(sendProcs, tmpIntData, recvProcs, recvLens) );

    for(unsigned i=0; i<recvLens.size(); ++i) {
      recvData[i]->resize(recvLens[i]);
    }
  }
  }
  catch(std::runtime_error& exc) {
    FEI_CERR << exc.what() << FEI_ENDL;
    ERReturn(-1);
  }

  //launch Irecv's for recvData:

  size_t numRecvProcs = recvProcs.size();
  int req_offset = 0;
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    if (recvProcs[i] == localProc_) {--numRecvProcs; continue;}

    size_t len = recvData[i]->size();
    std::vector<T>& rbuf = *recvData[i];

    CHK_MPI( MPI_Irecv(&rbuf[0], (int)len, mpi_dtype,
                       recvProcs[i], tag, comm_, &mpiReqs[req_offset++]) );
  }

  //send the sendData:

  for(unsigned i=0; i<sendProcs.size(); ++i) {
    if (sendProcs[i] == localProc_) continue;

    std::vector<T>& sbuf = *sendData[i];
    CHK_MPI( MPI_Send(&sbuf[0], (int)sbuf.size(), mpi_dtype,
                      sendProcs[i], tag, comm_) );
  }

  //complete the Irecv's:
  for(unsigned i=0; i<numRecvProcs; ++i) {
    if (recvProcs[i] == localProc_) continue;
    int index;
    MPI_Status status;
    CHK_MPI( MPI_Waitany((int)numRecvProcs, &mpiReqs[0], &index, &status) );
  }

#endif
  return(0);
}

//----------------------------------------------------------------------------
template<class T>
int snl_fei::CommUtils<T>::exchange(MessageHandler<T>* msgHandler)
{
#ifdef FEI_SER
  (void)msgHandler;
#else
  if (numProcs_ < 2) return(0);

  std::vector<int>& sendProcs = msgHandler->getSendProcs();
  int numSendProcs = sendProcs.size();
  std::vector<int>& recvProcs = msgHandler->getRecvProcs();
  int numRecvProcs = recvProcs.size();
  int i;

  if (numSendProcs < 1 && numRecvProcs < 1) {
    return(0);
  }

  std::vector<int> sendMsgLengths(numSendProcs), recvMsgLengths(numRecvProcs);

  for(i=0; i<numSendProcs; ++i) {
    CHK_ERR( msgHandler->getSendMessageLength(sendProcs[i], sendMsgLengths[i]) );
  }

  std::vector<T> recvMsgs;

  std::vector<std::vector<T>* > sendMsgs(numSendProcs);
  for(i=0; i<numSendProcs; ++i) {
    sendMsgs[i] = new std::vector<T>;
    CHK_ERR( msgHandler->getSendMessage(sendProcs[i], *(sendMsgs[i])) );
  }

  CHK_ERR( exchangeData(sendProcs, sendMsgs,
                        recvProcs, recvMsgLengths, false, recvMsgs) );

  int offset = 0;
  T* msgDataPtr = &recvMsgs[0];
  for(i=0; i<numRecvProcs; ++i) {
    int msgLen = recvMsgLengths[i];
    T* mdPtr = &(msgDataPtr[offset]);
    std::vector<T> recvMsg(mdPtr, mdPtr+msgLen);
    CHK_ERR( msgHandler->processRecvMessage(recvProcs[i], recvMsg ) );
    offset += msgLen;
  }

  for(i=0; i<numSendProcs; ++i) {
    delete sendMsgs[i];
  }
#endif

  return(0);
}

#endif // _snl_fei_CommUtils_hpp_

