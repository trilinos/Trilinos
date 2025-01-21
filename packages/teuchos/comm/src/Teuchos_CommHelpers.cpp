// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_CommHelpers.hpp"
#ifdef HAVE_TEUCHOS_MPI
#  include "Teuchos_Details_MpiCommRequest.hpp"
#  include "Teuchos_Details_MpiTypeTraits.hpp"
#endif // HAVE_TEUCHOS_MPI
#ifdef HAVE_TEUCHOSCORE_CXX11
#  include <memory>
#endif

namespace Teuchos {

#ifdef HAVE_TEUCHOS_MPI
namespace Details {

std::string getMpiErrorString (const int errCode) {
  // Space for storing the error string returned by MPI.
  // Leave room for null termination, since I don't know if MPI does this.
  char errString [MPI_MAX_ERROR_STRING+1];
  int errStringLen = MPI_MAX_ERROR_STRING; // output argument
  (void) MPI_Error_string (errCode, errString, &errStringLen);
  // errStringLen on output is the number of characters written.
  // I'm not sure (the MPI 3.0 Standard doesn't say) if this
  // includes the '\0', so I'll make sure.  We reserved space for
  // the extra '\0' if needed.
  if (errString[errStringLen-1] != '\0') {
    errString[errStringLen] = '\0';
  }
  return std::string (errString); // This copies the original string.
}

} // namespace Details
#endif // HAVE_TEUCHOS_MPI

namespace { // (anonymous)

/// \brief Generic implementation of reduceAll().
/// \tparam T The type of data on which to reduce.  The requirements
///   for this type are the same as for the template parameter T of
///   Teuchos::Details::MpiTypeTraits.
///
/// This generic implementation factors out common code among all full
/// specializations of reduceAll() in this file.
template<class T>
void
reduceAllImpl (const Comm<int>& comm,
               const EReductionType reductType,
               const int count,
               const T sendBuffer[],
               T globalReducts[])
{
#ifdef HAVE_TEUCHOS_MPI
  using Teuchos::Details::MpiTypeTraits;

  // mfh 17 Oct 2012: Even in an MPI build, Comm might be either a
  // SerialComm or an MpiComm.  If it's something else, we fall back
  // to the most general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    // Is it a SerialComm?
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) {
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
#ifdef HAVE_TEUCHOSCORE_CXX11
      std::unique_ptr<ValueTypeReductionOp<int, T> >
#else
      std::auto_ptr<ValueTypeReductionOp<int, T> >
#endif
          reductOp (createOp<int, T> (reductType));
      reduceAll (comm, *reductOp, count, sendBuffer, globalReducts);
    }
    else { // It's a SerialComm; there is only 1 process, so just copy.
      std::copy (sendBuffer, sendBuffer + count, globalReducts);
    }
  } else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Op rawMpiOp = ::Teuchos::Details::getMpiOpForEReductionType (reductType);
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    T t;
    MPI_Datatype rawMpiType = MpiTypeTraits<T>::getType (t);

    int err = MPI_SUCCESS;
    if (sendBuffer == globalReducts) {
      // NOTE (mfh 31 May 2017) This is only safe if the communicator
      // is NOT an intercomm.  The usual case is that communicators
      // are intracomms.
      err = MPI_Allreduce (MPI_IN_PLACE, globalReducts,
                           count, rawMpiType, rawMpiOp, rawMpiComm);
    }
    else {
      err = MPI_Allreduce (const_cast<T*> (sendBuffer), globalReducts,
                           count, rawMpiType, rawMpiOp, rawMpiComm);
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      err != MPI_SUCCESS,
      std::runtime_error,
      "MPI_Allreduce failed with the following error: "
      << ::Teuchos::Details::getMpiErrorString (err));
  }
#else
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  std::copy (sendBuffer, sendBuffer + count, globalReducts);
#endif // HAVE_TEUCHOS_MPI
}


/// \brief Generic implementation of gather().
/// \tparam T The type of data on which to reduce.  The requirements
///   for this type are the same as for the template parameter T of
///   Teuchos::Details::MpiTypeTraits.
///
/// This generic implementation factors out common code among all full
/// specializations of gather() in this file.
template<class T>
void
gatherImpl (const T sendBuf[],
            const int sendCount,
            T recvBuf[],
            const int recvCount,
            const int root,
            const Comm<int>& comm)
{
#ifdef HAVE_TEUCHOS_MPI
  using Teuchos::Details::MpiTypeTraits;

  // mfh 17 Oct 2012: Even in an MPI build, Comm might be either a
  // SerialComm or an MpiComm.  If it's something else, we fall back
  // to the most general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    // Is it a SerialComm?
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) {
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      gather<int, T> (sendBuf, sendCount, recvBuf, recvCount, root, comm);
    }
    else { // It's a SerialComm; there is only 1 process, so just copy.
      std::copy (sendBuf, sendBuf + sendCount, recvBuf);
    }
  } else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    T t;
    MPI_Datatype rawMpiType = MpiTypeTraits<T>::getType (t);
    const int err = MPI_Gather (const_cast<T*> (sendBuf), sendCount, rawMpiType,
                                recvBuf, recvCount, rawMpiType,
                                root, rawMpiComm);
    TEUCHOS_TEST_FOR_EXCEPTION(
      err != MPI_SUCCESS,
      std::runtime_error,
      "MPI_Gather failed with the following error: "
      << ::Teuchos::Details::getMpiErrorString (err));
  }
#else
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  std::copy (sendBuf, sendBuf + sendCount, recvBuf);
#endif // HAVE_TEUCHOS_MPI
}


/// \brief Generic implementation of scatter().
/// \tparam T The type of data on which to scatter.  The requirements
///   for this type are the same as for the template parameter T of
///   Teuchos::Details::MpiTypeTraits.
///
/// This generic implementation factors out common code among all full
/// specializations of scatter() in this file.
template<class T>
void
scatterImpl (const T sendBuf[],
             const int sendCount,
             T recvBuf[],
             const int recvCount,
             const int root,
             const Comm<int>& comm)
{
#ifdef HAVE_TEUCHOS_MPI
  using Teuchos::Details::MpiTypeTraits;

  // mfh 17 Oct 2012: Even in an MPI build, Comm might be either a
  // SerialComm or an MpiComm.  If it's something else, we fall back
  // to the most general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    // Is it a SerialComm?
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) {
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      scatter<int, T> (sendBuf, sendCount, recvBuf, recvCount, root, comm);
    }
    else { // It's a SerialComm; there is only 1 process, so just copy.
      std::copy (sendBuf, sendBuf + sendCount, recvBuf);
    }
  } else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    T t;
    MPI_Datatype rawMpiType = MpiTypeTraits<T>::getType (t);
    const int err =
      MPI_Scatter (const_cast<T*> (sendBuf), sendCount, rawMpiType,
                   recvBuf, recvCount, rawMpiType,
                   root, rawMpiComm);
    TEUCHOS_TEST_FOR_EXCEPTION
      (err != MPI_SUCCESS, std::runtime_error,
      "MPI_Scatter failed with the following error: "
      << ::Teuchos::Details::getMpiErrorString (err));
  }
#else
  // We've built without MPI, so just assume it's a SerialComm and
  // copy the data.
  std::copy (sendBuf, sendBuf + sendCount, recvBuf);
#endif // HAVE_TEUCHOS_MPI
}


/// \brief Generic implementation of reduce().
/// \tparam T The type of data on which to reduce.  The requirements
///   for this type are the same as for the template parameter T of
///   Teuchos::Details::MpiTypeTraits.
///
/// This generic implementation factors out common code among all full
/// specializations of reduce() in this file.
template<class T>
void
reduceImpl (const T sendBuf[],
            T recvBuf[],
            const int count,
            const EReductionType reductType,
            const int root,
            const Comm<int>& comm)
{
#ifdef HAVE_TEUCHOS_MPI
  using Teuchos::Details::MpiTypeTraits;

  // mfh 17 Oct 2012: Even in an MPI build, Comm might be either a
  // SerialComm or an MpiComm.  If it's something else, we fall back
  // to the most general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    // Is it a SerialComm?
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) {
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      reduce<int, T> (sendBuf, recvBuf, count, reductType, root, comm);
    }
    else { // It's a SerialComm; there is only 1 process, so just copy.
      std::copy (sendBuf, sendBuf + count, recvBuf);
    }
  } else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Op rawMpiOp = ::Teuchos::Details::getMpiOpForEReductionType (reductType);
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    T t;
    MPI_Datatype rawMpiType = MpiTypeTraits<T>::getType (t);
    const int err = MPI_Reduce (const_cast<T*> (sendBuf), recvBuf, count,
                                rawMpiType, rawMpiOp, root, rawMpiComm);
    TEUCHOS_TEST_FOR_EXCEPTION
      (err != MPI_SUCCESS, std::runtime_error, "MPI_Reduce failed with the "
       "following error: " << ::Teuchos::Details::getMpiErrorString (err));
  }
#else
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  std::copy (sendBuf, sendBuf + count, recvBuf);
#endif // HAVE_TEUCHOS_MPI
}


/// \brief Generic implementation of gatherv().
/// \tparam T The type of data on which to gather.  The requirements
///   for this type are the same as for the template parameter T of
///   Teuchos::Details::MpiTypeTraits.
///
/// This generic implementation factors out common code among all full
/// specializations of gatherv() in this file.
template<class T>
void
gathervImpl (const T sendBuf[],
             const int sendCount,
             T recvBuf[],
             const int recvCounts[],
             const int displs[],
             const int root,
             const Comm<int>& comm)
{
#ifdef HAVE_TEUCHOS_MPI
  using Teuchos::Details::MpiTypeTraits;

  // mfh 17 Oct 2012: Even in an MPI build, Comm might be either a
  // SerialComm or an MpiComm.  If it's something else, we fall back
  // to the most general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    // Is it a SerialComm?
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) {
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      gatherv<int, T> (sendBuf, sendCount, recvBuf, recvCounts, displs, root, comm);
    }
    else { // It's a SerialComm; there is only 1 process, so just copy.
      TEUCHOS_TEST_FOR_EXCEPTION(
        recvCounts[0] > sendCount, std::invalid_argument,
        "Teuchos::gatherv: If the input communicator contains only one "
        "process, then you cannot receive more entries than you send.  "
        "You aim to receive " << recvCounts[0] << " entries, but to send "
        << sendCount << " entries.");
      // Serial communicator case: just copy.  recvCounts[0] is the
      // amount to receive, so it's the amount to copy.  Start writing
      // to recvbuf at the offset displs[0].
      std::copy (sendBuf, sendBuf + recvCounts[0], recvBuf + displs[0]);
    }
  } else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    T t;
    MPI_Datatype rawMpiType = MpiTypeTraits<T>::getType (t);
    const int err = MPI_Gatherv (const_cast<T*> (sendBuf),
                                 sendCount,
                                 rawMpiType,
                                 recvBuf,
                                 const_cast<int*> (recvCounts),
                                 const_cast<int*> (displs),
                                 rawMpiType,
                                 root,
                                 rawMpiComm);
    TEUCHOS_TEST_FOR_EXCEPTION(
      err != MPI_SUCCESS,
      std::runtime_error,
      "MPI_Gatherv failed with the following error: "
      << ::Teuchos::Details::getMpiErrorString (err));
  }
#else
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  TEUCHOS_TEST_FOR_EXCEPTION(
    recvCounts[0] > sendCount, std::invalid_argument,
    "Teuchos::gatherv: If the input communicator contains only one "
    "process, then you cannot receive more entries than you send.  "
    "You aim to receive " << recvCounts[0] << " entries, but to send "
    << sendCount << " entries.");
  // Serial communicator case: just copy.  recvCounts[0] is the
  // amount to receive, so it's the amount to copy.  Start writing
  // to recvbuf at the offset displs[0].
  std::copy (sendBuf, sendBuf + recvCounts[0], recvBuf + displs[0]);
#endif // HAVE_TEUCHOS_MPI
}

/// \brief Generic implementation of ireceive() for any Comm subclass.
/// \tparam Packet The type of data to receive.
///
/// ireceiveImpl() falls back to this function if the given Comm is
/// neither an MpiComm, nor a SerialComm.
template<typename Packet>
RCP<Teuchos::CommRequest<int> >
ireceiveGeneral(const Comm<int>& comm,
                const ArrayRCP<Packet> &recvBuffer,
                const int sourceRank)
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::ireceive<int, " << "," << TypeNameTraits<Packet>::name ()
    << "> ( value type )"
    );
  ValueTypeSerializationBuffer<int, Packet>
    charRecvBuffer (recvBuffer.size (), recvBuffer.getRawPtr ());
  RCP<CommRequest<int> > commRequest =
    comm.ireceive (charRecvBuffer.getCharBufferView (), sourceRank);
  set_extra_data (recvBuffer, "buffer", inOutArg (commRequest));
  return commRequest;
}

/// \brief Variant of ireceiveGeneral that takes a tag.
/// It also restores the correct order of arguments.
template<typename Packet>
RCP<Teuchos::CommRequest<int> >
ireceiveGeneral (const ArrayRCP<Packet> &recvBuffer,
                 const int sourceRank,
                 const int tag,
                 const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::ireceive<int, " << "," << TypeNameTraits<Packet>::name ()
    << "> ( value type )"
    );
  ValueTypeSerializationBuffer<int, Packet>
    charRecvBuffer (recvBuffer.size (), recvBuffer.getRawPtr ());
  RCP<CommRequest<int> > commRequest =
    comm.ireceive (charRecvBuffer.getCharBufferView (), sourceRank, tag);
  set_extra_data (recvBuffer, "buffer", inOutArg (commRequest));
  return commRequest;
}

/// \brief Generic implementation of ireceive() for MpiComm.
/// \tparam T The type of data to receive.  The requirements for this
///   type are the same as for the template parameter T of
///   Teuchos::Details::MpiTypeTraits.
///
/// This generic implementation factors out common code among all full
/// specializations of ireceive() in this file.
///
/// \warning If the given Comm is actually a SerialComm, this method
///   will throw std::logic_error.  This is because SerialComm does
///   not correctly implement the equivalent of MPI_Irecv of a process
///   to itself.
template<class T>
RCP<CommRequest<int> >
ireceiveImpl (const Comm<int>& comm,
              const ArrayRCP<T>& recvBuffer,
              const int sourceRank)
{
#ifdef HAVE_TEUCHOS_MPI
  using Teuchos::Details::MpiTypeTraits;

  // Even in an MPI build, Comm might be either a SerialComm or an
  // MpiComm.  If it's something else, we fall back to the most
  // general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    // Is it a SerialComm?
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) {
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      return ireceiveGeneral<T> (comm, recvBuffer, sourceRank);
    }
    else { // SerialComm doesn't implement ireceive anyway.
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        std::logic_error,
        "ireceiveImpl: Not implemented for a serial communicator.");
    }
  }
  else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Comm rawComm = * (mpiComm->getRawMpiComm ());
    T t;
    MPI_Datatype rawType = MpiTypeTraits<T>::getType (t);
    T* rawRecvBuf = recvBuffer.getRawPtr ();
    const int count = as<int> (recvBuffer.size ());
    const int tag = mpiComm->getTag ();
    MPI_Request rawRequest = MPI_REQUEST_NULL;
    const int err = MPI_Irecv (rawRecvBuf, count, rawType, sourceRank, tag,
                               rawComm, &rawRequest);
    TEUCHOS_TEST_FOR_EXCEPTION(
      err != MPI_SUCCESS, std::runtime_error,
      "MPI_Irecv failed with the following error: "
      << ::Teuchos::Details::getMpiErrorString (err));

    ArrayRCP<const char> buf =
      arcp_const_cast<const char> (arcp_reinterpret_cast<char> (recvBuffer));
    RCP<Details::MpiCommRequest> req (new Details::MpiCommRequest (rawRequest, buf));
    return rcp_implicit_cast<CommRequest<int> > (req);
  }
#else
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,
    std::logic_error,
    "ireceiveImpl: Not implemented for a serial communicator.");

  // NOTE (mfh 15 Sep 2014): Most compilers have figured out that the
  // return statement below is unreachable.  Some older compilers
  // might not realize this.  That's why the return statement was put
  // there, so that those compilers don't warn that this function
  // doesn't return a value.  If it's a choice between one warning and
  // another, I would prefer the choice that produces less code and
  // doesn't have unreachable code (which never gets tested).

  //return null; // Guard to avoid compiler warning about not returning a value.
#endif // HAVE_TEUCHOS_MPI
}

/// \brief Variant of ireceiveImpl that takes a tag.
/// It also restores the correct order of arguments.
template<class T>
RCP<CommRequest<int> >
ireceiveImpl (const ArrayRCP<T>& recvBuffer,
              const int sourceRank,
              const int tag,
              const Comm<int>& comm)
{
#ifdef HAVE_TEUCHOS_MPI
  using Teuchos::Details::MpiTypeTraits;

  // Even in an MPI build, Comm might be either a SerialComm or an
  // MpiComm.  If it's something else, we fall back to the most
  // general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    // Is it a SerialComm?
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) {
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      return ireceiveGeneral<T> (recvBuffer, sourceRank, tag, comm);
    }
    else { // SerialComm doesn't implement ireceive anyway.
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        std::logic_error,
        "ireceiveImpl: Not implemented for a serial communicator.");
    }
  }
  else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Comm rawComm = * (mpiComm->getRawMpiComm ());
    T t;
    MPI_Datatype rawType = MpiTypeTraits<T>::getType (t);
    T* rawRecvBuf = recvBuffer.getRawPtr ();
    const int count = as<int> (recvBuffer.size ());
    MPI_Request rawRequest = MPI_REQUEST_NULL;
    const int err = MPI_Irecv (rawRecvBuf, count, rawType, sourceRank, tag,
                               rawComm, &rawRequest);
    TEUCHOS_TEST_FOR_EXCEPTION(
      err != MPI_SUCCESS, std::runtime_error,
      "MPI_Irecv failed with the following error: "
      << ::Teuchos::Details::getMpiErrorString (err));

    ArrayRCP<const char> buf =
      arcp_const_cast<const char> (arcp_reinterpret_cast<char> (recvBuffer));
    RCP<Details::MpiCommRequest> req (new Details::MpiCommRequest (rawRequest, buf));
    return rcp_implicit_cast<CommRequest<int> > (req);
  }
#else
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,
    std::logic_error,
    "ireceiveImpl: Not implemented for a serial communicator.");

  return null; // Guard to avoid compiler warning about not returning a value.
#endif // HAVE_TEUCHOS_MPI
}

/// \brief Generic implementation of send() for any Comm subclass.
/// \tparam T The type of data to send.
///
/// sendImpl() falls back to this function if the given Comm is
/// neither an MpiComm, nor a SerialComm.
template<class T>
void
sendGeneral (const Comm<int>& comm,
             const int count,
             const T sendBuffer[],
             const int destRank)
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::send<int, " << TypeNameTraits<T>::name () << ">");
  ConstValueTypeSerializationBuffer<int,T> charSendBuffer (count, sendBuffer);
  comm.send (charSendBuffer.getBytes (),
             charSendBuffer.getCharBuffer (),
             destRank);
}

/// \brief Variant of sendGeneral that takes a tag.
/// It also restores the correct order of arguments.
template<class T>
void
sendGeneral (const T sendBuffer[],
             const int count,
             const int destRank,
             const int tag,
             const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::send<int, " << TypeNameTraits<T>::name () << ">");
  ConstValueTypeSerializationBuffer<int,T> charSendBuffer (count, sendBuffer);
  comm.send (charSendBuffer.getBytes (),
             charSendBuffer.getCharBuffer (),
             destRank, tag);
}

/// \brief Generic implementation of send() for MpiComm.
/// \tparam T The type of data to send.  The requirements for this
///   type are the same as for the template parameter T of
///   Teuchos::Details::MpiTypeTraits.
///
/// This generic implementation factors out common code among all full
/// specializations of send() in this file.
///
/// \warning If the given Comm is actually a SerialComm, this method
///   will throw std::logic_error.  This is because SerialComm does
///   not correctly implement the equivalent of MPI_Send of a process
///   to itself.
template<class T>
void
sendImpl (const Comm<int>& comm,
          const int count,
          const T sendBuffer[],
          const int destRank)
{
#ifdef HAVE_TEUCHOS_MPI
  using Teuchos::Details::MpiTypeTraits;

  // Even in an MPI build, Comm might be either a SerialComm or an
  // MpiComm.  If it's something else, we fall back to the most
  // general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    // Is it a SerialComm?
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) {
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      sendGeneral<T> (comm, count, sendBuffer, destRank);
    }
    else { // SerialComm doesn't implement send correctly anyway.
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        std::logic_error,
        "sendImpl: Not implemented for a serial communicator.");
    }
  }
  else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Comm rawComm = * (mpiComm->getRawMpiComm ());
    T t;
    MPI_Datatype rawType = MpiTypeTraits<T>::getType (t);
    T* rawBuf = const_cast<T*> (sendBuffer);
    const int tag = mpiComm->getTag ();
    const int err = MPI_Send (rawBuf, count, rawType, destRank, tag, rawComm);
    TEUCHOS_TEST_FOR_EXCEPTION(
      err != MPI_SUCCESS,
      std::runtime_error,
      "MPI_Send failed with the following error: "
      << ::Teuchos::Details::getMpiErrorString (err));
  }
#else
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,
    std::logic_error,
    "sendImpl: Not implemented for a serial communicator.");
#endif // HAVE_TEUCHOS_MPI
}

/// \brief Variant of sendImpl that takes a tag.
/// It also restores the correct order of arguments.
template<class T>
void
sendImpl (const T sendBuffer[],
          const int count,
          const int destRank,
          const int tag,
          const Comm<int>& comm)
{
#ifdef HAVE_TEUCHOS_MPI
  using Teuchos::Details::MpiTypeTraits;

  // Even in an MPI build, Comm might be either a SerialComm or an
  // MpiComm.  If it's something else, we fall back to the most
  // general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    // Is it a SerialComm?
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) {
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      sendGeneral<T> (sendBuffer, count, destRank, tag, comm);
    }
    else { // SerialComm doesn't implement send correctly anyway.
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,
        std::logic_error,
        "sendImpl: Not implemented for a serial communicator.");
    }
  }
  else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Comm rawComm = * (mpiComm->getRawMpiComm ());
    T t;
    MPI_Datatype rawType = MpiTypeTraits<T>::getType (t);
    T* rawBuf = const_cast<T*> (sendBuffer);
    const int err = MPI_Send (rawBuf, count, rawType, destRank, tag, rawComm);
    TEUCHOS_TEST_FOR_EXCEPTION(
      err != MPI_SUCCESS,
      std::runtime_error,
      "MPI_Send failed with the following error: "
      << ::Teuchos::Details::getMpiErrorString (err));
  }
#else
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,
    std::logic_error,
    "sendImpl: Not implemented for a serial communicator.");
#endif // HAVE_TEUCHOS_MPI
}

/// \brief Generic implementation of isend() for any Comm subclass.
/// \tparam T The type of data to send.
///
/// isendImpl() falls back to this function if the given Comm is
/// neither an MpiComm, nor a SerialComm.
template<class T>
RCP<CommRequest<int> >
isendGeneral (const Comm<int>& comm,
              const ArrayRCP<const T>& sendBuffer,
              const int destRank)
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::isend<int," << TypeNameTraits<T>::name () << ">");
  ConstValueTypeSerializationBuffer<int, T>
    charSendBuffer (sendBuffer.size (), sendBuffer.getRawPtr ());
  RCP<CommRequest<int> > commRequest =
    comm.isend (charSendBuffer.getCharBufferView (), destRank);
  set_extra_data (sendBuffer, "buffer", inOutArg (commRequest));
  return commRequest;
}

/// \brief Generic implementation of isend() (with tag) for any Comm subclass.
/// \tparam T The type of data to send.
///
/// The version of isendImpl() that takes a tag falls back to this
/// function if the given Comm is neither an MpiComm, nor a
/// SerialComm.
template<class T>
RCP<CommRequest<int> >
isendGeneral (const ArrayRCP<const T>& sendBuffer,
              const int destRank,
              const int tag,
              const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::isend<int," << TypeNameTraits<T>::name () << ">");
  ConstValueTypeSerializationBuffer<int, T>
    charSendBuffer (sendBuffer.size (), sendBuffer.getRawPtr ());
  RCP<CommRequest<int> > commRequest =
    comm.isend (charSendBuffer.getCharBufferView (), destRank, tag);
  set_extra_data (sendBuffer, "buffer", inOutArg (commRequest));
  return commRequest;
}

/// \brief Variant of isendImpl() that takes a tag.
/// It also restores the correct order of arguments.
template<class T>
RCP<Teuchos::CommRequest<int> >
isendImpl (const ArrayRCP<const T>& sendBuffer,
           const int destRank,
           const int tag,
           const Comm<int>& comm)
{
#ifdef HAVE_TEUCHOS_MPI
  using Teuchos::Details::MpiTypeTraits;

  // Even in an MPI build, Comm might be either a SerialComm or an
  // MpiComm.  If it's something else, we fall back to the most
  // general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    // Is it a SerialComm?
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) {
      // We don't know what kind of Comm we have, so fall back to the
      // most general implementation.
      return isendGeneral<T> (sendBuffer, destRank, tag, comm);
    }
    else { // SerialComm doesn't implement send correctly anyway.
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        "isendImpl: Not implemented for a serial communicator.");
    }
  }
  else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Comm rawComm = * (mpiComm->getRawMpiComm ());
    T t;
    MPI_Datatype rawType = MpiTypeTraits<T>::getType (t);
    // MPI promises not to modify the send buffer; the const_cast
    // merely ensures compatibilty with C89, which does not have a
    // "const" keyword.
    T* rawSendBuf = const_cast<T*> (sendBuffer.getRawPtr ());
    const int count = as<int> (sendBuffer.size ());
    MPI_Request rawRequest = MPI_REQUEST_NULL;
    const int err = MPI_Isend (rawSendBuf, count, rawType, destRank, tag,
                               rawComm, &rawRequest);
    TEUCHOS_TEST_FOR_EXCEPTION(
      err != MPI_SUCCESS,
      std::runtime_error,
      "MPI_Isend failed with the following error: "
      << ::Teuchos::Details::getMpiErrorString (err));

    ArrayRCP<const char> buf = arcp_reinterpret_cast<const char> (sendBuffer);
    RCP<Details::MpiCommRequest> req (new Details::MpiCommRequest (rawRequest, buf));
    return rcp_implicit_cast<CommRequest<int> > (req);
  }
#else
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,
    std::logic_error,
    "isendImpl: Not implemented for a serial communicator.");
#endif // HAVE_TEUCHOS_MPI
}

} // namespace (anonymous)


// mfh 18 Oct 2012: Note on full template specializations
//
// To make Windows builds happy, declarations of full template
// specializations (as found in Teuchos_CommHelpers.hpp) must use the
// TEUCHOSCOMM_LIB_DLL_EXPORT macro.  However, _definitions_ of the
// specializations (as found in this file) must _not_ use the macro.
// That's why we don't use that macro here.

// amb See note in .hpp file.
#if 0
#ifdef HAVE_TEUCHOS_COMPLEX
// Specialization for Ordinal=int and Packet=std::complex<double>.
template<>
void
reduceAll<int, std::complex<double> > (const Comm<int>& comm,
                                       const EReductionType reductType,
                                       const int count,
                                       const std::complex<double> sendBuffer[],
                                       std::complex<double> globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, std::complex<double> > (" << count << ", "
    << toString (reductType) << ")"
    );
  reduceAllImpl<std::complex<double> > (comm, reductType, count, sendBuffer, globalReducts);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, std::complex<double> > (const Comm<int>& comm,
                                      const ArrayRCP<std::complex<double> >& recvBuffer,
                                      const int sourceRank)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, std::complex<double> >");
  return ireceiveImpl<std::complex<double> > (comm, recvBuffer, sourceRank);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, std::complex<double> > (const ArrayRCP<std::complex<double> >& recvBuffer,
                                      const int sourceRank,
                                      const int tag,
                                      const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, std::complex<double> >");
  return ireceiveImpl<std::complex<double> > (recvBuffer, sourceRank, tag, comm);
}

template<>
void
send<int, std::complex<double> > (const Comm<int>& comm,
                                  const int count,
                                  const std::complex<double> sendBuffer[],
                                  const int destRank)
{
  sendImpl<std::complex<double> > (comm, count, sendBuffer, destRank);
}

template<>
void
send<int, std::complex<double> > (const std::complex<double> sendBuffer[],
                                  const int count,
                                  const int destRank,
                                  const int tag,
                                  const Comm<int>& comm)
{
  sendImpl<std::complex<double> > (sendBuffer, count, destRank, tag, comm);
}

template<>
RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const std::complex<double> >& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<int>& comm)
{
  return isendImpl<std::complex<double> > (sendBuffer, destRank, tag, comm);
}

// Specialization for Ordinal=int and Packet=std::complex<float>.
template<>
void
reduceAll<int, std::complex<float> > (const Comm<int>& comm,
                                      const EReductionType reductType,
                                      const int count,
                                      const std::complex<float> sendBuffer[],
                                      std::complex<float> globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, std::complex<float> > (" << count << ", "
    << toString (reductType) << ")"
    );
  reduceAllImpl<std::complex<float> > (comm, reductType, count, sendBuffer, globalReducts);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, std::complex<float> > (const Comm<int>& comm,
                                     const ArrayRCP<std::complex<float> >& recvBuffer,
                                     const int sourceRank)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, std::complex<float> >");
  return ireceiveImpl<std::complex<float> > (comm, recvBuffer, sourceRank);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, std::complex<float> > (const ArrayRCP<std::complex<float> >& recvBuffer,
                                     const int sourceRank,
                                     const int tag,
                                     const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, std::complex<float> >");
  return ireceiveImpl<std::complex<float> > (recvBuffer, sourceRank, tag, comm);
}

template<>
void
send<int, std::complex<float> > (const Comm<int>& comm,
                                 const int count,
                                 const std::complex<float> sendBuffer[],
                                 const int destRank)
{
  return sendImpl<std::complex<float> > (comm, count, sendBuffer, destRank);
}

template<>
void
send<int, std::complex<float> > (const std::complex<float> sendBuffer[],
                                 const int count,
                                 const int destRank,
                                 const int tag,
                                 const Comm<int>& comm)
{
  return sendImpl<std::complex<float> > (sendBuffer, count, destRank, tag, comm);
}

template<>
RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const std::complex<float> >& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<int>& comm)
{
  return isendImpl<std::complex<float> > (sendBuffer, destRank, tag, comm);
}
#endif // HAVE_TEUCHOS_COMPLEX
#endif // if 0


// Specialization for Ordinal=int and Packet=double.
template<>
void
reduceAll<int, double> (const Comm<int>& comm,
                        const EReductionType reductType,
                        const int count,
                        const double sendBuffer[],
                        double globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, double> (" << count << ", "
    << toString (reductType) << ")"
    );
  reduceAllImpl<double> (comm, reductType, count, sendBuffer, globalReducts);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, double> (const Comm<int>& comm,
                       const ArrayRCP<double>& recvBuffer,
                       const int sourceRank)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, double>");
  return ireceiveImpl<double> (comm, recvBuffer, sourceRank);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, double> (const ArrayRCP<double>& recvBuffer,
                       const int sourceRank,
                       const int tag,
                       const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, double>");
  return ireceiveImpl<double> (recvBuffer, sourceRank, tag, comm);
}

template<>
void
send<int, double> (const Comm<int>& comm,
                   const int count,
                   const double sendBuffer[],
                   const int destRank)
{
  return sendImpl<double> (comm, count, sendBuffer, destRank);
}

template<>
void
send<int, double> (const double sendBuffer[],
                   const int count,
                   const int destRank,
                   const int tag,
                   const Comm<int>& comm)
{
  return sendImpl<double> (sendBuffer, count, destRank, tag, comm);
}

template<>
RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const double>& sendBuffer,
        const int destRank,
        const int tag,
        const Comm<int>& comm)
{
  return isendImpl<double> (sendBuffer, destRank, tag, comm);
}

// Specialization for Ordinal=int and Packet=float.
template<>
void
reduceAll<int, float> (const Comm<int>& comm,
                       const EReductionType reductType,
                       const int count,
                       const float sendBuffer[],
                       float globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, float> (" << count << ", "
    << toString (reductType) << ")"
    );
  reduceAllImpl<float> (comm, reductType, count, sendBuffer, globalReducts);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, float> (const Comm<int>& comm,
                      const ArrayRCP<float>& recvBuffer,
                      const int sourceRank)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, float>");
  return ireceiveImpl<float> (comm, recvBuffer, sourceRank);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, float> (const ArrayRCP<float>& recvBuffer,
                      const int sourceRank,
                      const int tag,
                      const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, float>");
  return ireceiveImpl<float> (recvBuffer, sourceRank, tag, comm);
}

template<>
void
send<int, float> (const Comm<int>& comm,
                  const int count,
                  const float sendBuffer[],
                  const int destRank)
{
  return sendImpl<float> (comm, count, sendBuffer, destRank);
}

template<>
void
send<int, float> (const float sendBuffer[],
                  const int count,
                  const int destRank,
                  const int tag,
                  const Comm<int>& comm)
{
  return sendImpl<float> (sendBuffer, count, destRank, tag, comm);
}

template<>
RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const float>& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<int>& comm)
{
  return isendImpl<float> (sendBuffer, destRank, tag, comm);
}


// Specialization for Ordinal=int and Packet=long long.
template<>
void
gather<int, long long> (const long long sendBuf[],
                        const int sendCount,
                        long long recvBuf[],
                        const int recvCount,
                        const int root,
                        const Comm<int>& comm)
{
  gatherImpl<long long> (sendBuf, sendCount, recvBuf, recvCount, root, comm);
}

template<>
void
gatherv<int, long long> (const long long sendBuf[],
                         const int sendCount,
                         long long recvBuf[],
                         const int recvCounts[],
                         const int displs[],
                         const int root,
                         const Comm<int>& comm)
{
  gathervImpl<long long> (sendBuf, sendCount, recvBuf, recvCounts, displs, root, comm);
}

template<>
void
reduceAll<int, long long> (const Comm<int>& comm,
                           const EReductionType reductType,
                           const int count,
                           const long long sendBuffer[],
                           long long globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, long long> (" << count << ", "
    << toString (reductType) << ")"
    );
  reduceAllImpl<long long> (comm, reductType, count, sendBuffer, globalReducts);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, long long> (const Comm<int>& comm,
                          const ArrayRCP<long long>& recvBuffer,
                          const int sourceRank)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, long long>");
  return ireceiveImpl<long long> (comm, recvBuffer, sourceRank);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, long long> (const ArrayRCP<long long>& recvBuffer,
                          const int sourceRank,
                          const int tag,
                          const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, long long>");
  return ireceiveImpl<long long> (recvBuffer, sourceRank, tag, comm);
}

template<>
void
send<int, long long> (const Comm<int>& comm,
                      const int count,
                      const long long sendBuffer[],
                      const int destRank)
{
  return sendImpl<long long> (comm, count, sendBuffer, destRank);
}

template<>
void
send<int, long long> (const long long sendBuffer[],
                      const int count,
                      const int destRank,
                      const int tag,
                      const Comm<int>& comm)
{
  return sendImpl<long long> (sendBuffer, count, destRank, tag, comm);
}

template<>
RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const long long>& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<int>& comm)
{
  return isendImpl<long long> (sendBuffer, destRank, tag, comm);
}

// Specialization for Ordinal=int and Packet=unsigned long long.
template<>
void
gather<int, unsigned long long> (const unsigned long long sendBuf[],
                                 const int sendCount,
                                 unsigned long long recvBuf[],
                                 const int recvCount,
                                 const int root,
                                 const Comm<int>& comm)
{
  gatherImpl<unsigned long long> (sendBuf, sendCount, recvBuf, recvCount, root, comm);
}

template<>
void
gatherv<int, unsigned long long> (const unsigned long long sendBuf[],
                                  const int sendCount,
                                  unsigned long long recvBuf[],
                                  const int recvCounts[],
                                  const int displs[],
                                  const int root,
                                  const Comm<int>& comm)
{
  gathervImpl<unsigned long long> (sendBuf, sendCount, recvBuf, recvCounts, displs, root, comm);
}

template<>
void
reduceAll<int, unsigned long long> (const Comm<int>& comm,
                                    const EReductionType reductType,
                                    const int count,
                                    const unsigned long long sendBuffer[],
                                    unsigned long long globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, unsigned long long> (" << count << ", "
    << toString (reductType) << ")"
    );
  reduceAllImpl<unsigned long long> (comm, reductType, count, sendBuffer, globalReducts);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, unsigned long long> (const Comm<int>& comm,
                                   const ArrayRCP<unsigned long long>& recvBuffer,
                                   const int sourceRank)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, unsigned long long>");
  return ireceiveImpl<unsigned long long> (comm, recvBuffer, sourceRank);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, unsigned long long> (const ArrayRCP<unsigned long long>& recvBuffer,
                                   const int sourceRank,
                                   const int tag,
                                   const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, unsigned long long>");
  return ireceiveImpl<unsigned long long> (recvBuffer, sourceRank, tag, comm);
}

template<>
void
send<int, unsigned long long> (const Comm<int>& comm,
                               const int count,
                               const unsigned long long sendBuffer[],
                               const int destRank)
{
  return sendImpl<unsigned long long> (comm, count, sendBuffer, destRank);
}

template<>
void
send<int, unsigned long long> (const unsigned long long sendBuffer[],
                               const int count,
                               const int destRank,
                               const int tag,
                               const Comm<int>& comm)
{
  return sendImpl<unsigned long long> (sendBuffer, count, destRank, tag, comm);
}

template<>
RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const unsigned long long>& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<int>& comm)
{
  return isendImpl<unsigned long long> (sendBuffer, destRank, tag, comm);
}


// Specialization for Ordinal=int and Packet=long.
template<>
void
gather<int, long> (const long sendBuf[],
                   const int sendCount,
                   long recvBuf[],
                   const int recvCount,
                   const int root,
                   const Comm<int>& comm)
{
  gatherImpl<long> (sendBuf, sendCount, recvBuf, recvCount, root, comm);
}

template<>
void
gatherv<int, long> (const long sendBuf[],
                    const int sendCount,
                    long recvBuf[],
                    const int recvCounts[],
                    const int displs[],
                    const int root,
                    const Comm<int>& comm)
{
  gathervImpl<long> (sendBuf, sendCount, recvBuf, recvCounts, displs, root, comm);
}

template<>
void
reduceAll<int, long> (const Comm<int>& comm,
                      const EReductionType reductType,
                      const int count,
                      const long sendBuffer[],
                      long globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, long> (" << count << ", "
    << toString (reductType) << ")"
    );
  reduceAllImpl<long> (comm, reductType, count, sendBuffer, globalReducts);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, long> (const Comm<int>& comm,
                     const ArrayRCP<long>& recvBuffer,
                     const int sourceRank)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, long>");
  return ireceiveImpl<long> (comm, recvBuffer, sourceRank);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, long> (const ArrayRCP<long>& recvBuffer,
                     const int sourceRank,
                     const int tag,
                     const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, long>");
  return ireceiveImpl<long> (recvBuffer, sourceRank, tag, comm);
}

template<>
void
send<int, long> (const Comm<int>& comm,
                 const int count,
                 const long sendBuffer[],
                 const int destRank)
{
  return sendImpl<long> (comm, count, sendBuffer, destRank);
}

template<>
void
send<int, long> (const long sendBuffer[],
                 const int count,
                 const int destRank,
                 const int tag,
                 const Comm<int>& comm)
{
  return sendImpl<long> (sendBuffer, count, destRank, tag, comm);
}

template<>
RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const long>& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<int>& comm)
{
  return isendImpl<long> (sendBuffer, destRank, tag, comm);
}


// Specialization for Ordinal=int and Packet=unsigned long.
template<>
void
gather<int, unsigned long> (const unsigned long sendBuf[],
                            const int sendCount,
                            unsigned long recvBuf[],
                            const int recvCount,
                            const int root,
                            const Comm<int>& comm)
{
  gatherImpl<unsigned long> (sendBuf, sendCount, recvBuf, recvCount, root, comm);
}

template<>
void
gatherv<int, unsigned long> (const unsigned long sendBuf[],
                             const int sendCount,
                             unsigned long recvBuf[],
                             const int recvCounts[],
                             const int displs[],
                             const int root,
                             const Comm<int>& comm)
{
  gathervImpl<unsigned long> (sendBuf, sendCount, recvBuf, recvCounts, displs, root, comm);
}

template<>
void
reduceAll<int, unsigned long> (const Comm<int>& comm,
                               const EReductionType reductType,
                               const int count,
                               const unsigned long sendBuffer[],
                               unsigned long globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, unsigned long> (" << count << ", "
    << toString (reductType) << ")"
    );
  reduceAllImpl<unsigned long> (comm, reductType, count, sendBuffer, globalReducts);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, unsigned long> (const Comm<int>& comm,
                              const ArrayRCP<unsigned long>& recvBuffer,
                              const int sourceRank)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, unsigned long>");
  return ireceiveImpl<unsigned long> (comm, recvBuffer, sourceRank);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, unsigned long> (const ArrayRCP<unsigned long>& recvBuffer,
                              const int sourceRank,
                              const int tag,
                              const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, unsigned long>");
  return ireceiveImpl<unsigned long> (recvBuffer, sourceRank, tag, comm);
}

template<>
void
send<int, unsigned long> (const Comm<int>& comm,
                          const int count,
                          const unsigned long sendBuffer[],
                          const int destRank)
{
  return sendImpl<unsigned long> (comm, count, sendBuffer, destRank);
}

template<>
void
send<int, unsigned long> (const unsigned long sendBuffer[],
                 const int count,
                 const int destRank,
                 const int tag,
                 const Comm<int>& comm)
{
  return sendImpl<unsigned long> (sendBuffer, count, destRank, tag, comm);
}

template<>
RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const unsigned long>& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<int>& comm)
{
  return isendImpl<unsigned long> (sendBuffer, destRank, tag, comm);
}

// Specialization for Ordinal=int and Packet=int.
template<>
void
gather<int, int> (const int sendBuf[],
                  const int sendCount,
                  int recvBuf[],
                  const int recvCount,
                  const int root,
                  const Comm<int>& comm)
{
  gatherImpl<int> (sendBuf, sendCount, recvBuf, recvCount, root, comm);
}

template<>
void
gatherv<int, int> (const int sendBuf[],
                   const int sendCount,
                   int recvBuf[],
                   const int recvCounts[],
                   const int displs[],
                   const int root,
                   const Comm<int>& comm)
{
  gathervImpl<int> (sendBuf, sendCount, recvBuf, recvCounts, displs, root, comm);
}

template<>
void
scatter<int, int> (const int sendBuf[],
                   const int sendCount,
                   int recvBuf[],
                   const int recvCount,
                   const int root,
                   const Comm<int>& comm)
{
  scatterImpl<int> (sendBuf, sendCount, recvBuf, recvCount, root, comm);
}

template<>
void
reduce<int, int> (const int sendBuf[],
                  int recvBuf[],
                  const int count,
                  const EReductionType reductType,
                  const int root,
                  const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR
    ("Teuchos::reduce<int, int> (" << count << ", " << toString (reductType)
     << ")");
  reduceImpl<int> (sendBuf, recvBuf, count, reductType, root, comm);
}
template<>
void
reduce<int, long> (const long sendBuf[],
                   long recvBuf[],
                   const int count,
                   const EReductionType reductType,
                   const int root,
                   const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR
    ("Teuchos::reduce<int, int> (" << count << ", " << toString (reductType)
     << ")");
  reduceImpl<long> (sendBuf, recvBuf, count, reductType, root, comm);
}

template<>
void
reduce<int, unsigned long> (const unsigned long sendBuf[],
                            unsigned long recvBuf[],
                            const int count,
                            const EReductionType reductType,
                            const int root,
                            const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR
    ("Teuchos::reduce<int, int> (" << count << ", " << toString (reductType)
     << ")");
  reduceImpl<unsigned long> (sendBuf, recvBuf, count, reductType, root, comm);
}

template<>
void
reduce<int, unsigned long long > (const unsigned long long sendBuf[],
                                  unsigned long long recvBuf[],
                                  const int count,
                                  const EReductionType reductType,
                                  const int root,
                                  const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR
    ("Teuchos::reduce<int, int> (" << count << ", " << toString (reductType)
     << ")");
  reduceImpl<unsigned long long> (sendBuf, recvBuf, count, reductType, root, comm);
}

template<>
void
reduce<int, double> (const double sendBuf[],
                     double recvBuf[],
                     const int count,
                     const EReductionType reductType,
                     const int root,
                     const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR
    ("Teuchos::reduce<int, int> (" << count << ", " << toString (reductType)
     << ")");
  reduceImpl<double> (sendBuf, recvBuf, count, reductType, root, comm);
}
template<>
void
reduceAll<int, int> (const Comm<int>& comm,
                     const EReductionType reductType,
                     const int count,
                     const int sendBuffer[],
                     int globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, int> (" << count << ", "
    << toString (reductType) << ")"
    );
  reduceAllImpl<int> (comm, reductType, count, sendBuffer, globalReducts);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, int> (const Comm<int>& comm,
                    const ArrayRCP<int>& recvBuffer,
                    const int sourceRank)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, int>");
  return ireceiveImpl<int> (comm, recvBuffer, sourceRank);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, int> (const ArrayRCP<int>& recvBuffer,
                    const int sourceRank,
                    const int tag,
                    const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, int>");
  return ireceiveImpl<int> (recvBuffer, sourceRank, tag, comm);
}

template<>
void
send<int, int> (const Comm<int>& comm,
                const int count,
                const int sendBuffer[],
                const int destRank)
{
  return sendImpl<int> (comm, count, sendBuffer, destRank);
}

template<>
void
send<int, int> (const int sendBuffer[],
                const int count,
                const int destRank,
                const int tag,
                const Comm<int>& comm)
{
  return sendImpl<int> (sendBuffer, count, destRank, tag, comm);
}

template<>
RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const int>& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<int>& comm)
{
  return isendImpl<int> (sendBuffer, destRank, tag, comm);
}

// Specialization for Ordinal=int and Packet=unsigned int.
template<>
void
gather<int, unsigned int> (const unsigned int sendBuf[],
                            const int sendCount,
                            unsigned int recvBuf[],
                            const int recvCount,
                            const int root,
                            const Comm<int>& comm)
{
  gatherImpl<unsigned int> (sendBuf, sendCount, recvBuf, recvCount, root, comm);
}

template<>
void
gatherv<int, unsigned int> (const unsigned int sendBuf[],
                             const int sendCount,
                             unsigned int recvBuf[],
                             const int recvCounts[],
                             const int displs[],
                             const int root,
                             const Comm<int>& comm)
{
  gathervImpl<unsigned int> (sendBuf, sendCount, recvBuf, recvCounts, displs, root, comm);
}

template<>
void
reduceAll<int, unsigned int> (const Comm<int>& comm,
                              const EReductionType reductType,
                              const int count,
                              const unsigned int sendBuffer[],
                              unsigned int globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, unsigned int> (" << count << ", "
    << toString (reductType) << ")"
    );
  reduceAllImpl<unsigned int> (comm, reductType, count, sendBuffer, globalReducts);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, unsigned int> (const Comm<int>& comm,
                             const ArrayRCP<unsigned int>& recvBuffer,
                             const int sourceRank)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, unsigned int>");
  return ireceiveImpl<unsigned int> (comm, recvBuffer, sourceRank);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, unsigned int> (const ArrayRCP<unsigned int>& recvBuffer,
                             const int sourceRank,
                             const int tag,
                             const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, unsigned int>");
  return ireceiveImpl<unsigned int> (recvBuffer, sourceRank, tag, comm);
}

template<>
void
send<int, unsigned int> (const Comm<int>& comm,
                         const int count,
                         const unsigned int sendBuffer[],
                         const int destRank)
{
  return sendImpl<unsigned int> (comm, count, sendBuffer, destRank);
}

template<>
void
send<int, unsigned int> (const unsigned int sendBuffer[],
                         const int count,
                         const int destRank,
                         const int tag,
                         const Comm<int>& comm)
{
  return sendImpl<unsigned int> (sendBuffer, count, destRank, tag, comm);
}

template<>
RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const unsigned int>& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<int>& comm)
{
  return isendImpl<unsigned int> (sendBuffer, destRank, tag, comm);
}


// Specialization for Ordinal=int and Packet=short.
template<>
void
gather<int, short> (const short sendBuf[],
                    const int sendCount,
                    short recvBuf[],
                    const int recvCount,
                    const int root,
                    const Comm<int>& comm)
{
  gatherImpl<short> (sendBuf, sendCount, recvBuf, recvCount, root, comm);
}

template<>
void
gatherv<int, short> (const short sendBuf[],
                     const int sendCount,
                     short recvBuf[],
                     const int recvCounts[],
                     const int displs[],
                     const int root,
                     const Comm<int>& comm)
{
  gathervImpl<short> (sendBuf, sendCount, recvBuf, recvCounts, displs, root, comm);
}

template<>
void
reduceAll<int, short> (const Comm<int>& comm,
                       const EReductionType reductType,
                       const int count,
                       const short sendBuffer[],
                       short globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, short> (" << count << ", "
    << toString (reductType) << ")"
    );
  reduceAllImpl<short> (comm, reductType, count, sendBuffer, globalReducts);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, short> (const Comm<int>& comm,
                      const ArrayRCP<short>& recvBuffer,
                      const int sourceRank)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, short>");
  return ireceiveImpl<short> (comm, recvBuffer, sourceRank);
}

template<>
RCP<Teuchos::CommRequest<int> >
ireceive<int, short> (const ArrayRCP<short>& recvBuffer,
                      const int sourceRank,
                      const int tag,
                      const Comm<int>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR("ireceive<int, short>");
  return ireceiveImpl<short> (recvBuffer, sourceRank, tag, comm);
}

template<>
void
send<int, short> (const Comm<int>& comm,
                  const int count,
                  const short sendBuffer[],
                  const int destRank)
{
  return sendImpl<short> (comm, count, sendBuffer, destRank);
}

template<>
void
send<int, short> (const short sendBuffer[],
                  const int count,
                  const int destRank,
                  const int tag,
                  const Comm<int>& comm)
{
  return sendImpl<short> (sendBuffer, count, destRank, tag, comm);
}

template<>
RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const short>& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<int>& comm)
{
  return isendImpl<short> (sendBuffer, destRank, tag, comm);
}

// mfh 18 Oct 2012: The specialization for Packet=char seems to be
// causing problems such as the following:
//
// http://testing.sandia.gov/cdash/testDetails.php?test=9909246&build=747699
//
// I am disabling it for now.  This should revert back to the old
// behavior for Packet=char.  That should fix the Tpetra errors, since
// many Tpetra objects inherit from DistObject<char, ...>.
#if 0
// Specialization for Ordinal=int and Packet=char.
template<>
void
reduceAll<int, char> (const Comm<int>& comm,
                      const EReductionType reductType,
                      const int count,
                      const char sendBuffer[],
                      char globalReducts[])
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::reduceAll<int, char> (" << count << ", "
    << toString (reductType) << ")"
    );
  reduceAllImpl<char> (comm, reductType, count, sendBuffer, globalReducts);
}
#endif // 0

} // namespace Teuchos
