// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_CommHelpers.hpp"
#ifdef HAVE_MPI
#  include "Teuchos_Details_MpiCommRequest.hpp"
#endif // HAVE_MPI

namespace Teuchos {
namespace { // (anonymous)

#ifdef HAVE_MPI
//! Get the raw MPI_Op corresponding to the given reduction type enum value.
MPI_Op getMpiOpForEReductionType (const enum EReductionType reductionType) {
  switch (reductionType) {
  case REDUCE_SUM: return MPI_SUM;
  case REDUCE_MIN: return MPI_MIN;
  case REDUCE_MAX: return MPI_MAX;
  case REDUCE_AND: return MPI_LAND; // logical AND, not bitwise AND
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      "The given EReductionType value is invalid.");
  }
}

/// \brief MPI's error string corresponding to the given integer error code.
///
/// \param errCode [in] Integer error code returned by MPI functions.
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

/// \class MpiTypeTraits
/// \brief Traits class mapping from type T to its MPI_Datatype
/// \tparam T The type being sent or received.  T must be default
///   constructible.  It must also be either one of C++'s built-in
///   types (like \c int or \c double), or a struct or "struct-like"
///   type like <tt>std::complex<double</tt>, for which sizeof(T)
///   correctly conveys the amount of data to send or receive.
template<class T>
class MpiTypeTraits {
public:
  /// \brief The MPI_Datatype corresponding to the given T instance.
  ///
  /// For more generality, this method requires passing in a T
  /// instance.  The method may or may not ignore this instance,
  /// depending on the type T.
  static MPI_Datatype getType (const T&);
};

#ifdef TEUCHOS_HAVE_COMPLEX
template<>
class MpiTypeTraits<std::complex<double> > {
public:
  static MPI_Datatype getType (const T&) {
    return MPI_C_DOUBLE_COMPLEX;
  }
};

template<>
class MpiTypeTraits<std::complex<float> > {
public:
  static MPI_Datatype getType (const T&) {
    return MPI_C_FLOAT_COMPLEX;
  }
};
#endif // TEUCHOS_HAVE_COMPLEX

template<>
class MpiTypeTraits<double> {
public:
  static MPI_Datatype getType (const double&) {
    return MPI_DOUBLE;
  }
};

template<>
class MpiTypeTraits<float> {
public:
  static MPI_Datatype getType (const float&) {
    return MPI_FLOAT;
  }
};

#ifdef TEUCHOS_HAVE_LONG_LONG_INT
template<>
class MpiTypeTraits<long long> {
public:
  static MPI_Datatype getType (const long long&) {
    return MPI_LONG_LONG;
  }
};
#endif // TEUCHOS_HAVE_LONG_LONG_INT

template<>
class MpiTypeTraits<long> {
public:
  static MPI_Datatype getType (const long&) {
    return MPI_LONG;
  }
};

template<>
class MpiTypeTraits<int> {
public:
  static MPI_Datatype getType (const int&) {
    return MPI_INT;
  }
};

template<>
class MpiTypeTraits<short> {
public:
  static MPI_Datatype getType (const short&) {
    return MPI_SHORT;
  }
};
#endif // HAVE_MPI


/// \brief Generic implementation of reduceAll().
/// \tparam T The type of data on which to reduce.  The requirements
///   for this type are the same as for the template parameter T of
///   MpiTypeTraits.
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
#ifdef HAVE_MPI
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
      std::auto_ptr<ValueTypeReductionOp<int, T> > reductOp (createOp<int, T> (reductType));
      reduceAll (comm, *reductOp, count, sendBuffer, globalReducts);
    }
    else { // It's a SerialComm; there is only 1 process, so just copy.
      std::copy (sendBuffer, sendBuffer + count, globalReducts);
    }
  } else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Op rawMpiOp = getMpiOpForEReductionType (reductType);
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    T t;
    MPI_Datatype rawMpiType = MpiTypeTraits<T>::getType (t);
    const int err = MPI_Allreduce (const_cast<T*> (sendBuffer),
      globalReducts, count, rawMpiType, rawMpiOp, rawMpiComm);
    TEUCHOS_TEST_FOR_EXCEPTION(
      err != MPI_SUCCESS,
      std::runtime_error,
      "MPI_Allreduce failed with the following error: "
      << getMpiErrorString (err));
  }
#else
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  std::copy (sendBuffer, sendBuffer + count, globalReducts);
#endif // HAVE_MPI
}


/// \brief Generic implementation of gather().
/// \tparam T The type of data on which to reduce.  The requirements
///   for this type are the same as for the template parameter T of
///   MpiTypeTraits.
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
#ifdef HAVE_MPI
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
      << getMpiErrorString (err));
  }
#else
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  std::copy (sendBuf, sendBuf + sendCount, recvBuf);
#endif // HAVE_MPI
}


/// \brief Generic implementation of gatherv().
/// \tparam T The type of data on which to reduce.  The requirements
///   for this type are the same as for the template parameter T of
///   MpiTypeTraits.
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
#ifdef HAVE_MPI
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
      << getMpiErrorString (err));
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
#endif // HAVE_MPI
}


/// \brief Generic implementation of reduceAllAndScatter().
/// \tparam T The type of data on which to reduce.  The requirements
///   for this type are the same as for the template parameter T of
///   MpiTypeTraits.
///
/// This generic implementation factors out common code among all full
/// specializations of reduceAllAndScatter() in this file.
template<class T>
void
reduceAllAndScatterImpl (const Comm<int>& comm,
                         const EReductionType reductType,
                         const int sendCount,
                         const T sendBuffer[],
                         const int recvCounts[],
                         T myGlobalReducts[])
{
#ifdef HAVE_MPI
  // mfh 17 Oct 2012: Even in an MPI build, Comm might be either a
  // SerialComm or an MpiComm.  If it's something else, we fall back
  // to the most general implementation.
  const MpiComm<int>* mpiComm = dynamic_cast<const MpiComm<int>* > (&comm);
  if (mpiComm == NULL) {
    // Is it a SerialComm?
    const SerialComm<int>* serialComm = dynamic_cast<const SerialComm<int>* > (&comm);
    if (serialComm == NULL) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "reduceAllAndScatter: "
        "Not implemented for Comm not MpiComm or SerialComm.");
    }
    else { // It's a SerialComm; there is only 1 process, so just copy.
      //
      // FIXME (mfh 09 Apr 2013) Is this right?  This is what
      // SerialComm does.
      //
      std::copy (sendBuffer, sendBuffer + sendCount, myGlobalReducts);
    }
  } else { // It's an MpiComm.  Invoke MPI directly.
    MPI_Op rawMpiOp = getMpiOpForEReductionType (reductType);
    MPI_Comm rawMpiComm = * (mpiComm->getRawMpiComm ());
    T t;
    MPI_Datatype rawMpiType = MpiTypeTraits<T>::getType (t);
    const int err = MPI_Reduce_scatter (const_cast<T*> (sendBuffer),
                                        myGlobalReducts,
                                        const_cast<int*> (recvCounts),
                                        rawMpiType,
                                        rawMpiOp,
                                        rawMpiComm);
    TEUCHOS_TEST_FOR_EXCEPTION(
      err != MPI_SUCCESS,
      std::runtime_error,
      "MPI_Reduce_scatter failed with the following error: "
      << getMpiErrorString (err));
  }
#else
  // We've built without MPI, so just assume it's a SerialComm and copy the data.
  //
  // FIXME (mfh 09 Apr 2013) Is this right?  This is what SerialComm does.
  //
  std::copy (sendBuffer, sendBuffer + sendCount, myGlobalReducts);
#endif // HAVE_MPI
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
///   MpiTypeTraits.
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
#ifdef HAVE_MPI
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
      << getMpiErrorString (err));

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
#endif // HAVE_MPI
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
#ifdef HAVE_MPI
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
      << getMpiErrorString (err));

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
#endif // HAVE_MPI
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
///   MpiTypeTraits.
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
#ifdef HAVE_MPI
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
      << getMpiErrorString (err));
  }
#else
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,
    std::logic_error,
    "sendImpl: Not implemented for a serial communicator.");
#endif // HAVE_MPI
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
#ifdef HAVE_MPI
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
      << getMpiErrorString (err));
  }
#else
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,
    std::logic_error,
    "sendImpl: Not implemented for a serial communicator.");
#endif // HAVE_MPI
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
#ifdef HAVE_MPI
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
      << getMpiErrorString (err));

    ArrayRCP<const char> buf = arcp_reinterpret_cast<const char> (sendBuffer);
    RCP<Details::MpiCommRequest> req (new Details::MpiCommRequest (rawRequest, buf));
    return rcp_implicit_cast<CommRequest<int> > (req);
  }
#else
  TEUCHOS_TEST_FOR_EXCEPTION(
    true,
    std::logic_error,
    "isendImpl: Not implemented for a serial communicator.");
#endif // HAVE_MPI
}

} // namespace (anonymous)

// mfh 18 Oct 2012: Note on full template specializations
//
// To make Windows builds happy, declarations of full template
// specializations (as found in Teuchos_CommHelpers.hpp) must use the
// TEUCHOSCOMM_LIB_DLL_EXPORT macro.  However, _definitions_ of the
// specializations (as found in this file) must _not_ use the macro.
// That's why we don't use that macro here.

#ifdef TEUCHOS_HAVE_COMPLEX
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
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, std::complex<double> > (const Comm<int>& comm,
                                  const int count,
                                  const std::complex<double> sendBuffer[],
                                  const int destRank)
{
  sendImpl<std::complex<double> > (comm, count, sendBuffer, destRank);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, std::complex<double> > (const std::complex<double> sendBuffer[],
                                  const int count,
                                  const int destRank,
                                  const int tag,
                                  const Comm<int>& comm)
{
  sendImpl<std::complex<double> > (sendBuffer, count, destRank, tag, comm);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<Teuchos::CommRequest<int> >
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
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, std::complex<float> > (const Comm<int>& comm,
                                 const int count,
                                 const std::complex<float> sendBuffer[],
                                 const int destRank)
{
  return sendImpl<std::complex<float> > (comm, count, sendBuffer, destRank);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, std::complex<float> > (const std::complex<float> sendBuffer[],
                                 const int count,
                                 const int destRank,
                                 const int tag,
                                 const Comm<int>& comm)
{
  return sendImpl<std::complex<float> > (sendBuffer, count, destRank, tag, comm);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const std::complex<float> >& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<int>& comm)
{
  return isendImpl<std::complex<float> > (sendBuffer, destRank, tag, comm);
}
#endif // TEUCHOS_HAVE_COMPLEX


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
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, double> (const Comm<int>& comm,
                   const int count,
                   const double sendBuffer[],
                   const int destRank)
{
  return sendImpl<double> (comm, count, sendBuffer, destRank);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, double> (const double sendBuffer[],
                   const int count,
                   const int destRank,
                   const int tag,
                   const Comm<int>& comm)
{
  return sendImpl<double> (sendBuffer, count, destRank, tag, comm);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<Teuchos::CommRequest<int> >
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
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, float> (const Comm<int>& comm,
                  const int count,
                  const float sendBuffer[],
                  const int destRank)
{
  return sendImpl<float> (comm, count, sendBuffer, destRank);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, float> (const float sendBuffer[],
                  const int count,
                  const int destRank,
                  const int tag,
                  const Comm<int>& comm)
{
  return sendImpl<float> (sendBuffer, count, destRank, tag, comm);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const float>& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<int>& comm)
{
  return isendImpl<float> (sendBuffer, destRank, tag, comm);
}


#ifdef TEUCHOS_HAVE_LONG_LONG_INT
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
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, long long> (const Comm<int>& comm,
                      const int count,
                      const long long sendBuffer[],
                      const int destRank)
{
  return sendImpl<long long> (comm, count, sendBuffer, destRank);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, long long> (const long long sendBuffer[],
                      const int count,
                      const int destRank,
                      const int tag,
                      const Comm<int>& comm)
{
  return sendImpl<long long> (sendBuffer, count, destRank, tag, comm);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const long long>& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<int>& comm)
{
  return isendImpl<long long> (sendBuffer, destRank, tag, comm);
}
#endif // TEUCHOS_HAVE_LONG_LONG_INT


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
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, long> (const Comm<int>& comm,
                 const int count,
                 const long sendBuffer[],
                 const int destRank)
{
  return sendImpl<long> (comm, count, sendBuffer, destRank);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, long> (const long sendBuffer[],
                 const int count,
                 const int destRank,
                 const int tag,
                 const Comm<int>& comm)
{
  return sendImpl<long> (sendBuffer, count, destRank, tag, comm);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const long>& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<int>& comm)
{
  return isendImpl<long> (sendBuffer, destRank, tag, comm);
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
void
reduceAllAndScatter<int, int> (const Comm<int>& comm,
                               const EReductionType reductType,
                               const int sendCount,
                               const int sendBuffer[],
                               const int recvCounts[],
                               int myGlobalReducts[])
{
  reduceAllAndScatterImpl<int> (comm, reductType, sendCount, sendBuffer,
                                recvCounts, myGlobalReducts);
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
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, int> (const Comm<int>& comm,
                const int count,
                const int sendBuffer[],
                const int destRank)
{
  return sendImpl<int> (comm, count, sendBuffer, destRank);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, int> (const int sendBuffer[],
                const int count,
                const int destRank,
                const int tag,
                const Comm<int>& comm)
{
  return sendImpl<int> (sendBuffer, count, destRank, tag, comm);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<Teuchos::CommRequest<int> >
isend (const ArrayRCP<const int>& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<int>& comm)
{
  return isendImpl<int> (sendBuffer, destRank, tag, comm);
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
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, short> (const Comm<int>& comm,
                  const int count,
                  const short sendBuffer[],
                  const int destRank)
{
  return sendImpl<short> (comm, count, sendBuffer, destRank);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, short> (const short sendBuffer[],
                  const int count,
                  const int destRank,
                  const int tag,
                  const Comm<int>& comm)
{
  return sendImpl<short> (sendBuffer, count, destRank, tag, comm);
}

template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<Teuchos::CommRequest<int> >
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
