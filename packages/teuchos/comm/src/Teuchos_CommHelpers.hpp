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

#ifndef TEUCHOS_COMM_HELPERS_HPP
#define TEUCHOS_COMM_HELPERS_HPP

#include "Teuchos_Comm.hpp"
#include "Teuchos_CommUtilities.hpp"
#include "Teuchos_SerializationTraitsHelpers.hpp"
#include "Teuchos_ReductionOpHelpers.hpp"
#include "Teuchos_SerializerHelpers.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_as.hpp"

#ifdef HAVE_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#endif // HAVE_MPI
#include "Teuchos_DefaultSerialComm.hpp"


namespace Teuchos {


//
// Teuchos::Comm Helper Functions
//

/** \brief Enumeration for selecting from a set of pre-defined reduction
 * operations.
 *
 * \relates Comm
 */
enum EReductionType {
  REDUCE_SUM, ///< Sum
  REDUCE_MIN, ///< Min
  REDUCE_MAX, ///< Max
  REDUCE_AND ///< Logical AND
};

/** \brief Convert to std::string representation.
 *
 * \relates EReductionType
 */
inline
const char* toString( const EReductionType reductType )
{
  switch(reductType) {
    case REDUCE_SUM: return "REDUCE_SUM";
    case REDUCE_MIN: return "REDUCE_MIN";
    case REDUCE_MAX: return "REDUCE_MAX";
    case REDUCE_AND: return "REDUCE_AND";
    default: TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  return 0; // Will never be called
}

/** \brief Get the process rank.
 *
 * \relates Comm
 */
template<typename Ordinal>
int rank(const Comm<Ordinal>& comm);

/** \brief Get the number of processes in the communicator.
 *
 * \relates Comm
 */
template<typename Ordinal>
int size(const Comm<Ordinal>& comm);

/** \brief Barrier.
 *
 * \relates Comm
 */
template<typename Ordinal>
void barrier(const Comm<Ordinal>& comm);

/** \brief Broadcast array of objects that use value semantics.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void broadcast(
  const Comm<Ordinal>& comm,
  const int rootRank,
  const Ordinal count, Packet buffer[]
  );

/** \brief Broadcast array of objects that use value semantics.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void broadcast(
  const Comm<Ordinal>& comm,
  const int rootRank,
  const ArrayView<Packet> &buffer
  );

/** \brief Broadcast single object that use value semantics.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void broadcast(
  const Comm<Ordinal>& comm,
  const int rootRank, Packet *object
  );

/** \brief Broadcast single object that use value semantics.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void broadcast(
  const Comm<Ordinal>& comm,
  const int rootRank, const Ptr<Packet> &object
  );

/** \brief Broadcast array of objects that use reference semantics.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void broadcast(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const int rootRank, const Ordinal count, Packet*const buffer[]
  );

/** \brief Broadcast array of objects that use reference semantics.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void broadcast(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const int rootRank, const ArrayView<const Ptr<Packet> > &buffer
  );

/** \brief Broadcast array of objects that use value semantics using
 * customized serializer.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet, typename Serializer>
void broadcast(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const int rootRank,
  const Ordinal count, Packet buffer[]
  );

/** \brief Gather array of objects that use value semantics from every process
 * to every process.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void gatherAll(
  const Comm<Ordinal>& comm,
  const Ordinal sendCount, const Packet sendBuffer[],
  const Ordinal recvCount, Packet recvBuffer[]
  );

/** \brief Gather array of objects that use reference semantics from every
 * process to every process.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void gatherAll(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const Ordinal sendCount, const Packet*const sendBuffer[],
  const Ordinal recvCount, Packet*const recvBuffer[]
  );

/** \brief Gather array of objects that use value semantics from every process
 * to every process using customized serializer.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet, typename Serializer>
void gatherAll(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const Ordinal sendCount, const Packet sendBuffer[],
  const Ordinal recvCount, Packet recvBuffer[]
  );

/** \brief Collective reduce all of array of objects using value semantics
 * using a user-defined reduction operator.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void reduceAll(
  const Comm<Ordinal>& comm, const ValueTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal count, const Packet sendBuffer[], Packet globalReducts[]
  );

/** \brief Collective reduce all of array of objects using value semantics
 * using a pre-defined reduction type.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void reduceAll(
  const Comm<Ordinal>& comm, const EReductionType reductType,
  const Ordinal count, const Packet sendBuffer[], Packet globalReducts[]
  );

/** \brief Collective reduce all for single object using value semantics using
 * a pre-defined reduction type.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void reduceAll(
  const Comm<Ordinal>& comm, const EReductionType reductType,
  const Packet &send, const Ptr<Packet> &globalReduct
  );

/** \brief Deprecated . */
template<typename Ordinal, typename Packet>
TEUCHOS_DEPRECATED void reduceAll(
  const Comm<Ordinal>& comm, const EReductionType reductType,
  const Packet &send, Packet *globalReduct
  )
{
  reduceAll<Ordinal,Packet>(comm, reductType, send, ptr(globalReduct));
}

/** \brief Collective reduce all for array of objects using reference
 * semantics.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void reduceAll(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const ReferenceTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal count, const Packet*const sendBuffer[], Packet*const globalReducts[]
  );

/** \brief Collective reduce all of array of objects using value semantics
 * using a user-defined reduction operator and customized serializer.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet, typename Serializer>
void reduceAll(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const ValueTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal count, const Packet sendBuffer[], Packet globalReducts[]
  );

/** \brief Collective reduce all of array of objects using value semantics
 * using a pre-defined reduction type and customized serializer.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet, typename Serializer>
void reduceAll(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const EReductionType reductType,
  const Ordinal count, const Packet sendBuffer[], Packet globalReducts[]
  );

/** \brief Reduce and Scatter array of objects that use value semantics using
 * a user-defined reduction object.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void reduceAllAndScatter(
  const Comm<Ordinal>& comm, const ValueTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal sendCount, const Packet sendBuffer[] ,
  const Ordinal recvCounts[], Packet myGlobalReducts[]
  );

/** \brief Reduce and Scatter array of objects that use value semantics using
 * a a pre-defined reduction type.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void reduceAllAndScatter(
  const Comm<Ordinal>& comm, const EReductionType reductType,
  const Ordinal sendCount, const Packet sendBuffer[] ,
  const Ordinal recvCounts[], Packet myGlobalReducts[]
  );

/** \brief Reduce and Scatter array of objects that use reference semantics
 * using a user-defined reduction object.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void reduceAllAndScatter(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const ReferenceTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal sendCount, const Packet*const sendBuffer[] ,
  const Ordinal recvCounts[], Packet*const myGlobalReducts[]
  );

/** \brief Scan/Reduce array of objects that use value semantics using a
 * user-defined reduction operator.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void scan(
  const Comm<Ordinal>& comm, const ValueTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal count, const Packet sendBuffer[], Packet scanReducts[]
  );

/** \brief Reduce and Scatter array of objects that use value semantics using
 * a user-defined reduction object and customized serializer.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet, typename Serializer>
void reduceAllAndScatter(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const ValueTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal sendCount, const Packet sendBuffer[],
  const Ordinal recvCounts[], Packet myGlobalReducts[]
  );

/** \brief Reduce and Scatter array of objects that use value semantics using
 * a a pre-defined reduction type and customized serializer.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet, typename Serializer>
void reduceAllAndScatter(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const EReductionType reductType,
  const Ordinal sendCount, const Packet sendBuffer[] ,
  const Ordinal recvCounts[], Packet myGlobalReducts[]
  );

/** \brief Scan/Reduce array of objects using value semantics using a
 * predefined reduction type.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void scan(
  const Comm<Ordinal>& comm, const EReductionType reductType,
  const Ordinal count, const Packet sendBuffer[], Packet scanReducts[]
  );

/** \brief Scan/Reduce single object using value semantics using a predefined
 * reduction type.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void scan(
  const Comm<Ordinal>& comm, const EReductionType reductType,
  const Packet &send, const Ptr<Packet> &scanReduct
  );

/** \brief Deprecated. */
template<typename Ordinal, typename Packet>
TEUCHOS_DEPRECATED void scan(
  const Comm<Ordinal>& comm, const EReductionType reductType,
  const Packet &send, Packet *scanReduct
  )
{
  scan(comm, reductType, send, ptr(scanReduct));
}

/** \brief Scan/Reduce array of objects that use reference semantics using a
 * user-defined reduction operator.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void scan(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const ReferenceTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal count, const Packet*const sendBuffer[], Packet*const scanReducts[]
  );

/** \brief Scan/Reduce array of objects that use value semantics using a
 * user-defined reduction operator and customized serializer.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet, typename Serializer>
void scan(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const ValueTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal count, const Packet sendBuffer[], Packet scanReducts[]
  );

/** \brief Scan/Reduce array of objects using value semantics using a
 * predefined reduction type and customized serializer.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet, typename Serializer>
void scan(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const EReductionType reductType,
  const Ordinal count, const Packet sendBuffer[], Packet scanReducts[]
  );

/** \brief Send objects that use values semantics to another process.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void send(
  const Comm<Ordinal>& comm,
  const Ordinal count, const Packet sendBuffer[], const int destRank
  );

//! Variant of send() that takes a tag (and restores the correct order of arguments).
template<typename Ordinal, typename Packet>
void
send (const Packet sendBuffer[],
      const Ordinal count,
      const int destRank,
      const int tag,
      const Comm<Ordinal>& comm);

/** \brief Synchronously send objects that use values semantics to another process.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void ssend(
  const Comm<Ordinal>& comm,
  const Ordinal count, const Packet sendBuffer[], const int destRank
  );

/** \brief Send a single object that use values semantics to another process.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void send(
  const Comm<Ordinal>& comm,
  const Packet &send, const int destRank
  );

/** \brief Synchronously send a single object that use values semantics to another process.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void ssend(
  const Comm<Ordinal>& comm,
  const Packet &send, const int destRank
  );

/** \brief Send objects that use reference semantics to another process.
 *
 * \relates Comm
 *
 * NOTE: Not implemented yet!
 */
template<typename Ordinal, typename Packet>
void send(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const Ordinal count, const Packet*const sendBuffer[], const int destRank
  );

/** \brief Send objects that use values semantics to another process
 * using customized serializer.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet, typename Serializer>
void send(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const Ordinal count, const Packet sendBuffer[], const int destRank
  );

/** \brief Receive objects that use values semantics from another process.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
int receive(
  const Comm<Ordinal>& comm,
  const int sourceRank, const Ordinal count, Packet recvBuffer[]
  );

/** \brief Receive a single object that use values semantics from another process.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
int receive(
  const Comm<Ordinal>& comm,
  const int sourceRank, Packet *recv
  );

/** \brief Receive objects that use reference semantics from another process.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
int receive(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const int sourceRank, const Ordinal count, Packet*const recvBuffer[]
  );

/** \brief Receive objects that use values semantics from another process
 * using customized serializer.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet, typename Serializer>
int receive(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const int sourceRank, const Ordinal count, Packet recvBuffer[]
  );

/** \brief Ready-Send an array of objects that use values semantics to another
 * process.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void readySend(
  const Comm<Ordinal>& comm,
  const ArrayView<const Packet> &sendBuffer,
  const int destRank
  );

/** \brief Ready-Send a single object that use values semantics to another process.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
void readySend(
  const Comm<Ordinal>& comm,
  const Packet &send,
  const int destRank
  );

/** \brief Ready-Send an array of objects that use values semantics to another
 * process using customized serializer.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet, typename Serializer>
void readySend(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const ArrayView<const Packet> &sendBuffer,
  const int destRank
  );

/** \brief Send objects that use values semantics to another process.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
RCP<CommRequest<Ordinal> > isend(
  const Comm<Ordinal>& comm,
  const ArrayRCP<const Packet> &sendBuffer,
  const int destRank
  );

//! Variant of isend() that takes a tag (and restores the correct order of arguments).
template<typename Ordinal, typename Packet>
RCP<CommRequest<Ordinal> >
isend (const ArrayRCP<const Packet>& sendBuffer,
       const int destRank,
       const int tag,
       const Comm<Ordinal>& comm);

/** \brief Send a single object that use values semantics to another process.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
RCP<CommRequest<Ordinal> > isend(
  const Comm<Ordinal>& comm,
  const RCP<const Packet> &send,
  const int destRank
  );

/** \brief Send objects that use values semantics to another process
 * using customized serializer.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet, typename Serializer>
RCP<CommRequest<Ordinal> > isend(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const ArrayRCP<const Packet> &sendBuffer,
  const int destRank
  );


// 2008/07/29: rabartl: ToDo: Add reference semantics version of isend!


/// \brief Receive one or more objects (that use values semantics) from another process.
/// \relates Comm
///
/// \param comm [in] The communicator.
/// \param recvBuffer [out] The buffer into which to receive the data.
/// \param sourceRank [in] The rank of the sending process.  A
///   negative source rank means that this will accept an incoming
///   message from any process on the given communicator.  (This is
///   the equivalent of MPI_ANY_SOURCE.)
template<typename Ordinal, typename Packet>
RCP<CommRequest<Ordinal> > ireceive(
  const Comm<Ordinal>& comm,
  const ArrayRCP<Packet> &recvBuffer,
  const int sourceRank
  );

//! Variant of ireceive that takes a tag argument (and restores the correct order of arguments).
template<typename Ordinal, typename Packet>
RCP<CommRequest<Ordinal> >
ireceive (const ArrayRCP<Packet> &recvBuffer,
          const int sourceRank,
          const int tag,
          const Comm<Ordinal>& comm);

/// \brief Receive one object (that uses values semantics) from another process.
/// \relates Comm
///
/// \param comm [in] The communicator.
/// \param recv [out] The buffer into which to receive the object.
/// \param sourceRank [in] The rank of the sending process.  A
///   negative source rank means that this will accept an incoming
///   message from any process on the given communicator.
///
/// \note To implementers: A negative source rank is the equivalent of
///   MPI_ANY_SOURCE, if the given Comm is an MpiComm.
template<typename Ordinal, typename Packet>
RCP<CommRequest<Ordinal> > ireceive(
  const Comm<Ordinal>& comm,
  const RCP<Packet> &recv,
  const int sourceRank
  );

/** \brief Send objects that use values semantics to another process
 * using customized serializer.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet, typename Serializer>
RCP<CommRequest<Ordinal> > ireceive(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const ArrayRCP<Packet> &recvBuffer,
  const int sourceRank
  );


// 2008/07/29: rabartl: ToDo: Add reference semantics version of ireceive!


/** \brief Wait for an array of Teuchos::CommRequest objects.
 *
 * Blocks until all communication operations associated with the CommRequest
 * objects have completed.
 *
 * \relates Comm
 */
template<typename Ordinal>
void waitAll(
  const Comm<Ordinal>& comm,
  const ArrayView<RCP<CommRequest<Ordinal> > > &requests
  );

/// \brief Wait on one or more communication requests, and return their statuses.
/// \relates Comm
///
/// This function will block until all the communication operations
/// complete.  You are responsible for ordering messages in a way that
/// avoids deadlock.  An erroneous ordering of messages may cause this
/// function to wait forever.
///
/// \pre <tt>requests.size() == statuses.size()</tt>
///
/// \pre For i in 0, 1, ..., <tt>requests.size() - 1</tt>,
///   <tt>requests[i]</tt> is either null or it was the return value
///   of a nonblocking communication request.
///
/// \post For i in 0, 1, ..., <tt>requests.size() - 1</tt>,
///   <tt>requests[i]</tt> is null.  If <tt>requests[i]</tt> were
///   nonnull on input, then the corresponding nonblocking
///   communication operation completed successfully.
///
/// \param comm [in] The communicator on which the communication
///   requests were made.
/// \param requests [in/out] On input: the requests on which to wait.
///   Null elements of the array are ignored.  On output: all elements
///   of the array are set to null.
/// \param statuses [out] CommStatus instances representing the
///   results of completing the requests.  You may query each for
///   information about its corresponding communication operation.
///   Any element of the array may be null, for example if its
///   corresponding CommRequest was null on input.
template<typename Ordinal>
void
waitAll (const Comm<Ordinal>& comm,
         const ArrayView<RCP<CommRequest<Ordinal> > >& requests,
         const ArrayView<RCP<CommStatus<Ordinal> > >& statuses);

/// \brief Wait on a single communication request, and return its status.
/// \relates Comm
///
/// This function will block until the communication operation
/// completes.  You are responsible for ordering messages in a way
/// that avoids deadlock.  An erroneous ordering of messages may cause
/// this function to wait forever.
///
/// \param request [in/out] A nonnull pointer to an RCP.  On input,
///   the RCP is either null (in which case this function does
///   nothing) or it points to a valid CommRequest representing an
///   outstanding communication request.  On output, this function
///   sets <tt>*request</tt> to null, which indicates that the request
///   completed successfully.  (This function will not return unless
///   the request completes.)
///
/// \return A CommStatus instance representing the result of
///   completing the request.  You may query it for information about
///   the completed communication operation.  This may be null, for
///   example if <tt>*request</tt> was null on input.
///
/// \pre <tt>!is_null(request)</tt> (that is, the Ptr is not null).
/// \post <tt>is_null(*request)</tt> (that is, the RCP is null).
template<typename Ordinal>
RCP<CommStatus<Ordinal> >
wait (const Comm<Ordinal>& comm, const Ptr<RCP<CommRequest<Ordinal> > >& request);

//
// Standard reduction subclasses for objects that use value semantics
//


/** \brief Standard summation operator for types with value semantics.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
class SumValueReductionOp : public ValueTypeReductionOp<Ordinal,Packet>
{
public:
  /** \brief . */
  void reduce(
    const Ordinal count,
    const Packet inBuffer[],
    Packet inoutBuffer[]
    ) const;
};


/** \brief Standard min operator for types with value semantics.
 *
 * Note, this class object will throw an std::exception when used with a packet
 * type where <tt>ScalarTraits<Packet>::isComparable==false</tt> but it will
 * still compile.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
class MinValueReductionOp : public ValueTypeReductionOp<Ordinal,Packet>
{
public:
  /** \brief . */
  void reduce(
    const Ordinal count,
    const Packet inBuffer[],
    Packet inoutBuffer[]
    ) const;
};


/** \brief Standard Max operator for types with value semantics.
 *
 * Note, this class object will throw an std::exception when used with a packet
 * type where <tt>ScalarTraits<Packet>::isComparable==false</tt> but it will
 * still compile.
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
class MaxValueReductionOp : public ValueTypeReductionOp<Ordinal,Packet>
{
public:
  /** \brief . */
  void reduce(
    const Ordinal count,
    const Packet inBuffer[],
    Packet inoutBuffer[]
    ) const;
};


/** \brief Standard logical AND operator for booleans
 *
 * \relates Comm
 */
template<typename Ordinal, typename Packet>
class ANDValueReductionOp : public ValueTypeReductionOp<Ordinal,Packet>
{
public:
  /** \brief . */
  void reduce(
    const Ordinal count,
    const Packet inBuffer[],
    Packet inoutBuffer[]
    ) const;
};


// ////////////////////////////////////////////////////////////
// Implementation details (not for geneal users to mess with)


//
// ReductionOp Utilities
//


namespace MixMaxUtilities {


template<bool isComparable, typename Ordinal, typename Packet>
class Min {};


template<typename Ordinal, typename Packet>
class Min<true,Ordinal,Packet> {
public:
  static void min(
    const Ordinal count,
    const Packet inBuffer[],
    Packet inoutBuffer[]
    )
    {
      for( int i = 0; i < count; ++i )
        inoutBuffer[i] = TEUCHOS_MIN(inoutBuffer[i],inBuffer[i]);
    }
};


template<typename Ordinal, typename Packet>
class Min<false,Ordinal,Packet> {
public:
  static void min(
    const Ordinal,
    const Packet[],
    Packet[]
    )
    {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,std::logic_error,
        "Error, the type "<<TypeNameTraits<Packet>::name()
        <<" does not support comparison operations!"
        );
    }
};


template<bool isComparable, typename Ordinal, typename Packet>
class Max {};


template<typename Ordinal, typename Packet>
class Max<true,Ordinal,Packet> {
public:
  static void max(
    const Ordinal count,
    const Packet inBuffer[],
    Packet inoutBuffer[]
    )
    {
      for( int i = 0; i < count; ++i )
        inoutBuffer[i] = TEUCHOS_MAX(inoutBuffer[i],inBuffer[i]);
    }
};


template<typename Ordinal, typename Packet>
class Max<false,Ordinal,Packet> {
public:
  static void max(
    const Ordinal,
    const Packet[],
    Packet[]
    )
    {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,std::logic_error,
        "Error, the type "<<TypeNameTraits<Packet>::name()
        <<" does not support comparison operations!"
        );
    }
};


template<bool isComparable, typename Ordinal, typename Packet>
class AND {};


template<typename Ordinal, typename Packet>
class AND<true,Ordinal,Packet> {
public:
  static void andOp(
    const Ordinal count,
    const Packet inBuffer[],
    Packet inoutBuffer[]
    )
    {
      for( int i = 0; i < count; ++i )
        inoutBuffer[i] = inoutBuffer[i] && inBuffer[i];
    }
};


template<typename Ordinal, typename Packet>
class AND<false,Ordinal,Packet> {
public:
  static void andOp(
    const Ordinal,
    const Packet[],
    Packet[]
    )
    {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true,std::logic_error,
        "Error, the type "<<TypeNameTraits<Packet>::name()
        <<" does not support logical AND operations!"
        );
    }
};


} // namespace MixMaxUtilities


template<typename Ordinal, typename Packet>
void SumValueReductionOp<Ordinal,Packet>::reduce(
  const Ordinal count,
  const Packet inBuffer[],
  Packet inoutBuffer[]
  ) const
{
  for( int i = 0; i < count; ++i )
    inoutBuffer[i] += inBuffer[i];
}


template<typename Ordinal, typename Packet>
void MinValueReductionOp<Ordinal,Packet>::reduce(
  const Ordinal count,
  const Packet inBuffer[],
  Packet inoutBuffer[]
  ) const
{
  typedef MixMaxUtilities::Min<ScalarTraits<Packet>::isComparable, Ordinal, Packet> min_type;
  min_type::min (count, inBuffer, inoutBuffer);
}


template<typename Ordinal, typename Packet>
void MaxValueReductionOp<Ordinal,Packet>::reduce(
  const Ordinal count,
  const Packet inBuffer[],
  Packet inoutBuffer[]
  ) const
{
  typedef MixMaxUtilities::Max<ScalarTraits<Packet>::isComparable, Ordinal, Packet> max_type;
  max_type::max (count,inBuffer,inoutBuffer);
}


template<typename Ordinal, typename Packet>
void ANDValueReductionOp<Ordinal,Packet>::reduce(
  const Ordinal count,
  const Packet inBuffer[],
  Packet inoutBuffer[]
  ) const
{
  typedef MixMaxUtilities::AND<ScalarTraits<Packet>::isComparable, Ordinal, Packet> and_type;
  and_type::andOp (count, inBuffer, inoutBuffer);
}


} // namespace Teuchos


// //////////////////////////
// Template implemenations


//
// ReductionOp utilities
//


namespace Teuchos {


// Not for the general user to use! I am returning a raw ReducionOp* pointer
// to avoid the overhead of using RCP. However, given the use case
// this is just fine since I can just use std::auto_ptr to make sure things
// are deleted correctly.
template<typename Ordinal, typename Packet>
ValueTypeReductionOp<Ordinal,Packet>*
createOp( const EReductionType reductType )
{
  typedef ScalarTraits<Packet> ST;
  switch(reductType) {
    case REDUCE_SUM: {
      return new SumValueReductionOp<Ordinal,Packet>();
      break;
    }
    case REDUCE_MIN: {
      TEUCHOS_TEST_FOR_EXCEPT(!ST::isComparable);
      return new MinValueReductionOp<Ordinal,Packet>();
      break;
    }
    case REDUCE_MAX: {
      TEUCHOS_TEST_FOR_EXCEPT(!ST::isComparable);
      return new MaxValueReductionOp<Ordinal,Packet>();
      break;
    }
    case REDUCE_AND: {
      return new ANDValueReductionOp<Ordinal, Packet>();
      break;
    }
    default:
      TEUCHOS_TEST_FOR_EXCEPT(true);
  }
  return 0; // Will never be called!
}


} // namespace Teuchos


//
// Teuchos::Comm wrapper functions
//


template<typename Ordinal>
int Teuchos::rank(const Comm<Ordinal>& comm)
{
  return comm.getRank();
}


template<typename Ordinal>
int Teuchos::size(const Comm<Ordinal>& comm)
{
  return comm.getSize();
}


template<typename Ordinal>
void Teuchos::barrier(const Comm<Ordinal>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: barrier<"
    <<OrdinalTraits<Ordinal>::name()
    <<">()"
    );
  comm.barrier();
}


template<typename Ordinal, typename Packet>
void Teuchos::broadcast(
  const Comm<Ordinal>& comm,
  const int rootRank, const Ordinal count, Packet buffer[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: broadcast<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type )"
    );
  ValueTypeSerializationBuffer<Ordinal,Packet>
    charBuffer(count,buffer);
  comm.broadcast(
    rootRank,charBuffer.getBytes(),charBuffer.getCharBuffer()
    );
}


template<typename Ordinal, typename Packet>
void Teuchos::broadcast(
  const Comm<Ordinal>& comm,
  const int rootRank,
  const ArrayView<Packet> &buffer
  )
{
  broadcast<Ordinal, Packet>(comm, rootRank, buffer.size(), buffer.getRawPtr() );
}


template<typename Ordinal, typename Packet>
void Teuchos::broadcast(
  const Comm<Ordinal>& comm,
  const int rootRank, Packet *object
  )
{
  broadcast<Ordinal,Packet>(comm,rootRank,1,object);
}


template<typename Ordinal, typename Packet>
void Teuchos::broadcast(
  const Comm<Ordinal>& comm,
  const int rootRank, const Ptr<Packet> &object
  )
{
  broadcast<Ordinal,Packet>(comm,rootRank,1,object.getRawPtr());
}


template<typename Ordinal, typename Packet>
void Teuchos::broadcast(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const int rootRank, const Ordinal count, Packet*const buffer[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: broadcast<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( reference type )"
    );
  ReferenceTypeSerializationBuffer<Ordinal,Packet>
    charBuffer(serializer, count, buffer);
  comm.broadcast(
    rootRank,charBuffer.getBytes(),charBuffer.getCharBuffer()
    );
}


template<typename Ordinal, typename Packet>
void Teuchos::broadcast(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const int rootRank, const ArrayView<const Ptr<Packet> > &buffer
  )
{
  Array<Packet*> bufferPtrArray;
  for (int i = 0; i < buffer.size(); ++i) {
    bufferPtrArray.push_back(buffer[i].getRawPtr());
  }
  broadcast<Ordinal,Packet>(comm, serializer, rootRank,
    buffer.size(), bufferPtrArray.getRawPtr());
}

template<typename Ordinal, typename Packet, typename Serializer>
void Teuchos::broadcast(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const int rootRank, const Ordinal count, Packet buffer[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: broadcast<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type )"
    );
  ValueTypeSerializationBuffer<Ordinal,Packet,Serializer>
    charBuffer(count,buffer,rcp(&serializer,false));
  comm.broadcast(
    rootRank,charBuffer.getBytes(),charBuffer.getCharBuffer()
    );
}


template<typename Ordinal, typename Packet>
void Teuchos::gatherAll(
  const Comm<Ordinal>& comm,
  const Ordinal sendCount, const Packet sendBuffer[],
  const Ordinal recvCount, Packet recvBuffer[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: gatherAll<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type )"
    );
  ConstValueTypeSerializationBuffer<Ordinal,Packet>
    charSendBuffer(sendCount,sendBuffer);
  ValueTypeSerializationBuffer<Ordinal,Packet>
    charRecvBuffer(recvCount,recvBuffer);
  comm.gatherAll(
    charSendBuffer.getBytes(),charSendBuffer.getCharBuffer()
    ,charRecvBuffer.getBytes(),charRecvBuffer.getCharBuffer()
    );
}


template<typename Ordinal, typename Packet>
void Teuchos::gatherAll(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const Ordinal sendCount, const Packet*const sendBuffer[],
  const Ordinal recvCount, Packet*const recvBuffer[]
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement and test when needed!
}

template<typename Ordinal, typename Packet, typename Serializer>
void Teuchos::gatherAll(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const Ordinal sendCount, const Packet sendBuffer[],
  const Ordinal recvCount, Packet recvBuffer[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: gatherAll<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type )"
    );
  ConstValueTypeSerializationBuffer<Ordinal,Packet,Serializer>
    charSendBuffer(sendCount,sendBuffer,rcp(&serializer,false));
  ValueTypeSerializationBuffer<Ordinal,Packet,Serializer>
    charRecvBuffer(recvCount,recvBuffer,rcp(&serializer,false));
  comm.gatherAll(
    charSendBuffer.getBytes(),charSendBuffer.getCharBuffer()
    ,charRecvBuffer.getBytes(),charRecvBuffer.getCharBuffer()
    );
}


template<typename Ordinal, typename Packet>
void Teuchos::reduceAll(
  const Comm<Ordinal>& comm, const ValueTypeReductionOp<Ordinal,Packet> &reductOp
  ,const Ordinal count, const Packet sendBuffer[], Packet globalReducts[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: reduceAll<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type, user-defined op )"
    );
  ConstValueTypeSerializationBuffer<Ordinal,Packet>
    charSendBuffer(count,sendBuffer);
  ValueTypeSerializationBuffer<Ordinal,Packet>
    charGlobalReducts(count,globalReducts);
  CharToValueTypeReductionOp<Ordinal,Packet>
    charReductOp(rcp(&reductOp,false));
  comm.reduceAll(
    charReductOp,charSendBuffer.getBytes(),charSendBuffer.getCharBuffer()
    ,charGlobalReducts.getCharBuffer()
    );
}


template<typename Ordinal, typename Packet>
void Teuchos::reduceAll(
  const Comm<Ordinal>& comm, const EReductionType reductType,
  const Ordinal count, const Packet sendBuffer[], Packet globalReducts[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: reduceAll<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type, "<<toString(reductType)<<" )"
    );
  std::auto_ptr<ValueTypeReductionOp<Ordinal,Packet> >
    reductOp(createOp<Ordinal,Packet>(reductType));
  reduceAll(comm,*reductOp,count,sendBuffer,globalReducts);
}


namespace Teuchos {

#ifdef TEUCHOS_HAVE_COMPLEX
// Specialization for Ordinal=int and Packet=std::complex<double>.
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
reduceAll<int, std::complex<double> > (const Comm<int>& comm,
                                       const EReductionType reductType,
                                       const int count,
                                       const std::complex<double> sendBuffer[],
                                       std::complex<double> globalReducts[]);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, std::complex<double> > (const Comm<int>& comm,
                                      const ArrayRCP<std::complex<double> >& recvBuffer,
                                      const int sourceRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, std::complex<double> > (const ArrayRCP<std::complex<double> > &recvBuffer,
                                      const int sourceRank,
                                      const int tag,
                                      const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, std::complex<double> > (const Comm<int>& comm,
                                  const int count,
                                  const std::complex<double> sendBuffer[],
                                  const int destRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, std::complex<double> > (const std::complex<double> sendBuffer[],
                                  const int count,
                                  const int destRank,
                                  const int tag,
                                  const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
isend<int, std::complex<double> > (const ArrayRCP<const std::complex<double> >& sendBuffer,
                                   const int destRank,
                                   const int tag,
                                   const Comm<int>& comm);

// Specialization for Ordinal=int and Packet=std::complex<float>.
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
reduceAll<int, std::complex<float> > (const Comm<int>& comm,
                                      const EReductionType reductType,
                                      const int count,
                                      const std::complex<float> sendBuffer[],
                                      std::complex<float> globalReducts[]);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, std::complex<float> > (const Comm<int>& comm,
                                     const ArrayRCP<std::complex<float> >& recvBuffer,
                                     const int sourceRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, std::complex<float> > (const ArrayRCP<std::complex<float> > &recvBuffer,
                                     const int sourceRank,
                                     const int tag,
                                     const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, std::complex<float> > (const Comm<int>& comm,
                                 const int count,
                                 const std::complex<float> sendBuffer[],
                                 const int destRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, std::complex<float> > (const std::complex<float> sendBuffer[],
                                 const int count,
                                 const int destRank,
                                 const int tag,
                                 const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
isend<int, std::complex<float> > (const ArrayRCP<const std::complex<float> >& sendBuffer,
                                  const int destRank,
                                  const int tag,
                                  const Comm<int>& comm);
#endif // TEUCHOS_HAVE_COMPLEX

// Specialization for Ordinal=int and Packet=double.
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
reduceAll<int, double> (const Comm<int>& comm,
                        const EReductionType reductType,
                        const int count,
                        const double sendBuffer[],
                        double globalReducts[]);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, double> (const Comm<int>& comm,
                       const ArrayRCP<double>& recvBuffer,
                       const int sourceRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, double> (const ArrayRCP<double> &recvBuffer,
                       const int sourceRank,
                       const int tag,
                       const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, double> (const Comm<int>& comm,
                   const int count,
                   const double sendBuffer[],
                   const int destRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, double> (const double sendBuffer[],
                   const int count,
                   const int destRank,
                   const int tag,
                   const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
isend<int, double> (const ArrayRCP<const double>& sendBuffer,
                    const int destRank,
                    const int tag,
                    const Comm<int>& comm);

// Specialization for Ordinal=int and Packet=float.
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
reduceAll<int, float> (const Comm<int>& comm,
                       const EReductionType reductType,
                       const int count,
                       const float sendBuffer[],
                       float globalReducts[]);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, float> (const Comm<int>& comm,
                      const ArrayRCP<float>& recvBuffer,
                      const int sourceRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, float> (const ArrayRCP<float> &recvBuffer,
                      const int sourceRank,
                      const int tag,
                      const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, float> (const Comm<int>& comm,
                  const int count,
                  const float sendBuffer[],
                  const int destRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, float> (const float sendBuffer[],
                  const int count,
                  const int destRank,
                  const int tag,
                  const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
isend<int, float> (const ArrayRCP<const float>& sendBuffer,
                   const int destRank,
                   const int tag,
                   const Comm<int>& comm);

#ifdef TEUCHOS_HAVE_LONG_LONG_INT
// Specialization for Ordinal=int and Packet=long long.
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
reduceAll<int, long long> (const Comm<int>& comm,
                           const EReductionType reductType,
                           const int count,
                           const long long sendBuffer[],
                           long long globalReducts[]);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, long long> (const Comm<int>& comm,
                          const ArrayRCP<long long>& recvBuffer,
                          const int sourceRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, long long> (const ArrayRCP<long long> &recvBuffer,
                          const int sourceRank,
                          const int tag,
                          const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, long long> (const Comm<int>& comm,
                      const int count,
                      const long long sendBuffer[],
                      const int destRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, long long> (const long long sendBuffer[],
                      const int count,
                      const int destRank,
                      const int tag,
                      const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
isend<int, long long> (const ArrayRCP<const long long>& sendBuffer,
                       const int destRank,
                       const int tag,
                       const Comm<int>& comm);
#endif // TEUCHOS_HAVE_LONG_LONG_INT

// Specialization for Ordinal=int and Packet=long.
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
reduceAll<int, long> (const Comm<int>& comm,
                      const EReductionType reductType,
                      const int count,
                      const long sendBuffer[],
                      long globalReducts[]);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, long> (const Comm<int>& comm,
                     const ArrayRCP<long>& recvBuffer,
                     const int sourceRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, long> (const ArrayRCP<long> &recvBuffer,
                     const int sourceRank,
                     const int tag,
                     const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, long> (const Comm<int>& comm,
                 const int count,
                 const long sendBuffer[],
                 const int destRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, long> (const long sendBuffer[],
                 const int count,
                 const int destRank,
                 const int tag,
                 const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
isend<int, long> (const ArrayRCP<const long>& sendBuffer,
                  const int destRank,
                  const int tag,
                  const Comm<int>& comm);

// Specialization for Ordinal=int and Packet=int.
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
reduceAll<int, int> (const Comm<int>& comm,
                     const EReductionType reductType,
                     const int count,
                     const int sendBuffer[],
                     int globalReducts[]);
// mfh 09 Apr 2013: We provide a full specialization of
// reduceAllAndScatter, because it is an important part of
// Tpetra::Distributor initialization.
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
reduceAllAndScatter<int, int> (const Comm<int>& comm,
                               const EReductionType reductType,
                               const int sendCount,
                               const int sendBuffer[],
                               const int recvCounts[],
                               int myGlobalReducts[]);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, int> (const Comm<int>& comm,
                    const ArrayRCP<int>& recvBuffer,
                    const int sourceRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, int> (const ArrayRCP<int> &recvBuffer,
                    const int sourceRank,
                    const int tag,
                    const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, int> (const Comm<int>& comm,
                const int count,
                const int sendBuffer[],
                const int destRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, int> (const int sendBuffer[],
                const int count,
                const int destRank,
                const int tag,
                const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
isend<int, int> (const ArrayRCP<const int>& sendBuffer,
                 const int destRank,
                 const int tag,
                 const Comm<int>& comm);

// Specialization for Ordinal=int and Packet=short.
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
reduceAll<int, short> (const Comm<int>& comm,
                       const EReductionType reductType,
                       const int count,
                       const short sendBuffer[],
                       short globalReducts[]);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, short> (const Comm<int>& comm,
                      const ArrayRCP<short>& recvBuffer,
                      const int sourceRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
ireceive<int, short> (const ArrayRCP<short> &recvBuffer,
                      const int sourceRank,
                      const int tag,
                      const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, short> (const Comm<int>& comm,
                  const int count,
                  const short sendBuffer[],
                  const int destRank);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT void
send<int, short> (const short sendBuffer[],
                  const int count,
                  const int destRank,
                  const int tag,
                  const Comm<int>& comm);
template<>
TEUCHOSCOMM_LIB_DLL_EXPORT RCP<CommRequest<int> >
isend<int, short> (const ArrayRCP<const short>& sendBuffer,
                   const int destRank,
                   const int tag,
                   const Comm<int>& comm);

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
TEUCHOSCOMM_LIB_DLL_EXPORT void
reduceAll<int, char> (const Comm<int>& comm,
                      const EReductionType reductType,
                      const int count,
                      const char sendBuffer[],
                      char globalReducts[]);
#endif // 0
} // namespace Teuchos


template<typename Ordinal, typename Packet>
void Teuchos::reduceAll(
  const Comm<Ordinal>& comm, const EReductionType reductType
  ,const Packet &send, const Ptr<Packet> &globalReduct
  )
{
  // mfh 17 Oct 2012: This will invoke the above specializations for
  // general count, so we don't need to specialize this function.
  reduceAll<Ordinal,Packet>(comm, reductType, 1, &send, &*globalReduct);
}


template<typename Ordinal, typename Packet>
void Teuchos::reduceAll(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const ReferenceTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal count, const Packet*const sendBuffer[], Packet*const globalReducts[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: reduceAll<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( reference type )"
    );
  ConstReferenceTypeSerializationBuffer<Ordinal,Packet>
    charSendBuffer(serializer,count,sendBuffer);
  ReferenceTypeSerializationBuffer<Ordinal,Packet>
    charGlobalReducts(serializer,count,globalReducts);
  CharToReferenceTypeReductionOp<Ordinal,Packet>
    charReductOp(rcp(&serializer,false),rcp(&reductOp,false));
  comm.reduceAll(
    charReductOp,charSendBuffer.getBytes(),charSendBuffer.getCharBuffer()
    ,charGlobalReducts.getCharBuffer()
    );
}

template<typename Ordinal, typename Packet, typename Serializer>
void Teuchos::reduceAll(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const ValueTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal count, const Packet sendBuffer[], Packet globalReducts[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: reduceAll<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type, user-defined op )"
    );
  ConstValueTypeSerializationBuffer<Ordinal,Packet,Serializer>
    charSendBuffer(count,sendBuffer,rcp(&serializer,false));
  ValueTypeSerializationBuffer<Ordinal,Packet,Serializer>
    charGlobalReducts(count,globalReducts,rcp(&serializer,false));
  CharToValueTypeReductionOp<Ordinal,Packet,Serializer>
    charReductOp(rcp(&reductOp,false),rcp(&serializer,false));
  comm.reduceAll(
    charReductOp,charSendBuffer.getBytes(),charSendBuffer.getCharBuffer()
    ,charGlobalReducts.getCharBuffer()
    );
}


template<typename Ordinal, typename Packet, typename Serializer>
void Teuchos::reduceAll(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const EReductionType reductType,
  const Ordinal count, const Packet sendBuffer[], Packet globalReducts[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: reduceAll<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type, "<<toString(reductType)<<" )"
    );
  std::auto_ptr<ValueTypeReductionOp<Ordinal,Packet> >
    reductOp(createOp<Ordinal,Packet>(reductType));
  reduceAll(comm,serializer,*reductOp,count,sendBuffer,globalReducts);
}


template<typename Ordinal, typename Packet>
void Teuchos::reduceAllAndScatter(
  const Comm<Ordinal>& comm, const ValueTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal sendCount, const Packet sendBuffer[] ,
  const Ordinal recvCounts[], Packet myGlobalReducts[]
  )
{

  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: reduceAllAndScatter<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type, user-defined op )"
    );

  const Ordinal size = Teuchos::size(comm);
  const Ordinal rank = Teuchos::rank(comm);

#ifdef TEUCHOS_DEBUG
  Ordinal sumRecvCounts = 0;
  for( Ordinal i = 0; i < size; ++i )
    sumRecvCounts += recvCounts[static_cast<ptrdiff_t>(i)];
  TEUCHOS_TEST_FOR_EXCEPT(!(sumRecvCounts==sendCount));
#endif

  ConstValueTypeSerializationBuffer<Ordinal,Packet>
    charSendBuffer(sendCount,sendBuffer);
  ValueTypeSerializationBuffer<Ordinal,Packet>
    charMyGlobalReducts(recvCounts[static_cast<ptrdiff_t>(rank)], myGlobalReducts);
  CharToValueTypeReductionOp<Ordinal,Packet>
    charReductOp(rcp(&reductOp,false));

  const Ordinal
    packetSize = charSendBuffer.getBytes()/sendCount;

  WorkspaceStore* wss = get_default_workspace_store().get();
  Workspace<Ordinal> charRecvCounts(wss, size);
  for (Ordinal k = 0; k < size; ++k) {
    charRecvCounts[k] = as<Ordinal>(recvCounts[static_cast<ptrdiff_t>(k)] * packetSize);
  }

  comm.reduceAllAndScatter(
    charReductOp, charSendBuffer.getBytes(), charSendBuffer.getCharBuffer(),
    &charRecvCounts[0], charMyGlobalReducts.getCharBuffer()
    );

}


template<typename Ordinal, typename Packet>
void Teuchos::reduceAllAndScatter(
  const Comm<Ordinal>& comm, const EReductionType reductType,
  const Ordinal sendCount, const Packet sendBuffer[] ,
  const Ordinal recvCounts[], Packet myGlobalReducts[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: reduceAllAndScatter<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type, "<<toString(reductType)<<" )"
    );
  std::auto_ptr<ValueTypeReductionOp<Ordinal,Packet> >
    reductOp(createOp<Ordinal,Packet>(reductType));
  reduceAllAndScatter(
    comm, *reductOp, sendCount, sendBuffer, recvCounts, myGlobalReducts
    );
}


template<typename Ordinal, typename Packet>
void Teuchos::reduceAllAndScatter(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const ReferenceTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal sendCount, const Packet*const sendBuffer[],
  const Ordinal recvCounts[], Packet*const myGlobalReducts[]
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement and test when needed!
}

template<typename Ordinal, typename Packet, typename Serializer>
void Teuchos::reduceAllAndScatter(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const ValueTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal sendCount, const Packet sendBuffer[] ,
  const Ordinal recvCounts[], Packet myGlobalReducts[]
  )
{

  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: reduceAllAndScatter<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type, user-defined op )"
    );

  const Ordinal size = Teuchos::size(comm);
  const Ordinal rank = Teuchos::rank(comm);

#ifdef TEUCHOS_DEBUG
  Ordinal sumRecvCounts = 0;
  for( Ordinal i = 0; i < size; ++i )
    sumRecvCounts += recvCounts[static_cast<ptrdiff_t>(i)];
  TEUCHOS_TEST_FOR_EXCEPT(!(sumRecvCounts==sendCount));
#endif

  ConstValueTypeSerializationBuffer<Ordinal,Packet,Serializer>
    charSendBuffer(sendCount,sendBuffer,rcp(&serializer,false));
  ValueTypeSerializationBuffer<Ordinal,Packet,Serializer>
    charMyGlobalReducts(recvCounts[static_cast<ptrdiff_t>(rank)], myGlobalReducts,rcp(&serializer,false));
  CharToValueTypeReductionOp<Ordinal,Packet,Serializer>
    charReductOp(rcp(&reductOp,false),rcp(&serializer,false));

  const Ordinal
    packetSize = charSendBuffer.getBytes()/sendCount;

  WorkspaceStore* wss = get_default_workspace_store().get();
  Workspace<Ordinal> charRecvCounts(wss, size);
  for (Ordinal k = 0; k < size; ++k) {
    charRecvCounts[k] = as<Ordinal>(recvCounts[static_cast<ptrdiff_t>(k)] * packetSize);
  }

  comm.reduceAllAndScatter(
    charReductOp, charSendBuffer.getBytes(), charSendBuffer.getCharBuffer(),
    &charRecvCounts[0], charMyGlobalReducts.getCharBuffer()
    );

}

template<typename Ordinal, typename Packet, typename Serializer>
void Teuchos::reduceAllAndScatter(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const EReductionType reductType,
  const Ordinal sendCount, const Packet sendBuffer[] ,
  const Ordinal recvCounts[], Packet myGlobalReducts[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: reduceAllAndScatter<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type, "<<toString(reductType)<<" )"
    );
  std::auto_ptr<ValueTypeReductionOp<Ordinal,Packet> >
    reductOp(createOp<Ordinal,Packet>(reductType));
  reduceAllAndScatter(
    comm, serializer, *reductOp, sendCount, sendBuffer, recvCounts,
    myGlobalReducts
    );
}


template<typename Ordinal, typename Packet>
void Teuchos::scan(
  const Comm<Ordinal>& comm, const ValueTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal count, const Packet sendBuffer[], Packet scanReducts[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: scan<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type, user-defined op )"
    );
  ConstValueTypeSerializationBuffer<Ordinal,Packet>
    charSendBuffer(count,sendBuffer);
  ValueTypeSerializationBuffer<Ordinal,Packet>
    charScanReducts(count,scanReducts);
  CharToValueTypeReductionOp<Ordinal,Packet>
    charReductOp(rcp(&reductOp,false));
  comm.scan(
    charReductOp,charSendBuffer.getBytes(),charSendBuffer.getCharBuffer()
    ,charScanReducts.getCharBuffer()
    );
}


template<typename Ordinal, typename Packet>
void Teuchos::scan(
  const Comm<Ordinal>& comm, const EReductionType reductType,
  const Ordinal count, const Packet sendBuffer[], Packet scanReducts[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: scan<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type, "<<toString(reductType)<<" )"
    );
  std::auto_ptr<ValueTypeReductionOp<Ordinal,Packet> >
    reductOp(createOp<Ordinal,Packet>(reductType));
  scan(comm,*reductOp,count,sendBuffer,scanReducts);
}


template<typename Ordinal, typename Packet>
void Teuchos::scan(
  const Comm<Ordinal>& comm, const EReductionType reductType,
  const Packet &send, const Ptr<Packet> &scanReduct
  )
{
  scan<Ordinal,Packet>(comm, reductType, 1, &send, &*scanReduct);
}


template<typename Ordinal, typename Packet>
void Teuchos::scan(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const ReferenceTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal count, const Packet*const sendBuffer[], Packet*const scanReducts[]
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement and test when needed!
}

template<typename Ordinal, typename Packet, typename Serializer>
void Teuchos::scan(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const ValueTypeReductionOp<Ordinal,Packet> &reductOp,
  const Ordinal count, const Packet sendBuffer[], Packet scanReducts[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: scan<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type, user-defined op )"
    );
  ConstValueTypeSerializationBuffer<Ordinal,Packet,Serializer>
    charSendBuffer(count,sendBuffer,rcp(&serializer,false));
  ValueTypeSerializationBuffer<Ordinal,Packet,Serializer>
    charScanReducts(count,scanReducts,rcp(&serializer,false));
  CharToValueTypeReductionOp<Ordinal,Packet,Serializer>
    charReductOp(rcp(&reductOp,false),rcp(&serializer,false));
  comm.scan(
    charReductOp,charSendBuffer.getBytes(),charSendBuffer.getCharBuffer()
    ,charScanReducts.getCharBuffer()
    );
}


template<typename Ordinal, typename Packet, typename Serializer>
void Teuchos::scan(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const EReductionType reductType,
  const Ordinal count, const Packet sendBuffer[], Packet scanReducts[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: scan<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type, "<<toString(reductType)<<" )"
    );
  std::auto_ptr<ValueTypeReductionOp<Ordinal,Packet> >
    reductOp(createOp<Ordinal,Packet>(reductType));
  scan(comm,serializer,*reductOp,count,sendBuffer,scanReducts);
}

template<typename Ordinal, typename Packet>
void Teuchos::send(
  const Comm<Ordinal>& comm,
  const Ordinal count, const Packet sendBuffer[], const int destRank
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: send<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type )"
    );
  ConstValueTypeSerializationBuffer<Ordinal,Packet>
    charSendBuffer(count,sendBuffer);
  comm.send(
    charSendBuffer.getBytes(),charSendBuffer.getCharBuffer()
    ,destRank
    );
}

template<typename Ordinal, typename Packet>
void
Teuchos::send (const Packet sendBuffer[],
               const Ordinal count,
               const int destRank,
               const int tag,
               const Comm<Ordinal>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: send<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type )"
    );
  ConstValueTypeSerializationBuffer<Ordinal,Packet> charSendBuffer (count, sendBuffer);
  comm.send (charSendBuffer.getBytes (), charSendBuffer.getCharBuffer (), destRank, tag);
}

template<typename Ordinal, typename Packet>
void Teuchos::ssend(
  const Comm<Ordinal>& comm,
  const Ordinal count, const Packet sendBuffer[], const int destRank
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: ssend<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type )"
    );
  ConstValueTypeSerializationBuffer<Ordinal,Packet>
    charSendBuffer(count,sendBuffer);
  comm.ssend(
    charSendBuffer.getBytes(),charSendBuffer.getCharBuffer()
    ,destRank
    );
}

template<typename Ordinal, typename Packet>
void Teuchos::send(
  const Comm<Ordinal>& comm,
  const Packet &send, const int destRank
  )
{
  Teuchos::send<Ordinal,Packet>(comm,1,&send,destRank);
}

template<typename Ordinal, typename Packet>
void Teuchos::ssend(
  const Comm<Ordinal>& comm,
  const Packet &send, const int destRank
  )
{
  Teuchos::ssend<Ordinal,Packet>(comm,1,&send,destRank);
}

template<typename Ordinal, typename Packet>
void Teuchos::send(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const Ordinal count, const Packet*const sendBuffer[], const int destRank
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement and test when needed!
}

template<typename Ordinal, typename Packet, typename Serializer>
void Teuchos::send(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const Ordinal count, const Packet sendBuffer[], const int destRank
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: send<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type )"
    );
  ConstValueTypeSerializationBuffer<Ordinal,Packet,Serializer>
    charSendBuffer(count,sendBuffer,rcp(&serializer,false));
  comm.send(
    charSendBuffer.getBytes(),charSendBuffer.getCharBuffer()
    ,destRank
    );
}

template<typename Ordinal, typename Packet>
int Teuchos::receive(
  const Comm<Ordinal>& comm,
  const int sourceRank, const Ordinal count, Packet recvBuffer[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: receive<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type )"
    );
  ValueTypeSerializationBuffer<Ordinal,Packet>
    charRecvBuffer(count,recvBuffer);
  return comm.receive(
    sourceRank
    ,charRecvBuffer.getBytes(),charRecvBuffer.getCharBuffer()
    );
}


template<typename Ordinal, typename Packet>
int Teuchos::receive(
  const Comm<Ordinal>& comm,
  const int sourceRank, Packet *recv
  )
{
  return Teuchos::receive<Ordinal,Packet>(comm,sourceRank,1,recv);
}


template<typename Ordinal, typename Packet>
int Teuchos::receive(
  const Comm<Ordinal>& comm, const Serializer<Ordinal,Packet> &serializer,
  const int sourceRank, const Ordinal count, Packet*const recvBuffer[]
  )
{
  TEUCHOS_TEST_FOR_EXCEPT(true); // ToDo: Implement and test when needed!
}

template<typename Ordinal, typename Packet, typename Serializer>
int Teuchos::receive(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const int sourceRank, const Ordinal count, Packet recvBuffer[]
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: receive<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type )"
    );
  ValueTypeSerializationBuffer<Ordinal,Packet,Serializer>
    charRecvBuffer(count,recvBuffer,rcp(&serializer,false));
  return comm.receive(
    sourceRank
    ,charRecvBuffer.getBytes(),charRecvBuffer.getCharBuffer()
    );
}

template<typename Ordinal, typename Packet>
void Teuchos::readySend(
  const Comm<Ordinal>& comm,
  const ArrayView<const Packet> &sendBuffer,
  const int destRank
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: readySend<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type )"
    );
  ConstValueTypeSerializationBuffer<Ordinal,Packet>
    charSendBuffer(sendBuffer.size(), sendBuffer.getRawPtr());
  comm.readySend( charSendBuffer.getCharBufferView(), destRank );
}


template<typename Ordinal, typename Packet>
void Teuchos::readySend(
  const Comm<Ordinal>& comm,
  const Packet &send,
  const int destRank
  )
{
  readySend<Ordinal, Packet>( comm, arrayView(&send,1), destRank );
}

template<typename Ordinal, typename Packet, typename Serializer>
void Teuchos::readySend(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const ArrayView<const Packet> &sendBuffer,
  const int destRank
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: readySend<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type )"
    );
  ConstValueTypeSerializationBuffer<Ordinal,Packet,Serializer>
    charSendBuffer(sendBuffer.size(), sendBuffer.getRawPtr(), serializer);
  comm.readySend( charSendBuffer.getCharBufferView(), destRank );
}

template<typename Ordinal, typename Packet>
Teuchos::RCP<Teuchos::CommRequest<Ordinal> >
Teuchos::isend(
  const Comm<Ordinal>& comm,
  const ArrayRCP<const Packet> &sendBuffer,
  const int destRank
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: isend<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type )"
    );
  ConstValueTypeSerializationBuffer<Ordinal,Packet>
    charSendBuffer(sendBuffer.size(), sendBuffer.getRawPtr());
  RCP<CommRequest<Ordinal> > commRequest = comm.isend(
    charSendBuffer.getCharBufferView(), destRank );
  set_extra_data( sendBuffer, "buffer", inOutArg(commRequest) );
  return commRequest;
}

template<typename Ordinal, typename Packet>
Teuchos::RCP<Teuchos::CommRequest<Ordinal> >
Teuchos::isend (const ArrayRCP<const Packet> &sendBuffer,
                const int destRank,
                const int tag,
                const Comm<Ordinal>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::isend<" << OrdinalTraits<Ordinal>::name () << ","
    << TypeNameTraits<Packet>::name () << ">");
  ConstValueTypeSerializationBuffer<Ordinal,Packet>
    charSendBuffer (sendBuffer.size (), sendBuffer.getRawPtr ());
  RCP<CommRequest<Ordinal> > commRequest =
    comm.isend (charSendBuffer.getCharBufferView (), destRank, tag);
  set_extra_data (sendBuffer, "buffer", inOutArg (commRequest));
  return commRequest;
}

template<typename Ordinal, typename Packet>
Teuchos::RCP<Teuchos::CommRequest<Ordinal> >
Teuchos::isend(
  const Comm<Ordinal>& comm,
  const RCP<const Packet> &send,
  const int destRank
  )
{
  const ArrayRCP<const Packet> sendBuffer =
    arcpWithEmbeddedObj( send.get(), 0, 1, send, false );
  // 2008/07/29: rabartl: Above: I need to write a helper function to create
  // new ArrayRCP object given a single object to copy.
  return isend<Ordinal, Packet>( comm, sendBuffer, destRank );
}

template<typename Ordinal, typename Packet, typename Serializer>
Teuchos::RCP<Teuchos::CommRequest<Ordinal> >
Teuchos::isend(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const ArrayRCP<const Packet> &sendBuffer,
  const int destRank
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: isend<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type )"
    );
  ConstValueTypeSerializationBuffer<Ordinal,Packet,Serializer>
    charSendBuffer(sendBuffer.size(), sendBuffer.getRawPtr(), serializer);
  RCP<CommRequest<Ordinal> > commRequest = comm.isend(
    charSendBuffer.getCharBufferView(), destRank );
  set_extra_data( sendBuffer, "buffer", inOutArg(commRequest) );
  return commRequest;
}

template<typename Ordinal, typename Packet>
Teuchos::RCP<Teuchos::CommRequest<Ordinal> >
Teuchos::ireceive(
  const Comm<Ordinal>& comm,
  const ArrayRCP<Packet> &recvBuffer,
  const int sourceRank
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::ireceive<int, " << "," << TypeNameTraits<Packet>::name () << ">");
  ValueTypeSerializationBuffer<Ordinal,Packet>
    charRecvBuffer(recvBuffer.size(), recvBuffer.getRawPtr());
  RCP<CommRequest<Ordinal> > commRequest = comm.ireceive(
    charRecvBuffer.getCharBufferView(), sourceRank );
  set_extra_data( recvBuffer, "buffer", inOutArg(commRequest) );
  return commRequest;
}

template<typename Ordinal, typename Packet>
Teuchos::RCP<Teuchos::CommRequest<Ordinal> >
Teuchos::ireceive (const Teuchos::ArrayRCP<Packet> &recvBuffer,
                   const int sourceRank,
                   const int tag,
                   const Teuchos::Comm<Ordinal>& comm)
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::ireceive<int, " << "," << TypeNameTraits<Packet>::name () << ">");
  ValueTypeSerializationBuffer<int, Packet>
    charRecvBuffer (recvBuffer.size (), recvBuffer.getRawPtr ());
  RCP<CommRequest<int> > commRequest =
    comm.ireceive (charRecvBuffer.getCharBufferView (), sourceRank, tag);
  set_extra_data (recvBuffer, "buffer", inOutArg (commRequest));
  return commRequest;
}

template<typename Ordinal, typename Packet>
Teuchos::RCP<Teuchos::CommRequest<Ordinal> >
Teuchos::ireceive(
  const Comm<Ordinal>& comm,
  const RCP<Packet> &recv,
  const int sourceRank
  )
{
  const ArrayRCP<Packet> recvBuffer =
    arcpWithEmbeddedObj( recv.get(), 0, 1, recv, false );
  // 2008/07/29: rabartl: Above: I need to write a helper function to create
  // new ArrayRCP object given a single object to copy.
  return ireceive<Ordinal, Packet>( comm, recvBuffer, sourceRank );
}

template<typename Ordinal, typename Packet, typename Serializer>
Teuchos::RCP<Teuchos::CommRequest<Ordinal> >
Teuchos::ireceive(
  const Comm<Ordinal>& comm,
  const Serializer& serializer,
  const ArrayRCP<Packet> &recvBuffer,
  const int sourceRank
  )
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::CommHelpers: ireceive<"
    <<OrdinalTraits<Ordinal>::name()<<","<<TypeNameTraits<Packet>::name()
    <<">( value type )"
    );
  ValueTypeSerializationBuffer<Ordinal,Packet,Serializer>
    charRecvBuffer(recvBuffer.size(), recvBuffer.getRawPtr(), serializer);
  RCP<CommRequest<Ordinal> > commRequest = comm.ireceive(
    charRecvBuffer.getCharBufferView(), sourceRank );
  set_extra_data( recvBuffer, "buffer", inOutArg(commRequest) );
  return commRequest;
}

template<typename Ordinal>
void Teuchos::waitAll(
  const Comm<Ordinal>& comm,
  const ArrayView<RCP<CommRequest<Ordinal> > > &requests
  )
{
  comm.waitAll(requests);
}


template<typename Ordinal>
void
Teuchos::waitAll (const Comm<Ordinal>& comm,
                  const ArrayView<RCP<CommRequest<Ordinal> > >& requests,
                  const ArrayView<RCP<CommStatus<Ordinal> > >& statuses)
{
  comm.waitAll (requests, statuses);
}


template<typename Ordinal>
Teuchos::RCP<Teuchos::CommStatus<Ordinal> >
Teuchos::wait (const Comm<Ordinal>& comm,
               const Ptr<RCP<CommRequest<Ordinal> > > &request)
{
  return comm.wait (request);
}


#endif // TEUCHOS_COMM_HELPERS_HPP
