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

#ifndef TEUCHOS_MPI_COMM_HPP
#define TEUCHOS_MPI_COMM_HPP


#include "Teuchos_Comm.hpp"
#include "Teuchos_CommUtilities.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_OpaqueWrapper.hpp"
#include "Teuchos_MpiReductionOpSetter.hpp"
#include "Teuchos_SerializationTraitsHelpers.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_Assert.hpp"
#include "mpi.h"


// This must be defined globally for the whole program!
//#define TEUCHOS_MPI_COMM_DUMP


#ifdef TEUCHOS_MPI_COMM_DUMP
#  include "Teuchos_VerboseObject.hpp"
#endif


namespace Teuchos {

//! Human-readable string version of the given MPI error code.
std::string
mpiErrorCodeToString (const int err);

#ifdef TEUCHOS_MPI_COMM_DUMP
template<typename Ordinal, typename T>
void dumpBuffer(
  const std::string &funcName, const std::string &buffName
  ,const Ordinal bytes, const T buff[]
  )
{
  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::OSTab tab(out);
  *out
    << "\n" << funcName << "::" << buffName << ":\n";
  tab.incrTab();
  for( Ordinal i = 0; i < bytes; ++i ) {
    *out << buffName << "[" << i << "] = '" << buff[i] << "'\n";
  }
  *out << "\n";
}
#endif // TEUCHOS_MPI_COMM_DUMP


/// \class MpiCommRequest
/// \brief MPI-specific implementation of \c CommRequest.
///
/// Users would not normally create an instance of this class.  Calls
/// to \c ireceive() and \c isend() return a CommRequest, which for an
/// MPI implementation of \c Comm is an MpiCommRequest.  Users might
/// wish to create an MpiCommRequest directly if they want to
/// encapsulate an MPI_Request returned by an external library or by
/// their own code, and pass it into one of our wrapper functions like
/// \c wait() or \c waitAll().
class MpiCommRequest : public CommRequest {
public:
  //! Constructor (from a raw MPI_Request).
  MpiCommRequest (MPI_Request rawMpiRequest) : 
    rawMpiRequest_ (rawMpiRequest) {}

  /// \brief Return and relinquish ownership of the raw MPI_Request.
  ///
  /// "Relinquish ownership" means that this object sets its raw
  /// MPI_Request to MPI_REQUEST_NULL, but returns the original
  /// MPI_Request.  This effectively gives the caller ownership of the
  /// raw MPI_Request.  This prevents hanging requests.
  MPI_Request releaseRawMpiRequest()
  {
    MPI_Request tmp_rawMpiRequest = rawMpiRequest_;
    rawMpiRequest_ = MPI_REQUEST_NULL;
    return tmp_rawMpiRequest;
  }
private:
  //! The raw request (an opaque object).
  MPI_Request rawMpiRequest_;
  MpiCommRequest(); // Not defined
};

/// \fn mpiCommRequest
/// \brief Nonmember constructor for \c MpiCommRequest.
/// \relates MpiCommRequest
inline RCP<MpiCommRequest>
mpiCommRequest (MPI_Request rawMpiRequest)
{
  return rcp (new MpiCommRequest (rawMpiRequest));
}

/// \class MpiCommStatus
/// \brief MPI-specific implementation of \c CommStatus.
///
/// Users would not normally create an instance of this class.  The
/// only time they might wish to do so is to encapsulate an MPI_Status
/// returned by an external library or by their own code, and pass it
/// into one of our functions like \c wait() or \c waitAll().
///
/// \tparam OrdinalType The same template parameter as \c Comm.  Only
///   use \c int here.  We only make this a template class for
///   compatibility with \c Comm.
template<class OrdinalType>
class MpiCommStatus : public CommStatus<OrdinalType> {
public:
  MpiCommStatus (MPI_Status status) : status_ (status) {}

  //! Destructor (declared virtual for memory safety)
  virtual ~MpiCommStatus() {}

  //! The source rank that sent the message.
  OrdinalType getSourceRank () { return status_.MPI_SOURCE; }

  //! The tag of the received message.
  OrdinalType getTag () { return status_.MPI_TAG; }

  //! The error code of the received message.
  OrdinalType getError () { return status_.MPI_ERROR; }

private:
  //! We forbid default construction syntactically.
  MpiCommStatus ();

  //! The raw MPI_Status struct that this class encapsulates.
  MPI_Status status_;
};

/// \fn mpiCommStatus
/// \brief Nonmember constructor for \c MpiCommStatus.
/// \relates MpiCommStatus
template<class OrdinalType>
inline RCP<MpiCommStatus<OrdinalType> >
mpiCommStatus (MPI_Status rawMpiStatus)
{
  return rcp (new MpiCommStatus<OrdinalType> (rawMpiStatus));
}

/** \brief Concrete communicator subclass based on MPI.
 *
 * <b>Assertions:</b><ul>
 * <li><tt>getRawMpiComm().get()!=NULL && *getRawMpiComm()!=MPI_COMM_NULL</tt>
 * <li><tt>getSize() > 0</tt>
 * <li><tt>0 <= getRank() && getRank() < getSize()</tt>
 * </ul>
 *
 * ToDo: Finish documentation!
 */
template<typename Ordinal>
class MpiComm : public Comm<Ordinal> {
public:

  //! @name Constructors 
  //@{

  /** \brief Construct given a wrapped MPI_Comm oqaque object.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>rawMpiComm.get()!=NULL && *rawMpiComm != MPI_COMM_NULL</tt>
   * </ul>
   */
  MpiComm(
    const RCP<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
    );

  /**
   * \brief Construct a communicator with a new context with the same properties
   * as the original.
   *
   * The newly constructed communicator will have a duplicate communication
   * space that has the same properties (e.g. processes, attributes,
   * topologies) as the input communicator.
   *
   * \param other The communicator to copy from.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>
   * other.getRawMpiComm().get() != NULL && *other.getRawMpiComm() != NULL
   * </tt></li>
   * </ul>
   */
  MpiComm(const MpiComm<Ordinal>& other);

  /** \brief Return the embedded wrapped opaque <tt>MPI_Comm</tt> object. */
  RCP<const OpaqueWrapper<MPI_Comm> > getRawMpiComm() const
  {return rawMpiComm_;}

  /// \brief Set the MPI error handler for this communicator.
  ///
  /// MPI lets you set an error handler function specific to each
  /// communicator.  MpiComm wraps this functionality.  Create an
  /// error handler using \c MPI_Errhandler_create(), or use one of
  /// the default error handlers that the MPI standard or your MPI
  /// implementation provides.  You will need to wrap the MPI error
  /// handler in an OpaqueWrapper first.  (See the documentation of
  /// OpaqueWrapper for the rationale behind not using MPI's opaque
  /// objects directly.)
  ///
  /// MpiComm will not attempt to call MPI_Errhandler_free() on the
  /// error handler you provide.  You can always set the RCP's custom
  /// "deallocator" function to free the error handler, if you want it
  /// taken care of automatically.
  ///
  /// \param errHandler [in] The error handler to set.  If null, do
  ///   nothing.
  void setErrorHandler (const RCP<const OpaqueWrapper<MPI_Errhandler> >& errHandler);

  //@}

  //! @name Overridden from Comm 
  //@{

  /** \brief . */
  virtual int getRank() const;
  /** \brief . */
  virtual int getSize() const;
  /** \brief . */
  virtual void barrier() const;
  /** \brief . */
  virtual void broadcast(
    const int rootRank, const Ordinal bytes, char buffer[]
    ) const;
  /** \brief . */
  virtual void gatherAll(
    const Ordinal sendBytes, const char sendBuffer[]
    ,const Ordinal recvBytes, char recvBuffer[]
    ) const;
  /** \brief . */
  virtual void reduceAll(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char globalReducts[]
    ) const;
  /** \brief . */
  virtual void reduceAllAndScatter(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal sendBytes, const char sendBuffer[]
    ,const Ordinal recvCounts[], char myGlobalReducts[]
    ) const;
  /** \brief . */
	virtual void scan(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char scanReducts[]
    ) const;
  /** \brief . */
  virtual void send(
    const Ordinal bytes, const char sendBuffer[], const int destRank
    ) const;
  /** \brief . */
  virtual int receive(
    const int sourceRank, const Ordinal bytes, char recvBuffer[]
    ) const;
  /** \brief . */
  virtual void readySend(
    const ArrayView<const char> &sendBuffer,
    const int destRank
    ) const;
  /** \brief . */
  virtual RCP<CommRequest> isend(
    const ArrayView<const char> &sendBuffer,
    const int destRank
    ) const;
  /** \brief . */
  virtual RCP<CommRequest> ireceive(
    const ArrayView<char> &Buffer,
    const int sourceRank
    ) const;
  /** \brief . */
  virtual void waitAll(
    const ArrayView<RCP<CommRequest> > &requests
    ) const;
  /** \brief . */
  virtual void 
  waitAll (const ArrayView<RCP<CommRequest> >& requests,
	   const ArrayView<RCP<CommStatus<Ordinal> > >& statuses) const;
  /** \brief . */
  virtual RCP<CommStatus<Ordinal> > 
  wait (const Ptr<RCP<CommRequest> >& request) const;
  /** \brief . */
  virtual RCP< Comm<Ordinal> > duplicate() const;
  /** \brief . */
  virtual RCP< Comm<Ordinal> > split(const int color, const int key) const;
  /** \brief . */
  virtual RCP< Comm<Ordinal> > createSubcommunicator(
    const ArrayView<const int>& ranks) const;
  //@}

  //! @name Overridden from Describable 
  //@{

  /** \brief . */
  std::string description() const;

  //@}

  // These should be private but the PGI compiler requires them be public

  static int const minTag_ = 26000; // These came from Teuchos::MpiComm???
  static int const maxTag_ = 26099; // ""

private:

  // Set internal data members once the rawMpiComm_ data member is valid.
  void setupMembersFromComm();
  static int tagCounter_;

  RCP<const OpaqueWrapper<MPI_Comm> > rawMpiComm_;
  int rank_;
  int size_;
  int tag_;

  //! MPI error handler.  If null, MPI uses the default error handler.
  RCP<const OpaqueWrapper<MPI_Errhandler> > customErrorHandler_;

  void assertRank(const int rank, const std::string &rankName) const;

  // Not defined and not to be called!
  MpiComm();

#ifdef TEUCHOS_MPI_COMM_DUMP
public:
  static bool show_dump;
#endif // TEUCHOS_MPI_COMM_DUMP
	
};


/** \brief Helper function that creates a dynamically allocated
 * <tt>MpiComm</tt> object or returns <tt>Teuchos::null</tt> to correctly
 * represent a null communicator.
 *
 * <b>Postconditions:</b></ul>
 * <li>[<tt>rawMpiComm.get()!=NULL && *rawMpiComm!=MPI_COMM_NULL</tt>]
 *     <tt>return.get()!=NULL</tt>
 * <li>[<tt>rawMpiComm.get()==NULL || *rawMpiComm==MPI_COMM_NULL</tt>]
 *     <tt>return.get()==NULL</tt>
 * </ul>
 *
 * \relates MpiComm
 */
template<typename Ordinal>
RCP<MpiComm<Ordinal> >
createMpiComm(
  const RCP<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
  );


// ////////////////////////
// Implementations


// Static members


template<typename Ordinal>
int MpiComm<Ordinal>::tagCounter_ = MpiComm<Ordinal>::minTag_;


// Constructors


template<typename Ordinal>
MpiComm<Ordinal>::MpiComm(
  const RCP<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
  )
{
  TEUCHOS_TEST_FOR_EXCEPT( rawMpiComm.get()==NULL );
  TEUCHOS_TEST_FOR_EXCEPT( *rawMpiComm == MPI_COMM_NULL );
  rawMpiComm_ = rawMpiComm;
  setupMembersFromComm();
}


template<typename Ordinal>
MpiComm<Ordinal>::MpiComm(const MpiComm<Ordinal>& other)
{
  TEUCHOS_TEST_FOR_EXCEPT(other.getRawMpiComm().get() == NULL);
  TEUCHOS_TEST_FOR_EXCEPT(*other.getRawMpiComm() == MPI_COMM_NULL);
  MPI_Comm newComm;
  MPI_Comm_dup(*other.getRawMpiComm(), &newComm);
  rawMpiComm_ = opaqueWrapper(newComm);
  setupMembersFromComm();
}


template<typename Ordinal>
void MpiComm<Ordinal>::setupMembersFromComm()
{
  MPI_Comm_size(*rawMpiComm_, &size_);
  MPI_Comm_rank(*rawMpiComm_, &rank_);
  if(tagCounter_ > maxTag_)
    tagCounter_ = minTag_;
  tag_ = tagCounter_++;
}


template<typename Ordinal>
void 
MpiComm<Ordinal>::
setErrorHandler (const RCP<const OpaqueWrapper<MPI_Errhandler> >& errHandler)
{
  if (! is_null (errHandler)) {
    const int err = MPI_Comm_set_errhandler (*getRawMpiComm(), *errHandler);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
      "Teuchos::MpiComm::setErrorHandler: MPI_Comm_set_errhandler() failed with "
      "error \"" << mpiErrorCodeToString (err) << "\".");
  }
  // Wait to set this until the end, in case MPI_Errhandler_set()
  // doesn't succeed.
  customErrorHandler_ = errHandler;
}



// Overridden from Comm

  
template<typename Ordinal>
int MpiComm<Ordinal>::getRank() const
{
  return rank_;
}

  
template<typename Ordinal>
int MpiComm<Ordinal>::getSize() const
{
  return size_;
}

  
template<typename Ordinal>
void MpiComm<Ordinal>::barrier() const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::barrier()"
    );
  MPI_Barrier(*rawMpiComm_);
}

  
template<typename Ordinal>
void MpiComm<Ordinal>::broadcast(
  const int rootRank, const Ordinal bytes, char buffer[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::broadcast(...)"
    );
  MPI_Bcast(buffer,bytes,MPI_CHAR,rootRank,*rawMpiComm_);
}

  
template<typename Ordinal>
void MpiComm<Ordinal>::gatherAll(
  const Ordinal sendBytes, const char sendBuffer[],
  const Ordinal recvBytes, char recvBuffer[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::gatherAll(...)"
    );
  TEUCHOS_ASSERT_EQUALITY((sendBytes*size_), recvBytes );
  MPI_Allgather(
    const_cast<char *>(sendBuffer), sendBytes, MPI_CHAR,
    recvBuffer, sendBytes, MPI_CHAR,
    *rawMpiComm_
    );
  // NOTE: 'sendBytes' is being sent above for the MPI arg recvcount (which is
  // very confusing in the MPI documentation) for MPI_Allgether(...).
}

  
template<typename Ordinal>
void MpiComm<Ordinal>::reduceAll(
  const ValueTypeReductionOp<Ordinal,char> &reductOp
  ,const Ordinal bytes, const char sendBuffer[], char globalReducts[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::reduceAll(...)"
    );
  MpiReductionOpSetter op(mpiReductionOp(rcp(&reductOp,false)));
  MPI_Datatype char_block;
  MPI_Type_contiguous(bytes, MPI_CHAR, &char_block);
  MPI_Type_commit(&char_block);
  MPI_Allreduce(
    const_cast<char*>(sendBuffer),globalReducts,1,char_block,op.mpi_op()
    ,*rawMpiComm_
    );
  MPI_Type_free(&char_block);
}


template<typename Ordinal>
void MpiComm<Ordinal>::reduceAllAndScatter(
  const ValueTypeReductionOp<Ordinal,char> &reductOp
  ,const Ordinal sendBytes, const char sendBuffer[]
  ,const Ordinal recvCounts[], char myGlobalReducts[]
  ) const
{

  (void)sendBytes; // Ignore if not in debug mode

  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::reduceAllAndScatter(...)"
    );

#ifdef TEUCHOS_DEBUG
  Ordinal sumRecvBytes = 0;
  for( Ordinal i = 0; i < size_; ++i ) {
    sumRecvBytes += recvCounts[i];
  }
  TEUCHOS_TEST_FOR_EXCEPT(!(sumRecvBytes==sendBytes));
#endif // TEUCHOS_DEBUG

#ifdef TEUCHOS_MPI_COMM_DUMP
  if(show_dump) {
    dumpBuffer<Ordinal,char>(
      "Teuchos::MpiComm<Ordinal>::reduceAllAndScatter(...)",
      "sendBuffer", sendBytes, sendBuffer );
    dumpBuffer<Ordinal,Ordinal>(
      "Teuchos::MpiComm<Ordinal>::reduceAllAndScatter(...)",
      "recvCounts", as<Ordinal>(size_), recvCounts );
    dumpBuffer<Ordinal,char>(
      "Teuchos::MpiComm<Ordinal>::reduceAllAndScatter(...)",
      "myGlobalReducts", as<char>(recvCounts[rank_]), myGlobalReducts );
  }
#endif // TEUCHOS_MPI_COMM_DUMP

  // Create a new recvCount[] if Ordinal!=int
  WorkspaceStore* wss = get_default_workspace_store().get();
  const bool Ordinal_is_int = typeid(int)==typeid(Ordinal);
  Workspace<int> ws_int_recvCounts(wss,Ordinal_is_int?0:size_);
  const int *int_recvCounts = 0;
  if(Ordinal_is_int) {
    int_recvCounts = reinterpret_cast<const int*>(recvCounts);
    // Note: We must do an reinterpet cast since this must
    // compile even if it is not executed.  I could implement
    // code that would not need to do this using template
    // conditionals but I don't want to bother.
  }
  else {
    std::copy(recvCounts, recvCounts+size_, &ws_int_recvCounts[0]);
    int_recvCounts = &ws_int_recvCounts[0];
  }

  // Perform the operation
  MpiReductionOpSetter op(mpiReductionOp(rcp(&reductOp, false)));
  MPI_Reduce_scatter(
    const_cast<char*>(sendBuffer), myGlobalReducts,
    const_cast<int*>(int_recvCounts),
    MPI_CHAR,
    op.mpi_op(),
    *rawMpiComm_
    );

}


template<typename Ordinal>
void MpiComm<Ordinal>::scan(
  const ValueTypeReductionOp<Ordinal,char> &reductOp
  ,const Ordinal bytes, const char sendBuffer[], char scanReducts[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::scan(...)" );

  MpiReductionOpSetter op(mpiReductionOp(rcp(&reductOp,false)));
  const int err = 
    MPI_Scan (const_cast<char*>(sendBuffer), scanReducts, bytes, MPI_CHAR, 
	      op.mpi_op(), *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::scan: MPI_Scan() failed with error code \"" 
    << mpiErrorCodeToString (err) << "\".");
}


template<typename Ordinal>
void 
MpiComm<Ordinal>::send (const Ordinal bytes, 
			const char sendBuffer[], 
			const int destRank) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::send(...)" );

#ifdef TEUCHOS_MPI_COMM_DUMP
  if(show_dump) {
    dumpBuffer<Ordinal,char>(
      "Teuchos::MpiComm<Ordinal>::send(...)"
      ,"sendBuffer", bytes, sendBuffer
      );
  }
#endif // TEUCHOS_MPI_COMM_DUMP

  const int err = MPI_Send (const_cast<char*>(sendBuffer), bytes, MPI_CHAR,
			    destRank, tag_, *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::send: MPI_Send() failed with error code \"" 
    << mpiErrorCodeToString (err) << "\".");
}


template<typename Ordinal>
void MpiComm<Ordinal>::readySend(
  const ArrayView<const char> &sendBuffer,
  const int destRank
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::readySend" );

#ifdef TEUCHOS_MPI_COMM_DUMP
  if(show_dump) {
    dumpBuffer<Ordinal,char>(
      "Teuchos::MpiComm<Ordinal>::readySend(...)"
      ,"sendBuffer", bytes, sendBuffer
      );
  }
#endif // TEUCHOS_MPI_COMM_DUMP

  const int err = 
    MPI_Rsend (const_cast<char*>(sendBuffer.getRawPtr()), sendBuffer.size(),
	       MPI_CHAR, destRank, tag_, *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::readySend: MPI_Rsend() failed with error code \"" 
    << mpiErrorCodeToString (err) << "\".");
}


template<typename Ordinal>
int 
MpiComm<Ordinal>::receive (const int sourceRank, 
			   const Ordinal bytes, 
			   char recvBuffer[]) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::receive(...)" );

  // A negative source rank indicates MPI_ANY_SOURCE, namely that we
  // will take an incoming message from any process, as long as the
  // tag matches.
  const int theSrcRank = (sourceRank < 0) ? MPI_ANY_SOURCE : sourceRank;

  MPI_Status status;
  const int err = MPI_Recv (recvBuffer, bytes, MPI_CHAR, theSrcRank, tag_, 
			    *rawMpiComm_, &status);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::receive: MPI_Recv() failed with error code \"" 
    << mpiErrorCodeToString (err) << "\".");

#ifdef TEUCHOS_MPI_COMM_DUMP
  if (show_dump) {
    dumpBuffer<Ordinal,char> ("Teuchos::MpiComm<Ordinal>::receive(...)", 
			      "recvBuffer", bytes, recvBuffer);
  }
#endif // TEUCHOS_MPI_COMM_DUMP

  // Returning the source rank is useful in the MPI_ANY_SOURCE case.
  return status.MPI_SOURCE;
}


template<typename Ordinal>
RCP<CommRequest> 
MpiComm<Ordinal>::isend (const ArrayView<const char> &sendBuffer,
			 const int destRank) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::isend(...)" );

  MPI_Request rawMpiRequest = MPI_REQUEST_NULL;
  const int err = 
    MPI_Isend (const_cast<char*>(sendBuffer.getRawPtr()), sendBuffer.size(), 
	       MPI_CHAR, destRank, tag_, *rawMpiComm_, &rawMpiRequest);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::isend: MPI_Isend() failed with error code \"" 
    << mpiErrorCodeToString (err) << "\".");

  return mpiCommRequest (rawMpiRequest);
}


template<typename Ordinal>
RCP<CommRequest> 
MpiComm<Ordinal>::ireceive (const ArrayView<char> &recvBuffer,
			    const int sourceRank) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::ireceive(...)" );

  // A negative source rank indicates MPI_ANY_SOURCE, namely that we
  // will take an incoming message from any process, as long as the
  // tag matches.
  const int theSrcRank = (sourceRank < 0) ? MPI_ANY_SOURCE : sourceRank;

  MPI_Request rawMpiRequest = MPI_REQUEST_NULL;
  const int err = 
    MPI_Irecv (const_cast<char*>(recvBuffer.getRawPtr()), recvBuffer.size(), 
	       MPI_CHAR, theSrcRank, tag_, *rawMpiComm_, &rawMpiRequest);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::ireceive: MPI_Irecv() failed with error code \"" 
    << mpiErrorCodeToString (err) << "\".");

  return mpiCommRequest (rawMpiRequest);
}


template<typename Ordinal>
void 
MpiComm<Ordinal>::
waitAll (const ArrayView<RCP<CommRequest> >& requests) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::waitAll(...)" );

  const int count = requests.size();
  if (count == 0) {
    return; // No requests on which to wait
  }
  
  // Pull out the raw MPI requests from the wrapped requests.
  Array<MPI_Request> rawMpiRequests (count, MPI_REQUEST_NULL);
  for (int i = 0; i < count; ++i) {
    RCP<CommRequest> request = requests[i];
    if (! is_null (request)) {
      RCP<MpiCommRequest> mpiCommRequest = 
	rcp_dynamic_cast<MpiCommRequest> (request);
      // This call makes this function wait() not satisfy the strong
      // exception guarantee, because releaseRawMpiRequest() modifies
      // MpiCommRequest.
      rawMpiRequests[i] = mpiCommRequest->releaseRawMpiRequest();
    }
    else { // Null requests map to MPI_REQUEST_NULL
      rawMpiRequests[i] = MPI_REQUEST_NULL;
    }
    requests[i] = null; // We've invalidated it above, so set to null.
  }

  Array<MPI_Status> rawMpiStatuses (count);
  const int err = MPI_Waitall (count, rawMpiRequests.getRawPtr(), 
			       rawMpiStatuses.getRawPtr());
  // The MPI standard doesn't say what happens to the rest of the
  // requests if one of them failed, in a multiple completion routine
  // like MPI_Waitall().  We conservatively abort in that case.  If
  // err == MPI_ERR_IN_STATUS, then we need to check each status'
  // MPI_ERROR field to find the error.
  if (err != MPI_SUCCESS) {
    if (err == MPI_ERR_IN_STATUS) {
      int firstErr = MPI_SUCCESS;
      int firstIndexFailed = 0;
      for (int i = 0; i < count; ++i) {
	if (rawMpiStatuses[i].MPI_ERROR != MPI_SUCCESS) {
	  firstErr = rawMpiStatuses[i].MPI_ERROR;
	  firstIndexFailed = i;
	  break;
	}
      }
      TEUCHOS_TEST_FOR_EXCEPTION(firstErr != MPI_SUCCESS, std::runtime_error, 
        "Teuchos::MpiComm::waitAll: MPI_Waitall() failed with error code "
        "\"MPI_ERR_IN_STATUS\".  Of the " << count << " request" 
        << (count != 1 ? "s" : "") << " given to MPI_Waitall(), the smallest "
        "0-based index that failed is " << firstIndexFailed << " and its error "
        "code is \"" << mpiErrorCodeToString (firstErr) << "\".");
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, 
        "Teuchos::MpiComm::waitAll: MPI_Waitall() failed with error code "
	<< mpiErrorCodeToString (err) << "\".");
    }
  }
}


template<typename Ordinal>
void 
MpiComm<Ordinal>::
waitAll (const ArrayView<RCP<CommRequest> >& requests,
	 const ArrayView<RCP<CommStatus<Ordinal> > >& statuses) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::waitAll(requests, statuses)" );

  const int count = requests.size();
  TEUCHOS_TEST_FOR_EXCEPTION(count != statuses.size(), 
    std::invalid_argument, "Teuchos::MpiComm::waitAll: requests.size() = " 
    << requests.size() << " != statuses.size() = " << statuses.size() << ".");
  if (count == 0) {
    return;
  }

  // MpiComm wraps MPI and can't expose any MPI structs or opaque
  // objects.  Thus, we have to unpack both requests and statuses into
  // separate arrays.  If that's too slow, then your code should just
  // call into MPI directly.
  Array<MPI_Request> rawMpiRequests (count, MPI_REQUEST_NULL);
  for (int i = 0; i < count; ++i) {
    if (! requests[i].is_null()) {
      RCP<MpiCommRequest> mpiRequest = 
	rcp_dynamic_cast<MpiCommRequest> (requests[i]);
      // This convinces the CommRequest to set its own raw MPI_Request
      // to MPI_REQUEST_NULL, so we won't have any hanging raw request
      // handles floating around afterwards.
      //
      // This call makes this function wait() not satisfy the strong
      // exception guarantee, because releaseRawMpiRequest() modifies
      // MpiCommRequest.
      rawMpiRequests[i] = mpiRequest->releaseRawMpiRequest ();
      requests[i] = null;
    }
  }

  Array<MPI_Status> rawMpiStatuses (count);
  // This is the part where we've finally peeled off the wrapper and
  // we can now interact with MPI directly.
  const int err = MPI_Waitall (count, rawMpiRequests.getRawPtr(), 
			       rawMpiStatuses.getRawPtr());

  // The MPI standard doesn't say what happens to the rest of the
  // requests if one of them failed, in a multiple completion routine
  // like MPI_Waitall().  We conservatively abort in that case.  If
  // err == MPI_ERR_IN_STATUS, then we need to check each status'
  // MPI_ERROR field to find the error.
  if (err != MPI_SUCCESS) {
    if (err == MPI_ERR_IN_STATUS) {
      int firstErr = MPI_SUCCESS;
      int firstIndexFailed = 0;
      for (int i = 0; i < count; ++i) {
	if (rawMpiStatuses[i].MPI_ERROR != MPI_SUCCESS) {
	  firstErr = rawMpiStatuses[i].MPI_ERROR;
	  firstIndexFailed = i;
	  break;
	}
      }
      TEUCHOS_TEST_FOR_EXCEPTION(firstErr != MPI_SUCCESS, std::runtime_error, 
        "Teuchos::MpiComm::waitAll: MPI_Waitall() failed with error code "
        "\"MPI_ERR_IN_STATUS\".  Of the " << count << " request" 
        << (count != 1 ? "s" : "") << " given to MPI_Waitall(), the smallest "
        "0-based index that failed is " << firstIndexFailed << " and its error "
        "code is \"" << mpiErrorCodeToString (firstErr) << "\".");
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, 
        "Teuchos::MpiComm::waitAll: MPI_Waitall() failed with error code "
	<< mpiErrorCodeToString (err) << "\".");
    }
  }

  // Repackage the raw MPI_Status structs into the CommStatus wrappers.
  for (int i = 0; i < count; ++i) {
    statuses[i] = mpiCommStatus<Ordinal> (rawMpiStatuses[i]);
  }
}


template<typename Ordinal>
RCP<CommStatus<Ordinal> > 
MpiComm<Ordinal>::wait (const Ptr<RCP<CommRequest> >& request) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::wait(...)" );

  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::wait(...)" );

  if (is_null(*request)) {
    return null; // Nothing to wait on ...
  }
  const RCP<MpiCommRequest> mpiCommRequest =
    rcp_dynamic_cast<MpiCommRequest>(*request);
  // This function doesn't satisfy the strong exception guarantee,
  // because releaseRawMpiRequest() modifies the MpiCommRequest.
  MPI_Request rawMpiRequest = mpiCommRequest->releaseRawMpiRequest();
  MPI_Status status;
  const int err = MPI_Wait (&rawMpiRequest, &status);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, 
    "Teuchos::MpiComm::wait: MPI_Wait() failed with error code "
    << mpiErrorCodeToString (err) << "\".");

  *request = null;
  return rcp (new MpiCommStatus<Ordinal> (status));
}


template<typename Ordinal>
RCP< Comm<Ordinal> >
MpiComm<Ordinal>::duplicate() const
{
  return rcp (new MpiComm<Ordinal> (*this));
}


template<typename Ordinal>
RCP< Comm<Ordinal> >
MpiComm<Ordinal>::split(const int color, const int key) const
{
  MPI_Comm newComm;
  int splitReturn = MPI_Comm_split(
    *rawMpiComm_,
    color < 0 ? MPI_UNDEFINED : color,
    key,
    &newComm);
  TEUCHOS_TEST_FOR_EXCEPTION(
    splitReturn != MPI_SUCCESS,
    std::logic_error,
    "Failed to create communicator with color " << color <<
    "and key " << key << ".");
  if (newComm == MPI_COMM_NULL) {
    return RCP< Comm<Ordinal> >();
  } else {
    return rcp(new MpiComm<Ordinal>(
      rcp_implicit_cast<const OpaqueWrapper<MPI_Comm> >( opaqueWrapper(newComm) )
      ) );
  }
}


template<typename Ordinal>
RCP< Comm<Ordinal> >
MpiComm<Ordinal>::createSubcommunicator(const ArrayView<const int> &ranks) const
{
  int mpiReturn;

  // Get the group that this communicator is in.
  MPI_Group thisGroup;
  mpiReturn = MPI_Comm_group(*rawMpiComm_, &thisGroup);
  TEUCHOS_TEST_FOR_EXCEPTION(mpiReturn != MPI_SUCCESS, std::logic_error,
                     "Failed to obtain group.");
  // Create a new group with the specified members.
  MPI_Group newGroup;
  mpiReturn = MPI_Group_incl(
    thisGroup, ranks.size(), const_cast<int *>(&ranks[0]), &newGroup);
  TEUCHOS_TEST_FOR_EXCEPTION(mpiReturn != MPI_SUCCESS, std::logic_error,
                     "Failed to create subgroup.");
  // Create a new communicator from the new group.
  MPI_Comm newComm;
  mpiReturn = MPI_Comm_create(*rawMpiComm_, newGroup, &newComm);
  TEUCHOS_TEST_FOR_EXCEPTION(mpiReturn != MPI_SUCCESS, std::logic_error,
                     "Failed to create subcommunicator.");
  if (newComm == MPI_COMM_NULL) {
    return RCP< Comm<Ordinal> >();
  } else {
    return rcp(new MpiComm<Ordinal>(
      rcp_implicit_cast<const OpaqueWrapper<MPI_Comm> >( opaqueWrapper(newComm) )
      ));
  }
}


// Overridden from Describable


template<typename Ordinal>
std::string MpiComm<Ordinal>::description() const
{
  std::ostringstream oss;
  oss
    << typeName(*this)
    << "{"
    << "size="<<size_
    << ",rank="<<rank_
    << ",rawMpiComm="<<static_cast<MPI_Comm>(*rawMpiComm_)
    <<"}";
  return oss.str();
}


#ifdef TEUCHOS_MPI_COMM_DUMP
template<typename Ordinal>
bool MpiComm<Ordinal>::show_dump = false;
#endif


// private


template<typename Ordinal>
void MpiComm<Ordinal>::assertRank(const int rank, const std::string &rankName) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! ( 0 <= rank && rank < size_ ), std::logic_error
    ,"Error, "<<rankName<<" = " << rank << " is not < 0 or is not"
    " in the range [0,"<<size_-1<<"]!"
    );
}


} // namespace Teuchos


template<typename Ordinal>
Teuchos::RCP<Teuchos::MpiComm<Ordinal> >
Teuchos::createMpiComm(
  const RCP<const OpaqueWrapper<MPI_Comm> > &rawMpiComm
  )
{
  if( rawMpiComm.get()!=NULL && *rawMpiComm != MPI_COMM_NULL )
    return rcp(new MpiComm<Ordinal>(rawMpiComm));
  return Teuchos::null;
}


#endif // TEUCHOS_MPI_COMM_HPP
