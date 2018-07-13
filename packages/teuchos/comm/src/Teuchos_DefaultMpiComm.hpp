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

/// \file Teuchos_DefaultMpiComm.hpp
/// \brief Implementation of Teuchos wrappers for MPI.
///
/// \warning It only makes sense to include this file if MPI is enabled.

#include <Teuchos_ConfigDefs.hpp>

// If MPI is not enabled, disable the contents of this file.
#ifdef HAVE_TEUCHOS_MPI

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
#include <mpi.h>
#include <iterator>

// This must be defined globally for the whole program!
//#define TEUCHOS_MPI_COMM_DUMP

#ifdef TEUCHOS_MPI_COMM_DUMP
#  include "Teuchos_VerboseObject.hpp"
#endif

namespace Teuchos {

//! Human-readable string version of the given MPI error code.
std::string
mpiErrorCodeToString (const int err);

namespace details {
  /// \brief Give \c comm to MPI_Comm_free, if MPI is not yet finalized.
  ///
  /// This function "frees" the given communicator by giving it to
  /// MPI_Comm_free.  It only does so if MPI_Finalize has not yet been
  /// called.  If MPI_Finalize has been called, this function does
  /// nothing, since it is not legal to call most MPI functions after
  /// MPI_Finalize has been called.  This function also ignores any
  /// errors returned by MPI_Finalize, which makes it suitable for use
  /// in a destructor.
  ///
  /// \note This function may allow a memory leak in your program, if
  ///   you have allowed the MPI_Comm to persist after MPI_Finalize
  ///   has been called.
  void safeCommFree (MPI_Comm* comm);

  /// Set the given communicator's error handler to \c handler.
  ///
  /// If the MPI version is >= 2, this calls MPI_Comm_set_handler().
  /// If the MPI version is 1, this calls MPI_Errhandler_set().
  int setCommErrhandler (MPI_Comm comm, MPI_Errhandler handler);

} // namespace details

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

/// \class MpiCommStatus
/// \brief MPI-specific implementation of CommStatus.
///
/// Users would not normally create an instance of this class.  The
/// only time they might wish to do so is to encapsulate an MPI_Status
/// returned by an external library or by their own code, and pass it
/// into one of our functions like wait() or waitAll().
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
/// \brief Nonmember constructor for MpiCommStatus.
/// \relates MpiCommStatus
template<class OrdinalType>
inline RCP<MpiCommStatus<OrdinalType> >
mpiCommStatus (MPI_Status rawMpiStatus)
{
  return rcp (new MpiCommStatus<OrdinalType> (rawMpiStatus));
}

/// \class MpiCommRequestBase
/// \brief Base class MPI implementation of CommRequest.
/// \tparam OrdinalType Same as the template parameter of Comm.
///
/// This class wraps MPI_Request, which is MPI's reification of a
/// nonblocking communication operation.
///
/// Users would not normally create an instance of this class.  Calls
/// to nonblocking communication operations (such as ireceive() or
/// isend()) return a pointer to a CommRequest.  If the Comm is an
/// MpiComm, then the returned CommRequest is an MpiCommRequest.
///
/// Users might wish to create an MpiCommRequest directly if they want
/// to encapsulate an MPI_Request returned by an external library or
/// by their own code.
template<class OrdinalType>
class MpiCommRequestBase : public CommRequest<OrdinalType> {
public:
  //! Default constructor.
  MpiCommRequestBase () :
    rawMpiRequest_ (MPI_REQUEST_NULL)
  {}

  //! Constructor (from a raw MPI_Request).
  MpiCommRequestBase (MPI_Request rawMpiRequest) :
    rawMpiRequest_ (rawMpiRequest)
  {}

  /// \brief Return and relinquish ownership of the raw MPI_Request.
  ///
  /// "Relinquish ownership" means that this object sets its raw
  /// MPI_Request to <tt>MPI_REQUEST_NULL</tt>, but returns the
  /// original MPI_Request.  This effectively gives the caller
  /// ownership of the raw MPI_Request.  This prevents hanging
  /// requests.
  MPI_Request releaseRawMpiRequest()
  {
    MPI_Request tmp_rawMpiRequest = rawMpiRequest_;
    rawMpiRequest_ = MPI_REQUEST_NULL;
    return tmp_rawMpiRequest;
  }

  //! Whether the raw MPI_Request is <tt>MPI_REQUEST_NULL</tt>.
  bool isNull() const {
    return rawMpiRequest_ == MPI_REQUEST_NULL;
  }

  /// \brief Wait on this communication request to complete.
  ///
  /// This is a blocking operation.  The user is responsible for
  /// avoiding deadlock.  (For example, a receive must have a matching
  /// send, otherwise a wait on the receive will wait forever.)
  RCP<CommStatus<OrdinalType> > wait () {
    MPI_Status rawMpiStatus;
    // Whether this function satisfies the strong exception guarantee
    // depends on whether MPI_Wait modifies its input request on error.
    const int err = MPI_Wait (&rawMpiRequest_, &rawMpiStatus);
    TEUCHOS_TEST_FOR_EXCEPTION(
      err != MPI_SUCCESS, std::runtime_error,
      "Teuchos: MPI_Wait() failed with error \""
      << mpiErrorCodeToString (err));
    // MPI_Wait sets the MPI_Request to MPI_REQUEST_NULL on success.
    return mpiCommStatus<OrdinalType> (rawMpiStatus);
  }

  /// \brief Cancel the communication request, and return its status.
  ///
  /// If this request is invalid or has already been invalidated, this
  /// method returns null.
  RCP<CommStatus<OrdinalType> > cancel () {
    if (rawMpiRequest_ == MPI_REQUEST_NULL) {
      return null;
    }
    else {
      int err = MPI_Cancel (&rawMpiRequest_);
      TEUCHOS_TEST_FOR_EXCEPTION(
        err != MPI_SUCCESS, std::runtime_error,
        "Teuchos: MPI_Cancel failed with the following error: "
        << mpiErrorCodeToString (err));

      // Wait on the request.  If successful, MPI_Wait will set the
      // MPI_Request to MPI_REQUEST_NULL.  The returned status may
      // still be useful; for example, one may call MPI_Test_cancelled
      // to test an MPI_Status from a nonblocking send.
      MPI_Status status;
      err = MPI_Wait (&rawMpiRequest_, &status);
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
        "Teuchos::MpiCommStatus::cancel: MPI_Wait failed with the following "
        "error: " << mpiErrorCodeToString (err));
      return mpiCommStatus<OrdinalType> (status);
    }
  }

  //! Destructor; cancels the request if it is still pending.
  virtual ~MpiCommRequestBase () {
    if (rawMpiRequest_ != MPI_REQUEST_NULL) {
      // We're in a destructor, so don't throw errors.  However, if
      // MPI_Cancel fails, it's probably a bad idea to call MPI_Wait.
      const int err = MPI_Cancel (&rawMpiRequest_);
      if (err == MPI_SUCCESS) {
        // The MPI_Cancel succeeded.  Now wait on the request.  Ignore
        // any reported error, since we can't do anything about those
        // in the destructor (other than kill the program).  If
        // successful, MPI_Wait will set the MPI_Request to
        // MPI_REQUEST_NULL.  We ignore the returned MPI_Status, since
        // if the user let the request fall out of scope, she must not
        // care about the status.
        //
        // mfh 21 Oct 2012: The MPI standard requires completing a
        // canceled request by calling a function like MPI_Wait,
        // MPI_Test, or MPI_Request_free.  MPI_Wait on a canceled
        // request behaves like a local operation (it does not
        // communicate or block waiting for communication).  One could
        // also call MPI_Request_free instead of MPI_Wait, but
        // MPI_Request_free is intended more for persistent requests
        // (created with functions like MPI_Recv_init).
        (void) MPI_Wait (&rawMpiRequest_, MPI_STATUS_IGNORE);
      }
    }
  }

private:
  //! The raw MPI request (an opaque object).
  MPI_Request rawMpiRequest_;
};

/// \class MpiCommRequest
/// \brief MPI implementation of CommRequest.
/// \tparam OrdinalType Same as the template parameter of Comm.
///
/// This class wraps MPI_Request, which is MPI's reification of a
/// nonblocking communication operation.
///
/// Users would not normally create an instance of this class.  Calls
/// to nonblocking communication operations (such as \c ireceive() or
/// \c isend()) return a pointer to a CommRequest.  If the Comm is an
/// MpiComm, then the returned CommRequest is an MpiCommRequest.
///
/// Users might wish to create an MpiCommRequest directly if they want
/// to encapsulate an MPI_Request returned by an external library or
/// by their own code.
template<class OrdinalType>
class MpiCommRequest : public MpiCommRequestBase<OrdinalType> {
public:
  //! Default constructor.
  MpiCommRequest () :
    MpiCommRequestBase<OrdinalType> (MPI_REQUEST_NULL),
    numBytes_ (0)
  {}

  //! Constructor (from a raw MPI_Request).
  MpiCommRequest (MPI_Request rawMpiRequest,
                  const ArrayView<char>::size_type numBytesInMessage) :
    MpiCommRequestBase<OrdinalType> (rawMpiRequest),
    numBytes_ (numBytesInMessage)
  {}

  /// \brief Number of bytes in the nonblocking send or receive request.
  ///
  /// Remembering this is inexpensive, and is also useful for
  /// debugging (e.g., for detecting whether the send and receive have
  /// matching message lengths).
  ArrayView<char>::size_type numBytes () const {
    return numBytes_;
  }

  //! Destructor; cancels the request if it is still pending.
  virtual ~MpiCommRequest () {}

private:
  //! Number of bytes in the nonblocking send or receive request.
  ArrayView<char>::size_type numBytes_;
};

/// \fn mpiCommRequest
/// \brief Nonmember constructor for MpiCommRequest.
/// \tparam OrdinalType Same as the template parameter of MpiCommRequest.
/// \relates MpiCommRequest
///
/// \param rawMpiRequest [in] The raw MPI_Request opaque object.
/// \param numBytes [in] The number of bytes in the nonblocking
///   send or receive request.
template<class OrdinalType>
inline RCP<MpiCommRequest<OrdinalType> >
mpiCommRequest (MPI_Request rawMpiRequest,
                const ArrayView<char>::size_type numBytes)
{
  return rcp (new MpiCommRequest<OrdinalType> (rawMpiRequest, numBytes));
}

/// \class MpiComm
/// \brief Implementation of Comm that uses MPI for communication.
/// \tparam Ordinal The index type for communication; same as the
///   template parameter of Comm.
///
/// This class uses MPI (the Message Passing Interface) to implement
/// the Comm interface.  It includes constructors that take an
/// MPI_Comm from the application.
///
/// Assertions:
/// - <tt>getRawMpiComm().get() != NULL<tt>
/// - <tt>*getRawMpiComm() != MPI_COMM_NULL</tt>
/// - <tt>getSize() > 0</tt>
/// - <tt>0 <= getRank() && getRank() < getSize()</tt>
///
template<typename Ordinal>
class MpiComm : public Comm<Ordinal> {
public:
  //! @name Constructors
  //@{

  /// \brief Construct an MpiComm with an MPI_Comm.
  ///
  /// This constructs an MpiComm that uses the given "raw" MPI
  /// communicator underneath.  The MPI_Comm must be valid for the
  /// lifetime of this MpiComm.  You are responsible for freeing the
  /// MPI_Comm (using MPI_Comm_free) if necessary.
  ///
  /// This constructor should be used only in one of two cases:
  /// 1. When the given MPI_Comm is one of the predefined
  ///    communicators that need not and must not be freed after use,
  ///    like MPI_COMM_WORLD or MPI_COMM_SELF.
  /// 2. When you know that the given MPI_Comm will not be freed until
  ///    after the lifetime of this MpiComm.
  ///
  /// If you cannot be sure of either of these two conditions, you
  /// should use the version of the constructor that takes an
  /// <tt>RCP<const OpaqueWrapper<MPI_Comm> ></tt>.
  ///
  /// Precondition:
  /// - <tt>rawMpiComm != MPI_COMM_NULL</tt>
  explicit MpiComm (MPI_Comm rawMpiComm);

  /// \brief Construct an MpiComm with the given wrapped MPI_Comm.
  ///
  /// This constructs an MpiComm that uses the given "raw" MPI
  /// communicator underneath.  This version of the constructor
  /// accepts the MPI_Comm wrapped in an OpaqueWrapper, which along
  /// with the RCP has the option to free the MPI_Comm (via
  /// MPI_Comm_free) after use if necessary.  You are responsible when
  /// creating the OpaqueWrapper for supplying a "free" function if
  /// needed.  We recommend using details::safeCommFree for the "free"
  /// function, if one is needed.
  ///
  /// Preconditions:
  /// - <tt>rawMpiComm.get() != NULL</tt>
  /// - <tt>*rawMpiComm != MPI_COMM_NULL</tt>
  MpiComm (const RCP<const OpaqueWrapper<MPI_Comm> >& rawMpiComm);

  /// \brief Construct an MpiComm with a wrapped MPI_Comm and a default tag.
  ///
  /// This constructor has the same meaning as the one-argument
  /// constructor that takes RCP<const OpaqueWrapper<MPI_Comm> >,
  /// except that it sets the default message tag on all processes to
  /// \c defaultTag.  This avoids the MPI_Bcast that the other two
  /// constructors do.
  ///
  /// This constructor is declared private for now, because it is an
  /// implementation detail of duplicate().  We may choose to expose
  /// it in the future.
  ///
  /// Preconditions:
  ///   - <tt>rawMpiComm.get() != NULL</tt>
  ///   - <tt>*rawMpiComm != MPI_COMM_NULL</tt>
  ///   - \c defaultTag is the same on all processes in the given
  ///     communicator
  MpiComm (const RCP<const OpaqueWrapper<MPI_Comm> >& rawMpiComm,
           const int defaultTag);

  /**
   * \brief Construct a communicator with a new context with the same
   *   properties as the original.
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
  MpiComm (const MpiComm<Ordinal>& other);

  /** \brief Return the embedded wrapped opaque <tt>MPI_Comm</tt> object. */
  RCP<const OpaqueWrapper<MPI_Comm> > getRawMpiComm () const {
    return rawMpiComm_;
  }

  /// \brief Set the MPI error handler for this communicator.
  ///
  /// \param errHandler [in] The error handler to set.  If null, do
  ///   nothing.
  ///
  /// MPI lets you set an error handler function specific to each
  /// communicator.  (See Section 8.3 of the MPI 3.0 Standard.)
  /// MpiComm wraps this functionality using this method.  You must
  /// first either create an error handler using
  /// MPI_Comm_create_errhandler() (or MPI_Errhandler_create() if you
  /// are stuck with an MPI 1 implementation), or use one of the
  /// default error handlers that the MPI standard or your MPI
  /// implementation provides.  You will need to wrap the MPI error
  /// handler in an OpaqueWrapper.  (See the documentation of
  /// OpaqueWrapper for the rationale behind not using MPI's opaque
  /// objects directly.)
  ///
  /// MpiComm will not attempt to call MPI_Errhandler_free() on the
  /// error handler you provide.  You are responsible for arranging
  /// that this be done.  Note that MPI_Comm_create_errhandler()
  /// (which creates an error handler, given a function pointer) does
  /// not attach the error handler to an MPI_Comm, so the lifetime of
  /// the error handler is not tied to the MPI_Comm to which it is
  /// assigned.  An error handler can be assigned to more than one
  /// MPI_Comm, in fact.  You just need to guarantee that if you
  /// create a custom error handler, then that handler gets freed at
  /// some point.  "The call to <tt>MPI_FINALIZE</tt> does not free
  /// objects created by MPI calls; these objects are freed using
  /// <tt>MPI_xxx_FREE</tt> calls" (Section 8.7, MPI 3.0 Standard).
  /// Note that it is legitimate to call MPI_Errhandler_free() right
  /// after setting the MPI_Comm's error handler; see Section 8.3.4 of
  /// the MPI 3.0 Standard ("The error handler [given to
  /// MPI_Errhandler_free] will be deallocated after all the objects
  /// associated with it (communicator, window, or file) have been
  /// deallocated").  You might instead attach your error handler as a
  /// attribute to <tt>MPI_COMM_SELF</tt>, in such a way that
  /// MPI_Errhandler_free() will be called when <tt>MPI_COMM_SELF</tt>
  /// is freed (which MPI_Finalize() does automatically).  We do not
  /// take responsibility for doing any of these things; you are
  /// responsible for freeing the error handler.
  ///
  /// Here is an example showing how to change an MpiComm's error
  /// handler.  The default error handler for any <tt>MPI_Comm</tt> is
  /// <tt>MPI_ERRORS_ARE_FATAL</tt>.  This handler immediately aborts
  /// if MPI encounters an error, without returning an error code from
  /// the MPI function.  Suppose that instead you would like MPI
  /// functions to return an error code if MPI should encounter an
  /// error.  (In that case, Teuchos' MPI wrappers will detect the
  /// returned error code and throw an exception with an appropriate
  /// error message.  If MPI aborts immediately on error, Teuchos
  /// won't have the chance to detect and report the error.)  If so,
  /// you may set the error handler to MPI_ERRORS_RETURN, one of MPI's
  /// built-in error handlers.  Here is how you may do this for an
  /// MpiComm:
  /// \code
  /// // Suppose that you've already created this MpiComm.
  /// RCP<const MpiComm<int> > comm = ...;
  ///
  /// // Wrap the error handler.
  /// RCP<const OpaqueWrapper<MPI_Errhandler> > errHandler =
  ///   rcp (new OpaqueWrapper<MPI_Errhandler> (MPI_ERRORS_RETURN));
  /// // Set the MpiComm's error handler.
  /// comm->setErrorHandler (errHandler);
  /// \endcode
  void setErrorHandler (const RCP<const OpaqueWrapper<MPI_Errhandler> >& errHandler);

  //@}
  //! @name Implementation of Comm interface
  //@{

  //! The calling process' rank.
  virtual int getRank() const;

  //! The number of processes in the communicator.
  virtual int getSize() const;

  //! Execute a barrier; must be called collectively.
  virtual void barrier() const;

  /** \brief . */
  virtual void broadcast(
    const int rootRank, const Ordinal bytes, char buffer[]
    ) const;

  //! Gather values from all processes to the root process.
  virtual void
  gather (const Ordinal sendBytes, const char sendBuffer[],
          const Ordinal recvBytes, char recvBuffer[],
          const int root) const;
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
  virtual void scan(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char scanReducts[]
    ) const;
  /** \brief . */
  virtual void send(
    const Ordinal bytes, const char sendBuffer[], const int destRank
    ) const;
  /** \brief . */
  virtual void
  send (const Ordinal bytes,
        const char sendBuffer[],
        const int destRank,
        const int tag) const;
  /** \brief . */
  virtual void ssend(
    const Ordinal bytes, const char sendBuffer[], const int destRank
    ) const;
  //! Variant of ssend() that takes a message tag.
  virtual void
  ssend (const Ordinal bytes,
         const char sendBuffer[],
         const int destRank,
         const int tag) const;
  /** \brief . */
  virtual int receive(
    const int sourceRank, const Ordinal bytes, char recvBuffer[]
    ) const;
  /** \brief . */
  virtual void readySend(
    const ArrayView<const char> &sendBuffer,
    const int destRank
    ) const;
  //! Variant of readySend() that accepts a message tag.
  virtual void
  readySend (const Ordinal bytes,
             const char sendBuffer[],
             const int destRank,
             const int tag) const;
  /** \brief . */
  virtual RCP<CommRequest<Ordinal> > isend(
    const ArrayView<const char> &sendBuffer,
    const int destRank
    ) const;
  //! Variant of isend() that takes a tag.
  virtual RCP<CommRequest<Ordinal> >
  isend (const ArrayView<const char> &sendBuffer,
         const int destRank,
         const int tag) const;
  /** \brief . */
  virtual RCP<CommRequest<Ordinal> > ireceive(
    const ArrayView<char> &Buffer,
    const int sourceRank
    ) const;
  /** \brief . */
  virtual RCP<CommRequest<Ordinal> >
  ireceive (const ArrayView<char> &Buffer,
            const int sourceRank,
            const int tag) const;
  /** \brief . */
  virtual void waitAll(
    const ArrayView<RCP<CommRequest<Ordinal> > > &requests
    ) const;
  /** \brief . */
  virtual void
  waitAll (const ArrayView<RCP<CommRequest<Ordinal> > >& requests,
           const ArrayView<RCP<CommStatus<Ordinal> > >& statuses) const;
  /** \brief . */
  virtual RCP<CommStatus<Ordinal> >
  wait (const Ptr<RCP<CommRequest<Ordinal> > >& request) const;
  /** \brief . */
  virtual RCP< Comm<Ordinal> > duplicate() const;
  /** \brief . */
  virtual RCP< Comm<Ordinal> > split(const int color, const int key) const;
  /** \brief . */
  virtual RCP< Comm<Ordinal> > createSubcommunicator(
    const ArrayView<const int>& ranks) const;

  //@}
  //! @name Implementation of Describable interface
  //@{

  /** \brief . */
  std::string description() const;

  //@}

  // These should be private but the PGI compiler requires them be public

  static int const minTag_ = 26000; // These came from Teuchos::MpiComm???
  static int const maxTag_ = 26099; // ""

  /// \brief The current tag.
  ///
  /// \warning This method is ONLY for use by Teuchos developers.
  ///   Users should not depend on the interface of this method.
  ///   It may change or disappear at any time without warning.
  int getTag () const { return tag_; }

private:

  /// \brief Set internal data members once the rawMpiComm_ data member is valid.
  ///
  /// This method should only be called from MpiComm's constructor.
  void setupMembersFromComm();
  static int tagCounter_;

  /// \brief The "raw" MPI_Comm (communicator).
  ///
  /// We wrap the MPI_Comm so that it is freed automatically when its
  /// reference count goes to zero, if it does need to be freed after
  /// use by calling MPI_Comm_free.  (Predefined MPI_Comm, which
  /// include but are not limited to MPI_COMM_WORLD and MPI_COMM_SELF,
  /// need not and must not be freed after use.)
  RCP<const OpaqueWrapper<MPI_Comm> > rawMpiComm_;

  //! The rank of the calling process.
  int rank_;

  //! The number of processes in the communicator.
  int size_;

  /// \brief The current tag, to use for all MPI functions that need it.
  ///
  /// Each MpiComm instance always uses the same tag.  Different
  /// MpiComm instances use different tags.  The tag is set in
  /// MpiComm's constructor.  Please refer to
  /// <a href="https://software.sandia.gov/bugzilla/show_bug.cgi?id=5740">Bug 5740</a>
  /// for further discussion.
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


/** \brief Helper function that extracts a raw <tt>MPI_Comm</tt> object out of
 * a <tt>Teuchos::MpiComm</tt> wrapper object.
 *
 * <b>Preconditions:</b></ul>
 * <li><tt>dynamic_cast<const MpiComm<Ordinal>*>(&comm) != 0</tt>
 * </ul>
 *
 * If the underlying type is not an <tt>MpiComm<Ordinal></tt> object, then the
 * function with throw an exception which contains the type information as for
 * why it failed.
 *
 * <b>WARNING:</b> The lifetime of the returned <tt>MPI_Comm</tt> object is
 * controlled by the owning <tt>RCP<OpaqueWrapper<MPI_Comm> ></tt> object and
 * is not guaranteed to live the entire life of the program.  Therefore, only
 * use this function to grab and use the underlying <tt>MPI_Comm</tt> object
 * in a vary narrow scope and then forget it.  If you need it again, get it
 * off of the <tt>comm</tt> object each time.
 *
 * If you want a memory safe <tt><tt>RCP<OpaqueWrapper<MPI_Comm> ></tt> to the
 * raw <tt>MPI_Comm</tt> object, then call:
 *
 * \verbatim
 * dyn_cast<const MpiComm<Ordinal> >(comm).getRawMpiComm()
 * \endverbatim
 *
 * \relates MpiComm
 */
template<typename Ordinal>
MPI_Comm
getRawMpiComm(const Comm<Ordinal> &comm);


// ////////////////////////
// Implementations


// Static members


template<typename Ordinal>
int MpiComm<Ordinal>::tagCounter_ = MpiComm<Ordinal>::minTag_;


// Constructors


template<typename Ordinal>
MpiComm<Ordinal>::
MpiComm (const RCP<const OpaqueWrapper<MPI_Comm> >& rawMpiComm)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    rawMpiComm.get () == NULL, std::invalid_argument,
    "Teuchos::MpiComm constructor: The input RCP is null.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    *rawMpiComm == MPI_COMM_NULL, std::invalid_argument,
    "Teuchos::MpiComm constructor: The given MPI_Comm is MPI_COMM_NULL.");

  rawMpiComm_ = rawMpiComm;

  // mfh 09 Jul 2013: Please resist the temptation to modify the given
  // MPI communicator's error handler here.  See Bug 5943.  Note that
  // an MPI communicator's default error handler is
  // MPI_ERRORS_ARE_FATAL, which immediately aborts on error (without
  // returning an error code from the MPI function).  Users who want
  // MPI functions instead to return an error code if they encounter
  // an error, should set the error handler to MPI_ERRORS_RETURN.  DO
  // NOT SET THE ERROR HANDLER HERE!!!  Teuchos' MPI wrappers should
  // always check the error code returned by an MPI function,
  // regardless of the error handler.  Users who want to set the error
  // handler on an MpiComm may call its setErrorHandler method.

  setupMembersFromComm ();
}


template<typename Ordinal>
MpiComm<Ordinal>::
MpiComm (const RCP<const OpaqueWrapper<MPI_Comm> >& rawMpiComm,
         const int defaultTag)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    rawMpiComm.get () == NULL, std::invalid_argument,
    "Teuchos::MpiComm constructor: The input RCP is null.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    *rawMpiComm == MPI_COMM_NULL, std::invalid_argument,
    "Teuchos::MpiComm constructor: The given MPI_Comm is MPI_COMM_NULL.");

  rawMpiComm_ = rawMpiComm;
  // Set size_ (the number of processes in the communicator).
  int err = MPI_Comm_size (*rawMpiComm_, &size_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm constructor: MPI_Comm_size failed with "
    "error \"" << mpiErrorCodeToString (err) << "\".");
  // Set rank_ (the calling process' rank).
  err = MPI_Comm_rank (*rawMpiComm_, &rank_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm constructor: MPI_Comm_rank failed with "
    "error \"" << mpiErrorCodeToString (err) << "\".");
  tag_ = defaultTag; // set the default message tag
}


template<typename Ordinal>
MpiComm<Ordinal>::MpiComm (MPI_Comm rawMpiComm)
{
  TEUCHOS_TEST_FOR_EXCEPTION(rawMpiComm == MPI_COMM_NULL,
    std::invalid_argument, "Teuchos::MpiComm constructor: The given MPI_Comm "
    "is MPI_COMM_NULL.");
  // We don't supply a "free" function here, since this version of the
  // constructor makes the caller responsible for freeing rawMpiComm
  // after use if necessary.
  rawMpiComm_ = opaqueWrapper<MPI_Comm> (rawMpiComm);

  // mfh 09 Jul 2013: Please resist the temptation to modify the given
  // MPI communicator's error handler here.  See Bug 5943.  Note that
  // an MPI communicator's default error handler is
  // MPI_ERRORS_ARE_FATAL, which immediately aborts on error (without
  // returning an error code from the MPI function).  Users who want
  // MPI functions instead to return an error code if they encounter
  // an error, should set the error handler to MPI_ERRORS_RETURN.  DO
  // NOT SET THE ERROR HANDLER HERE!!!  Teuchos' MPI wrappers should
  // always check the error code returned by an MPI function,
  // regardless of the error handler.  Users who want to set the error
  // handler on an MpiComm may call its setErrorHandler method.

  setupMembersFromComm ();
}


template<typename Ordinal>
MpiComm<Ordinal>::MpiComm (const MpiComm<Ordinal>& other) :
  rawMpiComm_ (opaqueWrapper<MPI_Comm> (MPI_COMM_NULL)) // <- This will be set below
{
  // These are logic errors, since they violate MpiComm's invariants.
  RCP<const OpaqueWrapper<MPI_Comm> > origCommPtr = other.getRawMpiComm ();
  TEUCHOS_TEST_FOR_EXCEPTION(origCommPtr == null, std::logic_error,
    "Teuchos::MpiComm copy constructor: "
    "The input's getRawMpiComm() method returns null.");
  MPI_Comm origComm = *origCommPtr;
  TEUCHOS_TEST_FOR_EXCEPTION(origComm == MPI_COMM_NULL, std::logic_error,
    "Teuchos::MpiComm copy constructor: "
    "The input's raw MPI_Comm is MPI_COMM_NULL.");

  // mfh 19 Oct 2012: Don't change the behavior of MpiComm's copy
  // constructor for now.  Later, we'll switch to the version that
  // calls MPI_Comm_dup.  For now, we just copy other's handle over.
  // Note that the new MpiComm's tag is still different than the input
  // MpiComm's tag.  See Bug 5740.
  if (true) {
    rawMpiComm_ = origCommPtr;
  }
  else { // false (not run)
    MPI_Comm newComm;
    const int err = MPI_Comm_dup (origComm, &newComm);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
      "Teuchos::MpiComm copy constructor: MPI_Comm_dup failed with "
      "the following error: " << mpiErrorCodeToString (err));
    // No side effects until after everything has succeeded.
    rawMpiComm_ = opaqueWrapper (newComm, details::safeCommFree);
  }

  setupMembersFromComm ();
}


template<typename Ordinal>
void MpiComm<Ordinal>::setupMembersFromComm ()
{
  int err = MPI_Comm_size (*rawMpiComm_, &size_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm constructor: MPI_Comm_size failed with "
    "error \"" << mpiErrorCodeToString (err) << "\".");
  err = MPI_Comm_rank (*rawMpiComm_, &rank_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm constructor: MPI_Comm_rank failed with "
    "error \"" << mpiErrorCodeToString (err) << "\".");

  // Set the default tag to make unique across all communicators
  if (tagCounter_ > maxTag_) {
    tagCounter_ = minTag_;
  }
  tag_ = tagCounter_++;
  // Ensure that the same tag is used on all processes.
  //
  // FIXME (mfh 09 Jul 2013) This would not be necessary if MpiComm
  // were just to call MPI_Comm_dup (as every library should) when
  // given its communicator.  Of course, MPI_Comm_dup may also be
  // implemented as a collective, and may even be more expensive than
  // a broadcast.  If we do decide to use MPI_Comm_dup, we can get rid
  // of the broadcast below, and also get rid of tag_, tagCounter_,
  // minTag_, and maxTag_.
  MPI_Bcast (&tag_, 1, MPI_INT, 0, *rawMpiComm_);
}


template<typename Ordinal>
void
MpiComm<Ordinal>::
setErrorHandler (const RCP<const OpaqueWrapper<MPI_Errhandler> >& errHandler)
{
  if (! is_null (errHandler)) {
    const int err = details::setCommErrhandler (*getRawMpiComm (), *errHandler);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
      "Teuchos::MpiComm: Setting the MPI_Comm's error handler failed with "
      "error \"" << mpiErrorCodeToString (err) << "\".");
  }
  // Wait to set this until the end, in case setting the error handler
  // doesn't succeed.
  customErrorHandler_ = errHandler;
}

//
// Overridden from Comm
//

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
  const int err = MPI_Barrier (*rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::barrier: MPI_Barrier failed with error \""
    << mpiErrorCodeToString (err) << "\".");
}


template<typename Ordinal>
void MpiComm<Ordinal>::broadcast(
  const int rootRank, const Ordinal bytes, char buffer[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::broadcast(...)"
    );
  const int err = MPI_Bcast (buffer, bytes, MPI_CHAR, rootRank, *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::broadcast: MPI_Bcast failed with error \""
    << mpiErrorCodeToString (err) << "\".");
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
  const int err =
    MPI_Allgather (const_cast<char *>(sendBuffer), sendBytes, MPI_CHAR,
                   recvBuffer, sendBytes, MPI_CHAR, *rawMpiComm_);
  // NOTE: 'sendBytes' is being sent above for the MPI arg recvcount (which is
  // very confusing in the MPI documentation) for MPI_Allgether(...).

  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::gatherAll: MPI_Allgather failed with error \""
    << mpiErrorCodeToString (err) << "\".");
}


template<typename Ordinal>
void
MpiComm<Ordinal>::gather (const Ordinal sendBytes,
                          const char sendBuffer[],
                          const Ordinal recvBytes,
                          char recvBuffer[],
                          const int root) const
{
  (void) recvBytes; // silence compile warning for "unused parameter"

  TEUCHOS_COMM_TIME_MONITOR(
    "Teuchos::MpiComm<"<<OrdinalTraits<Ordinal>::name()<<">::gather(...)"
    );
  const int err =
    MPI_Gather (const_cast<char *> (sendBuffer), sendBytes, MPI_CHAR,
                recvBuffer, sendBytes, MPI_CHAR, root, *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::gather: MPI_Gather failed with error \""
    << mpiErrorCodeToString (err) << "\".");
}


template<typename Ordinal>
void
MpiComm<Ordinal>::
reduceAll (const ValueTypeReductionOp<Ordinal,char> &reductOp,
           const Ordinal bytes,
           const char sendBuffer[],
           char globalReducts[]) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::reduceAll(...)" );
  int err = MPI_SUCCESS;

  Details::MpiReductionOp<Ordinal> opWrap (reductOp);
  MPI_Op op = Details::setMpiReductionOp (opWrap);

  // FIXME (mfh 23 Nov 2014) Ross decided to mash every type into
  // char.  This can cause correctness issues if we're actually doing
  // a reduction over, say, double.  Thus, he creates a custom
  // MPI_Datatype here that represents a contiguous block of char, so
  // that MPI doesn't split up the reduction type and thus do the sum
  // wrong.  It's a hack but it works.

  MPI_Datatype char_block;
  err = MPI_Type_contiguous (bytes, MPI_CHAR, &char_block);
  TEUCHOS_TEST_FOR_EXCEPTION(
    err != MPI_SUCCESS, std::runtime_error, "Teuchos::reduceAll: "
    "MPI_Type_contiguous failed with error \"" << mpiErrorCodeToString (err)
    << "\".");
  err = MPI_Type_commit (&char_block);
  TEUCHOS_TEST_FOR_EXCEPTION(
    err != MPI_SUCCESS, std::runtime_error, "Teuchos::reduceAll: "
    "MPI_Type_commit failed with error \"" << mpiErrorCodeToString (err)
    << "\".");

  if (sendBuffer == globalReducts) {
    // NOTE (mfh 31 May 2017) This is only safe if the communicator is
    // NOT an intercomm.  The usual case is that communicators are
    // intracomms.
    err = MPI_Allreduce (MPI_IN_PLACE, globalReducts, 1,
                         char_block, op, *rawMpiComm_);
  }
  else {
    err = MPI_Allreduce (const_cast<char*> (sendBuffer), globalReducts, 1,
                         char_block, op, *rawMpiComm_);
  }
  if (err != MPI_SUCCESS) {
    // Don't throw until we release the type resources we allocated
    // above.  If freeing fails for some reason, let the memory leak
    // go; we already have more serious problems if MPI_Allreduce
    // doesn't work.
    (void) MPI_Type_free (&char_block);
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error, "Teuchos::reduceAll (MPI, custom op): "
      "MPI_Allreduce failed with error \"" << mpiErrorCodeToString (err)
      << "\".");
  }
  err = MPI_Type_free (&char_block);
  TEUCHOS_TEST_FOR_EXCEPTION(
    err != MPI_SUCCESS, std::runtime_error, "Teuchos::reduceAll: "
    "MPI_Type_free failed with error \"" << mpiErrorCodeToString (err)
    << "\".");
}


template<typename Ordinal>
void MpiComm<Ordinal>::scan(
  const ValueTypeReductionOp<Ordinal,char> &reductOp
  ,const Ordinal bytes, const char sendBuffer[], char scanReducts[]
  ) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::scan(...)" );

  Details::MpiReductionOp<Ordinal> opWrap (reductOp);
  MPI_Op op = Details::setMpiReductionOp (opWrap);
  const int err =
    MPI_Scan (const_cast<char*> (sendBuffer), scanReducts, bytes, MPI_CHAR,
              op, *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::scan: MPI_Scan() failed with error \""
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
    "Teuchos::MpiComm::send: MPI_Send() failed with error \""
    << mpiErrorCodeToString (err) << "\".");
}


template<typename Ordinal>
void
MpiComm<Ordinal>::send (const Ordinal bytes,
                        const char sendBuffer[],
                        const int destRank,
                        const int tag) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::send(...)" );
  const int err = MPI_Send (const_cast<char*> (sendBuffer), bytes, MPI_CHAR,
                            destRank, tag, *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::send: MPI_Send() failed with error \""
    << mpiErrorCodeToString (err) << "\".");
}


template<typename Ordinal>
void
MpiComm<Ordinal>::ssend (const Ordinal bytes,
                         const char sendBuffer[],
                         const int destRank) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::ssend(...)" );

#ifdef TEUCHOS_MPI_COMM_DUMP
  if(show_dump) {
    dumpBuffer<Ordinal,char>(
      "Teuchos::MpiComm<Ordinal>::send(...)"
      ,"sendBuffer", bytes, sendBuffer
      );
  }
#endif // TEUCHOS_MPI_COMM_DUMP

  const int err = MPI_Ssend (const_cast<char*>(sendBuffer), bytes, MPI_CHAR,
                             destRank, tag_, *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::send: MPI_Ssend() failed with error \""
    << mpiErrorCodeToString (err) << "\".");
}

template<typename Ordinal>
void
MpiComm<Ordinal>::ssend (const Ordinal bytes,
                         const char sendBuffer[],
                         const int destRank,
                         const int tag) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::ssend(...)" );
  const int err =
    MPI_Ssend (const_cast<char*>(sendBuffer), bytes, MPI_CHAR,
               destRank, tag, *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::send: MPI_Ssend() failed with error \""
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
    MPI_Rsend (const_cast<char*>(sendBuffer.getRawPtr()), static_cast<int>(sendBuffer.size()),
               MPI_CHAR, destRank, tag_, *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::readySend: MPI_Rsend() failed with error \""
    << mpiErrorCodeToString (err) << "\".");
}


template<typename Ordinal>
void MpiComm<Ordinal>::
readySend (const Ordinal bytes,
           const char sendBuffer[],
           const int destRank,
           const int tag) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::readySend" );
  const int err =
    MPI_Rsend (const_cast<char*> (sendBuffer), bytes,
               MPI_CHAR, destRank, tag, *rawMpiComm_);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::readySend: MPI_Rsend() failed with error \""
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
    "Teuchos::MpiComm::receive: MPI_Recv() failed with error \""
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
RCP<CommRequest<Ordinal> >
MpiComm<Ordinal>::isend (const ArrayView<const char> &sendBuffer,
                         const int destRank) const
{
  using Teuchos::as;
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::isend(...)" );

  MPI_Request rawMpiRequest = MPI_REQUEST_NULL;
  const int err =
    MPI_Isend (const_cast<char*> (sendBuffer.getRawPtr ()),
               as<Ordinal> (sendBuffer.size ()), MPI_CHAR,
               destRank, tag_, *rawMpiComm_, &rawMpiRequest);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::isend: MPI_Isend() failed with error \""
    << mpiErrorCodeToString (err) << "\".");

  return mpiCommRequest<Ordinal> (rawMpiRequest, sendBuffer.size ());
}


template<typename Ordinal>
RCP<CommRequest<Ordinal> >
MpiComm<Ordinal>::
isend (const ArrayView<const char> &sendBuffer,
       const int destRank,
       const int tag) const
{
  using Teuchos::as;
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::isend(...)" );

  MPI_Request rawMpiRequest = MPI_REQUEST_NULL;
  const int err =
    MPI_Isend (const_cast<char*> (sendBuffer.getRawPtr ()),
               as<Ordinal> (sendBuffer.size ()), MPI_CHAR,
               destRank, tag, *rawMpiComm_, &rawMpiRequest);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::isend: MPI_Isend() failed with error \""
    << mpiErrorCodeToString (err) << "\".");

  return mpiCommRequest<Ordinal> (rawMpiRequest, sendBuffer.size ());
}


template<typename Ordinal>
RCP<CommRequest<Ordinal> >
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
    "Teuchos::MpiComm::ireceive: MPI_Irecv() failed with error \""
    << mpiErrorCodeToString (err) << "\".");

  return mpiCommRequest<Ordinal> (rawMpiRequest, recvBuffer.size());
}

template<typename Ordinal>
RCP<CommRequest<Ordinal> >
MpiComm<Ordinal>::ireceive (const ArrayView<char> &recvBuffer,
                            const int sourceRank,
                            const int tag) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::ireceive(...)" );

  // A negative source rank indicates MPI_ANY_SOURCE, namely that we
  // will take an incoming message from any process, as long as the
  // tag matches.
  const int theSrcRank = (sourceRank < 0) ? MPI_ANY_SOURCE : sourceRank;

  MPI_Request rawMpiRequest = MPI_REQUEST_NULL;
  const int err =
    MPI_Irecv (const_cast<char*> (recvBuffer.getRawPtr ()), recvBuffer.size (),
               MPI_CHAR, theSrcRank, tag, *rawMpiComm_, &rawMpiRequest);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error,
    "Teuchos::MpiComm::ireceive: MPI_Irecv() failed with error \""
    << mpiErrorCodeToString (err) << "\".");

  return mpiCommRequest<Ordinal> (rawMpiRequest, recvBuffer.size ());
}

namespace {
  // Called by the two-argument MpiComm::waitAll() variant.
  template<typename Ordinal>
  void
  waitAllImpl (const ArrayView<RCP<CommRequest<Ordinal> > >& requests,
               const ArrayView<MPI_Status>& rawMpiStatuses)
  {
    typedef typename ArrayView<RCP<CommRequest<Ordinal> > >::size_type size_type;
    const size_type count = requests.size();
    // waitAllImpl() is not meant to be called by users, so it's a bug
    // for the two views to have different lengths.
    TEUCHOS_TEST_FOR_EXCEPTION(rawMpiStatuses.size() != count,
      std::logic_error, "Teuchos::MpiComm's waitAllImpl: rawMpiStatus.size() = "
      << rawMpiStatuses.size() << " != requests.size() = " << requests.size()
      << ".  Please report this bug to the Tpetra developers.");
    if (count == 0) {
      return; // No requests on which to wait
    }

    // MpiComm wraps MPI and can't expose any MPI structs or opaque
    // objects.  Thus, we have to unpack requests into a separate array.
    // If that's too slow, then your code should just call into MPI
    // directly.
    //
    // Pull out the raw MPI requests from the wrapped requests.
    // MPI_Waitall should not fail if a request is MPI_REQUEST_NULL, but
    // we keep track just to inform the user.
    bool someNullRequests = false;
    Array<MPI_Request> rawMpiRequests (count, MPI_REQUEST_NULL);
    for (int i = 0; i < count; ++i) {
      RCP<CommRequest<Ordinal> > request = requests[i];
      if (! is_null (request)) {
        RCP<MpiCommRequestBase<Ordinal> > mpiRequest =
          rcp_dynamic_cast<MpiCommRequestBase<Ordinal> > (request);
        // releaseRawMpiRequest() sets the MpiCommRequest's raw
        // MPI_Request to MPI_REQUEST_NULL.  This makes waitAll() not
        // satisfy the strong exception guarantee.  That's OK because
        // MPI_Waitall() doesn't promise that it satisfies the strong
        // exception guarantee, and we would rather conservatively
        // invalidate the handles than leave dangling requests around
        // and risk users trying to wait on the same request twice.
        rawMpiRequests[i] = mpiRequest->releaseRawMpiRequest();
      }
      else { // Null requests map to MPI_REQUEST_NULL
        rawMpiRequests[i] = MPI_REQUEST_NULL;
        someNullRequests = true;
      }
    }

    // This is the part where we've finally peeled off the wrapper and
    // we can now interact with MPI directly.
    //
    // One option in the one-argument version of waitAll() is to ignore
    // the statuses completely.  MPI lets you pass in the named constant
    // MPI_STATUSES_IGNORE for the MPI_Status array output argument in
    // MPI_Waitall(), which would tell MPI not to bother with the
    // statuses.  However, we want the statuses because we can use them
    // for detailed error diagnostics in case something goes wrong.
    const int err = MPI_Waitall (count, rawMpiRequests.getRawPtr(),
                                 rawMpiStatuses.getRawPtr());

    // In MPI_Waitall(), an error indicates that one or more requests
    // failed.  In that case, there could be requests that completed
    // (their MPI_Status' error field is MPI_SUCCESS), and other
    // requests that have not completed yet but have not necessarily
    // failed (MPI_PENDING).  We make no attempt here to wait on the
    // pending requests.  It doesn't make sense for us to do so, because
    // in general Teuchos::Comm doesn't attempt to provide robust
    // recovery from failed messages.
    if (err != MPI_SUCCESS) {
      if (err == MPI_ERR_IN_STATUS) {
        //
        // When MPI_Waitall returns MPI_ERR_IN_STATUS (a standard error
        // class), it's telling us to check the error codes in the
        // returned statuses.  In that case, we do so and generate a
        // detailed exception message.
        //
        // Figure out which of the requests failed.
        Array<std::pair<size_type, int> > errorLocationsAndCodes;
        for (size_type k = 0; k < rawMpiStatuses.size(); ++k) {
          const int curErr = rawMpiStatuses[k].MPI_ERROR;
          if (curErr != MPI_SUCCESS) {
            errorLocationsAndCodes.push_back (std::make_pair (k, curErr));
          }
        }
        const size_type numErrs = errorLocationsAndCodes.size();
        if (numErrs > 0) {
          // There was at least one error.  Assemble a detailed
          // exception message reporting which requests failed,
          // their error codes, and their source
          std::ostringstream os;
          os << "Teuchos::MpiComm::waitAll: MPI_Waitall() failed with error \""
             << mpiErrorCodeToString (err) << "\".  Of the " << count
             << " total request" << (count != 1 ? "s" : "") << ", " << numErrs
             << " failed.  Here are the indices of the failed requests, and the "
            "error codes extracted from their returned MPI_Status objects:"
             << std::endl;
          for (size_type k = 0; k < numErrs; ++k) {
            const size_type errInd = errorLocationsAndCodes[k].first;
            os << "Request " << errInd << ": MPI_ERROR = "
               << mpiErrorCodeToString (rawMpiStatuses[errInd].MPI_ERROR)
               << std::endl;
          }
          if (someNullRequests) {
            os << "  On input to MPI_Waitall, there was at least one MPI_"
              "Request that was MPI_REQUEST_NULL.  MPI_Waitall should not "
              "normally fail in that case, but we thought we should let you know "
              "regardless.";
          }
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, os.str());
        }
        // If there were no actual errors in the returned statuses,
        // well, then I guess everything is OK.  Just keep going.
      }
      else {
        std::ostringstream os;
        os << "Teuchos::MpiComm::waitAll: MPI_Waitall() failed with error \""
           << mpiErrorCodeToString (err) << "\".";
        if (someNullRequests) {
          os << "  On input to MPI_Waitall, there was at least one MPI_Request "
            "that was MPI_REQUEST_NULL.  MPI_Waitall should not normally fail in "
            "that case, but we thought we should let you know regardless.";
        }
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, os.str());
      }
    }

    // Invalidate the input array of requests by setting all entries
    // to null.
    std::fill (requests.begin(), requests.end(), null);
  }



  // Called by the one-argument MpiComm::waitAll() variant.
  template<typename Ordinal>
  void
  waitAllImpl (const ArrayView<RCP<CommRequest<Ordinal> > >& requests)
  {
    typedef typename ArrayView<RCP<CommRequest<Ordinal> > >::size_type size_type;
    const size_type count = requests.size ();
    if (count == 0) {
      return; // No requests on which to wait
    }

    // MpiComm wraps MPI and can't expose any MPI structs or opaque
    // objects.  Thus, we have to unpack requests into a separate
    // array.  If that's too slow, then your code should just call
    // into MPI directly.
    //
    // Pull out the raw MPI requests from the wrapped requests.
    // MPI_Waitall should not fail if a request is MPI_REQUEST_NULL,
    // but we keep track just to inform the user.
    bool someNullRequests = false;
    Array<MPI_Request> rawMpiRequests (count, MPI_REQUEST_NULL);
    for (int i = 0; i < count; ++i) {
      RCP<CommRequest<Ordinal> > request = requests[i];
      if (! request.is_null ()) {
        RCP<MpiCommRequestBase<Ordinal> > mpiRequest =
          rcp_dynamic_cast<MpiCommRequestBase<Ordinal> > (request);
        // releaseRawMpiRequest() sets the MpiCommRequest's raw
        // MPI_Request to MPI_REQUEST_NULL.  This makes waitAll() not
        // satisfy the strong exception guarantee.  That's OK because
        // MPI_Waitall() doesn't promise that it satisfies the strong
        // exception guarantee, and we would rather conservatively
        // invalidate the handles than leave dangling requests around
        // and risk users trying to wait on the same request twice.
        rawMpiRequests[i] = mpiRequest->releaseRawMpiRequest ();
      }
      else { // Null requests map to MPI_REQUEST_NULL
        rawMpiRequests[i] = MPI_REQUEST_NULL;
        someNullRequests = true;
      }
    }

    // This is the part where we've finally peeled off the wrapper and
    // we can now interact with MPI directly.
    //
    // MPI lets us pass in the named constant MPI_STATUSES_IGNORE for
    // the MPI_Status array output argument in MPI_Waitall(), which
    // tells MPI not to bother writing out the statuses.
    const int err = MPI_Waitall (count, rawMpiRequests.getRawPtr(),
                                 MPI_STATUSES_IGNORE);

    // In MPI_Waitall(), an error indicates that one or more requests
    // failed.  In that case, there could be requests that completed
    // (their MPI_Status' error field is MPI_SUCCESS), and other
    // requests that have not completed yet but have not necessarily
    // failed (MPI_PENDING).  We make no attempt here to wait on the
    // pending requests.  It doesn't make sense for us to do so,
    // because in general Teuchos::Comm doesn't attempt to provide
    // robust recovery from failed messages.
    if (err != MPI_SUCCESS) {
      std::ostringstream os;
      os << "Teuchos::MpiComm::waitAll: MPI_Waitall() failed with error \""
         << mpiErrorCodeToString (err) << "\".";
      if (someNullRequests) {
        os << std::endl << "On input to MPI_Waitall, there was at least one "
          "MPI_Request that was MPI_REQUEST_NULL.  MPI_Waitall should not "
          "normally fail in that case, but we thought we should let you know "
          "regardless.";
      }
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, os.str());
    }

    // Invalidate the input array of requests by setting all entries
    // to null.  We delay this until the end, since some
    // implementations of CommRequest might hold the only reference to
    // the communication buffer, and we don't want that to go away
    // until we've waited on the communication operation.
    std::fill (requests.begin(), requests.end(), null);
  }

} // namespace (anonymous)



template<typename Ordinal>
void
MpiComm<Ordinal>::
waitAll (const ArrayView<RCP<CommRequest<Ordinal> > >& requests) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::waitAll(requests)" );
  // Call the one-argument version of waitAllImpl, to avoid overhead
  // of handling statuses (which the user didn't want anyway).
  waitAllImpl<Ordinal> (requests);
}


template<typename Ordinal>
void
MpiComm<Ordinal>::
waitAll (const ArrayView<RCP<CommRequest<Ordinal> > >& requests,
         const ArrayView<RCP<CommStatus<Ordinal> > >& statuses) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::waitAll(requests, statuses)" );

  typedef typename ArrayView<RCP<CommRequest<Ordinal> > >::size_type size_type;
  const size_type count = requests.size();

  TEUCHOS_TEST_FOR_EXCEPTION(count != statuses.size(),
    std::invalid_argument, "Teuchos::MpiComm::waitAll: requests.size() = "
    << count << " != statuses.size() = " << statuses.size() << ".");

  Array<MPI_Status> rawMpiStatuses (count);
  waitAllImpl<Ordinal> (requests, rawMpiStatuses());

  // Repackage the raw MPI_Status structs into the wrappers.
  for (size_type i = 0; i < count; ++i) {
    statuses[i] = mpiCommStatus<Ordinal> (rawMpiStatuses[i]);
  }
}


template<typename Ordinal>
RCP<CommStatus<Ordinal> >
MpiComm<Ordinal>::wait (const Ptr<RCP<CommRequest<Ordinal> > >& request) const
{
  TEUCHOS_COMM_TIME_MONITOR( "Teuchos::MpiComm::wait(...)" );

  if (is_null (*request)) {
    return null; // Nothing to wait on ...
  }
  else {
    RCP<CommStatus<Ordinal> > status = (*request)->wait ();
    // mfh 22 Oct 2012: The unit tests expect waiting on the
    // CommRequest to invalidate it by setting it to null.
    *request = null;
    return status;
  }
}

template<typename Ordinal>
RCP< Comm<Ordinal> >
MpiComm<Ordinal>::duplicate() const
{
  MPI_Comm origRawComm = *rawMpiComm_;
  MPI_Comm newRawComm = MPI_COMM_NULL;
  const int err = MPI_Comm_dup (origRawComm, &newRawComm);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, "Teuchos"
    "::MpiComm::duplicate: MPI_Comm_dup failed with the following error: "
    << mpiErrorCodeToString (err));

  // Wrap the raw communicator, and pass the (const) wrapped
  // communicator to MpiComm's constructor.  We created the raw comm,
  // so we have to supply a function that frees it after use.
  RCP<OpaqueWrapper<MPI_Comm> > wrapped =
    opaqueWrapper<MPI_Comm> (newRawComm, details::safeCommFree);
  // Since newComm's raw MPI_Comm is the result of an MPI_Comm_dup,
  // its messages cannot collide with those of any other MpiComm.
  // This means we can assign its tag without an MPI_Bcast.
  RCP<MpiComm<Ordinal> > newComm =
    rcp (new MpiComm<Ordinal> (wrapped.getConst (), minTag_));
  return rcp_implicit_cast<Comm<Ordinal> > (newComm);
}


template<typename Ordinal>
RCP< Comm<Ordinal> >
MpiComm<Ordinal>::split(const int color, const int key) const
{
  MPI_Comm newComm;
  const int splitReturn =
    MPI_Comm_split (*rawMpiComm_,
                    color < 0 ? MPI_UNDEFINED : color,
                    key,
                    &newComm);
  TEUCHOS_TEST_FOR_EXCEPTION(
    splitReturn != MPI_SUCCESS,
    std::logic_error,
    "Teuchos::MpiComm::split: Failed to create communicator with color "
    << color << "and key " << key << ".  MPI_Comm_split failed with error \""
    << mpiErrorCodeToString (splitReturn) << "\".");
  if (newComm == MPI_COMM_NULL) {
    return RCP< Comm<Ordinal> >();
  } else {
    RCP<const OpaqueWrapper<MPI_Comm> > wrapped =
      opaqueWrapper<MPI_Comm> (newComm, details::safeCommFree);
    // Since newComm's raw MPI_Comm is the result of an
    // MPI_Comm_split, its messages cannot collide with those of any
    // other MpiComm.  This means we can assign its tag without an
    // MPI_Bcast.
    return rcp (new MpiComm<Ordinal> (wrapped, minTag_));
  }
}


template<typename Ordinal>
RCP< Comm<Ordinal> >
MpiComm<Ordinal>::createSubcommunicator(const ArrayView<const int> &ranks) const
{
  int err = MPI_SUCCESS; // For error codes returned by MPI functions

  // Get the group that this communicator is in.
  MPI_Group thisGroup;
  err = MPI_Comm_group (*rawMpiComm_, &thisGroup);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::logic_error,
    "Failed to obtain the current communicator's group.  "
    "MPI_Comm_group failed with error \""
    << mpiErrorCodeToString (err) << "\".");

  // Create a new group with the specified members.
  MPI_Group newGroup;
  // It's rude to cast away const, but MPI functions demand it.
  //
  // NOTE (mfh 14 Aug 2012) Please don't ask for &ranks[0] unless you
  // know that ranks.size() > 0.  That's why I'm using getRawPtr().
  err = MPI_Group_incl (thisGroup, ranks.size(),
                        const_cast<int*> (ranks.getRawPtr ()), &newGroup);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::logic_error,
    "Failed to create subgroup.  MPI_Group_incl failed with error \""
    << mpiErrorCodeToString (err) << "\".");

  // Create a new communicator from the new group.
  MPI_Comm newComm;
  try {
    err = MPI_Comm_create (*rawMpiComm_, newGroup, &newComm);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::logic_error,
      "Failed to create subcommunicator.  MPI_Comm_create failed with error \""
      << mpiErrorCodeToString (err) << "\".");
  } catch (...) {
    // Attempt to free the new group before rethrowing.  If
    // successful, this will prevent a memory leak due to the "lost"
    // group that was allocated successfully above.  Since we're
    // throwing std::logic_error anyway, we can only promise
    // best-effort recovery; thus, we don't check the error code.
    (void) MPI_Group_free (&newGroup);
    (void) MPI_Group_free (&thisGroup);
    throw;
  }

  // We don't need the group any more, so free it.
  err = MPI_Group_free (&newGroup);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::logic_error,
    "Failed to free subgroup.  MPI_Group_free failed with error \""
    << mpiErrorCodeToString (err) << "\".");
  err = MPI_Group_free (&thisGroup);
  TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::logic_error,
    "Failed to free subgroup.  MPI_Group_free failed with error \""
    << mpiErrorCodeToString (err) << "\".");

  if (newComm == MPI_COMM_NULL) {
    return RCP<Comm<Ordinal> > ();
  } else {
    using Teuchos::details::safeCommFree;
    typedef OpaqueWrapper<MPI_Comm> ow_type;
    RCP<const ow_type> wrapper =
      rcp_implicit_cast<const ow_type> (opaqueWrapper (newComm, safeCommFree));
    // Since newComm's raw MPI_Comm is the result of an
    // MPI_Comm_create, its messages cannot collide with those of any
    // other MpiComm.  This means we can assign its tag without an
    // MPI_Bcast.
    return rcp (new MpiComm<Ordinal> (wrapper, minTag_));
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


template<typename Ordinal>
MPI_Comm
Teuchos::getRawMpiComm(const Comm<Ordinal> &comm)
{
  return *(
    dyn_cast<const MpiComm<Ordinal> >(comm).getRawMpiComm()
    );
}


#endif // HAVE_TEUCHOS_MPI
#endif // TEUCHOS_MPI_COMM_HPP

