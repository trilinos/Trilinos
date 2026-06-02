// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_DEFAULTMPICOMM_DECL_HPP
#define TEUCHOS_DEFAULTMPICOMM_DECL_HPP

#include "Teuchos_Comm.hpp"

namespace Teuchos {

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
  MpiCommStatus (MPI_Status status);

  //! Destructor (declared virtual for memory safety)
  virtual ~MpiCommStatus();

  //! The source rank that sent the message.
  OrdinalType getSourceRank ();

  //! The tag of the received message.
  OrdinalType getTag ();

  //! The error code of the received message.
  OrdinalType getError ();

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
RCP<MpiCommStatus<OrdinalType> >
mpiCommStatus (MPI_Status rawMpiStatus);

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
  MpiCommRequestBase ();

  //! Constructor (from a raw MPI_Request).
  MpiCommRequestBase (MPI_Request rawMpiRequest);

  /// \brief Return and relinquish ownership of the raw MPI_Request.
  ///
  /// "Relinquish ownership" means that this object sets its raw
  /// MPI_Request to <tt>MPI_REQUEST_NULL</tt>, but returns the
  /// original MPI_Request.  This effectively gives the caller
  /// ownership of the raw MPI_Request.  This prevents hanging
  /// requests.
  MPI_Request releaseRawMpiRequest();

  //! Whether the raw MPI_Request is <tt>MPI_REQUEST_NULL</tt>.
  bool isNull() const;

  bool isReady();

  /// \brief Wait on this communication request to complete.
  ///
  /// This is a blocking operation.  The user is responsible for
  /// avoiding deadlock.  (For example, a receive must have a matching
  /// send, otherwise a wait on the receive will wait forever.)
  RCP<CommStatus<OrdinalType> > wait ();

  /// \brief Cancel the communication request, and return its status.
  ///
  /// If this request is invalid or has already been invalidated, this
  /// method returns null.
  RCP<CommStatus<OrdinalType> > cancel ();

  //! Destructor; cancels the request if it is still pending.
  virtual ~MpiCommRequestBase ();

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
  MpiCommRequest ();

  //! Constructor (from a raw MPI_Request).
  MpiCommRequest (MPI_Request rawMpiRequest,
                  const ArrayView<char>::size_type numBytesInMessage);

  /// \brief Number of bytes in the nonblocking send or receive request.
  ///
  /// Remembering this is inexpensive, and is also useful for
  /// debugging (e.g., for detecting whether the send and receive have
  /// matching message lengths).
  ArrayView<char>::size_type numBytes () const;

  //! Destructor; cancels the request if it is still pending.
  virtual ~MpiCommRequest ();

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
RCP<MpiCommRequest<OrdinalType> >
mpiCommRequest (MPI_Request rawMpiRequest,
                const ArrayView<char>::size_type numBytes);

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

  int incrementTag() {
    ++tag_;
    if (tag_ == std::numeric_limits<int>::max())
      tag_ = 0;
    return tag_;
  }

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
              const RCP<const OpaqueWrapper<MPI_Comm> > &rawMpiComm,
              const int defaultTag
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

}

#endif
