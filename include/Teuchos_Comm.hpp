// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_COMM_HPP
#define TEUCHOS_COMM_HPP

#include "Teuchos_ReductionOp.hpp"
#include "Teuchos_ArrayRCP.hpp"


namespace Teuchos {

/// \class CommStatus
/// \brief Encapsulation of the result of a receive (blocking or nonblocking).
///
/// An instance of this class encapsulates the result of a receive.
/// (An MPI implementation would wrap MPI_Status.)  You can query it
/// for information like the rank of the process that sent you the
/// message.  (This is useful if your receive specified a negative
/// source rank, indicating that you would accept a message from any
/// process in the communicator.)
///
/// \tparam OrdinalType The same template parameter as Comm.  Only use
///   \c int here.  We only make this a template class for
///   compatibility with Comm.
///
/// \note For now, this class only exposes the rank of the process
///   that sent the message (the "source rank") and its tag.  Later,
///   we might expose other fields of MPI_Status in this interface.
///   For now, you can attempt a dynamic cast to MpiCommStatus to
///   access all three fields (MPI_SOURCE, MPI_TAG, and MPI_ERROR).
template<class OrdinalType>
class CommStatus {
public:
  //! Destructor (declared virtual for memory safety)
  virtual ~CommStatus() {}

  //! The source rank that sent the message.
  virtual OrdinalType getSourceRank () = 0;

  //! The tag of the received message.
  virtual OrdinalType getTag () = 0;
};

// Forward declaration for CommRequest::wait.
template<class OrdinalType>
class Comm;

/// \class CommRequest
/// \brief Encapsulation of a pending nonblocking communication operation.
/// \tparam OrdinalType Same as the template parameter of Comm.
///
/// An instance of (a subclass of) this class represents a nonblocking
/// communication operation, such as a nonblocking send, receive, or
/// collective.  To wait on the communication operation, you may give
/// the CommRequest to functions like wait() or waitAll() (which may
/// be found in Teuchos_CommHelpers.hpp).  Here is an example of how
/// to use wait().
/// \code
/// const int sourceRank = ...; // Rank of the sending process.
/// RCP<const Comm<int> > comm = ...; // The communicator.
/// ArrayRCP<double> buf (...); // Buffer for incoming data.
/// RCP<CommRequest<int> > req = ireceive (comm, buf, sourceRank);
///
/// // ... Do some other things ...
///
/// // Wait on the request.  This blocks on the sending process.
/// // When it finishes, it invalidates the req reference, and
/// // returns a status (which wraps MPI_Status in an MPI
/// // implementation).
/// RCP<CommStatus<int> > status = wait (comm, ptr (&req));
/// \endcode
///
/// This object's destructor cancels the request without
/// communication.  If you wish, you may rely on this behavior for
/// speculative communication.  For example:
/// \code
/// const int sourceRank = ...; // Rank of the sending process.
/// RCP<const Comm<int> > comm = ...; // The communicator.
/// ArrayRCP<double> buf (...); // Buffer for incoming data.
/// RCP<CommRequest<int> > req = ireceive (comm, buf, sourceRank);
///
/// // ... Do some other things ...
/// // ... Find out we didn't need to receive data ...
///
/// // This cancels the request.  We could also just let
/// // the one reference to the request fall out of scope.
/// req = null;
/// \endcode
///
/// \note To implementers: The MPI implementation of this class
///   (MpiCommRequest) wraps MPI_Request.  The MPI version of
///   waitAll() will need to unpack the array of wrapped requests, and
///   then pack up the resulting MPI_Request after waiting on them.
///   It would be preferable to have a class \c CommRequests that
///   encapsulates a set of requests, so that you can avoid this
///   unpacking and packing.
template<class OrdinalType>
class CommRequest : public Teuchos::Describable {
public:
  /// \brief Destructor; cancels the request if it is still pending.
  ///
  /// Canceling a communication request must always be a local
  /// operation.  An MPI implementation may achieve this by first
  /// calling MPI_Cancel to cancel the request, then calling MPI_Wait
  /// (which behaves as a local operation for a canceled request) to
  /// complete the canceled request (as required by the MPI standard).
  virtual ~CommRequest() {}

  virtual bool isReady() = 0;

  /// Wait on this request (a blocking operation).
  virtual RCP<CommStatus<OrdinalType> > wait () = 0;
};

/// \class Comm
/// \brief Abstract interface for distributed-memory communication.
/// \tparam Ordinal Type of indices used for communication.
///
/// \section Teuchos_Comm_What What is Comm?
///
/// This class is Teuchos' interface to distributed-memory
/// communication between one or more parallel processes.  It presents
/// an interface very much like that of MPI (the Message Passing
/// Interface).  Teuchos provides two implementations of Comm:
/// - An MPI (Message Passing Interface) implementation, MpiComm
/// - A "serial" implementation, SerialComm, that only has one process
///
/// Comm is an abstract interface.  You cannot create a Comm directly.
/// You have to create one of the subclasses.  The normal way to
/// handle a Comm is to pass it around using RCP (a reference-counted
/// "smart" pointer).  For example:
///
/// \code
/// // Make a Comm.  This one happens to wrap MPI_COMM_WORLD.
/// RCP<const Comm<int> > comm = rcp (new MpiComm (MPI_COMM_WORLD));
/// // Equivalent of MPI_Comm_rank
/// const int myRank = comm->getRank ();
/// // Equivalent of MPI_Comm_size
/// const int numProcs = comm->getSize ();
/// // Equivalent of MPI_Comm_barrier
/// comm->barrier ();
/// \endcode
///
/// Comm's communication methods that actually send or receive data
/// accept that data as an array of \c char.  You should never call
/// these methods directly.  Instead, you should use the nonmember
/// "helper" functions in Teuchos_CommHelpers.hpp.  These methods are
/// templated on the \c Packet type, that is, the type of data you
/// want to send or receive.  See the example below.
///
/// \section Teuchos_Comm_Handle Treat <tt>RCP<const Comm<int> ></tt> like an opaque handle
///
/// You should consider an <tt>RCP<const Comm<int> ></tt> as
/// equivalent to the MPI_Comm opaque handle, except that the RCP also
/// does reference counting to ensure memory safety when using the
/// same communicator in different parts of the code.  That is,
/// copying the RCP does not create a new communicator; the following
/// two codes do about the same thing, except with a different syntax
/// (and reference counting in the second code).
///
/// Raw MPI_Comm handles:
/// \code
/// MPI_Comm comm = ...;
/// // sameComm is THE SAME COMMUNICATOR as comm.
/// MPI_Comm sameComm = comm;
/// \endcode
///
/// Reference-counted pointers to Comm:
/// \code
/// RCP<const Comm<int> > comm = ...;
/// // *sameComm is THE SAME COMMUNICATOR as *comm.
/// RCP<const Comm<int> > sameComm = comm;
/// \endcode
///
/// If you want to make a "new communicator" rather than just "copy
/// the handle," you should call the duplicate() method.  This has the
/// same behavior as MPI_Comm_dup (which see).
///
/// The "reference counting" feature means that the subclass of Comm
/// will take care of freeing the underlying MPI_Comm (and any other
/// data structures it may use) by calling MPI_Comm_free if necessary,
/// once the reference count of the RCP goes to zero.
///
/// \warning Do <i>not</i> pass around subclasses of Comm by value!
///   Comm or its subclasses by themselves do not have handle
///   semantics.  Their copy constructors likely do not behave as you
///   would expect if the classes had handle semantics.
///
/// \section Teuchos_Comm_How How do I make a Comm?
///
/// Comm works whether or not you have build Trilinos with MPI
/// support.  If you want to make a "default" Comm that is the
/// equivalent of MPI_COMM_WORLD, but you don't know if your Trilinos
/// with MPI enabled, you may use GlobalMPISession to call MPI_Init if
/// necessary, and DefaultComm to "get a default communicator."  For
/// example:
/// \code
/// int main (int argc, char* argv[]) {
///   using Teuchos::Comm;
///   using Teuchos::DefaultComm;
///   using Teuchos::RCP;
///
///   // This replaces the call to MPI_Init.  If you didn't
///   // build with MPI, this doesn't call MPI functions.
///   Teuchos::GlobalMPISesssion session (&argc, &argv, NULL);
///   // comm is the equivalent of MPI_COMM_WORLD.
///   RCP<const Comm<int> > comm = DefaultComm<int>::getComm ();
///
///   // ... use comm in your code as you would use MPI_COMM_WORLD ...
///
///   // We don't need to call MPI_Finalize, since the
///   // destructor of GlobalMPISession does that for us.
///   return EXIT_SUCCESS;
/// }
/// \endcode
/// This code works whether or not you built Trilinos with MPI
/// support.  It is not necessary to use GlobalMPISession, but it's
/// useful so you don't have to remember to call MPI_Finalize.  If you
/// don't want to use GlobalMPISession, you can still call
/// <tt>DefaultComm<int>::getComm()</tt>, though you must have called
/// MPI_Init first if you build Trilinos with MPI support.
/// Furthermore, if you know MPI is present, you don't need to use
/// DefaultComm.  You may simply pass MPI_COMM_WORLD directly to
/// MpiComm, like this:
/// \code
/// RCP<const Comm<int> > comm = rcp (new MpiComm (MPI_COMM_WORLD));
/// \endcode
/// You may also pass an arbitrary MPI_Comm directly into MpiComm's
/// constructor, though you are responsible for freeing it after use
/// (via MPI_Comm_free) if necessary.  You may automate the freeing
/// of your MPI_Comm by using OpaqueWrapper (which see).
///
/// \section Teuchos_Comm_Use How do I use Comm?
///
/// As we mentioned above, for communication of data with Comm, you
/// you should use the nonmember "helper" functions in
/// Teuchos_CommHelpers.hpp.  These methods are templated on the
/// <tt>Packet</tt> type, that is, the type of data you want to send
/// or receive.  For example, suppose you have two processes (with
/// ranks 0 and 1, respectively), and you want to send an array of
/// 10 <tt>double</tt> from Process 0 to Process 1.  Both processes have
/// defined <tt>RCP<const Comm<int> > comm</tt> as above.  Here is the
/// code on Process 0:
/// \code
/// const int count = 10; // Send 10 doubles
/// double values[10] = ...;
/// const int destinationRank = 1; // Send to Process 1
/// // You may be able to omit the template arguments of 'send' here.
/// Teuchos::send<int, double> (*comm, 10, values, destinationRank);
/// \endcode
/// Here is the code on Process 1:
/// \code
/// const int count = 10; // Receive 10 doubles
/// double values[10]; // Will be overwritten by receive
/// const int sourceRank = 0; // Receive from Process 0
/// // You may be able to omit the template arguments of 'receive' here.
/// Teuchos::receive<int, double> (*comm, sourceRank, 10, values);
/// \endcode
/// Please refer to the documentation in Teuchos_CommHelpers.hpp for
/// more details.
///
/// \section Teuchos_Comm_Former Former documentation
///
/// This interface is templated on the ordinal type but only deals with buffers
/// of untyped data represented as arrays <tt>char</tt> type. All reduction
/// operations that are initiated by the concreate communicator object are
/// performed by user-defined <tt>ReductOpBase</tt> objects.  It is the
/// responsibility of the <tt>ReductOpBase</tt> object to know what the currect
/// data type is, to perform casts or serializations/unserializations to and
/// from <tt>char[]</tt> buffers, and to know how to reduce the objects
/// correctly.  It is strictly up to the client to correctly convert data types
/// to <tt>char[]</tt> arrays but there is a great deal of helper code to make
/// this easy and safe.
template<typename Ordinal>
class Comm : virtual public Describable {
public:
  /// \brief The current tag.
  ///
  /// \warning This method is ONLY for use by Teuchos developers.
  ///   Users should not depend on the interface of this method.
  ///   It may change or disappear at any time without warning.
  virtual int getTag () const = 0;

  //! @name Destructor
  //@{

  //! Destructor, declared virtual for safety of derived classes.
  virtual ~Comm() {}
  //@}

  //! @name Query functions
  //@{

  /** \brief Returns the rank of this process.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>0 <= return && return < this->getSize()</tt>
   * </ul>
   */
  virtual int getRank() const = 0;

  /** \brief Returns the number of processes that make up this communicator.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>return > 0</tt>
   * </ul>
   */
  virtual int getSize() const = 0;

  //@}

  //! @name Collective Operations
  //@{

  /** \brief Pause every process in <tt>*this</tt> communicator until all the
   * processes reach this point.
   */
  virtual void barrier() const = 0;

  /** \brief Broadcast values from the root process to the slave processes.
   *
   * \param rootRank [in] The rank of the root process.
   *
   * \param count [in] The number of bytes in <tt>buffer[]</tt>.
   *
   * \param buffer [in/out] Array (length <tt>bytes</tt>) of packed data.
   * Must be set on input on the root processes with rank <tt>root</tt>.  On
   * output, each processs, including the root process contains the data.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= rootRank && rootRank < this->getSize()</tt>
   * </ul>
   */
  virtual void broadcast(
    const int rootRank, const Ordinal bytes, char buffer[]
    ) const = 0;

  //! Gather values from all processes to the root process.
  virtual void
  gather (const Ordinal sendBytes, const char sendBuffer[],
          const Ordinal recvBytes, char recvBuffer[],
          const int root) const = 0;

  /** \brief Gather values from each process to collect on all processes.
   *
   * \param sendBytes [in] Number of entires in <tt>sendBuffer[]</tt> on
   * input.
   *
   * \param sendBuffer [in] Array (length <tt>sendBytes</tt>) of data being
   * sent from each process.
   *
   * \param recvBytes [in] Number of entires in <tt>recvBuffer[]</tt> which
   * must be equal to <tt>sendBytes*this->getSize()</tt>.  This field is just
   * here for debug checking.
   *
   * \param recvBuffer [out] Array (length <tt>recvBytes</tt>) of all of the
   * entires sent from each processes.  Specifically,
   * <tt>recvBuffer[sendBytes*j+i]</tt>, for <tt>j=0...this->getSize()-1</tt>
   * and <tt>i=0...sendBytes-1</tt>, is the entry <tt>sendBuffer[i]</tt> from
   * process with rank <tt>j</tt>.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>recvBytes==sendBytes*this->getSize()</tt>
   * </ul>
   */
  virtual void gatherAll(
    const Ordinal sendBytes, const char sendBuffer[]
    ,const Ordinal recvBytes, char recvBuffer[]
    ) const = 0;

  /** \brief Global reduction.
   *
   * \param reductOp [in] The user-defined reduction operation
   *
   * \param bytes [in] The length of the buffers <tt>sendBuffer[]</tt> and
   * <tt>globalReducts[]</tt>.
   *
   * \param sendBuffer [in] Array (length <tt>bytes</tt>) of the data
   * contributed from each process.
   *
   * \param globalReducts [out] Array (length <tt>bytes</tt>) of the global
   * reduction from each process.
   */
  virtual void reduceAll(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char globalReducts[]
    ) const = 0;

  /** \brief Scan reduction.
   *
   * \param reductOp [in] The user-defined reduction operation
   *
   * \param bytes [in] The length of the buffers <tt>sendBuffer[]</tt> and
   * <tt>scanReducts[]</tt>.
   *
   * \param sendBuffer [in] Array (length <tt>bytes</tt>) of the data
   * contributed from each process.
   *
   * \param scanReducts [out] Array (length <tt>bytes</tt>) of the reduction
   * up to and including this process.
   */
        virtual void scan(
    const ValueTypeReductionOp<Ordinal,char> &reductOp
    ,const Ordinal bytes, const char sendBuffer[], char scanReducts[]
    ) const = 0;

  //! @name Blocking Point-to-Point Operations
  //@{

  /** \brief Possibly blocking send of data from this process to another process.
   *
   * This routine does not return until you can reuse the send buffer.
   * Whether this routine blocks depends on whether the MPI
   * implementation buffers.
   *
   * \param bytes [in] The number of bytes of data being passed between
   * processes.
   *
   * \param sendBuffer [in] Array (length <tt>bytes</tt>) of data being sent
   * from this process.  This buffer can be immediately destroyed or reused as
   * soon as the function exits (that is why this function is "blocking").
   *
   * \param destRank [in] The rank of the process to receive the data.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= destRank && destRank < this->getSize()</tt>
   * <li><tt>destRank != this->getRank()</tt>
   * </ul>
   */
  virtual void send(
    const Ordinal bytes, const char sendBuffer[], const int destRank
    ) const = 0;

  //! Variant of send() that takes a tag.
  virtual void
  send (const Ordinal bytes,
        const char sendBuffer[],
        const int destRank,
        const int tag) const = 0;

  /** \brief Always blocking send of data from this process to another process.
   *
   * This routine blocks until the matching receive posts.  After it
   * returns, you are allowed to reuse the send buffer.
   *
   * \param bytes [in] The number of bytes of data being passed between
   * processes.
   *
   * \param sendBuffer [in] Array (length <tt>bytes</tt>) of data being sent
   * from this process.  This buffer can be immediately destroyed or reused as
   * soon as the function exits (that is why this function is "blocking").
   *
   * \param destRank [in] The rank of the process to receive the data.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= destRank && destRank < this->getSize()</tt>
   * <li><tt>destRank != this->getRank()</tt>
   * </ul>
   */
  virtual void ssend(
    const Ordinal bytes, const char sendBuffer[], const int destRank
    ) const = 0;

  //! Variant of ssend() that takes a message tag.
  virtual void
  ssend (const Ordinal bytes,
         const char sendBuffer[],
         const int destRank,
         const int tag) const = 0;

  /** \brief Blocking receive of data from this process to another process.
   *
   * \param sourceRank [in] The rank of the process to receive the data from.
   * If <tt>sourceRank < 0</tt> then data will be received from any process.
   *
   * \param bytes [in] The number of bytes of data being passed between
   * processes.
   *
   * \param recvBuffer [out] Array (length <tt>bytes</tt>) of data being
   * received from this process.  This buffer can be immediately used to
   * access the data as soon as the function exits (that is why this function
   * is "blocking").
   *
   * <b>Preconditions:</b><ul>
   * <li>[<tt>sourceRank >= 0] <tt>sourceRank < this->getSize()</tt>
   * <li><tt>sourceRank != this->getRank()</tt>
   * </ul>
   *
   * \return Returns the senders rank.
   */
  virtual int receive(
    const int sourceRank, const Ordinal bytes, char recvBuffer[]
    ) const = 0;


  /** \brief Ready send of data from this process to another process.
   *
   * \param sendBuffer [in] The data to be sent.
   *
   * \param destRank [in] The rank of the process to receive the data.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= destRank && destRank < this->getSize()</tt>
   * <li><tt>destRank != this->getRank()</tt>
   * </ul>
   */
  virtual void readySend(
    const ArrayView<const char> &sendBuffer,
    const int destRank
    ) const = 0;

  //! Variant of readySend() that accepts a message tag.
  virtual void
  readySend (const Ordinal bytes,
             const char sendBuffer[],
             const int destRank,
             const int tag) const = 0;

  //@}
  //! @name Non-blocking Point-to-Point Operations
  //@{

  /** \brief Non-blocking send.
   *
   * \param sendBuffer [in] The data buffer to be sent.
   *
   * \param destRank [in] The rank of the process to receive the data.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>0 <= destRank && destRank < this->getSize()</tt>
   * <li><tt>destRank != this->getRank()</tt>
   * </ul>
   */
  virtual RCP<CommRequest<Ordinal> > isend(
    const ArrayView<const char> &sendBuffer,
    const int destRank
    ) const = 0;

  //! Variant of isend() that takes a tag.
  virtual RCP<CommRequest<Ordinal> >
  isend (const ArrayView<const char> &sendBuffer,
         const int destRank,
         const int tag) const = 0;

  /** \brief Non-blocking receive.
   *
   * \param recvBuffer [out] The location for storing the received data.
   *
   * \param sourceRank [in] The rank of the process to receive the data from.
   * If <tt>sourceRank < 0</tt> then data will be received from any process.
   *
   * <b>Preconditions:</b><ul>
   * <li>[<tt>sourceRank >= 0] <tt>sourceRank < this->getSize()</tt>
   * <li><tt>sourceRank != this->getRank()</tt>
   * </ul>
   *
   * \return Returns the senders rank.
   */
  virtual RCP<CommRequest<Ordinal> > ireceive(
    const ArrayView<char> &recvBuffer,
    const int sourceRank
    ) const = 0;

  //! Variant of ireceive that takes a tag.
  virtual RCP<CommRequest<Ordinal> >
  ireceive (const ArrayView<char> &recvBuffer,
            const int sourceRank,
            const int tag) const = 0;

  /** \brief Wait on a set of communication requests.
   *
   * <b>Preconditions:</b><ul>
   * <li> <tt>requests.size() > 0</tt>
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li> <tt>is_null(request[i]))</tt> for <tt>i=0...requests.size()-1</tt>
   * </ul>
   */
  virtual void waitAll(
    const ArrayView<RCP<CommRequest<Ordinal> > > &requests
    ) const = 0;

  /// \brief Wait on communication requests, and return their statuses.
  ///
  /// \pre requests.size() == statuses.size()
  ///
  /// \pre For i in 0, 1, ..., requests.size()-1, requests[i] is
  ///   either null or requests[i] was returned by an ireceive() or
  ///   isend().
  ///
  /// \post For i in 0, 1, ..., requests.size()-1,
  ///   requests[i].is_null() is true.
  ///
  /// \param requests [in/out] On input: the requests on which to
  ///   wait.  On output: all set to null.
  ///
  /// \param statuses [out] The status results of waiting on the
  ///   requests.
  virtual void
  waitAll (const ArrayView<RCP<CommRequest<Ordinal> > >& requests,
           const ArrayView<RCP<CommStatus<Ordinal> > >& statuses) const = 0;

  /// \brief Wait on a single communication request, and return its status.
  ///
  /// \param request [in/out] On input: request is not null, and
  /// <tt>*request</tt> is either null (in which case this function
  /// does nothing and returns null) or an RCP of a valid CommRequest
  /// instance representing an outstanding communication request.  On
  /// output: If the communication request completed successfully, we
  /// set <tt>*request</tt> to null, indicating that the request has
  /// completed.  (This helps prevent common bugs like trying to
  /// complete the same request twice.)
  ///
  /// \return If *request is null, this method returns null.
  /// Otherwise this method returns a CommStatus instance representing
  /// the result of completing the request.  In the case of a
  /// nonblocking receive request, you can query the CommStatus
  /// instance for the process ID of the sending process.  (This is
  /// useful for receiving from any process via \c MPI_ANY_SOURCE.)
  ///
  /// \pre <tt>!is_null(request)</tt> (that is, the Ptr is not null).
  /// \post <tt>is_null(*request)</tt> (that is, the RCP is null).
  ///
  /// This function blocks until the communication operation
  /// associated with the CommRequest object has completed.
  virtual RCP<CommStatus<Ordinal> >
  wait (const Ptr<RCP<CommRequest<Ordinal> > >& request) const = 0;

  //@}

  //! @name Subcommunicator Operations
  //@{

  /**
   * \brief Duplicate this communicator.
   *
   * Make a copy of this communicator with a duplicate communication
   * space.  Note that the returned communicator has the same
   * properties (including process ranks, attributes and topologies)
   * as this communicator, but is distinct from the original.
   * "Distinct" means that if you send a message on the original
   * communicator, you can't receive it on the new one, and vice
   * versa.  The new communicator represents a separate message space.
   * This has the same semantics as MPI_Comm_dup.  (In fact, the
   * subclass MpiComm implements this using MPI_Comm_dup.)
   *
   * Most users don't want to do this.  The duplicate() method returns
   * a <i>new communicator</i>.  In MPI terms, it is a <i>different
   * MPI_Comm</i>.  If you want a shallow copy of the handle, you
   * should pass the <tt>Comm<Ordinal><tt> around by const pointer,
   * like this:
   * \code
   * RCP<const Comm<int> > comm = ...; // my original communicator
   * // ... do some stuff with comm ...
   * // Make a shallow copy.
   * RCP<const Comm<int> > diffHandleSameComm = comm;
   * // ... do some stuff with diffHandleSameComm ...
   * \endcode
   * This behaves the same as the following "raw MPI" code:
   * \code
   * MPI_Comm comm = ...; // my original communicator
   * // ... do some stuff with comm ...
   * // Make a shallow copy.
   * MPI_Comm diffHandleSameComm = comm;
   * // ... do some stuff with diffHandleSameComm ...
   * \endcode
   * The subclass of Comm ensures that the "raw" MPI handle is freed
   * only after the last reference to it by a subclass instance
   * disappears.  (It does reference counting underneath.)
   *
   * Please, please do not invoke the copy constructor or assignment
   * operator of Comm.  Of course it's not legal to do that anyway,
   * because Comm is pure virtual.  However, even if you could do it,
   * you must never do this!  For example, do <i>not</i> do this:
   * \code
   * RCP<const Comm<int> > comm = ...; // my original communicator
   * // ... do some stuff with comm ...
   * // DO NOT DO THIS, EVER!!!  THIS IS VERY BAD!!!
   * RCP<const Comm<int> > badComm (new Comm<int> (*comm));
   * \endcode
   * and do <i>not</i> do this:
   * \code
   * RCP<const Comm<int> > comm = ...; // my original communicator
   * // ... do some stuff with comm ...
   * // DO NOT DO THIS, EITHER!!!  THIS IS JUST AS BAD!!!
   * RCP<const Comm<int> > badComm = rcp (new Comm<int> (*comm));
   * \endcode
   * This is bad because it ignores the subclass' data.  Depending on
   * the subclass of Comm that you are actually using, it may be
   * appropriate to invoke the copy constructor or assignment operator
   * of the specific subclass, but <i>never</i> those of Comm itself.
   *
   * Users are not responsible for freeing the returned communicator.
   * The destructor of the subclass of Comm handles that itself.
   *
   * In an MPI implementation, the returned communicator is created
   * using MPI_Comm_dup, with the resulting semantic implications.
   *
   * \return A new communicator.
   */
  virtual RCP< Comm > duplicate() const = 0;

  /**
   * \brief Split a communicator into subcommunicators based on color
   * and key.
   *
   * Partition this communicator into multiple disjoint groups, and
   * return the communicator corresponding to the group to which this
   * process belongs.  There will be as many groups as there are
   * globally many distinct values for the <tt>color</tt> parameter.
   * Within each subset of the partition, the ranks will be ordered
   * according to the key value each process passed for the
   * <tt>key</tt> parameter. If multiple processes pass the same value
   * for <tt>key</tt>, then they will be ordered according to their
   * rank in the original communicator.  To return a valid
   * communicator, this function requires a nonnegative value for
   * <tt>color</tt>.  If <tt>color</tt> is negative, this method will
   * return a null communicator.
   *
   * This method must be called as a collective on all processes in
   * this communicator.  That is, if this method is called at all, it
   * must be called on all processes in the communicator.
   *
   * Users are not responsible for freeing the returned communicator.
   * The destructor of the subclass of Comm handles that itself.
   *
   * In an MPI implementation, the returned communicator is created
   * using MPI_Comm_split, with the resulting semantic implications.
   *
   * \param color [in] An integer representing the color for the local
   *   rank.  In the MPI implementation, if this is negative,
   *   MPI_Comm_split gets <tt>MPI_UNDEFINED</tt> as the color.
   *
   * \param key [in] A key value to order processes of the same color.
   *   In the MPI implementation, this is passed directly to
   *   MPI_Comm_split.
   *
   * \return A partitioned communicator.
   */
  virtual RCP<Comm> split (const int color, const int key) const = 0;

  /**
   * \brief Create a subcommunicator containing the specified processes.
   *
   * Create and return a subcommunicator of this communicator.  The
   * subcommunicator contains the processes in this communicator with
   * the given ranks, in which they are listed in the input vector.
   * Processes whose ranks are not included in the input vector will
   * be given a null communicator.
   *
   * This method must be called as a collective on all processes in
   * this communicator.  That is, if this method is called at all, it
   * must be called on all processes in the communicator.
   *
   * Users are not responsible for freeing the returned communicator.
   * The destructor of the subclass of Comm handles that itself.
   *
   * In an MPI implementation, the subcommunicator is created using
   * MPI_Comm_create, with the resulting semantic implications.
   *
   * \param ranks The ranks of the processes to include in the subcommunicator.
   * \return The subcommunicator.
   */
  virtual RCP<Comm>
  createSubcommunicator (const ArrayView<const int>& ranks) const = 0;
  //@}

}; // class Comm

} // namespace Teuchos

#endif // TEUCHOS_COMM_HPP
