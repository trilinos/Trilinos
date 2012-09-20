/// \class MpiTypeTraits
/// \brief Traits class mapping a given \c Packet type to its MPI_Datatype.
/// \tparam Packet The type of data being sent and received.
///
/// \section Teuchos_MpiTypeTraits_Summary Summary
///
/// For a given C++ type Packet, this traits class tells you:
/// 1. The corresponding MPI_Datatype and count (i.e., number of
///    contiguous instances of that MPI_Datatype).  These tell MPI how
///    to send and receive objects of type \c Packet.
/// 2. Whether the MPI_Datatype is a <i>basic</i> or <i>derived</i>
///    type.  Derived types must be freed after use, by calling
///    MPI_Type_free().  Basic types need not and must not be freed
///    after use.
/// 3. Whether the MPI_Datatype and count are valid for all instances
///    of Packet on an MPI process, at all times in the program, on
///    all MPI processes.  (These are three separate conditions.  See
///    Section \ref Teuchos_MpiTypeTraits_Three ("Three consistency
///    conditions") for details.)
///
/// This class is meant mainly for Trilinos packages such as Tpetra
/// and Zoltan2, that use Comm (Teuchos' MPI wrapper) to communicate
/// data.  It is also appropriate for any direct user of Comm.  If you
/// do not use Comm directly, you need not worry about this class.
///
/// This class works by specialization.  If there is no specialization
/// of MpiTypeTraits for a type \c Packet, then expressions involving
/// <tt> MpiTypeTraits<Packet> </tt> will not compile.  We provide
/// specializations for many commonly used C++ types, including all
/// standard built-in integer and floating-point types.  However, you
/// must write a specialization for any \c Packet types we do not
/// cover.  Be sure to check this header file for the list of
/// specializations that we provide.  Duplicating an existing
/// specialization will result in an error either at compile time or
/// link time.
///
/// \section Teuchos_MpiTypeTraits_Prereq Prerequisites
///
/// In order to understand this documentation, you must first be
/// familiar with MPI (the Message Passing Interface).  In particular,
/// you should be aware that MPI uses an MPI_Datatype and a count to
/// refer to the type and amount of data being sent, received, or
/// accumulated in a communication operation.  You should understand
/// that some MPI_Datatypes ("derived" ones) need to be freed after
/// use by calling MPI_Type_free(), and others ("basic" ones, which
/// come predefined by the MPI standard) need not and must not be
/// freed in this way.  Finally, you should know what an opaque handle
/// is in the context of MPI, that MPI_Datatype is an example of an
/// opaque handle, and that you can use Teuchos' OpaqueWrapper to
/// handle opaque handles safely.
///
/// \section Teuchos_MpiTypeTraits_What What is an MPI_Datatype?
///
/// MPI_Datatype is an <i>opaque handle</i> type.  It is typically
/// implemented either as a pointer to a struct (whose definition you
/// can't see or shouldn't rely upon), or as an integer index into a
/// table.  For more details on opaque handles and how to manage them,
/// please see the documentation of OpaqueHandle.  We recommend using
/// OpaqueHandle to manage the MPI_Datatype returned by makeType().
///
/// \subsection Teuchos_MpiTypeTraits_What_DontIgnore Ignore MPI_Datatype at your own peril
///
/// Let me say this again: <i>Ignore MPI_Datatype at your own
/// peril</i>.  There are at least three reasons why you should always
/// use the right MPI_Datatype for your data.  Two of them relate to
/// correctness and availability of functionality, and one relates to
/// performance.
///
/// 1. The MPI standard only requires that one-sided communication
///    (e.g., \c MPI_Put, \c MPI_Get, and \c MPI_Accumulate) be
///    defined for basic MPI_Datatype instances.  They need not
///    necessarily work for derived MPI_Datatypes (those which are not
///    built in, which must be freed after use by calling
///    MPI_Type_free()).
/// 2. The MPI standard allows certain reduction operations (e.g.,
///    MPI_Allreduce) to produce incorrect results when using the
///    wrong datatype (e.g., \c MPI_CHAR instead of MPI_DOUBLE for the
///    C++ built-in type \c double), <i>even if</i> the total length
///    of the data is correct.  (This can happen in practice, and in
///    fact was the motivation for redesigning Teuchos' MPI
///    interface.)
/// 3. Using the right MPI_Datatype exposes optimizations to the
///    underlying communication hardware.  For example, reductions
///    with certain basic MPI_Datatypes may be performed in the
///    network hardware, without involving the CPU at all.  This may
///    significantly reduce latency of certain collective operations.
///
/// \subsection Teuchos_MpiTypeTraits_What_NotSer Not just for serialization
///
/// Note that MPI_Datatype is not merely a description of how to
/// serialize an object.  MPI does provide serialization (\c MPI_Pack)
/// and deserialization (\c MPI_Unpack) functions which use
/// MPI_Datatype in this way, but an MPI implementation need not
/// necessarily serialize data in order to send or receive it.  For
/// example, sending an array of \c double between two processes on
/// the same shared-memory node need only copy memory.
///
/// \section Teuchos_MpiTypeTraits_Three Three consistency conditions
///
/// The MPI_Datatype and count corresponding to many built-in C++
/// <tt>Packet</tt> types, like \c int or \c double, have the pleasant
/// property that they are valid for
/// 1. all <i>instances</i> of Packet on an MPI process,
/// 2. at all <i>times</i> in an execution of the program,
/// 3. on all <i>MPI processes</i>.
///
/// These conditions do not hold for all \c Packet types.  This
/// matters because in order for MPI to send, receive, or perform
/// collectives correctly on objects of type \c Packet, the
/// MPI_Datatypes and counts must be consistent on all participating
/// processes.  (Note that we didn't say "the same," only
/// "consistent."  It's possible to have different MPI_Datatypes on
/// the sender and receiver.  This is useful for doing things like
/// transposing a matrix.  Read an MPI book for more examples.)
///
/// The problem is that the receiving process needs to know how many
/// data there are, and how to receive them.  If the MPI_Datatype and
/// count are consistent on both sender and receiver, the sender and
/// receiver need not coordinate, other than to do the actual send or
/// receive.  Otherwise, correct communication requires coordination.
/// For example, if the counts differ on each end but the MPI_Datatype
/// instances on each end are consistent, then one of the following
/// must take place.  Either
/// - the sending process sends two messages, first the count and then
///   the actual object, or
/// - the sender sends without coordination, the receiving process
///   probes (via MPI_Probe) to figure out the count, and then it
///   receives the message.
///
/// The same problem can happen if the MPI_Datatypes are not
/// consistent, even if the counts are.  Both of these coordination
/// options cost convenience and performance.  Thus, it's important to
/// know in advance whether coordination is necessary.  Furthermore,
/// coordination affects the ability to do nonblocking communication.
/// For example, since an MPI_Irecv can't start without the message
/// count, the receiver must block on the message containing the count
/// before it can start the MPI_Irecv.
///
/// \subsection Teuchos_MpiTypeTraits_Three_Ex Examples violating the consistency conditions
///
/// This section describes several examples of useful \c Packet types
/// which violate one or more of the consistency conditions.
///
/// The C++ Standard Library class <tt>std::string</tt> violates all
/// three consistency conditions.  Different strings may have
/// different lengths (and thus different counts of <tt>MPI_CHAR</tt>)
/// (#1), the length of an individual <tt>std::string</tt> instance
/// may change (#2), and therefore of course, different strings may
/// differ in length on different processes (#3).  We see that in
/// general, violation of either #1 or #2 implies violation of #3.
///
/// Each instance of the arbitrary-precision real type provided by the
/// MPFR library may have a different precision.  The precision is set
/// in the MPFR object's constructor.  Therefore, each MPFR instance
/// may have a different size, which violates #1, #2, and therefore
/// #3.  Users could, however, choose to make all MPFR instances the
/// same precision (#1) at all times (#2) on all processes (#3).  This
/// assertion would remove the need for coordination between sending
/// and receiving processes.  The right place to make this assertion
/// is not here, in the traits class, nor is it in Teuchos MPI wrapper
/// Comm or Comm's nonmember helper functions.  It's wherever those
/// MPI wrappers get used to send and receive \c Packet data, for
/// example in Tpetra::Distributor.
///
/// The automatic differentiation types provided by the Sacado package
/// in Trilinos may violate all three conditions, just like MPFR
/// arbitrary-precision types.
///
/// The ARPREC library's arbitrary-precision real type obeys different
/// rules than MPFR.  Each ARPREC instance on a process has the same
/// precision.  However, this precision is a global per-process
/// parameter which may be changed at any time, and which may be
/// different on each MPI process.  This violates Conditions 2 and 3,
/// but not necessarily Condition 1.  Furthermore, users could promise
/// not to change the precision, which would allow the same
/// optimization as mentioned above for MPFR.
///
/// \subsection Teuchos_MpiTypeTraits_Three_NoPath We forbid pathological heterogeneity
///
/// The above discussion excludes certain pathological cases relating
/// to heterogeneous processors or compilers.  For example, our
/// provided specializations of MpiTypeTraits assume that on all MPI
/// processes, for all built-in C++ \c Packet types,
/// <tt>sizeof(Packet)</tt> has the same value.  (This excludes
/// built-in C++ types which have a size guaranteed by the language
/// standard, such as \c char, \c float, and \c double.)
///
/// It is technically possible to run on clusters with heterogeneous
/// nodes, such that <tt>sizeof(Packet)</tt> for types such as \c
/// bool, \c int, or \c size_t may have different values on different
/// nodes.  Some MPI implementations, such as OpenMPI, even support
/// this case correctly via the XDR back-end communication standard.
/// However, Teuchos and its client packages probably do not work
/// correctly for this use case, so we do not support it here.
///
/// We make no effort to check whether different processes have
/// different values of <tt>sizeof(Packet)</tt> for all built-in C++
/// \c Packet types.  If you find yourself running on a platform with
/// heterogeneous sizes of some built-in types, you will have to
/// rewrite the specializations of MpiTypeTraits for the affected
/// types.  Their \c globallyConsistent and \c globallyConsistentType
/// fields must be \c false, and you must reimplement their packing
/// methods (e.g., pack() and unpack()) to handle this case.  You must
/// also make sure that your MPI implementation can handle it as well.
///
/// \subsection Teuchos_MpiTypeTraits_Three_Use How to use consistency conditions to write generic code
///
/// This traits class tells you enough about which of the above three
/// consistency conditions are true that you may write correct,
/// efficient generic code.  <i>Generic</i> means that you do not have
/// to know what type \c Packet is.  We have distilled the consistency
/// conditions into two Boolean fields:
/// - <tt>globallyConsistent</tt>: If true, then \c Packet satisfies
///   all three consistency conditions.
/// - <tt>globallyConsistentType</tt>: If true, then all
///   <tt>Packet</tt> instances, at all times during execution of a
///   program, on all MPI processes, have a consistent MPI_Datatype.
///   Different \c Packet instances may still have different sizes,
///   depending on whether \c globallyConsistent is true.
///
/// If \c globallyConsistent is true, then the MPI_Datatype and count
/// returned by makeType() may be used to communicate \c Packet
/// instances, without prior coordination of count information.  The
/// sending process need not send the count first in a separate
/// message, and the receiving process need not probe the size.
///
/// If \c globallyConsistentType is true, then MPI processes need only
/// coordinate on the count, not the MPI_Datatype.  Thus, for sending
/// or receiving a \c Packet array, it suffices to use the
/// MPI_Datatype of the first instance in the array.  If users know
/// that the count for a \c Packet instance won't change, they should
/// save the count whenever possible for reuse, so as to avoid the
/// extra coordination step.
///
/// If \c globallyConsistentType is false, then different entries of
/// an array of \c Packet may have inconsistent MPI_Datatypes.  Thus
/// you cannot send the array directly using a single MPI_Send call,
/// or receive it using a single MPI_Recv call.  A good way to handle
/// this case would be to use MPI_Pack on the sending process to pack
/// each entry in turn into a send buffer of MPI_PACKED, and and
/// MPI_Unpack on the receiving process to unpack the entries.
///
/// We further subdivide the <tt>globallyConsistentType == false</tt>
/// case into the following dichotomy:
/// 1. Knowing the MPI_Datatype and count before a receive suffices
///    for preparing the \c Packet (e.g., resizing) for the receive.
/// 2. Knowing the MPI_Datatype and count before a receive does not
///    suffice.  The received data actually has to come in before 
///    the \c Packet can be prepared.
/// 
/// Case 2 is unusual for \c Packet types that most users want to send
/// or receive.  It might happen if the data must be unpacked using
/// MPI_Unpack(), and the initially unpacked data determine the
/// <tt>Packet</tt>'s contents.  For example, one could encode
/// <tt>Teuchos::any</tt> by packing it (using MPI_Pack()) with an
/// initial header containing type info, which the receiving process
/// must unpack and read before it knows what C++ type lives inside
/// the <tt>Teuchos::any</tt>.  (We don't think
/// <tt>Packet=Teuchos::any</tt> this is a good idea, but it is
/// technically possible.)  Case 1 is by far the most common for the
/// kinds of C++ types that Trilinos users might want to communicate
/// between MPI processes.
///
/// \section Teuchos_MpiTypeTraits_Specialize How to write a specialization
///
/// This header file (Teuchos_MpiTypeTraits.hpp) and the source file
/// Teuchos_MpiTypeTraits.cpp provide many examples of how to
/// specialize this class for specific \c Packet types.  Examples
/// include both built-in (basic) MPI_Datatypes, and derived
/// MPI_Datatypes.  You might find yourself writing a specialization
/// that produces a basic MPI_Datatype if your MPI implementation
/// provides useful basic MPI_Datatypes beyond what the standard
/// guarantees, or if later versions of the MPI standard add new basic
/// MPI_Datatypes.
///
/// The values of \c needFree, \c globallyConsistent, and
/// <tt>globallyConsistentType</tt> in your specialization must be
/// compile-time constants whose values are visible in the header
/// file.  This lets them be used as template parameters.  It also
/// ensures that their values are the same on all MPI processes.
///
/// Here is an example specialization for <tt>std::string</tt>:
/// \code
/// template<>
/// class MpiTypeTraits<std::string> {
///   // Don't use MPI_Type_free(); MPI_CHAR is a basic MPI_Datatype.
///   static const bool needFree = false; 
///   // Different strings may have different lengths.
///   // We have to ask each string for its length at a particular time.
///   static const bool globallyConsistent = false; 
///   // All strings at all times on all MPI processes use MPI_CHAR.
///   static const bool globallyConsistentType = true; 
///   std::pair<MPI_Datatype, size_t> makeType (const std::string& s) {
///     return std::make_pair (MPI_CHAR, s.size ());
///   }
///   void 
///   resizeForReceive (std::string& packet, 
///                     const std::pair<MPI_Datatype, size_t> typeInfo) 
///   {
///     packet.resize (typeInfo.second);
///   }
/// };
/// \endcode
///
/// Note that specializing this class does not make sense for pointer
/// types.  These include "raw" pointers (<tt>T*</tt> for a C++ type
/// \c T), "smart" pointers such as RCP, iterators, and elements of
/// linked data structures with possible cycles (such as vertices in a
/// graph).
///
/// \section Teuchos_MpiTypeTraits_Design Design discussion
///
/// This section is likely only of interest to implementers of
/// specializations of MpiTypeTraits, or to those who have questions
/// about the design choices made by the authors of MpiTypeTraits.
///
/// \subsection Teuchos_MpiTypeTraits_Design_Datatypes Notes on derived datatypes
///
/// The <a
/// href="http://www.mpi-forum.org/docs/mpi22-report/node78.htm#Node78">MPI
/// 2.2 standard</a> says that MPI_Type_free "[m]arks the datatype
/// object associated with datatype for deallocation and sets datatype
/// to MPI_DATATYPE_NULL. Any communication that is currently using
/// this datatype will complete normally.  Freeing a datatype does not
/// affect any other datatype that was built from the freed
/// datatype. The system behaves as if input datatype arguments to
/// derived datatype constructors are passed by value."  We may
/// conclude the following:
/// 1. We may call MPI_Type_free on any MPI_Datatype, even if a
///    nonblocking operation (e.g., an MPI_Irecv) with that datatype
///    hasn't yet completed.
/// 2. We may safely create composite datatypes (e.g.,
///    <tt>std::pair<X,Y></tt>) from previously created derived
///    datatypes.
///
/// The first statement means that we may free an MPI_Datatype any
/// time after using it in an MPI function call.  We need not track
/// MPI_Irecv or any other nonblocking operation (such as the new
/// nonblocking collectives in MPI 3.0).
///
/// The second statement means that a composite MPI_Datatype will work
/// correctly even if its component MPI_Datatypes are freed later.
/// Furthermore, we may free the composite MPI_Datatype either before
/// or after its component MPI_Datatypes are freed.  Thus, we may free
/// an MPI_Datatype any time after using it.  Otherwise, we could
/// never just create an MPI_Datatype; we would have to maintain a
/// table of all MPI_Datatypes we ever needed to create, and manage
/// use of them by wrapping each in an
/// <tt>RCP<OpaqueWrapper<MPI_Datatype> ></tt> with an elaborate tree
/// of references, so that each parent MPI_Datatype persists at least
/// as long as its children.
///
/// \subsection Teuchos_MpiTypeTraits_Design_Raw Why return a "raw" MPI_Datatype?
///
/// Users might wonder why makeType() returns a "raw" MPI_Datatype
/// rather than a "wrapped" <tt>RCP<OpaqueWrapper<MPI_Datatype> ></tt>
/// with a suitable deallocator.  That would be helpful for creating
/// composite datatypes; for example, it would make it easier to write
/// exception-safe code.  We chose to return the "raw" MPI_Datatype
/// first, for the following reasons:
/// 1. Wrapping it adds overhead for the overwhelmingly common case of
///    basic (not derived) datatypes.
/// 2. We want to give clients of this class a choice about their
///    deallocation mechanism.
///
/// Users who want their MPI_Datatypes wrapped in this way may use the
/// nonmember template function makeMpiType().
///
/// \subsection Teuchos_MpiTypeTraits_Design_Instance Why does makeType() need a \c Packet input?
///
/// In general, for maximum portability, one may need a Packet in
/// order to compute its MPI_Datatype, even if \c
/// <tt>globallyConsistentType == true</tt>.  This is because C++ does
/// not provide a way to perform introspection on the layout of a
/// class or struct, without actually having an instance of it.
///
/// The makeType() function could just assume that \c Packet is
/// default constructible, and compute the MPI_Datatype and size of a
/// default constructed instance.  However, this is wrong for at least
/// two reasons:
/// 1. If different instances of \c Packet may differ in size, then
///    the size of a default constructed object may not be meaningful.
///    A good example is <tt>std::string</tt>.
/// 2. It requires that \c Packet be default constructible.  
/// 
/// \subsection Teuchos_MpiTypeTraits_Design_Assumptions What this class assumes about \c Packet
///
/// 1. <tt>Packet</tt>'s assignment operator (<tt>operator=()</tt>)
///    does the right thing.  We don't need to include an extra
///    <tt>copy()</tt> function here.
/// 2. \c Packet is not immutable.  \c Packet instances may be modified
///    after their construction.  (MPI assumes this anyway.)
template<class Packet>
class MpiTypeTraits {
 public:
  //! The type of the data being received and sent.
  typedef Packet packet_type;
  
  /// \brief Whether the MPI_Datatype returned by makeType() must be freed after use.
  ///
  /// Derived MPI_Datatype instances (which were created using one of
  /// MPI's type constructor functions) must be freed after use by
  /// calling MPI_Type_free().  Basic MPI_Datatype instances (which
  /// come predefined by MPI) must <i>not</i> be freed.  See the
  /// <a href="http://www.mcs.anl.gov/research/projects/mpi/www/www3/MPI_Type_free.html">documentation of MPI_Type_free()</a> 
  /// for details.
  ///
  /// \note This is a compile-time constant so that you may use it as
  ///   a template parameter.  It is visible in the header fie so that
  ///   users can be sure that its value is the same on all MPI
  ///   processes.
  static const bool needFree = false;

  /// Whether the MPI_Datatype and size returned by makeType() are the
  /// same for all \c Packet instances, for all MPI processes, at all
  /// times in the program's execution.
  ///
  /// See the class documentation for a detailed explanation.
  static const bool globallyConsistent = false;

  /// Whether all Packet instances on all MPI processes have a
  /// consistent MPI_Datatype (but possibly different sizes) at all
  /// times in the program's execution.
  ///
  /// See the class documentation for a detailed explanation.
  ///
  /// \note To Trilinos developers: If implementing a cache for
  ///   MPI_Datatype instances, you should not cache those
  ///   orresponding to \c Packet types for which this is false.
  static const bool globallyConsistentType = false;

  //! Return the MPI_Datatype and count corresponding to the Packet type.
  static std::pair<MPI_Datatype, size_t> makeType (const Packet& example) {
    (void) example; // Silence compiler warning for unused input.

    // Raise a compile-time error in case no specialization of this
    // traits class has been defined for the given Packet type.
    Packet::noSpecializationForThisType();

    // Return a valid value, so the compiler doesn't complain.
    return std::make_pair (MPI_TYPE_NULL, static_cast<size_t> (1));
  }

  /// \brief Resize \c packet, if necessary, for a receive.
  ///
  /// \param packet [in/out] The Packet instance to resize.
  /// \param typeInfo [in] The MPI_Datatype and size for the receive.
  ///
  /// The default implementation of this function does nothing.  This
  /// suffices for any \c Packet type for which \c globallyConsistent
  /// is true.  However, some \c Packet types require resizing in
  /// order to fit an incoming receive.  A good example is
  /// <tt>std::string</tt>.
  ///
  /// A type-generic receive implementation need only call this
  /// function if <tt>globallyConsistent</tt> is false.  However, this
  /// function must always be defined so that it is correct to call it
  /// if <tt>globallyConsistent</tt> is true.
  ///
  /// Note that \c typeInfo may not necessarily be the return value of
  /// <tt>makeType(packet)</tt> on input.  We have not received the
  /// packet yet, so we do not know yet what its size should be.
  /// We may not even know the incoming MPI_Datatype, if
  /// <tt>globallyConsistentType == false</tt>.  If indeed
  /// <tt>globallyConsistentType</tt> is false, then implementers may
  /// have to work backwards from the given MPI_Datatype and size to
  /// figure out how to resize the given packet.  For example, they
  /// might have to call MPI_Type_get_envelope() or
  /// MPI_Type_get_contents() to deconstruct the MPI_Datatype, detect
  /// and reject invalid MPI_Datatype input, and match valid
  /// MPI_Datatype input to \c Packet.  We consider this a case to be
  /// avoided.  Whenever possible, specializations of MpiTypeTraits
  /// should prefer a globally consistent MPI_Datatype whenever
  /// possible.
  static void resizeForReceive (Packet& packet, const std::pair<MPI_Datatype, size_t> typeInfo) {
    (void) packet; // Silence compiler warnings for unused arguments.
    (void) typeInfo;
  }
};

/// \fn makeMpiType
/// \brief Get MPI_Datatype (safely wrapped in OpaqueHandle) and size
///   for a \c Packet instance.
/// \tparam Packet The type of data being sent and received; same as
///   the \c Packet template parameter of MpiTypeTraits.
/// 
/// Please refer to the documentation of MpiTypeTraits for details.
template<class Packet>
std::pair<RCP<OpaqueHandle<MPI_Datatype> >, size_t>
makeMpiType (const Packet& p)
{
  std::pair<MPI_Datatype, size_t> result = MpiTypeTraits<Packet>::makeType (p);
  RCP<OpaqueHandle<MPI_Datatype> > handle;  
  if (MpiTypeTraits<Packet>::needFree) {
    handle = opaqueHandle (result.first, MPI_Type_free);
  } else { // It's a basic MPI_Datatype; don't call MPI_Type_free after use.
    handle = opaqueHandle (result.first);
  }
  return std::make_pair (handle, result.second);
}
