/// \class MpiTypeTraits
/// \brief Traits class mapping a given Packet type to its MPI_Datatype.
/// \tparam Packet The type of data being sent and received.
///
/// For a given C++ type Packet, this traits class tells you:
/// 1. The corresponding MPI_Datatype and count (i.e., number of
///    contiguous instances of that MPI_Datatype).  These tell MPI how
///    to send and receive objects of type Packet.
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
/// of MpiTypeTraits for a type Packet, then expressions involving
/// <tt> MpiTypeTraits<Packet> </tt> will not compile.  We provide
/// specializations for many commonly used C++ types, including all
/// standard built-in integer and floating-point types.  However, you
/// must write a specialization for any Packet types we do not
/// cover.  Be sure to check this header file for the list of
/// specializations that we provide.  Duplicating an existing
/// specialization will result in an error either at compile time or
/// link time.
template<class Packet>
class MpiTypeTraits {
 public:
  //! The type of the data being received and sent.
  typedef Packet packet_type;
  
  /// \brief Whether the MPI_Datatype returned by getType() must be freed after use.
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
  static const bool mustFreeDatatype = false;

  /// \brief Whether all instances of Packet have the same MPI_Datatype.
  ///
  /// If this is true, then all Packet instances have the same
  /// MPI_Datatype, at all times, on all processes.  This does not
  /// distinguish between the local and global cases.
  static const bool sameDatatype = false;

  /// \brief Whether all local instances (on the calling process) of Packet have the same count.
  ///
  /// If sameGlobalCount is true, then this must be true.  Vice versa
  /// is not necessarily true (ARPREC's arbitrary-precision real type
  /// is an example).
  static const bool sameLocalCount = false;

  /// \brief Whether all instances of Packet on all processes have the same count.
  ///
  /// If this is true, then the count is the same for all Packet
  /// instances, for all MPI processes, at all times in the program's
  /// execution.  In that case, sameLocalCount must also be true.
  /// Vice versa is not necessarily true (ARPREC's arbitrary-precision
  /// real type is an example).
  static const bool sameGlobalCount = false;

  /// \brief Whether a pointer to Packet is a valid MPI input or output buffer.
  ///
  /// You should always use getPtr() to get a pointer to the right
  /// data in a Packet instance.  The MPI_Datatype returned by
  /// getType() applies only to the pointer returned by getPtr().
  ///
  /// We say that Packet is "indirect" if direct is false.  An
  /// instance of an indirect types for which mustSerialize is false
  /// may be communicated directly, by giving the raw pointer returned
  /// by getPtr() directly to MPI as an input or output buffer.
  static const bool direct = false;

  /// \brief Whether you must serialize individual Packet instances before communicating them.
  /// 
  /// You must always serialize arrays of indirect types, as well as
  /// arrays of arrays of any type, regardless of the value of this
  /// parameter.
  static const bool mustSerialize = false;

  //! Pointer to give to MPI as an input or output buffer.
  static void* getPtr (const Packet& x) {
    (void) x; // Silence compiler warning for unused input.

    // Raise a compile-time error in case no specialization of this
    // traits class has been defined for the given Packet type.
    Packet::noSpecializationForThisType();

    // Return a legitimate pointer, so the compiler doesn't complain.
    return reinterpret_cast<void*> (const_cast<Packet*> (&x));
  }    

  /// \brief Return the count for the given Packet instance.
  ///
  /// If sameLocalCount is false, then different instances of Packet
  /// on the calling process may have different counts.  In that case,
  /// you must supply the specific Packet you want to send (or
  /// receive, if its layout has been coordinated with the sender), in
  /// order to get the correct count.  If sameLocalCount is true, you
  /// can use any Packet instance.  The implementation may or may not
  /// read it.  (For example: every \c double has MPI_Datatype
  /// MPI_DOUBLE.  An implementation of this when Packet is a simple
  /// struct might use MPI_Type_struct to build the MPI_Datatype, in
  /// which case it must take the address of different fields of the
  /// Packet.)
  static size_t getCount (const Packet& x) {
    (void) x; // Silence compiler warning for unused input.

    // Raise a compile-time error in case no specialization of this
    // traits class has been defined for the given Packet type.
    Packet::noSpecializationForThisType();

    // Return a legitimate count, so the compiler doesn't complain.
    return 1;
  }

  /// \brief Return the MPI_Datatype for the given Packet.
  ///
  /// This interface does allow different instances of Packet to have
  /// different datatypes.  For example, if you have a C++ type
  /// representing a structured grid, it would be appropriate to use
  /// an MPI_Datatype to represent an exterior ghost region being
  /// received, or an interior ghost region being sent.  
  ///
  /// It is possible that different instances of a scalar (a "single
  /// value" rather than an array, mesh, vector, matrix, or other
  /// container of values) Packet type might have different
  /// MPI_Datatype.  However, we recommend that you use serialization
  /// for such types, since that makes coordination of size and layout
  /// information much easier (just send the count).  If you use MPI's
  /// serialization, the corresponding MPI_Datatype is MPI_PACKED.
  ///
  /// If you are implementing makeType() for your own specialization
  /// of MpiTypeTraits, and if you plan not to read from the Packet
  /// input, you should take the usual measures to prevent compiler
  /// warnings for unused variables.
  static MPI_Datatype getType (const Packet& example) {
    (void) x; // Silence compiler warning for unused input.

    // Raise a compile-time error in case no specialization of this
    // traits class has been defined for the given Packet type.
    Packet::noSpecializationForThisType();

    // Return a legitimate MPI_Datatype, so the compiler doesn't complain.
    return MPI_TYPE_NULL;
  }
};

/// \brief Free an MPI_Datatype safely.
///
/// This function checks if MPI_Finalize has been called.  If it has,
/// it just returns.  Otherwise, it calls MPI_Type_free on the given
/// MPI_Datatype.  This ensures that no MPI functions will be called
/// after MPI_Finalize has been called.  
///
/// The current implementation of this function conceals a memory leak
/// in your code.  If you've allowed a derived MPI_Datatype to persist
/// past MPI_Finalize, then you are leaking memory.  We could fix this
/// by using MPI's caching mechanism to make the MPI_Datatype an
/// attribute of MPI_COMM_SELF, and adding a free function to the
/// attribute.  MPI_Finalize invokes the free functions of all the
/// attributes of MPI_COMM_SELF (as specified in the MPI 3.0
/// Standard), so this would free the MPI_Datatype on calling
/// MPI_Finalize.  We don't do this currently, because we would need a
/// name for the attribute that would not collide with any possible
/// user-specified attribute.
void TEUCHOS_LIB_DLL_EXPORT safeMpiDatatypeFree (MPI_Datatype datatype);

/// \fn getMpiDatatype
/// \brief Get MPI_Datatype (safely wrapped in OpaqueWrapper) for a Packet instance.
/// \tparam Packet The type of data being sent and received; same as
///   the Packet template parameter of MpiTypeTraits.
///
/// This function uses MpiTypeTraits<Packet> to get the MPI_Datatype
/// for communicating the given Packet instance.  Please refer to the
/// documentation of MpiTypeTraits for details.  We recommend using
/// this function if you do not want to manage freeing the
/// MPI_Datatype after use, as you would have to do if you used
/// MpiTypeTraits<Packet>::getType().
///
/// If MpiTypeTraits<Packet>::sameLocalCount is false, then different
/// Packet instances on the calling process may have different counts.
/// This means you need to give this function the specific Packet
/// instance you want to send or receive.  Otherwise, all Packet
/// instances on the calling process have the same count and
/// MPI_Datatype, so you can use any Packet instance.
///
/// \param p [in] The Packet instance for which to compute the MPI_Datatype.
///
/// \return The corresponding MPI_Datatype, wrapped in an RCP of
///   OpaqueWrapper.  The RCP will automatically calls MPI_Type_free if
///   necessary (and only if necessary) when the RCP's reference count
///   goes to zero.
template<class Packet>
RCP<OpaqueWrapper<MPI_Datatype> >
getMpiDatatype (const Packet& p)
{
  MPI_Datatype rawDatatype = MpiTypeTraits<Packet>::getType (p);
  RCP<OpaqueWrapper<MPI_Datatype> > handle;  
  if (MpiTypeTraits<Packet>::mustFreeDatatype) { 
    handle = opaqueWrapper (rawDatatype, safeMpiDatatypeFree);
  } else { // It's a basic MPI_Datatype; don't use MPI_Type_free.
    handle = opaqueWrapper (rawDatatype);
  }
  return handle;
}


/// \fn getMpiCount
/// \brief Get the count for communicating a Packet instance with MPI.
/// \tparam Packet The type of data being sent and received; same as
///   the Packet template parameter of MpiTypeTraits.
///
/// This function uses MpiTypeTraits<Packet> to get the count for
/// communicating the given Packet instance using MPI.  Please refer
/// to the documentation of MpiTypeTraits for details.  
///
/// If MpiTypeTraits<Packet>::sameLocalCount is false, then different
/// Packet instances on the calling process may have different counts.
/// This means you need to give this function the specific Packet
/// instance you want to send or receive.  Otherwise, all Packet
/// instances on the calling process have the same count, so you can
/// use any Packet instance.
///
/// \param p [in] The Packet instance for which to compute the count.
template<class Packet>
size_t getMpiCount (const Packet& p)
{
  return MpiTypeTraits<Packet>::getCount (p);
}
