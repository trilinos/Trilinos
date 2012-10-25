/// \brief Partial specialization for an std::pair of two types X and Y.
///
/// The implementation assumes that X and Y may be different types.
/// It makes no attempt to check whether they are the same type.  You
/// should provide specializations for this case yourself.
template<class X, class Y>
class MpiTypeTraits<std::pair<X, Y> > {
 public:
  typedef std::pair<X, Y> packet_type;
  static const bool mustFree = true;
  static const bool sameDatatype = 
    MpiTypeTraits<X>::sameDatatype && MpiTypeTraits<Y>::sameDatatype;
  static const bool sameLocalCount = 
    MpiTypeTraits<X>::sameLocalCount && MpiTypeTraits<Y>::sameLocalCount;
  static const bool sameGlobalCount = 
    MpiTypeTraits<X>::sameGlobalCount && MpiTypeTraits<Y>::sameGlobalCount;
  static const bool direct = 
    MpiTypeTraits<X>::direct && MpiTypeTraits<Y>::direct;
  static const bool mustSerialize = 
    MpiTypeTraits<X>::mustSerialize && MpiTypeTraits<Y>::mustSerialize;

  static void* getPtr (const std::pair<X, Y>& packet) {
    // If either X or Y are indirect types, then there is no way to
    // get the pointer to the pair.  (A pair is really an array with a
    // constant size of two elements, so this is the generalization of
    // an array of indirect type.)  Thus, we just return the address.
    return reinterpret_cast<void*> (const_cast<Packet*> (&packet));
  }

  static size_t getCount (const std::pair<X, Y>& ) {
    return static_cast<size_t> (1);
  }

  static MPI_Datatype getType (const std::pair<X, Y>& packet) {
    // NOTE (mfh 30 Aug 2012) If we call getType() for X and Y, we
    // might be wasting some effort if X and Y are derived (instead of
    // basic) types.  However, this behavior is still correct; it does
    // not cause problems if X's or Y's type gets freed but the type
    // of their pair is not.  This comes from the MPI 2.2 standard,
    // in particular
    // <a href="http://www.mpi-forum.org/docs/mpi22-report/node78.htm#Node78">this section</a>.
    //
    // "Freeing a datatype does not affect any other datatype that was
    // built from the freed datatype.  The system behaves as if input
    // datatype arguments to derived datatype constructors are passed by
    // value."
    //
    // It's unlikely that our users will be constructing any composite
    // datatypes from datatypes that themselves are expensive to
    // compute.  If that does become the case, however, we could
    // refactor type creation to use a datatype cache that returns
    // RCP<OpaqueWrapper<MPI_Datatype> > for the datatype
    // corresponding to a given type T.  This could even be optional,
    // so that getType() could take in a datatype cache if one is
    // available, and otherwise construct datatypes from scratch.

    // "Subtypes" for X and Y.
    MPI_Datatype types[2] = {MPI_TYPE_NULL, MPI_TYPE_NULL};
    MPI_Datatype t = MPI_TYPE_NULL; // The return value.
    int err = 0; // Error code returned by MPI functions.

    try {
      int blkLens[2];
      blkLens[0] = static_cast<int> (countX);
      blkLens[1] = static_cast<int> (countY);
      MPI_Aint disps[2] = {&(packet.first) - &packet, &(packet.second) - &packet};

      // The operations that follow might raise an exception.  If they
      // do, we'll have to roll back by freeing the datatypes if needed.
      types[0] = MpiTypeTraits<X>::getType (packet.first);
      types[1] = MpiTypeTraits<Y>::getType (packet.second);
      const size_t countX = MpiTypeTraits<X>::getCount (packet.first);
      const size_t countY = MpiTypeTraits<Y>::getCount (packet.second);

      err = MPI_Type_create_struct (2, blkLens, disps, types, &t);
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, 
        "MpiTypeTraits<std::pair<" << TypeNameTraits<X>::name() << ", " 
        << TypeNameTraits<Y>::name() << ">: MPI_Type_create_struct failed.  "
        "This may mean either Teuchos programmer error or something wrong "
        "with the current state of MPI.");

      err = MPI_Type_commit (&t);
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, 
        "MpiTypeTraits<std::pair<" << TypeNameTraits<X>::name() << ", " 
        << TypeNameTraits<Y>::name() << ">: MPI_Type_commit failed.  "
        "This may mean either Teuchos programmer error or something wrong "
        "with the current state of MPI.");
    } 
    catch (...) {
      // Free the datatypes for X and Y if necessary.  If the above
      // code threw an exception, then MPI_Type_commit on t didn't
      // finish, so we don't have to free t.
      for (int i = 0; i < 2; ++i) {
	if (types[i] != MPI_TYPE_NULL) {
	  // Don't bother checking for errors.  We'll just try our best,
	  // since we are responding to an error anyway.
	  (void) MPI_Type_free (&types[i]);
	}
      }
      throw;
    }
    return t;
  }
};
