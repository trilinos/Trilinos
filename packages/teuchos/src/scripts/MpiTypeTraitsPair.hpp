/// \brief Partial specialization for an std::pair of two types X and Y.
///
/// The implementation assumes that X and Y may be different types.
/// It makes no attempt to check whether they are the same type.  You
/// should provide specializations for this case yourself.
template<class X, class Y>
class MpiTypeTraits<std::pair<X, Y> > {
 public:
  typedef std::pair<X, Y> packet_type;
  static const bool needFree = true;
  static const bool globallyConsistent = 
    MpiTypeTraits<X>::globallyConsistent && 
    MpiTypeTraits<Y>::globallyConsistent;
  static const bool globallyConsistentType = 
    MpiTypeTraits<X>::globallyConsistentType && 
    MpiTypeTraits<Y>::globallyConsistentType;

  static std::pair<MPI_Datatype, size_t> 
  makeType (const std::pair<X, Y>& packet) 
  {
    int err = 0;
    MPI_Datatype t; // the return value
    int blkLens[2] = {1, 1}; // one each of X and Y
    MPI_Aint disps[2] = {&(packet.first) - &packet, &(packet.second) - &packet};

    // NOTE (mfh 30 Aug 2012) If we call makeType() for X and Y, we
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
    // so that makeType() could take in a datatype cache if one is
    // available, and otherwise construct datatypes from scratch.
    std::pair<MPI_Datatype, size_t> typeInfoX = MpiTypeTraits<X>::makeType ();
    RCP<OpaqueWrapper<MPI_Datatype> > dataTypeX;
    if (MpiTypeTraits<X>::needFree) {
      dataTypeX = opaqueWrapper (typeInfoX.first, MPI_Type_free);
    } else {
      dataTypeX = opaqueWrapper (typeInfoX.first);
    }
    std::pair<MPI_Datatype, size_t> typeInfoY = MpiTypeTraits<Y>::makeType ();
    RCP<OpaqueWrapper<MPI_Datatype> > dataTypeY;
    if (MpiTypeTraits<X>::needFree) {
      dataTypeY = opaqueWrapper (typeInfoY.first, MPI_Type_free);
    } else {
      dataTypeY = opaqueWrapper (typeInfoY.first);
    }

    MPI_Datatype types[2];
    types[0] = typeInfoX.first;
    types[1] = typeInfoY.first;

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
    return std::make_pair (t, static_cast<size_t> (1));
  }
};
