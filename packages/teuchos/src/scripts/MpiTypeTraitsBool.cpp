// C++ bool is problematic to MPI because MPI's C++ bindings are
// deprecated (as of MPI 2.2).  This means that we can't rely on a
// native standard MPI_Datatype (in this case, MPI::BOOL) for bool.
// It also means that we don't get standard reduction operators for
// bool; we have to define our own custom operators, or hope that one
// of the standard operators for an integer type of the same size
// works correctly.  Finally, std::vector<bool> is not at all an array
// of bool; it's a bit set.  Users who attempt to pass the address of
// the first entry of an std::vector<bool> into an MPI function will
// experience undefined and likely undesirable behavior.
//
// Note that MPI_C_BOOL is for C99 '_Bool', which is not guaranteed to
// have the same representation (or even the same sizeof value) as C++
// 'bool'.  The 2003 version of the C++ standard (Section 5.3.3/1)
// says that sizeof(bool) is implementation defined, whereas
// sizeof(char) == sizeof(signed char) == sizeof(unsigned char) == 1.
// (There are legitimate performance reasons to want bool to be as big
// as a word; some CPUs are inefficient at extracting and using data
// smaller than a word.) 
//
// We combine two approaches in order to ensure that MpiTypeTraits 
// is correct on all platforms. 
//
// 1. Find the integer type of the same length as bool.
//
// It's reasonable to assume that bool has the same length as some
// kind of standard integer type.  If it doesn't, we fall back to the
// second approach:
//
// 2. Create a derived datatype, using sizeof(bool) contiguous MPI_CHAR.
//
// This works because bool carries information (so sizeof(bool) > 0),
// and sizeof(T) for any type T must be a nonnegative integer.  Thus,
// sizeof(bool) must be a multiple of sizeof(char) (which == 1 by the
// 2003 version of the C++ standard).

namespace { // anonymous
  //! Fall-back for returning a derived MPI_Datatype for Packet=bool.
  static MPI_Datatype makeDerivedTypeForBool () {
    int err = 0;
    MPI_Datatype t;
    // This is always guaranteed to work.  sizeof(T) for any type T
    // must be a nonnegative integer, and thus must be a multiple of
    // sizeof(char) (which == 1 by the C++03 standard).
    err = MPI_Type_contiguous (sizeof(bool), MPI_CHAR, &t);
    if (err != MPI_SUCCESS) {
      TEUCHOS_TEST_FOR_EXCEPTION(err == MPI_ERR_INTERN, std::runtime_error, 
        "MpiTypeTraits<bool>: MPI_Type_contiguous returned MPI_ERR_INTERN.  "
        "This means that MPI failed to acquire memory for allocating the "
        "custom MPI_Datatype.");
      TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::logic_error, 
        "MpiTypeTraits<bool>: MPI_Type_contiguous failed, "
        "probably due to Teuchos programmer error.  "
	"Please report this to the Teuchos developers.");
    } 
    err = MPI_Type_commit (&t);
    TEUCHOS_TEST_FOR_EXCEPTION(err != MPI_SUCCESS, std::runtime_error, 
      "MpiTypeTraits<bool>: MPI_Type_commit failed.  This may mean either "
      "Teuchos programmer error or something wrong with the current state of "
      "MPI.");
    return t;
  }
} // namespace (anonymous)

template<>
MPI_Datatype MpiTypeTraits<bool>::getType () {
  // The compiler should be able to optimize away the false branches.
  if (sizeof(bool) == sizeof(char)) {
    return MPI_CHAR;
  } else if (sizeof(bool) == sizeof(short)) {
    return MPI_SHORT;
  } else if (sizeof(bool) == sizeof(int)) {
    return MPI_INT;
  } else if (sizeof(bool) == sizeof(long)) {
    return MPI_LONG;
  } else { // Fall-back if bool isn't the same size as a common integer type.
    return makeDerivedTypeForBool ();
  }
}

template<>
size_t MpiTypeTraits<bool>::getCount () {
  return static_cast<size_t> (1);
}
