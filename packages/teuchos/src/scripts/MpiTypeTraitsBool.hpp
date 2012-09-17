/// \\brief Specialization of MpiTypeTraits for Packet type \c bool.
///
/// \\note For various reasons, we cannot promise that makeType() for
///   the Packet=bool specialization returns a basic datatype.  Even
///   if it does, we cannot promise that the MPI_Op reduction
///   functions for typical Boolean operations are not user-defined.
///   This means that bool should not in general be used in one-sided
///   communication operations, since we cannot guarantee this is
///   possible on all platforms.  In general, you should prefer using
///   C++ built-in integer types to encode Boolean values.  In most
///   cases, \\c int should suffice.  For large arrays of Boolean 
///   values, you might consider using \\c char to save space.
///
/// \warning std::vector<bool> is not an array of bool; it is a bit
///   set with a different representation than std::vector<T> for
///   other types T.  Users who attempt to pass the address of the
///   first entry of an std::vector<bool> into an MPI function will
///   experience undefined and likely undesirable behavior.
///
/// \warning We assume that on all MPI processes, sizeof(bool) has
///   the same value.  It is technically possible to run on clusters
///   with heterogeneous nodes, such that sizeof(bool) may have
///   different values on different nodes.  However, we do not support
///   this (currently uncommon) use case.
template<>
class MpiTypeTraits<bool> {
 public:
  typedef bool packet_type;
  // The right-hand side of this expression is a compile-time constant.
  static const bool needFree = sizeof(bool) != sizeof(char) &&
                               sizeof(bool) != sizeof(short) &&
                               sizeof(bool) != sizeof(int) &&
                               sizeof(bool) != sizeof(long);
  // We assume that every bool on every MPI process has the same size.
  static const bool globallyConsistent = true;
  static const bool globallyConsistentType = true;

  static std::pair<MPI_Datatype, size_t> makeType (const bool& x);
};
