// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Teuchos_Details_Allocator.hpp
/// \brief Declaration of Teuchos::Details::Allocator, a tracking and
///   logging implementation of the C++ Standard Library's Allocator
///   concept.
/// \warning This header file is an implementation detail of Teuchos.
///   We make no promises of backwards compatibility for this header
///   file or its contents.  They may change or disappear at any time.

#ifndef TEUCHOS_DETAILS_ALLOCATOR
#define TEUCHOS_DETAILS_ALLOCATOR

#include <Teuchos_ConfigDefs.hpp>
#include <iostream>
#include <limits>
#include <type_traits>
#include <typeinfo>

namespace Teuchos {
namespace Details {

/// \brief Logging implementation used by Allocator (see below).
///
/// \warning This is an implementation detail of Teuchos.  We make no
///   promises of backwards compatibility for this header file or this
///   class.  They may change or disappear at any time.
///
/// Please refer to
/// <tt> teuchos/core/test/HashTable/Allocator_atexit.cpp </tt>
/// for an example of how to register an atexit() hook that reports
/// current and maximum memory usage at exit from main().  It would
/// be easy to adapt this example to register an MPI_Finalize() hook,
/// using the MPI standard idiom of attaching an attribute to
/// MPI_COMM_SELF.
class AllocationLogger {
public:
  /// \brief Type of the size of an allocation or deallocation.
  ///
  /// This must match Allocator::size_type (see below).
  typedef std::size_t size_type;

  /// \brief Log an allocation.
  ///
  /// This is useful for Allocator<T>, but not so useful for users,
  /// unless they are writing a custom Allocator and want it use this
  /// class for logging.
  ///
  /// \param out [out] Output stream to which to write logging
  ///   information, if \c verbose is true.
  ///
  /// \param numEntries [in] Number of entries (of type \c T) in the
  ///   allocated array (1 if just allocating one instance of \c T).
  ///
  /// \param numBytes [in] Number of bytes in the allocated array.
  ///   This would normally be <tt>numEntries*sizeof(T)</tt>, at least
  ///   for types \c T which are "plain old data" (POD) or structs
  ///   thereof.  We make you compute this, so that the logger doesn't
  ///   have to know about (e.g., be templated on) the type \c T.
  ///
  /// \param typeName [in] Human-readable name of the type \c T, for
  ///   which you are logging an allocation.
  ///
  /// \param verbose [in] Whether to print logging information to the
  ///   given output stream \c out.
  static void
  logAllocation (std::ostream& out,
                 const size_type numEntries,
                 const size_type numBytes,
                 const char typeName[],
                 const bool verbose);

  /// \brief Log a deallocation, that was previously logged using
  ///   logAllocation().
  ///
  /// This is useful for Allocator<T>, but not so useful for users,
  /// unless they are writing a custom Allocator and want it use this
  /// class for logging.
  ///
  /// \param out [out] Output stream to which to write logging
  ///   information, if \c verbose is true.
  ///
  /// \param numEntries [in] Number of entries (of type \c T) in the
  ///   deallocated array (1 if just allocating one instance of \c T).
  ///   Note that the allocate() function in the C++ Standard
  ///   Library's Allocator concept gets this information, unlike
  ///   free() or <tt>operator delete[]</tt>.
  ///
  /// \param numBytes [in] Number of bytes in the deallocated array.
  ///   This would normally be <tt>numEntries*sizeof(T)</tt>, at least
  ///   for types \c T which are "plain old data" (POD) or structs
  ///   thereof.  We make you compute this, so that the logger doesn't
  ///   have to know about (e.g., be templated on) the type \c T.
  ///
  /// \param typeName [in] Human-readable name of the type \c T, for
  ///   which you are logging a deallocation.
  ///
  /// \param verbose [in] Whether to print logging information to the
  ///   given output stream \c out.
  static void
  logDeallocation (std::ostream& out,
                   const size_type numEntries,
                   const size_type numBytes,
                   const char typeName[],
                   const bool verbose);

  /// \brief Current total allocation in bytes.
  ///
  /// This count includes allocations using Allocator<T> for all \c T.
  static size_type curAllocInBytes ();

  /// \brief Max total allocation ("high water mark") in bytes.
  ///
  /// This count includes allocations using Allocator<T> for all \c T.
  static size_type maxAllocInBytes ();

  /// \brief Reset the current and max total allocation numbers to zero.
  ///
  /// This is mainly helpful just for tests.
  static void resetAllocationCounts ();

private:
  //! Current total allocation in bytes.
  static size_type curAllocInBytes_;

  //! Max total allocation ("high water mark") in bytes.
  static size_type maxAllocInBytes_;
};

/// \brief Optional tracking allocator for Teuchos Memory Management classes.
/// \tparam T The type of the entries in arrays to allocate.
///
/// \warning This is an implementation detail of Teuchos.  We make no
///   promises of backwards compatibility for this header file or this
///   class.  They may change or disappear at any time.
///
/// This class implements the C++ Standard Library's "Allocator"
/// concept.  It does "high water mark" tracking, and has the option
/// to print to \c stderr on allocations and deallocations.
///
/// \note To Teuchos developers: Making this class be the "allocator"
///   (second) template parameter of the std::vector type that
///   Teuchos::Array uses, would make Teuchos::Array track memory.
///   (This would require changing Teuchos::Array.)  It should also be
///   possible to change Teuchos::ArrayRCP to use this class as well,
///   but only to track memory that it allocates (vs. memory given it
///   by the user, possibly with a custom deallocator functor).
template<class T>
class Allocator {
private:
  /// \brief Internal enum, identifying whether an operation is an
  ///   allocation or a deallocation.
  enum EAllocatorOp {
    ALLOCATOR_ALLOCATE,
    ALLOCATOR_DEALLOCATE
  };

  //! Whether this Allocator logs.
  bool tracking () const { return track_; }

  //! Whether this allocator prints verbosely.
  bool verbose () const { return verbose_; }

  // This lets tracking() and verbose() stay private,
  // without breaking the templated copy constructor.
  template<class U>
  friend class Allocator;

public:
  //! Type of the template parameter of this class.
  typedef T value_type;

  /// \brief Type of a pointer to T.
  ///
  /// This is only needed for C++98, not for C++11 and newer.
  typedef T* pointer;
  /// \brief Type of a pointer to const T.
  ///
  /// This is only needed for C++98, not for C++11 and newer.
  typedef const T* const_pointer;
  /// \brief Type of a reference to T.
  ///
  /// This is only needed for C++98, not for C++11 and newer.
  typedef T& reference;
  /// \brief Type of a reference to const T.
  ///
  /// This is only needed for C++98, not for C++11 and newer.
  typedef const T& const_reference;

  /// \brief Type of the size of an allocation or deallocation.
  ///
  /// This must match AllocationLogger::size_type (see above).  If we
  /// didn't need this to have the same type as that, then we wouldn't
  /// need this typedef.
  typedef AllocationLogger::size_type size_type;

  /// \typedef difference_type
  /// \brief Integer type representing the difference between two pointers.
  ///
  /// This is only needed for C++98, not for C++11 and newer.
  /// However, we want this typedef to exist in both cases.  Thus, if
  /// C++11 is enabled, we use size_type above to compute this, in
  /// order to ensure consistency.
#ifdef HAVE_TEUCHOSCORE_CXX11
  typedef std::make_signed<size_type>::type difference_type;
#else
  typedef std::ptrdiff_t difference_type;
#endif // HAVE_TEUCHOSCORE_CXX11

  //! Default constructor.
  Allocator () :
    track_ (true), verbose_ (false)
  {}

  /// \brief Constructor.
  ///
  /// \param trackMemory [in] Whether to track memory usage at all.
  ///   If false, verboseOutput is ignored.
  /// \param verboseOutput [in] Whether to print on every allocation
  ///   and deallocation.  Even if this is true, trackMemory must also
  ///   be true in order for this to work.
  Allocator (const bool trackMemory,
             const bool verboseOutput) :
    track_ (trackMemory), verbose_ (verboseOutput)
  {}

  //! Copy constructor that takes an Allocator<U> for any U.
  template<class U>
  Allocator (const Allocator<U>& src) :
    track_ (src.tracking ()), verbose_ (src.verbose ())
  {}

  /// \brief Mapping to an Allocator for a different type \c U.
  ///
  /// Most users do not need this.
  ///
  /// \note To Teuchos developers: This is only optional if the
  ///   Allocator has the form Allocator<T, Args>, where Args is zero
  ///   or more additional template parameters.  I don't want to make
  ///   Allocator require C++11, so it does not have this form.  Thus,
  ///   the rebind struct is required.
  template<class U>
  struct rebind { typedef Allocator<U> other; };

  /// \brief Upper bound (possibly loose) on maximum allocation size.
  ///
  /// Implementations of the Allocator concept should not NEED this
  /// method, but it makes the Clang compiler happy.
  size_type max_size() const {
    return std::numeric_limits<size_type>::max();
  }

  /// \brief Allocate an array of n instances of value_type.
  ///
  /// \param n [in] Number of entries in the array.
  ///
  /// The optional second argument provides an "allocation hint."
  /// This implementation ignores the hint.
  ///
  /// \return The allocated but uninitialized array.
  value_type* allocate (const size_type& n, const void* = 0) {
    if (tracking ()) {
      AllocationLogger::logAllocation (std::cerr, n, n * sizeof (value_type),
                                       typeid (value_type).name (), verbose_);
    }
    return (value_type*) (::operator new (n * sizeof (T)));
  }

  /// \brief Deallocate n instances of value_type.
  ///
  /// \param p [in] Pointer to the array to deallocate.
  /// \param n [in] Number of entries in the array p.
  void deallocate (value_type* p, const size_type& n) {
    if (tracking ()) {
      // Thankfully, this method accepts the array size.  Thus, we don't
      // have to do tricks like allocating extra space and stashing the
      // size in the array.
      AllocationLogger::logDeallocation (std::cerr, n, n * sizeof (value_type),
                                         typeid (value_type).name (), verbose_);
    }
    ::operator delete ((void*) p);
  }

  //! Current total allocation in bytes, over all Allocator<U>.
  size_type curAllocInBytes () {
    return AllocationLogger::curAllocInBytes ();
  }

  //! Max total allocation ("high water mark") in bytes, over all Allocator<U>.
  size_type maxAllocInBytes () {
    return AllocationLogger::maxAllocInBytes ();
  }

#ifndef HAVE_TEUCHOSCORE_CXX11
  /// \brief Invoke the constructor of an instance of \c T, without
  ///   allocating.
  ///
  /// \warning This variant only exists if C++11 is OFF.  C++11
  ///   requires a method with a different signature, but it also
  ///   supplies a good default implementation if this method is
  ///   missing.
  ///
  /// \param p [in] Pointer to an area of memory in which to construct
  ///   a T instance, using placement new.
  /// \param val [in] Argument to T's (copy) constructor.
  void construct (pointer p, const_reference val) {
    new ((void*) p) T (val);
  }
#endif // HAVE_TEUCHOSCORE_CXX11

#ifndef HAVE_TEUCHOSCORE_CXX11
  /// \brief Invoke the destructor of an instance of \c T, without
  ///   deallocating.
  ///
  /// \warning This variant only exists if C++11 is OFF.  C++11
  ///   requires a method with a different signature, but it also
  ///   supplies a good default implementation if this method is
  ///   missing.
  ///
  /// \param p [in] Pointer to an instance of \c T to destroy.
  void destroy (pointer p) {
    ((T*)p)->~T ();
  }
#endif // HAVE_TEUCHOSCORE_CXX11

private:
  bool track_;
  bool verbose_;
};


/// \brief Return true if and only if the two input Allocator
///   instances are interchangeable.
///
/// All instances of our Allocator class use the same memory space,
/// memory pool, and allocation mechanism.  Thus, by this function
/// returning true, we state that all instances of our Allocator class
/// are interchangeable.
template<class T, class U>
bool operator== (const Allocator<T>&, const Allocator<U>&) {
  return true;
}

//! Return <tt> ! (a_t == a_u) </tt> (see above).
template<class T, class U>
bool operator!= (const Allocator<T>& a_t, const Allocator<U>& a_u) {
  return ! (a_t == a_u);
}

} // namespace Details
} // namespace Teuchos

#endif // TEUCHOS_DETAILS_ALLOCATOR
