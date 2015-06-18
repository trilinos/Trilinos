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

#ifndef TEUCHOS_DETAILS_ALLOCATOR
#define TEUCHOS_DETAILS_ALLOCATOR

#include <memory>
#include <iostream>

namespace Teuchos {
namespace Details {

/// \brief Logging implementation used by Allocator (see below).
///
/// \warning This is an implementation detail of Allocator, which in
///   turn is an implementation detail of Teuchos.  DO NOT USE THIS
///   CLASS DIRECTLY IN YOUR CODE.
class AllocationLogger {
public:
  /// \brief Type of the size of an allocation or deallocation.
  ///
  /// This must match Allocation::size_type (see below).
  typedef std::size_t size_type;

  static void
  logAllocation (std::ostream& out,
                 const size_type numEntries,
                 const size_type numBytes,
                 const char typeName[],
                 const bool verbose);

  static void
  logDeallocation (std::ostream& out,
                   const size_type numEntries,
                   const size_type numBytes,
                   const char typeName[],
                   const bool verbose);

  //! Current total allocation in bytes.
  static size_type curAllocInBytes ();

  //! Max total allocation ("high water mark") in bytes.
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

/// \brief Tracking allocator for Teuchos Memory Management classes.
/// \tparam T The type of the entries in arrays.
template<class T>
class Allocator {
private:
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
  //! Type of a pointer to T.
  typedef T* pointer;
  //! Type of a pointer to const T.
  typedef const T* const_pointer;
  //! Type of a reference to T.
  typedef T& reference;
  //! Type of a reference to const T.
  typedef const T& const_reference;

  /// \brief Type of the size of an allocation or deallocation.
  ///
  /// This must match AllocationLogger::size_type (see above).
  typedef std::size_t size_type;

  //! Signed integer type representing the difference of two pointers.
  typedef std::ptrdiff_t difference_type;

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
      AllocationLogger::logAllocation (std::cerr, n, n * sizeof (value_type), typeid (value_type).name (), verbose_);
    }
    return (pointer) (::operator new (n * sizeof (T)));
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
      AllocationLogger::logDeallocation (std::cerr, n, n * sizeof (value_type), typeid (value_type).name (), verbose_);
    }
    ::operator delete ((void*) p);
  }

  /// \brief Invoke the constructor of an instance of U, without allocating.
  /// \tparam U The type of the object to construct.
  ///
  /// The "class..." thing and std::forward are C++11 syntax for a
  /// variable number of arguments corresponding to possibly different
  /// template parameters.  That lets this method deal with any kind
  /// of constructor that U might have.
  template<class U, class... Args>
  void construct (U* p, Args&&... args) {
    new ((void*) p) U (std::forward<Args> (args)...);
  }

  /// \brief Invoke the destructor of an instance of U, without deallocating.
  /// \tparam U The type of the object to destroy.
  template<class U>
  void destroy (U* p) {
    p->~U ();
  }

  template<class U>
  struct rebind { typedef Allocator<U> other; };

  //! Current total allocation in bytes, over all Allocator<U>.
  size_type curAllocInBytes () {
    return AllocationLogger::curAllocInBytes ();
  }

  //! Max total allocation ("high water mark") in bytes, over all Allocator<U>.
  size_type maxAllocInBytes () {
    return AllocationLogger::maxAllocInBytes ();
  }

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
