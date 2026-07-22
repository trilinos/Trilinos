// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_Details_Allocator.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include "Teuchos_UnitTestHarness.hpp"
#include <string>
#include <vector>

namespace { // (anonymous)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Allocator, Test1, T )
{
  using Teuchos::Details::AllocationLogger;
  using Teuchos::TypeNameTraits;
  using std::endl;

  Teuchos::OSTab tab0 (out);
  out << "Test Teuchos::Details::Allocator for T = "
      << TypeNameTraits<T>::name () << endl;
  AllocationLogger::resetAllocationCounts ();

  typedef Teuchos::Details::Allocator<T> alloc_type;
  alloc_type alloc;

  typedef typename alloc_type::size_type size_type;

  // At this point, we haven't allocated anything yet.  The allocator
  // does not track whatever memory it uses in its implementation.
  TEST_EQUALITY_CONST( alloc.curAllocInBytes (), static_cast<size_type> (0) );

  // We'll use this below.
  size_type oldMaxAlloc = 0;

  // Put the test in an inner scope, so that the std::vector gets
  // deallocated before this test finishes.  This lets us print
  // whether deallocation succeeded.
  {
    const size_type numEntries = 10;

    typedef std::vector<T, alloc_type> vec_type;
    // C++14 defines a two-argument std::vector constructor (count,
    // alloc), but C++11 only has a three-argument constructor that
    // takes a count and the allocator.  Thus, we need some default
    // value T.  I choose 22 because it is one plus the sum of
    // integers from 1 to 10, inclusive.  It's not a default value,
    // like zero, and it's positive, so it works if T is unsigned.
    // It also fits exactly in float and double.
    T val = static_cast<T> (22);
    vec_type vec (numEntries, val, alloc);

    TEST_EQUALITY( vec.size (), numEntries );
    TEST_EQUALITY_CONST( vec.capacity () >= vec.size (), true );

    oldMaxAlloc = alloc.maxAllocInBytes ();
    const size_type curAlloc = alloc.curAllocInBytes ();
    const size_type expectedCurAlloc = numEntries * sizeof (T);

    // We don't need strict equality, because the allocator may choose
    // to allocate more memory than necessary (e.g., to stash
    // additional information in each allocation).
    TEST_EQUALITY_CONST( curAlloc >= expectedCurAlloc, true );
    TEST_EQUALITY_CONST( oldMaxAlloc >= expectedCurAlloc, true );

    // Make sure that the std::vector's constructor correctly filled
    // it using val.  We have to test this because std::vector defers
    // to alloc_type::construct for this.
    for (size_type k = 0; k < numEntries; ++k) {
      TEST_EQUALITY( vec[k], val );
    }
  }

  // At this point, alloc.curAlloc() should be zero, and
  // alloc.maxAlloc() should not have changed.
  const size_type newMaxAlloc = alloc.maxAllocInBytes ();
  TEST_EQUALITY( oldMaxAlloc, newMaxAlloc );
  TEST_EQUALITY_CONST( alloc.curAllocInBytes (), static_cast<size_type> (0) );

  out << "Done with test!" << endl;
}

//
// Repeat Test1, but with verbose logging on.
//
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Allocator, Test2, T )
{
  using Teuchos::Details::AllocationLogger;
  using Teuchos::TypeNameTraits;
  using std::endl;

  Teuchos::OSTab tab0 (out);
  out << "Test Teuchos::Details::Allocator for T = "
      << TypeNameTraits<T>::name () << ", with verbose logging on" << endl;
  AllocationLogger::resetAllocationCounts ();

  typedef Teuchos::Details::Allocator<T> alloc_type;
  // Tell the Allocator to track memory.
  alloc_type alloc (true, true);

  typedef typename alloc_type::size_type size_type;

  // At this point, we haven't allocated anything yet.  The allocator
  // does not track whatever memory it uses in its implementation.
  TEST_EQUALITY_CONST( alloc.curAllocInBytes (), static_cast<size_type> (0) );

  // We'll use this below.
  size_type oldMaxAlloc = 0;

  // Put the test in an inner scope, so that the std::vector gets
  // deallocated before this test finishes.  This lets us print
  // whether deallocation succeeded.
  {
    const size_type numEntries = 10;

    typedef std::vector<T, alloc_type> vec_type;
    // C++14 defines a two-argument std::vector constructor (count,
    // alloc), but C++11 only has a three-argument constructor that
    // takes a count and the allocator.  Thus, we need some default
    // value T.  I choose 22 because it is one plus the sum of
    // integers from 1 to 10, inclusive.  It's not a default value,
    // like zero, and it's positive, so it works if T is unsigned.
    // It also fits exactly in float and double.
    T val = static_cast<T> (22);
    vec_type vec (numEntries, val, alloc);

    TEST_EQUALITY( vec.size (), numEntries );
    TEST_EQUALITY_CONST( vec.capacity () >= vec.size (), true );

    oldMaxAlloc = alloc.maxAllocInBytes ();
    const size_type curAlloc = alloc.curAllocInBytes ();
    const size_type expectedCurAlloc = numEntries * sizeof (T);

    // We don't need strict equality, because the allocator may choose
    // to allocate more memory than necessary (e.g., to stash
    // additional information in each allocation).
    TEST_EQUALITY_CONST( curAlloc >= expectedCurAlloc, true );
    TEST_EQUALITY_CONST( oldMaxAlloc >= expectedCurAlloc, true );

    // Make sure that the std::vector's constructor correctly filled
    // it using val.  We have to test this because std::vector defers
    // to alloc_type::construct for this.
    for (size_type k = 0; k < numEntries; ++k) {
      TEST_EQUALITY( vec[k], val );
    }
  }

  // At this point, alloc.curAlloc() should be zero, and
  // alloc.maxAlloc() should not have changed.
  const size_type newMaxAlloc = alloc.maxAllocInBytes ();
  TEST_EQUALITY( oldMaxAlloc, newMaxAlloc );
  TEST_EQUALITY_CONST( alloc.curAllocInBytes (), static_cast<size_type> (0) );

  out << "Done with test!" << endl;
}

//
// Make sure that mixing std::vector<T> instances for different T
// still gives the right current and max allocation numbers.
//
TEUCHOS_UNIT_TEST( Allocator, TestMixedTypes )
{
  using Teuchos::Details::AllocationLogger;
  using std::endl;

  Teuchos::OSTab tab0 (out);
  out << "Test Teuchos::Details::Allocator<T> for mixed T" << endl;
  AllocationLogger::resetAllocationCounts ();

  typedef Teuchos::Details::Allocator<int> int_alloc_type;
  typedef int_alloc_type::size_type size_type;
  const size_type numEntries = 10;
  const size_type expectedMaxAlloc = numEntries * sizeof (int) + numEntries * sizeof (double);

  // At this point, we haven't allocated anything yet.
  TEST_EQUALITY_CONST( AllocationLogger::curAllocInBytes (), static_cast<size_type> (0) );

  {
    std::vector<int, int_alloc_type> intVec (numEntries);

    typedef Teuchos::Details::Allocator<double> double_alloc_type;
    std::vector<double, double_alloc_type> dblVec (numEntries);

    // Both std::vector types must report the same current and max total
    // allocation sizes, since they share the same allocation mechanism.
    TEST_EQUALITY( intVec.get_allocator ().curAllocInBytes (), dblVec.get_allocator ().curAllocInBytes () );
    TEST_EQUALITY( intVec.get_allocator ().maxAllocInBytes (), dblVec.get_allocator ().maxAllocInBytes () );


    TEST_EQUALITY_CONST( intVec.get_allocator ().curAllocInBytes () >= expectedMaxAlloc, true );
    TEST_EQUALITY_CONST( intVec.get_allocator ().maxAllocInBytes () >= expectedMaxAlloc, true );
  }

  TEST_EQUALITY_CONST( AllocationLogger::curAllocInBytes (), static_cast<size_type> (0) );
  TEST_EQUALITY_CONST( AllocationLogger::maxAllocInBytes () >= expectedMaxAlloc, true );

  out << "Done with test!" << endl;
}

//
// Make sure that the Allocator works for types T that do run-time
// allocation.  std::string is a good example.
//
// This is the test that shows why you CANNOT use "return new T[n]" to
// implement allocate(), and "delete [] p" to implement deallocate().
// (Try it on a Mac and watch your debug malloc complain that you're
// trying to free something it never malloc'd.)
//
TEUCHOS_UNIT_TEST( Allocator, TestString )
{
  using Teuchos::Details::AllocationLogger;
  using std::endl;

  Teuchos::OSTab tab0 (out);
  out << "Test Teuchos::Details::Allocator<std::string>" << endl;
  AllocationLogger::resetAllocationCounts ();

  typedef Teuchos::Details::Allocator<std::string> string_alloc_type;
  typedef string_alloc_type::size_type size_type;
  const size_type numEntries = 10;
  // Even though std::string itself does run-time allocation inside,
  // this is still the correct max allocation for an array of
  // std::string.
  const size_type expectedMaxAlloc = numEntries * sizeof (std::string);

  // At this point, we haven't allocated anything yet.
  TEST_EQUALITY_CONST( AllocationLogger::curAllocInBytes (), static_cast<size_type> (0) );

  // First, try it without setting any of the strings.
  {
    std::vector<std::string, string_alloc_type> vec (numEntries);

    TEST_EQUALITY_CONST( vec.get_allocator ().curAllocInBytes () >= expectedMaxAlloc, true );
    TEST_EQUALITY_CONST( vec.get_allocator ().maxAllocInBytes () >= expectedMaxAlloc, true );
  }

  TEST_EQUALITY_CONST( AllocationLogger::curAllocInBytes (), static_cast<size_type> (0) );
  TEST_EQUALITY_CONST( AllocationLogger::maxAllocInBytes () >= expectedMaxAlloc, true );

  // Next, construct the std::vector, setting all entries to a string
  // of nonzero length.
  {
    string_alloc_type alloc;
    std::string val ("I'm a little teapot, short and stout.");
    std::vector<std::string, string_alloc_type> vec (numEntries, val, alloc);

    TEST_EQUALITY_CONST( vec.get_allocator ().curAllocInBytes () >= expectedMaxAlloc, true );
    TEST_EQUALITY_CONST( vec.get_allocator ().maxAllocInBytes () >= expectedMaxAlloc, true );
  }

  TEST_EQUALITY_CONST( AllocationLogger::curAllocInBytes (), static_cast<size_type> (0) );
  TEST_EQUALITY_CONST( AllocationLogger::maxAllocInBytes () >= expectedMaxAlloc, true );

  // Next, construct the std::vector without setting any of the
  // strings, then iterate through it and set the strings one by one
  // to different values (that circumvents possible reference counting
  // in std::string).
  {
    string_alloc_type alloc;
    std::vector<std::string, string_alloc_type> vec (numEntries);

    for (size_type k = 0; k < numEntries; ++k) {
      std::ostringstream os;
      os << "Current index: " << k;
      vec[k] = os.str ();
    }

    TEST_EQUALITY_CONST( vec.get_allocator ().curAllocInBytes () >= expectedMaxAlloc, true );
    TEST_EQUALITY_CONST( vec.get_allocator ().maxAllocInBytes () >= expectedMaxAlloc, true );
  }

  TEST_EQUALITY_CONST( AllocationLogger::curAllocInBytes (), static_cast<size_type> (0) );
  TEST_EQUALITY_CONST( AllocationLogger::maxAllocInBytes () >= expectedMaxAlloc, true );

  out << "Done with test!" << endl;
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Allocator, Test1, int )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Allocator, Test1, double )


} // namespace (anonymous)
