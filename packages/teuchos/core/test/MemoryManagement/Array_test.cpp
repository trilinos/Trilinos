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

#include "Teuchos_Array.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_LocalTestingHelpers.hpp"

#include "TestClasses.hpp"


//
// Main templated array test function
//


template<class T>
bool testArray( const int n, Teuchos::FancyOStream &out )
{

  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::outArg;
  using Teuchos::getConst;
  using Teuchos::NullIteratorTraits;
  using Teuchos::TypeNameTraits;
  using Teuchos::as;
  using Teuchos::tuple;
  typedef typename Array<T>::size_type size_type;

  bool success = true;

  out
    << "\n***"
    << "\n*** Testing "<<TypeNameTraits<Array<T> >::name()<<" of size = "<<n
    << "\n***\n";

  Teuchos::OSTab tab(out);

  //
  out << "\nA) Initial setup ...\n\n";
  //

  // Tests construction using size

  Array<T> a(n);

  TEST_EQUALITY_CONST( a.empty(), false );
  TEST_EQUALITY( a.length(), n );
  TEST_EQUALITY( as<int>(a.size()), n );
  TEST_EQUALITY( a.getRawPtr(), &a[0] );
  TEST_EQUALITY( getConst(a).getRawPtr(), &getConst(a)[0] );
  TEST_COMPARE( a.max_size(), >=, as<size_type>(n) );
  TEST_COMPARE( as<int>(a.capacity()), >=, n );

  {
    out << "\nInitializing data ...\n";
    for( int i = 0; i < n; ++i )
      a[i] = as<T>(i); // tests non-const operator[](i)
  }

  {
    out << "\nTest that a[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( a, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    out << "\nTest that a.at(i) == i ...\n";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEUCHOS_TEST_EQUALITY( a.at(i), as<T>(i), out, local_success );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  //
  out << "\nB) Test constructors, assignment operators etc ...\n";
  //

  {
    out << "\nTest default constructor ...\n";
    Array<T> a2;
    TEST_EQUALITY_CONST( as<int>(a2.size()), 0 );
    TEST_EQUALITY_CONST( as<bool>(a2.empty()), true );
    TEST_EQUALITY_CONST( a2.getRawPtr(), 0 );
    TEST_EQUALITY_CONST( getConst(a2).getRawPtr(), 0 );
  }

  {
    out << "\nTest copy conversion to and from Teuchos::Array and std::vector ...\n";
    std::vector<T> v2 = createVector(a);
    Array<T> a2(v2);
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest assignment operator taking an std::vector ...\n";
    std::vector<T> v2 = createVector(a);
    Array<T> a2;
    a2 = v2;
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest construction using iterators ...\n";
    std::vector<T> v2 = createVector(a);
    Array<T> a2(a.begin(),a.end());
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest copy construction ...\n";
    Array<T> a2(a);
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest array assignment operator ...\n";
    Array<T> a2;
    a2 = a;
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest array assign(...) ...\n";
    Array<T> a2;
    a2.assign(a.begin(),a.end());
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest iterator access and then resize ...\n";
    Array<T> a2(a);
    const Array<T> &ca2 = a2;
    Array<T> a3(ca2.begin(),ca2.end());
    TEST_COMPARE_ARRAYS( a3, a );
    TEST_NOTHROW(a2.resize(0)); // This used to throw exception!
  }

  //
  out << "\nC) Test element access ...\n";
  //

  TEST_EQUALITY_CONST( a.front(), as<T>(0) );
  TEST_EQUALITY( a.back(), as<T>(n-1) );
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW( a[-1], Teuchos::RangeError );
  TEST_THROW( a[n], Teuchos::RangeError );
  TEST_THROW( a.at(-1), Teuchos::RangeError );
  TEST_THROW( a.at(n), Teuchos::RangeError );
#else //HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW( a.at(-1), std::out_of_range );
  TEST_THROW( a.at(n), std::out_of_range );
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  //
  out << "\nD) Test iterator access ...\n";
  //

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  {
    out << "\nTesting functions that should throw for empty container ...\n";
    Array<T> a2;
    TEST_THROW( *a2.begin(), Teuchos::NullReferenceError );
    TEST_THROW( a2.front(), Teuchos::NullReferenceError );
    TEST_THROW( a2.back(), Teuchos::NullReferenceError );
    TEST_THROW( getConst(a2).front(), Teuchos::NullReferenceError );
    TEST_THROW( getConst(a2).back(), Teuchos::NullReferenceError );
    TEST_THROW( a2.erase(a2.begin()), Teuchos::NullReferenceError );
  }

#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  {
    out << "\nTest that a2.begin() == a2.end() for empty a2 ...\n";
    Array<T> a2;
    TEST_ITER_EQUALITY( a2.begin(), a2.end() );
  }

  {
    out << "\nTest nonconst forward iterator access ... ";
    bool local_success = true;
    typedef typename Array<T>::iterator iter_t;
    iter_t iter = a.begin();
    for ( int i = 0; i < n; ++i, ++iter )
      TEST_ARRAY_ELE_EQUALITY( a, i, *iter );
    iter = NullIteratorTraits<iter_t>::getNull();
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    out << "\nTest const forward iterator access ... ";
    bool local_success = true;
    typedef typename Array<T>::const_iterator iter_t;
    iter_t iter = getConst(a).begin();
    for ( int i = 0; i < n; ++i, ++iter )
      TEST_ARRAY_ELE_EQUALITY( a, i, *iter );
    iter = NullIteratorTraits<iter_t>::getNull();
    if (local_success) out << "passed\n";
    else success = false;
  }

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  {
    out << "\nTest forward iterators dereferenced out of bounds ...\n";
    TEST_THROW( *(a.begin()-1), Teuchos::RangeError );
    TEST_THROW( *a.end(), Teuchos::RangeError );
    TEST_THROW( *(getConst(a).begin()-1), Teuchos::RangeError );
    TEST_THROW( *getConst(a).end(), Teuchos::RangeError );
  }

#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  {
    out << "\nTest that a2.rbegin() == a2.rend() for empty a2 ...\n";
    Array<T> a2;
    TEST_ITER_EQUALITY( a2.rbegin(), a2.rend() );
  }

  {
    out << "\nTest nonconst reverse iterator access ... ";
    bool local_success = true;
    typedef typename Array<T>::reverse_iterator iter_t;
    iter_t iter = a.rbegin();
    for ( int i = n-1; i >= 0; --i, ++iter )
      TEST_ARRAY_ELE_EQUALITY( a, i, *iter );
    iter = NullIteratorTraits<iter_t>::getNull();
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    out << "\nTest const reverse iterator access ... ";
    bool local_success = true;
    typedef typename Array<T>::const_reverse_iterator iter_t;
    iter_t iter = getConst(a).rbegin();
    for ( int i = n-1; i >= 0; --i, ++iter )
      TEST_ARRAY_ELE_EQUALITY( a, i, *iter );
    iter = NullIteratorTraits<iter_t>::getNull();
    if (local_success) out << "passed\n";
    else success = false;
  }

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  {
    out << "\nTest reverse iterators dereferenced out of bounds ...\n";
    TEST_THROW( *(a.rbegin()-1), Teuchos::RangeError );
    TEST_THROW( *a.rend(), Teuchos::RangeError );
    TEST_THROW( *(getConst(a).rbegin()-1), Teuchos::RangeError );
    TEST_THROW( *getConst(a).rend(), Teuchos::RangeError );
  }
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  {
    out << "\nTest that an iterator reference set to null does not throw ...\n";
    typedef typename Array<T>::iterator iter_t;
    iter_t iter = NullIteratorTraits<iter_t>::getNull();
    TEST_NOTHROW( Array<T> a2(n); iter = a2.begin();
      iter = NullIteratorTraits<iter_t>::getNull() );
  }

  {
    out << "\nTest that a dangling iterator reference throws exception ...\n";
    typedef typename Array<T>::iterator iter_t;
    iter_t iter = NullIteratorTraits<iter_t>::getNull();
    {
      Array<T> a2(n);
      iter = a2.begin();
    }
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    TEST_THROW(*iter=0, Teuchos::DanglingReferenceError );
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  }


  //
  out << "\nE) Test insertion and deletion functions ...\n";
  //

  {
    out << "\nTest push_back(x) ...\n";
    Array<T> a2;
    for ( int i = 0; i < n; ++i ) {
      a2.push_back(as<T>(i));
      TEST_EQUALITY_CONST(a2.front(),as<T>(0));
      TEST_EQUALITY_CONST(getConst(a2).front(),as<T>(0));
      TEST_EQUALITY(a2.back(),as<T>(i));
      TEST_EQUALITY(getConst(a2).back(),as<T>(i));
    }
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest pop_back() ...\n";
    Array<T> a2(a);
    for ( int i = n-1; i >= 0; --i ) {
      TEST_EQUALITY(a2.back(),as<T>(i));
      a2.pop_back();
    }
  }

  {
    out << "\nTest insert(iter,x) ...\n";
    Array<T> a2;
    for ( int i = 0; i < n; ++i ) {
      const typename Array<T>::iterator
        iter = a2.insert(a2.end(), as<T>(i));
      TEST_EQUALITY(*iter, as<T>(i));
    }
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest insert(iter,1,x) ...\n";
    Array<T> a2;
    for ( int i = 0; i < n; ++i )
      a2.insert(a2.end(),1,i);
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest insert(iter,first,last) ...\n";
    Array<T> a2;
    for ( int i = 0; i < n; ++i )
      a2.insert(a2.end(),a.begin()+i,a.begin()+i+1);
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest append(x) ...\n";
    Array<T> a2;
    for ( int i = 0; i < n; ++i )
      a2.append(as<T>(i));
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest erase(iter) ...\n";
    Array<T> a2(a);
    for ( int i = 0; i < n; ++i ) {
      TEST_EQUALITY( as<int>(a2.size()), n-i );
      TEST_EQUALITY( a2.front(), as<T>(i) );
      a2.erase(a2.begin());
    }
    TEST_EQUALITY_CONST( a2.empty(), true );
  }

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  {
    out << "\nTest trying to erase twice with the same iterator which should throw ...\n";
    Array<T> a2(a);
    const typename Array<T>::iterator iter = a2.begin();
    a2.erase(iter); // After this point, the iterator is no longer valid!
    // This is no longer a valid iterator and should throw!
    TEST_THROW( a2.erase(iter), Teuchos::IncompatibleIteratorsError );
  }

#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  // 2007/11/08: rabartl: ToDo: Above, I have tested one use case where the
  // iterator should be invalidated and this tests that it throws an exception
  // as it should.  However, currently, I don't have code written that will
  // catch the problem where the client would try to dereference the iterator
  // or something like that.  This is a big no-no.  I could do this by adding
  // an is_valid() property to RCP_node and then setting this to null when the
  // structure of the iterator changes.  Then, I would have to put asserts in
  // ArrayRCP to constantly check is_valid() (with an assert_is_valid()
  // function or something) on any call other than operator=(...) which would
  // reset this iterator.  Catching all of these user errors is a lot of work!

  {
    out << "\nTest remove(i) ...\n";
    Array<T> a2(a);
    for ( int i = 0; i < n; ++i ) {
      TEST_EQUALITY( as<int>(a2.size()), n-i );
      TEST_EQUALITY( a2.front(), as<T>(i) );
      a2.remove(0); // Always remove the "first" entry!
    }
    TEST_EQUALITY_CONST( a2.empty(), true );
  }

  {
    out << "\nTest erase(begin(),end()) ...\n";
    Array<T> a2(a);
    a2.erase(a2.begin(),a2.end());
    TEST_EQUALITY_CONST( a2.empty(), true );
  }

  {
    out << "\nTest member swap() ...\n";
    Array<T> a2(a);
    Array<T> a3(a);
    for ( int i = 0; i < n; ++i )
      a2[i] += as<T>(1);
    a2.swap(a3);
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest non-member swap() ...\n";
    Array<T> a2(a);
    Array<T> a3(a);
    for ( int i = 0; i < n; ++i )
      a2[i] += as<T>(1);
    swap(a2,a3);
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest clear() ...\n";
    Array<T> a2(a);
    a2.clear();
    TEST_EQUALITY_CONST( a2.empty(), true );
    TEST_EQUALITY_CONST( a2.size(), 0 );
  }

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  // mfh 28 Aug 2012: This test no longer passes, because we've
  // specialized ArrayView<T>::toString() for T = (const) float,
  // (const) double.  We've done the specialization to print float and
  // double in scientific notation.  That was a hack to fix a bug; it
  // would make more sense to provide a standard toString() for float
  // and double, and have the test (or even the
  // ArrayView<T>::toString() specialization) use that.
  //
  // {
  //   out << "\nTest to string ...\n";
  //   std::ostringstream o;
  //   o << "{";
  //   for ( int i = 0; i < n; ++i ) {
  //     o << as<T>(i) << ( i < n-1 ? ", " : "" );
  //   }
  //   o << "}";
  //   TEST_EQUALITY( o.str(), a.toString() );
  // }
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  {
    out << "\nTest hasArrayBoundsChecking() ... \n";
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    TEST_EQUALITY_CONST( a.hasBoundsChecking(), true );
#else
    TEST_EQUALITY_CONST( a.hasBoundsChecking(), false );
#endif
  }

  //
  out << "\nG) Test views ...\n";
  //

  {
    out << "\nTest full non-const subview ...\n";
    const ArrayView<T> av2 = a(0,n);
    TEST_COMPARE_ARRAYS( av2, a );
  }

  {
    out << "\nTest full shorthand non-const subview ...\n";
    const ArrayView<T> av2 = a();
    TEST_COMPARE_ARRAYS( av2, a );
  }

  {
    out << "\nTest full const subview ...\n";
    const ArrayView<const T> cav2 = getConst(a)(0, n);
    TEST_COMPARE_ARRAYS( cav2, a );
  }

  {
    out << "\nTest full non-const to const subview ...\n";
    const ArrayView<const T> cav2 = a(0, n);
    TEST_COMPARE_ARRAYS( cav2, a );
  }

  {
    out << "\nTest full short-hand const subview ...\n";
    const ArrayView<const T> cav2 = getConst(a)();
    TEST_COMPARE_ARRAYS( cav2, a );
  }

  {
    out << "\nTest non-const initial range view ...\n";
    Array<T> a2(n,as<T>(-1));
    const ArrayView<T> av2 = a2; // Tests implicit conversion!
    const ArrayView<T> av2_end = av2(0,n-1);
    TEST_EQUALITY( av2_end.size(), n-1 );
    av2_end.assign( a(0,n-1) );
    av2.back() = as<T>(n-1);
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest non-const middle range view ...\n";
    Array<T> a2(n,as<T>(-1));
    const ArrayView<T> av2 = a2; // Tests implicit conversion!
    const ArrayView<T> av2_middle = av2(1,n-2);
    TEST_EQUALITY( av2_middle.size(), n-2 );
    av2_middle.assign( a(1,n-2) );
    av2.front() = as<T>(0);
    av2.back() = as<T>(n-1);
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest const view ... ";
    const ArrayView<const T> av2 = a; // Tests implicit conversion to const!
    const ArrayView<const T> av2_middle = av2(1,n-2);
    TEST_EQUALITY( av2_middle.size(), n-2 );
    bool local_success = true;
    for ( int i = 0; i < n-2; ++i )
      TEST_ARRAY_ELE_EQUALITY( av2_middle, i, as<T>(i+1) );
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    out << "\nTest constructing Array<T> from ArrayView<T> ...\n";
    const ArrayView<T> av2 = a;
    Array<T> a2(av2);
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest constructing Array<T> from ArrayView<const T> ...\n";
    const ArrayView<const T> av2 = a;
    Array<T> a2(av2);
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest comparison operators ...\n";
    Array<T> a2(a);
    TEST_EQUALITY_CONST( (a2==a), true );
    TEST_EQUALITY_CONST( (a2!=a), false );
    TEST_EQUALITY_CONST( (a2<=a), true );
    TEST_EQUALITY_CONST( (a2>=a), true );
    TEST_EQUALITY_CONST( (a2<a), false );
    TEST_EQUALITY_CONST( (a2>a), false );
  }

  //
  out << "\nH) Test tuple(...) construction ...\n";
  //

  {
    const size_type m = 1;
    out << "\nTest Array<T> = tuple(0,...,"<<m-1<<")\n";
    Array<T> am = tuple<T>(0);
    TEST_EQUALITY_CONST(am.size(), m);
    out << "Test that am[i] == i ... ";
    bool local_success = true;
    for( size_type i = 0; i < m; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( am, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const size_type m = 2;
    out << "\nTest Array<T> = tuple(0,...,"<<m-1<<")\n";
    Array<T> am = tuple<T>(0,1);
    TEST_EQUALITY_CONST(am.size(),m);
    out << "Test that am[i] == i ... ";
    bool local_success = true;
    for( size_type i = 0; i < m; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( am, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const size_type m = 3;
    out << "\nTest Array<T> = tuple(0,...,"<<m-1<<")\n";
    Array<T> am = tuple<T>(0,1,2);
    TEST_EQUALITY_CONST(am.size(),m);
    out << "Test that am[i] == i ... ";
    bool local_success = true;
    for( size_type i = 0; i < m; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( am, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const size_type m = 4;
    out << "\nTest Array<T> = tuple(0,...,"<<m-1<<")\n";
    Array<T> am = tuple<T>(0,1,2,3);
    TEST_EQUALITY_CONST(am.size(),m);
    out << "Test that am[i] == i ... ";
    bool local_success = true;
    for( size_type i = 0; i < m; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( am, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const size_type m = 5;
    out << "\nTest Array<T> = tuple(0,...,"<<m-1<<")\n";
    Array<T> am = tuple<T>(0,1,2,3,4);
    TEST_EQUALITY_CONST(am.size(),m);
    out << "Test that am[i] == i ... ";
    bool local_success = true;
    for( size_type i = 0; i < m; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( am, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const size_type m = 6;
    out << "\nTest Array<T> = tuple(0,...,"<<m-1<<")\n";
    Array<T> am = tuple<T>(0,1,2,3,4,5);
    TEST_EQUALITY_CONST(am.size(),m);
    out << "Test that am[i] == i ... ";
    bool local_success = true;
    for( size_type i = 0; i < m; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( am, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const size_type m = 7;
    out << "\nTest Array<T> = tuple(0,...,"<<m-1<<")\n";
    Array<T> am = tuple<T>(0,1,2,3,4,5,6);
    TEST_EQUALITY_CONST(am.size(),m);
    out << "Test that am[i] == i ... ";
    bool local_success = true;
    for( size_type i = 0; i < m; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( am, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const size_type m = 8;
    out << "\nTest Array<T> = tuple(0,...,"<<m-1<<")\n";
    Array<T> am = tuple<T>(0,1,2,3,4,5,6,7);
    TEST_EQUALITY_CONST(am.size(),m);
    out << "Test that am[i] == i ... ";
    bool local_success = true;
    for( size_type i = 0; i < m; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( am, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const size_type m = 9;
    out << "\nTest Array<T> = tuple(0,...,"<<m-1<<")\n";
    Array<T> am = tuple<T>(0,1,2,3,4,5,6,7,8);
    TEST_EQUALITY_CONST(am.size(),m);
    out << "Test that am[i] == i ... ";
    bool local_success = true;
    for( size_type i = 0; i < m; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( am, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const size_type m = 10;
    out << "\nTest Array<T> = tuple(0,...,"<<m-1<<")\n";
    Array<T> am = tuple<T>(0,1,2,3,4,5,6,7,8,9);
    TEST_EQUALITY_CONST(am.size(),m);
    out << "Test that am[i] == i ... ";
    bool local_success = true;
    for( size_type i = 0; i < m; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( am, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    out << "\nTest taking an empty view ...\n";
    const ArrayView<T> av = a(0,0);
    TEST_EQUALITY_CONST( av.size(), 0 );
  }

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  {
    out << "\nTest taking views outside of valid range ...\n";
    TEST_THROW( const ArrayView<T> av = a(-1,n), Teuchos::RangeError );
    TEST_THROW( const ArrayView<T> av = a(0,n+1), Teuchos::RangeError );
    TEST_THROW( const ArrayView<T> av = a(0,-1), Teuchos::RangeError );
  }
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  return success;

}


template<class T>
bool testArrayOpaqueWithoutTNT( const std::string &T_name, const int n,
  const T &someValue, Teuchos::FancyOStream &out )
{

  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::TypeNameTraits;
  using Teuchos::as;
  typedef typename Array<T>::size_type size_type;

  bool success = true;

  out
    << "\n***"
    << "\n*** Testing Array<"<<T_name<<"> for opaque type without TNT of size = "<<n
    << "\n***\n";

  Teuchos::OSTab tab(out);

  //
  out << "\nA) Initial setup ...\n\n";
  //

  // Tests construction using size

  Array<T> a(n);

  TEST_EQUALITY_CONST( a.empty(), false );
  TEST_EQUALITY( a.length(), n );
  TEST_EQUALITY( as<int>(a.size()), n );
  TEST_EQUALITY( a.getRawPtr(), &a[0] );
  TEST_EQUALITY( getConst(a).getRawPtr(), &getConst(a)[0] );
  TEST_COMPARE( a.max_size(), >=, as<size_type>(n) );
  TEST_COMPARE( as<int>(a.capacity()), >=, n );

  {
    out << "\nInitializing data ...\n";
    for( int i = 0; i < n; ++i )
      a[i] = someValue; // tests non-const operator[](i)
  }

  {
    out << "\nTest that a[i] == "<<someValue<<" ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( a, i, someValue );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

#ifndef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  {
    out << "\nTest taking a view of the array ...\n";
    const ArrayView<T> av = a();
    TEST_COMPARE_ARRAYS( av, a );
  }
  // 2008/08/01: rabartl: Above: We can not create an array view of an
  // undefined type in debug mode without a specialization of TypeNameTraits.
#endif // not HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  // ToDo: Do we need to be testing other things for opaque objects?

  return success;

}


template<class T>
bool testArrayOpaqueWithTNT( const int n, const T &someValue, Teuchos::FancyOStream &out )
{

  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::TypeNameTraits;
  using Teuchos::as;
  typedef typename Array<T>::size_type size_type;

  bool success = true;

  out
    << "\n***"
    << "\n*** Testing "<<TypeNameTraits<Array<T> >::name()<<" for opaque type with TNT of size = "<<n
    << "\n***\n";

  Teuchos::OSTab tab(out);

  //
  out << "\nA) Initial setup ...\n\n";
  //

  // Tests construction using size

  Array<T> a(n);

  TEST_EQUALITY_CONST( a.empty(), false );
  TEST_EQUALITY( a.length(), n );
  TEST_EQUALITY( as<int>(a.size()), n );
  TEST_EQUALITY( a.getRawPtr(), &a[0] );
  TEST_EQUALITY( getConst(a).getRawPtr(), &getConst(a)[0] );
  TEST_COMPARE( a.max_size(), >=, as<size_type>(n) );
  TEST_COMPARE( as<int>(a.capacity()), >=, n );

  {
    out << "\nInitializing data ...\n";
    for( int i = 0; i < n; ++i )
      a[i] = someValue; // tests non-const operator[](i)
  }

  {
    out << "\nTest that a[i] == "<<someValue<<" ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( a, i, someValue );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    out << "\nTest taking a view of the array ...\n";
    const ArrayView<T> av = a();
    TEST_COMPARE_ARRAYS( av, a );
  }

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  {
    out << "\nTest taking views outside of valid range ...\n";
    TEST_THROW( const ArrayView<T> av = a(-1,n), Teuchos::RangeError );
    TEST_THROW( const ArrayView<T> av = a(0,n+1), Teuchos::RangeError );
    TEST_THROW( const ArrayView<T> av = a(0,-1), Teuchos::RangeError );
  }
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  // 2008/08/01: rabartl: Above: We can create ArrayViews and any other thing
  // that we would like since we have defined a TypeNameTraits class for the
  // undefined type.

  // ToDo: Do we need to be testing other things for opaque objects?

  return success;

}


//
// Main testing program
//

int main( int argc, char* argv[] ) {

  using Teuchos::CommandLineProcessor;
  using Teuchos::Array;

        bool success = true;
  bool result;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  //const int procRank = Teuchos::GlobalMPISession::getRank();

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

        try {

    //
                // Read options from the commandline
    //

    CommandLineProcessor clp(false); // Don't throw exceptions

    int n = 4;
    clp.setOption( "n", &n, "Number of elements in the array" );

                CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);

                if ( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
                        *out << "\nEnd Result: TEST FAILED" << std::endl;
                        return parse_return;
                }

    *out << std::endl << Teuchos::Teuchos_Version() << std::endl;

    result = testArray<int>(n,*out);
    if (!result) success = false;

    result = testArray<float>(n,*out);
    if (!result) success = false;

    result = testArray<double>(n,*out);
    if (!result) success = false;

    //result = testArray<std::complex<double> >(n,*out);
    //if (!result) success = false;
    // 2007/12/03: rabartl: Commented this out so I can test comparison operators

    result = testArrayOpaqueWithoutTNT<Opaque_handle>("Opaque_handle", n,
      OPAQUE_HANDLE_NULL, *out);
    if (!result) success = false;

    result = testArrayOpaqueWithTNT<Opaque2_handle>(n, OPAQUE2_HANDLE_NULL, *out);
    if (!result) success = false;

    result = testArrayOpaqueWithTNT<Opaque3_handle>(n, OPAQUE3_HANDLE_NULL, *out);
    if (!result) success = false;

    // ToDo: Fill in the rest of the types!

        }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);

  if (success)
    *out << "\nEnd Result: TEST PASSED" << std::endl;
  else
    *out << "\nEnd Result: TEST FAILED" << std::endl;

  return ( success ? 0 : 1 );

}
