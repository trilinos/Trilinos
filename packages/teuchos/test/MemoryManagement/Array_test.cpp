// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_Array.new.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_TestingHelpers.hpp"


//
// Define local macros to make defining tests easier for this particular test
// code.
//
// Note, macros with these types of names should only exist in a *.cpp file
// after all #includes are done!
//


#define TEST_EQUALITY_CONST( v1, v2 ) \
  TEUCHOS_TEST_EQUALITY_CONST( v1, v2, out, success )

#define TEST_EQUALITY( v1, v2 ) \
  TEUCHOS_TEST_EQUALITY( v1, v2, out, success )

#define TEST_ITER_EQUALITY( iter1, iter2 ) \
  TEUCHOS_TEST_ITER_EQUALITY( iter1, iter2, out, success )

#define TEST_ARRAY_ELE_EQUALITY( a, i, val ) \
   TEUCHOS_TEST_ARRAY_ELE_EQUALITY( a, i, val, false, out, local_success )

#define TEST_COMPARE( v1, comp, v2 ) \
  TEUCHOS_TEST_COMPARE( v1, comp, v2, out, success )

#define TEST_COMPARE_ARRAYS( a1, a2 ) \
  { \
    const bool result = compareArrays(a1,#a1,a2,#a2,out); \
    if (!result) success = false; \
  }

#define TEST_THROW( code, ExceptType  ) \
  TEUCHOS_TEST_THROW( code, ExceptType, out, success  )

#define TEST_NOTHROW( code  ) \
  TEUCHOS_TEST_NOTHROW( code, out, success  )


//
// Main templated array test function
//


template<class T>
bool testArray( const int n, Teuchos::FancyOStream &out )
{
  
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::setToNull;
  using Teuchos::outArg;
  using Teuchos::getConst;
  using Teuchos::NullIteratorTraits;
  using Teuchos::TypeNameTraits;
  using Teuchos::as;
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
  TEST_COMPARE( a.max_size(), >=, as<size_type>(n) );
  TEST_COMPARE( as<int>(a.capacity()), >=, n );
 
  {
    out << "\nInitializing data ...\n";
    for( int i = 0; i < n; ++i )
      a[i] = i; // tests non-const operator[](i)
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
    TEST_EQUALITY_CONST( as<int>(a2.empty()), true );
  }

  {
    out << "\nTest copy conversion to and from Teuchos::Array and std::vector ...\n";
    std::vector<T> v2 = a.toVector();
    Array<T> a2(v2);
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest assignment operator taking an std::vector ...\n";
    std::vector<T> v2 = a.toVector();
    Array<T> a2;
    a2 = v2;
    TEST_COMPARE_ARRAYS( a2, a );
  }

  {
    out << "\nTest construction using iterators ...\n";
    std::vector<T> v2 = a.toVector();
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

  {
    out << "\nTest that an iterator reference set to null does not throw ...\n";
    typedef typename Array<T>::iterator iter_t;
    iter_t iter = NullIteratorTraits<iter_t>::getNull();
    TEST_NOTHROW( Array<T> a2(n); iter = a2.begin(); setToNull(outArg(iter)) );
  }

  {
    out << "\nTest that a dangling iterator reference throws exception ...\n";
    typedef typename Array<T>::iterator iter_t;
    iter_t iter = NullIteratorTraits<iter_t>::getNull();
    TEST_THROW( { Array<T> a2(n); iter = a2.begin(); }, Teuchos::DanglingReferenceError );
  }

#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  //
  out << "\nE) Test insertion and deletion functions ...\n";
  //

  {
    out << "\nTest push_back(x) ...\n";
    Array<T> a2;
    for ( int i = 0; i < n; ++i ) {
      a2.push_back(i);
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
      const typename Array<T>::iterator iter = a2.insert(a2.end(),i);
      TEST_EQUALITY(*iter,as<T>(i));
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
      a2.append(i);
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
    out << "\nTest swap() ...\n";
    Array<T> a2(a);
    Array<T> a3(a);
    for ( int i = 0; i < n; ++i )
      a2[i] += 1;
    a2.swap(a3);
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
  {
    out << "\nTest to string ...\n";
    std::ostringstream o;
    o << "{";
    for ( int i = 0; i < n; ++i ) {
      o << as<T>(i) << ( i < n-1 ? ", " : "" );
    }
    o << "}";
    TEST_EQUALITY( o.str(), a.toString() );
  }
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
    const ArrayView<const T> cav2 = getConst(a)(0,n);
    TEST_COMPARE_ARRAYS( cav2, a );
  }

  {
    out << "\nTest full non-const to const subview ...\n";
    const ArrayView<const T> cav2 = a(0,n);
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

  // ToDo: Add more tests!

  // ToDo: Test requesting views outside of valid range!
  
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
 
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
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

    result = testArray<std::complex<double> >(n,*out);
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
