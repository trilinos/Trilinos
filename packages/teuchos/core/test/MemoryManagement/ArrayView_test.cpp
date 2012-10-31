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

#include "Teuchos_ArrayView.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Version.hpp"
#include "Teuchos_getConst.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_TestingHelpers.hpp"


// Uncomment to show compile errors from invalid usage
//#define SHOW_INVALID_COPY_CONSTRUCTION
//#define SHOW_INVALID_CONST_ASSIGN
//#define SHOW_INVALID_CONST_ITER_MODIFICATION

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
bool testArrayView( const int n, Teuchos::FancyOStream &out )
{
  
  using Teuchos::ArrayView;
  using Teuchos::arrayView;
  using Teuchos::arrayViewFromVector;
  using Teuchos::outArg;
  using Teuchos::NullIteratorTraits;
  using Teuchos::TypeNameTraits;
  using Teuchos::getConst;
  using Teuchos::as;
  typedef typename ArrayView<T>::size_type size_type;

  bool success = true;
 
  out
    << "\n***"
    << "\n*** Testing "<<TypeNameTraits<ArrayView<T> >::name()<<" of size = "<<n
    << "\n***\n";
  
  Teuchos::OSTab tab(out);

  //
  out << "\nA) Initial setup testing ...\n\n";
  //
  
  {
    out << "\nTesting basic null construction!\n\n";
    ArrayView<T> av2 = Teuchos::null;
    TEST_EQUALITY_CONST(is_null(av2),true);
    TEST_EQUALITY_CONST(av2.size(),0);
    TEST_EQUALITY_CONST(av2.getRawPtr(),0);
    TEST_ITER_EQUALITY(av2.begin(),av2.end());
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    TEST_THROW(av2[0],Teuchos::NullReferenceError);
    TEST_THROW(*av2.begin(), Teuchos::NullReferenceError);
    TEST_THROW(*av2.end(), Teuchos::NullReferenceError);
    TEST_THROW(av2.assign(av2), Teuchos::NullReferenceError);
    TEST_THROW(av2.front(), Teuchos::NullReferenceError);
    TEST_THROW(av2.back(), Teuchos::NullReferenceError);
#endif
    ArrayView<const T> cav2(av2); // Tests copy constructor and implicit conversion operator! 
    TEST_EQUALITY_CONST(cav2.size(),0);
    TEST_EQUALITY_CONST(cav2.getRawPtr(),0);
    TEST_ITER_EQUALITY(cav2.begin(),av2.end());
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    TEST_THROW(cav2[0],Teuchos::NullReferenceError);
    TEST_THROW(*cav2.begin(), Teuchos::NullReferenceError);
    TEST_THROW(*cav2.end(), Teuchos::NullReferenceError);
    TEST_THROW(cav2.back(), Teuchos::NullReferenceError);
#endif
#ifdef SHOW_INVALID_CONST_ASSIGN
    TEST_NOTHROW(cav2.assign(av2)); // Should not compile!
#endif
  }

  std::vector<T> v(n);

  const ArrayView<T> av = arrayViewFromVector(v);
  TEST_EQUALITY_CONST(is_null(av), false);
  TEST_EQUALITY( as<int>(av.size()), n );

  const ArrayView<const T> cav = av;
 
  {
    out << "\nInitializing data for std::vector v through view av ...\n";
    for( int i = 0; i < n; ++i )
      av[i] = i; // tests non-const operator[](i)
  }

  {
    out << "\nTest that v[i] == i through ArrayView<const T> ...\n";
    const ArrayView<const T> cav2 = cav;
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( cav2, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    out << "\nTest explicit copy to std::vector from non-const array view ...\n";
    std::vector<T> v2 = Teuchos::createVector(av);
    TEST_COMPARE_ARRAYS( v2, v );
  }

  {
    out << "\nTest explicit copy to std::vector from const array view ...\n";
    std::vector<T> v2 = Teuchos::createVector(cav);
    TEST_COMPARE_ARRAYS( v2, v );
  }

  {
    out << "\nTest shallow implicit conversion from ArrayView<T> to ArrayView<T> ...\n";
    ArrayView<T> av2(av);
    TEST_COMPARE_ARRAYS( av2, av );
  }

  {
    out << "\nTest shallow implicit conversion from ArrayView<const T> to ArrayView<const T> ...\n";
    ArrayView<const T> cav2(cav);
    TEST_COMPARE_ARRAYS( cav2, cav );
  }

  {
    out << "\nTest shallow implicit conversion from ArrayView<const T> to ArrayView<T> ...\n";
    ArrayView<const T> cav2(av);
    TEST_COMPARE_ARRAYS( cav2, av );
  }

  {
    out << "\nTest shallow implicit conversion from std::vector<T> to ArrayView<T> ...\n";
    std::vector<T> v2 = Teuchos::createVector(cav);
    ArrayView<T> cav2(v2);
    TEST_COMPARE_ARRAYS( cav2, av );
  }

  {
    out << "\nTest shallow implicit conversion from const std::vector<T> to ArrayView<const T> ...\n";
    const std::vector<T> v2 = Teuchos::createVector(cav);
    ArrayView<const T> cav2(v2);
    TEST_COMPARE_ARRAYS( cav2, av );
  }

  {
    // Try to copy construct from ArrayView<const T> to ArrayView<T> ..
#ifdef SHOW_INVALID_COPY_CONSTRUCTION
    ArrayView<T> cav2(cav); // should not compile!
#endif
  }

  {
    out << "\ntest assign(...) ...\n";
    std::vector<T> v2(n);
    ArrayView<T> av2(v2);
    av2.assign(av);
    TEST_COMPARE_ARRAYS( v2, v );
  }

  //
  out << "\nB) Test element access ...\n";
  //


  TEST_EQUALITY_CONST( av.front(), as<T>(0) );
  TEST_EQUALITY( av.back(), as<T>(n-1) );
  TEST_EQUALITY_CONST( cav.front(), as<T>(0) );
  TEST_EQUALITY( cav.back(), as<T>(n-1) );
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_THROW( av[-1], Teuchos::RangeError );
  TEST_THROW( av[n], Teuchos::RangeError );
  TEST_THROW( cav[-1], Teuchos::RangeError );
  TEST_THROW( cav[n], Teuchos::RangeError );
#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  //
  out << "\nC) Test iterator access ...\n";
  //


  {
    out << "\nTest non-const forward iterator access ...\n";
    std::vector<T> v2(n);
    ArrayView<T> av2(v2);
    typedef typename ArrayView<T>::iterator iter_t;
    iter_t iter = av2.begin();
    for ( int i = 0; iter != av2.end(); ++i )
      *iter++ = i;
    TEST_COMPARE_ARRAYS( v2, v );
  }

  {
    out << "\nTest const forward iterator access ... ";
    bool local_success = true;
    typedef typename ArrayView<const T>::iterator iter_t;
    const ArrayView<const T> cav2 = av.getConst();
    iter_t iter = cav2.begin();
    for ( int i = 0; i < n; ++i, ++iter ) {
      TEST_ARRAY_ELE_EQUALITY( av, i, *iter );

#ifdef SHOW_INVALID_CONST_ITER_MODIFICATION
      *iter = as<T>(i); // Should not compile!
#endif
    }
    iter = NullIteratorTraits<iter_t>::getNull();
    if (local_success) out << "passed\n";
    else success = false;
  }

  //
  out << "\nD) Test sub-views ...\n";
  //

  {
    out << "\nTest full non-const subview ...\n";
    const ArrayView<T> av2 = av(0,n);
    TEST_COMPARE_ARRAYS( av2, av );
  }

  {
    out << "\nTest full shorthand non-const subview ...\n";
    const ArrayView<T> av2 = av();
    TEST_COMPARE_ARRAYS( av2, av );
  }

  {
    out << "\nTest full const subview ...\n";
    const ArrayView<const T> cav2 = cav(0,n);
    TEST_COMPARE_ARRAYS( cav2, cav );
  }

  {
    out << "\nTest full non-const to const subview ...\n";
    const ArrayView<const T> cav2 = av(0,n);
    TEST_COMPARE_ARRAYS( cav2, cav );
  }

  {
    out << "\nTest full short-hand const subview ...\n";
    const ArrayView<const T> cav2 = cav();
    TEST_COMPARE_ARRAYS( cav2, cav );
  }

  {
    out << "\nTest non-const initial range view ...\n";
    std::vector<T> v2(n,as<T>(-1));
    const ArrayView<T> av2(v2);
    const ArrayView<T> av2_init = av2(0,n-1);
    TEST_EQUALITY( av2_init.size(), n-1 );
    av2_init.assign( av(0,n-1) );
    av2.back() = as<T>(n-1);
    TEST_COMPARE_ARRAYS( v2, v );
  }

  {
    out << "\nTest non-const final range view ...\n";
    std::vector<T> v2(n,as<T>(-1));
    const ArrayView<T> av2(v2);
    const ArrayView<T> av2_init = av2(1,n-1);
    TEST_EQUALITY( av2_init.size(), n-1 );
    av2_init.assign( av(1,n-1) );
    av2.front() = as<T>(0);
    TEST_COMPARE_ARRAYS( v2, v );
  }

  {
    out << "\nTest non-const middle range view ...\n";
    std::vector<T> v2(n,as<T>(-1));
    const ArrayView<T> av2(v2);
    const ArrayView<T> av2_init = av2(1,n-2);
    TEST_EQUALITY( av2_init.size(), n-2 );
    av2_init.assign( av(1,n-2) );
    av2.front() = as<T>(0);
    av2.back() = as<T>(n-1);
    TEST_COMPARE_ARRAYS( v2, v );
  }

  // ToDo: Test requesting views outside of valid range!

  return success;

}


//
// Main testing program
//

int main( int argc, char* argv[] )
{

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  using Teuchos::CommandLineProcessor;
	
	bool success = true;
  bool result;
 
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
 
    result = testArrayView<int>(n,*out);
    if (!result) success = false;

    result = testArrayView<float>(n,*out);
    if (!result) success = false;

    result = testArrayView<double>(n,*out);
    if (!result) success = false;

    result = testArrayView<std::complex<double> >(n,*out);
    if (!result) success = false;
 
	}
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);
 
  if (success)
    *out << "\nEnd Result: TEST PASSED" << std::endl;
  else
    *out << "\nEnd Result: TEST FAILED" << std::endl;
 
  return ( success ? 0 : 1 );
 
}
