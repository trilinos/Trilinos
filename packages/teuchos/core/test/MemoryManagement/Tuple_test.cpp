// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_Tuple.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
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


#define TEST_EQUALITY_CONST( v1, v2 )                   \
  TEUCHOS_TEST_EQUALITY_CONST( v1, v2, out, success )

#define TEST_EQUALITY( v1, v2 )                 \
  TEUCHOS_TEST_EQUALITY( v1, v2, out, success )

#define TEST_ITER_EQUALITY( iter1, iter2 )                      \
  TEUCHOS_TEST_ITER_EQUALITY( iter1, iter2, out, success )

#define TEST_ARRAY_ELE_EQUALITY( a, i, val )                            \
  TEUCHOS_TEST_ARRAY_ELE_EQUALITY( a, i, val, false, out, local_success )

#define TEST_COMPARE( v1, comp, v2 )                    \
  TEUCHOS_TEST_COMPARE( v1, comp, v2, out, success )

#define TEST_COMPARE_ARRAYS( a1, a2 )                           \
  {                                                             \
    const bool result = compareArrays(a1,#a1,a2,#a2,out);       \
    if (!result) success = false;                               \
  }

#define TEST_THROW( code, ExceptType  )                 \
  TEUCHOS_TEST_THROW( code, ExceptType, out, success  )

#define TEST_NOTHROW( code  )                   \
  TEUCHOS_TEST_NOTHROW( code, out, success  )


//
// Main templated array test function
//


template<class T, int N>
bool testTuple( Teuchos::FancyOStream &out )
{

  using Teuchos::Tuple;
  using Teuchos::tuple;
  using Teuchos::ArrayView;
  using Teuchos::arrayView;
  using Teuchos::arrayViewFromVector;
  using Teuchos::outArg;
  using Teuchos::NullIteratorTraits;
  using Teuchos::TypeNameTraits;
  using Teuchos::getConst;
  using Teuchos::as;
  typedef typename ArrayView<T>::size_type size_type;
  // mfh 03 Apr 2014: The point of the above line of code is to ensure
  // that ArrayView<T> has a public size_type typedef.  However, the
  // above line of code in isolation causes some compilers to warn
  // about a declared but unused typedef.  We deal with this by
  // declaring a variable (actually, the oxymoron "const variable") of
  // type size_type, then using the "cast to void" trick to forestall
  // compiler warnings for the declared but unused variable.  (Fun
  // fact: "oxymoron" means "sharp dull" and is itself an oxymoron.)
  // The "cast to void" trick doesn't always work, but if it doesn't,
  // it's easy to make it go away by printing it to the output stream
  // 'out'.
  const size_type arbitrarySizeTypeValue = 0;
  (void) arbitrarySizeTypeValue;

  bool success = true;

  out
    << "\n***"
    << "\n*** Testing "<<TypeNameTraits<Tuple<T,N> >::name()<<" of size = "<<N
    << "\n***\n";

  Teuchos::OSTab tab(out);

  //
  out << "\nA) Initial setup testing ...\n\n";
  //

  Tuple<T,N> t;
  TEST_EQUALITY_CONST(t.size(),N);
  for( int i = 0; i < N; ++i )
    t[i] = i; // tests non-const operator[](i)

  {
    out << "\nTest that t[i] == i ... ";
    const ArrayView<const T> cav2 = t;
    bool local_success = true;
    for( int i = 0; i < N; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( cav2, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const int n = 1;
    out << "\nTest Tuple<T,"<<n<<"> = tuple(...)\n";
    Tuple<T,n> tn = tuple<T>(0);
    TEST_EQUALITY_CONST(tn.size(),n);
    out << "Test that tn[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( tn, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const int n = 2;
    out << "\nTest Tuple<T,"<<n<<"> = tuple(...)\n";
    Tuple<T,n> tn = tuple<T>(0,1);
    TEST_EQUALITY_CONST(tn.size(),n);
    out << "Test that tn[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( tn, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const int n = 3;
    out << "\nTest Tuple<T,"<<n<<"> = tuple(...)\n";
    Tuple<T,n> tn = tuple<T>(0,1,2);
    TEST_EQUALITY_CONST(tn.size(),n);
    out << "Test that tn[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( tn, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const int n = 4;
    out << "\nTest Tuple<T,"<<n<<"> = tuple(...)\n";
    Tuple<T,n> tn = tuple<T>(0,1,2,3);
    TEST_EQUALITY_CONST(tn.size(),n);
    out << "Test that tn[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( tn, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const int n = 5;
    out << "\nTest Tuple<T,"<<n<<"> = tuple(...)\n";
    Tuple<T,n> tn = tuple<T>(0,1,2,3,4);
    TEST_EQUALITY_CONST(tn.size(),n);
    out << "Test that tn[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( tn, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const int n = 6;
    out << "\nTest Tuple<T,"<<n<<"> = tuple(...)\n";
    Tuple<T,n> tn = tuple<T>(0,1,2,3,4,5);
    TEST_EQUALITY_CONST(tn.size(),n);
    out << "Test that tn[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( tn, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const int n = 7;
    out << "\nTest Tuple<T,"<<n<<"> = tuple(...)\n";
    Tuple<T,n> tn = tuple<T>(0,1,2,3,4,5,6);
    TEST_EQUALITY_CONST(tn.size(),n);
    out << "Test that tn[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( tn, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const int n = 8;
    out << "\nTest Tuple<T,"<<n<<"> = tuple(...)\n";
    Tuple<T,n> tn = tuple<T>(0,1,2,3,4,5,6,7);
    TEST_EQUALITY_CONST(tn.size(),n);
    out << "Test that tn[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( tn, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const int n = 9;
    out << "\nTest Tuple<T,"<<n<<"> = tuple(...)\n";
    Tuple<T,n> tn = tuple<T>(0,1,2,3,4,5,6,7,8);
    TEST_EQUALITY_CONST(tn.size(),n);
    out << "Test that tn[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( tn, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const int n = 10;
    out << "\nTest Tuple<T,"<<n<<"> = tuple(...)\n";
    Tuple<T,n> tn = tuple<T>(0,1,2,3,4,5,6,7,8,9);
    TEST_EQUALITY_CONST(tn.size(),n);
    out << "Test that tn[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( tn, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const int n = 11;
    out << "\nTest Tuple<T,"<<n<<"> = tuple(...)\n";
    Tuple<T,n> tn = tuple<T>(0,1,2,3,4,5,6,7,8,9,10);
    TEST_EQUALITY_CONST(tn.size(),n);
    out << "Test that tn[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( tn, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const int n = 12;
    out << "\nTest Tuple<T,"<<n<<"> = tuple(...)\n";
    Tuple<T,n> tn = tuple<T>(0,1,2,3,4,5,6,7,8,9,10,11);
    TEST_EQUALITY_CONST(tn.size(),n);
    out << "Test that tn[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( tn, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const int n = 13;
    out << "\nTest Tuple<T,"<<n<<"> = tuple(...)\n";
    Tuple<T,n> tn = tuple<T>(0,1,2,3,4,5,6,7,8,9,10,11,12);
    TEST_EQUALITY_CONST(tn.size(),n);
    out << "Test that tn[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( tn, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const int n = 14;
    out << "\nTest Tuple<T,"<<n<<"> = tuple(...)\n";
    Tuple<T,n> tn = tuple<T>(0,1,2,3,4,5,6,7,8,9,10,11,12,13);
    TEST_EQUALITY_CONST(tn.size(),n);
    out << "Test that tn[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( tn, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    const int n = 15;
    out << "\nTest Tuple<T,"<<n<<"> = tuple(...)\n";
    Tuple<T,n> tn = tuple<T>(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14);
    TEST_EQUALITY_CONST(tn.size(),n);
    out << "Test that tn[i] == i ... ";
    bool local_success = true;
    for( int i = 0; i < n; ++i ) {
      TEST_ARRAY_ELE_EQUALITY( tn, i, as<T>(i) );
    }
    if (local_success) out << "passed\n";
    else success = false;
  }

  {
    out << "\nTest constructing Array<const T> from Tuple<T,N> ...\n";
    const ArrayView<const T> av2 = t;
    TEST_COMPARE_ARRAYS( av2, t );
  }

  // ToDo: Add more tests!

  return success;

}


//
// Main testing program
//

int main( int argc, char* argv[] ) {

  using Teuchos::CommandLineProcessor;

  bool success = true;
  bool result;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

  try {

    //
    // Read options from the commandline
    //

    CommandLineProcessor clp(false); // Don't throw exceptions

    CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);

    if ( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
      *out << "\nEnd Result: TEST FAILED" << std::endl;
      return parse_return;
    }

    *out << std::endl << Teuchos::Teuchos_Version() << std::endl;

    const int N = 8;

    result = testTuple<int,N>(*out);
    if (!result) success = false;

    result = testTuple<float,N>(*out);
    if (!result) success = false;

    result = testTuple<double,N>(*out);
    if (!result) success = false;

    result = testTuple<std::complex<double> ,N>(*out);
    if (!result) success = false;

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);

  if (success)
    *out << "\nEnd Result: TEST PASSED" << std::endl;
  else
    *out << "\nEnd Result: TEST FAILED" << std::endl;

  return ( success ? 0 : 1 );

}
