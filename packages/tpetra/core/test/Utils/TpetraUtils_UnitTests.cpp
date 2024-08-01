// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_Util.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_as.hpp>

#include <algorithm>
#include <iterator>

using Teuchos::Array;
using Teuchos::as;
using Teuchos::CommandLineProcessor;
using Tpetra::sort2;
using Tpetra::sort3;
using std::ostream_iterator;
using std::endl;

namespace Teuchos {
  template<class Array1, class Array2>
  bool compareArraysNeg(
      const Array1 &a1, const std::string &a1_name,
      const Array2 &a2, const std::string &a2_name,
      FancyOStream &out
      )
  {
    bool success = true;
    out << "Comparing " << a1_name << " == -(" << a2_name << ") ... ";
    const int n = a1.size();
    // Compare sizes
    if (as<int>(a2.size()) != n) {
      out << "\nError, "<<a1_name<<".size() = "<<a1.size()<<" == "
        << a2_name<<".size() = "<<a2.size()<<" : failed!\n";
      return false;
    }
    // Compare elements
    for( int i = 0; i < n; ++i ) {
      const bool result = ( a1[i] == -a2[i] ); // Tests C::operator[](i) const
      if (!result) {
        out << "\nError, "<<a1_name<<"["<<i<<"] = "<<a1[i]<<" == -("
          << a2_name<<"["<<i<<"]) = -("<<a2[i]<<"): failed!\n";
        success = false;
      }
    }
    if (success) {
      out << "passed\n";
    }
    return success;
  }
}

#define TEST_COMPARE_ARRAYS_NEG( a1, a2 ) \
  { \
    const bool l_result = Teuchos::compareArraysNeg(a1,#a1,a2,#a2,out); \
    if (!l_result) success = false; \
  }

template <typename ForwardIterator>
bool tpetra_is_sorted(ForwardIterator first, ForwardIterator last)
{
  if (first == last)
    return true;
  ForwardIterator next = first;
  for (++next; next != last; first = next, ++next) {
    if (*next < *first)
      return false;
  }
  return true;
}

namespace {

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
  }

  //
  // UNIT TESTS
  //


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( TpetraUtils, Sort2, T1, T2 )
  {
    const int NUMENTRIES = 20;
    // typedef typename Array<T1>::iterator IT1;
    // typedef typename Array<T2>::iterator IT2;
    Array<T1> arr1(NUMENTRIES);
    Array<T2> arr2(NUMENTRIES);
    for (Teuchos_Ordinal i=0; i<NUMENTRIES; ++i) {
      arr1[i] = -Teuchos::as<T1>(i);
      arr2[i] =  Teuchos::as<T2>(i);
    }
    // currently, arr1 goes from high to low, arr2 from low to high
    // sort according to arr1
    out << "before (arr1): "; std::copy( arr1.begin(), arr1.end(), ostream_iterator<T1>(out, " ")); out << endl;
    out << "       (arr2): "; std::copy( arr2.begin(), arr2.end(), ostream_iterator<T2>(out, " ")); out << endl;
    sort2( arr1.begin(), arr1.end(), arr2.begin() );
    out << "after  (arr1): "; std::copy( arr1.begin(), arr1.end(), ostream_iterator<T1>(out, " ")); out << endl;
    out << "       (arr2): "; std::copy( arr2.begin(), arr2.end(), ostream_iterator<T2>(out, " ")); out << endl;
    TEST_EQUALITY_CONST( tpetra_is_sorted(arr1.begin(), arr1.end()), true );
    TEST_EQUALITY_CONST( tpetra_is_sorted(arr2.begin(), arr2.end()), false );
    TEST_COMPARE_ARRAYS_NEG( arr1, arr2 )
    // currently, arr2 goes from high to low, arr1 from low to high
    // sort according to arr2
    out << "before (arr1): "; std::copy( arr1.begin(), arr1.end(), ostream_iterator<T1>(out, " ")); out << endl;
    out << "       (arr2): "; std::copy( arr2.begin(), arr2.end(), ostream_iterator<T2>(out, " ")); out << endl;
    sort2( arr2.begin(), arr2.end(), arr1.begin() );
    out << "after  (arr1): "; std::copy( arr1.begin(), arr1.end(), ostream_iterator<T1>(out, " ")); out << endl;
    out << "       (arr2): "; std::copy( arr2.begin(), arr2.end(), ostream_iterator<T2>(out, " ")); out << endl;
    TEST_EQUALITY_CONST( tpetra_is_sorted(arr2.begin(), arr2.end()), true );
    TEST_EQUALITY_CONST( tpetra_is_sorted(arr1.begin(), arr1.end()), false );
    TEST_COMPARE_ARRAYS_NEG( arr1, arr2 )
  }


  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( TpetraUtils, Sort3, T1, T2 )
  {
    const int NUMENTRIES = 20;
    // typedef typename Array<T1>::iterator IT1;
    // typedef typename Array<T2>::iterator IT2;
    Array<T1> arr1(NUMENTRIES),
              arr3(NUMENTRIES);
    Array<T2> arr2(NUMENTRIES);
    for (Teuchos_Ordinal i=0; i<NUMENTRIES; ++i) {
      arr1[i] = -Teuchos::as<T1>(i);
      arr2[i] =  Teuchos::as<T2>(i);
      arr3[i] = -Teuchos::as<T1>(i);
    }
    // currently, arr1/arr3 goes from high to low, arr2 from low to high
    // sort according to arr1
    out << "before (arr1): "; std::copy( arr1.begin(), arr1.end(), ostream_iterator<T1>(out, " ")); out << endl;
    out << "       (arr2): "; std::copy( arr2.begin(), arr2.end(), ostream_iterator<T2>(out, " ")); out << endl;
    out << "       (arr3): "; std::copy( arr3.begin(), arr3.end(), ostream_iterator<T1>(out, " ")); out << endl;
    sort3( arr1.begin(), arr1.end(), arr2.begin(), arr3.begin() );
    out << "after  (arr1): "; std::copy( arr1.begin(), arr1.end(), ostream_iterator<T1>(out, " ")); out << endl;
    out << "       (arr2): "; std::copy( arr2.begin(), arr2.end(), ostream_iterator<T2>(out, " ")); out << endl;
    out << "       (arr3): "; std::copy( arr3.begin(), arr3.end(), ostream_iterator<T1>(out, " ")); out << endl;
    TEST_EQUALITY_CONST( tpetra_is_sorted(arr1.begin(), arr1.end()), true );
    TEST_EQUALITY_CONST( tpetra_is_sorted(arr2.begin(), arr2.end()), false );
    TEST_EQUALITY_CONST( tpetra_is_sorted(arr3.begin(), arr3.end()), true );
    TEST_COMPARE_ARRAYS( arr1, arr3 );
    TEST_COMPARE_ARRAYS_NEG( arr1, arr2 )
    // currently, arr2 goes from high to low, arr1/arr3 from low to high
    // sort according to arr2
    out << "before (arr1): "; std::copy( arr1.begin(), arr1.end(), ostream_iterator<T1>(out, " ")); out << endl;
    out << "       (arr2): "; std::copy( arr2.begin(), arr2.end(), ostream_iterator<T2>(out, " ")); out << endl;
    out << "       (arr3): "; std::copy( arr3.begin(), arr3.end(), ostream_iterator<T1>(out, " ")); out << endl;
    sort3( arr2.begin(), arr2.end(), arr1.begin(), arr3.begin() );
    out << "after  (arr1): "; std::copy( arr1.begin(), arr1.end(), ostream_iterator<T1>(out, " ")); out << endl;
    out << "       (arr2): "; std::copy( arr2.begin(), arr2.end(), ostream_iterator<T2>(out, " ")); out << endl;
    out << "       (arr3): "; std::copy( arr3.begin(), arr3.end(), ostream_iterator<T1>(out, " ")); out << endl;
    TEST_EQUALITY_CONST( tpetra_is_sorted(arr1.begin(), arr1.end()), false );
    TEST_EQUALITY_CONST( tpetra_is_sorted(arr2.begin(), arr2.end()), true );
    TEST_EQUALITY_CONST( tpetra_is_sorted(arr3.begin(), arr3.end()), false );
    TEST_COMPARE_ARRAYS( arr1, arr3 );
    TEST_COMPARE_ARRAYS_NEG( arr1, arr2 )
    // currently, arr1/arr3 goes from high to low, arr2 from low to high
    // sort according to arr3
    out << "before (arr1): "; std::copy( arr1.begin(), arr1.end(), ostream_iterator<T1>(out, " ")); out << endl;
    out << "       (arr2): "; std::copy( arr2.begin(), arr2.end(), ostream_iterator<T2>(out, " ")); out << endl;
    out << "       (arr3): "; std::copy( arr3.begin(), arr3.end(), ostream_iterator<T1>(out, " ")); out << endl;
    sort3( arr3.begin(), arr3.end(), arr2.begin(), arr1.begin() );
    out << "after  (arr1): "; std::copy( arr1.begin(), arr1.end(), ostream_iterator<T1>(out, " ")); out << endl;
    out << "       (arr2): "; std::copy( arr2.begin(), arr2.end(), ostream_iterator<T2>(out, " ")); out << endl;
    out << "       (arr3): "; std::copy( arr3.begin(), arr3.end(), ostream_iterator<T1>(out, " ")); out << endl;
    TEST_EQUALITY_CONST( tpetra_is_sorted(arr1.begin(), arr1.end()), true );
    TEST_EQUALITY_CONST( tpetra_is_sorted(arr2.begin(), arr2.end()), false );
    TEST_EQUALITY_CONST( tpetra_is_sorted(arr3.begin(), arr3.end()), true );
    TEST_COMPARE_ARRAYS( arr1, arr3 );
    TEST_COMPARE_ARRAYS_NEG( arr1, arr2 )
  }


  //
  // INSTANTIATIONS
  //

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_T1_T2( T1, T2 ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( TpetraUtils, Sort2, T1, T2 ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( TpetraUtils, Sort3, T1, T2 )

# ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
     UNIT_TEST_GROUP_T1_T2(int,double)
# else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD
     UNIT_TEST_GROUP_T1_T2(double,int)
     UNIT_TEST_GROUP_T1_T2(int,double)
     UNIT_TEST_GROUP_T1_T2(double,double)
# endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

}


