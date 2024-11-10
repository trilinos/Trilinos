// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

// This test requires C++11 (for static_assert), so why not use the
// standard type traits
#include <type_traits>
#include <utility>
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TestingHelpers.hpp"

#include "Sacado_mpl_vector.hpp"
#include "Sacado_mpl_find.hpp"
#include "Sacado_mpl_size.hpp"

// These tests are all compile-time tests, so if the test compiles,
// it passes...

TEUCHOS_UNIT_TEST( MPL, Vector )
{
  using std::is_same;
  using Sacado::mpl::vector;
  using Sacado::mpl::size;
  using Sacado::mpl::push_back;
  using Sacado::mpl::at;
  using Sacado::mpl::find;

  // Some types to stick into a vector -- make some of these non-class types
  struct T1 {};
  typedef int T2;
  struct T3 {};
  struct T4 {};
  typedef char T5;

  // length-0 vector
  typedef vector<> vec0;
  static_assert( size<vec0>::value == 0, "" );

  // length-1 vector -- test push_back, at, size
  typedef typename push_back<vec0,T1>::type vec1;
  static_assert( is_same< vec1, vector<T1> >::value, "");
  static_assert( is_same< at<vec1,0>::type, T1>::value, "");
  static_assert( size<vec1>::value == 1, "" );

  // length-2 vector -- test push_back, at, size
  typedef typename push_back<vec1,T2>::type vec2;
  static_assert( is_same< vec2, vector<T1,T2> >::value, "");
  static_assert( is_same< at<vec2,0>::type, T1>::value, "");
  static_assert( is_same< at<vec2,1>::type, T2>::value, "");
  static_assert( size<vec2>::value == 2, "" );

  // length-3 vector -- test push_back, at, size
  typedef typename push_back<vec2,T3>::type vec3;
  static_assert( is_same< vec3, vector<T1,T2,T3> >::value, "");
  static_assert( is_same< at<vec3,0>::type, T1>::value, "");
  static_assert( is_same< at<vec3,1>::type, T2>::value, "");
  static_assert( is_same< at<vec3,2>::type, T3>::value, "");
  static_assert( size<vec3>::value == 3, "" );

  // length-4 vector -- test push_back, at, size
  typedef typename push_back<vec3,T4>::type vec4;
  static_assert( is_same< vec4, vector<T1,T2,T3,T4> >::value, "");
  static_assert( is_same< at<vec4,0>::type, T1>::value, "");
  static_assert( is_same< at<vec4,1>::type, T2>::value, "");
  static_assert( is_same< at<vec4,2>::type, T3>::value, "");
  static_assert( is_same< at<vec4,3>::type, T4>::value, "");
  static_assert( size<vec4>::value == 4, "" );

  // length-5 vector -- test push_back, at, size
  typedef typename push_back<vec4,T5>::type vec5;
  static_assert( is_same< vec5, vector<T1,T2,T3,T4,T5> >::value, "");
  static_assert( is_same< at<vec5,0>::type, T1>::value, "");
  static_assert( is_same< at<vec5,1>::type, T2>::value, "");
  static_assert( is_same< at<vec5,2>::type, T3>::value, "");
  static_assert( is_same< at<vec5,3>::type, T4>::value, "");
  static_assert( is_same< at<vec5,4>::type, T5>::value, "");
  static_assert( size<vec5>::value == 5, "" );

  // length-6 vector -- test push_back, at, size
  typedef typename push_back<vec5,T5>::type vec6;
  static_assert( is_same< vec6, vector<T1,T2,T3,T4,T5,T5> >::value, "");
  static_assert( is_same< at<vec6,0>::type, T1>::value, "");
  static_assert( is_same< at<vec6,1>::type, T2>::value, "");
  static_assert( is_same< at<vec6,2>::type, T3>::value, "");
  static_assert( is_same< at<vec6,3>::type, T4>::value, "");
  static_assert( is_same< at<vec6,4>::type, T5>::value, "");
  static_assert( is_same< at<vec6,5>::type, T5>::value, "");
  static_assert( size<vec6>::value == 6, "" );

  // length-7 vector -- test push_back, at, size
  typedef typename push_back<vec6,T4>::type vec7;
  static_assert( is_same< vec7, vector<T1,T2,T3,T4,T5,T5,T4> >::value, "");
  static_assert( is_same< at<vec7,0>::type, T1>::value, "");
  static_assert( is_same< at<vec7,1>::type, T2>::value, "");
  static_assert( is_same< at<vec7,2>::type, T3>::value, "");
  static_assert( is_same< at<vec7,3>::type, T4>::value, "");
  static_assert( is_same< at<vec7,4>::type, T5>::value, "");
  static_assert( is_same< at<vec7,5>::type, T5>::value, "");
  static_assert( is_same< at<vec7,6>::type, T4>::value, "");
  static_assert( size<vec7>::value == 7, "" );

  // length-8 vector -- test push_back, at, size
  typedef typename push_back<vec7,T3>::type vec8;
  static_assert( is_same< vec8, vector<T1,T2,T3,T4,T5,T5,T4,T3> >::value, "");
  static_assert( is_same< vec7, vector<T1,T2,T3,T4,T5,T5,T4> >::value, "");
  static_assert( is_same< at<vec8,0>::type, T1>::value, "");
  static_assert( is_same< at<vec8,1>::type, T2>::value, "");
  static_assert( is_same< at<vec8,2>::type, T3>::value, "");
  static_assert( is_same< at<vec8,3>::type, T4>::value, "");
  static_assert( is_same< at<vec8,4>::type, T5>::value, "");
  static_assert( is_same< at<vec8,5>::type, T5>::value, "");
  static_assert( is_same< at<vec8,6>::type, T4>::value, "");
  static_assert( is_same< at<vec8,7>::type, T3>::value, "");
  static_assert( size<vec8>::value == 8, "" );

  // length-9 vector -- test push_back, at, size
  typedef typename push_back<vec8,T2>::type vec9;
  static_assert( is_same< vec9, vector<T1,T2,T3,T4,T5,T5,T4,T3,T2> >::value, "");
  static_assert( is_same< at<vec9,0>::type, T1>::value, "");
  static_assert( is_same< at<vec9,1>::type, T2>::value, "");
  static_assert( is_same< at<vec9,2>::type, T3>::value, "");
  static_assert( is_same< at<vec9,3>::type, T4>::value, "");
  static_assert( is_same< at<vec9,4>::type, T5>::value, "");
  static_assert( is_same< at<vec9,5>::type, T5>::value, "");
  static_assert( is_same< at<vec9,6>::type, T4>::value, "");
  static_assert( is_same< at<vec9,7>::type, T3>::value, "");
  static_assert( is_same< at<vec9,8>::type, T2>::value, "");
  static_assert( size<vec9>::value == 9, "" );

  // length-10 vector -- test push_back, at, size
  typedef typename push_back<vec9,T1>::type vec10;
  static_assert( is_same< vec10, vector<T1,T2,T3,T4,T5,T5,T4,T3,T2,T1> >::value, "");
  static_assert( is_same< at<vec10,0>::type, T1>::value, "");
  static_assert( is_same< at<vec10,1>::type, T2>::value, "");
  static_assert( is_same< at<vec10,2>::type, T3>::value, "");
  static_assert( is_same< at<vec10,3>::type, T4>::value, "");
  static_assert( is_same< at<vec10,4>::type, T5>::value, "");
  static_assert( is_same< at<vec10,5>::type, T5>::value, "");
  static_assert( is_same< at<vec10,6>::type, T4>::value, "");
  static_assert( is_same< at<vec10,7>::type, T3>::value, "");
  static_assert( is_same< at<vec10,8>::type, T2>::value, "");
  static_assert( is_same< at<vec10,9>::type, T1>::value, "");
  static_assert( size<vec10>::value == 10, "" );

  // length-11 vector -- test push_back, at, size
  typedef typename push_back<vec10,T1>::type vec11;
  static_assert( is_same< vec11, vector<T1,T2,T3,T4,T5,T5,T4,T3,T2,T1,T1> >::value, "");
  static_assert( is_same< at<vec11,0>::type, T1>::value, "");
  static_assert( is_same< at<vec11,1>::type, T2>::value, "");
  static_assert( is_same< at<vec11,2>::type, T3>::value, "");
  static_assert( is_same< at<vec11,3>::type, T4>::value, "");
  static_assert( is_same< at<vec11,4>::type, T5>::value, "");
  static_assert( is_same< at<vec11,5>::type, T5>::value, "");
  static_assert( is_same< at<vec11,6>::type, T4>::value, "");
  static_assert( is_same< at<vec11,7>::type, T3>::value, "");
  static_assert( is_same< at<vec11,8>::type, T2>::value, "");
  static_assert( is_same< at<vec11,9>::type, T1>::value, "");
  static_assert( is_same< at<vec11,10>::type, T1>::value, "");
  static_assert( size<vec11>::value == 11, "" );

  // length-12 vector -- test push_back, at, size
  typedef typename push_back<vec11,T2>::type vec12;
  static_assert( is_same< vec12, vector<T1,T2,T3,T4,T5,T5,T4,T3,T2,T1,T1,T2> >::value, "");
  static_assert( is_same< at<vec12,0>::type, T1>::value, "");
  static_assert( is_same< at<vec12,1>::type, T2>::value, "");
  static_assert( is_same< at<vec12,2>::type, T3>::value, "");
  static_assert( is_same< at<vec12,3>::type, T4>::value, "");
  static_assert( is_same< at<vec12,4>::type, T5>::value, "");
  static_assert( is_same< at<vec12,5>::type, T5>::value, "");
  static_assert( is_same< at<vec12,6>::type, T4>::value, "");
  static_assert( is_same< at<vec12,7>::type, T3>::value, "");
  static_assert( is_same< at<vec12,8>::type, T2>::value, "");
  static_assert( is_same< at<vec12,9>::type, T1>::value, "");
  static_assert( is_same< at<vec12,10>::type, T1>::value, "");
  static_assert( is_same< at<vec12,11>::type, T2>::value, "");
  static_assert( size<vec12>::value == 12, "" );

  // length-13 vector -- test push_back, at, size
  typedef typename push_back<vec12,T3>::type vec13;
  static_assert( is_same< vec13, vector<T1,T2,T3,T4,T5,T5,T4,T3,T2,T1,T1,T2,T3> >::value, "");
  static_assert( is_same< at<vec13,0>::type, T1>::value, "");
  static_assert( is_same< at<vec13,1>::type, T2>::value, "");
  static_assert( is_same< at<vec13,2>::type, T3>::value, "");
  static_assert( is_same< at<vec13,3>::type, T4>::value, "");
  static_assert( is_same< at<vec13,4>::type, T5>::value, "");
  static_assert( is_same< at<vec13,5>::type, T5>::value, "");
  static_assert( is_same< at<vec13,6>::type, T4>::value, "");
  static_assert( is_same< at<vec13,7>::type, T3>::value, "");
  static_assert( is_same< at<vec13,8>::type, T2>::value, "");
  static_assert( is_same< at<vec13,9>::type, T1>::value, "");
  static_assert( is_same< at<vec13,10>::type, T1>::value, "");
  static_assert( is_same< at<vec13,11>::type, T2>::value, "");
  static_assert( is_same< at<vec13,12>::type, T3>::value, "");
  static_assert( size<vec13>::value == 13, "" );

  // length-14 vector -- test push_back, at, size
  typedef typename push_back<vec13,T4>::type vec14;
  static_assert( is_same< vec14, vector<T1,T2,T3,T4,T5,T5,T4,T3,T2,T1,T1,T2,T3,T4> >::value, "");
  static_assert( is_same< at<vec14,0>::type, T1>::value, "");
  static_assert( is_same< at<vec14,1>::type, T2>::value, "");
  static_assert( is_same< at<vec14,2>::type, T3>::value, "");
  static_assert( is_same< at<vec14,3>::type, T4>::value, "");
  static_assert( is_same< at<vec14,4>::type, T5>::value, "");
  static_assert( is_same< at<vec14,5>::type, T5>::value, "");
  static_assert( is_same< at<vec14,6>::type, T4>::value, "");
  static_assert( is_same< at<vec14,7>::type, T3>::value, "");
  static_assert( is_same< at<vec14,8>::type, T2>::value, "");
  static_assert( is_same< at<vec14,9>::type, T1>::value, "");
  static_assert( is_same< at<vec14,10>::type, T1>::value, "");
  static_assert( is_same< at<vec14,11>::type, T2>::value, "");
  static_assert( is_same< at<vec14,12>::type, T3>::value, "");
  static_assert( is_same< at<vec14,13>::type, T4>::value, "");
  static_assert( size<vec14>::value == 14, "" );

  // length-15 vector -- test push_back, at, size
  typedef typename push_back<vec14,T5>::type vec15;
  static_assert( is_same< vec15, vector<T1,T2,T3,T4,T5,T5,T4,T3,T2,T1,T1,T2,T3,T4,T5> >::value, "");
  static_assert( is_same< at<vec15,0>::type, T1>::value, "");
  static_assert( is_same< at<vec15,1>::type, T2>::value, "");
  static_assert( is_same< at<vec15,2>::type, T3>::value, "");
  static_assert( is_same< at<vec15,3>::type, T4>::value, "");
  static_assert( is_same< at<vec15,4>::type, T5>::value, "");
  static_assert( is_same< at<vec15,5>::type, T5>::value, "");
  static_assert( is_same< at<vec15,6>::type, T4>::value, "");
  static_assert( is_same< at<vec15,7>::type, T3>::value, "");
  static_assert( is_same< at<vec15,8>::type, T2>::value, "");
  static_assert( is_same< at<vec15,9>::type, T1>::value, "");
  static_assert( is_same< at<vec15,10>::type, T1>::value, "");
  static_assert( is_same< at<vec15,11>::type, T2>::value, "");
  static_assert( is_same< at<vec15,12>::type, T3>::value, "");
  static_assert( is_same< at<vec15,13>::type, T4>::value, "");
  static_assert( is_same< at<vec15,14>::type, T5>::value, "");
  static_assert( size<vec15>::value == 15, "" );

  // The implementation now uses variatic templates, so there
  // is no hard limit on the length of mpl::vector.  However
  // sizes up to 15 still seems like a good length to test.

  // Check find
  static_assert( Sacado::mpl::find<vec5, T1>::value == 0, "" );
  static_assert( Sacado::mpl::find<vec5, T2>::value == 1, "" );
  static_assert( Sacado::mpl::find<vec5, T3>::value == 2, "" );
  static_assert( Sacado::mpl::find<vec5, T4>::value == 3, "" );
  static_assert( Sacado::mpl::find<vec5, T5>::value == 4, "" );
  static_assert( Sacado::mpl::find<vec10, T1>::value == 0, "" );
  static_assert( Sacado::mpl::find<vec10, T2>::value == 1, "" );
  static_assert( Sacado::mpl::find<vec10, T3>::value == 2, "" );
  static_assert( Sacado::mpl::find<vec10, T4>::value == 3, "" );
  static_assert( Sacado::mpl::find<vec10, T5>::value == 4, "" );
  static_assert( Sacado::mpl::find<vec15, T1>::value == 0, "" );
  static_assert( Sacado::mpl::find<vec15, T2>::value == 1, "" );
  static_assert( Sacado::mpl::find<vec15, T3>::value == 2, "" );
  static_assert( Sacado::mpl::find<vec15, T4>::value == 3, "" );
  static_assert( Sacado::mpl::find<vec15, T5>::value == 4, "" );

  success = true;
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
