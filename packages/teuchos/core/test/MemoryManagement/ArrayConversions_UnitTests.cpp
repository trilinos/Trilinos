// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#include "Teuchos_ArrayConversions.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Array_UnitTest_helpers.hpp"
#include "ArrayConversions_UnitTest_helpers.hpp"
#include "TestClasses.hpp"

namespace {


using ArrayUnitTestHelpers::n;
using ArrayConversionsUnitTestHelpers::generateArrayRcp;
using ArrayConversionsUnitTestHelpers::generateArrayRcpGen;
using ArrayConversionsUnitTestHelpers::testArrayViewInput;
using ArrayConversionsUnitTestHelpers::testArrayViewOutput;
using Teuchos::arrayPtrConv;
using Teuchos::arrayRcpConv;
using Teuchos::arrayViewPtrConv;
using Teuchos::arrayViewRcpConv;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::Ptr;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::as;
using Teuchos::implicit_ptr_cast;


// Verify generateArrayRcp works correctly
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, generateArrayRcp, T )
{
  const Array<RCP<T> > a_in = generateArrayRcp<T>(n);
  TEST_EQUALITY_CONST( as<Teuchos_Ordinal>(a_in.size()), n );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY_CONST( Teuchos::is_null(a_in[i]), false );
    TEST_EQUALITY_CONST( *a_in[i], as<T>(i) );
  }
}


// Verify testArrayViewInput works correctly
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, testArrayViewInput, T )
{
  //typedef Teuchos::ScalarTraits<T> ST; // unused
  const Array<RCP<T> > a_data = generateArrayRcp<T>(n);
  Array<Ptr<const T> > a_in(n);
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    a_in[i] = a_data[i].ptr();
  }
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY_CONST( Teuchos::is_null(a_in[i]), false );
  }
  T a_out = testArrayViewInput<T>(a_in);
  TEST_EQUALITY_CONST( a_out, as<T>(n*(n-1)/2) );
}


// Verify testArrayViewOutput works correctly
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, testArrayViewOutput, T )
{
  typedef Teuchos::ScalarTraits<T> ST;
  const Array<RCP<T> > a_data = generateArrayRcp<T>(n);
  Array<Ptr<T> > a_out;
  a_out.reserve(n);
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    *a_data[i] = ST::zero();
    a_out.push_back( a_data[i].ptr() );
  }
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY_CONST( Teuchos::is_null(a_out[i]), false );
  }
  testArrayViewOutput<T>(a_out);
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY_CONST( *a_out[i], as<T>(i) );
  }

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayPtrConv_RcpNonconst_to_PtrConst, T )
{
  const Array<RCP<T> > a_in = generateArrayRcp<T>(n);
  const Array<Ptr<const T> > a_out = arrayPtrConv<const T>(a_in);
  TEST_EQUALITY( a_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( a_out[i].get(), a_in[i].get() );
    TEST_EQUALITY( *a_out[i], *a_in[i] );
  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayPtrConv_RcpNonconst_to_PtrNonconst, T )
{
  const Array<RCP<T> > a_in = generateArrayRcp<T>(n);
  const Array<Ptr<T> > a_out = arrayPtrConv<T>(a_in);
  TEST_EQUALITY( a_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( a_out[i].get(), a_in[i].get() );
    TEST_EQUALITY( *a_out[i], *a_in[i] );
  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayPtrConv_PtrNonconst_to_PtrNonconst, T )
{
  const Array<RCP<T> > a_in = generateArrayRcp<T>(n);
  const Array<Ptr<T> > a1_out = arrayPtrConv<T>(a_in);
  const Array<Ptr<T> > a2_out = arrayPtrConv<T>(a1_out());
  TEST_COMPARE_ARRAYS( a2_out, a1_out );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayPtrConv_PtrNonconst_to_PtrConst, T )
{
  const Array<RCP<T> > a_in = generateArrayRcp<T>(n);
  const Array<Ptr<T> > a1_out = arrayPtrConv<T>(a_in);
  const Array<Ptr<const T> > a2_out = arrayPtrConv<const T>(a1_out());
  TEST_COMPARE_ARRAYS( a2_out, a1_out );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayPtrConv_PtrConst_to_PtrConst, T )
{
  const Array<RCP<T> > a_in = generateArrayRcp<T>(n);
  const Array<Ptr<const T> > a1_out = arrayPtrConv<const T>(a_in);
  const Array<Ptr<const T> > a2_out = arrayPtrConv<const T>(a1_out());
  TEST_COMPARE_ARRAYS( a2_out, a1_out );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayPtrConv_PassConst, T )
{
  Array<RCP<T> > a_in = generateArrayRcp<T>(n);
  T a = testArrayViewInput<T>(arrayPtrConv<const T>(a_in));
  T a_exact = as<T>(n*(n-1)/2);
  TEST_EQUALITY( a, a_exact );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayPtrConv_PassNonconst, T )
{
  typedef Teuchos::ScalarTraits<T> ST;
  Array<RCP<T> > a_out = generateArrayRcp<T>(n);
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    *a_out[i] = ST::zero();
  }
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY_CONST( *a_out[i], ST::zero() );
  }
  testArrayViewOutput<T>(arrayPtrConv<T>(a_out));
  TEST_EQUALITY_CONST( as<Teuchos_Ordinal>(a_out.size()), n );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY_CONST( *a_out[i], as<T>(i) );
  }
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayRcpConv_RcpNonconst_to_RcpNonconst, T )
{
  const Array<RCP<T> > a_in = generateArrayRcp<T>(n);
  const Array<RCP<T> > a1_out = arrayRcpConv<T>(a_in);
  TEST_COMPARE_ARRAYS( a1_out, a_in );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayRcpConv_RcpNonconst_to_RcpConst, T )
{
  const Array<RCP<T> > a_in = generateArrayRcp<T>(n);
  const Array<RCP<const T> > a1_out = arrayRcpConv<const T>(a_in);
  TEST_COMPARE_ARRAYS( a1_out, a_in );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayRcpConv_RcpConst_to_RcpConst, T )
{
  const Array<RCP<T> > a_in = generateArrayRcp<T>(n);
  const Array<RCP<const T> > a1_out = arrayRcpConv<const T>(a_in);
  const Array<RCP<const T> > a2_out = arrayRcpConv<const T>(a1_out);
  TEST_COMPARE_ARRAYS( a2_out, a1_out );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayConstPtrConstCast_nonconst_to_const,
  T )
{
  Array<RCP<T> > a_rcp = generateArrayRcp<T>(n);
  Array<Ptr<T> > a_ptr = arrayPtrConv<T>(a_rcp);
  const ArrayView<const Ptr<T> > av_ptr_nonconst = a_ptr();
  const ArrayView<const Ptr<const T> > av_ptr_const =
    Teuchos::arrayConstPtrConstCast(av_ptr_nonconst);
  TEST_COMPARE_ARRAYS( av_ptr_nonconst, av_ptr_const );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayConstPtrConstCast_const_to_const,
  T )
{
  Array<RCP<T> > a_rcp = generateArrayRcp<T>(n);
  Array<Ptr<const T> > a_ptr = arrayPtrConv<const T>(a_rcp);
  const ArrayView<const Ptr<const T> > av_ptr_const1 = a_ptr();
  const ArrayView<const Ptr<const T> > av_ptr_const2 =
    Teuchos::arrayConstPtrConstCast(av_ptr_const1);
  TEST_COMPARE_ARRAYS( av_ptr_const1, av_ptr_const2 );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayConstRcpConstCast_nonconst_to_const,
  T )
{
  Array<RCP<T> > a_rcp = generateArrayRcp<T>(n);
  const ArrayView<const RCP<T> > av_rcp_nonconst = a_rcp();
  const ArrayView<const RCP<const T> > av_rcp_const =
    Teuchos::arrayConstRcpConstCast(av_rcp_nonconst);
  TEST_COMPARE_ARRAYS( av_rcp_nonconst, av_rcp_const );
}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ArrayConversions, arrayConstRcpConstCast_const_to_const,
  T )
{
  Array<RCP<T> > a_rcp_orig = generateArrayRcp<T>(n);
  Array<RCP<const T> > a_rcp = arrayRcpConv<const T>(a_rcp_orig);
  const ArrayView<const RCP<const T> > av_rcp_const1 = a_rcp();
  const ArrayView<const RCP<const T> > av_rcp_const2 =
    Teuchos::arrayConstRcpConstCast(av_rcp_const1);
  TEST_COMPARE_ARRAYS( av_rcp_const1, av_rcp_const2 );
}


#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#  define DEBUG_UNIT_TEST_GROUP( T )

#else // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#  define DEBUG_UNIT_TEST_GROUP( T )

#endif // HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

#define UNIT_TEST_GROUP( T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, generateArrayRcp, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, testArrayViewInput, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, testArrayViewOutput, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayPtrConv_RcpNonconst_to_PtrConst, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayPtrConv_RcpNonconst_to_PtrNonconst, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayPtrConv_PtrNonconst_to_PtrNonconst, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayPtrConv_PtrNonconst_to_PtrConst, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayPtrConv_PtrConst_to_PtrConst, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayRcpConv_RcpConst_to_RcpConst, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayPtrConv_PassConst, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayPtrConv_PassNonconst, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayRcpConv_RcpNonconst_to_RcpNonconst, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayRcpConv_RcpNonconst_to_RcpConst, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayConstPtrConstCast_nonconst_to_const, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayConstPtrConstCast_const_to_const, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayConstRcpConstCast_nonconst_to_const, T ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( ArrayConversions, arrayConstRcpConstCast_const_to_const, T ) \
  DEBUG_UNIT_TEST_GROUP( T )


UNIT_TEST_GROUP(Teuchos_Ordinal)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)



TEUCHOS_UNIT_TEST( ArrayConversions, arrayPtrConv_RcpNonconstDerived_to_PtrNonconstBase)
{
  const Array<RCP<C> > a_in = generateArrayRcpGen<C>(n);
  const Array<Ptr<A> > a_out = arrayPtrConv<A>(a_in);
  TEST_EQUALITY( a_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( &*a_out[i], implicit_ptr_cast<A>(&*a_in[i]) );
  }
}


TEUCHOS_UNIT_TEST( ArrayConversions, arrayPtrConv_RcpNonconstDerived_to_PtrConstBase)
{
  const Array<RCP<C> > a_in = generateArrayRcpGen<C>(n);
  const Array<Ptr<const A> > a_out = arrayPtrConv<const A>(a_in);
  TEST_EQUALITY( a_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( &*a_out[i], implicit_ptr_cast<const A>(&*a_in[i]) );
  }
}


TEUCHOS_UNIT_TEST( ArrayConversions, arrayPtrConv_RcpConstDerived_to_PtrConstBase)
{
  const Array<RCP<C> > a_in = generateArrayRcpGen<C>(n);
  const Array<RCP<const C> > a1_out = arrayRcpConv<const C>(a_in);
  const Array<Ptr<const A> > a2_out = arrayPtrConv<const A>(a1_out);
  TEST_EQUALITY( a2_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( &*a2_out[i], implicit_ptr_cast<const A>(&*a_in[i]) );
  }
}


TEUCHOS_UNIT_TEST( ArrayConversions, arrayRcpConv_RcpNonconstDerived_to_RcpNonconstBase)
{
  const Array<RCP<C> > a_in = generateArrayRcpGen<C>(n);
  const Array<RCP<A> > a_out = arrayRcpConv<A>(a_in);
  TEST_EQUALITY( a_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( &*a_out[i], implicit_ptr_cast<A>(&*a_in[i]) );
  }
}


TEUCHOS_UNIT_TEST( ArrayConversions, arrayRcpConv_RcpNonconstDerived_to_RcpConstBase)
{
  const Array<RCP<C> > a_in = generateArrayRcpGen<C>(n);
  const Array<RCP<const A> > a_out = arrayRcpConv<const A>(a_in);
  TEST_EQUALITY( a_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( &*a_out[i], implicit_ptr_cast<A>(&*a_in[i]) );
  }
}


TEUCHOS_UNIT_TEST( ArrayConversions, arrayRcpConv_RcpConstDerived_to_RcpConstBase)
{
  const Array<RCP<C> > a_in = generateArrayRcpGen<C>(n);
  const Array<RCP<const C> > a1_out = arrayRcpConv<const C>(a_in);
  const Array<RCP<const A> > a2_out = arrayRcpConv<const A>(a1_out);
  TEST_EQUALITY( a2_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( &*a2_out[i], implicit_ptr_cast<A>(&*a_in[i]) );
  }
}


TEUCHOS_UNIT_TEST( ArrayConversions, arrayViewPtrConv_RcpNonconstDerived_to_PtrNonconstBase)
{
  const Array<RCP<C> > a_in = generateArrayRcpGen<C>(n);
  Array<Ptr<A> > a_out(n);
  arrayViewPtrConv(a_in, a_out());
  TEST_EQUALITY( a_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( &*a_out[i], implicit_ptr_cast<A>(&*a_in[i]) );
  }
}


TEUCHOS_UNIT_TEST( ArrayConversions, arrayViewPtrConv_RcpNonconstDerived_to_PtrConstBase)
{
  const Array<RCP<C> > a_in = generateArrayRcpGen<C>(n);
  Array<Ptr<const A> > a_out(n);
  arrayViewPtrConv(a_in, a_out());
  TEST_EQUALITY( a_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( &*a_out[i], implicit_ptr_cast<const A>(&*a_in[i]) );
  }
}


TEUCHOS_UNIT_TEST( ArrayConversions, arrayViewPtrConv_RcpConstDerived_to_PtrConstBase)
{
  const Array<RCP<C> > a_in = generateArrayRcpGen<C>(n);
  Array<RCP<const C> > a1_out(n);
  arrayViewRcpConv(a_in, a1_out());
  Array<Ptr<const A> > a2_out(n);
  arrayViewPtrConv(a1_out, a2_out());
  TEST_EQUALITY( a2_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( &*a2_out[i], implicit_ptr_cast<const A>(&*a_in[i]) );
  }
}


TEUCHOS_UNIT_TEST( ArrayConversions, arrayViewRcpConv_RcpNonconstDerived_to_RcpNonconstBase)
{
  const Array<RCP<C> > a_in = generateArrayRcpGen<C>(n);
  Array<RCP<A> > a_out(n);
  arrayViewRcpConv(a_in, a_out());
  TEST_EQUALITY( a_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( &*a_out[i], implicit_ptr_cast<A>(&*a_in[i]) );
  }
}


TEUCHOS_UNIT_TEST( ArrayConversions, arrayViewRcpConv_RcpNonconstDerived_to_RcpConstBase)
{
  const Array<RCP<C> > a_in = generateArrayRcpGen<C>(n);
  Array<RCP<const A> > a_out(n);
  arrayViewRcpConv(a_in, a_out());
  TEST_EQUALITY( a_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( &*a_out[i], implicit_ptr_cast<A>(&*a_in[i]) );
  }
}


TEUCHOS_UNIT_TEST( ArrayConversions, arrayViewRcpConv_RcpConstDerived_to_RcpConstBase)
{
  const Array<RCP<C> > a_in = generateArrayRcpGen<C>(n);
  Array<RCP<const C> > a1_out(n);
  arrayViewRcpConv(a_in, a1_out());
  Array<RCP<const A> > a2_out(n);
  arrayViewRcpConv(a1_out, a2_out());
  TEST_EQUALITY( a2_out.size(), a_in.size() );
  for (Teuchos_Ordinal i=0 ; i<n ; ++i) {
    TEST_EQUALITY( &*a2_out[i], implicit_ptr_cast<A>(&*a_in[i]) );
  }
}


} // namespace
