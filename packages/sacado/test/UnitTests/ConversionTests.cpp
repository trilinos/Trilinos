// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include <type_traits>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TestingHelpers.hpp"

#include "Sacado_No_Kokkos.hpp"
#include "Sacado_Fad_SimpleFad.hpp"
#include "Sacado_Tay_CacheTaylor.hpp"
#include "Sacado_mpl_apply.hpp"

// Some classes for testing std::is_convertible<From,To>
struct A {};
struct B {
  B() {}
  B(const A&) {}
};
struct C : public A {};

// Size used for all Fad types
const int global_fad_size = 10;

// Test is_convertible<From,To> behaves as expected
TEUCHOS_UNIT_TEST( Conversion, IsConvertible )
{
  const bool is_b_a = std::is_convertible<B,A>::value;
  const bool is_a_b = std::is_convertible<A,B>::value;
  const bool is_c_a = std::is_convertible<C,A>::value;
  const bool is_int_double = std::is_convertible<int,double>::value;
  const bool is_double_int = std::is_convertible<double,int>::value;
  const bool is_double_a = std::is_convertible<double,A>::value;
  TEST_EQUALITY( is_b_a, false );
  TEST_EQUALITY( is_a_b, true );
  TEST_EQUALITY( is_c_a, true );
  TEST_EQUALITY( is_int_double, true );
  TEST_EQUALITY( is_double_int, true );
  TEST_EQUALITY( is_double_a, false );
}

template <typename ad_type>
bool test_ad_conversions(Teuchos::FancyOStream& out)
{
  bool success = true;
  typedef typename Sacado::ValueType<ad_type>::type value_type;
  typedef typename Sacado::ScalarType<ad_type>::type scalar_type;

  const bool is_value_ad =
    std::is_convertible<value_type,ad_type>::value;
  const bool is_ad_value =
    std::is_convertible<ad_type,value_type>::value;
  const bool is_scalar_ad =
    std::is_convertible<scalar_type,ad_type>::value;
  const bool is_ad_scalar =
    std::is_convertible<ad_type,scalar_type>::value;
  const bool is_not_view = ! Sacado::IsView<ad_type>::value;

  const bool is_int_ad =
    std::is_convertible<value_type,ad_type>::value;

  TEST_EQUALITY( is_value_ad, is_not_view );
  TEST_EQUALITY_CONST( is_ad_value, false );
  TEST_EQUALITY( is_scalar_ad, is_not_view );
  TEST_EQUALITY_CONST( is_ad_scalar, false );
  TEST_EQUALITY( is_int_ad, is_not_view );

  // Get the type of the result of the expression 'ad_type * ad_type'
  // The use of declval gets around actually instantiation objects of type
  // ad_type.
  typedef decltype(std::declval<ad_type>()*std::declval<ad_type>()) ad_expr_type;
  typedef decltype(std::declval<value_type>()*std::declval<value_type>()) val_expr_type;

  const bool is_ad_expr_ad =
    std::is_convertible<ad_expr_type,ad_type>::value;
  const bool is_val_expr_ad =
    std::is_convertible<val_expr_type,ad_type>::value;

  TEST_EQUALITY( is_ad_expr_ad, is_not_view );
  TEST_EQUALITY( is_val_expr_ad, is_not_view );

  // typedef typename ad_expr_type::value_type ad_expr_value_type;
  // std::cout << typeid(ad_expr_value_type).name() << std::endl;

  return success;
}

// Check various AD conversions work as expected
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Conversion, ADConversions, AD )
{
  typedef AD ad_type;
  typedef typename ad_type::value_type value_type;
  typedef typename Sacado::mpl::apply< ad_type, ad_type >::type ad_ad_type;

  success = true;
  success = success && test_ad_conversions<ad_type>(out);
  success = success && test_ad_conversions<ad_ad_type>(out);

  // Check value-type expression to nested fad-fad works
  ad_type x(global_fad_size, value_type(1.5));
  for (int i=0; i<global_fad_size; ++i)
    x.fastAccessDx(i) = 2.0;
  ad_ad_type y = x + x;
  TEST_EQUALITY_CONST( y.val().val(), 3.0 );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.val().dx(i), 4.0 );
    TEST_EQUALITY_CONST( y.dx(i).val(), 0.0 );
    for (int j=0; j<global_fad_size; ++j)
      TEST_EQUALITY_CONST( y.dx(i).dx(j), 0.0 );
  }

  // Check mixed value-type/Fad expression with nested fad-fad works
  ad_ad_type z = (x + x) + y;
  TEST_EQUALITY_CONST( z.val().val(), 6.0 );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( z.val().dx(i), 8.0 );
    TEST_EQUALITY_CONST( z.dx(i).val(), 0.0 );
    for (int j=0; j<global_fad_size; ++j)
      TEST_EQUALITY_CONST( z.dx(i).dx(j), 0.0 );
  }

  // Check mix-arithmetic with int's works
  y += 1;
  TEST_EQUALITY_CONST( y.val().val(), 4.0 );
  for (int i=0; i<global_fad_size; ++i) {
    TEST_EQUALITY_CONST( y.val().dx(i), 4.0 );
    TEST_EQUALITY_CONST( y.dx(i).val(), 0.0 );
    for (int j=0; j<global_fad_size; ++j)
      TEST_EQUALITY_CONST( y.dx(i).dx(j), 0.0 );
  }
}

// Check various view conversions work as expected
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Conversion, ViewConversions, AD )
{
  typedef AD ad_type;
  typedef typename Sacado::mpl::apply< ad_type, ad_type >::type ad_ad_type;

  success = true;
  success = success && test_ad_conversions<ad_type>(out);
  success = success && test_ad_conversions<ad_ad_type>(out);

  // ad_ad_type x;
  // ad_ad_type y = x*x;
}

// Check various other conversions work as expected
// These are for types that aren't expected to be nested, but may be nesteed
// inside other Fad types
TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Conversion, OtherConversions, AD )
{
  typedef AD ad_type;
  typedef Sacado::Fad::DFad<ad_type> fad_ad_type;

  success = true;
  success = success && test_ad_conversions<ad_type>(out);
  success = success && test_ad_conversions<fad_ad_type>(out);
}

typedef Sacado::Fad::DFad<double> Fad_DFadType;
typedef Sacado::Fad::SLFad<double,global_fad_size> Fad_SLFadType;
typedef Sacado::Fad::SFad<double,global_fad_size> Fad_SFadType;
typedef Sacado::Fad::DVFad<double> Fad_DVFadType;
typedef Sacado::Fad::SimpleFad<double> Fad_SimpleFadType;
typedef Sacado::Fad::ViewFad<double,global_fad_size,1,Fad_DFadType> Fad_VFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ADConversions, Fad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ADConversions, Fad_SLFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ADConversions, Fad_SFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ADConversions, Fad_DVFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ADConversions, Fad_SimpleFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ViewConversions, Fad_VFadType )

typedef Sacado::ELRFad::DFad<double> ELRFad_DFadType;
typedef Sacado::ELRFad::SLFad<double,global_fad_size> ELRFad_SLFadType;
typedef Sacado::ELRFad::SFad<double,global_fad_size> ELRFad_SFadType;
typedef Sacado::ELRFad::ViewFad<double,global_fad_size,1,ELRFad_DFadType> ELRFad_VFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ADConversions, ELRFad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ADConversions, ELRFad_SLFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ADConversions, ELRFad_SFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ViewConversions, ELRFad_VFadType )

typedef Sacado::CacheFad::DFad<double> CacheFad_DFadType;
typedef Sacado::CacheFad::SLFad<double,global_fad_size> CacheFad_SLFadType;
typedef Sacado::CacheFad::SFad<double,global_fad_size> CacheFad_SFadType;
typedef Sacado::CacheFad::ViewFad<double,global_fad_size,1,CacheFad_DFadType> CacheFad_VFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ADConversions, CacheFad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ADConversions, CacheFad_SLFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ADConversions, CacheFad_SFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ViewConversions, CacheFad_VFadType )

typedef Sacado::ELRCacheFad::DFad<double> ELRCacheFad_DFadType;
typedef Sacado::ELRCacheFad::SLFad<double,global_fad_size> ELRCacheFad_SLFadType;
typedef Sacado::ELRCacheFad::SFad<double,global_fad_size> ELRCacheFad_SFadType;
typedef Sacado::ELRCacheFad::ViewFad<double,global_fad_size,1,ELRCacheFad_DFadType> ELRCacheFad_VFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ADConversions, ELRCacheFad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ADConversions, ELRCacheFad_SLFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ADConversions, ELRCacheFad_SFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, ViewConversions, ELRCacheFad_VFadType )

// The dx() tests in ADConversions don't make sense for these types.
// They also need more than the default constructor, and aren't designed to be
// nested.
typedef Sacado::LFad::LogicalSparse<double,bool> LFadType;
typedef Sacado::FlopCounterPack::ScalarFlopCounter<double> SFCType;
typedef Sacado::Tay::Taylor<double> TaylorType;
typedef Sacado::Tay::CacheTaylor<double> CacheTaylorType;
typedef Sacado::Rad::ADvar<double> RadType;
typedef Sacado::Rad2::ADvar<double> Rad2Type;
typedef Sacado::RadVec::ADvar<double> RadVecType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, OtherConversions, LFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, OtherConversions, SFCType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, OtherConversions, TaylorType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, OtherConversions, CacheTaylorType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, OtherConversions, RadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, OtherConversions, Rad2Type )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Conversion, OtherConversions, RadVecType )

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
