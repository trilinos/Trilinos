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

#include "Sacado_No_Kokkos.hpp"
#include "Sacado_Fad_SimpleFad.hpp"
#include "Sacado_Tay_CacheTaylor.hpp"
#include "Sacado_mpl_apply.hpp"

template <typename ad_type>
bool testADPromote() {
  using Sacado::Promote;
  using std::is_same;

  typedef typename Sacado::ValueType<ad_type>::type value_type;
  typedef typename Sacado::ScalarType<ad_type>::type scalar_type;

  // Get the type of the result of the expression '- ad_type'
  // The use of declval gets around actually instantiation objects of type
  // ad_type.
  // We use a unary expression to catch special-case problems with Promote
  // since the AD type may be convertible to the expression type
  typedef decltype(-std::declval<ad_type>()) expr_type;

  static_assert(
    is_same<typename Promote<ad_type,ad_type>::type, ad_type >::value,
    "Promote<ad_type,ad_type>::type != ad_type");

  static_assert(
    is_same<typename Promote<ad_type,value_type>::type, ad_type >::value,
    "Promote<ad_type,value_type>::type != ad_type");

  static_assert(
    is_same<typename Promote<value_type,ad_type>::type, ad_type >::value,
    "Promote<value_type,ad_type>::type != ad_type");

  static_assert(
    is_same<typename Promote<ad_type,scalar_type>::type, ad_type >::value,
    "Promote<ad_type,scalar_type>::type != ad_type");

  static_assert(
    is_same<typename Promote<scalar_type,ad_type>::type, ad_type >::value,
    "Promote<scalar_type,ad_type>::type != ad_type");

  static_assert(
    is_same<typename Promote<ad_type,expr_type>::type, ad_type >::value,
    "Promote<ad_type,expr_type>::type != ad_type");

  static_assert(
    is_same<typename Promote<expr_type,ad_type>::type, ad_type >::value,
    "Promote<expr_type,ad_type>::type != ad_type");

  static_assert(
    is_same<typename Promote<expr_type,value_type>::type, ad_type >::value,
    "Promote<expr_type,value_type>::type != ad_type");

  static_assert(
    is_same<typename Promote<value_type,expr_type>::type, ad_type >::value,
    "Promote<value_type,expr_type>::type != ad_type");

  static_assert(
    is_same<typename Promote<expr_type,scalar_type>::type, ad_type >::value,
    "Promote<expr_type,scalar_type>::type != ad_type");

  static_assert(
    is_same<typename Promote<scalar_type,expr_type>::type, ad_type >::value,
    "Promote<scalar_type,expr_type>::type != ad_type");

  // These tests are all compile-time tests, so if the test compiles,
  // it passes...
  return true;
}

template <typename view_type>
bool testViewPromote() {
  using Sacado::Promote;
  using std::is_same;

  typedef typename Sacado::ValueType<view_type>::type value_type;
  typedef typename Sacado::ScalarType<view_type>::type scalar_type;
  typedef typename view_type::base_fad_type base_fad_type;

  // Get the type of the result of the expression '- view_type'
  // The use of declval gets around actually instantiation objects of type
  // view_type.
  // We use a unary expression to catch special-case problems with Promote
  // since the AD type may be convertible to the expression type
  typedef decltype(-std::declval<view_type>()) expr_type;

  static_assert(
    is_same<typename Promote<view_type,view_type>::type, base_fad_type >::value,
    "Promote<view_type,view_type>::type != base_fad_type");

  static_assert(
    is_same<typename Promote<view_type,value_type>::type, base_fad_type >::value,
    "Promote<view_type,value_type>::type != base_fad_type");

  static_assert(
    is_same<typename Promote<value_type,view_type>::type, base_fad_type >::value,
    "Promote<value_type,view_type>::type != base_fad_type");

  static_assert(
    is_same<typename Promote<view_type,scalar_type>::type, base_fad_type >::value,
    "Promote<view_type,scalar_type>::type != base_fad_type");

  static_assert(
    is_same<typename Promote<scalar_type,view_type>::type, base_fad_type >::value,
    "Promote<scalar_type,view_type>::type != base_fad_type");

  static_assert(
    is_same<typename Promote<view_type,expr_type>::type, base_fad_type >::value,
    "Promote<view_type,expr_type>::type != base_fad_type");

  static_assert(
    is_same<typename Promote<expr_type,view_type>::type, base_fad_type >::value,
    "Promote<expr_type,view_type>::type != base_fad_type");

  static_assert(
    is_same<typename Promote<expr_type,value_type>::type, base_fad_type >::value,
    "Promote<expr_type,value_type>::type != base_fad_type");

  static_assert(
    is_same<typename Promote<value_type,expr_type>::type, base_fad_type >::value,
    "Promote<value_type,expr_type>::type != base_fad_type");

  static_assert(
    is_same<typename Promote<expr_type,scalar_type>::type, base_fad_type >::value,
    "Promote<expr_type,scalar_type>::type != base_fad_type");

  static_assert(
    is_same<typename Promote<scalar_type,expr_type>::type, base_fad_type >::value,
    "Promote<scalar_type,expr_type>::type != base_fad_type");

  // These tests are all compile-time tests, so if the test compiles,
  // it passes...
  return true;
}

template <typename fad_type>
bool testFadPromote() {
  using Sacado::Promote;
  using std::is_same;

  typedef typename Sacado::ViewFadType<fad_type,0,1>::type view_fad_type;

  typedef typename Sacado::mpl::apply< fad_type, fad_type >::type fad_fad_type;
  typedef typename Sacado::mpl::apply< view_fad_type, fad_type >::type view_fad_fad_type;
  typedef typename Sacado::mpl::apply< view_fad_type, view_fad_type >::type view_view_fad_type;

  typedef typename view_fad_type::base_fad_type base_fad_type;
  typedef typename view_fad_fad_type::base_fad_type base_fad_fad_type;

  testADPromote<fad_type>();
  testADPromote<fad_fad_type>();
  testViewPromote<view_fad_type>();
  testViewPromote<view_fad_fad_type>();
  testViewPromote<view_view_fad_type>();

  static_assert(
    is_same<typename Promote<view_fad_type,fad_type>::type, fad_type >::value,
    "Promote<view_fad_type,fad_type>::type != fad_type");

  static_assert(
    is_same<typename Promote<fad_type,view_fad_type>::type, fad_type >::value,
    "Promote<fad_type,view_fad_type>::type != fad_type");

  static_assert(
    is_same<typename Promote<view_fad_fad_type,fad_fad_type>::type, fad_fad_type >::value,
    "Promote<view_fad_fad_type,fad_fad_type>::type != fad_fad_type");

  static_assert(
    is_same<typename Promote<fad_fad_type,view_fad_fad_type>::type, fad_fad_type >::value,
    "Promote<fad_fad_type,view_fad_fad_type>::type != fad_fad_type");

  typedef decltype(std::declval<fad_type>()*std::declval<fad_type>()) fad_expr_type;
  typedef decltype(std::declval<view_fad_type>()*std::declval<view_fad_type>()) view_fad_expr_type;

  static_assert(
    is_same<typename Promote<view_fad_type,fad_expr_type>::type, base_fad_type >::value,
    "Promote<view_fad_type,fad_expr_type>::type != base_fad_type");

  static_assert(
    is_same<typename Promote<fad_expr_type,view_fad_type>::type, base_fad_type >::value,
    "Promote<fad_expr_type,view_fad_type>::type != base_fad_type");

  static_assert(
    is_same<typename Promote<fad_type,view_fad_expr_type>::type, fad_type >::value,
    "Promote<fad_type,view_fad_expr_type>::type != fad_type");

  static_assert(
    is_same<typename Promote<view_fad_expr_type,fad_type>::type, fad_type >::value,
    "Promote<view_fad_expr_type,fad_type>::type != fad_type");

  typedef decltype(-std::declval<fad_fad_type>()) fad_fad_expr_type;
  typedef decltype(-std::declval<view_fad_fad_type>()) view_fad_fad_expr_type;
   typedef decltype(-std::declval<view_view_fad_type>()) view_view_fad_expr_type;

  static_assert(
    is_same<typename Promote<view_fad_type,fad_fad_expr_type>::type, fad_fad_type >::value,
    "Promote<view_fad_type,fad_fad_expr_type>::type != fad_fad_type");

  static_assert(
    is_same<typename Promote<fad_fad_expr_type,view_fad_type>::type, fad_fad_type >::value,
    "Promote<fad_fad_expr_type,view_fad_type>::type != fad_fad_type");

  static_assert(
    is_same<typename Promote<view_fad_type,view_fad_fad_expr_type>::type, base_fad_fad_type >::value,
    "Promote<view_fad_type,view_fad_fad_expr_type>::type != base_fad_fad_type");

  static_assert(
    is_same<typename Promote<view_fad_fad_expr_type,view_fad_type>::type, base_fad_fad_type >::value,
    "Promote<view_fad_fad_expr_type,view_fad_type>::type != base_fad_fad_type");

  static_assert(
    is_same<typename Promote<fad_type,view_view_fad_expr_type>::type, base_fad_fad_type >::value,
    "Promote<fad_type,view_fad_fad_expr_type>::type != base_fad_fad_type");

  static_assert(
    is_same<typename Promote<view_view_fad_expr_type,fad_type>::type, base_fad_fad_type >::value,
    "Promote<view_fad_fad_expr_type,fad_type>::type != base_fad_fad_type");

  static_assert(
    is_same<typename Promote<fad_expr_type,fad_fad_expr_type>::type, fad_fad_type >::value,
    "Promote<fad_expr_type,fad_fad_expr_type>::type != fad_fad_type");

   static_assert(
    is_same<typename Promote<fad_fad_expr_type,fad_expr_type>::type, fad_fad_type >::value,
    "Promote<fad_fad_expr_type,fad_expr_type>::type != fad_fad_type");

   static_assert(
    is_same<typename Promote<view_fad_expr_type,fad_fad_expr_type>::type, fad_fad_type >::value,
    "Promote<view_fad_expr_type,fad_fad_expr_type>::type != fad_fad_type");

   static_assert(
    is_same<typename Promote<fad_fad_expr_type,view_fad_expr_type>::type, fad_fad_type >::value,
    "Promote<fad_fad_expr_type,view_fad_expr_type>::type != fad_fad_type");

   static_assert(
    is_same<typename Promote<fad_expr_type,view_fad_fad_expr_type>::type, base_fad_fad_type >::value,
    "Promote<fad_expr_type,view_fad_fad_expr_type>::type != base_fad_fad_type");

   static_assert(
    is_same<typename Promote<view_fad_fad_expr_type,fad_expr_type>::type, base_fad_fad_type >::value,
    "Promote<view_fad_fad_expr_type,fad_expr_type>::type != base_fad_fad_type");

   static_assert(
    is_same<typename Promote<view_fad_expr_type,view_fad_fad_expr_type>::type, base_fad_fad_type >::value,
    "Promote<view_fad_expr_type,view_fad_fad_expr_type>::type != base_fad_fad_type");

   static_assert(
    is_same<typename Promote<view_fad_fad_expr_type,view_fad_expr_type>::type, base_fad_fad_type >::value,
    "Promote<view_fad_fad_expr_type,view_fad_expr_type>::type != base_fad_fad_type");

   static_assert(
    is_same<typename Promote<fad_expr_type,view_view_fad_expr_type>::type, base_fad_fad_type >::value,
    "Promote<fad_expr_type,view_view_fad_expr_type>::type != base_fad_fad_type");

   static_assert(
    is_same<typename Promote<view_view_fad_expr_type,fad_expr_type>::type, base_fad_fad_type >::value,
    "Promote<view_view_fad_expr_type,fad_expr_type>::type != base_fad_fad_type");

   static_assert(
    is_same<typename Promote<view_fad_expr_type,view_view_fad_expr_type>::type, base_fad_fad_type >::value,
    "Promote<view_fad_expr_type,view_view_fad_expr_type>::type != base_fad_fad_type");

   static_assert(
    is_same<typename Promote<view_view_fad_expr_type,view_fad_expr_type>::type, base_fad_fad_type >::value,
    "Promote<view_view_fad_expr_type,view_fad_expr_type>::type != base_fad_fad_type");

  // These tests are all compile-time tests, so if the test compiles,
  // it passes...
  return true;
}

template <typename scalar_type>
bool testPromote() {
  typedef Sacado::Fad::DFad<scalar_type> fad_scalar_type;

  testADPromote<scalar_type>();
  testADPromote<fad_scalar_type>();

  // These tests are all compile-time tests, so if the test compiles,
  // it passes...
  return true;
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Promote, Fad, FAD )
{
  success = testFadPromote<FAD>();
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Promote, Other, AD )
{
  success = testPromote<AD>();
}

const int global_fad_size = 10;

typedef Sacado::Fad::DFad<double> Fad_DFadType;
typedef Sacado::Fad::SLFad<double,2*global_fad_size> Fad_SLFadType;
typedef Sacado::Fad::SFad<double,global_fad_size> Fad_SFadType;
typedef Sacado::Fad::DVFad<double> Fad_DVFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Fad, Fad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Fad, Fad_SFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Fad, Fad_SLFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Fad, Fad_DVFadType )

typedef Sacado::ELRFad::DFad<double> ELRFad_DFadType;
typedef Sacado::ELRFad::SLFad<double,2*global_fad_size> ELRFad_SLFadType;
typedef Sacado::ELRFad::SFad<double,global_fad_size> ELRFad_SFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Fad, ELRFad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Fad, ELRFad_SFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Fad, ELRFad_SLFadType )

typedef Sacado::CacheFad::DFad<double> CacheFad_DFadType;
typedef Sacado::CacheFad::SLFad<double,2*global_fad_size> CacheFad_SLFadType;
typedef Sacado::CacheFad::SFad<double,global_fad_size> CacheFad_SFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Fad, CacheFad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Fad, CacheFad_SFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Fad, CacheFad_SLFadType )

typedef Sacado::ELRCacheFad::DFad<double> ELRCacheFad_DFadType;
typedef Sacado::ELRCacheFad::SLFad<double,2*global_fad_size> ELRCacheFad_SLFadType;
typedef Sacado::ELRCacheFad::SFad<double,global_fad_size> ELRCacheFad_SFadType;
typedef Sacado::ELRCacheFad::ViewFad<double,global_fad_size,1,ELRCacheFad_DFadType> ELRCacheFad_VFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Fad, ELRCacheFad_DFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Fad, ELRCacheFad_SFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Fad, ELRCacheFad_SLFadType )

typedef Sacado::Fad::SimpleFad<double> SimpleFadType;
typedef Sacado::LFad::LogicalSparse<double,bool> LFadType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Other, SimpleFadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Other, LFadType )

typedef Sacado::FlopCounterPack::ScalarFlopCounter<double> SFCType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Other, SFCType )

typedef Sacado::Tay::Taylor<double> TaylorType;
typedef Sacado::Tay::CacheTaylor<double> CacheTaylorType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Other, TaylorType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Other, CacheTaylorType )

typedef Sacado::Rad::ADvar<double> RadType;
typedef Sacado::Rad2::ADvar<double> Rad2Type;
typedef Sacado::RadVec::ADvar<double> RadVecType;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Other, RadType )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Other, Rad2Type )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Other, RadVecType )

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
