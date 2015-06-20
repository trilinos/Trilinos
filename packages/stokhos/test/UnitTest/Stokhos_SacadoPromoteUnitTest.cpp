// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

// This test requires C++11 (for static_assert), so why not use the
// standard type traits
#include <type_traits>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TestingHelpers.hpp"

#include "Sacado.hpp"
#include "Stokhos_Sacado.hpp"
#include "Stokhos_Sacado_Kokkos.hpp"
#include "Sacado_mpl_apply.hpp"

template <typename uq_type, typename expr_type>
bool testUQExprPromote() {
  using Sacado::Promote;
  using std::is_same;

  typedef typename Sacado::ValueType<uq_type>::type value_type;
  typedef typename Sacado::ScalarType<uq_type>::type scalar_type;

  static_assert(
    is_same<typename Promote<uq_type,uq_type>::type, uq_type >::value,
    "Promote<uq_type,uq_type>::type != uq_type");

  static_assert(
    is_same<typename Promote<uq_type,value_type>::type, uq_type >::value,
    "Promote<uq_type,value_type>::type != uq_type");

  static_assert(
    is_same<typename Promote<value_type,uq_type>::type, uq_type >::value,
    "Promote<value_type,uq_type>::type != uq_type");

  static_assert(
    is_same<typename Promote<uq_type,scalar_type>::type, uq_type >::value,
    "Promote<uq_type,scalar_type>::type != uq_type");

  static_assert(
    is_same<typename Promote<scalar_type,uq_type>::type, uq_type >::value,
    "Promote<scalar_type,uq_type>::type != uq_type");

  static_assert(
    is_same<typename Promote<uq_type,expr_type>::type, uq_type >::value,
    "Promote<expr_type,uq_type>::type != uq_type");

  static_assert(
    is_same<typename Promote<expr_type,uq_type>::type, uq_type >::value,
    "Promote<expr_type,uq_type>::type != uq_type");

  static_assert(
    is_same<typename Promote<scalar_type,expr_type>::type, uq_type >::value,
    "Promote<scalar_type,uq_type>::type != uq_type");

  static_assert(
    is_same<typename Promote<expr_type,scalar_type>::type, uq_type >::value,
    "Promote<expr_type,scalar_type>::type != uq_type");

   static_assert(
    is_same<typename Promote<value_type,expr_type>::type, uq_type >::value,
    "Promote<value_type,uq_type>::type != uq_type");

  static_assert(
    is_same<typename Promote<expr_type,value_type>::type, uq_type >::value,
    "Promote<expr_type,value_type>::type != uq_type");

  // These tests are all compile-time tests, so if the test compiles,
  // it passes...
  return true;
}

template <typename uq_type, typename expr1_type, typename expr2_type>
bool testUQExprPromote2() {
  using Sacado::Promote;
  using std::is_same;

  static_assert(
    is_same<typename Promote<expr1_type,expr2_type>::type, uq_type >::value,
    "Promote<expr1_type,expr2_type>::type != uq_type");

  // These tests are all compile-time tests, so if the test compiles,
  // it passes...
  return true;
}

template <typename uq_type>
bool testUQPromote() {

  // Get the type of the result of the expression 'ad_type * ad_type'
  // The use of declval gets around actually instantiation objects of type
  // ad_type.
  typedef decltype(std::declval<uq_type>()*std::declval<uq_type>()) bi_expr_type;
  bool res1 = testUQExprPromote<uq_type,bi_expr_type>();

  // Get the type of the result of the expression '-ad_type'
  // This provides additional testing for Promote specializations as unary
  // expressions are convertible to/from the scalar type
  typedef decltype(-std::declval<uq_type>()) un_expr_type;
  bool res2 = testUQExprPromote<uq_type,un_expr_type>();

  bool res3 = testUQExprPromote2<uq_type,bi_expr_type,un_expr_type>();

  return res1 && res2 && res3;
}

template <typename uq_type>
bool testPromote() {
  using Sacado::Promote;
  using std::is_same;
  typedef Sacado::Fad::DFad<uq_type> fad_uq_type;

  bool res1 = testUQPromote<uq_type>();
  bool res2 = testUQPromote<fad_uq_type>();

  typedef decltype(std::declval<uq_type>()*std::declval<uq_type>()) uq_expr_type;
  typedef decltype(-std::declval<fad_uq_type>()) fad_uq_expr_type;
  static_assert(
    is_same<typename Promote<uq_expr_type,fad_uq_expr_type>::type, fad_uq_type >::value,
    "Promote<uq_expr_type,fad_uq_expr_type>::type != fad_uq_type");

  return res1 && res2;
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Promote, Promote, UQ )
{
  success = testPromote<UQ>();
}

//
// Sacado::ETV::Vector
//
typedef Stokhos::StandardStorage<int,double> storage_type;
typedef Sacado::ETV::Vector<double,storage_type> vector_type;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Promote, vector_type )

//
// Sacado::ETV::Vector2
//
typedef Sacado::ETV::Vector2<double,storage_type> vector2_type;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Promote, vector2_type )

//
// Sacado::PCE::OrthogPoly
//
typedef Sacado::PCE::OrthogPoly<double,storage_type> orthog_poly_type;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Promote, orthog_poly_type )

//
// Sacado::ETPCE::OrthogPoly
//
typedef Sacado::ETPCE::OrthogPoly<double,storage_type> et_orthog_poly_type;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Promote, et_orthog_poly_type )

//
// Sacado::MP::Vector
//
typedef Kokkos::DefaultExecutionSpace device;
typedef Stokhos::StaticFixedStorage<int,double,32,device> static_storage_type;
typedef Sacado::MP::Vector<static_storage_type> kokkos_mp_type;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Promote, kokkos_mp_type )

//
// Sacado::UQ::PCE
//
typedef Stokhos::DynamicStorage<int,double,device> dynamic_storage_type;
typedef Sacado::UQ::PCE<dynamic_storage_type> kokkos_pce_type;
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Promote, Promote, kokkos_pce_type )

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
