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
#ifdef HAVE_STOKHOS_KOKKOSCORE
#include "Stokhos_Sacado_Kokkos.hpp"
#endif
#include "Sacado_mpl_apply.hpp"

template <typename uq_type>
bool testUQPromote() {
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

  // These tests are all compile-time tests, so if the test compiles,
  // it passes...
  return true;
}

template <typename uq_type>
bool testPromote() {
  typedef Sacado::Fad::DFad<uq_type> fad_uq_type;

  testUQPromote<uq_type>();
  testUQPromote<fad_uq_type>();

  // These tests are all compile-time tests, so if the test compiles,
  // it passes...
  return true;
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

#ifdef HAVE_STOKHOS_KOKKOSCORE
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
#endif

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
