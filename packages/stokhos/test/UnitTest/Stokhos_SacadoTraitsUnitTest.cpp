// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <type_traits>

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos_UnitTestHelpers.hpp"

#include "Stokhos_Sacado.hpp"
#include "Stokhos_Sacado_Kokkos.hpp"

#include <Kokkos_Core.hpp>

#if defined( KOKKOS_ENABLE_OPENMP )
typedef Kokkos::OpenMP node_type;
#elif defined( KOKKOS_ENABLE_THREADS )
typedef Kokkos::Threads node_type;
#else
typedef Kokkos::Serial node_type;
#endif

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Traits, ScalarType, ad_type, scalar_type )
{
  const bool is_same =
    std::is_same< typename Sacado::ScalarType<ad_type>::type, scalar_type >::value;

  TEUCHOS_TEST_EQUALITY(is_same, true, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Traits, ValueType, ad_type, value_type )
{
  const bool is_same =
    std::is_same< typename Sacado::ValueType<ad_type>::type, value_type >::value;
  TEUCHOS_TEST_EQUALITY(is_same, true, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Traits, IsADType, ad_type )
{
  const bool is_ad = Sacado::IsADType<ad_type>::value;
  TEUCHOS_TEST_EQUALITY(is_ad, true, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( Traits, IsScalarType, ad_type )
{
  const bool is_scalar = Sacado::IsScalarType<ad_type>::value;
  TEUCHOS_TEST_EQUALITY(is_scalar, false, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Traits, Value, ad_type, value_type )
{
  value_type v(1.0);
  ad_type a(v);
  TEUCHOS_TEST_EQUALITY(Sacado::Value<ad_type>::eval(a), v, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Traits, ScalarValue, ad_type, scalar_type )
{
  scalar_type v(1.0);
  ad_type a(v);
  TEUCHOS_TEST_EQUALITY(Sacado::ScalarValue<ad_type>::eval(a), v, out, success);
}

#ifdef HAVE_STOKHOS_PCE_SCALAR_TYPE
//
// Sacado::ETPCE::OrthogPoly
//
typedef Stokhos::StandardStorage<int,double> storage_type;
typedef Sacado::ETPCE::OrthogPoly<double,storage_type> pce_type;
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ScalarType, pce_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ValueType, pce_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Traits, IsADType, pce_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, Value, pce_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ScalarValue, pce_type, double )

typedef Sacado::ETPCE::Expr< Sacado::ETPCE::OrthogPolyImpl<double,storage_type> > pce_expr_type;
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ScalarType, pce_expr_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ValueType, pce_expr_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Traits, IsADType, pce_expr_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, Value, pce_expr_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ScalarValue, pce_expr_type, double )
#endif

#ifdef HAVE_STOKHOS_ENSEMBLE_SCALAR_TYPE
//
// Sacado::MP::Vector
//
typedef Stokhos::DynamicStorage<int,double,node_type> kokkos_storage_type;
typedef Sacado::MP::Vector<kokkos_storage_type> mp_type;
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ScalarType, mp_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ValueType, mp_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Traits, IsADType, mp_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, Value, mp_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ScalarValue, mp_type, double )
#endif

// Can't test Expr because there is no way to actually create one

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
