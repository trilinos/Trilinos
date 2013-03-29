// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Stokhos_UnitTestHelpers.hpp"

#include "Stokhos_Sacado.hpp"
#ifdef HAVE_STOKHOS_KOKKOSARRAY
#include "Stokhos_Sacado_Kokkos.hpp"
#endif

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Traits, ScalarType, ad_type, scalar_type )
{
  const bool is_same =
    Sacado::mpl::is_same< typename Sacado::ScalarType<ad_type>::type, scalar_type >::value;

  TEUCHOS_TEST_EQUALITY(is_same, true, out, success);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Traits, ValueType, ad_type, value_type )
{
  const bool is_same =
    Sacado::mpl::is_same< typename Sacado::ValueType<ad_type>::type, value_type >::value;
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

//
// Sacado::ETV::Vector
//
typedef Stokhos::StandardStorage<int,double> storage_type;
typedef Sacado::ETV::Vector<double,storage_type> vector_type;
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ScalarType, vector_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ValueType, vector_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Traits, IsADType, vector_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, Value, vector_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ScalarValue, vector_type, double )

typedef Sacado::ETV::Expr< Sacado::ETV::VectorImpl<double,storage_type> > vector_expr_type;
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ScalarType, vector_expr_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ValueType, vector_expr_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Traits, IsADType, vector_expr_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, Value, vector_expr_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ScalarValue, vector_expr_type, double )

//
// Sacado::ETV::Vector2 -- no need to test Expr as it is the same as Vector
//
typedef Sacado::ETV::Vector2<double,storage_type> vector2_type;
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ScalarType, vector2_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ValueType, vector2_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Traits, IsADType, vector2_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, Value, vector2_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ScalarValue, vector2_type, double )

// 
// Sacado::ETPCE::OrthogPoly
//
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

//
// Sacado::MP::Vector
//
#ifdef HAVE_STOKHOS_KOKKOSARRAY
typedef KokkosArray::Host node_type;
typedef Stokhos::DynamicStorage<int,double,node_type> kokkos_storage_type;
typedef Sacado::MP::Vector<kokkos_storage_type,node_type> mp_type;
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ScalarType, mp_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ValueType, mp_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( Traits, IsADType, mp_type )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, Value, mp_type, double )
TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( Traits, ScalarValue, mp_type, double )

// Can't test Expr because there is no way to actually create one
#endif

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
