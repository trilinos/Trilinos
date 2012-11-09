// $Id$ 
// $Source$ 
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
#include "Teuchos_ArrayView.hpp"

#include "Stokhos.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

namespace SmolyakBasisUtilsUnitTest {

  // Common setup for unit tests
  template <typename ordinal_type, typename value_type>
  struct UnitTestSetup {
    value_type rtol, atol;
    value_type apply_rtol, apply_atol;
    ordinal_type p,d;
    
    UnitTestSetup() {
      rtol = 1e-12; 
      atol = 1e-12; 
      apply_rtol = 1e-2;
      apply_atol = 1e-3;
      p = 4;
      d = 2;
    }
    
  };

  typedef int ordinal_type;
  typedef double value_type;
  typedef UnitTestSetup<ordinal_type,value_type> setup_type;
  setup_type setup;

  // Function for testing quadratures
  value_type quad_func1(const Teuchos::ArrayView<const value_type>& x) {
    value_type val = 0.0;
    for (int i=0; i<x.size(); i++)
      val += x[i];
    return std::exp(val);
  }

  // Function for testing quadratures
  value_type quad_func2(const Teuchos::ArrayView<const value_type>& x) {
    value_type val = 0.0;
    for (int i=0; i<x.size(); i++)
      val += x[i];
    return std::sin(val);
  }

  TEUCHOS_UNIT_TEST( Direct, Linear ) {
    // Build isotropic Smolyak basis of dimension d and order p with 
    // linear growth 
     Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(setup.p, true));
    Stokhos::TotalOrderIndexSet<ordinal_type> coeff_index_set(setup.d, setup.p);
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set, coeff_growth));

    // Build corresponding pseudospectral operator
    Teuchos::Array< Stokhos::EvenGrowthRule<ordinal_type> > point_growth(
      setup.d);
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, point_growth, false);

     // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      smolyak_basis, setup.p, 1e-12, Pecos::SLOW_RESTRICTED_GROWTH);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *smolyak_basis, quad);

    // Test grid
    success = Stokhos::testPseudoSpectralPoints(
      sm_op, quad_op, setup.rtol, setup.atol, out);

    // Test apply
    success = success && 
      Stokhos::testPseudoSpectralApply(
	sm_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);

    // Test transpose apply
    success = success && 
      Stokhos::testPseudoSpectralApplyTrans(
	sm_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);

    // Test discrete orthogonality
    success = success && 
      Stokhos::testPseudoSpectralDiscreteOrthogonality(
	*smolyak_basis, sm_op, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Direct, ModerateLinear ) {
    // Build isotropic Smolyak basis of dimension d and order p with 
    // moderate linear growth 
     Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(setup.p, true));
    Stokhos::TotalOrderIndexSet<ordinal_type> coeff_index_set(setup.d, setup.p);
    Teuchos::Array< Stokhos::LinearGrowthRule<ordinal_type> > coeff_growth(
      setup.d, Stokhos::LinearGrowthRule<ordinal_type>(2,0));
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set, coeff_growth));

    // Build corresponding pseudospectral operator
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > point_growth(
      setup.d);
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, point_growth, false);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      smolyak_basis, setup.p, 1e-12, Pecos::MODERATE_RESTRICTED_GROWTH);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *smolyak_basis, quad);
    
    // Test grid
    success = Stokhos::testPseudoSpectralPoints(
      sm_op, quad_op, setup.rtol, setup.atol, out);

    // Test apply
    success = success && 
      Stokhos::testPseudoSpectralApply(
	sm_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);

    // Test transpose apply
    success = success && 
      Stokhos::testPseudoSpectralApplyTrans(
	sm_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);

    // Direct apply will not satisfy discrete orthogonality
  }

  TEUCHOS_UNIT_TEST( Smolyak, Linear ) {
    // Build isotropic Smolyak basis of dimension d and order p with 
    // linear growth 
     Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(setup.p, true));
    Stokhos::TotalOrderIndexSet<ordinal_type> coeff_index_set(setup.d, setup.p);
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > coeff_growth(
      setup.d);
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set, coeff_growth));

    // Build corresponding pseudospectral operator
    Teuchos::Array< Stokhos::EvenGrowthRule<ordinal_type> > point_growth(
      setup.d);
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, point_growth, true);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      smolyak_basis, setup.p, 1e-12, Pecos::SLOW_RESTRICTED_GROWTH);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *smolyak_basis, quad);

    Stokhos::TensorProductQuadrature<ordinal_type,value_type> tp_quad(
      smolyak_basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> tp_quad_op(*smolyak_basis, tp_quad);

    // Test grid
    success = Stokhos::testPseudoSpectralPoints(
      sm_op, quad_op, setup.rtol, setup.atol, out);

    // Test apply
    success = success && 
      Stokhos::testPseudoSpectralApply(
	sm_op, tp_quad_op, quad_func1, quad_func2, setup.apply_rtol, setup.apply_atol, out);

    // Test transpose apply
    success = success && 
      Stokhos::testPseudoSpectralApplyTrans(
	sm_op, tp_quad_op, quad_func1, quad_func2, setup.apply_rtol, setup.apply_atol, out);

    // Test discrete orthogonality
    success = success && 
      Stokhos::testPseudoSpectralDiscreteOrthogonality(
	*smolyak_basis, sm_op, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Smolyak, ModerateLinear ) {
    // Build isotropic Smolyak basis of dimension d and order p with 
    // moderate linear growth 
     Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(setup.p, true));
    Stokhos::TotalOrderIndexSet<ordinal_type> coeff_index_set(setup.d, setup.p);
    Teuchos::Array< Stokhos::LinearGrowthRule<ordinal_type> > coeff_growth(
      setup.d, Stokhos::LinearGrowthRule<ordinal_type>(2,0));
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set, coeff_growth));

    // Build corresponding pseudospectral operator
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > point_growth(
      setup.d);
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, point_growth, true);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      smolyak_basis, setup.p, 1e-12, Pecos::MODERATE_RESTRICTED_GROWTH);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *smolyak_basis, quad);

    Stokhos::TensorProductQuadrature<ordinal_type,value_type> tp_quad(
      smolyak_basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> tp_quad_op(*smolyak_basis, tp_quad);
    
    // Test grid
    success = Stokhos::testPseudoSpectralPoints(
      sm_op, quad_op, setup.rtol, setup.atol, out);

    // Test apply
    success = success && 
      Stokhos::testPseudoSpectralApply(
	sm_op, tp_quad_op, quad_func1, quad_func2, setup.apply_rtol, setup.apply_atol, out);

    // Test transpose apply
    success = success && 
      Stokhos::testPseudoSpectralApplyTrans(
	sm_op, tp_quad_op, quad_func1, quad_func2, setup.apply_rtol, setup.apply_atol, out);

    // Test discrete orthogonality
    success = success && 
      Stokhos::testPseudoSpectralDiscreteOrthogonality(
	*smolyak_basis, sm_op, setup.rtol, setup.atol, out);
  }

#ifdef HAVE_STOKHOS_DAKOTA 

  TEUCHOS_UNIT_TEST( Direct, ClenshawCurtis ) {
    // Build isotropic Smolyak basis of dimension d and order p with 
    // exponential growth 
     Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::ClenshawCurtisLegendreBasis<ordinal_type,value_type>(setup.p, true));
    Stokhos::TotalOrderIndexSet<ordinal_type> coeff_index_set(setup.d, setup.p);
    Teuchos::Array< Stokhos::ClenshawCurtisExponentialGrowthRule<ordinal_type> > coeff_growth(setup.d);
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set, coeff_growth));

    // Build corresponding pseudospectral operator
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > point_growth(
      setup.d);
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, point_growth, false);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      smolyak_basis, setup.p, 1e-12, Pecos::UNRESTRICTED_GROWTH);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *smolyak_basis, quad);
    
    // Test grid
    success = Stokhos::testPseudoSpectralPoints(
      sm_op, quad_op, setup.rtol, setup.atol, out);

    // Test apply
    success = success && 
      Stokhos::testPseudoSpectralApply(
	sm_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);

    // Test transpose apply
    success = success && 
      Stokhos::testPseudoSpectralApplyTrans(
	sm_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);

    // Direct apply will not satisfy discrete orthogonality
  }

  TEUCHOS_UNIT_TEST( Direct, GaussPatterson ) {
    // Build isotropic Smolyak basis of dimension d and order p with 
    // exponential Gauss-Patterson growth 
     Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::GaussPattersonLegendreBasis<ordinal_type,value_type>(setup.p, true));
    Stokhos::TotalOrderIndexSet<ordinal_type> coeff_index_set(setup.d, setup.p);
    Teuchos::Array< Stokhos::GaussPattersonExponentialGrowthRule<ordinal_type> > coeff_growth(setup.d);
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set, coeff_growth));

    // Build corresponding pseudospectral operator
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > point_growth(
      setup.d);
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, point_growth, false);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      smolyak_basis, setup.p, 1e-12, Pecos::UNRESTRICTED_GROWTH);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *smolyak_basis, quad);
    
    // Test grid
    success = Stokhos::testPseudoSpectralPoints(
      sm_op, quad_op, setup.rtol, setup.atol, out);

    // Test apply
    success = success && 
      Stokhos::testPseudoSpectralApply(
	sm_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);

    // Test transpose apply
    success = success && 
      Stokhos::testPseudoSpectralApplyTrans(
	sm_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);

    // Direct apply will not satisfy discrete orthogonality
  }

  TEUCHOS_UNIT_TEST( Smolyak, ClenshawCurtis ) {
    // Build isotropic Smolyak basis of dimension d and order p with 
    // exponential growth 
     Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::ClenshawCurtisLegendreBasis<ordinal_type,value_type>(setup.p, true));
    Stokhos::TotalOrderIndexSet<ordinal_type> coeff_index_set(setup.d, setup.p);
    Teuchos::Array< Stokhos::ClenshawCurtisExponentialGrowthRule<ordinal_type> > coeff_growth(setup.d);
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set, coeff_growth));

    // Build corresponding pseudospectral operator
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > point_growth(
      setup.d);
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, point_growth, true);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      smolyak_basis, setup.p, 1e-12, Pecos::UNRESTRICTED_GROWTH);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *smolyak_basis, quad);

    Stokhos::TensorProductQuadrature<ordinal_type,value_type> tp_quad(
      smolyak_basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> tp_quad_op(*smolyak_basis, tp_quad);
    
    // Test grid
    success = Stokhos::testPseudoSpectralPoints(
      sm_op, quad_op, setup.rtol, setup.atol, out);

    // Test apply
    success = success && 
      Stokhos::testPseudoSpectralApply(
	sm_op, tp_quad_op, quad_func1, quad_func2, setup.apply_rtol, setup.apply_atol, out);

    // Test transpose apply
    success = success && 
      Stokhos::testPseudoSpectralApplyTrans(
	sm_op, tp_quad_op, quad_func1, quad_func2, setup.apply_rtol, setup.apply_atol, out);

    // Test discrete orthogonality
    success = success && 
      Stokhos::testPseudoSpectralDiscreteOrthogonality(
	*smolyak_basis, sm_op, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Smolyak, GaussPatterson ) {
    // Build isotropic Smolyak basis of dimension d and order p with 
    // exponential Gauss-Patterson growth 
     Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::GaussPattersonLegendreBasis<ordinal_type,value_type>(setup.p, true));
    Stokhos::TotalOrderIndexSet<ordinal_type> coeff_index_set(setup.d, setup.p);
    Teuchos::Array< Stokhos::GaussPattersonExponentialGrowthRule<ordinal_type> > coeff_growth(setup.d);
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set, coeff_growth));

    // Build corresponding pseudospectral operator
    Teuchos::Array< Stokhos::IdentityGrowthRule<ordinal_type> > point_growth(
      setup.d);
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, point_growth, true);

    // Generate sparse grids using original approach
    Stokhos::SparseGridQuadrature<ordinal_type,value_type> quad(
      smolyak_basis, setup.p, 1e-12, Pecos::UNRESTRICTED_GROWTH);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *smolyak_basis, quad);

    Stokhos::TensorProductQuadrature<ordinal_type,value_type> tp_quad(
      smolyak_basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> tp_quad_op(*smolyak_basis, tp_quad);
    
    // Test grid
    success = Stokhos::testPseudoSpectralPoints(
      sm_op, quad_op, setup.rtol, setup.atol, out);

    // Test apply
    success = success && 
      Stokhos::testPseudoSpectralApply(
	sm_op, tp_quad_op, quad_func1, quad_func2, setup.apply_rtol, setup.apply_atol, out);

    // Test transpose apply
    success = success && 
      Stokhos::testPseudoSpectralApplyTrans(
	sm_op, tp_quad_op, quad_func1, quad_func2, setup.apply_rtol, setup.apply_atol, out);

    // Test discrete orthogonality
    success = success && 
      Stokhos::testPseudoSpectralDiscreteOrthogonality(
	*smolyak_basis, sm_op, setup.rtol, setup.atol, out);
  }

#endif

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
