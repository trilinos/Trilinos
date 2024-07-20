// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set));

    // Build corresponding pseudospectral operator
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, false);

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
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(setup.p, true, Stokhos::MODERATE_GROWTH));
    Stokhos::TotalOrderIndexSet<ordinal_type> coeff_index_set(setup.d, setup.p);
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set));

    // Build corresponding pseudospectral operator
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, false);

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
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set));

    // Build corresponding pseudospectral operator
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, true);

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
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(setup.p, true, Stokhos::MODERATE_GROWTH));
    Stokhos::TotalOrderIndexSet<ordinal_type> coeff_index_set(setup.d, setup.p);
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set));

    // Build corresponding pseudospectral operator
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, true);

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
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set));

    // Build corresponding pseudospectral operator
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, false);

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
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set));

    // Build corresponding pseudospectral operator
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, false);

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
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set));

    // Build corresponding pseudospectral operator
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, true);

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
    Teuchos::RCP< Stokhos::SmolyakBasis<ordinal_type,value_type> > smolyak_basis = Teuchos::rcp(new Stokhos::SmolyakBasis<ordinal_type,value_type>(bases, coeff_index_set));

    // Build corresponding pseudospectral operator
    Stokhos::SmolyakPseudoSpectralOperator<ordinal_type,value_type> sm_op(
      *smolyak_basis, true);

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
