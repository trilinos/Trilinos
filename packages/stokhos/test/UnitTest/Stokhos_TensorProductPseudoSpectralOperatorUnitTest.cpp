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

namespace TensorProductUnitTest {

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    ValueType rtol, atol;
    OrdinalType p,d;
    
    UnitTestSetup() {
      rtol = 1e-12;
      atol = 1e-12;
      d = 3;
      p = 5;
    }
    
  };

  typedef int ordinal_type;
  typedef double value_type;
  UnitTestSetup<ordinal_type,value_type> setup;

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

  TEUCHOS_UNIT_TEST( TensorProduct, IsotropicPoints ) {
    success = true;

    // Build tensor product basis of dimension d and order p
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(setup.p, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TensorProductBasis<ordinal_type,value_type>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis);
    Stokhos::TensorProductQuadrature<ordinal_type,value_type> quad(basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *basis, quad);

    success = Stokhos::testPseudoSpectralPoints(
      tp_op, quad_op, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( TensorProduct, AnisotropicPoints ) {
    success = true;

    // Build aniostropic tensor product basis of dimension d
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TensorProductBasis<ordinal_type,value_type>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis);
    Stokhos::TensorProductQuadrature<ordinal_type,value_type> quad(basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *basis, quad);

    success = Stokhos::testPseudoSpectralPoints(
      tp_op, quad_op, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( TensorProduct, Apply ) {
    success = true;

    // Build anisotropic tensor product basis of dimension d
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TensorProductBasis<ordinal_type,value_type>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis);
    Stokhos::TensorProductQuadrature<ordinal_type,value_type> quad(basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *basis, quad);

    success = Stokhos::testPseudoSpectralApply(
      tp_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( TensorProduct, Apply_Trans ) {
    success = true;

    // Build anisotropic tensor product basis of dimension d
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TensorProductBasis<ordinal_type,value_type>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis);
    Stokhos::TensorProductQuadrature<ordinal_type,value_type> quad(basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *basis, quad);

    success = Stokhos::testPseudoSpectralApplyTrans(
      tp_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( TensorProduct, Apply_PST_Lexo ) {
    success = true;

    // Build anisotropic tensor product basis of dimension d
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    typedef Stokhos::LexographicLess<Stokhos::MultiIndex<ordinal_type> > coeff_compare;
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TensorProductBasis<ordinal_type,value_type,coeff_compare>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis, true);
    Stokhos::TensorProductQuadrature<ordinal_type,value_type> quad(basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *basis, quad);

    success = Stokhos::testPseudoSpectralApply(
      tp_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( TensorProduct, Apply_PST_TotalOrder ) {
    success = true;

    // Build anisotropic tensor product basis of dimension d
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TensorProductBasis<ordinal_type,value_type>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis, true);
    Stokhos::TensorProductQuadrature<ordinal_type,value_type> quad(basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *basis, quad);

    success = Stokhos::testPseudoSpectralApply(
      tp_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( TensorProduct, Apply_PST_Trans_Lexo ) {
    success = true;

    // Build anisotropic tensor product basis of dimension d
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    typedef Stokhos::LexographicLess<Stokhos::MultiIndex<ordinal_type> > coeff_compare;
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TensorProductBasis<ordinal_type,value_type,coeff_compare>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis, true);
    Stokhos::TensorProductQuadrature<ordinal_type,value_type> quad(basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *basis, quad);

    success = Stokhos::testPseudoSpectralApplyTrans(
      tp_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( TensorProduct, Apply_PST_Trans_TotalOrder ) {
    success = true;

    // Build anisotropic tensor product basis of dimension d
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TensorProductBasis<ordinal_type,value_type>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis, true);
    Stokhos::TensorProductQuadrature<ordinal_type,value_type> quad(basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *basis, quad);

    success = Stokhos::testPseudoSpectralApplyTrans(
      tp_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( TensorProduct, IsotropicDiscreteOrthogonality ) {
    success = true;

    // Build tensor product basis of dimension d and order p
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(setup.p, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TensorProductBasis<ordinal_type,value_type>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis, true);

    success = Stokhos::testPseudoSpectralDiscreteOrthogonality(
      *basis, tp_op, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( TensorProduct, AnisotropicDiscreteOrthogonality ) {
    success = true;

    // Build tensor product basis of dimension d and order p
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TensorProductBasis<ordinal_type,value_type>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis, true);

    success = Stokhos::testPseudoSpectralDiscreteOrthogonality(
      *basis, tp_op, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( TotalOrder, IsotropicPoints ) {
    success = true;

    // Build total order basis of dimension d and order p
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(setup.p, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TotalOrderBasis<ordinal_type,value_type>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis);
    Stokhos::TensorProductQuadrature<ordinal_type,value_type> quad(basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *basis, quad);

    success = Stokhos::testPseudoSpectralPoints(
      tp_op, quad_op, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( TotalOrder, AnisotropicPoints ) {
    success = true;

    // Build aniostropic tensor product basis of dimension d
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TotalOrderBasis<ordinal_type,value_type>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis);
    Stokhos::TensorProductQuadrature<ordinal_type,value_type> quad(basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *basis, quad);

    success = Stokhos::testPseudoSpectralPoints(
      tp_op, quad_op, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( TotalOrder, Apply ) {
    success = true;

    // Build anisotropic tensor product basis of dimension d
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TotalOrderBasis<ordinal_type,value_type>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis);
    Stokhos::TensorProductQuadrature<ordinal_type,value_type> quad(basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *basis, quad);

    success = Stokhos::testPseudoSpectralApply(
      tp_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( TotalOrder, ApplyTrans ) {
    success = true;

    // Build anisotropic tensor product basis of dimension d
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TotalOrderBasis<ordinal_type,value_type>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis);
    Stokhos::TensorProductQuadrature<ordinal_type,value_type> quad(basis);
    Stokhos::QuadraturePseudoSpectralOperator<ordinal_type,value_type> quad_op(
      *basis, quad);

    success = Stokhos::testPseudoSpectralApplyTrans(
      tp_op, quad_op, quad_func1, quad_func2, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( TotalOrder, IsotropicDiscreteOrthogonality ) {
    success = true;

    // Build tensor product basis of dimension d and order p
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(setup.p, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TotalOrderBasis<ordinal_type,value_type>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis);

    success = Stokhos::testPseudoSpectralDiscreteOrthogonality(
      *basis, tp_op, setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( TotalOrder, AnisotropicDiscreteOrthogonality ) {
    success = true;

    // Build tensor product basis of dimension d and order p
    Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> > > bases(setup.d);
    for (ordinal_type i=0; i<setup.d; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<ordinal_type,value_type>(i+1, true));
    Teuchos::RCP<const Stokhos::ProductBasis<ordinal_type,value_type> > basis = Teuchos::rcp(new Stokhos::TotalOrderBasis<ordinal_type,value_type>(bases));
    Stokhos::TensorProductPseudoSpectralOperator<ordinal_type,value_type> tp_op(
      *basis);

    success = Stokhos::testPseudoSpectralDiscreteOrthogonality(
      *basis, tp_op, setup.rtol, setup.atol, out);
  }

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  int res = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  Teuchos::TimeMonitor::summarize(std::cout);
  return res;
}
