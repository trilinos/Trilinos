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

#include "Stokhos.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

namespace SparseGridQuadratureUnitTest {

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<OrdinalType,ValueType> > basis;
    const OrdinalType d;
    const OrdinalType p;

    UnitTestSetup() : d(2), p(5) {

      // Create product basis
      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<OrdinalType,ValueType> > > bases(d);
      for (OrdinalType i=0; i<d; i++)
        bases[i] =
          Teuchos::rcp(new Stokhos::LegendreBasis<OrdinalType,ValueType>(
                         p, true, Stokhos::MODERATE_GROWTH));
      basis =
        Teuchos::rcp(new Stokhos::CompletePolynomialBasis<OrdinalType,ValueType>(bases));
    }

  };

  UnitTestSetup<int,double> setup;

#ifdef HAVE_STOKHOS_DAKOTA

  TEUCHOS_UNIT_TEST( Stokhos_SparseGridQuadrature, NumPoints ) {
    const Stokhos::SparseGridQuadrature<int,double> quad(
      setup.basis, setup.p, 1e-12, Pecos::MODERATE_RESTRICTED_GROWTH);
    const Teuchos::Array<double>& weights = quad.getQuadWeights();
    int nqp = weights.size();
    int nqp_gold = 181;

    if (nqp == nqp_gold)
      success = true;
    else
      success = false;

    out << std::endl
        << "Check: quad_weight.size() = " << nqp << " == " << nqp_gold
        << " : ";
    if (success) out << "Passed.";
    else
      out << "Failed!";
    out << std::endl;
  }

#endif

  TEUCHOS_UNIT_TEST( Stokhos_SmolyakSparseGridQuadrature, NumPoints ) {
    const Stokhos::TotalOrderIndexSet<int> index_set(setup.d, setup.p);
    const Stokhos::SmolyakSparseGridQuadrature<int,double> quad(
      setup.basis, index_set, 1e-12);
    const Teuchos::Array<double>& weights = quad.getQuadWeights();
    int nqp = weights.size();
    int nqp_gold = 181;

    if (nqp == nqp_gold)
      success = true;
    else
      success = false;

    out << std::endl
        << "Check: quad_weight.size() = " << nqp << " == " << nqp_gold
        << " : ";
    if (success) out << "Passed.";
    else
      out << "Failed!";
    out << std::endl;
  }

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
