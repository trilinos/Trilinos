// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
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
