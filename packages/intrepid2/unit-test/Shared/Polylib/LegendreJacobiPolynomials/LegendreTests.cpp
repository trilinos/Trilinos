// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov),
//                    Mauro Perego  (mperego@sandia.gov), or
//                    Nate Roberts  (nvrober@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   LegendreTests.cpp
    \brief  Tests to verify various implementations of Legendre polynomials.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_Polylib.hpp"     // provides JacobiPolynomial implementation used for verification
#include "Intrepid2_Polynomials.hpp" // provides the Legendre polynomials under test
#include "Intrepid2_Types.hpp"

#include "Intrepid2_TestUtils.hpp"

#include "Kokkos_Core.hpp"

namespace
{
  using namespace Intrepid2;
  
  void testLegendreDerivatives(const double x, const int maxPolyOrder, const int derivativeOrder, const double tol,
                               Teuchos::FancyOStream &out, bool &success)
  {
    Kokkos::View<double*> legendreDerivatives("legendre derivatives",maxPolyOrder+1);
    
    using Intrepid2::Polynomials::legendreDerivativeValues;
    // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int dummy_index)
    {
      legendreDerivativeValues(legendreDerivatives, maxPolyOrder, x, derivativeOrder);
    });
    
    auto legendreDerivativesHost = getHostCopy(legendreDerivatives);
    
    // Polylib Jacobi computes a single function at potentially multiple points
    const int numPoints = 1;
    Kokkos::View<double*> inputPoints("x point", numPoints);
    Kokkos::deep_copy(inputPoints,x);
    
    Kokkos::View<double*> jacobiValue("jacobi value from Polynomials::JacobiPolynomial",numPoints);
    for (int i=0; i<=maxPolyOrder; i++)
    {
      double expectedValue;
      if (i <= derivativeOrder)
      {
        expectedValue = 0.0;
      }
      else
      {
        const Kokkos::View<double*> null;
        const double alpha = 0;
        const double beta = 0;
        double scaleFactor = 1.0;
        for (int j=1; j<=derivativeOrder; j++)
        {
          scaleFactor *= 0.5 * (j+alpha+beta+i);
        }
        
        // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
        Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int dummy_index)
        {
          Intrepid2::Polylib::Serial::JacobiPolynomial(numPoints, inputPoints, jacobiValue, null, i-derivativeOrder, alpha+derivativeOrder, beta+derivativeOrder);
        });
        
        auto jacobiValueHost = getHostCopy(jacobiValue);
        
        expectedValue = scaleFactor * jacobiValueHost(0);
        bool valuesAreBothSmall = valuesAreSmall(expectedValue, legendreDerivativesHost(i), tol);
        if (! valuesAreBothSmall)
        {
          if (! approximatelyEqual(expectedValue, legendreDerivativesHost(i), tol) )
          {
            double diff = abs(expectedValue - legendreDerivativesHost(i));
            success = false;
            out << "FAILURE: testLegendreDerivative -- expectedValue was " << expectedValue;
            out << "; actual was " << legendreDerivativesHost(i) << " (diff: " << diff << ") for derivative " << derivativeOrder << " of i = " << i;
            out << " at x = " << x << std::endl;
          }
        }
      }
    }
  }
  
  TEUCHOS_UNIT_TEST( Legendre, DerivativesMatchJacobi )
  {
    // test confirms that legendreDerivativeValues() and JacobiPolynomial() values agree
    const double tol = TEST_TOLERANCE_TIGHT;
    
    std::vector<double> x_values = {{-1.0,-0.5,-1.0/3.0,0.0,1.0/3.0,0.50,1.0}};
    int polyOrderMax = 10;
    
    for (int derivativeOrder=1; derivativeOrder<=10; derivativeOrder++)
    {
      for (auto x : x_values)
      {
        testLegendreDerivatives(x, polyOrderMax, derivativeOrder, tol, out, success);
      }
    }
  }
} // namespace
