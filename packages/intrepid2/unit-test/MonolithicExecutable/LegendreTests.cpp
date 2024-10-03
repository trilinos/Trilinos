// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  
  using DeviceType = DefaultTestDeviceType;
  using ExecutionSpace = typename DeviceType::execution_space;

  void testLegendreDerivatives(const double x, const int maxPolyOrder, const int derivativeOrder, const double tol,
                               Teuchos::FancyOStream &out, bool &success)
  {
    Kokkos::View<double*,DeviceType> legendreDerivatives("legendre derivatives",maxPolyOrder+1);
    
    using Intrepid2::Polynomials::legendreDerivativeValues;
    // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
    auto policy = Kokkos::RangePolicy<>(ExecutionSpace(),0,1);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int dummy_index)
    {
      legendreDerivativeValues(legendreDerivatives, maxPolyOrder, x, derivativeOrder);
    });
    
    auto legendreDerivativesHost = getHostCopy(legendreDerivatives);
    
    // Polylib Jacobi computes a single function at potentially multiple points
    const int numPoints = 1;
    Kokkos::View<double*,DeviceType> inputPoints("x point", numPoints);
    Kokkos::deep_copy(inputPoints,x);
    
    Kokkos::View<double*,DeviceType> jacobiValue("jacobi value from Polynomials::JacobiPolynomial",numPoints);
    for (int i=0; i<=maxPolyOrder; i++)
    {
      double expectedValue;
      if (i <= derivativeOrder)
      {
        expectedValue = 0.0;
      }
      else
      {
        const Kokkos::View<double*,DeviceType> null;
        const double alpha = 0;
        const double beta = 0;
        double scaleFactor = 1.0;
        for (int j=1; j<=derivativeOrder; j++)
        {
          scaleFactor *= 0.5 * (j+alpha+beta+i);
        }
        
        // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
        auto policy = Kokkos::RangePolicy<>(ExecutionSpace(),0,1);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int dummy_index)
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
