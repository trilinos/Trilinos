// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   JacobiTests.cpp
    \brief  Tests to verify that Jacobi implementations in Intrepid2_Polynomials.hpp agree with those in Intrepid2_Polylib.hpp.
    \author Created by N.V. Roberts.
 */
 
#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_Polylib.hpp"     // provides JacobiPolynomial implementation used for verification
#include "Intrepid2_Polynomials.hpp" // provides the Legendre polynomials under test
#include "Intrepid2_TestUtils.hpp"
#include "Intrepid2_Types.hpp"

#include "Kokkos_Core.hpp"

namespace
{
  using namespace Intrepid2;
  std::vector<double> alpha_values = {{0.0, 0.2, 0.4, 0.6, 0.8, 1.0}};
  std::vector<double> x_values = {{-1.0,-0.5,-1.0/3.0,0.0,1.0/3.0,0.50,1.0}};
  std::vector<double> t_values = {{0.0, 0.2, 0.4, 0.6, 0.8, 1.0}};

  using DeviceType = DefaultTestDeviceType;
  using ExecutionSpace = typename DeviceType::execution_space;
  
  void testAgreementWithPolylib(const int polyOrder, const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    // Polylib Jacobi computes a single function at potentially multiple points
    const int numPoints = 1;
    Kokkos::View<double*,DeviceType> inputPoints("x point", numPoints);
    Kokkos::View<double*,DeviceType> jacobiValues("jacobi values from Intrepid2::Polynomials",polyOrder+1);
    Kokkos::View<double*,DeviceType> polylibJacobiValues("jacobi values from Intrepid2::Polylib",polyOrder+1);
    const Kokkos::View<double*,DeviceType> null; // placeholder View to indicate we don't want the derivatives from Polylib's Jacobi implementation
    
    const double t = 1.0; // this is a scaling value used by shiftedScaledJacobiValues; the 1.0 value corresponds to what's done in Polylib
    
    for (auto alpha : alpha_values)
    {
      const double beta = 0.0; // implementation in Intrepid2::Polynomials assumes beta = 0; the one in Polylib is more general
      for (auto x : x_values)
      {
        Kokkos::deep_copy(inputPoints, x);
        double x_shiftedScaled = (x + 1.0) / 2.0; // shiftedScaledJacobiValues computes on a [0,1] domain, compared with a [-1,1] domain

        using Intrepid2::Polynomials::shiftedScaledJacobiValues;

        auto singlePolicy = Kokkos::RangePolicy<>(ExecutionSpace(),0,1);
        for (int p=0; p<polyOrder+1; p++)
        {
          auto polylibJacobiValue = Kokkos::subview(polylibJacobiValues, std::make_pair(p,p+1));
          // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
          Kokkos::parallel_for(singlePolicy, KOKKOS_LAMBDA(const int dummy_index)
          {
            Intrepid2::Polylib::Serial::JacobiPolynomial(numPoints, inputPoints, polylibJacobiValue, null, p, alpha, beta);
          });
        }
        
        // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
        Kokkos::parallel_for(singlePolicy, KOKKOS_LAMBDA(const int dummy_index)
        {
          shiftedScaledJacobiValues(jacobiValues, alpha, polyOrder, x_shiftedScaled, t);
        });
        
        auto jacobiValuesHost        = getHostCopy(jacobiValues);
        auto polylibJacobiValuesHost = getHostCopy(polylibJacobiValues);
        
        for (int i=0; i<=polyOrder; i++)
        {
          bool valuesAreBothSmall = valuesAreSmall(jacobiValuesHost(i), polylibJacobiValuesHost(i), tol);
          
          if (!valuesAreBothSmall)
          {
            if (! approximatelyEqual(jacobiValuesHost(i), polylibJacobiValuesHost(i), tol) )
            {
              out << "for polyOrder " << i << ", alpha = " << alpha << ", x = " << x << ", t = " << t << ": ";
              out << jacobiValuesHost(i) << " != " << polylibJacobiValuesHost(i);
              out << " (diff = " << abs(jacobiValuesHost(i) - polylibJacobiValuesHost(i));
              out << "; tol = " << tol << ")\n";
              success = false;
            }
          }
        }
      }
    }
  }
  
  void testAgreementWithAnalyticLegendreForP2(const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    const double alpha  = 0;
    const int polyOrder = 2;
    Kokkos::View<double*,DeviceType> jacobiValues("jacobi values from Intrepid2::Polynomials",polyOrder+1);
    
    for (auto t : t_values)
    {
      for (auto x_onMinusOneToOne : x_values)
      {
        double x = (x_onMinusOneToOne + 1.0) / 2.0; // shiftedScaledJacobiValues computes on a [0,1] domain, compared with a [-1,1] domain
        
        using Intrepid2::Polynomials::shiftedScaledJacobiValues;
        
        // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
        auto policy = Kokkos::RangePolicy<>(ExecutionSpace(),0,1);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int dummy_index)
        {
          shiftedScaledJacobiValues(jacobiValues, alpha, polyOrder, x, t);
        });
        
        auto jacobiValuesHost        = getHostCopy(jacobiValues);
        
        const double expectedValue = 6.0 * x * x - 6.0 * x * t + t*t;
        
        const int i = 2;
        bool valuesAreBothSmall = valuesAreSmall(jacobiValuesHost(i), expectedValue, tol);
        
        if (!valuesAreBothSmall)
        {
          if (! approximatelyEqual(jacobiValuesHost(i), expectedValue, tol) )
          {
            out << "for polyOrder " << i << ", alpha = " << alpha << ", x = " << x << ", t = " << t << ": ";
            out << jacobiValuesHost(i) << " != " << expectedValue;
            out << " (diff = " << abs(jacobiValuesHost(i) - expectedValue);
            out << "; tol = " << tol << ")\n";
            success = false;
          }
        }
      }
    }
  }
  
  TEUCHOS_UNIT_TEST( Jacobi, AgreesWithLegendreForP2 )
  {
    const double tol = TEST_TOLERANCE_TIGHT;
    testAgreementWithAnalyticLegendreForP2(tol, out, success);
  }
  
  TEUCHOS_UNIT_TEST( Jacobi, AgreesWithPolylib )
  {
    const int polyOrderMax = 1;
    const double tol = TEST_TOLERANCE_TIGHT;
    testAgreementWithPolylib(polyOrderMax, tol, out, success);
  }
  
} // namespace
