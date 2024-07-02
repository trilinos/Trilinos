// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   IntegratedJacobiTests.cpp
    \brief  Tests to verify that implementations of integrated Jacobi polynomials satisfy some basic consistency requirements.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_Types.hpp"

#include "Intrepid2_TestUtils.hpp"

#include "Kokkos_Core.hpp"

namespace
{
  using namespace Intrepid2;
  
  // x, t, alpha values used across these tests
  std::vector<double> alpha_values = {{0.0, 0.2, 0.4, 0.6, 0.8, 1.0}};
  std::vector<double> t_values = {{0.0, 0.2, 0.4, 0.6, 0.8, 1.0}};
  std::vector<double> x_values = {{-1.0,-0.5,-1.0/3.0,0.0,1.0/3.0,0.50,1.0}};

  using DeviceType = DefaultTestDeviceType;
  using ExecutionSpace = typename DeviceType::execution_space;
  
  void testIntegratedJacobiIsZeroAtZero(const int polyOrder, const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    // for all alpha, t, integrated jacobi should evaluate to 0 at 0.
    Kokkos::View<double*,DeviceType> integratedJacobiView("integrated jacobi values",polyOrder+1);
    
    using Intrepid2::Polynomials::shiftedScaledIntegratedJacobiValues;
    using Intrepid2::Polynomials::shiftedScaledJacobiValues;

    const double x = 0.0;
    
    for (auto alpha : alpha_values)
    {
      for (auto t : t_values)
      {
        // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
        auto policy = Kokkos::RangePolicy<>(ExecutionSpace(),0,1);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int dummy_index)
        {
          shiftedScaledIntegratedJacobiValues(integratedJacobiView, alpha, polyOrder, x, t);
        });
        
        Kokkos::fence();
        auto integratedJacobiViewHost = getHostCopy(integratedJacobiView);
        
        for (int i=1; i<=polyOrder; i++)
        {
          if ( abs(integratedJacobiViewHost(i)) > tol)
          {
            success = false;
            out << "for alpha = " << alpha << ", t = " << t << ", integrated Jacobi for polyOrder " << i;
            out << " at x=0 is not zero (it is " << integratedJacobiViewHost(i) << ")\n";
          }
        }
      }
    }
  }
  
  void testIntegratedJacobiAnalyticAlphaZeroJ2(const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    // for alpha=0, integrated jacobi j=2 should be x^2 - x*t
    const int polyOrder = 2;
    Kokkos::View<double*,DeviceType> integratedJacobiView("integrated jacobi values",polyOrder+1);
    
    using Intrepid2::Polynomials::shiftedScaledIntegratedJacobiValues;
    
    const double alpha = 0.0;
    
    for (auto x : x_values)
    {
      for (auto t : t_values)
      {
        double expected_value = x * x - x * t;
        // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
        auto policy = Kokkos::RangePolicy<>(ExecutionSpace(),0,1);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int dummy_index)
        {
          shiftedScaledIntegratedJacobiValues(integratedJacobiView, alpha, polyOrder, x, t);
        });
        
        Kokkos::fence();
        
        auto integratedJacobiViewHost = getHostCopy(integratedJacobiView);
        const int i = 2;
        double diff = integratedJacobiViewHost(i) - expected_value;
        if ( abs(diff) > tol)
        {
          success = false;
          out << "for alpha = " << alpha << ", t = " << t << ", integrated Jacobi for polyOrder " << i;
          out << " at x=" << x << " is not x^2 - x * t (" << expected_value << "); instead, it is " << integratedJacobiViewHost(i);
          out << ", a difference of " << abs(diff) << ")\n";
        }
      }
    }
  }
  
  void testIntegratedJacobiTwoPathsMatch(const int polyOrder, const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    Kokkos::View<double*,DeviceType> jacobiView("jacobi values",polyOrder+1);
    Kokkos::View<double*,DeviceType> integratedJacobiViewSecondPath("integrated jacobi values (second path)",polyOrder+1);
    Kokkos::View<double*,DeviceType> integratedJacobiView("integrated jacobi values",polyOrder+1);

    for (auto alpha : alpha_values)
    {
      for (auto x : x_values)
      {
        for (auto t : t_values)
        {
          using Intrepid2::Polynomials::shiftedScaledIntegratedJacobiValues;
          using Intrepid2::Polynomials::shiftedScaledJacobiValues;
          
          // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
          auto policy = Kokkos::RangePolicy<>(ExecutionSpace(),0,1);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int dummy_index)
          {
            shiftedScaledJacobiValues(jacobiView, alpha, polyOrder, x, t);
            shiftedScaledIntegratedJacobiValues(integratedJacobiView, alpha, polyOrder, x, t);
            shiftedScaledIntegratedJacobiValues(integratedJacobiViewSecondPath, jacobiView, alpha, polyOrder, x, t);
          });
          
          auto integratedJacobiViewSecondPathHost = getHostCopy(integratedJacobiViewSecondPath);
          auto integratedJacobiViewHost           = getHostCopy(integratedJacobiView);
          
          for (int i=0; i<=polyOrder; i++)
          {
            bool valuesAreBothSmall = valuesAreSmall(integratedJacobiViewSecondPathHost(i), integratedJacobiViewHost(i), tol);
            if (!valuesAreBothSmall)
            {
              if (! approximatelyEqual(integratedJacobiViewSecondPathHost(i), integratedJacobiViewHost(i), tol) )
              {
                out << "for polyOrder " << i << ", alpha = " << alpha << ", x = " << x << ", t = " << t << ": ";
                out << integratedJacobiViewSecondPathHost(i) << " != " << integratedJacobiViewHost(i);
                out << " (diff = " << abs(integratedJacobiViewSecondPathHost(i) - integratedJacobiViewHost(i));
                out << "; tol = " << tol << ")\n";
                success = false;
              }
            }
          }
        }
      }
    }
  }
  
  void testIntegratedJacobi_dtTwoPathsMatch(const int polyOrder, const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    Kokkos::View<double*,DeviceType> jacobiView("jacobi values",polyOrder+1);
    Kokkos::View<double*,DeviceType> integratedJacobiView_dt("d/dt(integrated jacobi) values",polyOrder+1);
    Kokkos::View<double*,DeviceType> secondPathIntegratedJacobiView_dt("d/dt(integrated jacobi) values (second path)",polyOrder+1);
    
    for (auto alpha : alpha_values)
    {
      for (auto x : x_values)
      {
        for (auto t : t_values)
        {
          using Intrepid2::Polynomials::shiftedScaledIntegratedJacobiValues_dt;
          using Intrepid2::Polynomials::shiftedScaledJacobiValues;
          
          // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
          auto policy = Kokkos::RangePolicy<>(ExecutionSpace(),0,1);
          Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int dummy_index)
          {
            shiftedScaledJacobiValues(jacobiView, alpha, polyOrder, x, t);
            shiftedScaledIntegratedJacobiValues_dt(integratedJacobiView_dt, alpha, polyOrder, x, t);
            shiftedScaledIntegratedJacobiValues_dt(secondPathIntegratedJacobiView_dt, jacobiView, alpha, polyOrder, x, t);
          });
          
          auto secondPathIntegratedJacobiView_dt_host = getHostCopy(secondPathIntegratedJacobiView_dt);
          auto integratedJacobiView_dt_host           = getHostCopy(integratedJacobiView_dt);
          
          for (int i=0; i<=polyOrder; i++)
          {
            bool valuesAreBothSmall = valuesAreSmall(secondPathIntegratedJacobiView_dt_host(i), integratedJacobiView_dt_host(i), tol);
            if (!valuesAreBothSmall)
            {
              if (! approximatelyEqual(secondPathIntegratedJacobiView_dt_host(i), integratedJacobiView_dt_host(i), tol) )
              {
                out << "for polyOrder " << i << ", alpha = " << alpha << ", x = " << x << ", t = " << t << ": ";
                out << secondPathIntegratedJacobiView_dt_host(i) << " != " << integratedJacobiView_dt_host(i);
                out << " (diff = " << abs(secondPathIntegratedJacobiView_dt_host(i) - integratedJacobiView_dt_host(i));
                out << "; tol = " << tol << ")\n";
                success = false;
              }
            }
          }
        }
      }
    }
  }
  
  TEUCHOS_UNIT_TEST( IntegratedJacobi, TwoPathsMatch )
  {
    const int polyOrderMax = 10;
    const double tol = TEST_TOLERANCE_TIGHT * 1.0e2; // 10th order fails for a single test case on some CUDA platforms with TEST_TOLERANCE_TIGHT, so we relax this a bit.
    testIntegratedJacobiTwoPathsMatch(polyOrderMax, tol, out, success);
  }
  
  TEUCHOS_UNIT_TEST( IntegratedJacobi, dtTwoPathsMatch )
  {
    const int polyOrderMax = 10;
    const double tol = TEST_TOLERANCE_TIGHT * 1.0e2; // prefer to keep the tolerances the same for these two tests.
    testIntegratedJacobi_dtTwoPathsMatch(polyOrderMax, tol, out, success);
  }
  
  TEUCHOS_UNIT_TEST( IntegratedJacobi, ZeroAtZero)
  {
    const int polyOrderMax = 10;
    const double tol = TEST_TOLERANCE_TIGHT;
    testIntegratedJacobiIsZeroAtZero(polyOrderMax, tol, out, success);
  }
  
  TEUCHOS_UNIT_TEST( IntegratedJacobi, AnalyticAlphaZeroJ2)
  {
    // analytically, we can determine that with alpha=0, the second integrated Jacobi polynomial
    // should be x^2 - x*t
    const double tol = TEST_TOLERANCE_TIGHT;
    testIntegratedJacobiAnalyticAlphaZeroJ2(tol, out, success);
  }
} // namespace
