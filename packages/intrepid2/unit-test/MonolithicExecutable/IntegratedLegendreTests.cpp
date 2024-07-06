// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   IntegratedLegendreTests.cpp
    \brief  Tests to verify various implementations of integrated Legendre polynomials.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_Polynomials.hpp"

#include "Intrepid2_TestUtils.hpp"

#include "Kokkos_Core.hpp"

namespace
{
  using namespace Intrepid2;
  
  // x and t values used across these tests
  std::vector<double> t_values = {{0.0, 0.2, 0.4, 0.6, 0.8, 1.0}};
  std::vector<double> x_values = {{-1.0,-0.5,-1.0/3.0,0.0,1.0/3.0,0.50,1.0}};

  using DeviceType = DefaultTestDeviceType;
  using ExecutionSpace = typename DeviceType::execution_space;
  
  void testIntegratedLegendreTwoPathsMatch(const int polyOrder, const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    Kokkos::View<double*,DeviceType> integrated_legendre_values("integrated legendre values",polyOrder+1);
    Kokkos::View<double*,DeviceType> shifted_scaled_legendre_values("shifted scaled legendre values",polyOrder+1);
    Kokkos::View<double*,DeviceType> integrated_legendre_values_second_path("integrated legendre values (second path)",polyOrder+1);
    
    for (auto x : x_values)
    {
      for (auto t : t_values)
      {
        using Intrepid2::Polynomials::shiftedScaledIntegratedLegendreValues;
        using Intrepid2::Polynomials::shiftedScaledLegendreValues;
        
        // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
        auto policy = Kokkos::RangePolicy<>(ExecutionSpace(),0,1);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int dummy_index)
        {
          shiftedScaledIntegratedLegendreValues(integrated_legendre_values, polyOrder, x, t);
          shiftedScaledLegendreValues(shifted_scaled_legendre_values, polyOrder, x, t);
          shiftedScaledIntegratedLegendreValues(integrated_legendre_values_second_path, shifted_scaled_legendre_values, polyOrder, x, t);
        });

        auto integrated_legendre_values_host             = getHostCopy(integrated_legendre_values);
        auto shifted_scaled_legendre_values_host         = getHostCopy(shifted_scaled_legendre_values);
        auto integrated_legendre_values_second_path_host = getHostCopy(integrated_legendre_values_second_path);
        
        for (int i=0; i<=polyOrder; i++)
        {
          bool valuesAreBothSmall = valuesAreSmall(integrated_legendre_values_second_path_host(i), integrated_legendre_values_host(i), tol);
          if (! valuesAreBothSmall)
          {
            if (! approximatelyEqual(integrated_legendre_values_second_path_host(i), integrated_legendre_values_host(i), tol) )
            {
              out << "for polyOrder " << i << ", x = " << x << ", t = " << t << ": ";
              out << integrated_legendre_values_second_path_host(i) << " != " << integrated_legendre_values_host(i);
              out << " (diff = " << abs(integrated_legendre_values_second_path_host(i) - integrated_legendre_values_host(i));
              out << "; tol = " << tol << ")\n";
              success = false;
            }
          }
        }
      }
    }
  }
  
  void testIntegratedLegendre_dtTwoPathsMatch(const int polyOrder, const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    Kokkos::View<double*,DeviceType> integrated_legendre_values_dt("d/dt(integrated legendre) values",polyOrder+1);
    Kokkos::View<double*,DeviceType> shifted_scaled_legendre_values("shifted scaled legendre values",polyOrder+1);
    Kokkos::View<double*,DeviceType> integrated_legendre_values_dt_second_path("d/dt(integrated legendre) values (second path)",polyOrder+1);
    
    for (auto x : x_values)
    {
      for (auto t : t_values)
      {
        using Intrepid2::Polynomials::shiftedScaledIntegratedLegendreValues_dt;
        using Intrepid2::Polynomials::shiftedScaledLegendreValues;
        
        // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
        auto policy = Kokkos::RangePolicy<>(ExecutionSpace(),0,1);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int dummy_index)
        {
          shiftedScaledIntegratedLegendreValues_dt(integrated_legendre_values_dt, polyOrder, x, t);
          shiftedScaledLegendreValues(shifted_scaled_legendre_values, polyOrder, x, t);
          shiftedScaledIntegratedLegendreValues_dt(integrated_legendre_values_dt_second_path, shifted_scaled_legendre_values, polyOrder, x, t);
        });
        
        auto integrated_legendre_values_dt_host             = getHostCopy(integrated_legendre_values_dt);
        auto integrated_legendre_values_dt_second_path_host = getHostCopy(integrated_legendre_values_dt_second_path);
        
        for (int i=0; i<=polyOrder; i++)
        {
          bool valuesAreBothSmall = valuesAreSmall(integrated_legendre_values_dt_second_path_host(i), integrated_legendre_values_dt_host(i), tol);
          if (!valuesAreBothSmall)
          {
            if (! approximatelyEqual(integrated_legendre_values_dt_second_path_host(i), integrated_legendre_values_dt_host(i), tol) )
            {
              out << "for polyOrder " << i << ", x = " << x << ", t = " << t << ": ";
              out << integrated_legendre_values_dt_second_path_host(i) << " != " << integrated_legendre_values_dt_host(i);
              out << " (diff = " << abs(integrated_legendre_values_dt_second_path_host(i) - integrated_legendre_values_dt_host(i));
              out << "; tol = " << tol << ")\n";
              success = false;
            }
          }
        }
      }
    }
  }
  
  void testIntegratedLegendreMatchesFormula(const int polyOrder, const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    Kokkos::View<double*,DeviceType> integrated_legendre_values("integrated legendre values",polyOrder+1);
    Kokkos::View<double*,DeviceType> shifted_scaled_legendre_values("shifted scaled legendre values",polyOrder+1);
    
    for (auto x : x_values)
    {
      for (auto t : t_values)
      {
        using Intrepid2::Polynomials::shiftedScaledIntegratedLegendreValues;
        using Intrepid2::Polynomials::shiftedScaledLegendreValues;
        
        // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
        auto policy = Kokkos::RangePolicy<>(ExecutionSpace(),0,1);
        Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int dummy_index)
        {
          shiftedScaledIntegratedLegendreValues(integrated_legendre_values, polyOrder, x, t);
          shiftedScaledLegendreValues(shifted_scaled_legendre_values, polyOrder, x, t);
        });
        
        auto integrated_legendre_values_host     = getHostCopy(integrated_legendre_values);
        auto shifted_scaled_legendre_values_host = getHostCopy(shifted_scaled_legendre_values);
        
        for (int i=0; i<=polyOrder; i++)
        {
          double L_i_actual = integrated_legendre_values_host(i);
          double L_i_expected;
          if (i == 0)
          {
            L_i_expected = 1.0;
          }
          else if (i == 1)
          {
            L_i_expected = x;
          }
          else
          {
            double i_factor = 2. * (2.* i -1.);
            double P_i = shifted_scaled_legendre_values_host(i);
            double P_i_minus_2 = shifted_scaled_legendre_values_host(i-2);
            L_i_expected = (P_i - t * t * P_i_minus_2) / i_factor;
          }
          
          bool valuesMatch = true;
          bool valuesAreBothSmall = valuesAreSmall(L_i_actual, L_i_expected, tol);
          if (!valuesAreBothSmall)
          {
            TEUCHOS_TEST_FLOATING_EQUALITY(L_i_actual, L_i_expected, tol, out, valuesMatch);
          }
      
          if ( !valuesMatch )
          {
            out << "for polyOrder " << i << ", x = " << x << ", t = " << t << ": ";
            out << L_i_actual << " != " << L_i_expected;
            out << " (diff = " << abs(L_i_actual - L_i_expected);
            out << "; tol = " << tol << ")\n";
            success = false;
          }
        }
      }
    }
  }
  
  // test disabled.  See note by commented-out integratedLegendreValues() in Intrepid_Polynomials.hpp.
  //bool testIntegratedLegendreAgreesWithShiftedScaled(const double x, const int polyOrder, const double tol)
  //{
  //  bool success = true;
  //
  //  Kokkos::View<double*,DeviceType> integrated_legendre_values("integrated legendre values",polyOrder+1);
  //
  //  using Intrepid2::PolyLib::integratedLegendreValues;
  //  using Intrepid2::PolyLib::shiftedScaledIntegratedLegendreValues;
  //  using Intrepid2::PolyLib::shiftedScaledLegendreValues;
  //
  //  integratedLegendreValues(integrated_legendre_values, polyOrder, x);
  //
  //  Kokkos::View<double*,DeviceType> shifted_scaled_integrated_legendre_values("shifted scaled integrated legendre values",polyOrder+1);
  //  using Intrepid2::PolyLib::shiftedScaledIntegratedLegendreValues;
  //
  //  const double x_scaled = (x + 1.0) / 2.0;
  //  const double t = 1.0;
  //  shiftedScaledIntegratedLegendreValues(shifted_scaled_integrated_legendre_values, polyOrder, x_scaled, t);
  //
  //  for (int i=0; i<=polyOrder; i++)
  //  {
  //    if (! approximatelyEqual(shifted_scaled_integrated_legendre_values(i), integrated_legendre_values(i), tol) )
  //    {
  //      success = false;
  //      double diff = std::abs(shifted_scaled_integrated_legendre_values(i) - integrated_legendre_values(i));
  //      std::cout << "FAILURE: testIntegratedLegendreAgreesWithShiftedScaled expected ";
  //      std::cout << shifted_scaled_integrated_legendre_values(i) << " for P_" << i << ", but got ";
  //      std::cout << integrated_legendre_values(i) << " (diff: " << diff << ")" << std::endl;
  //    }
  //  }
  //  return success;
  //}
  
  TEUCHOS_UNIT_TEST( IntegratedLegendre, TwoPathsMatch )
  {
    const int polyOrderMax = 10;
    const double tol = TEST_TOLERANCE_TIGHT;
    testIntegratedLegendreTwoPathsMatch(polyOrderMax, tol, out, success);
  }
  
  TEUCHOS_UNIT_TEST( IntegratedLegendre, dtTwoPathsMatch )
  {
    const int polyOrderMax = 10;
    const double tol = TEST_TOLERANCE_TIGHT;
    testIntegratedLegendre_dtTwoPathsMatch(polyOrderMax, tol, out, success);
  }
  
  TEUCHOS_UNIT_TEST( IntegratedLegendre, MatchesFormula )
  {
    const int polyOrderMax = 10;
    const double tol = TEST_TOLERANCE_TIGHT;
    testIntegratedLegendreMatchesFormula(polyOrderMax, tol, out, success);
  }
} // namespace
