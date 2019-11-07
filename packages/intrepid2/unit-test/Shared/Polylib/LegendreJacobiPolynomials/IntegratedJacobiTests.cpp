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
  
  void testIntegratedJacobiTwoPathsMatch(const int polyOrder, const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    Kokkos::View<double*> jacobiView("jacobi values",polyOrder+1);
    Kokkos::View<double*> integratedJacobiViewSecondPath("integrated jacobi values (second path)",polyOrder+1);
    Kokkos::View<double*> integratedJacobiView("integrated jacobi values",polyOrder+1);

    for (auto alpha : alpha_values)
    {
      for (auto x : x_values)
      {
        for (auto t : t_values)
        {
          using Intrepid2::Polynomials::integratedJacobiValues;
          using Intrepid2::Polynomials::shiftedScaledJacobiValues;
          
          // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
          Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int dummy_index)
          {
            shiftedScaledJacobiValues(jacobiView, alpha, polyOrder, x, t);
            integratedJacobiValues(integratedJacobiView, alpha, polyOrder, x, t);
            integratedJacobiValues(integratedJacobiViewSecondPath, jacobiView, alpha, polyOrder, x, t);
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
    Kokkos::View<double*> jacobiView("jacobi values",polyOrder+1);
    Kokkos::View<double*> integratedJacobiView_dt("d/dt(integrated jacobi) values",polyOrder+1);
    Kokkos::View<double*> secondPathIntegratedJacobiView_dt("d/dt(integrated jacobi) values (second path)",polyOrder+1);
    
    for (auto alpha : alpha_values)
    {
      for (auto x : x_values)
      {
        for (auto t : t_values)
        {
          using Intrepid2::Polynomials::integratedJacobiValues_dt;
          using Intrepid2::Polynomials::shiftedScaledJacobiValues;
          
          // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
          Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int dummy_index)
          {
            shiftedScaledJacobiValues(jacobiView, alpha, polyOrder, x, t);
            integratedJacobiValues_dt(integratedJacobiView_dt, alpha, polyOrder, x, t);
            integratedJacobiValues_dt(secondPathIntegratedJacobiView_dt, jacobiView, alpha, polyOrder, x, t);
          });
          
          auto secondPathIntegratedJacobiView_dt_host = getHostCopy(secondPathIntegratedJacobiView_dt);
          auto integratedJacobiView_dt_host           = getHostCopy(integratedJacobiView_dt);
          
          for (int i=0; i<=polyOrder; i++)
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
  
} // namespace
