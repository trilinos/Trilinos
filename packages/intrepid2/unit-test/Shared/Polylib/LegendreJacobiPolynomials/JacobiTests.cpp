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
  
  void testAgreementWithPolylib(const int polyOrder, const double tol, Teuchos::FancyOStream &out, bool &success)
  {
    // Polylib Jacobi computes a single function at potentially multiple points
    const int numPoints = 1;
    Kokkos::View<double*> inputPoints("x point", numPoints);
    Kokkos::View<double*> jacobiValues("jacobi values from Intrepid2::Polynomials",polyOrder+1);
    Kokkos::View<double*> polylibJacobiValues("jacobi values from Intrepid2::Polylib",polyOrder+1);
    const Kokkos::View<double*> null; // placeholder View to indicate we don't want the derivatives from Polylib's Jacobi implementation
    
    const double t = 1.0; // this is a scaling value used by shiftedScaledJacobiValues; the 1.0 value corresponds to what's done in Polylib
    
    for (auto alpha : alpha_values)
    {
      const double beta = 0.0; // implementation in Intrepid2::Polynomials assumes beta = 0; the one in Polylib is more general
      for (auto x : x_values)
      {
        Kokkos::deep_copy(inputPoints, x);
        double x_shiftedScaled = (x + 1.0) / 2.0; // shiftedScaledJacobiValues computes on a [0,1] domain, compared with a [-1,1] domain

        using Intrepid2::Polynomials::shiftedScaledJacobiValues;

        for (int p=0; p<polyOrder+1; p++)
        {
          auto polylibJacobiValue = Kokkos::subview(polylibJacobiValues, std::make_pair(p,p+1));
          // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
          Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int dummyIndex)
          {
            Intrepid2::Polylib::Serial::JacobiPolynomial(numPoints, inputPoints, polylibJacobiValue, null, p, alpha, beta);
          });
        }
        
        // wrap invocation in parallel_for just to ensure execution on device (for CUDA)
        Kokkos::parallel_for(1, KOKKOS_LAMBDA(const int dummyIndex)
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
  
  TEUCHOS_UNIT_TEST( Jacobi, AgreesWithPolylib )
  {
    const int polyOrderMax = 1;
    const double tol = TEST_TOLERANCE_TIGHT;
    testAgreementWithPolylib(polyOrderMax, tol, out, success);
  }
  
} // namespace
