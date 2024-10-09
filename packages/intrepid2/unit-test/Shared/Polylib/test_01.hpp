// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
///////////////////////////////////////////////////////////////////////////////
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description:
// This file is redistributed with the Intrepid package. It should be used
// in accordance with the above MIT license, at the request of the authors.
// This file is NOT covered by the usual Intrepid/Trilinos LGPL license.
//
// Origin: Nektar++ library, http://www.nektar.info, downloaded on
//         March 10, 2009.
//
///////////////////////////////////////////////////////////////////////////////


/** \file   test_01.hpp
    \brief  Test file for a set of functions providing orthogonal polynomial
            polynomial calculus and interpolation.
    \author Created by Spencer Sherwin, Aeronautics, Imperial College London,
            modified and redistributed by D. Ridzal.
            modified and Kokkorized by Kyungjoo Kim
*/

/* -------------------------------------------------------------------
   This is a routine to test the integration, differentiation and
   interpolation routines in the polylib.c.
   
   First, it performs the integral
   
   /1      alpha   beta  alpha,beta
   |   (1-x)   (1+x)     P (x)       dx  = 0
   /-1                    n
   
   for all   -0.5 <= alpha <= 5   (increments of 0.5)
   -0.5 <= beta  <= 5   (increments of 0.5)
   
   using np points where
   NPLOWER <= np <= NPUPPER
   2     <= n  <= 2*np - delta
   
   delta = 1 (gauss), 2(radau), 3(lobatto).
   The integral is evaluated and if it is larger that EPS then the
   value of alpha,beta,np,n and the integral is printed to the screen.
   
   After every alpha value the statement
   "finished checking all beta values for alpha = #"
   is printed
   
   The routine then evaluates the derivate of
   
   d   n      n-1
   -- x  = n x
   dx
   
   for all   -0.5 <= alpha <= 5   (increments of 0.5)
   -0.5 <= beta  <= 5   (increments of 0.5)
   
   using np points where
   NPLOWER <= np <= NPUPPER
   2     <= n  <= np - 1
   
   The error is check in a pointwise sense and if it is larger than
   EPS then the value of alpha,beta,np,n and the error is printed to
   the screen. After every alpha value the statement
   "finished checking all beta values for alpha = #"
   is printed
   
   Finally the routine  evaluates the interpolation of
   
   n      n
   z  to  x
   
   where z are the quadrature zeros and x are the equispaced points
   
   2*i
   x    =   -----   - 1.0    (0 <= i <= np-1)
   i       (np-1)
   
   
   for all   -0.5 <= alpha <= 5   (increments of 0.5)
   -0.5 <= beta  <= 5   (increments of 0.5)
   
   using np points where
   NPLOWER <= np <= NPUPPER
   2     <= n  <= np - 1
   
   The error is check in a pointwise sense and if it is larger than
   EPS then the value of alpha,beta,np,n and the error is printed to
   the screen. After every alpha value the statement
   "finished checking all beta values for alpha = #"
   is printed
   
   The above checks are performed for all the Gauss, Gauss-Radau and
   Gauss-Lobatto points. If you want to disable any routine then set
   GAUSS_INT, GAUSS_RADAU_INT, GAUSS_LOBATTO_INT = 0
   for the integration rouintes
   GAUSS_DIFF,GAUSS_RADAU_DIFF, GAUSS_LOBATTO_DIFF = 0
   for the differentiation routines
   GAUSS_INTERP,GAUSS_RADAU_INTERP, GAUSS_LOBATTO_INTERP = 0
   for the interpolation routines.
   ------------------------------------------------------------------*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Polylib.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {
  
  namespace Test {

#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    }                                                                   \
    catch (std::logic_error &err) {                                      \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }

    template<typename ValueType, 
             typename xViewType,
             typename yViewType>
    ValueType ddot(const ordinal_type n, 
                   const xViewType x,
                   const yViewType y) {
      ValueType sum = 0.;
      for (auto i=0;i<n;++i)
        sum += x(i)*y(i);
      return sum;
    }

    template<typename ValueType, typename DeviceType>
    int Polylib_Test01(const bool verbose) {

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      *outStream
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                       Unit Test (IntrepidPolylib)                           |\n"
        << "|                                                                             |\n"
        << "|     1) the original Polylib tests                                           |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n"
        << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n";

      int errorFlag = 0;

      *outStream
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 1: exceptions                                                          |\n"\
        << "===============================================================================\n";

      try{
        ordinal_type nthrow = 0, ncatch = 0;
        // GammaFunction cannot throw exception
        // INTREPID2_TEST_ERROR_EXPECTED( Polylib::GammaFunction((ValueType)0.33) );
        if (nthrow != ncatch) {
          errorFlag++;
          *outStream << std::setw(70) << "^^^^----FAILURE!" << "\n";
          *outStream << "# of catch ("<< ncatch << ") is different from # of throw (" << ncatch << ")\n";
        }
      } catch (std::logic_error &err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      }

      *outStream
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 2: correctness of math operations                                      |\n"\
        << "===============================================================================\n";

      outStream->precision(5);

      const ordinal_type npLower = 5, npUpper = Polylib::MaxPolylibPoint; // npUpper: 31 right now
      const ValueType tol = 1000.0 * tolerence();
      const double  lowOrderTol = tol;
      const double highOrderTol = tol * 100;

      try {
        Kokkos::View<ValueType*,Kokkos::HostSpace> 
          z("z", npUpper), 
          w("w", npUpper), 
          p("p", npUpper),
          null;
        
        Kokkos::View<ValueType**,Kokkos::HostSpace> 
          d("d", npUpper, npUpper);

        const EPolyType ps[] = {
          POLYTYPE_GAUSS,
          POLYTYPE_GAUSS_RADAU_LEFT,
          POLYTYPE_GAUSS_RADAU_RIGHT,
          POLYTYPE_GAUSS_LOBATTO,
          POLYTYPE_MAX
        };

        const ordinal_type offs[] = { 1, 2, 2, 3 };

        for (auto h=0;ps[h]!=POLYTYPE_MAX;++h) {
          const auto poly = ps[h];
          const auto off  = offs[h];

          *outStream << "Begin checking integration: " << EPolyTypeToString(poly) << "\n";
          {
            ValueType alpha = -0.5;
            while (alpha <= 5.0) {
              ValueType beta = -0.5;
              while (beta <= 5.0) {
                for (auto np = npLower; np <= npUpper; ++np){
                  const double localTol = (np > 20) ? highOrderTol : lowOrderTol;
                  Polylib::Serial::getCubature(z, w, np, alpha, beta, poly);

                  for (auto n = 2; n < 2*np-off; ++n){
                    Polylib::Serial::JacobiPolynomial(np, z, p, null, n, alpha, beta);
                    const ValueType sum = ddot<ValueType>(np, w, p);
                    if (std::isnan(sum) || std::abs(sum) > localTol) {
                      errorFlag = -1000;
                      *outStream << "ERROR:  alpha = " << alpha << ", beta = " << beta <<
                        ", np = " << np << ", n = " << n << "  integral was " << sum << "\n";
                    }
                  }
                }
                beta += 0.5;
              }
              *outStream << "finished checking all beta values for alpha = " << alpha << "\n";
              alpha += 0.5;
            }
          }
          *outStream << "Finished checking integration: " << EPolyTypeToString(poly) << "\n\n";

          *outStream << "Begin checking differentiation: " << EPolyTypeToString(poly) << "\n";
          {
            ValueType alpha = -0.5;
            while (alpha <= 5.0) {
              ValueType beta = -0.5;
              while (beta <= 5.0) {
                for (auto np = npLower; np <= npUpper; ++np) {
                  Polylib::Serial::getCubature(z, w, np, alpha, beta, poly);
                  const double localTol = (np > 20) ? highOrderTol : lowOrderTol;

                  for (auto n = 2; n < np-1; ++n) {
                    Polylib::Serial::getDerivative(d, z, np, alpha, beta, poly);
                    
                    for (auto i = 0; i < np; ++i) 
                      p(i) = std::pow(z(i), n);

                    ValueType sum = 0.0;
                    for (auto i = 0; i < np; ++i)
                      sum += std::abs(ddot<ValueType>(np, Kokkos::subview(d, i, Kokkos::ALL()), p) - n*std::pow(z(i),n-1));
                    sum /= (ValueType)np;
                    if (std::abs(sum)>localTol) {
                      errorFlag = -1000;
                      *outStream << "ERROR:  alpha = " << alpha << ", beta = " << beta <<
                        ", np = " << np << ", n = " << n << "  difference " << sum << "\n";
                    }
                  }
                }
                beta += 0.5;
              }
              *outStream << "finished checking all beta values for alpha = " << alpha << "\n";
              alpha += 0.5;
            }
          }
          *outStream << "Finished checking differentiation: " << EPolyTypeToString(poly) << "\n\n";
          
          *outStream << "Begin checking interpolation: " << EPolyTypeToString(poly) << "\n";
          {
            ValueType alpha = -0.5;
            while (alpha <= 5.0) {
              ValueType beta = -0.5;
              while (beta <= 5.0) {
                
                for (auto np = npLower; np <= npUpper; ++np) {
                  const double localTol = (np > 20) ? highOrderTol : lowOrderTol;
                  Polylib::Serial::getCubature(z, w, np, alpha, beta, poly);

                  for (auto n = 2; n < np-1; ++n) {
                    for (auto i = 0; i < np; ++i) {
                      w(i) = 2.0*i/(ValueType)(np-1)-1.0;
                      p(i) = std::pow(z(i),n);
                    }
                    Polylib::Serial::getInterpolationOperator(d, z, w, np, np, alpha, beta, poly);

                    ValueType sum = 0;
                    for (auto i = 0; i < np; ++i)
                      sum += std::abs(ddot<ValueType>(np, Kokkos::subview(d, i, Kokkos::ALL()), p) - std::pow(w(i),n));
                    sum /= (ValueType)np;
                    if (std::abs(sum)>localTol) {
                      errorFlag = -1000;
                      *outStream << "ERROR:  alpha = " << alpha << ", beta = " << beta <<
                        ", np = " << np << ", n = " << n << "  difference " << sum << "\n";
                    }
                  }
                }
                beta += 0.5;
              }
              *outStream << "finished checking all beta values for alpha = " << alpha << "\n";
              alpha += 0.5;
            }
            *outStream << "Finished checking interpolation: " << EPolyTypeToString(poly) << "\n\n";
          }
        }          
        /******************************************/
        *outStream << "\n";
      } catch (std::logic_error &err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      }

      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";
      
      // reset format state of std::cout
      std::cout.copyfmt(oldFormatState);

      return errorFlag;
    }
  }
}
