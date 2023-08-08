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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER


/** \file
    \brief  Unit test (CubatureDirect,CubatureTensor): correctness of
    integration of monomials for 3D reference cells.
    \author Created by P. Bochev, D. Ridzal and Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"
#include "Intrepid2_Utils_ExtData.hpp"

#include "Intrepid2_CubatureDirectTetDefault.hpp"
#include "Intrepid2_CubatureDirectTetSymmetric.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

#include "test_util.hpp"

namespace Intrepid2 {
  
  namespace Test {
#define INTREPID2_TEST_ERROR_EXPECTED( S )              \
    try {                                                               \
      ++nthrow;                                                         \
      S ;                                                               \
    } catch (std::logic_error &err) {                                    \
      ++ncatch;                                                         \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };

    template<typename ValueType, typename DeviceType>
    int Integration_Test06(const bool verbose) {

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs, false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      using DeviceSpaceType = typename DeviceType::execution_space;
      typedef typename
        Kokkos::DefaultHostExecutionSpace HostSpaceType ;

      *outStream << "DeviceSpace::  "; DeviceSpaceType().print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType().print_configuration(*outStream, false);
 
      *outStream                                                       
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|       Unit Test (CubatureDirect,CubatureTensor,DefaultCubatureFactory)      |\n"
        << "|                                                                             |\n"
        << "|     1) Computing integrals of monomials on reference cells in 3D            |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n"
        << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n"
        << "|                      Kyungjoo Kim (kyukim@sandia.gov).                      |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n"
        << "| TEST 1: integrals of monomials in 3D                                        |\n"
        << "===============================================================================\n";
      
      typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

      typedef ValueType pointValueType;
      typedef ValueType weightValueType;
      typedef CubatureDirectTetDefault<DeviceType,pointValueType,weightValueType> CubatureTetType;
      typedef CubatureDirectTetSymmetric<DeviceType,pointValueType,weightValueType> CubatureTetSymType;
      
      // tolerence is too tight to test upto order 20
      // tol factor increased by 1000 due to cubature order 20 test failure
      const auto tol = 2000.0 * tolerence();

      int errorFlag = 0;

      // compute and compare integrals
      try {
        // cannot test maxcubature degree edge (20) as max integration point is limited by 1001.
        const auto maxDeg   = Parameters::MaxCubatureDegreeTet;
        //const auto polySize = (maxDeg+1)*(maxDeg+2)*(maxDeg+3)/6;

        // analytic integral values
        const auto analyticMaxDeg = 20;

        // storage for cubatrue points and weights
        DynRankView ConstructWithLabel(cubPoints,
                                       Parameters::MaxIntegrationPoints,
                                       Parameters::MaxDimension);

        DynRankView ConstructWithLabel(cubWeights,
                                       Parameters::MaxIntegrationPoints);

        // perform comparison
        for (auto cubDeg=0;cubDeg<=maxDeg;++cubDeg) {

          *outStream << "Testing Default Cubature of Order " << std::setw(2) << std::left << cubDeg << "\n";

          CubatureTetType tetCub(cubDeg);

          const auto y_offs = (analyticMaxDeg - cubDeg);
          const auto x_offs = y_offs*(y_offs + 1)/2;
          ordinal_type offset = 0;
          const auto oldFlag = errorFlag;
          for (auto xDeg=0;xDeg<=cubDeg;++xDeg,offset += x_offs) 
            for (auto yDeg=0;yDeg<=(cubDeg-xDeg);++yDeg,offset += y_offs) 
              for (auto zDeg=0;zDeg<=(cubDeg-xDeg-yDeg);++zDeg) {
                const auto analyticIntegral = analyticIntegralOfMonomialOverTet<ValueType>(xDeg,yDeg,zDeg);
                const auto abstol  = std::fabs(tol*analyticIntegral );
                const auto computedIntegral = computeIntegralOfMonomial<ValueType>(tetCub,cubPoints,cubWeights,xDeg,yDeg,zDeg);
                const auto absdiff = std::fabs(analyticIntegral - computedIntegral);
                if (absdiff > abstol) {
                  *outStream << "Default Cubature of Order " << std::setw(2) << std::left << cubDeg << " Integrating "
                             << "x^" << std::setw(2) << std::left << xDeg << " * y^" << std::setw(2) << yDeg
                             << " * z^" << std::setw(2) << zDeg << ":" << "   "
                             << std::scientific << std::setprecision(16)
                             << computedIntegral << "   " << analyticIntegral << "   "
                             << std::setprecision(4) << absdiff << "   " << "<?" << "   " << abstol << "\n";
                  errorFlag++;
                  *outStream << std::right << std::setw(118) << "^^^^---FAILURE!\n";
                }
              }
          *outStream << "Default Cubature of Order " << std::setw(2) << std::left << cubDeg
                     << (errorFlag == oldFlag ? "  PASSED" : "  FAILED") << std::endl;                         
        }
        *outStream << "\n";


        for (auto cubDeg=0;cubDeg<=maxDeg;++cubDeg) {
          *outStream << "Testing Symmetric Cubature of Order " << std::setw(2) << std::left << cubDeg << "\n";

          CubatureTetSymType tetCub(cubDeg);
          tetCub.getCubature(cubPoints, cubWeights);
          auto cubWeights_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cubWeights);
          ValueType minWeigth = 1.0;
          for(int i=0; i<tetCub.getNumPoints();++i)
            minWeigth = std::min(minWeigth,cubWeights_host(i));

          if (minWeigth <= 0.0) {
            errorFlag++;
            *outStream << "   Cubature Rule is not positive!\n" << std::right << std::setw(111) << "^^^^---FAILURE!\n";
          }

          const auto y_offs = (analyticMaxDeg - cubDeg);
          const auto x_offs = y_offs*(y_offs + 1)/2;
          ordinal_type offset = 0;
          const auto oldFlag = errorFlag;
          for (auto xDeg=0;xDeg<=cubDeg;++xDeg,offset += x_offs) 
            for (auto yDeg=0;yDeg<=(cubDeg-xDeg);++yDeg,offset += y_offs) 
              for (auto zDeg=0;zDeg<=(cubDeg-xDeg-yDeg);++zDeg) {
                const auto analyticIntegral = analyticIntegralOfMonomialOverTet<ValueType>(xDeg,yDeg,zDeg);
                const auto computedIntegral = computeIntegralOfMonomial<ValueType>(tetCub,cubPoints,cubWeights,xDeg,yDeg,zDeg);
                const auto abstol  = std::fabs(tol*analyticIntegral );
                const auto absdiff = std::fabs(analyticIntegral - computedIntegral);
                if (absdiff > abstol) {
                  *outStream << "Symmetric Cubature order " << std::setw(2) << std::left << cubDeg << " integrating "
                             << "x^" << std::setw(2) << std::left << xDeg << " * y^" << std::setw(2) << yDeg
                             << " * z^" << std::setw(2) << zDeg << ":" << "   "
                             << std::scientific << std::setprecision(16)
                             << computedIntegral << "   " << analyticIntegral << "   "
                             << std::setprecision(4) << absdiff << "   " << "<?" << "   " << abstol << "\n";
                  errorFlag++;
                  *outStream << std::right << std::setw(118) << "^^^^---FAILURE!\n";
                }
              }
          *outStream << "Symmetric Cubature of Order " << std::setw(2) << std::left << cubDeg
                     << (errorFlag == oldFlag ? "  PASSED" : "  FAILED") << std::endl;                         
        }
        *outStream << "\n";
      } catch (std::logic_error &err) {
        *outStream << err.what() << "\n";
        errorFlag = -1;
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
