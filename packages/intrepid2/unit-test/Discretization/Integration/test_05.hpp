// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Unit test (CubatureDirect,CubatureTensor): correctness of
    integration of monomials for 2D reference cells.
    \author Created by P. Bochev, D. Ridzal and Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"
#include "Intrepid2_Utils_ExtData.hpp"

#include "Intrepid2_CubatureDirectTriDefault.hpp"
#include "Intrepid2_CubatureDirectTriSymmetric.hpp"

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
    int Integration_Test05(const bool verbose) {

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
        << "|     Unit Test (CubatureDirect,CubatureTensor)                               |\n"
        << "|                                                                             |\n"
        << "|     1) Computing integrals of monomials on reference cells in 2D            |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n"
        << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n"
        << "|                      Kyungjoo Kim (kyukim@sandia.gov).                      |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n"
        << "| TEST 1: integrals of monomials in 2D                                        |\n"
        << "===============================================================================\n";
      
      typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

      typedef ValueType pointValueType;
      typedef ValueType weightValueType;
      typedef CubatureDirectTriDefault<DeviceType,pointValueType,weightValueType> CubatureTriType;
      typedef CubatureDirectTriSymmetric<DeviceType,pointValueType,weightValueType> CubatureTriSymType;

      const auto tol = 10.0 * tolerence();

      int errorFlag = 0;

      // compute and compare integrals
      try {
        const auto maxDeg   = Parameters::MaxCubatureDegreeTri;
        // storage for cubatrue points and weights
        DynRankView ConstructWithLabel(cubPoints,
                                       Parameters::MaxIntegrationPoints,
                                       Parameters::MaxDimension);
        
        DynRankView ConstructWithLabel(cubWeights,
                                       Parameters::MaxIntegrationPoints);
        
        // compute integrals
        for (auto cubDeg=0;cubDeg<=20;++cubDeg) {
          CubatureTriType triCub(cubDeg);

          *outStream << "Default Cubature of rder " << std::setw(2) << std::left << cubDeg << " Testing\n";   
          ordinal_type cnt = 0;
          for (auto xDeg=0;xDeg<=cubDeg;++xDeg) {
            for (auto yDeg=0;yDeg<=(cubDeg-xDeg);++yDeg,++cnt) {
              auto computedIntegral = computeIntegralOfMonomial<ValueType>(triCub,
                                                                          cubPoints,
                                                                          cubWeights,
                                                                          xDeg, 
                                                                          yDeg);

              auto anlyticIntegral = analyticIntegralOfMonomialOverTri<ValueType>(xDeg,yDeg);
              const auto abstol  = std::fabs(tol*anlyticIntegral );
              const auto absdiff = std::fabs(anlyticIntegral - computedIntegral);
              if (absdiff > abstol) {
                *outStream << "Default Cubature of order " << std::setw(2) << std::left << cubDeg << " integrating "
                        << "x^" << std::setw(2) << std::left << xDeg << " * y^" << std::setw(2) << yDeg << ":" << "   "
                        << std::scientific << std::setprecision(16) << computedIntegral << "   " << anlyticIntegral << "   "
                        << std::setprecision(4) << absdiff << "   " << "<?" << "   " << abstol << "\n";
                *outStream << std::right << std::setw(111) << "^^^^---FAILURE!\n";
                errorFlag++;
              }
            }
          }
        }

        for (auto cubDeg=0;cubDeg<=maxDeg;++cubDeg) {

          *outStream << "Symmetric Cubature of order " << std::setw(2) << std::left << cubDeg << "  Testing\n";   

          CubatureTriSymType triCub(cubDeg);
          triCub.getCubature(cubPoints, cubWeights);
          bool isInvariant = IsQuadratureInvariantToOrientation(triCub, cubPoints, cubWeights, shards::Triangle<3>::key);

          if (!isInvariant) {
            errorFlag++;
            *outStream << "   Cubature Rule is not invariant to rotations!\n" << std::right << std::setw(111) << "^^^^---FAILURE!\n";
          }

          auto cubWeights_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), cubWeights);
          ValueType minWeigth = 1.0;
          for(int i=0; i<triCub.getNumPoints();++i)
            minWeigth = std::min(minWeigth,cubWeights_host(i));

          if (minWeigth <= 0.0) {
            errorFlag++;
            *outStream << "   Cubature Rule is not positive!\n" << std::right << std::setw(111) << "^^^^---FAILURE!\n";
          }

          
          ordinal_type cnt = 0;
          for (auto xDeg=0;xDeg<=cubDeg;++xDeg) {
            for (auto yDeg=0;yDeg<=(cubDeg-xDeg);++yDeg,++cnt) {
              auto computedIntegral = computeIntegralOfMonomial<ValueType>(triCub,
                                                                          cubPoints,
                                                                          cubWeights,
                                                                          xDeg, 
                                                                          yDeg);
              auto anlyticIntegral = analyticIntegralOfMonomialOverTri<ValueType>(xDeg,yDeg);
              const auto abstol  = std::fabs(tol*anlyticIntegral );
              const auto absdiff = std::fabs(anlyticIntegral - computedIntegral);
              if (absdiff > abstol) {
                *outStream << "Symmetric Cubature order " << std::setw(2) << std::left << cubDeg << " integrating "
                        << "x^" << std::setw(2) << std::left << xDeg << " * y^" << std::setw(2) << yDeg << ":" << "   "
                        << std::scientific << std::setprecision(16) << computedIntegral << "   " << anlyticIntegral << "   "
                        << std::setprecision(4) << absdiff << "   " << "<?" << "   " << abstol << "\n";
                *outStream << std::right << std::setw(111) << "^^^^---FAILURE!\n";
                errorFlag++;
              }
            }
          }
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
