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

#include "Intrepid2_CubatureDirectLineGauss.hpp"
#include "Intrepid2_CubatureTensor.hpp"

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
    int Integration_Test04(const bool verbose) {

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
      
      typedef Kokkos::DynRankView<ValueType,Kokkos::HostSpace> DynRankViewHost;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

      typedef ValueType pointValueType;
      typedef ValueType weightValueType;
      typedef CubatureDirectLineGauss<DeviceType,pointValueType,weightValueType> CubatureLineType;
      typedef CubatureTensor<DeviceType,pointValueType,weightValueType> CubatureTensorType;
      
      const auto tol = 10.0 * tolerence();

      int errorFlag = 0;

      // get names of files with analytic values
      std::string basedir = "./data";
      std::stringstream namestream;
      std::string filename;
      namestream << basedir << "/HEX_integrals" << ".dat";
      namestream >> filename;
      *outStream << "filename = " << filename << std::endl;
      std::ifstream filecompare(filename);

      // compute and compare integrals
      try {
        // cannot test maxcubature degree edge (20) as max integration point is limited by 1001.
        /// 19 is too large
        const auto maxDeg   = 10; // 19; //Parameters::MaxCubatureDegreeEdge;
        const auto polySize = (maxDeg+1)*(maxDeg+2)*(maxDeg+3)/6;

        // test integral values
        DynRankViewHost ConstructWithLabel(testInt, maxDeg+1, polySize);

        // analytic integral values
        const auto analyticMaxDeg = 61;
        const auto analyticPolySize = (analyticMaxDeg+1)*(analyticMaxDeg+2)*(analyticMaxDeg+3)/6;

        DynRankViewHost ConstructWithLabel(analyticInt, analyticPolySize, 1);

        // compute integrals
        for (auto cubDeg=0;cubDeg<=maxDeg;++cubDeg) {
          CubatureLineType line(cubDeg);
          CubatureTensorType hexCub( line, line, line );
          *outStream << "Cubature order " << std::setw(2) << std::left << cubDeg << "  Testing\n";
          
          Kokkos::Array<int,3> degrees;
          
          auto cubPoints  = hexCub.allocateCubaturePoints();
          auto cubWeights = hexCub.allocateCubatureWeights();
          hexCub.getCubature(cubPoints, cubWeights);
          
          ordinal_type cnt = 0;
          for (auto xDeg=0;xDeg<=cubDeg;++xDeg)
          {
            degrees[0] = xDeg;
            for (auto yDeg=0;yDeg<=(cubDeg-xDeg);++yDeg)
            {
              degrees[1] = yDeg;
              for (auto zDeg=0;zDeg<=(cubDeg-xDeg-yDeg);++zDeg,++cnt) {
                degrees[2] = zDeg;
                testInt(cubDeg, cnt) = computeIntegralOfMonomial<ValueType>(cubPoints,
                                                                            cubWeights,
                                                                            degrees);
              }
            }
          }
        }


        // get analytic values
        if (filecompare.is_open()) {
          getAnalytic(analyticInt, filecompare);
          filecompare.close();
        }

        // perform comparison
        for (auto cubDeg=0;cubDeg<=maxDeg;++cubDeg) {

          const auto y_offs = (analyticMaxDeg - cubDeg);
          const auto x_offs = y_offs*(y_offs + 1)/2;

          ordinal_type offset = 0, cnt = 0;
          const auto oldFlag = errorFlag;
          for (auto xDeg=0;xDeg<=cubDeg;++xDeg,offset += x_offs) 
            for (auto yDeg=0;yDeg<=(cubDeg-xDeg);++yDeg,offset += y_offs) 
              for (auto zDeg=0;zDeg<=(cubDeg-xDeg-yDeg);++zDeg,++cnt) {
                const auto loc = cnt + offset;
                const auto abstol  = ( analyticInt(loc,0) == 0.0 ? tol : std::fabs(tol*analyticInt(loc,0)) );
                const auto absdiff = std::fabs(analyticInt(loc,0) - testInt(cubDeg,cnt));
                if (absdiff > abstol) {
                  *outStream << "Cubature order " << std::setw(2) << std::left << cubDeg << " integrating "
                             << "x^" << std::setw(2) << std::left << xDeg << " * y^" << std::setw(2) << yDeg
                             << " * z^" << std::setw(2) << zDeg << ":" << "   "
                             << std::scientific << std::setprecision(16)
                             << testInt(cubDeg,cnt) << "   " << analyticInt(loc,0) << "   "
                             << std::setprecision(4) << absdiff << "   " << "<?" << "   " << abstol << "\n";
                  errorFlag++;
                  *outStream << std::right << std::setw(118) << "^^^^---FAILURE!\n";
                }
              }
          *outStream << "Cubature order " << std::setw(2) << std::left << cubDeg
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
