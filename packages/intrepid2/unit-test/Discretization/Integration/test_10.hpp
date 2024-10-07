// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file
    \brief  Unit test (CubaturePolylib): correctness of
    integration of monomials for 1D reference cells.
    \author Created by P. Bochev, D. Ridzal and Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"
#include "Intrepid2_Utils_ExtData.hpp"

#include "Intrepid2_CubaturePolylib.hpp"

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

    template<typename ValueType, typename DeviceSpaceType>
    int Integration_Test10(const bool verbose) {

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs, false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      *outStream
        << "===============================================================================\n"
        << "|                                                                             |\n"
        << "|                 Unit Test (CubaturePolylib)                                 |\n"
        << "|                                                                             |\n"
        << "|     1) Computing integrals of monomials on reference cells in 1D            |\n"
        << "|                                                                             |\n"
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n"
        << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n"
        << "|                      Kyungjoo Kim (kyukim@sandia.gov).                      |\n"
        << "|                                                                             |\n"
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n"
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n"
        << "|                                                                             |\n"
        << "===============================================================================\n"
        << "| TEST 1: integrals of monomials in 1D                                        |\n"
        << "===============================================================================\n";

      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

      typedef ValueType pointValueType;
      typedef ValueType weightValueType;
      typedef CubaturePolylib<DeviceSpaceType,pointValueType,weightValueType> CubatureLineType;

      const auto tol = 100.0 * tolerence(); // relaxing this as we support higher-order cubature

      int errorFlag = 0;

      // open file with analytic values
      std::string basedir = "./data";
      std::stringstream namestream;
      std::string filename;
      namestream <<  basedir << "/EDGE_integrals" << ".dat";
      namestream >> filename;
      *outStream << "filename = "  << filename << std::endl;
      std::ifstream filecompare(filename);

      *outStream << "\n-> Integrals of monomials on a reference line (edge):\n";

      // compute and compare integrals
      try {
        const auto maxDeg   = Parameters::MaxCubatureDegreeEdge;
        const auto polySize = maxDeg + 1;

        // test inegral values
        DynRankView ConstructWithLabel(testInt, maxDeg+1, polySize);
        
        // analytic integral values
        DynRankView ConstructWithLabel(analyticInt, maxDeg+1, polySize);

        // get analytic values
        if (filecompare.is_open()) {
          getAnalytic(analyticInt, filecompare);
          filecompare.close();
        } else {
          INTREPID2_TEST_FOR_EXCEPTION( true, std::runtime_error,
                                        ">>> ERROR (Integration::Test02): Cannot open analytic solution file" );
        }
        
        // compute integrals
        EPolyType polyType[4] = { POLYTYPE_GAUSS,
                                  POLYTYPE_GAUSS_RADAU_LEFT,
                                  POLYTYPE_GAUSS_RADAU_RIGHT,
                                  POLYTYPE_GAUSS_LOBATTO };
        for (ordinal_type pid=0;pid<4;++pid) {
          const auto ptype = polyType[pid];
          *outStream << "\n  -> Testing poly type " << EPolyTypeToString(ptype) << "\n";

          for (ordinal_type cubDeg=0;cubDeg<=maxDeg;++cubDeg) {
            CubatureLineType lineCub(cubDeg, ptype);
            auto cubPoints  = lineCub.allocateCubaturePoints();
            auto cubWeights = lineCub.allocateCubatureWeights();
            lineCub.getCubature(cubPoints,cubWeights);
            Kokkos::Array<int,1> degrees;
            for (auto polyDeg=0;polyDeg<=cubDeg;++polyDeg)
            {
              degrees[0] = polyDeg;
              testInt(cubDeg, polyDeg) = computeIntegralOfMonomial<ValueType>(cubPoints,
                                                                              cubWeights,
                                                                              degrees);
            }
          }
                    
          // perform comparison
          for (ordinal_type cubDeg=0;cubDeg<=maxDeg;++cubDeg) {
            for (ordinal_type polyDeg=0;polyDeg<=cubDeg;++polyDeg) {
              const auto abstol  = ( analyticInt(polyDeg,0) == 0 ? tol : std::fabs(tol*analyticInt(polyDeg,0)) );
              const auto absdiff = std::fabs(analyticInt(polyDeg,0) - testInt(cubDeg,polyDeg));
              *outStream << "Cubature order " << std::setw(2) << std::left << cubDeg << " integrating "
                         << "x^" << std::setw(2) << std::left << polyDeg <<  ":" << "   "
                         << std::scientific << std::setprecision(16) << testInt(cubDeg,polyDeg) << "   " << analyticInt(polyDeg,0) << "   "
                         << std::setprecision(4) << absdiff << "   " << "<?" << "   " << abstol << "\n";
              if (absdiff > abstol) {
                errorFlag++;
                *outStream << std::right << std::setw(104) << "^^^^---FAILURE!\n";
              }
            }
            *outStream << "\n";
          }
        } 
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







