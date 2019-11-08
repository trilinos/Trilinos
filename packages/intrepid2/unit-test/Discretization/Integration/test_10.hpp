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

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

      *outStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);
 
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

      const auto tol = 10.0 * tolerence();

      int errorFlag = 0;

      // open file with analytic values
      std::string basedir = "../data";
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
        
        // storage for cubatrue points and weights
        DynRankView ConstructWithLabel(cubPoints, 
                                       Parameters::MaxIntegrationPoints, 
                                       Parameters::MaxDimension);

        DynRankView ConstructWithLabel(cubWeights, 
                                       Parameters::MaxIntegrationPoints);
        
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
            for (auto polyDeg=0;polyDeg<=cubDeg;++polyDeg) 
              testInt(cubDeg, polyDeg) = computeIntegralOfMonomial<ValueType>(lineCub,
                                                                              cubPoints,
                                                                              cubWeights,
                                                                              polyDeg);
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







