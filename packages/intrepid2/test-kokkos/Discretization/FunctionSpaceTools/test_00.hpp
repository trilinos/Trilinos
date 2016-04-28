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

/** \file test_00.cpp
    \brief  Unit test for the FunctionSpaceTools class.
    \author Created by D. Ridzal, P. Bochev, K. Peterson and Kyungjoo Kim
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

namespace Intrepid2 {
  
  namespace Test {
#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    try {                                                               \
      S ;                                                               \
    }                                                                   \
    catch (std::logic_error err) {                                      \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };                                                                  
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)
    
    template<typename ValueType, typename DeviceSpaceType>
    int FunctionSpaceTools_Test00(const bool verbose) {
      
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

      *outStream                                                        \
        << "===============================================================================\n" \
        << "|                                                                             |\n" \
        << "|                      Unit Test (FunctionSpaceTools)                         |\n" \
        << "|                                                                             |\n" \
        << "|     1) exception tests for general functionalities                          |\n" \
        << "|                                                                             |\n" \
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
        << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
        << "|                      Kyungjoo Kim (kyukim@sandia.gov).                      |\n" \
        << "|                                                                             |\n" \
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
        << "|                                                                             |\n" \
        << "===============================================================================\n";

      typedef FunctionSpaceTools<DeviceSpaceType> fst;
      typedef Kokkos::DynRankView<ValueType,DeviceSpaceType> DynRankView;
      
      int errorFlag = 0;

      *outStream                                \
        << "\n"
        << "===============================================================================\n" \
        << "| TEST 1: exceptions                                                          |\n" \
        << "===============================================================================\n";
      
      try{
#ifdef HAVE_INTREPID2_DEBUG
        DynRankView ConstructWithLabel( a_2, 2 );
        DynRankView ConstructWithLabel( a_2_2, 2, 2 );
        DynRankView ConstructWithLabel( a_2_3, 2, 3 );
        DynRankView ConstructWithLabel( a_3_2, 3, 2 );
        DynRankView ConstructWithLabel( a_2_2_3, 2, 2, 3 );
        DynRankView ConstructWithLabel( a_2_2_3_3, 2, 2, 3, 3 );
        DynRankView ConstructWithLabel( a_2_2_2, 2, 2, 2 );
        DynRankView ConstructWithLabel( a_2_2_2_3_3, 2, 2, 2, 3, 3 );
        DynRankView ConstructWithLabel( a_2_2_2_2_2, 2, 2, 2, 2, 2 );
        DynRankView ConstructWithLabel( a_2_2_2_2, 2, 2, 2, 2 );
        DynRankView ConstructWithLabel( a_3_2_2_2, 3, 2, 2, 2 );
        DynRankView ConstructWithLabel( a_2_3_2_2, 2, 3, 2, 2 );
        DynRankView ConstructWithLabel( a_2_2_3_2, 2, 2, 3, 2 );
        DynRankView ConstructWithLabel( a_2_2_2_3, 2, 2, 2, 3 );

        // scratch space
        DynRankView ConstructWithLabel( scratch, 256 );

        // hex is a place holder to fill the input argument (not specifically testing only for hex)
        const auto cellTopo = shards::getCellTopologyData< shards::Hexahedron<> >();
        
        *outStream << "-> Errors are expected due to rank mismatches:\n";
        
        *outStream << "-> computeCellMeasure:\n";
        
        INTREPID2_TEST_ERROR_EXPECTED( fst::computeCellMeasure(a_2_2, a_2, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::computeCellMeasure(a_2_2, a_2_2, a_2) );
        
        *outStream << "-> computeFaceMeasure:\n";
        
        INTREPID2_TEST_ERROR_EXPECTED( fst::computeFaceMeasure(a_2_2, a_2, a_2, 0, cellTopo, scratch) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::computeFaceMeasure(a_2_2, a_2_2_3_3, a_2, 0, cellTopo, scratch) );
        
        *outStream << "-> computeEdgeMeasure:\n";

        INTREPID2_TEST_ERROR_EXPECTED( fst::computeEdgeMeasure(a_2_2, a_2, a_2, 0, cellTopo, scratch) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::computeEdgeMeasure(a_2_2, a_2_2_2_2, a_2, 0, cellTopo, scratch) );

        *outStream << "-> integrate:\n";

        INTREPID2_TEST_ERROR_EXPECTED( fst::integrate(a_2_2_2_2, a_2_2_2, a_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::integrate(a_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::integrate(a_2, a_2_2_3, a_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::integrate(a_2, a_2_2_3_3, a_2_2_3_3) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::integrate(a_2_2, a_2_2, a_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::integrate(a_2_2, a_2_2_3, a_2_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::integrate(a_2_2, a_2_2_3_3, a_2_2_2_3_3) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::integrate(a_2_2_2, a_2_2_2, a_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::integrate(a_2_2_2, a_2_2_2_3, a_2_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::integrate(a_2_2_2, a_2_2_2_3_3, a_2_2_2_3_3) );
        
        *outStream << "-> applyLeftFieldSigns:\n";

        INTREPID2_TEST_ERROR_EXPECTED( fst::applyLeftFieldSigns(a_2, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::applyLeftFieldSigns(a_2_2_2, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::applyLeftFieldSigns(a_2_2_2, a_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::applyLeftFieldSigns(a_2_2_2, a_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::applyLeftFieldSigns(a_2_2_2, a_2_2) );
        
        *outStream << "-> applyRightFieldSigns:\n";

        INTREPID2_TEST_ERROR_EXPECTED( fst::applyRightFieldSigns(a_2, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::applyRightFieldSigns(a_2_2_2, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::applyRightFieldSigns(a_2_2_2, a_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::applyRightFieldSigns(a_2_2_2, a_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::applyRightFieldSigns(a_2_2_2, a_2_2) );

        *outStream << "-> applyFieldSigns:\n";

        INTREPID2_TEST_ERROR_EXPECTED( fst::applyFieldSigns(a_2, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::applyFieldSigns(a_2_2, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::applyFieldSigns(a_2_2, a_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::applyFieldSigns(a_2_2, a_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::applyFieldSigns(a_2_2_2_3_3, a_2_2) );

        *outStream << "-> evaluate:\n";

        INTREPID2_TEST_ERROR_EXPECTED( fst::evaluate(a_2, a_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::evaluate(a_2, a_2, a_2_2_2_3_3) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::evaluate(a_2, a_2_2, a_2_2_2_3_3) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::evaluate(a_2_2_3_3, a_3_2, a_2_2_2_3_3) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::evaluate(a_2_2_3_3, a_2_3, a_2_2_2_3_3) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::evaluate(a_3_2_2_2, a_2_2, a_2_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::evaluate(a_2_3_2_2, a_2_2, a_2_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::evaluate(a_2_2_3_2, a_2_2, a_2_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::evaluate(a_2_2_2_3, a_2_2, a_2_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( fst::evaluate(a_2_2_2_2, a_2_2, a_2_2_2_2_2) );
#endif
      } catch (std::logic_error err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };

      if (errorFlag != 0)
        std::cout << "End Result: TEST FAILED\n";
      else
        std::cout << "End Result: TEST PASSED\n";
      
      // reset format state of std::cout
      std::cout.copyfmt(oldFormatState);
      
      return errorFlag;
    }
  } // end of namespace TEST
} // end of namespace Intrepid2
