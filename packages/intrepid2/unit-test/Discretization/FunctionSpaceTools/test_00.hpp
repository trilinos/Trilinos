// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    catch (std::logic_error &err) {                                      \
      *outStream << "Expected Error ----------------------------------------------------------------\n"; \
      *outStream << err.what() << '\n';                                 \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    };                                                                  
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)
    
    template<typename ValueType, typename DeviceType>
    int FunctionSpaceTools_Test00(const bool verbose) {
      using ExecSpaceType = typename DeviceType::execution_space;
      
      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose)
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs, false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef typename
        Kokkos::DefaultHostExecutionSpace HostSpaceType ;

      *outStream << "DeviceSpace::  ";   ExecSpaceType().print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType().print_configuration(*outStream, false);

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

      int errorFlag = 0;

      *outStream                                \
        << "\n"
        << "===============================================================================\n" \
        << "| TEST 1: exceptions                                                          |\n" \
        << "===============================================================================\n";
      
      try{
#ifdef HAVE_INTREPID2_DEBUG
        typedef FunctionSpaceTools<DeviceType> fst;
        typedef Kokkos::DynRankView<ValueType,DeviceType> DynRankView;
        #define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

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
      } catch (std::logic_error &err) {
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
