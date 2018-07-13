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

/** \file test_04.hpp
\brief  Unit test for the dot multiply operations of the ArrayTools class.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_RealSpaceTools.hpp"
#include "Intrepid2_ArrayTools.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_GlobalMPISession.hpp"

namespace Intrepid2 {

  namespace Test {

#define INTREPID2_TEST_ERROR_EXPECTED( S )                              \
    {                                                                   \
      try {                                                             \
        S ;                                                             \
      }                                                                 \
      catch (std::logic_error err) {                                    \
        *outStream << "Expected Error ----------------------------------------------------------------\n"; \
        *outStream << err.what() << '\n';                               \
        *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
      };                                                                \
    }
    
    template<typename ValueType, typename DeviceSpaceType>
    int ArrayTools_Test04(const bool verbose) {
      
      typedef ValueType value_type;

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose) 
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      typedef typename
        Kokkos::Impl::is_space<DeviceSpaceType>::host_mirror_space::execution_space HostSpaceType ;

      *outStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(std::cout, false);
      *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(std::cout, false);

      *outStream                                                        \
        << "===============================================================================\n" \
        << "|                                                                             |\n" \
        << "|                       Unit Test (ArrayTools)                                |\n" \
        << "|                                                                             |\n" \
        << "|     1) Array operations: multiplication, contraction                        |\n" \
        << "|                                                                             |\n" \
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
        << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
        << "|                                                                             |\n" \
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
        << "|                                                                             |\n" \
        << "===============================================================================\n";      
      
      typedef RealSpaceTools<DeviceSpaceType> rst;
      typedef ArrayTools<DeviceSpaceType> art; 
      typedef Kokkos::DynRankView<value_type,DeviceSpaceType> DynRankView;

#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)
      
      const value_type tol = tolerence()*10000.0;
      int errorFlag = 0;
      
      
      /************************************************************************************************
        *                                                                                             *
        *  Exception tests: should only run when compiled in DEBUG mode                               *
        *                                                                                             *
        ***********************************************************************************************/
      int C  = 12,  C1 = 10;
      int P  = 8,   P1 = 16;
      int F  = 4,   F1 = 8;
      int D1 = 1,   D2 = 2,   D3 = 3,   D4 = 4;
  
      DynRankView ConstructWithLabel(fc_P, P);
      DynRankView ConstructWithLabel(fc_P_D1, P, D1);
      DynRankView ConstructWithLabel(fc_P_D2, P, D2);
      DynRankView ConstructWithLabel(fc_P_D3, P, D3);
      DynRankView ConstructWithLabel(fc_P_D2_D2, P, D2, D2);
      DynRankView ConstructWithLabel(fc_P1_D3, P1, D3);
      DynRankView ConstructWithLabel(fc_P1_D2_D2, P1, D2, D2);
      DynRankView ConstructWithLabel(fc_P1_D3_D3, P1, D3, D3);

      DynRankView ConstructWithLabel(fc_C, C);
      DynRankView ConstructWithLabel(fc_C1_P, C1,P);
      DynRankView ConstructWithLabel(fc_C_P1, C, P1);
      DynRankView ConstructWithLabel(fc_C_P, C, P);
      DynRankView ConstructWithLabel(fc_C_P_D1, C, P, D1);
      DynRankView ConstructWithLabel(fc_C_P_D2, C, P, D2);
      DynRankView ConstructWithLabel(fc_C_P_D3, C, P, D3);
      DynRankView ConstructWithLabel(fc_C_P_D4, C, P, D4);
      DynRankView ConstructWithLabel(fc_C1_P_D2, C1,P, D2);
      DynRankView ConstructWithLabel(fc_C1_P_D3, C1,P, D3);
      DynRankView ConstructWithLabel(fc_C_P1_D1, C, P1,D1);
      DynRankView ConstructWithLabel(fc_C_P1_D2, C, P1,D2);
      DynRankView ConstructWithLabel(fc_C_P1_D3, C, P1,D3);
      DynRankView ConstructWithLabel(fc_C_P_D1_D1, C, P, D1, D1);
      DynRankView ConstructWithLabel(fc_C_P_D1_D2, C, P, D1, D2);
      DynRankView ConstructWithLabel(fc_C_P_D1_D3, C, P, D1, D3);
      DynRankView ConstructWithLabel(fc_C_P_D2_D2, C, P, D2, D2);
      DynRankView ConstructWithLabel(fc_C_P_D3_D1, C, P, D3, D1);
      DynRankView ConstructWithLabel(fc_C_P_D3_D2, C, P, D3, D2);
      DynRankView ConstructWithLabel(fc_C_P_D2_D3, C, P, D2, D3);
      DynRankView ConstructWithLabel(fc_C_P_D3_D3, C, P, D3, D3);
      DynRankView ConstructWithLabel(fc_C1_P_D3_D3, C1,P, D3, D3);
      DynRankView ConstructWithLabel(fc_C1_P_D2_D2, C1,P, D2, D2);
      DynRankView ConstructWithLabel(fc_C_P1_D2_D2, C, P1,D2, D2);
      DynRankView ConstructWithLabel(fc_C_P1_D3_D3, C, P1,D3, D3);
      DynRankView ConstructWithLabel(fc_C_P_D3_D3_D3, C, P, D3, D3, D3);

      DynRankView ConstructWithLabel(fc_F_P, F, P);
      DynRankView ConstructWithLabel(fc_F_P_D1, F, P, D1);
      DynRankView ConstructWithLabel(fc_F_P_D2, F, P, D2);
      DynRankView ConstructWithLabel(fc_F_P1_D2, F, P1, D2);
      DynRankView ConstructWithLabel(fc_F_P_D3, F, P, D3);
      DynRankView ConstructWithLabel(fc_F_P_D3_D3, F, P, D3, D3);
      DynRankView ConstructWithLabel(fc_F1_P_D2, F1,P, D2);
      DynRankView ConstructWithLabel(fc_F1_P_D3, F1,P, D3);
      DynRankView ConstructWithLabel(fc_F1_P_D3_D3, F1,P, D3, D3);
      DynRankView ConstructWithLabel(fc_F_P1_D3, F, P1,D3);
      DynRankView ConstructWithLabel(fc_F_P1_D3_D3, F, P1,D3, D3);
      DynRankView ConstructWithLabel(fc_F_P_D2_D2, F, P, D2, D2);
      DynRankView ConstructWithLabel(fc_F_P1_D2_D2, F, P1,D2, D2);
      DynRankView ConstructWithLabel(fc_C_F_P, C, F, P);
      DynRankView ConstructWithLabel(fc_C1_F_P, C1, F, P);
      DynRankView ConstructWithLabel(fc_C_F1_P, C, F1,P);
      DynRankView ConstructWithLabel(fc_C_F_P1, C, F, P1);
      DynRankView ConstructWithLabel(fc_C_F_P_D1, C, F, P, D1);
      DynRankView ConstructWithLabel(fc_C_F_P_D2, C, F, P, D2);
      DynRankView ConstructWithLabel(fc_C_F_P_D3, C, F, P, D3);
      DynRankView ConstructWithLabel(fc_C1_F_P_D2, C1, F, P,D2);
      DynRankView ConstructWithLabel(fc_C1_F_P_D3, C1, F, P,D3);
      DynRankView ConstructWithLabel(fc_C_F1_P_D2, C, F1,P, D2);
      DynRankView ConstructWithLabel(fc_C_F1_P_D3, C, F1,P, D3);
      DynRankView ConstructWithLabel(fc_C_F1_P_D3_D3, C, F1,P, D3, D3);
      DynRankView ConstructWithLabel(fc_C_F_P1_D2, C, F, P1,D2);
      DynRankView ConstructWithLabel(fc_C_F_P1_D3, C, F, P1,D3);
      DynRankView ConstructWithLabel(fc_C_F_P_D1_D1, C, F, P, D1, D1);
      DynRankView ConstructWithLabel(fc_C_F_P_D2_D2, C, F, P, D2, D2);
      DynRankView ConstructWithLabel(fc_C_F_P_D1_D3, C, F, P, D1, D3);
      DynRankView ConstructWithLabel(fc_C_F_P_D2_D3, C, F, P, D2, D3);
      DynRankView ConstructWithLabel(fc_C_F_P_D3_D1, C, F, P, D3, D1);
      DynRankView ConstructWithLabel(fc_C_F_P_D3_D2, C, F, P, D3, D2);
      DynRankView ConstructWithLabel(fc_C_F_P_D3_D3, C, F, P, D3, D3);
      DynRankView ConstructWithLabel(fc_C_F_P1_D2_D2, C, F, P1,D2, D2);
      DynRankView ConstructWithLabel(fc_C_F_P1_D3_D3, C, F, P1,D3, D3);
      DynRankView ConstructWithLabel(fc_C1_F_P_D2_D2, C1,F, P, D2, D2);
      DynRankView ConstructWithLabel(fc_C1_F_P1_D2_D2, C1,F, P1,D2, D2);
      DynRankView ConstructWithLabel(fc_C1_F_P_D3_D3, C1,F, P, D3, D3);
      DynRankView ConstructWithLabel(fc_C_F_P_D3_D3_D3, C, F, P, D3, D3, D3);
  
  
      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 1: crossProductDataField exceptions                                    |\n"\
        << "===============================================================================\n";
    #ifdef HAVE_INTREPID2_DEBUG
      try{
        // 39 exceptions
        // Test rank and D dimension of inputData
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3, fc_C_P,       fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3, fc_C_P_D2_D2, fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3, fc_C_P_D1,    fc_C_F_P_D3) );

        // Test rank and D dimension of inputFields
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3,  fc_C_P_D3,  fc_F_P) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3,  fc_C_P_D3,  fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3,  fc_C_P_D3,  fc_C_F_P_D1) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3,  fc_C_P_D3,  fc_F_P_D1) );

        // Test rank of outputFields
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3,  fc_C_P_D2,  fc_C_F_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P,     fc_C_P_D3,  fc_C_F_P_D3) );
        
        // Dimension cross-check: (1) inputData    vs. inputFields
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3,  fc_C_P_D3, fc_C1_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3,  fc_C_P_D3, fc_C_F_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3,  fc_C_P_D3, fc_C_F_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P,     fc_C_P_D2, fc_C1_F_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P,     fc_C_P_D2, fc_C_F_P1_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P,     fc_C_P_D2, fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3,  fc_C_P_D3, fc_F_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3,  fc_C_P_D3, fc_F_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P,     fc_C_P_D2, fc_F_P1_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P,     fc_C_P_D2, fc_F_P_D3) );

        // Dimension cross-check: (2) outputFields vs. inputData
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C1_F_P,  fc_C_P_D2, fc_C_F_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P1,  fc_C_P_D2, fc_C_F_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C1_F_P,  fc_C_P_D2, fc_F_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P1,  fc_C_P_D2, fc_F_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C1_F_P_D3, fc_C_P_D3, fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P1_D3, fc_C_P_D3, fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D2,  fc_C_P_D3, fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C1_F_P_D3, fc_C_P_D3, fc_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P1_D3, fc_C_P_D3, fc_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D2,  fc_C_P_D3, fc_F_P_D3) );

         // Dimension cross-check: (3) outputFields vs. inputFields
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P,  fc_C1_P_D2, fc_C1_F_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P,  fc_C_P1_D2, fc_C_F_P1_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P,  fc_C_P_D2, fc_C_F1_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P,  fc_C_P1_D2, fc_F_P1_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P,  fc_C_P_D2, fc_F1_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3, fc_C1_P_D3,  fc_C1_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3, fc_C_P1_D3,  fc_C_F_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3, fc_C_P_D3,   fc_C_F1_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3, fc_C_P1_D3,  fc_F_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataField(fc_C_F_P_D3, fc_C_P_D3,   fc_F1_P_D3) );

        *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 2: crossProductDataData exceptions                                     |\n"\
        << "===============================================================================\n";

        // 18 exceptions
        // inputDataL is (C, P, D) and 2 <= D <= 3 is required
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P_D3, fc_C_P,       fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P_D3, fc_C_P_D2_D2, fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P_D3, fc_C_P_D1,    fc_C_P_D3) );

        // inputDataRight is (C, P, D) or (P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required.
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P_D3, fc_C_P_D3,    fc_C) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P_D3, fc_C_P_D3,    fc_C_P_D2_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P_D3, fc_C_P_D3,    fc_C_P_D1) );

        // outputData is (C,P,D) in 3D and (C,P) in 2D => rank = inputDataLeft.extent(2)
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P_D3, fc_C_P_D2,    fc_C_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P,    fc_C_P_D3,    fc_C_P_D3) );

        // Dimension cross-check (1):
        // inputDataLeft(C,P,D) vs. inputDataRight(C,P,D): C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P_D3, fc_C1_P_D3,   fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P_D3, fc_C_P1_D3,   fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P_D3, fc_C_P_D3,    fc_C_P_D2) );
        // inputDataLeft(C, P,D) vs. inputDataRight(P,D):  P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P_D3, fc_C_P1_D3,   fc_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P_D3, fc_C_P_D3,    fc_P_D2) );

        // Dimension cross-check (2):
        // in 2D: outputData(C,P) vs. inputDataLeft(C,P,D): dimensions C, P must match
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C1_P,    fc_C_P_D2,   fc_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P1,    fc_C_P_D2,   fc_P_D2) );
        // in 3D: outputData(C,P,D) vs. inputDataLeft(C,P,D): all dimensions C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C1_P_D3, fc_C_P_D3,   fc_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P1_D3, fc_C_P_D3,   fc_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::crossProductDataData(fc_C_P_D2,  fc_C_P_D3,   fc_P_D3) );

        *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 3: outerProductDataField exceptions                                    |\n"\
        << "===============================================================================\n";
        // 28 exceptions
        // Test rank and D dimension: inputData(C, P, D) and 2 <= D <= 3 is required
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P,       fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P_D2_D2, fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P_D1,    fc_C_F_P_D3) );

        // Test rank and D dimension: inputFields(C,F,P,D)/(F,P,D) and 2 <= (D=dimension(rank - 1)) <= 3 is required.
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3,    fc_F_P) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3,    fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3,    fc_C_F_P_D1) );

        //  Test rank and D dimension: outputFields(C,F,P,D,D)
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P,         fc_C_P_D3,    fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3,      fc_C_P_D3,    fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3_D3,fc_C_P_D3,    fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D1_D3,   fc_C_P_D3,    fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D1,   fc_C_P_D3,    fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D1_D1,   fc_C_P_D3,    fc_C_F_P_D3) );

        // Cross-check (2): outputFields(C,F,P,D,D) vs. inputData(C,P,D): dimensions C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C1_P_D3,   fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P1_D3,   fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P_D2,    fc_C_F_P_D2) );

        // Cross-check (1): inputData(C,P,D) vs. inputFields(C,F,P,D):  dimensions  C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3,   fc_C1_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3,   fc_C_F_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3,   fc_C_F_P_D2) );
        // Cross-check (1): inputData(C,P,D) vs. inputFields(F,P,D): dimensions  P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3,   fc_F_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3,   fc_F_P_D2) );

        // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(C,F,P,D): dimensions C, F, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C1_P_D3,   fc_C1_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3,    fc_C_F1_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P1_D3,   fc_C_F_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P_D2,    fc_C_F_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D2_D3, fc_C_P_D2,    fc_C_F_P_D2) );
        // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(F,P,D): dimensions F, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3,    fc_F1_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P1_D3,   fc_F_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataField(fc_C_F_P_D3_D3, fc_C_P_D2,    fc_F_P_D2) );
        *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 4: outerProductDataData exceptions                                     |\n"\
        << "===============================================================================\n";
        // 26 exceptions
        // (1) inputDataLeft is (C, P, D) and 2 <= D <= 3 is required
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P,         fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,   fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P_D1,      fc_C_P_D3) );

        // (2) inputDataRight is (C, P, D) or (P, D) and 2 <= (D=dimension(rank - 1)) <= 3 is required.
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P_D3,      fc_C) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P_D3,      fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P_D3,      fc_C_P_D1) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P_D3,      fc_P_D1) );

        // (3) outputData is (C,P,D,D)
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3,      fc_C_P_D3,      fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3_D3,fc_C_P_D3,      fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D1,   fc_C_P_D3,      fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D1_D2,   fc_C_P_D3,      fc_C_P_D3) );

        // Cross-check (2): outputData(C,P,D,D) vs. inputDataLeft(C,P,D): dimensions C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C1_P_D3,       fc_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P1_D3,       fc_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P_D2,        fc_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D2_D3, fc_C_P_D2,        fc_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D2, fc_C_P_D2,        fc_P_D3) );

        // Cross-check (1): inputDataLeft(C,P,D) vs. inputDataRight(C,P,D):  all dimensions  C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P_D3,       fc_C1_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P_D3,       fc_C_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P_D3,       fc_C_P_D2) );
        // Cross-check (1): inputDataLeft(C,P,D) vs. inputDataRight(P,D): dimensions  P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P_D3,       fc_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P_D3,       fc_P_D2) );

        // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(C,P,D): dimensions C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C1_P_D3,       fc_C1_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P1_D3,       fc_C_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P_D2,        fc_C_P_D2) );
        // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(P,D): dimensions P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P1_D3,      fc_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::outerProductDataData(fc_C_P_D3_D3, fc_C_P_D2,       fc_P_D2) );

        *outStream \
         << "\n"
         << "===============================================================================\n"\
         << "| TEST 5: matvecProductDataField exceptions                                   |\n"\
         << "===============================================================================\n";
        // 34 exceptions
        // (1) inputData is (C,P), (C,P,D) or (C, P, D, D) and 1 <= D <= 3 is required
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C,              fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D4,         fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3_D3_D3,   fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3_D1,      fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D1_D3,      fc_C_F_P_D3) );

        // (2) inputFields is (C, F, P, D) or (F, P, D) and 1 <= (D=dimension(rank - 1)) <= 3 is required.
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3_D3,    fc_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3_D3,    fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3_D3,    fc_C_F_P_D1) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3_D3,    fc_F_P_D1) );
        // (3) outputFields is (C,F,P,D)
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,    fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P,       fc_C_P_D3_D3,    fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D1,    fc_C_P_D3_D3,    fc_C_F_P_D3) );

        // Cross-check (2): outputFields(C,F,P,D) vs. inputData(C,P,D) and (C,P,D,D): dimensions C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C1_P_D3_D3,    fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P1_D3_D3,    fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D2_D2,     fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C1_P,          fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P1,          fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C1_P_D3,       fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P1_D3,       fc_C_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D2,        fc_C_F_P_D3) );

        // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(C,F,P,D):  dimensions  C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3_D3,    fc_C1_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3_D3,    fc_C_F_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3_D3,    fc_C_F_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3,       fc_C_F_P_D2) );

        // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(F,P,D): dimensions  P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3_D3,    fc_F_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3_D3,    fc_F_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3,       fc_F_P_D2) );

        // Cross-check (3): outputFields(C,F,P,D) vs. inputFields(C,F,P,D): all dimensions C, F, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3_D3,    fc_C1_F_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3_D3,    fc_C_F1_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P1_D3_D3,   fc_C_F_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D2, fc_C_P_D2_D2,    fc_C_F_P_D3) );

        // Cross-check (3): outputFields(C,F,P,D) vs. inputFields(F,P,D): dimensions F, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3_D3,    fc_C_F1_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P1_D3_D3,   fc_C_F_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataField(fc_C_F_P_D3, fc_C_P_D3_D3,    fc_C_F_P_D2) );

        *outStream \
         << "\n"
         << "===============================================================================\n"\
         << "| TEST 6: matvecProductDataData exceptions                                    |\n"\
         << "===============================================================================\n";
        // 37 exceptions
        // (1) inputDataLeft is (C,P), (C,P,D) or (C,P,D,D) and 1 <= D <= 3 is required
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C,             fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3_D3_D3,  fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3_D1,     fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D1_D3,     fc_C_P_D3) );

        // (2) inputDataRight is (C, P, D) or (P, D) and 1 <= (D=dimension(rank - 1)) <= 3 is required.
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3_D3,    fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3_D3,    fc_P) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3_D3,    fc_P_D1) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3_D3,    fc_C_P_D1) );

        // (3) outputData is (C,P,D)
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P,       fc_C_P_D3_D3,    fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,    fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D1,    fc_C_P_D3_D3,    fc_C_P_D3) );

        // Cross-check (2): outputData(C,P,D) vs. inputDataLeft(C,P), (C,P,D), (C,P,D,D): dimensions C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C1_P_D3, fc_C_P_D3_D3,  fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P1_D3, fc_C_P_D3_D3,  fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D2,  fc_C_P_D3_D3,  fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C1_P_D3, fc_C_P,        fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P1_D3, fc_C_P,        fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C1_P_D3, fc_C_P_D3,     fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P1_D3, fc_C_P_D3,     fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D2,  fc_C_P_D3,     fc_C_P_D3) );

        // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(C,P,D):  dimensions  C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3_D3,    fc_C1_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3_D3,    fc_C_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3_D3,    fc_C_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P,          fc_C1_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P,          fc_C_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3,       fc_C1_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3,       fc_C_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3,       fc_C_P_D2) );

        // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(P,D): dimensions  P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3_D3,    fc_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3_D3,    fc_P_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P,          fc_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3,       fc_P1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3,       fc_P_D2) );

        // Cross-check (3): outputData(C,P,D) vs. inputDataRight(C,P,D): all dimensions C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C1_P_D3, fc_C1_P_D3_D3,    fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P1_D3, fc_C_P1_D3_D3,    fc_C_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D2,  fc_C_P_D3_D3,     fc_C_P_D3) );

        // Cross-check (3): outputData(C,P,D) vs. inputDataRight(P,D): dimensions P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P1_D3, fc_C_P1_D3_D3,    fc_P_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matvecProductDataData(fc_C_P_D3, fc_C_P_D3_D3,    fc_P_D2) );

        *outStream \
         << "\n"
         << "===============================================================================\n"\
         << "| TEST 7: matmatProductDataField exceptions                                   |\n"\
         << "===============================================================================\n";
        // 46 exceptions
        // (1) inputData is (C,P), (C,P,D), or (C,P,D,D) and 1 <= D <= 3 is required
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C,             fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3_D3,  fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D1_D3,     fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D1,     fc_C_F_P_D3_D3) );

        // (2) inputFields is (C,F,P,D,D) or (F,P,D,D) and 1 <= (dimension(rank-1), (rank-2)) <= 3 is required.
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_C_F_P) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_C_F_P_D3_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_C_F_P_D1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_C_F_P_D3_D1) );

        // (3) outputFields is (C,F,P,D,D)
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3,       fc_C_P_D3_D3,  fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3_D3, fc_C_P_D3_D3,  fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D1_D3,    fc_C_P_D3_D3,  fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D1,    fc_C_P_D3_D3,  fc_C_F_P_D3_D3) );

        // Cross-check (2): outputFields(C,F,P,D,D) vs. inputData(C,P), (C,P,D), or (C,P,D,D): dimensions C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C1_P_D3_D3,  fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P1_D3_D3,  fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D2_D2,   fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P1_D2_D2,  fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D2_D2, fc_C_P_D3_D3,   fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C1_P,  fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P1,  fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C1_P_D3,  fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P1_D3,  fc_C_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D1,   fc_C_F_P_D3_D3) );

        // Cross-check (1): inputData(C,P), (C,P,D) or (C,P,D,D) vs. inputFields(C,F,P,D,D):  dimensions  C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_C1_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_C_F_P1_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_C_F_P_D2_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_C_F_P1_D2_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_C1_F_P_D2_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_C1_F_P1_D2_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P,  fc_C1_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P,  fc_C_F_P1_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3,  fc_C1_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3,  fc_C_F_P1_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3,  fc_C_F_P_D2_D2) );

        // Cross-check (1): inputData(C,P), (C,P,D), or (C,P,D,D) vs. inputFields(F,P,D,D): dimensions  P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_F_P1_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_F_P_D2_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_F_P1_D2_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P,    fc_F_P1_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3, fc_F_P1_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3, fc_F_P_D2_D2) );

        // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(C,F,P,D,D): all dimensions C, F, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_C1_F_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_C_F1_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_C_F_P1_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_C_F_P_D2_D2) );

        // Cross-check (3): outputFields(C,F,P,D,D) vs. inputFields(F,P,D,D): dimensions F, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_F1_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_F_P1_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_F_P_D3_D3, fc_C_P_D3_D3,  fc_F_P_D2_D2) );
        *outStream \
         << "\n"
         << "===============================================================================\n"\
         << "| TEST 8: matmatProductDataData exceptions                                    |\n"\
         << "===============================================================================\n";
        // 45 exceptions
        // (1) inputDataLeft is (C,P), (C,P,D), or (C,P,D,D) and 2 <= D <= 3 is required
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C,             fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3_D3,  fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D1_D3,     fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D1,     fc_C_P_D3_D3) );

        // (2) inputDataRight is (C,P,D,D) or (P,D,D) and 1 <= (dimension(rank-1), (rank-2)) <= 3 is required.
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P,        fc_C_P) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3,     fc_C_P_D3_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_C_P_D1_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_C_P_D3_D1) );

        // (3) outputData is (C,P,D,D)
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3,       fc_C_P,        fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3_D3, fc_C_P_D3,     fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D1_D3,    fc_C_P_D3_D3,  fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D1,    fc_C_P_D3_D3,  fc_C_P_D3_D3) );

        // Cross-check (2): outputData(C,P,D,D) vs. inputDataLeft(C,P), (C,P,D), or (C,P,D,D): dimensions C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C1_P_D3_D3,  fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P1_D3_D3,  fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D2_D2,   fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P1_D2_D2,  fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D2_D2, fc_C_P_D3_D3,   fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C1_P,        fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P1,        fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C1_P_D3,     fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P1_D3,     fc_C_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D1,      fc_C_P_D3_D3) );

        // Cross-check (1): inputDataLeft(C,P), (C,P,D) or (C,P,D,D) vs. inputDataRight(C,P,D,D):  dimensions  C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_C1_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_C_P1_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_C_P_D2_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_C_P1_D2_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_C1_P_D2_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_C_P1_D2_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P,  fc_C1_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P,  fc_C_P1_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3,  fc_C1_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3,  fc_C_P1_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3,  fc_C_P_D2_D2) );

        // Cross-check (1): inputDataLeft(C,P), (C,P,D), or (C,P,D,D) vs. inputDataRight(P,D,D): dimensions  P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_P1_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_P_D2_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_P1_D2_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_P_D3_D3, fc_C_P,    fc_P1_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_P_D3_D3, fc_C_P_D3, fc_P1_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataField(fc_C_P_D3_D3, fc_C_P_D3, fc_P_D2_D2) );

        // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(C,P,D,D): all dimensions C, P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_C1_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_C1_P_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_C_P1_D3_D3) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_C_P_D2_D2) );

        // Cross-check (3): outputData(C,P,D,D) vs. inputDataRight(P,D,D): dimensions P, D must match
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_P_D2_D2) );
        INTREPID2_TEST_ERROR_EXPECTED(art::matmatProductDataData(fc_C_P_D3_D3, fc_C_P_D3_D3,  fc_P1_D3_D3) );
        }

        catch (std::logic_error err) {
          *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
          *outStream << err.what() << '\n';
          *outStream << "-------------------------------------------------------------------------------" << "\n\n";
          errorFlag = -1000;
        };

        #endif

      /************************************************************************************************
        *                                                                                             *
        *                                 Operation corectness tests                                  *
        *                                                                                             *
        ***********************************************************************************************/

      try{

      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 1.a: 3D crossProductDataField operations: (C,P,D) and (C,F,P,D)        |\n"\
        << "===============================================================================\n";
      /*
       *        ijkData_1a      ijkFields_1a       Expected result in (C,F,P,D) array
       *          i,i          i,i  j,j,  k,k           0, 0    k, k  -j,-j
       *          j,j          i,i  j,j,  k,k          -k,-k    0, 0   i, i
       *          k,k          i,i  j,j,  k,k           j, j   -i,-i   0, 0
       */

      // input data is (C,P,D)
      DynRankView ConstructWithLabel(ijkData_1a, 3, 2, 3);
      // C=0 contains i
      ijkData_1a(0, 0, 0) = 1.0;   ijkData_1a(0, 0, 1) = 0.0;   ijkData_1a(0, 0, 2) = 0.0;
      ijkData_1a(0, 1, 0) = 1.0;   ijkData_1a(0, 1, 1) = 0.0;   ijkData_1a(0, 1, 2) = 0.0;
      // C=1 contains j
      ijkData_1a(1, 0, 0) = 0.0;   ijkData_1a(1, 0, 1) = 1.0;   ijkData_1a(1, 0, 2) = 0.0;
      ijkData_1a(1, 1, 0) = 0.0;   ijkData_1a(1, 1, 1) = 1.0;   ijkData_1a(1, 1, 2) = 0.0;
      // C=2 contains k
      ijkData_1a(2, 0, 0) = 0.0;   ijkData_1a(2, 0, 1) = 0.0;   ijkData_1a(2, 0, 2) = 1.0;
      ijkData_1a(2, 1, 0) = 0.0;   ijkData_1a(2, 1, 1) = 0.0;   ijkData_1a(2, 1, 2) = 1.0;


      DynRankView ConstructWithLabel(ijkFields_1a,3, 3, 2, 3);
      // C=0, F=0 is i
      ijkFields_1a(0, 0, 0, 0) = 1.0; ijkFields_1a(0, 0, 0, 1) = 0.0; ijkFields_1a(0, 0, 0, 2) = 0.0;
      ijkFields_1a(0, 0, 1, 0) = 1.0; ijkFields_1a(0, 0, 1, 1) = 0.0; ijkFields_1a(0, 0, 1, 2) = 0.0;
      // C=0, F=1 is j
      ijkFields_1a(0, 1, 0, 0) = 0.0; ijkFields_1a(0, 1, 0, 1) = 1.0; ijkFields_1a(0, 1, 0, 2) = 0.0;
      ijkFields_1a(0, 1, 1, 0) = 0.0; ijkFields_1a(0, 1, 1, 1) = 1.0; ijkFields_1a(0, 1, 1, 2) = 0.0;
      // C=0, F=2 is k
      ijkFields_1a(0, 2, 0, 0) = 0.0; ijkFields_1a(0, 2, 0, 1) = 0.0; ijkFields_1a(0, 2, 0, 2) = 1.0;
      ijkFields_1a(0, 2, 1, 0) = 0.0; ijkFields_1a(0, 2, 1, 1) = 0.0; ijkFields_1a(0, 2, 1, 2) = 1.0;

      // C=1, F=0 is i
      ijkFields_1a(1, 0, 0, 0) = 1.0; ijkFields_1a(1, 0, 0, 1) = 0.0; ijkFields_1a(1, 0, 0, 2) = 0.0;
      ijkFields_1a(1, 0, 1, 0) = 1.0; ijkFields_1a(1, 0, 1, 1) = 0.0; ijkFields_1a(1, 0, 1, 2) = 0.0;
      // C=1, F=1 is j
      ijkFields_1a(1, 1, 0, 0) = 0.0; ijkFields_1a(1, 1, 0, 1) = 1.0; ijkFields_1a(1, 1, 0, 2) = 0.0;
      ijkFields_1a(1, 1, 1, 0) = 0.0; ijkFields_1a(1, 1, 1, 1) = 1.0; ijkFields_1a(1, 1, 1, 2) = 0.0;
      // C=1, F=2 is k
      ijkFields_1a(1, 2, 0, 0) = 0.0; ijkFields_1a(1, 2, 0, 1) = 0.0; ijkFields_1a(1, 2, 0, 2) = 1.0;
      ijkFields_1a(1, 2, 1, 0) = 0.0; ijkFields_1a(1, 2, 1, 1) = 0.0; ijkFields_1a(1, 2, 1, 2) = 1.0;

      // C=2, F=0 is i
      ijkFields_1a(2, 0, 0, 0) = 1.0; ijkFields_1a(2, 0, 0, 1) = 0.0; ijkFields_1a(2, 0, 0, 2) = 0.0;
      ijkFields_1a(2, 0, 1, 0) = 1.0; ijkFields_1a(2, 0, 1, 1) = 0.0; ijkFields_1a(2, 0, 1, 2) = 0.0;
      // C=2, F=1 is j
      ijkFields_1a(2, 1, 0, 0) = 0.0; ijkFields_1a(2, 1, 0, 1) = 1.0; ijkFields_1a(2, 1, 0, 2) = 0.0;
      ijkFields_1a(2, 1, 1, 0) = 0.0; ijkFields_1a(2, 1, 1, 1) = 1.0; ijkFields_1a(2, 1, 1, 2) = 0.0;
      // C=2, F=2 is k
      ijkFields_1a(2, 2, 0, 0) = 0.0; ijkFields_1a(2, 2, 0, 1) = 0.0; ijkFields_1a(2, 2, 0, 2) = 1.0;
      ijkFields_1a(2, 2, 1, 0) = 0.0; ijkFields_1a(2, 2, 1, 1) = 0.0; ijkFields_1a(2, 2, 1, 2) = 1.0;


      DynRankView ConstructWithLabel(outFields, 3, 3, 2, 3);
      art::crossProductDataField(outFields, ijkData_1a, ijkFields_1a);

      // checks for C = 0
      if( !(outFields(0,0,0,0)==0.0 && outFields(0,0,0,1)==0.0 && outFields(0,0,0,2)==0.0 &&
            outFields(0,0,1,0)==0.0 && outFields(0,0,1,1)==0.0 && outFields(0,0,1,2)==0.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (1): i x i != 0; ";
        errorFlag = -1000;
      }
      if( !(outFields(0,1,0,0)==0.0 && outFields(0,1,0,1)==0.0 && outFields(0,1,0,2)==1.0 &&
            outFields(0,1,1,0)==0.0 && outFields(0,1,1,1)==0.0 && outFields(0,1,1,2)==1.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (2): i x j != k; ";
        errorFlag = -1000;
      }
      if( !(outFields(0,2,0,0)==0.0 && outFields(0,2,0,1)==-1.0 && outFields(0,2,0,2)==0.0 &&
            outFields(0,2,1,0)==0.0 && outFields(0,2,1,1)==-1.0 && outFields(0,2,1,2)==0.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (3): i x k != -j; ";
        errorFlag = -1000;
      }

      // checks for C = 1
      if( !(outFields(1,0,0,0)==0.0 && outFields(1,0,0,1)==0.0 && outFields(1,0,0,2)==-1.0 &&
            outFields(1,0,1,0)==0.0 && outFields(1,0,1,1)==0.0 && outFields(1,0,1,2)==-1.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (4): j x i != -k; ";
        errorFlag = -1000;
      }
      if( !(outFields(1,1,0,0)==0.0 && outFields(1,1,0,1)==0.0 && outFields(1,1,0,2)==0.0 &&
            outFields(1,1,1,0)==0.0 && outFields(1,1,1,1)==0.0 && outFields(1,1,1,2)==0.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (5): j x j != 0; ";
        errorFlag = -1000;
      }
      if( !(outFields(1,2,0,0)==1.0 && outFields(1,2,0,1)==0.0 && outFields(1,2,0,2)==0.0 &&
            outFields(1,2,1,0)==1.0 && outFields(1,2,1,1)==0.0 && outFields(1,2,1,2)==0.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (6): j x k != i; ";
        errorFlag = -1000;
      }

      // checks for C = 2
      if( !(outFields(2,0,0,0)==0.0 && outFields(2,0,0,1)==1.0 && outFields(2,0,0,2)==0.0 &&
            outFields(2,0,1,0)==0.0 && outFields(2,0,1,1)==1.0 && outFields(2,0,1,2)==0.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (7): k x i != j; ";
        errorFlag = -1000;
      }
      if( !(outFields(2,1,0,0)==-1.0 && outFields(2,1,0,1)==0.0 && outFields(2,1,0,2)==0.0 &&
            outFields(2,1,1,0)==-1.0 && outFields(2,1,1,1)==0.0 && outFields(2,1,1,2)==0.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (8): k x j != -i; ";
        errorFlag = -1000;
      }
      if( !(outFields(2,2,0,0)==0.0 && outFields(2,2,0,1)==0.0 && outFields(2,2,0,2)==0.0 &&
            outFields(2,2,1,0)==0.0 && outFields(2,2,1,1)==0.0 && outFields(2,2,1,2)==0.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (9): k x k != 0; ";
        errorFlag = -1000;
      }

      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 1.b: 3D crossProductDataField operations:  (C,P,D) and (F,P,D)         |\n"\
        << "===============================================================================\n";
      /*
       *        ijkData_1b         ijkFields_1b   expected result in (C,F,P,D) array
       *          i,i,i               i,j,k                 0, k,-j
       *          j,j,j                                    -k, 0, i
       *          k,k,k                                     j,-i, 0
       */

      // input data is (C,P,D)
      DynRankView ConstructWithLabel(ijkData_1b, 3, 3, 3);
      // C=0 contains i
      ijkData_1b(0, 0, 0) = 1.0;   ijkData_1b(0, 0, 1) = 0.0;   ijkData_1b(0, 0, 2) = 0.0;
      ijkData_1b(0, 1, 0) = 1.0;   ijkData_1b(0, 1, 1) = 0.0;   ijkData_1b(0, 1, 2) = 0.0;
      ijkData_1b(0, 2, 0) = 1.0;   ijkData_1b(0, 2, 1) = 0.0;   ijkData_1b(0, 2, 2) = 0.0;
      // C=1 contains j
      ijkData_1b(1, 0, 0) = 0.0;   ijkData_1b(1, 0, 1) = 1.0;   ijkData_1b(1, 0, 2) = 0.0;
      ijkData_1b(1, 1, 0) = 0.0;   ijkData_1b(1, 1, 1) = 1.0;   ijkData_1b(1, 1, 2) = 0.0;
      ijkData_1b(1, 2, 0) = 0.0;   ijkData_1b(1, 2, 1) = 1.0;   ijkData_1b(1, 2, 2) = 0.0;
      // C=2 contains k
      ijkData_1b(2, 0, 0) = 0.0;   ijkData_1b(2, 0, 1) = 0.0;   ijkData_1b(2, 0, 2) = 1.0;
      ijkData_1b(2, 1, 0) = 0.0;   ijkData_1b(2, 1, 1) = 0.0;   ijkData_1b(2, 1, 2) = 1.0;
      ijkData_1b(2, 2, 0) = 0.0;   ijkData_1b(2, 2, 1) = 0.0;   ijkData_1b(2, 2, 2) = 1.0;

      // input fields are (F,P,D)
      DynRankView ConstructWithLabel(ijkFields_1b, 1, 3, 3);
      // F=0 at 3 points is (i,j,k)
      ijkFields_1b(0, 0, 0) = 1.0; ijkFields_1b(0, 0, 1) = 0.0; ijkFields_1b(0, 0, 2) = 0.0;
      ijkFields_1b(0, 1, 0) = 0.0; ijkFields_1b(0, 1, 1) = 1.0; ijkFields_1b(0, 1, 2) = 0.0;
      ijkFields_1b(0, 2, 0) = 0.0; ijkFields_1b(0, 2, 1) = 0.0; ijkFields_1b(0, 2, 2) = 1.0;

      // Output array is (C,F,P,D)
      outFields = DynRankView("outFields", 3, 1, 3, 3);
      art::crossProductDataField(outFields, ijkData_1b, ijkFields_1b);

      // checks for C = 0
      if( !(outFields(0,0,0,0)==0.0 && outFields(0,0,0,1)==0.0 && outFields(0,0,0,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (10): i x i != 0; ";
        errorFlag = -1000;
      }
      if( !(outFields(0,0,1,0)==0.0 && outFields(0,0,1,1)==0.0 && outFields(0,0,1,2)==1.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (11): i x j != k; ";
        errorFlag = -1000;
      }
      if( !(outFields(0,0,2,0)==0.0 && outFields(0,0,2,1)==-1.0 && outFields(0,0,2,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (12): i x k != -j; ";
        errorFlag = -1000;
      }

      // checks for C = 1
      if( !(outFields(1,0,0,0)==0.0 && outFields(1,0,0,1)==0.0 && outFields(1,0,0,2)==-1.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (13): j x i != -k; ";
        errorFlag = -1000;
      }
      if( !(outFields(1,0,1,0)==0.0 && outFields(1,0,1,1)==0.0 && outFields(1,0,1,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (14): j x j != 0; ";
        errorFlag = -1000;
      }
      if( !(outFields(1,0,2,0)==1.0 && outFields(1,0,2,1)==0.0 && outFields(1,0,2,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (15): j x k != i; ";
        errorFlag = -1000;
      }

      // checks for C = 2
      if( !(outFields(2,0,0,0)==0.0 && outFields(2,0,0,1)==1.0 && outFields(2,0,0,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (16): k x i != j; ";
        errorFlag = -1000;
      }
      if( !(outFields(2,0,1,0)==-1.0 && outFields(2,0,1,1)==0.0 && outFields(2,0,1,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (17): k x j != -i; ";
        errorFlag = -1000;
      }
      if( !(outFields(2,0,2,0)==0.0 && outFields(2,0,2,1)==0.0 && outFields(2,0,2,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (18): k x k != 0; ";
        errorFlag = -1000;
      }

      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 1.c: 2D crossProductDataField operations: (C,P,D) and (C,F,P,D)        |\n"\
        << "===============================================================================\n";

      /*
       *        ijData_1c        ijFields_1c    expected result in (C,F,P) array
       *          i,i               i,i  j,j             0, 0   1, 1
       *          j,j               i,i  j,j            -1,-1   0, 0
       */
      // input data is (C,P,D)
      DynRankView ConstructWithLabel(ijData_1c, 2, 2, 2);
      // C=0 contains i
      ijData_1c(0, 0, 0) = 1.0;   ijData_1c(0, 0, 1) = 0.0;
      ijData_1c(0, 1, 0) = 1.0;   ijData_1c(0, 1, 1) = 0.0;
      // C=1 contains j
      ijData_1c(1, 0, 0) = 0.0;   ijData_1c(1, 0, 1) = 1.0;
      ijData_1c(1, 1, 0) = 0.0;   ijData_1c(1, 1, 1) = 1.0;


      DynRankView ConstructWithLabel(ijFields_1c, 2, 2, 2, 2);
      // C=0, F=0 is i
      ijFields_1c(0, 0, 0, 0) = 1.0; ijFields_1c(0, 0, 0, 1) = 0.0;
      ijFields_1c(0, 0, 1, 0) = 1.0; ijFields_1c(0, 0, 1, 1) = 0.0;
      // C=0, F=1 is j
      ijFields_1c(0, 1, 0, 0) = 0.0; ijFields_1c(0, 1, 0, 1) = 1.0;
      ijFields_1c(0, 1, 1, 0) = 0.0; ijFields_1c(0, 1, 1, 1) = 1.0;

      // C=1, F=0 is i
      ijFields_1c(1, 0, 0, 0) = 1.0; ijFields_1c(1, 0, 0, 1) = 0.0;
      ijFields_1c(1, 0, 1, 0) = 1.0; ijFields_1c(1, 0, 1, 1) = 0.0;
      // C=1, F=1 is j
      ijFields_1c(1, 1, 0, 0) = 0.0; ijFields_1c(1, 1, 0, 1) = 1.0;
      ijFields_1c(1, 1, 1, 0) = 0.0; ijFields_1c(1, 1, 1, 1) = 1.0;

      // Output array is (C,F,P)
      outFields = DynRankView("outFields", 2, 2, 2);
      art::crossProductDataField(outFields, ijData_1c, ijFields_1c);

      if( !(outFields(0,0,0)==0.0 && outFields(0,0,1)==0.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (19): i x i != 0; ";
        errorFlag = -1000;
      }
      if( !(outFields(0,1,0)==1.0 && outFields(0,1,1)==1.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (20): i x j != 1; ";
        errorFlag = -1000;
      }

      if( !(outFields(1,0,0)==-1.0 && outFields(1,0,1)==-1.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (21): j x i != -1; ";
        errorFlag = -1000;
      }
      if( !(outFields(1,1,0)==0.0 && outFields(1,1,1)==0.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (22): j x j != 0; ";
        errorFlag = -1000;
      }

      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 1.d: 2D crossProductDataField operations: (C,P,D) and (F,P,D)          |\n"\
        << "===============================================================================\n";
      /*
       *        ijData_1c      ijFields_1d           expected result in (C,F,P) array
       *          i,i               i,j                      0, 1
       *          j,j                                       -1, 0
       */
      // inputFields is (F,P,D)
      DynRankView ConstructWithLabel(ijFields_1d, 1, 2, 2);
      // F=0 at 2 points is i,j
      ijFields_1d(0, 0, 0) = 1.0; ijFields_1d(0, 0, 1) = 0.0;
      ijFields_1d(0, 1, 0) = 0.0; ijFields_1d(0, 1, 1) = 1.0;

      // Output array is (C,F,P)
      outFields = DynRankView("outFields", 2, 1, 2);
      art::crossProductDataField(outFields, ijData_1c, ijFields_1d);

      if( !(outFields(0,0,0)==0.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (23): i x i != 0; ";
        errorFlag = -1000;
      }
      if( !(outFields(0,0,1)==1.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (24): i x j != 1; ";
        errorFlag = -1000;
      }
      if( !(outFields(1,0,0)==-1.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (25): j x i != -1; ";
        errorFlag = -1000;
      }
      if( !(outFields(1,0,1)==0.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataField (26): j x j != 0; ";
        errorFlag = -1000;
      }


      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 2.a: 3D crossProductDataData operations: (C,P,D) and (C,P,D)           |\n"\
        << "===============================================================================\n";
      /*
       *     ijkData_1a    jkiData_2a   kijData_2a
       *        i,i          j,j          k,k
       *        j,j          k,k          i,i
       *        k,k          i,i          j,j
       */
      DynRankView ConstructWithLabel(jkiData_2a, 3, 2, 3);
      // C=0 contains j
      jkiData_2a(0, 0, 0) = 0.0;   jkiData_2a(0, 0, 1) = 1.0;   jkiData_2a(0, 0, 2) = 0.0;
      jkiData_2a(0, 1, 0) = 0.0;   jkiData_2a(0, 1, 1) = 1.0;   jkiData_2a(0, 1, 2) = 0.0;
      // C=1 contains k
      jkiData_2a(1, 0, 0) = 0.0;   jkiData_2a(1, 0, 1) = 0.0;   jkiData_2a(1, 0, 2) = 1.0;
      jkiData_2a(1, 1, 0) = 0.0;   jkiData_2a(1, 1, 1) = 0.0;   jkiData_2a(1, 1, 2) = 1.0;
      // C=2 contains i
      jkiData_2a(2, 0, 0) = 1.0;   jkiData_2a(2, 0, 1) = 0.0;   jkiData_2a(2, 0, 2) = 0.0;
      jkiData_2a(2, 1, 0) = 1.0;   jkiData_2a(2, 1, 1) = 0.0;   jkiData_2a(2, 1, 2) = 0.0;

      DynRankView ConstructWithLabel(kijData_2a, 3, 2, 3);
      // C=0 contains k
      kijData_2a(0, 0, 0) = 0.0;   kijData_2a(0, 0, 1) = 0.0;   kijData_2a(0, 0, 2) = 1.0;
      kijData_2a(0, 1, 0) = 0.0;   kijData_2a(0, 1, 1) = 0.0;   kijData_2a(0, 1, 2) = 1.0;
      // C=1 contains i
      kijData_2a(1, 0, 0) = 1.0;   kijData_2a(1, 0, 1) = 0.0;   kijData_2a(1, 0, 2) = 0.0;
      kijData_2a(1, 1, 0) = 1.0;   kijData_2a(1, 1, 1) = 0.0;   kijData_2a(1, 1, 2) = 0.0;
      // C=2 contains j
      kijData_2a(2, 0, 0) = 0.0;   kijData_2a(2, 0, 1) = 1.0;   kijData_2a(2, 0, 2) = 0.0;
      kijData_2a(2, 1, 0) = 0.0;   kijData_2a(2, 1, 1) = 1.0;   kijData_2a(2, 1, 2) = 0.0;


      // ijkData_1a x ijkData_1a: outData should contain ixi=0, jxj=0, kxk=0
      DynRankView ConstructWithLabel(outData, 3,2,3);
      art::crossProductDataData(outData, ijkData_1a, ijkData_1a);

      for(int i = 0; i < 3; i++)
        for(int j = 0; j < 2; j++)
          for(int k = 0; k < 3; k++){
            if(outData(i,j,k) != 0) {
              *outStream << "\n\nINCORRECT crossProductDataData (1): i x i, j x j, or k x k != 0; ";
              errorFlag = -1000;
            }
          }


      // ijkData_1a x jkiData_2a
      art::crossProductDataData(outData, ijkData_1a, jkiData_2a);

      // cell 0 should contain i x j = k
      if( !( outData(0,0,0)==0.0 && outData(0,0,1)==0.0 && outData(0,0,2)==1.0 &&
             outData(0,1,0)==0.0 && outData(0,1,1)==0.0 && outData(0,1,2)==1.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (2): i x j != k; ";
        errorFlag = -1000;
      }

      // cell 1 should contain j x k = i
      if( !( outData(1,0,0)==1.0 && outData(1,0,1)==0.0 && outData(1,0,2)==0.0 &&
             outData(1,1,0)==1.0 && outData(1,1,1)==0.0 && outData(1,1,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (3): j x k != i; ";
        errorFlag = -1000;
      }

      // cell 2 should contain k x i = j
      if( !( outData(2,0,0)==0.0 && outData(2,0,1)==1.0 && outData(2,0,2)==0.0 &&
             outData(2,1,0)==0.0 && outData(2,1,1)==1.0 && outData(2,1,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (4): k x i != j; ";
        errorFlag = -1000;
      }


      // ijkData_1a x kijData_2a
      art::crossProductDataData(outData, ijkData_1a, kijData_2a);

      // cell 0 should contain i x k = -j
      if( !( outData(0,0,0)==0.0 && outData(0,0,1)==-1.0 && outData(0,0,2)==0.0 &&
             outData(0,1,0)==0.0 && outData(0,1,1)==-1.0 && outData(0,1,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (5): i x k != -j; ";
        errorFlag = -1000;
      }

      // cell 1 should contain j x i = -k
      if( !( outData(1,0,0)==0.0 && outData(1,0,1)==0.0 && outData(1,0,2)==-1.0 &&
             outData(1,1,0)==0.0 && outData(1,1,1)==0.0 && outData(1,1,2)==-1.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (6): j x i != -k; ";
        errorFlag = -1000;
      }

      // cell 2 should contain k x j = -i
      if( !( outData(2,0,0)==-1.0 && outData(2,0,1)==0.0 && outData(2,0,2)==0.0 &&
             outData(2,1,0)==-1.0 && outData(2,1,1)==0.0 && outData(2,1,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (7): k x j != -i; ";
        errorFlag = -1000;
      }


      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 2.b: 3D crossProductDataData operations: (C,P,D) and (P,D)             |\n"\
        << "===============================================================================\n";
      /*
       *        ijkData_1b         ijkData_2b        expected result in (C,P,D) array
       *          i,i,i               i,j,k               0, k,-j
       *          j,j,j                                  -k, 0, i
       *          k,k,k                                   j,-i, 0
       */
      // input data is (P,D)
      DynRankView ConstructWithLabel(ijkData_2b, 3, 3);
      // F=0 at 3 points is (i,j,k)
      ijkData_2b(0, 0) = 1.0;   ijkData_2b(0, 1) = 0.0;   ijkData_2b(0, 2) = 0.0;
      ijkData_2b(1, 0) = 0.0;   ijkData_2b(1, 1) = 1.0;   ijkData_2b(1, 2) = 0.0;
      ijkData_2b(2, 0) = 0.0;   ijkData_2b(2, 1) = 0.0;   ijkData_2b(2, 2) = 1.0;

      // Output array is (C,P,D)
      outData = DynRankView("outData", 3, 3, 3);
      art::crossProductDataData(outData, ijkData_1b, ijkData_2b);
      // checks for C = 0
      if( !(outData(0,0,0)==0.0 && outData(0,0,1)==0.0 && outData(0,0,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (8): i x i != 0; ";
        errorFlag = -1000;
      }
      if( !(outData(0,1,0)==0.0 && outData(0,1,1)==0.0 && outData(0,1,2)==1.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (9): i x j != k; ";
        errorFlag = -1000;
      }
      if( !(outData(0,2,0)==0.0 && outData(0,2,1)==-1.0 && outData(0,2,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (10): i x k != -j; ";
        errorFlag = -1000;
      }

      // checks for C = 1
      if( !(outData(1,0,0)==0.0 && outData(1,0,1)==0.0 && outData(1,0,2)==-1.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (11): j x i != -k; ";
        errorFlag = -1000;
      }
      if( !(outData(1,1,0)==0.0 && outData(1,1,1)==0.0 && outData(1,1,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (12): j x j != 0; ";
        errorFlag = -1000;
      }
      if( !(outData(1,2,0)==1.0 && outData(1,2,1)==0.0 && outData(1,2,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (13): j x k != i; ";
        errorFlag = -1000;
      }

      // checks for C = 2
      if( !(outData(2,0,0)==0.0 && outData(2,0,1)==1.0 && outData(2,0,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (14): k x i != j; ";
        errorFlag = -1000;
      }
      if( !(outData(2,1,0)==-1.0 && outData(2,1,1)==0.0 && outData(2,1,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (15): k x j != -i; ";
        errorFlag = -1000;
      }
      if( !(outData(2,2,0)==0.0 && outData(2,2,1)==0.0 && outData(2,2,2)==0.0) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (16): k x k != 0; ";
        errorFlag = -1000;
      }


      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 2.c: 2D crossProductDataData operations: (C,P,D) and (C,P,D)           |\n"\
        << "===============================================================================\n";
      /*
       *     ijData_1c    jiData_2c
       *        i,i          j,j
       *        j,j          i,i
       */
      DynRankView ConstructWithLabel(jiData_2c, 2, 2, 2);
      // C=0 contains j
      jiData_2c(0, 0, 0) = 0.0;   jiData_2c(0, 0, 1) = 1.0;
      jiData_2c(0, 1, 0) = 0.0;   jiData_2c(0, 1, 1) = 1.0;
      // C=1 contains i
      jiData_2c(1, 0, 0) = 1.0;   jiData_2c(1, 0, 1) = 0.0;
      jiData_2c(1, 1, 0) = 1.0;   jiData_2c(1, 1, 1) = 0.0;


      // ijData_1c x ijData_1c: outData should contain ixi=0, jxj=0
      outData = DynRankView("outData", 2,2);
      art::crossProductDataData(outData, ijData_1c, ijData_1c);

      for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++){
          if(outData(i,j) != 0) {
            *outStream << "\n\nINCORRECT crossProductDataData (17): i x i or j x j != 0; ";
            errorFlag = -1000;
          }
        }

      // ijData_1c x jiData_1c: outData should contain ixi=0, jxj=0
      art::crossProductDataData(outData, ijData_1c, jiData_2c);

      if( !(outData(0,0)==1.0 && outData(0,1)==1.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (18): i x j != 1; ";
        errorFlag = -1000;
      }
      if( !(outData(1,0)==-1.0 && outData(1,1)==-1.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (19): j x i != -1; ";
        errorFlag = -1000;
      }

      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 2.d: 2D crossProductDataData operations: (C,P,D) and (P,D)             |\n"\
        << "===============================================================================\n";
      /*
       *        ijData_1c      ijData_2d        expected result in (C,P) array
       *          i,i             i,j                 0, 1
       *          j,j                                -1, 0
       */
      DynRankView ConstructWithLabel(ijData_2d, 2, 2);
      ijData_2d(0, 0) = 1.0;   ijData_2d(0, 1) = 0.0;
      ijData_2d(1, 0) = 0.0;   ijData_2d(1, 1) = 1.0;

      art::crossProductDataData(outData, ijData_1c, ijData_2d);

      if( !(outData(0,0)==0.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (20): i x i != 0; ";
        errorFlag = -1000;
      }
      if( !(outData(0,1)==1.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (21): i x j != 1; ";
        errorFlag = -1000;
      }
      if( !(outData(1,0)==-1.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (22): j x i != -1; ";
        errorFlag = -1000;
      }
      if( !(outData(1,1)==0.0 ) ) {
        *outStream << "\n\nINCORRECT crossProductDataData (23): j x j != 0; ";
        errorFlag = -1000;
      }


      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 3.a: outerProductDataField operations: (C,P,D) and (C,F,P,D)           |\n"\
        << "===============================================================================\n";
      /*
       *        ijkData_1a      ijkFields_1a       Expected result in (C,F,P,D,D) array:
       *          i,i          i,i  j,j,  k,k            (0,0) (0,1) (0,2)
       *          j,j          i,i  j,j,  k,k            (1,0) (1,1) (1,2)
       *          k,k          i,i  j,j,  k,k            (2,0) (2,1) (2,2)
       *   Indicates the only non-zero element of (C,F,P,*,*), specifically,
       *   element with row = cell and col = field should equal 1; all other should equal 0
       */

      outFields = DynRankView("outFields", 3, 3, 2, 3, 3);
      art::outerProductDataField(outFields, ijkData_1a, ijkFields_1a);

      for(size_type cell = 0; cell < ijkData_1a.extent(0); cell++){
        for(size_type field = 0; field < ijkFields_1a.extent(1); field++){
          for(size_type point = 0; point < ijkData_1a.extent(1); point++){
            for(size_type row = 0; row < ijkData_1a.extent(2); row++){
              for(size_type col = 0; col < ijkData_1a.extent(2); col++){

                // element with row = cell and col = field should equal 1; all other should equal 0
                if( (row == cell && col == field) ){
                  if(outFields(cell, field, point, row, col) != 1.0) {
                    *outStream << "\n\nINCORRECT outerProductDataField (1): computed value is "
                    << outFields(cell, field, point, row, col) << " whereas correct value is 1.0";
                    errorFlag = -1000;
                  }
                }
                else {
                  if(outFields(cell, field, point, row, col) != 0.0) {
                    *outStream << "\n\nINCORRECT outerProductDataField (2): computed value is "
                    << outFields(cell, field, point, row, col) << " whereas correct value is 0.0";
                    errorFlag = -1000;
                  }
                } // if
              }// col
            }// row
          }// point
        }// field
      }// cell

      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 3.b: outerProductDataField operations: (C,P,D) and (F,P,D)             |\n"\
        << "===============================================================================\n";
      /*
       *        ijkData_1b         ijkFields_1b   expected result in (C,F,P,D,D) array
       *          i,i,i               i,j,k            (0,0) (0,1) (0,2)
       *          j,j,j                                (1,0) (1,1) (1,2)
       *          k,k,k                                (2,0) (2,1) (2,2)
       *   Indicates the only non-zero element of (C,F,P,*,*), specifically,
       *   element with row = cell and col = point should equal 1; all other should equal 0
       */

      outFields = DynRankView("outFields", 3, 1, 3, 3, 3);
      art::outerProductDataField(outFields, ijkData_1b, ijkFields_1b);

      for(size_type cell = 0; cell < ijkData_1b.extent(0); cell++){
        for(size_type field = 0; field < ijkFields_1b.extent(0); field++){
          for(size_type point = 0; point < ijkData_1b.extent(1); point++){
            for(size_type row = 0; row < ijkData_1b.extent(2); row++){
              for(size_type col = 0; col < ijkData_1b.extent(2); col++){

                // element with row = cell and col = point should equal 1; all other should equal 0
                if( (row == cell && col == point) ){
                  if(outFields(cell, field, point, row, col) != 1.0) {
                    *outStream << "\n\nINCORRECT outerProductDataField (3): computed value is "
                    << outFields(cell, field, point, row, col) << " whereas correct value is 1.0";
                    errorFlag = -1000;

                  }
                }
                else {
                  if(outFields(cell, field, point, row, col) != 0.0) {
                    *outStream << "\n\nINCORRECT outerProductDataField (4): computed value is "
                    << outFields(cell, field, point, row, col) << " whereas correct value is 0.0";
                    errorFlag = -1000;
                  }
                } // if
              }// col
            }// row
          }// point
        }// field
      }// cell

      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 4.a: outerProductDataData operations: (C,P,D) and (C,P,D)              |\n"\
        << "===============================================================================\n";
      /*
       *     ijkData_1a    jkiData_2a   kijData_2a
       *        i,i          j,j          k,k
       *        j,j          k,k          i,i
       *        k,k          i,i          j,j
       *
       *     Expected results are stated with each test case.
       */
      outData = DynRankView("outData", 3, 2, 3, 3);

      art::outerProductDataData(outData, ijkData_1a, ijkData_1a);
      for(size_type cell = 0; cell < ijkData_1a.extent(0); cell++){
          for(size_type point = 0; point < ijkData_1a.extent(1); point++){
            for(size_type row = 0; row < ijkData_1a.extent(2); row++){
              for(size_type col = 0; col < ijkData_1a.extent(2); col++){

                // element with row = cell and col = cell should equal 1; all other should equal 0
                if( (row == cell && col == cell) ){
                  if(outData(cell, point, row, col) != 1.0) {
                    *outStream << "\n\nINCORRECT outerProductDataData (1): computed value is "
                    << outData(cell, point, row, col) << " whereas correct value is 1.0";
                    errorFlag = -1000;
                  }
                }
                else {
                  if(outData(cell, point, row, col) != 0.0) {
                    *outStream << "\n\nINCORRECT outerProductDataData (2): computed value is "
                    << outData(cell, point, row, col) << " whereas correct value is 0.0";
                    errorFlag = -1000;
                  }
                } // if
              }// col
            }// row
          }// point
      }// cell

      deep_copy(outData, 0.0);
      art::outerProductDataData(outData, ijkData_1a, jkiData_2a);
      for(size_type cell = 0; cell < ijkData_1a.extent(0); cell++){
        for(size_type point = 0; point < ijkData_1a.extent(1); point++){
          for(size_type row = 0; row < ijkData_1a.extent(2); row++){
            for(size_type col = 0; col < ijkData_1a.extent(2); col++){

              // element with row = cell and col = cell + 1 (mod 3) should equal 1; all other should equal 0
              if( (row == cell && col == (cell + 1) % 3) ){
                if(outData(cell, point, row, col) != 1.0) {
                  *outStream << "\n\nINCORRECT outerProductDataData (3): computed value is "
                  << outData(cell, point, row, col) << " whereas correct value is 1.0";
                  errorFlag = -1000;
                }
              }
              else {
                if(outData(cell, point, row, col) != 0.0) {
                  *outStream << "\n\nINCORRECT outerProductDataData (4): computed value is "
                  << outData(cell, point, row, col) << " whereas correct value is 0.0";
                  errorFlag = -1000;
                }
              } // if
            }// col
          }// row
        }// point
      }// cell


      deep_copy(outData, 0.0);
      art::outerProductDataData(outData, ijkData_1a, kijData_2a);
      for(size_type cell = 0; cell < ijkData_1a.extent(0); cell++){
        for(size_type point = 0; point < ijkData_1a.extent(1); point++){
          for(size_type row = 0; row < ijkData_1a.extent(2); row++){
            for(size_type col = 0; col < ijkData_1a.extent(2); col++){

              // element with row = cell and col = cell + 2 (mod 3) should equal 1; all other should equal 0
              if( (row == cell && col == (cell + 2) % 3) ){
                if(outData(cell, point, row, col) != 1.0) {
                  *outStream << "\n\nINCORRECT outerProductDataData (5): computed value is "
                  << outData(cell, point, row, col) << " whereas correct value is 1.0";
                  errorFlag = -1000;
                }
              }
              else {
                if(outData(cell, point, row, col) != 0.0) {
                  *outStream << "\n\nINCORRECT outerProductDataData (6): computed value is "
                  << outData(cell, point, row, col) << " whereas correct value is 0.0";
                  errorFlag = -1000;
                }
              } // if
            }// col
          }// row
        }// point
      }// cell


      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 4.b: outerProductDataData operations: (C,P,D) and (P,D)                |\n"\
        << "===============================================================================\n";
      /*
       *        ijkData_1b         ijkData_2b   expected result in (C,P,D,D) array
       *          i,i,i               i,j,k            (0,0) (0,1) (0,2)
       *          j,j,j                                (1,0) (1,1) (1,2)
       *          k,k,k                                (2,0) (2,1) (2,2)
       *   Indicates the only non-zero element of (C,P,*,*), specifically,
       *   element with row = cell and col = point should equal 1; all other should equal 0
       */
      outData = DynRankView("outData", 3,3,3,3);
      art::outerProductDataData(outData, ijkData_1b, ijkData_2b);
      for(size_type cell = 0; cell < ijkData_1b.extent(0); cell++){
        for(size_type point = 0; point < ijkData_1b.extent(1); point++){
          for(size_type row = 0; row < ijkData_1b.extent(2); row++){
            for(size_type col = 0; col < ijkData_1b.extent(2); col++){

              // element with row = cell and col = cell + 2 (mod 3) should equal 1; all other should equal 0
              if( (row == cell && col == point) ){
                if(outData(cell, point, row, col) != 1.0) {
                  *outStream << "\n\nINCORRECT outerProductDataData (7): computed value is "
                  << outData(cell, point, row, col) << " whereas correct value is 1.0";
                  errorFlag = -1000;
                }
              }
              else {
                if(outData(cell, point, row, col) != 0.0) {
                  *outStream << "\n\nINCORRECT outerProductDataData (8): computed value is "
                  << outData(cell, point, row, col) << " whereas correct value is 0.0";
                  errorFlag = -1000;
                }
              } // if
            }// col
          }// row
        }// point
      }// cell

      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 5.a: matvecProductDataField operations: (C,P,D,D) and (C,F,P,D)        |\n"\
        << "===============================================================================\n";
      /*
       *         inputMat              inputVecFields           outFields
       *     1  1  1     0  0  0     0   1       -1  -1      0   3      0   0
       *    -1  2 -1    -1 -2 -3     0   1       -1   1      0   0      6   2
       *     1  2  3    -2  6 -4     0   1       -1  -1      0   6      0  12
       */

      // (C,P,D,D)
      DynRankView ConstructWithLabel(inputMat, 2,1,3,3);
      // cell 0
      inputMat(0,0,0,0) = 1.0;  inputMat(0,0,0,1) = 1.0;  inputMat(0,0,0,2) = 1.0;
      inputMat(0,0,1,0) =-1.0;  inputMat(0,0,1,1) = 2.0;  inputMat(0,0,1,2) =-1.0;
      inputMat(0,0,2,0) = 1.0;  inputMat(0,0,2,1) = 2.0;  inputMat(0,0,2,2) = 3.0;
      // cell 1
      inputMat(1,0,0,0) = 0.0;  inputMat(1,0,0,1) = 0.0;  inputMat(1,0,0,2) = 0.0;
      inputMat(1,0,1,0) =-1.0;  inputMat(1,0,1,1) =-2.0;  inputMat(1,0,1,2) =-3.0;
      inputMat(1,0,2,0) =-2.0;  inputMat(1,0,2,1) = 6.0;  inputMat(1,0,2,2) =-4.0;

      // (C,F,P,D)
      DynRankView ConstructWithLabel(inputVecFields, 2,2,1,3);
      // cell 0; fields 0,1
      inputVecFields(0,0,0,0) = 0.0;  inputVecFields(0,0,0,1) = 0.0;  inputVecFields(0,0,0,2) = 0.0;
      inputVecFields(0,1,0,0) = 1.0;  inputVecFields(0,1,0,1) = 1.0;  inputVecFields(0,1,0,2) = 1.0;
      // cell 1; fields 0,1
      inputVecFields(1,0,0,0) =-1.0;  inputVecFields(1,0,0,1) =-1.0;  inputVecFields(1,0,0,2) =-1.0;
      inputVecFields(1,1,0,0) =-1.0;  inputVecFields(1,1,0,1) = 1.0;  inputVecFields(1,1,0,2) =-1.0;

      // (C,F,P,D) - true
      DynRankView ConstructWithLabel(outFieldsCorrect, 2,2,1,3);
      // cell 0; fields 0,1
      outFieldsCorrect(0,0,0,0) = 0.0;  outFieldsCorrect(0,0,0,1) = 0.0;  outFieldsCorrect(0,0,0,2) = 0.0;
      outFieldsCorrect(0,1,0,0) = 3.0;  outFieldsCorrect(0,1,0,1) = 0.0;  outFieldsCorrect(0,1,0,2) = 6.0;
      // cell 1; fields 0,1
      outFieldsCorrect(1,0,0,0) = 0.0;  outFieldsCorrect(1,0,0,1) = 6.0;  outFieldsCorrect(1,0,0,2) = 0.0;
      outFieldsCorrect(1,1,0,0) = 0.0;  outFieldsCorrect(1,1,0,1) = 2.0;  outFieldsCorrect(1,1,0,2) = 12.0;

      // (C,F,P,D)
      outFields = DynRankView("outFields", 2,2,1,3);
      art::matvecProductDataField(outFields, inputMat, inputVecFields);

      // test loop
      for(size_type cell = 0; cell < outFields.extent(0); cell++){
        for(size_type field = 0; field < outFields.extent(1); field++){
          for(size_type point = 0; point < outFields.extent(2); point++){
            for(size_type row = 0; row < outFields.extent(3); row++){
              if(outFields(cell, field, point, row) != outFieldsCorrect(cell, field, point, row)) {
                *outStream << "\n\nINCORRECT matvecProductDataField (1): \n value at multi-index ("
                << cell << "," << field << "," << point << "," << row << ") = "
                << outFields(cell, field, point, row) << " but correct value is "
                << outFieldsCorrect(cell, field, point, row) <<"\n";
                errorFlag = -1000;
              }
            }//row
          }// point
        }// field
      }// cell


      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 5.b: matvecProductDataField operations: (C,P,D,D) and (F,P,D)          |\n"\
        << "===============================================================================\n";
      /*
       *         inputMat              inputVecFields        outFields
       *     1  1  1     0  0  0       0   1  -1  -1       0  3 -3 -1     0  0  0  0
       *    -1  2 -1    -1 -2 -3       0   1  -1   1       0  0  0  4     0 -6  6  2
       *     1  2  3    -2  6 -4       0   1  -1  -1       0  6 -6 -2     0  0  0 12
       *     Use the same 4 vector fields as above but formatted as (F,P,D) array
       */
      // (C,F,P,D)
      inputVecFields= DynRankView("inputVecFields", 4,1,3);
      // fields 0,1,2,3
      inputVecFields(0,0,0) = 0.0;  inputVecFields(0,0,1) = 0.0;  inputVecFields(0,0,2) = 0.0;
      inputVecFields(1,0,0) = 1.0;  inputVecFields(1,0,1) = 1.0;  inputVecFields(1,0,2) = 1.0;
      inputVecFields(2,0,0) =-1.0;  inputVecFields(2,0,1) =-1.0;  inputVecFields(2,0,2) =-1.0;
      inputVecFields(3,0,0) =-1.0;  inputVecFields(3,0,1) = 1.0;  inputVecFields(3,0,2) =-1.0;

      // (C,F,P,D) - true
      outFieldsCorrect= DynRankView("outFieldsCorrect", 2,4,1,3);
      // cell 0; fields 0,1,2,3
      outFieldsCorrect(0,0,0,0) = 0.0;  outFieldsCorrect(0,0,0,1) = 0.0;  outFieldsCorrect(0,0,0,2) = 0.0;
      outFieldsCorrect(0,1,0,0) = 3.0;  outFieldsCorrect(0,1,0,1) = 0.0;  outFieldsCorrect(0,1,0,2) = 6.0;
      outFieldsCorrect(0,2,0,0) =-3.0;  outFieldsCorrect(0,2,0,1) = 0.0;  outFieldsCorrect(0,2,0,2) =-6.0;
      outFieldsCorrect(0,3,0,0) =-1.0;  outFieldsCorrect(0,3,0,1) = 4.0;  outFieldsCorrect(0,3,0,2) =-2.0;
      // cell 1; fields 0,1,2,3
      outFieldsCorrect(1,0,0,0) = 0.0;  outFieldsCorrect(1,0,0,1) = 0.0;  outFieldsCorrect(1,0,0,2) = 0.0;
      outFieldsCorrect(1,1,0,0) = 0.0;  outFieldsCorrect(1,1,0,1) =-6.0;  outFieldsCorrect(1,1,0,2) = 0.0;
      outFieldsCorrect(1,2,0,0) = 0.0;  outFieldsCorrect(1,2,0,1) = 6.0;  outFieldsCorrect(1,2,0,2) = 0.0;
      outFieldsCorrect(1,3,0,0) = 0.0;  outFieldsCorrect(1,3,0,1) = 2.0;  outFieldsCorrect(1,3,0,2) =12.0;

      // (C,F,P,D)
      outFields = DynRankView("outFields", 2,4,1,3);
      art::matvecProductDataField(outFields, inputMat, inputVecFields);

      // test loop
      for(size_type cell = 0; cell < outFields.extent(0); cell++){
        for(size_type field = 0; field < outFields.extent(1); field++){
          for(size_type point = 0; point < outFields.extent(2); point++){
            for(size_type row = 0; row < outFields.extent(3); row++){
              if(outFields(cell, field, point, row) != outFieldsCorrect(cell, field, point, row)) {
                *outStream << "\n\nINCORRECT matvecProductDataField (2): \n value at multi-index ("
                << cell << "," << field << "," << point << "," << row << ") = "
                << outFields(cell, field, point, row) << " but correct value is "
                << outFieldsCorrect(cell, field, point, row) <<"\n";
                errorFlag = -1000;
              }
            }//row
          }// point
        }// field
      }// cell


      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 5.c: matvecProductDataField random tests: branch inputFields(C,F,P,D)  |\n"\
        << "===============================================================================\n";
      /*
       *  d1 is spatial dimension and should be 1, 2 or 3. If d1>3, the rst function 'inverse' will fail
       */
      {// test 5.c scope
        int c=5, p=9, f=7, d1=3;


        DynRankView ConstructWithLabel(in_c_f_p_d, c, f, p, d1);
        DynRankView ConstructWithLabel(out_c_f_p_d, c, f, p, d1);
        DynRankView ConstructWithLabel(outi_c_f_p_d, c, f, p, d1);

        DynRankView ConstructWithLabel(data_c_p, c, p);
        DynRankView ConstructWithLabel(datainv_c_p, c, p);
        DynRankView ConstructWithLabel(data_c_1, c, 1);
        DynRankView ConstructWithLabel(datainv_c_1, c, 1);
        DynRankView ConstructWithLabel(data_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(data_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(data_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(data_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_1_d_d, c, 1, d1, d1);
        /***********************************************************************************************
          *                          Constant diagonal tensor: inputData(C,P)                          *
          **********************************************************************************************/

        // fill with random numbers
         for (auto i=0;i<c;++i) {
           for (auto j=0;j<f;++j)
             for (auto k=0;k<p;++k)
               for (auto m=0;m<d1;++m)
                 in_c_f_p_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();

           for (auto j=0;j<p;++j) {
             data_c_p(i,j)  = Teuchos::ScalarTraits<value_type>::random();
             datainv_c_p(i,j) = 1./data_c_p(i,j);
           }
           data_c_1(i,0)  = Teuchos::ScalarTraits<value_type>::random();
           datainv_c_1(i,0) = 1./data_c_1(i,0);
         }


        // Tensor values vary by point:
        art::matvecProductDataField(out_c_f_p_d, data_c_p, in_c_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_p, out_c_f_p_d);
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (3): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
    
        // Tensor values do not vary by point:
        art::matvecProductDataField(out_c_f_p_d, data_c_1, in_c_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_1, out_c_f_p_d);
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (4): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *                     Non-onstant diagonal tensor: inputData(C,P,D)                          *
          **********************************************************************************************/
    

        // fill with random numbers
         for (auto i=0;i<c;++i) {
           for (auto j=0;j<f;++j)
             for (auto k=0;k<p;++k)
               for (auto m=0;m<d1;++m)
                 in_c_f_p_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();

           for (auto j=0;j<p;++j)
             for (auto k=0;k<d1;++k){
               data_c_p_d(i,j,k)  = Teuchos::ScalarTraits<value_type>::random();
               datainv_c_p_d(i,j,k) = 1./data_c_p_d(i,j,k);
             }
           for (auto k=0;k<d1;++k) {
             data_c_1_d(i,0,k)  = Teuchos::ScalarTraits<value_type>::random();
             datainv_c_1_d(i,0,k) = 1./data_c_1_d(i,0,k);
           }
         }

        // Tensor values vary by point:
        art::matvecProductDataField(out_c_f_p_d, data_c_p_d, in_c_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_p_d, out_c_f_p_d);
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (5): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point
        art::matvecProductDataField(out_c_f_p_d, data_c_1_d, in_c_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_1_d, out_c_f_p_d);
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (6): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *                               Full tensor: inputData(C,P,D,D)                              *
          **********************************************************************************************/

        for (int i=0; i < c; i++) {
          for (int j=0; j < f; j++)
            for (int k=0; k < p; k++)
              for (int m=0; m<d1; m++)
                in_c_f_p_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
          }
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }
        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);

        // Tensor values vary by point: test "N" and "T" options (no transpose/transpose)
        art::matvecProductDataField(out_c_f_p_d, data_c_p_d_d, in_c_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_p_d_d, out_c_f_p_d);
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (7): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matvecProductDataField(out_c_f_p_d, data_c_p_d_d, in_c_f_p_d, 't');
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_p_d_d, out_c_f_p_d, 't');
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (8): check matrix inverse property, w/ double transpose\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point: test "N" and "T" options (no transpose/transpose)
        art::matvecProductDataField(out_c_f_p_d, data_c_1_d_d, in_c_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_1_d_d, out_c_f_p_d);
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (9): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matvecProductDataField(out_c_f_p_d, data_c_1_d_d, in_c_f_p_d, 't');
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_1_d_d, out_c_f_p_d, 't');
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (10): check matrix inverse property, w/ double transpose\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *             Full tensor: inputData(C,P,D,D) test inverse transpose                         *
          **********************************************************************************************/


        for (int i=0; i < c; i++) {
          for (int j=0; j < f; j++)
            for (int k=0; k < p; k++)
              for (int m=0; m<d1; m++)
                in_c_f_p_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
          }
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }



        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::transpose(datainvtrn_c_p_d_d, datainv_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);
        rst::transpose(datainvtrn_c_1_d_d, datainv_c_1_d_d);


        // Tensor values vary by point:
        art::matvecProductDataField(out_c_f_p_d, data_c_p_d_d, in_c_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainvtrn_c_p_d_d, out_c_f_p_d, 't');
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (11): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point
        art::matvecProductDataField(out_c_f_p_d, data_c_1_d_d, in_c_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainvtrn_c_1_d_d, out_c_f_p_d, 't');
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (12): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
      }// test 5.c scope


      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 5.d: matvecProductDataField random tests: branch inputFields(F,P,D)    |\n"\
        << "===============================================================================\n";
      /*
       *  d1 is the spatial dimension and should be 1, 2 or 3. If d1>3, the RealSpaceTools function 'inverse' will fail
       */
      {// test 5.d scope
        int c=5, p=9, f=7, d1=3;

        DynRankView ConstructWithLabel(in_f_p_d, f, p, d1);
        DynRankView ConstructWithLabel(in_c_f_p_d, c, f, p, d1);
        DynRankView ConstructWithLabel(data_c_p, c, p);
        DynRankView ConstructWithLabel(datainv_c_p, c, p);
        DynRankView ConstructWithLabel(data_c_1, c, 1);
        DynRankView ConstructWithLabel(datainv_c_1, c, 1);
        DynRankView ConstructWithLabel(data_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(data_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(data_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(data_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(data_c_p_one, c, p);
        DynRankView ConstructWithLabel(data_c_1_one, c, 1);
        DynRankView ConstructWithLabel(out_c_f_p_d, c, f, p, d1);
        DynRankView ConstructWithLabel(outi_c_f_p_d, c, f, p, d1);
        /***********************************************************************************************
          *                          Constant diagonal tensor: inputData(C,P)                          *
          **********************************************************************************************/

        for (int i=0; i < c; i++) {
          for (int j=0; j < f; j++)
            for (int k=0; k < p; k++)
              for (int m=0; m<d1; m++)
                in_c_f_p_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
          for (int j=0; j < p; j++) {
            data_c_p(i, j) = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_p(i, j) = 1.0 / data_c_p(i,j);
            data_c_p_one(i, j) = 1.0;
          }
        data_c_1(i, 0) = Teuchos::ScalarTraits<value_type>::random();
        datainv_c_1(i, 0) = 1.0 / data_c_1(i,0);
        }

        // Tensor values vary by point
        art::matvecProductDataField(in_c_f_p_d, data_c_p_one, in_f_p_d);
        art::matvecProductDataField(out_c_f_p_d, data_c_p, in_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_p, out_c_f_p_d);
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (13): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point
        art::matvecProductDataField(in_c_f_p_d, data_c_p_one, in_f_p_d);
        art::matvecProductDataField(out_c_f_p_d, data_c_1, in_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_1, out_c_f_p_d);
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (14): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *                       Non-constant diagonal tensor: inputData(C,P,D)                       *
          **********************************************************************************************/

    
        for (int i=0; i < c; i++) {
          for (int j=0; j < f; j++)
            for (int k=0; k < p; k++)
              for (int m=0; m<d1; m++)
                in_c_f_p_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
          for (int j=0; j < p; j++)
            for (int m=0; m<d1; m++) {
              data_c_p_d(i, j, m) = Teuchos::ScalarTraits<value_type>::random();
              datainv_c_p_d(i, j, m) = 1.0 / data_c_p_d(i,j,m);
          }
          for (int m=0; m<d1; m++) {
            data_c_1_d(i, 0, m) = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_1_d(i, 0, m) = 1.0 / data_c_1_d(i,0, m);
          }
        }

        // Tensor values vary by point:
        art::matvecProductDataField(in_c_f_p_d, data_c_p_one, in_f_p_d);
        art::matvecProductDataField(out_c_f_p_d, data_c_p_d, in_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_p_d, out_c_f_p_d);
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (15): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point:
        art::matvecProductDataField(in_c_f_p_d, data_c_p_one, in_f_p_d);
        art::matvecProductDataField(out_c_f_p_d, data_c_1_d, in_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_1_d, out_c_f_p_d);
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (16): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *                              Full tensor: inputData(C,P,D,D)                               *
          **********************************************************************************************/

        for (int i=0; i < f; i++)
          for (int j=0; j < p; j++)
            for (int m=0; m<d1; m++)
              in_f_p_d(i, j, m) = Teuchos::ScalarTraits<value_type>::random();
        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++)
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();

          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }

        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);

        // Tensor values vary by point: test "N" and "T" (no-transpose/transpose) options
        art::matvecProductDataField(in_c_f_p_d, data_c_p_one, in_f_p_d);
        art::matvecProductDataField(out_c_f_p_d, data_c_p_d_d, in_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_p_d_d, out_c_f_p_d);
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (17): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matvecProductDataField(in_c_f_p_d, data_c_p_one, in_f_p_d);
        art::matvecProductDataField(out_c_f_p_d, data_c_p_d_d, in_f_p_d, 't');
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_p_d_d, out_c_f_p_d, 't');
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (18): check matrix inverse property, w/ double tvalue_typese\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point: test "N" and "T" (no-transpose/transpose) options
        art::matvecProductDataField(in_c_f_p_d, data_c_p_one, in_f_p_d);
        art::matvecProductDataField(out_c_f_p_d, data_c_1_d_d, in_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_1_d_d, out_c_f_p_d);
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (19): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matvecProductDataField(in_c_f_p_d, data_c_p_one, in_f_p_d);
        art::matvecProductDataField(out_c_f_p_d, data_c_1_d_d, in_f_p_d, 't');
        art::matvecProductDataField(outi_c_f_p_d, datainv_c_1_d_d, out_c_f_p_d, 't');
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (20): check matrix inverse property, w/ double transpose\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *             Full tensor: inputData(C,P,D,D) test inverse transpose                         *
          **********************************************************************************************/

        for (int i=0; i < f; i++)
          for (int j=0; j < p; j++)
            for (int m=0; m<d1; m++)
              in_f_p_d(i,j,m) = Teuchos::ScalarTraits<value_type>::random();

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
          }
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }


        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::transpose(datainvtrn_c_p_d_d, datainv_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);
        rst::transpose(datainvtrn_c_1_d_d, datainv_c_1_d_d);



        // Tensor values vary by point:
        art::matvecProductDataField(in_c_f_p_d, data_c_p_one, in_f_p_d);
        art::matvecProductDataField(out_c_f_p_d, data_c_p_d_d, in_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainvtrn_c_p_d_d, out_c_f_p_d, 't');
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (21): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point:
        art::matvecProductDataField(in_c_f_p_d, data_c_p_one, in_f_p_d);
        art::matvecProductDataField(out_c_f_p_d, data_c_1_d_d, in_f_p_d);
        art::matvecProductDataField(outi_c_f_p_d, datainvtrn_c_1_d_d, out_c_f_p_d, 't');
        rst::subtract(outi_c_f_p_d, in_c_f_p_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataField (22): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
      }// test 5.d scope

      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 6.a: matvecProductDataData operations: (C,P,D,D) and (C,P,D)           |\n"\
        << "===============================================================================\n";
      /*
       *                     inputMat                     inputVecFields       outFields
       *      1  1  1   0  0  0    1  1  1    0  0  0     0   1  -1  -1      0   0  -3   0
       *     -1  2 -1  -1 -2 -3   -1  2 -1   -1 -2 -3     0   1  -1   1      0  -6   0   2
       *      1  2  3  -2  6 -4    1  2  3   -2  6 -4     0   1  -1  -1      0   0  -6  12
       */

      // (C,P,D,D)
      inputMat = DynRankView("inputMat", 4,1,3,3);
      // cell 0
      inputMat(0,0,0,0) = 1.0;  inputMat(0,0,0,1) = 1.0;  inputMat(0,0,0,2) = 1.0;
      inputMat(0,0,1,0) =-1.0;  inputMat(0,0,1,1) = 2.0;  inputMat(0,0,1,2) =-1.0;
      inputMat(0,0,2,0) = 1.0;  inputMat(0,0,2,1) = 2.0;  inputMat(0,0,2,2) = 3.0;
      // cell 1
      inputMat(1,0,0,0) = 0.0;  inputMat(1,0,0,1) = 0.0;  inputMat(1,0,0,2) = 0.0;
      inputMat(1,0,1,0) =-1.0;  inputMat(1,0,1,1) =-2.0;  inputMat(1,0,1,2) =-3.0;
      inputMat(1,0,2,0) =-2.0;  inputMat(1,0,2,1) = 6.0;  inputMat(1,0,2,2) =-4.0;
      // cell 2
      inputMat(2,0,0,0) = 1.0;  inputMat(2,0,0,1) = 1.0;  inputMat(2,0,0,2) = 1.0;
      inputMat(2,0,1,0) =-1.0;  inputMat(2,0,1,1) = 2.0;  inputMat(2,0,1,2) =-1.0;
      inputMat(2,0,2,0) = 1.0;  inputMat(2,0,2,1) = 2.0;  inputMat(2,0,2,2) = 3.0;
      // cell 3
      inputMat(3,0,0,0) = 0.0;  inputMat(3,0,0,1) = 0.0;  inputMat(3,0,0,2) = 0.0;
      inputMat(3,0,1,0) =-1.0;  inputMat(3,0,1,1) =-2.0;  inputMat(3,0,1,2) =-3.0;
      inputMat(3,0,2,0) =-2.0;  inputMat(3,0,2,1) = 6.0;  inputMat(3,0,2,2) =-4.0;

      // (C,P,D)
      inputVecFields= DynRankView("inputVecFields", 4,1,3);
      inputVecFields(0,0,0) = 0.0;  inputVecFields(0,0,1) = 0.0;  inputVecFields(0,0,2) = 0.0;
      inputVecFields(1,0,0) = 1.0;  inputVecFields(1,0,1) = 1.0;  inputVecFields(1,0,2) = 1.0;
      inputVecFields(2,0,0) =-1.0;  inputVecFields(2,0,1) =-1.0;  inputVecFields(2,0,2) =-1.0;
      inputVecFields(3,0,0) =-1.0;  inputVecFields(3,0,1) = 1.0;  inputVecFields(3,0,2) =-1.0;

      // (C,P,D) - true
      outFieldsCorrect= DynRankView("outFieldsCorrect", 4,1,3);
      outFieldsCorrect(0,0,0) = 0.0;  outFieldsCorrect(0,0,1) = 0.0;  outFieldsCorrect(0,0,2) = 0.0;
      outFieldsCorrect(1,0,0) = 0.0;  outFieldsCorrect(1,0,1) =-6.0;  outFieldsCorrect(1,0,2) = 0.0;
      outFieldsCorrect(2,0,0) =-3.0;  outFieldsCorrect(2,0,1) = 0.0;  outFieldsCorrect(2,0,2) =-6.0;
      outFieldsCorrect(3,0,0) = 0.0;  outFieldsCorrect(3,0,1) = 2.0;  outFieldsCorrect(3,0,2) = 12.0;

      // (C,P,D)
      outFields = DynRankView("outFields", 4,1,3);
      art::matvecProductDataData(outFields, inputMat, inputVecFields);

      // test loop
      for(size_type cell = 0; cell < outFields.extent(0); cell++){
        for(size_type point = 0; point < outFields.extent(1); point++){
          for(size_type row = 0; row < outFields.extent(2); row++){
            if(outFields(cell, point, row) != outFieldsCorrect(cell, point, row)) {
              *outStream << "\n\nINCORRECT matvecProductDataData (1): \n value at multi-index ("
              << cell << "," << point << "," << row << ") = "
              << outFields(cell, point, row) << " but correct value is "
              << outFieldsCorrect(cell, point, row) <<"\n";
              errorFlag = -1000;
            }
          }//row
        }// point
      }// cell


      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 6.b: matvecProductDataData operations: (C,P,D,D) and (P,D)             |\n"\
        << "===============================================================================\n";
      /*
       *                     inputMat                     inputVecFields       outFields
       *      1  1  1   0  0  0    1  1  1    0  0  0     0   1  -1  -1      0   0  -3   0
       *     -1  2 -1  -1 -2 -3   -1  2 -1   -1 -2 -3     0   1  -1   1      0  -6   0   2
       *      1  2  3  -2  6 -4    1  2  3   -2  6 -4     0   1  -1  -1      0   0  -6  12
       */
      // (C,P,D,D)
      inputMat = DynRankView("inputMat", 1,4,3,3);
      // point 0
      inputMat(0,0,0,0) = 1.0;  inputMat(0,0,0,1) = 1.0;  inputMat(0,0,0,2) = 1.0;
      inputMat(0,0,1,0) =-1.0;  inputMat(0,0,1,1) = 2.0;  inputMat(0,0,1,2) =-1.0;
      inputMat(0,0,2,0) = 1.0;  inputMat(0,0,2,1) = 2.0;  inputMat(0,0,2,2) = 3.0;
      // point 1
      inputMat(0,1,0,0) = 0.0;  inputMat(0,1,0,1) = 0.0;  inputMat(0,1,0,2) = 0.0;
      inputMat(0,1,1,0) =-1.0;  inputMat(0,1,1,1) =-2.0;  inputMat(0,1,1,2) =-3.0;
      inputMat(0,1,2,0) =-2.0;  inputMat(0,1,2,1) = 6.0;  inputMat(0,1,2,2) =-4.0;
      // point 2
      inputMat(0,2,0,0) = 1.0;  inputMat(0,2,0,1) = 1.0;  inputMat(0,2,0,2) = 1.0;
      inputMat(0,2,1,0) =-1.0;  inputMat(0,2,1,1) = 2.0;  inputMat(0,2,1,2) =-1.0;
      inputMat(0,2,2,0) = 1.0;  inputMat(0,2,2,1) = 2.0;  inputMat(0,2,2,2) = 3.0;
      // point 3
      inputMat(0,3,0,0) = 0.0;  inputMat(0,3,0,1) = 0.0;  inputMat(0,3,0,2) = 0.0;
      inputMat(0,3,1,0) =-1.0;  inputMat(0,3,1,1) =-2.0;  inputMat(0,3,1,2) =-3.0;
      inputMat(0,3,2,0) =-2.0;  inputMat(0,3,2,1) = 6.0;  inputMat(0,3,2,2) =-4.0;

      // (P,D)
      inputVecFields= DynRankView("inputVecFields", 4,3);
      //
      inputVecFields(0,0) = 0.0;  inputVecFields(0,1) = 0.0;  inputVecFields(0,2) = 0.0;
      inputVecFields(1,0) = 1.0;  inputVecFields(1,1) = 1.0;  inputVecFields(1,2) = 1.0;
      inputVecFields(2,0) =-1.0;  inputVecFields(2,1) =-1.0;  inputVecFields(2,2) =-1.0;
      inputVecFields(3,0) =-1.0;  inputVecFields(3,1) = 1.0;  inputVecFields(3,2) =-1.0;

      // (C,P,D) - true
      outFieldsCorrect= DynRankView("outFieldsCorrect", 1,4,3);
      outFieldsCorrect(0,0,0) = 0.0;  outFieldsCorrect(0,0,1) = 0.0;  outFieldsCorrect(0,0,2) = 0.0;
      outFieldsCorrect(0,1,0) = 0.0;  outFieldsCorrect(0,1,1) =-6.0;  outFieldsCorrect(0,1,2) = 0.0;
      outFieldsCorrect(0,2,0) =-3.0;  outFieldsCorrect(0,2,1) = 0.0;  outFieldsCorrect(0,2,2) =-6.0;
      outFieldsCorrect(0,3,0) = 0.0;  outFieldsCorrect(0,3,1) = 2.0;  outFieldsCorrect(0,3,2) = 12.0;

      // (C,P,D)
      outFields = DynRankView("outFields", 1,4,3);
      art::matvecProductDataData(outFields, inputMat, inputVecFields);

      // test loop
      for(size_type cell = 0; cell < outFields.extent(0); cell++){
        for(size_type point = 0; point < outFields.extent(1); point++){
          for(size_type row = 0; row < outFields.extent(2); row++){
            if(outFields(cell, point, row) != outFieldsCorrect(cell, point, row)) {
              *outStream << "\n\nINCORRECT matvecProductDataData (2): \n value at multi-index ("
              << cell << "," << point << "," << row << ") = "
              << outFields(cell, point, row) << " but correct value is "
              << outFieldsCorrect(cell, point, row) <<"\n";
              errorFlag = -1000;
            }
          }//row
        }// point
      }// cell



      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 6.c: matvecProductDataData random tests: branch inputDataRight(C,P,D)  |\n"\
        << "===============================================================================\n";
      /*
       * Test derived from Test 5.c
       */
      {// test 6.c scope
        int c=5, p=9, d1=3;

        DynRankView ConstructWithLabel(in_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(out_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(outi_c_p_d, c, p, d1);

        DynRankView ConstructWithLabel(data_c_p, c, p);
        DynRankView ConstructWithLabel(datainv_c_p, c, p);
        DynRankView ConstructWithLabel(data_c_1, c, 1);
        DynRankView ConstructWithLabel(datainv_c_1, c, 1);
        DynRankView ConstructWithLabel(data_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(data_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(data_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(data_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_1_d_d, c, 1, d1, d1);
        /***********************************************************************************************
          *                          Constant diagonal tensor: inputDataLeft(C,P)                      *
          **********************************************************************************************/

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++)
              in_c_p_d(i,j,m) = Teuchos::ScalarTraits<value_type>::random();

            data_c_p(i, j) = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_p(i, j) = 1.0/data_c_p(i, j);
          }
          data_c_1(i, 0) = Teuchos::ScalarTraits<value_type>::random();
          datainv_c_1(i, 0) = 1.0/data_c_1(i, 0);
        }



        // Tensor values vary by point:
        art::matvecProductDataData(out_c_p_d, data_c_p, in_c_p_d);
        art::matvecProductDataData(outi_c_p_d, datainv_c_p, out_c_p_d);
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (3): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
    
        // Tensor values do not vary by point:
        art::matvecProductDataData(out_c_p_d, data_c_1, in_c_p_d);
        art::matvecProductDataData(outi_c_p_d, datainv_c_1, out_c_p_d);
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (4): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *                     Non-onstant diagonal tensor: inputDataLeft(C,P,D)                      *
          **********************************************************************************************/

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++)
            for (int m=0; m<d1; m++) {
              in_c_p_d(i,j,m) = Teuchos::ScalarTraits<value_type>::random();
              data_c_p_d(i, j, m) = Teuchos::ScalarTraits<value_type>::random();
              datainv_c_p_d(i, j, m) = 1.0/data_c_p_d(i, j, m);
            }
          for (int m=0; m<d1; m++) {
            data_c_1_d(i, 0, m) = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_1_d(i, 0, m) = 1.0/data_c_1_d(i, 0, m);
          }
        }
    
        // Tensor values vary by point:
        art::matvecProductDataData(out_c_p_d, data_c_p_d, in_c_p_d);
        art::matvecProductDataData(outi_c_p_d, datainv_c_p_d, out_c_p_d);
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (5): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point
        art::matvecProductDataData(out_c_p_d, data_c_1_d, in_c_p_d);
        art::matvecProductDataData(outi_c_p_d, datainv_c_1_d, out_c_p_d);
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (6): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *                            Full tensor: inputDataLeft(C,P,D,D)                             *
          **********************************************************************************************/

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++) {
              in_c_p_d(i, j, m) = Teuchos::ScalarTraits<value_type>::random();
              for (int n=0; n<d1; n++)
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
            }
          }
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }
        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);

        // Tensor values vary by point: test "N" and "T" options (no transpose/transpose)
        art::matvecProductDataData(out_c_p_d, data_c_p_d_d, in_c_p_d);
        art::matvecProductDataData(outi_c_p_d, datainv_c_p_d_d, out_c_p_d);
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (7): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matvecProductDataData(out_c_p_d, data_c_p_d_d, in_c_p_d, 't');
        art::matvecProductDataData(outi_c_p_d, datainv_c_p_d_d, out_c_p_d, 't');
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (8): check matrix inverse property, w/ double transpose\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point: test "N" and "T" options (no transpose/transpose)
        art::matvecProductDataData(out_c_p_d, data_c_1_d_d, in_c_p_d);
        art::matvecProductDataData(outi_c_p_d, datainv_c_1_d_d, out_c_p_d);
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (9): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matvecProductDataData(out_c_p_d, data_c_1_d_d, in_c_p_d, 't');
        art::matvecProductDataData(outi_c_p_d, datainv_c_1_d_d, out_c_p_d, 't');
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (10): check matrix inverse property, w/ double transpose\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *             Full tensor: inputDataLeft(C,P,D,D) test inverse transpose                     *
          **********************************************************************************************/


        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++) {
              in_c_p_d(i, j, m) = Teuchos::ScalarTraits<value_type>::random();
              for (int n=0; n<d1; n++)
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
            }
          }
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }

        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::transpose(datainvtrn_c_p_d_d, datainv_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);
        rst::transpose(datainvtrn_c_1_d_d, datainv_c_1_d_d);

        // Tensor values vary by point:
        art::matvecProductDataData(out_c_p_d, data_c_p_d_d, in_c_p_d);
        art::matvecProductDataData(outi_c_p_d, datainvtrn_c_p_d_d, out_c_p_d, 't');
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (11): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point
        art::matvecProductDataData(out_c_p_d, data_c_1_d_d, in_c_p_d);
        art::matvecProductDataData(outi_c_p_d, datainvtrn_c_1_d_d, out_c_p_d, 't');
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (12): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
      }// test 6.c scope


      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 6.d: matvecProductDataData random tests: branch inputDataRight(P,D)    |\n"\
        << "===============================================================================\n";
      /*
       * Test derived from Test 5.d
       */
      {// test 6.d scope
        int c=5, p=9, d1=3;

        DynRankView ConstructWithLabel(in_p_d, p, d1);
        DynRankView ConstructWithLabel(in_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(data_c_p, c, p);
        DynRankView ConstructWithLabel(datainv_c_p, c, p);
        DynRankView ConstructWithLabel(data_c_1, c, 1);
        DynRankView ConstructWithLabel(datainv_c_1, c, 1);
        DynRankView ConstructWithLabel(data_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(data_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(data_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(data_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(data_c_p_one, c, p);
        DynRankView ConstructWithLabel(data_c_1_one, c, 1);
        DynRankView ConstructWithLabel(out_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(outi_c_p_d, c, p, d1);
        /***********************************************************************************************
          *                          Constant diagonal tensor: inputData(C,P)                          *
          **********************************************************************************************/


        for (int i=0; i < p; i++)
          for (int m=0; m<d1; m++)
            in_p_d(i, m) = Teuchos::ScalarTraits<value_type>::random();

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
              data_c_p(i, j) = Teuchos::ScalarTraits<value_type>::random();
              datainv_c_p(i, j) = 1.0/data_c_p(i, j);
              data_c_p_one(i, j) = 1.0;
          }
          data_c_1(i, 0) = Teuchos::ScalarTraits<value_type>::random();
          datainv_c_1(i, 0) = 1.0/data_c_1(i, 0);
        }

        // Tensor values vary by point
        art::matvecProductDataData(in_c_p_d, data_c_p_one, in_p_d);
        art::matvecProductDataData(out_c_p_d, data_c_p, in_p_d);
        art::matvecProductDataData(outi_c_p_d, datainv_c_p, out_c_p_d);
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (13): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point
        art::matvecProductDataData(in_c_p_d, data_c_p_one, in_p_d);
        art::matvecProductDataData(out_c_p_d, data_c_1, in_p_d);
        art::matvecProductDataData(outi_c_p_d, datainv_c_1, out_c_p_d);
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (14): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *                       Non-constant diagonal tensor: inputData(C,P,D)                       *
          **********************************************************************************************/
    
        for (int i=0; i < p; i++)
          for (int m=0; m<d1; m++)
            in_p_d(i, m) = Teuchos::ScalarTraits<value_type>::random();

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++)
            for (int m=0; m<d1; m++) {
              data_c_p_d(i, j, m) = Teuchos::ScalarTraits<value_type>::random();
              datainv_c_p_d(i, j, m) = 1.0/data_c_p_d(i, j, m);
          }
          for (int m=0; m<d1; m++) {
            data_c_1_d(i, 0, m) = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_1_d(i, 0, m) = 1.0/data_c_1_d(i, 0, m);
          }
        }

        // Tensor values vary by point:
        art::matvecProductDataData(in_c_p_d, data_c_p_one, in_p_d);
        art::matvecProductDataData(out_c_p_d, data_c_p_d, in_p_d);
        art::matvecProductDataData(outi_c_p_d, datainv_c_p_d, out_c_p_d);
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (15): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point:
        art::matvecProductDataData(in_c_p_d, data_c_p_one, in_p_d);
        art::matvecProductDataData(out_c_p_d, data_c_1_d, in_p_d);
        art::matvecProductDataData(outi_c_p_d, datainv_c_1_d, out_c_p_d);
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (16): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *                              Full tensor: inputData(C,P,D,D)                               *
          **********************************************************************************************/

        for (int i=0; i < p; i++)
            for (int m=0; m<d1; m++)
              in_p_d(i, m) = Teuchos::ScalarTraits<value_type>::random();

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
            }

          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }

        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);


        // Tensor values vary by point: test "N" and "T" (no-transpose/transpose) options
        art::matvecProductDataData(in_c_p_d, data_c_p_one, in_p_d);
        art::matvecProductDataData(out_c_p_d, data_c_p_d_d, in_p_d);
        art::matvecProductDataData(outi_c_p_d, datainv_c_p_d_d, out_c_p_d);
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (17): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matvecProductDataData(in_c_p_d, data_c_p_one, in_p_d);
        art::matvecProductDataData(out_c_p_d, data_c_p_d_d, in_p_d, 't');
        art::matvecProductDataData(outi_c_p_d, datainv_c_p_d_d, out_c_p_d, 't');
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (18): check matrix inverse property, w/ double transpose\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point: test "N" and "T" (no-transpose/transpose) options
        art::matvecProductDataData(in_c_p_d, data_c_p_one, in_p_d);
        art::matvecProductDataData(out_c_p_d, data_c_1_d_d, in_p_d);
        art::matvecProductDataData(outi_c_p_d, datainv_c_1_d_d, out_c_p_d);
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (19): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matvecProductDataData(in_c_p_d, data_c_p_one, in_p_d);
        art::matvecProductDataData(out_c_p_d, data_c_1_d_d, in_p_d, 't');
        art::matvecProductDataData(outi_c_p_d, datainv_c_1_d_d, out_c_p_d, 't');
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (20): check matrix inverse property, w/ double transpose\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *             Full tensor: inputData(C,P,D,D) test inverse transpose                         *
          **********************************************************************************************/


        for (int i=0; i < p; i++)
            for (int m=0; m<d1; m++)
              in_p_d(i, m) = Teuchos::ScalarTraits<value_type>::random();

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
            }
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }

        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::transpose(datainvtrn_c_p_d_d, datainv_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);
        rst::transpose(datainvtrn_c_1_d_d, datainv_c_1_d_d);

        // Tensor values vary by point:
        art::matvecProductDataData(in_c_p_d, data_c_p_one, in_p_d);
        art::matvecProductDataData(out_c_p_d, data_c_p_d_d, in_p_d);
        art::matvecProductDataData(outi_c_p_d, datainvtrn_c_p_d_d, out_c_p_d, 't');
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (21): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point:
        art::matvecProductDataData(in_c_p_d, data_c_p_one, in_p_d);
        art::matvecProductDataData(out_c_p_d, data_c_1_d_d, in_p_d);
        art::matvecProductDataData(outi_c_p_d, datainvtrn_c_1_d_d, out_c_p_d, 't');
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matvecProductDataData (22): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
      }// test 6.d scope



      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 7.a: matmatProductDataField random tests: branch inputFields(C,F,P,D,D)|\n"\
        << "===============================================================================\n";
      {// Test 7.a scope
        int c=5, p=9, f=7, d1=3;

        DynRankView ConstructWithLabel(in_c_f_p_d_d, c, f, p, d1, d1);
        DynRankView ConstructWithLabel(data_c_p, c, p);
        DynRankView ConstructWithLabel(datainv_c_p, c, p);
        DynRankView ConstructWithLabel(data_c_1, c, 1);
        DynRankView ConstructWithLabel(datainv_c_1, c, 1);
        DynRankView ConstructWithLabel(data_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(data_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(data_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(data_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(out_c_f_p_d_d, c, f, p, d1, d1);
        DynRankView ConstructWithLabel(outi_c_f_p_d_d, c, f, p, d1, d1);
        /***********************************************************************************************
         *                          Constant diagonal tensor: inputData(C,P)                          *
         **********************************************************************************************/




        for (int i=0; i < c; i++) {
          for (int j=0; j < f; j++)
            for (int k=0; k < p; k++)
              for (int m=0; m<d1; m++)
                for (int n=0; n<d1; n++)
                  in_c_f_p_d_d(i, j, k, m, n)  = Teuchos::ScalarTraits<value_type>::random();
          for (int j=0; j < p; j++) {
              data_c_p(i, j) = Teuchos::ScalarTraits<value_type>::random();
              datainv_c_p(i, j) = 1.0/data_c_p(i, j);
          }
          data_c_1(i, 0) = Teuchos::ScalarTraits<value_type>::random();
          datainv_c_1(i, 0) = 1.0/data_c_1(i, 0);
        }

        // Tensor values vary by point:
        art::matmatProductDataField(out_c_f_p_d_d, data_c_p, in_c_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d);
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (1): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point:
        art::matmatProductDataField(out_c_f_p_d_d, data_c_1, in_c_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d);
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (2): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
         *                     Non-onstant diagonal tensor: inputData(C,P,D)                          *
         **********************************************************************************************/

        for (int i=0; i < c; i++) {
          for (int j=0; j < f; j++)
            for (int k=0; k < p; k++)
              for (int m=0; m<d1; m++)
                for (int n=0; n<d1; n++)
                  in_c_f_p_d_d(i, j, k, m, n)  = Teuchos::ScalarTraits<value_type>::random();
          for (int j=0; j < p; j++)
            for (int m=0; m<d1; m++) {
              data_c_p_d(i, j, m) = Teuchos::ScalarTraits<value_type>::random();
              datainv_c_p_d(i, j, m) = 1.0/data_c_p_d(i, j, m);
          }
          for (int m=0; m<d1; m++) {
            data_c_1_d(i, 0, m) = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_1_d(i, 0, m) = 1.0/data_c_1_d(i, 0, m);
          }
        }

        // Tensor values vary by point:
        art::matmatProductDataField(out_c_f_p_d_d, data_c_p_d, in_c_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_p_d, out_c_f_p_d_d);
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (3): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point
        art::matmatProductDataField(out_c_f_p_d_d, data_c_1_d, in_c_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_1_d, out_c_f_p_d_d);
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (4): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
         *                               Full tensor: inputData(C,P,D,D)                              *
         **********************************************************************************************/

        for (int i=0; i < c; i++) {
          for (int j=0; j < f; j++)
            for (int k=0; k < p; k++)
              for (int m=0; m<d1; m++)
                for (int n=0; n<d1; n++)
                  in_c_f_p_d_d(i, j, k, m, n)  = Teuchos::ScalarTraits<value_type>::random();
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
            }
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }

        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);

        // Tensor values vary by point: test "N" and "T" options (no transpose/transpose)
        art::matmatProductDataField(out_c_f_p_d_d, data_c_p_d_d, in_c_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_p_d_d, out_c_f_p_d_d);
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (5): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matmatProductDataField(out_c_f_p_d_d, data_c_p_d_d, in_c_f_p_d_d, 't');
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_p_d_d, out_c_f_p_d_d, 't');
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (6): check matrix inverse property, w/ double transpose\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point: test "N" and "T" options (no transpose/transpose)
        art::matmatProductDataField(out_c_f_p_d_d, data_c_1_d_d, in_c_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_1_d_d, out_c_f_p_d_d);
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (7): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matmatProductDataField(out_c_f_p_d_d, data_c_1_d_d, in_c_f_p_d_d,'t');
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_1_d_d, out_c_f_p_d_d, 't');
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (8): check matrix inverse property, w/ double transpose\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
         *             Full tensor: inputData(C,P,D,D) test inverse transpose                         *
         **********************************************************************************************/
    
        for (int i=0; i < c; i++) {
          for (int j=0; j < f; j++)
            for (int k=0; k < p; k++)
              for (int m=0; m<d1; m++)
                for (int n=0; n<d1; n++)
                  in_c_f_p_d_d(i, j, k, m, n)  = Teuchos::ScalarTraits<value_type>::random();
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
            }
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }

        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::transpose(datainvtrn_c_p_d_d, datainv_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);
        rst::transpose(datainvtrn_c_1_d_d, datainv_c_1_d_d);

    
        // Tensor values vary by point:
        art::matmatProductDataField(out_c_f_p_d_d, data_c_p_d_d, in_c_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainvtrn_c_p_d_d, out_c_f_p_d_d, 't');
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (9): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point
        art::matmatProductDataField(out_c_f_p_d_d, data_c_1_d_d, in_c_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainvtrn_c_1_d_d, out_c_f_p_d_d, 't');
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (10): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
      }// test 7.a scope



      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 7.b: matmatProductDataField random tests: branch inputFields(F,P,D,D)  |\n"\
        << "===============================================================================\n";
      // kj
      {// Test 7.b scope
        int c=5, p=9, f=7, d1=3;

        DynRankView ConstructWithLabel(in_f_p_d_d, f, p, d1, d1);
        DynRankView ConstructWithLabel(in_c_f_p_d_d, c, f, p, d1, d1);
        DynRankView ConstructWithLabel(data_c_p, c, p);
        DynRankView ConstructWithLabel(datainv_c_p, c, p);
        DynRankView ConstructWithLabel(data_c_1, c, 1);
        DynRankView ConstructWithLabel(datainv_c_1, c, 1);
        DynRankView ConstructWithLabel(data_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(data_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(data_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(data_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(data_c_p_one, c, p);
        DynRankView ConstructWithLabel(data_c_1_one, c, 1);
        DynRankView ConstructWithLabel(out_c_f_p_d_d, c, f, p, d1, d1);
        DynRankView ConstructWithLabel(outi_c_f_p_d_d, c, f, p, d1, d1);
        /***********************************************************************************************
          *                          Constant diagonal tensor: inputData(C,P)                          *
          **********************************************************************************************/
        for (int i=0; i < f; i++)
          for (int j=0; j < p; j++)
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                in_f_p_d_d(i, j, m, n)  = Teuchos::ScalarTraits<value_type>::random();
        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
            data_c_p(i, j) = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_p(i, j) = 1.0/data_c_p(i, j);
            data_c_p_one(i, j) = 1.0;
          }
          data_c_1(i, 0) = Teuchos::ScalarTraits<value_type>::random();
          datainv_c_1(i, 0) = 1.0/data_c_1(i, 0);
        }


        // Tensor values vary by point
        art::matmatProductDataField(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
        art::matmatProductDataField(out_c_f_p_d_d, data_c_p, in_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d);
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (11): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point
        art::matmatProductDataField(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
        art::matmatProductDataField(out_c_f_p_d_d, data_c_1, in_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d);
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (12): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *                       Non-constant diagonal tensor: inputData(C,P,D)                       *
          **********************************************************************************************/
        for (int i=0; i < f; i++)
          for (int j=0; j < p; j++)
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                in_f_p_d_d(i, j, m, n)  = Teuchos::ScalarTraits<value_type>::random();
        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++)
            for (int m=0; m<d1; m++) {
              data_c_p_d(i, j, m) = Teuchos::ScalarTraits<value_type>::random();
              datainv_c_p_d(i, j, m) = 1.0/data_c_p_d(i, j, m);
          }
          for (int m=0; m<d1; m++) {
            data_c_1_d(i, 0, m) = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_1_d(i, 0, m) = 1.0/data_c_1_d(i, 0, m);
          }
        }
    
        // Tensor values vary by point:
        art::matmatProductDataField(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
        art::matmatProductDataField(out_c_f_p_d_d, data_c_p_d, in_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_p_d, out_c_f_p_d_d);
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (13): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point:
        art::matmatProductDataField(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
        art::matmatProductDataField(out_c_f_p_d_d, data_c_1_d, in_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_1_d, out_c_f_p_d_d);
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (14): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *                              Full tensor: inputData(C,P,D,D)                               *
          **********************************************************************************************/

        for (int i=0; i < f; i++)
          for (int j=0; j < p; j++)
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                in_f_p_d_d(i, j, m, n)  = Teuchos::ScalarTraits<value_type>::random();
        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
            }
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }

        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);
    
        // Tensor values vary by point: test "N" and "T" (no-transpose/transpose) options
        art::matmatProductDataField(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
        art::matmatProductDataField(out_c_f_p_d_d, data_c_p_d_d, in_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_p_d_d, out_c_f_p_d_d);
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (15): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matmatProductDataField(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
        art::matmatProductDataField(out_c_f_p_d_d, data_c_p_d_d, in_f_p_d_d, 't');
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_p_d_d, out_c_f_p_d_d, 't');
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (16): check matrix inverse property, w/ double transpose\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point: test "N" and "T" (no-transpose/transpose) options
        art::matmatProductDataField(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
        art::matmatProductDataField(out_c_f_p_d_d, data_c_1_d_d, in_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_1_d_d, out_c_f_p_d_d);
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (17): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matmatProductDataField(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
        art::matmatProductDataField(out_c_f_p_d_d, data_c_1_d_d, in_f_p_d_d, 't');
        art::matmatProductDataField(outi_c_f_p_d_d, datainv_c_1_d_d, out_c_f_p_d_d, 't');
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (18): check matrix inverse property, w/ double transpose\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *             Full tensor: inputData(C,P,D,D) test inverse transpose                         *
          **********************************************************************************************/

        for (int i=0; i < f; i++)
          for (int j=0; j < p; j++)
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                in_f_p_d_d(i, j, m, n)  = Teuchos::ScalarTraits<value_type>::random();

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
            }
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }

        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::transpose(datainvtrn_c_p_d_d, datainv_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);
        rst::transpose(datainvtrn_c_1_d_d, datainv_c_1_d_d);


        // Tensor values vary by point:
        art::matmatProductDataField(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
        art::matmatProductDataField(out_c_f_p_d_d, data_c_p_d_d, in_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainvtrn_c_p_d_d, out_c_f_p_d_d, 't');
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (19): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point:
        art::matmatProductDataField(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
        art::matmatProductDataField(out_c_f_p_d_d, data_c_1_d_d, in_f_p_d_d);
        art::matmatProductDataField(outi_c_f_p_d_d, datainvtrn_c_1_d_d, out_c_f_p_d_d, 't');
        rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataField (20): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
      }// test 7.b scope



      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 8.a: matmatProductDataData random tests: branch inputDataRight(C,P,D,D)|\n"\
        << "===============================================================================\n";
      /*
       * Test derived from Test 7.a
       */
      {// test 8.a scope
        int c=5, p=9, d1=3;


        DynRankView ConstructWithLabel(in_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(out_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(outi_c_p_d_d, c, p, d1, d1);

        DynRankView ConstructWithLabel(data_c_p, c, p);
        DynRankView ConstructWithLabel(datainv_c_p, c, p);
        DynRankView ConstructWithLabel(data_c_1, c, 1);
        DynRankView ConstructWithLabel(datainv_c_1, c, 1);
        DynRankView ConstructWithLabel(data_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(data_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(data_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(data_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_1_d_d, c, 1, d1, d1);
       /***********************************************************************************************
        *                          Constant diagonal tensor: inputDataLeft(C,P)                      *
        **********************************************************************************************/

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
            data_c_p(i, j) = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_p(i, j) = 1.0/data_c_p(i, j);
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                in_c_p_d_d(i, j, m, n)  = Teuchos::ScalarTraits<value_type>::random();
          }
          data_c_1(i, 0) = Teuchos::ScalarTraits<value_type>::random();
          datainv_c_1(i, 0) = 1.0/data_c_1(i, 0);
        }

        // Tensor values vary by point:
        art::matmatProductDataData(out_c_p_d_d, data_c_p, in_c_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_p, out_c_p_d_d);
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (1): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point:
        art::matmatProductDataData(out_c_p_d_d, data_c_1, in_c_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_1, out_c_p_d_d);
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (2): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *                     Non-onstant diagonal tensor: inputDataLeft(C,P,D)                      *
          **********************************************************************************************/

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++)
            for (int m=0; m<d1; m++) {
              data_c_p_d(i, j, m) = Teuchos::ScalarTraits<value_type>::random();
              datainv_c_p_d(i, j, m) = 1.0/data_c_p_d(i, j, m);
              for (int n=0; n<d1; n++)
                in_c_p_d_d(i, j, m, n)  = Teuchos::ScalarTraits<value_type>::random();
          }
          for (int m=0; m<d1; m++) {
            data_c_1_d(i, 0, m) = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_1_d(i, 0, m) = 1.0/data_c_1_d(i, 0, m);
          }
        }

        // Tensor values vary by point:
        art::matmatProductDataData(out_c_p_d_d, data_c_p_d, in_c_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_p_d, out_c_p_d_d);
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (3): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point
        art::matmatProductDataData(out_c_p_d_d, data_c_1_d, in_c_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_1_d, out_c_p_d_d);
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (4): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *                            Full tensor: inputDataLeft(C,P,D,D)                             *
          **********************************************************************************************/

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++) {
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
                in_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
              }
            }

          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }
    
        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);
    
        // Tensor values vary by point: test "N" and "T" options (no transpose/transpose)
        art::matmatProductDataData(out_c_p_d_d, data_c_p_d_d, in_c_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_p_d_d, out_c_p_d_d);
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (5): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matmatProductDataData(out_c_p_d_d, data_c_p_d_d, in_c_p_d_d, 't');
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_p_d_d, out_c_p_d_d, 't');
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (6): check matrix inverse property, w/ double transpose\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point: test "N" and "T" options (no transpose/transpose)
        art::matmatProductDataData(out_c_p_d_d, data_c_1_d_d, in_c_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_1_d_d, out_c_p_d_d);
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (7): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matmatProductDataData(out_c_p_d_d, data_c_1_d_d, in_c_p_d_d, 't');
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_1_d_d, out_c_p_d_d, 't');
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (8): check matrix inverse property, w/ double transpose\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *             Full tensor: inputDataLeft(C,P,D,D) test inverse transpose                     *
          **********************************************************************************************/

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++) {
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
                in_c_p_d_d(i, j, m, n)  = Teuchos::ScalarTraits<value_type>::random();
              }
            }
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }

        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::transpose(datainvtrn_c_p_d_d, datainv_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);
        rst::transpose(datainvtrn_c_1_d_d, datainv_c_1_d_d);

        // Tensor values vary by point:
        art::matmatProductDataData(out_c_p_d_d, data_c_p_d_d, in_c_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainvtrn_c_p_d_d, out_c_p_d_d, 't');
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (9): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point
        art::matmatProductDataData(out_c_p_d_d, data_c_1_d_d, in_c_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainvtrn_c_1_d_d, out_c_p_d_d, 't');
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (10): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
      }// test 8.a scope
      *outStream \
        << "\n"
        << "===============================================================================\n"\
        << "| TEST 8.b: matmatProductDataData random tests: branch inputDataRight(P,D,D)  |\n"\
        << "===============================================================================\n";
      /*
       * Test derived from Test 7.b
       */
      {// test 8.b scope
        int c=5, p=9, d1=3;


        DynRankView ConstructWithLabel(in_p_d_d, p, d1, d1);
        DynRankView ConstructWithLabel(in_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(out_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(outi_c_p_d_d, c, p, d1, d1);

        DynRankView ConstructWithLabel(data_c_p, c, p);
        DynRankView ConstructWithLabel(datainv_c_p, c, p);
        DynRankView ConstructWithLabel(data_c_1, c, 1);
        DynRankView ConstructWithLabel(datainv_c_1, c, 1);
        DynRankView ConstructWithLabel(data_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d, c, p, d1);
        DynRankView ConstructWithLabel(data_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d, c, 1, d1);
        DynRankView ConstructWithLabel(data_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_p_d_d, c, p, d1, d1);
        DynRankView ConstructWithLabel(data_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainv_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(datainvtrn_c_1_d_d, c, 1, d1, d1);
        DynRankView ConstructWithLabel(data_c_p_one, c, p);
        DynRankView ConstructWithLabel(data_c_1_one, c, 1);
        /*********************************************************************************************
        *                          Constant diagonal tensor: inputData(C,P)                          *
        **********************************************************************************************/

        for (int i=0; i < p; i++)
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              in_p_d_d(i, m, n)  = Teuchos::ScalarTraits<value_type>::random();

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
            data_c_p(i, j) = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_p(i, j) = 1.0/data_c_p(i, j);
            data_c_p_one(i, j) = 1.0;
          }
          data_c_1(i, 0) = Teuchos::ScalarTraits<value_type>::random();
          datainv_c_1(i, 0) = 1.0/data_c_1(i, 0);
        }
        // Tensor values vary by point
        art::matmatProductDataData(in_c_p_d_d, data_c_p_one, in_p_d_d);
        art::matmatProductDataData(out_c_p_d_d, data_c_p, in_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_p, out_c_p_d_d);
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (11): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point
        art::matmatProductDataData(in_c_p_d_d, data_c_p_one, in_p_d_d);
        art::matmatProductDataData(out_c_p_d_d, data_c_1, in_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_1, out_c_p_d_d);
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (12): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *                       Non-constant diagonal tensor: inputData(C,P,D)                       *
          **********************************************************************************************/

        for (int i=0; i < p; i++)
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              in_p_d_d(i, m, n)  = Teuchos::ScalarTraits<value_type>::random();

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++)
            for (int m=0; m<d1; m++) {
              data_c_p_d(i, j, m) = Teuchos::ScalarTraits<value_type>::random();
              datainv_c_p_d(i, j, m) = 1.0/data_c_p_d(i, j, m);
          }
          for (int m=0; m<d1; m++) {
            data_c_1_d(i, 0, m) = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_1_d(i, 0, m) = 1.0/data_c_1_d(i, 0, m);
          }
        }

        // Tensor values vary by point:
        art::matmatProductDataData(in_c_p_d_d, data_c_p_one, in_p_d_d);
        art::matmatProductDataData(out_c_p_d_d, data_c_p_d, in_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_p_d, out_c_p_d_d);
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (13): check scalar inverse property\n\n";
          errorFlag = -1000;
        }

        // Tensor values do not vary by point:
        art::matmatProductDataData(in_c_p_d_d, data_c_p_one, in_p_d_d);
        art::matmatProductDataData(out_c_p_d_d, data_c_1_d, in_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_1_d, out_c_p_d_d);
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (14): check scalar inverse property\n\n";
          errorFlag = -1000;
        }

        /***********************************************************************************************
          *                              Full tensor: inputData(C,P,D,D)                               *
          **********************************************************************************************/

        for (int i=0; i < p; i++)
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              in_p_d_d(i, m, n) = Teuchos::ScalarTraits<value_type>::random();

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
            }

          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }

        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);

        // Tensor values vary by point: test "N" and "T" (no-transpose/transpose) options
        art::matmatProductDataData(in_c_p_d_d, data_c_p_one, in_p_d_d);
        art::matmatProductDataData(out_c_p_d_d, data_c_p_d_d, in_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_p_d_d, out_c_p_d_d);
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (15): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matmatProductDataData(in_c_p_d_d, data_c_p_one, in_p_d_d);
        art::matmatProductDataData(out_c_p_d_d, data_c_p_d_d, in_p_d_d, 't');
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_p_d_d, out_c_p_d_d, 't');
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (16): check matrix inverse property, w/ double transpose\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point: test "N" and "T" (no-transpose/transpose) options
        art::matmatProductDataData(in_c_p_d_d, data_c_p_one, in_p_d_d);
        art::matmatProductDataData(out_c_p_d_d, data_c_1_d_d, in_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_1_d_d, out_c_p_d_d);
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (17): check matrix inverse property\n\n";
          errorFlag = -1000;
        }
        art::matmatProductDataData(in_c_p_d_d, data_c_p_one, in_p_d_d);
        art::matmatProductDataData(out_c_p_d_d, data_c_1_d_d, in_p_d_d, 't');
        art::matmatProductDataData(outi_c_p_d_d, datainv_c_1_d_d, out_c_p_d_d, 't');
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (18): check matrix inverse property, w/ double transpose\n\n";
          errorFlag = -1000;
        }
        /***********************************************************************************************
          *             Full tensor: inputData(C,P,D,D) test inverse transpose                         *
          **********************************************************************************************/

        for (int i=0; i < p; i++)
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              in_p_d_d(i, m, n)  = Teuchos::ScalarTraits<value_type>::random();

        for (int i=0; i < c; i++) {
          for (int j=0; j < p; j++) {
            for (int m=0; m<d1; m++)
              for (int n=0; n<d1; n++)
                data_c_p_d_d(i, j, m, n) = Teuchos::ScalarTraits<value_type>::random();
            }
          for (int m=0; m<d1; m++)
            for (int n=0; n<d1; n++)
              data_c_1_d_d(i, 0, m, n) = Teuchos::ScalarTraits<value_type>::random();
        }

        rst::inverse(datainv_c_p_d_d, data_c_p_d_d);
        rst::transpose(datainvtrn_c_p_d_d, datainv_c_p_d_d);
        rst::inverse(datainv_c_1_d_d, data_c_1_d_d);
        rst::transpose(datainvtrn_c_1_d_d, datainv_c_1_d_d);

        // Tensor values vary by point:
        art::matmatProductDataData(in_c_p_d_d, data_c_p_one, in_p_d_d);
        art::matmatProductDataData(out_c_p_d_d, data_c_p_d_d, in_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainvtrn_c_p_d_d, out_c_p_d_d, 't');
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (19): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
        // Tensor values do not vary by point:
        art::matmatProductDataData(in_c_p_d_d, data_c_p_one, in_p_d_d);
        art::matmatProductDataData(out_c_p_d_d, data_c_1_d_d, in_p_d_d);
        art::matmatProductDataData(outi_c_p_d_d, datainvtrn_c_1_d_d, out_c_p_d_d, 't');
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT matmatProductDataData (20): check matrix inverse transpose property\n\n";
          errorFlag = -1000;
        }
      }// test 8.b scope
      }//try

      /*************************************************************************************************
        *                                      Finish test                                             *
        ************************************************************************************************/

      catch (std::logic_error err) {
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
  }
}
