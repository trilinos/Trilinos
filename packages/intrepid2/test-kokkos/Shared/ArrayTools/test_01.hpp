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

/** \file test_01.cpp
\brief  Unit test for the ArrayTools class.
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
    int ArrayTools_Test01(const bool verbose) {
      
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

      *outStream      \
			  << "===============================================================================\n" \
			  << "|                                                                             |\n" \
			  << "|                       Unit Test (ArrayTools)                                |\n" \
			  << "|                                                                             |\n" \
			  << "|     1) Array operations: contractions                                       |\n" \
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
      #ifdef HAVE_INTREPID2_DEBUG
  			art atools;
			#endif
      #define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

      int errorFlag = 0;

      *outStream                                \
        << "\n"
        << "===============================================================================\n" \
        << "| TEST 1: exceptions                                                   |\n" \
        << "===============================================================================\n";




  try{

    #ifdef HAVE_INTREPID2_DEBUG
        DynRankView ConstructWithLabel(a_2_2, 2, 2);
        DynRankView ConstructWithLabel(a_10_2, 10, 2);
        DynRankView ConstructWithLabel(a_10_3, 10, 3);
        DynRankView ConstructWithLabel(a_10_2_2, 10, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_3, 10, 2, 3);        
        DynRankView ConstructWithLabel(a_10_3_2, 10, 3, 2);
        DynRankView ConstructWithLabel(a_9_2_2, 9, 2, 2);

        std::cout << "Line: " << __LINE__ <<std::endl;
    
        *outStream << "-> contractFieldFieldScalar:\n";
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldScalar(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldScalar(a_2_2, a_10_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldScalar(a_2_2, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldScalar(a_10_2_2, a_9_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldScalar(a_10_2_2, a_10_2_2, a_10_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldScalar(a_9_2_2, a_10_3_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldScalar(a_10_3_2, a_10_2_2, a_10_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldScalar(a_10_2_2, a_10_2_2, a_10_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldScalar(a_10_2_3, a_10_2_2, a_10_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldScalar(a_10_2_3, a_10_2_2, a_10_3_2) );

        std::cout << "Line: " << __LINE__ <<std::endl;


        DynRankView ConstructWithLabel(a_10_2_2_2, 10, 2, 2, 2);
        std::cout << "Line: " << __LINE__ <<std::endl;
        DynRankView ConstructWithLabel(a_9_2_2_2, 9, 2, 2, 2);
        std::cout << "Line: " << __LINE__ <<std::endl;
        DynRankView ConstructWithLabel(a_10_3_2_2, 10, 3, 2, 2);
        std::cout << "Line: " << __LINE__ <<std::endl;
        DynRankView ConstructWithLabel(a_10_2_3_2, 10, 2, 3, 2);
        std::cout << "Line: " << __LINE__ <<std::endl;
        DynRankView ConstructWithLabel(a_10_2_2_3, 10, 2, 2, 3);

        std::cout << "Line: " << __LINE__ <<std::endl;


        *outStream << "-> contractFieldFieldVector:\n";
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldVector(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldVector(a_2_2, a_10_2_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldVector(a_2_2, a_10_2_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldVector(a_10_2_2, a_9_2_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldVector(a_10_2_2, a_10_2_2_2, a_10_2_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldVector(a_10_2_2, a_10_2_2_2, a_10_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldVector(a_9_2_2, a_10_2_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldVector(a_10_3_2, a_10_2_2_2, a_10_3_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldVector(a_10_2_2, a_10_2_2_2, a_10_3_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldVector(a_10_2_3, a_10_2_2_2, a_10_3_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldVector(a_10_2_3, a_10_2_2_2, a_10_3_2_2) );

        std::cout << "Line: " << __LINE__ <<std::endl;


        DynRankView ConstructWithLabel(a_10_2_2_2_2, 10, 2, 2, 2, 2);
        DynRankView ConstructWithLabel(a_9_2_2_2_2, 9, 2, 2, 2, 2);
        DynRankView ConstructWithLabel(a_10_3_2_2_2, 10, 3, 2, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_3_2_2, 10, 2, 3, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_2_3_2, 10, 2, 2, 3, 2);
        DynRankView ConstructWithLabel(a_10_2_2_2_3, 10, 2, 2, 2, 3);

        std::cout << "Line: " << __LINE__ <<std::endl;

        *outStream << "-> contractFieldFieldTensor:\n";
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldTensor(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldTensor(a_2_2, a_10_2_2_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldTensor(a_2_2, a_10_2_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldTensor(a_10_2_2, a_9_2_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldTensor(a_10_2_2, a_10_2_2_2_2, a_10_2_3_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldTensor(a_10_2_2, a_10_2_2_2_2, a_10_2_2_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldTensor(a_10_2_2, a_10_2_2_2_2, a_10_2_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldTensor(a_9_2_2, a_10_2_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldTensor(a_10_3_2, a_10_2_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldTensor(a_10_2_2, a_10_2_2_2_2, a_10_3_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldTensor(a_10_2_3, a_10_2_2_2_2, a_10_3_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractFieldFieldTensor(a_10_2_3, a_10_2_2_2_2, a_10_3_2_2_2) );

        std::cout << "Line: " << __LINE__ <<std::endl;


        DynRankView ConstructWithLabel(a_9_2, 9, 2);
        DynRankView ConstructWithLabel(a_10_1, 10, 1);

        std::cout << "Line: " << __LINE__ <<std::endl;

        *outStream << "-> contractDataFieldScalar:\n";
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldScalar(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldScalar(a_2_2, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldScalar(a_10_2_2, a_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldScalar(a_2_2, a_10_2, a_9_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldScalar(a_2_2, a_10_3, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldScalar(a_9_2, a_10_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldScalar(a_10_2, a_10_2, a_10_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldScalar(a_10_2, a_10_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldScalar(a_10_2, a_10_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldScalar(a_10_2, a_10_1, a_10_2_2) );

        std::cout << "Line: " << __LINE__ <<std::endl;


        DynRankView ConstructWithLabel(a_10_1_2, 10, 1, 2);
        DynRankView ConstructWithLabel(a_10_1_3, 10, 1, 3);

        std::cout << "Line: " << __LINE__ <<std::endl;

        *outStream << "-> contractDataFieldVector:\n";
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldVector(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldVector(a_2_2, a_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldVector(a_10_2_2, a_10_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldVector(a_2_2, a_10_2_2, a_9_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldVector(a_2_2, a_10_2_2, a_10_2_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldVector(a_2_2, a_10_2_2, a_10_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldVector(a_9_2, a_10_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldVector(a_10_2, a_10_2_2, a_10_3_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldVector(a_10_2, a_10_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldVector(a_10_2, a_10_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldVector(a_10_2, a_10_1_2, a_10_2_2_2) );

        std::cout << "Line: " << __LINE__ <<std::endl;


        DynRankView ConstructWithLabel(a_10_1_2_2, 10, 1, 2, 2);

        std::cout << "Line: " << __LINE__ <<std::endl;

        *outStream << "-> contractDataFieldTensor:\n";
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldTensor(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldTensor(a_2_2, a_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldTensor(a_10_2_2, a_10_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldTensor(a_2_2, a_10_2_2_2, a_9_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldTensor(a_2_2, a_10_2_2_2, a_10_2_3_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldTensor(a_2_2, a_10_2_2_2, a_10_2_2_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldTensor(a_2_2, a_10_2_2_2, a_10_2_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldTensor(a_9_2, a_10_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldTensor(a_10_2, a_10_2_2_2, a_10_3_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldTensor(a_10_2, a_10_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldTensor(a_10_2, a_10_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataFieldTensor(a_10_2, a_10_1_2_2, a_10_2_2_2_2) );

        std::cout << "Line: " << __LINE__ <<std::endl;


        DynRankView ConstructWithLabel(a_2, 2);
        DynRankView ConstructWithLabel(a_10, 10);

        std::cout << "Line: " << __LINE__ <<std::endl;

        *outStream << "-> contractDataDataScalar:\n";
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataScalar(a_2_2, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataScalar(a_2_2, a_10_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataScalar(a_2_2, a_10_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataScalar(a_2, a_9_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataScalar(a_2, a_10_2, a_10_3) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataScalar(a_2, a_10_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataScalar(a_10, a_10_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataScalar(a_10, a_10_2, a_10_2) );

        std::cout << "Line: " << __LINE__ <<std::endl;


        *outStream << "-> contractDataDataVector:\n";
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataVector(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataVector(a_2_2, a_10_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataVector(a_2_2, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataVector(a_10, a_9_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataVector(a_10, a_10_3_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataVector(a_10, a_10_2_3, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataVector(a_2, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataVector(a_10, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataVector(a_10, a_10_2_2, a_10_2_2) );

        std::cout << "Line: " << __LINE__ <<std::endl;


        *outStream << "-> contractDataDataTensor:\n";
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataTensor(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataTensor(a_2_2, a_10_2_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataTensor(a_2_2, a_10_2_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataTensor(a_10, a_9_2_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataTensor(a_10, a_10_2_2_2, a_10_3_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataTensor(a_10, a_10_2_2_2, a_10_2_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataTensor(a_10, a_10_2_2_2, a_10_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataTensor(a_2, a_10_2_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataTensor(a_10, a_10_2_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( atools.contractDataDataTensor(a_10, a_10_2_2_2, a_10_2_2_2) );

        std::cout << "Line: " << __LINE__ <<std::endl;


    #endif

      }
      catch (std::logic_error err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      };


      *outStream \
      << "\n"
      << "===============================================================================\n"\
      << "| TEST 2: correctness of math operations                                      |\n"\
      << "===============================================================================\n";

      outStream->precision(20);

      try {
          { // start scope
          *outStream << "\n************ Checking contractFieldFieldScalar ************\n";

          int c=5, p=9, l=3, r=7;

          DynRankView ConstructWithLabel(in_c_l_p, c, l, p);
          DynRankView ConstructWithLabel(in_c_r_p, c, r, p);
          DynRankView ConstructWithLabel(out1_c_l_r, c, l, r);
          DynRankView ConstructWithLabel(out2_c_l_r, c, l, r);
          value_type zero = INTREPID_TOL*10000.0;

          // fill with random numbers
          for (auto i=0;i<c;++i) {
            for (auto j=0;j<l;++j)
              for (auto k=0;k<p;++k)
                  in_c_l_p(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
            for (auto j=0;j<r;++j)
               for (auto k=0;k<p;++k)
                   in_c_r_p(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
          }

          art::contractFieldFieldScalar(out1_c_l_r, in_c_l_p, in_c_r_p);
          art::contractFieldFieldScalar(out2_c_l_r, in_c_l_p, in_c_r_p);
          rst::subtract(out1_c_l_r, out2_c_l_r);
          if (rst::Serial::vectorNorm(out1_c_l_r, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractFieldFieldScalar (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_r, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          // with sumInto:

          //fill with 2.0
          for (auto i=0;i<c;++i)
            for (auto j=0;j<l;++j)
              for (auto k=0;k<r;++k)
                out1_c_l_r(i,j,k) = out2_c_l_r(i,j,k) = 2.0;

          art::contractFieldFieldScalar(out1_c_l_r, in_c_l_p, in_c_r_p, true);
          art::contractFieldFieldScalar(out2_c_l_r, in_c_l_p, in_c_r_p, true);
          rst::subtract(out1_c_l_r, out2_c_l_r);
          if (rst::Serial::vectorNorm(out1_c_l_r, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractFieldFieldScalar (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_r, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          } // end scope

          { // start scope
          *outStream << "\n************ Checking contractFieldFieldVector ************\n";

          int c=5, p=9, l=3, r=7, d=13;

          DynRankView ConstructWithLabel(in_c_l_p_d, c, l, p, d);
          DynRankView ConstructWithLabel(in_c_r_p_d, c, r, p, d);
          DynRankView ConstructWithLabel(out1_c_l_r, c, l, r);
          DynRankView ConstructWithLabel(out2_c_l_r, c, l, r);
          value_type zero = INTREPID_TOL*10000.0;

          // fill with random numbers
          for (auto i=0;i<c;++i) {
            for (auto j=0;j<l;++j)
              for (auto k=0;k<p;++k)
                for (auto m=0;m<d;++m)
                  in_c_l_p_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
            for (auto j=0;j<r;++j)
              for (auto k=0;k<p;++k)
                for (auto m=0;m<d;++m)
                  in_c_r_p_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
          }


          art::contractFieldFieldVector(out1_c_l_r, in_c_l_p_d, in_c_r_p_d);
          art::contractFieldFieldVector(out2_c_l_r, in_c_l_p_d, in_c_r_p_d);

          rst::subtract(out1_c_l_r, out2_c_l_r);
          if (rst::Serial::vectorNorm(out1_c_l_r, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractFieldFieldVector (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_r, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }

          // with sumInto:

          //fill with 2.0
          for (auto i=0;i<c;++i)
            for (auto j=0;j<l;++j)
              for (auto k=0;k<r;++k)
                out1_c_l_r(i,j,k) = out2_c_l_r(i,j,k) = 2.0;

          art::contractFieldFieldVector(out1_c_l_r, in_c_l_p_d, in_c_r_p_d, true);
          art::contractFieldFieldVector(out2_c_l_r, in_c_l_p_d, in_c_r_p_d, true);

          rst::subtract(out1_c_l_r, out2_c_l_r);
          if (rst::Serial::vectorNorm(out1_c_l_r, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractFieldFieldVector (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_r, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          } // end scope

          { // start scope
          *outStream << "\n************ Checking contractFieldFieldTensor ************\n";

          int c=5, p=9, l=3, r=7, d1=13, d2=5;

          DynRankView ConstructWithLabel(in_c_l_p_d_d, c, l, p, d1, d2);
          DynRankView ConstructWithLabel(in_c_r_p_d_d, c, r, p, d1, d2);
          DynRankView ConstructWithLabel(out1_c_l_r, c, l, r);
          DynRankView ConstructWithLabel(out2_c_l_r, c, l, r);
          value_type zero = INTREPID_TOL*10000.0;

          // fill with random numbers
          for (auto i=0;i<c;++i) {
            for (auto j=0;j<l;++j)
              for (auto k=0;k<p;++k)
                for (auto m=0;m<d1;++m)
                  for (auto n=0;n<d2;++n)
                    in_c_l_p_d_d(i,j,k,m,n) = Teuchos::ScalarTraits<value_type>::random();
            for (auto j=0;j<r;++j)
              for (auto k=0;k<p;++k)
                for (auto m=0;m<d1;++m)
                  for (auto n=0;n<d2;++n)
                    in_c_r_p_d_d(i,j,k,m,n) = Teuchos::ScalarTraits<value_type>::random();
          }

          art::contractFieldFieldTensor(out1_c_l_r, in_c_l_p_d_d, in_c_r_p_d_d);
          art::contractFieldFieldTensor(out2_c_l_r, in_c_l_p_d_d, in_c_r_p_d_d);

          rst::subtract(out1_c_l_r, out2_c_l_r);
          if (rst::Serial::vectorNorm(out1_c_l_r, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractFieldFieldTensor (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_r, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }

          // with sumInto:

          //fill with 2.0
          for (auto i=0;i<c;++i)
            for (auto j=0;j<l;++j)
              for (auto k=0;k<r;++k)
                out1_c_l_r(i,j,k) = out2_c_l_r(i,j,k) = 2.0;

          art::contractFieldFieldTensor(out1_c_l_r, in_c_l_p_d_d, in_c_r_p_d_d, true);
          art::contractFieldFieldTensor(out2_c_l_r, in_c_l_p_d_d, in_c_r_p_d_d, true);

          rst::subtract(out1_c_l_r, out2_c_l_r);
          if (rst::Serial::vectorNorm(out1_c_l_r, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractFieldFieldTensor (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_r, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          } // end scope

          { // start scope
          *outStream << "\n************ Checking contractDataFieldScalar ************\n";

          int c=5, p=9, l=7;

          DynRankView ConstructWithLabel(in_c_l_p, c, l, p);
          DynRankView ConstructWithLabel(data_c_p, c, p);
          DynRankView ConstructWithLabel(data_c_1, c, 1);
          DynRankView ConstructWithLabel(out1_c_l, c, l);
          DynRankView ConstructWithLabel(out2_c_l, c, l);
          value_type zero = INTREPID_TOL*10000.0;


          // fill with random numbers
          for (auto i=0;i<c;++i) {
            for (auto j=0;j<l;++j)
              for (auto k=0;k<p;++k)
                in_c_l_p(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
            for (auto j=0;j<p;++j)
              data_c_p(i,j) = Teuchos::ScalarTraits<value_type>::random();
            data_c_1(i,0) = Teuchos::ScalarTraits<value_type>::random();
          }

          // nonconstant data
          art::contractDataFieldScalar(out1_c_l, data_c_p, in_c_l_p);
          art::contractDataFieldScalar(out2_c_l, data_c_p, in_c_l_p);
          rst::subtract(out1_c_l, out2_c_l);
          if (rst::Serial::vectorNorm(out1_c_l, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractDataFieldScalar (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          // constant data
          art::contractDataFieldScalar(out1_c_l, data_c_1, in_c_l_p);
          art::contractDataFieldScalar(out2_c_l, data_c_1, in_c_l_p);
          rst::subtract(out1_c_l, out2_c_l);
          if (rst::Serial::vectorNorm(out1_c_l, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractDataFieldScalar (2): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          // nonconstant data with sumInto
          for (auto i=0;i<c;++i)
            for (auto j=0;j<l;++j)
              out1_c_l(i,j) = out2_c_l(i,j) = 2.0;

          art::contractDataFieldScalar(out1_c_l, data_c_p, in_c_l_p, true);
          art::contractDataFieldScalar(out2_c_l, data_c_p, in_c_l_p, true);
          rst::subtract(out1_c_l, out2_c_l);
          if (rst::Serial::vectorNorm(out1_c_l, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractDataFieldScalar (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          } // end scope

          { // start scope
          *outStream << "\n************ Checking contractDataFieldVector ************\n";

          int c=5, p=9, l=7, d=3;

          DynRankView ConstructWithLabel(in_c_l_p_d, c, l, p, d);
          DynRankView ConstructWithLabel(data_c_p_d, c, p, d);
          DynRankView ConstructWithLabel(data_c_1_d, c, 1, d);
          DynRankView ConstructWithLabel(out1_c_l, c, l);
          DynRankView ConstructWithLabel(out2_c_l, c, l);
          value_type zero = INTREPID_TOL*10000.0;

          // fill with random numbers
          for (auto i=0;i<c;++i) {
            for (auto j=0;j<l;++j)
              for (auto k=0;k<p;++k)
                for (auto m=0;m<d;++m)
                  in_c_l_p_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
            for (auto j=0;j<p;++j)
              for (auto k=0;k<d;++k)
                data_c_p_d(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
            for (auto k=0;k<d;++k)
              data_c_1_d(i,0,k) = Teuchos::ScalarTraits<value_type>::random();
          }


          // nonconstant data
          art::contractDataFieldVector(out1_c_l, data_c_p_d, in_c_l_p_d);
          art::contractDataFieldVector(out2_c_l, data_c_p_d, in_c_l_p_d);
          rst::subtract(out1_c_l, out2_c_l);
          if (rst::Serial::vectorNorm(out1_c_l, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractDataFieldVector (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          // constant data
          art::contractDataFieldVector(out1_c_l, data_c_1_d, in_c_l_p_d);
          art::contractDataFieldVector(out2_c_l, data_c_1_d, in_c_l_p_d);
          rst::subtract(out1_c_l, out2_c_l);
          if (rst::Serial::vectorNorm(out1_c_l, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractDataFieldVector (2): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          // nonconstant data with sumInto

          //fill with 2.0
          value_type two=2.0;
          deep_copy(out1_c_l, two);
          for (auto i=0;i<c;++i)
            for (auto j=0;j<l;++j)
              out1_c_l(i,j) = out2_c_l(i,j) = 2.0;

          art::contractDataFieldVector(out1_c_l, data_c_p_d, in_c_l_p_d, true);
          art::contractDataFieldVector(out2_c_l, data_c_p_d, in_c_l_p_d, true);
          rst::subtract(out1_c_l, out2_c_l);
          if (rst::Serial::vectorNorm(out1_c_l,NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractDataFieldVector (3): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          } // end scope

          { // start scope
          *outStream << "\n************ Checking contractDataFieldTensor ************\n";

          int c=5, p=9, l=7, d1=3, d2=13;

          DynRankView ConstructWithLabel(in_c_l_p_d_d, c, l, p, d1, d2);
          DynRankView ConstructWithLabel(data_c_p_d_d, c, p, d1, d2);
          DynRankView ConstructWithLabel(data_c_1_d_d, c, 1, d1, d2);
          DynRankView ConstructWithLabel(out1_c_l, c, l);
          DynRankView ConstructWithLabel(out2_c_l, c, l);
          value_type zero = INTREPID_TOL*10000.0;


          // fill with random numbers
          for (auto i=0;i<c;++i) {
            for (auto j=0;j<l;++j)
              for (auto k=0;k<p;++k)
                for (auto m=0;m<d1;++m)
                  for (auto n=0;n<d2;++n)
                    in_c_l_p_d_d(i,j,k,m,n) = Teuchos::ScalarTraits<value_type>::random();
            for (auto j=0;j<p;++j)
              for (auto k=0;k<d1;++k)
                for (auto m=0;m<d2;++m)
                  data_c_p_d_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
            for (auto k=0;k<d1;++k)
              for (auto m=0;m<d2;++m)
                data_c_1_d_d(i,0,k,m) = Teuchos::ScalarTraits<value_type>::random();
          }

          // nonconstant data
          art::contractDataFieldTensor(out1_c_l, data_c_p_d_d, in_c_l_p_d_d);
          art::contractDataFieldTensor(out2_c_l, data_c_p_d_d, in_c_l_p_d_d);
          rst::subtract(out1_c_l, out2_c_l);
          if (rst::Serial::vectorNorm(out1_c_l, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractDataFieldTensor (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }

          // constant data
          art::contractDataFieldTensor(out1_c_l, data_c_1_d_d, in_c_l_p_d_d);
          art::contractDataFieldTensor(out2_c_l, data_c_1_d_d, in_c_l_p_d_d);
          rst::subtract(out1_c_l, out2_c_l);
          if (rst::Serial::vectorNorm(out1_c_l, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractDataFieldTensor (2): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }

          // nonconstant data with sumInto

          //fill with 2.0
          for (auto i=0;i<c;++i)
            for (auto j=0;j<l;++j)
              out1_c_l(i,j) = out2_c_l(i,j) = 2.0;

          art::contractDataFieldTensor(out1_c_l, data_c_p_d_d, in_c_l_p_d_d, true);
          art::contractDataFieldTensor(out2_c_l, data_c_p_d_d, in_c_l_p_d_d, true);
          rst::subtract(out1_c_l, out2_c_l);
          if (rst::Serial::vectorNorm(out1_c_l, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractDataFieldTensor (3): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          } // end scope

          { // start scope
          *outStream << "\n************ Checking contractDataDataScalar ************\n";

          int c=5, p=9;

          DynRankView ConstructWithLabel(inl_c_p, c, p);
          DynRankView ConstructWithLabel(inr_c_p, c, p);
          DynRankView ConstructWithLabel(out1_c, c);
          DynRankView ConstructWithLabel(out2_c, c);
          value_type zero = INTREPID_TOL*10000.0;

          // fill with random numbers
          for (auto i=0;i<c;++i)
            for (auto j=0;j<p;++j) {
              inl_c_p(i,j) = Teuchos::ScalarTraits<value_type>::random();
              inr_c_p(i,j) = Teuchos::ScalarTraits<value_type>::random();
            }

          art::contractDataDataScalar(out1_c, inl_c_p, inr_c_p);
          art::contractDataDataScalar(out2_c, inl_c_p, inr_c_p);
          rst::subtract(out1_c, out2_c);
          if (rst::Serial::vectorNorm(out1_c, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractDataDataScalar (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          // with sumInto:

          //fill with 2.0
          for (auto i=0;i<c;++i)
            out1_c(i) = out2_c(i) = 2.0;

          art::contractDataDataScalar(out1_c, inl_c_p, inr_c_p, true);
          art::contractDataDataScalar(out2_c, inl_c_p, inr_c_p, true);
          rst::subtract(out1_c, out2_c);
          if (rst::Serial::vectorNorm(out1_c, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractDataDataScalar (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          } // end scope

          { // start scope
          *outStream << "\n************ Checking contractDataDataVector ************\n";

          int c=5, p=9, d=13;

          DynRankView ConstructWithLabel(inl_c_p_d, c, p, d);
          DynRankView ConstructWithLabel(inr_c_p_d, c, p, d);
          DynRankView ConstructWithLabel(out1_c, c);
          DynRankView ConstructWithLabel(out2_c, c);
          value_type zero = INTREPID_TOL*10000.0;

          // fill with random numbers
          for (auto i=0;i<c;++i)
            for (auto j=0;j<p;++j)
              for (auto k=0;k<d;++k) {
              inl_c_p_d(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
              inr_c_p_d(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
            }

          art::contractDataDataVector(out1_c, inl_c_p_d, inr_c_p_d);
          art::contractDataDataVector(out2_c, inl_c_p_d, inr_c_p_d);

          rst::subtract(out1_c, out2_c);
          if (rst::Serial::vectorNorm(out1_c, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractDataDataVector (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }

          // with sumInto:

          //fill with 2.0
          for (auto i=0;i<c;++i)
            out1_c(i) = out2_c(i) = 2.0;

          art::contractDataDataVector(out1_c, inl_c_p_d, inr_c_p_d, true);
          art::contractDataDataVector(out2_c, inl_c_p_d, inr_c_p_d, true);

          rst::subtract(out1_c, out2_c);
          if (rst::Serial::vectorNorm(out1_c, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractDataDataVector (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          } // end scope

          { // start scope
          *outStream << "\n************ Checking contractDataDataTensor ************\n";

          int c=5, p=9, d1=13, d2=5;

          DynRankView ConstructWithLabel(inl_c_p_d_d, c, p, d1, d2);
          DynRankView ConstructWithLabel(inr_c_p_d_d, c, p, d1, d2);
          DynRankView ConstructWithLabel(out1_c, c);
          DynRankView ConstructWithLabel(out2_c, c);
          value_type zero = INTREPID_TOL*10000.0;

          // fill with random numbers
          for (auto i=0;i<c;++i)
            for (auto j=0;j<p;++j)
              for (auto k=0;k<d1;++k)
                for (auto m=0;m<d2;++m) {
                  inl_c_p_d_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
                  inr_c_p_d_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
                }

          art::contractDataDataTensor(out1_c, inl_c_p_d_d, inr_c_p_d_d);
          art::contractDataDataTensor(out2_c, inl_c_p_d_d, inr_c_p_d_d);

          rst::subtract(out1_c, out2_c);
          if (rst::Serial::vectorNorm(out1_c, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractDataDataTensor (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }

          // with sumInto:

          //fill with 2.0
          for (auto i=0;i<c;++i)
            out1_c[i] = out2_c[i] = 2.0;

          art::contractDataDataTensor(out1_c, inl_c_p_d_d, inr_c_p_d_d, true);
          art::contractDataDataTensor(out2_c, inl_c_p_d_d, inr_c_p_d_d, true);

          rst::subtract(out1_c, out2_c);
          if (rst::Serial::vectorNorm(out1_c, NORM_ONE) > zero) {
            *outStream << "\n\nINCORRECT contractDataDataTensor (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          } // end scope

          /******************************************/
          *outStream << "\n";
      }
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
