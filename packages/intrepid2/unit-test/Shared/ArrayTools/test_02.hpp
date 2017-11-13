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

/** \file test_02.hpp
\brief  Unit test for the scalar multiply operations of the ArrayTools class.
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
    int ArrayTools_Test02(const bool verbose) {
      
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
        << "|     1) Array operations: scalar multiply                                    |\n" \
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

        DynRankView ConstructWithLabel(a_10_2_2_2, 10, 2, 2, 2);
        DynRankView ConstructWithLabel(a_9_2_2_2, 9, 2, 2, 2);
        DynRankView ConstructWithLabel(a_10_3_2_2, 10, 3, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_3_2, 10, 2, 3, 2);
        DynRankView ConstructWithLabel(a_10_2_2_3, 10, 2, 2, 3);

        DynRankView ConstructWithLabel(a_10_2_2_2_2, 10, 2, 2, 2, 2);
        DynRankView ConstructWithLabel(a_9_2_2_2_2, 9, 2, 2, 2, 2);
        DynRankView ConstructWithLabel(a_10_3_2_2_2, 10, 3, 2, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_3_2_2, 10, 2, 3, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_2_3_2, 10, 2, 2, 3, 2);
        DynRankView ConstructWithLabel(a_10_2_2_2_3, 10, 2, 2, 2, 3);

        DynRankView ConstructWithLabel(a_9_2, 9, 2);
        DynRankView ConstructWithLabel(a_10_1, 10, 1);

        DynRankView ConstructWithLabel(a_10_1_2, 10, 1, 2);
        DynRankView ConstructWithLabel(a_10_1_3, 10, 1, 3);

        DynRankView ConstructWithLabel(a_10_1_2_2, 10, 1, 2, 2);

        DynRankView ConstructWithLabel(a_2_3_2_2, 2, 3, 2, 2);
        DynRankView ConstructWithLabel(a_2_2_2_2, 2, 2, 2, 2);
        DynRankView ConstructWithLabel(a_2_10, 2, 10);
        DynRankView ConstructWithLabel(a_2, 2);

        *outStream << "-> scalarMultiplyDataField:\n";
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_2_2, a_10_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_2_2, a_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_2, a_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_2, a_10_3, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_9_2_2_2_2, a_10_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_3_2_2_2, a_10_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_3_2_2, a_10_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_2_3_2, a_10_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_2_2_3, a_10_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_2_2_2, a_10_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_2_2_2, a_10_1, a_10_2_2_2_2) );
        //
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_2_2, a_10_2_2, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_2_2, a_2_2, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_2_2, a_2_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_2, a_10_2, a_2_10) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_2_2_2, a_9_2, a_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_3_2_2_2, a_10_2, a_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_3_2_2, a_10_2, a_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_2_3_2, a_10_2, a_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_2_2_3, a_10_2, a_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_2_2_2, a_10_2, a_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataField(a_10_2_2_2_2, a_10_1, a_2_2_2_2) );


        DynRankView ConstructWithLabel(a_2_2_2, 2, 2, 2);

        *outStream << "-> scalarMultiplyDataData:\n";
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_2_2, a_10_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_2, a_2_2, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_2_2, a_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_2_2, a_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_2_2, a_10_3, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_9_2_2_2, a_10_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_3_2_2, a_10_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_2_3_2, a_10_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_2_2_3, a_10_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_2_2_2, a_10_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_2_2_2, a_10_1, a_10_2_2_2) );
        //
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_2_2, a_10_2_2, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_2_2_2_2, a_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_2_2_2, a_2_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_2_2, a_2_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_2_2_2, a_9_2, a_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_3_2_2, a_10_2, a_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_2_3_2, a_10_2, a_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_2_2_3, a_10_2, a_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_2_2_2, a_10_2, a_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::scalarMultiplyDataData(a_10_2_2_2, a_10_1, a_2_2_2) );
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
        *outStream << "\n************ Checking scalarMultiplyDataField, reciprocal=false (branch 1) ************\n";

        int c=5, p=9, f=7, d1=7, d2=13;

        DynRankView ConstructWithLabel(in_c_f_p, c, f, p);
        DynRankView ConstructWithLabel(in_c_f_p_d, c, f, p, d1);
        DynRankView ConstructWithLabel(in_c_f_p_d_d, c, f, p, d1, d2);
        DynRankView ConstructWithLabel(data_c_p, c, p);
        DynRankView ConstructWithLabel(datainv_c_p, c, p);
        DynRankView ConstructWithLabel(data_c_1, c, 1);
        DynRankView ConstructWithLabel(datainv_c_1, c, 1);
        DynRankView ConstructWithLabel(out_c_f_p, c, f, p);
        DynRankView ConstructWithLabel(outi_c_f_p, c, f, p);
        DynRankView ConstructWithLabel(out_c_f_p_d, c, f, p, d1);
        DynRankView ConstructWithLabel(outi_c_f_p_d, c, f, p, d1);
        DynRankView ConstructWithLabel(out_c_f_p_d_d, c, f, p, d1, d2);
        DynRankView ConstructWithLabel(outi_c_f_p_d_d, c, f, p, d1, d2);


        // fill with random numbers
        for (auto i=0;i<c;++i) {
          for (auto j=0;j<f;++j)
            for (auto k=0;k<p;++k) {
              in_c_f_p(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
              for (auto m=0;m<d1;++m) {
                in_c_f_p_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
                for (auto n=0;n<d2;++n)
                  in_c_f_p_d_d(i,j,k,m,n) = Teuchos::ScalarTraits<value_type>::random();
              }
            }
          for (auto j=0;j<p;++j) {
            data_c_p(i,j)  = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_p(i,j) = 1./data_c_p(i,j);
          }
          data_c_1(i,0)  = Teuchos::ScalarTraits<value_type>::random();
          datainv_c_1(i,0) = 1./data_c_1(i,0);
        }





      art::scalarMultiplyDataField(out_c_f_p, data_c_p, in_c_f_p);
      art::scalarMultiplyDataField(outi_c_f_p, datainv_c_p, out_c_f_p);
      rst::subtract(outi_c_f_p, in_c_f_p);
      if (rst::Serial::vectorNorm(outi_c_f_p, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      art::scalarMultiplyDataField(out_c_f_p_d, data_c_p, in_c_f_p_d);
      art::scalarMultiplyDataField(outi_c_f_p_d, datainv_c_p, out_c_f_p_d);
      rst::subtract(outi_c_f_p_d, in_c_f_p_d);
      if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (2): check scalar inverse property" << rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_p, in_c_f_p_d_d);
      art::scalarMultiplyDataField(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_1, in_c_f_p_d_d);
      art::scalarMultiplyDataField(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
      deep_copy(in_c_f_p_d_d, 1.0); deep_copy(data_c_p,5.0); deep_copy(data_c_1, 5.0);

      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_p, in_c_f_p_d_d);
      if (std::abs(rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) - data_c_p(0,0)*in_c_f_p_d_d(0,0,0,0,0)*c*p*f*d1*d2) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (5): check result: "
                   << rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) << " != "
                   << data_c_p(0,0)*in_c_f_p_d_d(0,0,0,0,0)*c*f*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_1, in_c_f_p_d_d);
      if (std::abs(rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) - data_c_1(0,0)*in_c_f_p_d_d(0,0,0,0,0)*c*p*f*d1*d2) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (6): check result: "
                   << rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) << " != "
                   << data_c_1(0,0)*in_c_f_p_d_d(0,0)*c*f*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking scalarMultiplyDataField, reciprocal=false (branch 2) ************\n";

      int c=5, p=9, f=7, d1=7, d2=13;

				DynRankView ConstructWithLabel(in_f_p, f, p);
				DynRankView ConstructWithLabel(in_f_p_d, f, p, d1);
				DynRankView ConstructWithLabel(in_f_p_d_d, f, p, d1, d2);
				DynRankView ConstructWithLabel(in_c_f_p, c, f, p);
				DynRankView ConstructWithLabel(in_c_f_p_d, c, f, p, d1);
				DynRankView ConstructWithLabel(in_c_f_p_d_d, c, f, p, d1, d2);
				DynRankView ConstructWithLabel(data_c_p, c, p);
				DynRankView ConstructWithLabel(datainv_c_p, c, p);
				DynRankView ConstructWithLabel(data_c_1, c, 1);
				DynRankView ConstructWithLabel(datainv_c_1, c, 1);
				DynRankView ConstructWithLabel(data_c_p_one, c, p);
				DynRankView ConstructWithLabel(data_c_1_one, c, 1);
				DynRankView ConstructWithLabel(out_c_f_p, c, f, p);
				DynRankView ConstructWithLabel(outi_c_f_p, c, f, p);
				DynRankView ConstructWithLabel(out_c_f_p_d, c, f, p, d1);
				DynRankView ConstructWithLabel(outi_c_f_p_d, c, f, p, d1);
				DynRankView ConstructWithLabel(out_c_f_p_d_d, c, f, p, d1, d2);
				DynRankView ConstructWithLabel(outi_c_f_p_d_d, c, f, p, d1, d2);

        // fill with random numbers
        for (auto i=0;i<f;++i)
          for (auto j=0;j<p;++j) {
              in_f_p(i,j) = Teuchos::ScalarTraits<value_type>::random();
              for (auto m=0;m<d1;++m) {
                in_f_p_d(i,j,m) = Teuchos::ScalarTraits<value_type>::random();
                for (auto n=0;n<d2;++n)
                  in_f_p_d_d(i,j,m,n) = Teuchos::ScalarTraits<value_type>::random();
              }
            }

        for (auto i=0;i<c;++i) {
          for (auto j=0;j<f;++j)
            for (auto k=0;k<p;++k){
              in_c_f_p(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
              for (auto m=0;m<d1;++m) {
                in_c_f_p_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
                for (auto n=0;n<d2;++n)
                  in_c_f_p_d_d(i,j,k,m,n) = Teuchos::ScalarTraits<value_type>::random();
              }
            }
          for (auto j=0;j<p;++j) {
            data_c_p(i,j)  = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_p(i,j) = 1./data_c_p(i,j);
            data_c_p_one(i,j) = 1.;
          }
          data_c_1(i,0)  = Teuchos::ScalarTraits<value_type>::random();
          datainv_c_1(i,0) = 1./data_c_1(i,0);
          data_c_1_one(i,0) = 1.;
        }

      art::scalarMultiplyDataField(out_c_f_p, data_c_p, in_f_p);
      art::scalarMultiplyDataField(outi_c_f_p, datainv_c_p, out_c_f_p);
      art::scalarMultiplyDataField(in_c_f_p, data_c_p_one, in_f_p);
      rst::subtract(outi_c_f_p, in_c_f_p);
      if (rst::Serial::vectorNorm(outi_c_f_p, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField(out_c_f_p_d, data_c_p, in_f_p_d);
      art::scalarMultiplyDataField(outi_c_f_p_d, datainv_c_p, out_c_f_p_d);
      art::scalarMultiplyDataField(in_c_f_p_d, data_c_p_one, in_f_p_d);
      rst::subtract(outi_c_f_p_d, in_c_f_p_d);
      if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_p, in_f_p_d_d);
      art::scalarMultiplyDataField(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d);
      art::scalarMultiplyDataField(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_1, in_f_p_d_d);
      art::scalarMultiplyDataField(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d);
      art::scalarMultiplyDataField(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
      deep_copy(in_f_p_d_d, 1.0); deep_copy(data_c_p, 5.0); deep_copy(data_c_1, 5.0);

      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_p, in_f_p_d_d);
      if (std::abs(rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) - data_c_p(0,0)*in_f_p_d_d(0,0,0,0)*c*p*f*d1*d2) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (5): check result: "
                   << rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) << " != "
                   << data_c_p(0,0)*in_f_p_d_d(0,0,0,0)*c*f*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_1, in_f_p_d_d);
      if (std::abs(rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) - data_c_1(0,0)*in_f_p_d_d(0,0,0,0)*c*p*f*d1*d2) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (6): check result: "
                   << rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) << " != "
                   << data_c_1(0,0)*in_f_p_d_d(0,0,0,0)*c*f*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking scalarMultiplyDataField, reciprocal=true, i.e. division (branch 1) ************\n";

      int c=5, p=9, f=7, d1=7, d2=13;

				DynRankView ConstructWithLabel(in_c_f_p, c, f, p);
				DynRankView ConstructWithLabel(in_c_f_p_d, c, f, p, d1);
				DynRankView ConstructWithLabel(in_c_f_p_d_d, c, f, p, d1, d2);
				DynRankView ConstructWithLabel(data_c_p, c, p);
				DynRankView ConstructWithLabel(datainv_c_p, c, p);
				DynRankView ConstructWithLabel(data_c_1, c, 1);
				DynRankView ConstructWithLabel(datainv_c_1, c, 1);
				DynRankView ConstructWithLabel(out_c_f_p, c, f, p);
				DynRankView ConstructWithLabel(outi_c_f_p, c, f, p);
				DynRankView ConstructWithLabel(out_c_f_p_d, c, f, p, d1);
				DynRankView ConstructWithLabel(outi_c_f_p_d, c, f, p, d1);
				DynRankView ConstructWithLabel(out_c_f_p_d_d, c, f, p, d1, d2);
				DynRankView ConstructWithLabel(outi_c_f_p_d_d, c, f, p, d1, d2);


        // fill with random numbers
        for (auto i=0;i<c;++i) {
          for (auto j=0;j<f;++j)
            for (auto k=0;k<p;++k){
              in_c_f_p(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
              for (auto m=0;m<d1;++m) {
                in_c_f_p_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
                for (auto n=0;n<d2;++n)
                  in_c_f_p_d_d(i,j,k,m,n) = Teuchos::ScalarTraits<value_type>::random();
              }
            }
          for (auto j=0;j<p;++j) {
            data_c_p(i,j)  = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_p(i,j) = 1./data_c_p(i,j);
          }
          data_c_1(i,0)  = Teuchos::ScalarTraits<value_type>::random();
          datainv_c_1(i,0) = 1./data_c_1(i,0);
        }



      art::scalarMultiplyDataField(out_c_f_p, data_c_p, in_c_f_p, true);
      art::scalarMultiplyDataField(outi_c_f_p, datainv_c_p, out_c_f_p, true);
      rst::subtract(outi_c_f_p, in_c_f_p);
      if (rst::Serial::vectorNorm(outi_c_f_p, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField(out_c_f_p_d, data_c_p, in_c_f_p_d, true);
      art::scalarMultiplyDataField(outi_c_f_p_d, datainv_c_p, out_c_f_p_d, true);
      rst::subtract(outi_c_f_p_d, in_c_f_p_d);
      if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_p, in_c_f_p_d_d, true);
      art::scalarMultiplyDataField(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d, true);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_1, in_c_f_p_d_d, true);
      art::scalarMultiplyDataField(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d, true);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
      deep_copy(in_c_f_p_d_d, 1.0); deep_copy(data_c_p, 5.0); deep_copy(data_c_1, 5.0);

      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_p, in_c_f_p_d_d, true);
      if (std::abs(rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) - (1.0/data_c_p(0,0))*in_c_f_p_d_d(0,0,0,0,0)*c*p*f*d1*d2)/rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (5): check result: "
                   << rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) << " != "
                   << (1.0/data_c_p(0,0))*in_c_f_p_d_d(0,0,0,0,0)*c*p*f*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_1, in_c_f_p_d_d, true);
      if (std::abs(rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) - (1.0/data_c_1(0,0))*in_c_f_p_d_d(0,0,0,0,0)*c*p*f*d1*d2)/rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (6): check result: "
                   << rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) << " != "
                   << (1.0/data_c_p(0,0))*in_c_f_p_d_d(0,0,0,0,0)*c*p*f*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking scalarMultiplyDataField, reciprocal=true, i.e. division (branch 2) ************\n";

      int c=5, p=9, f=7, d1=7, d2=13;

				DynRankView ConstructWithLabel(in_f_p, f, p);
				DynRankView ConstructWithLabel(in_f_p_d, f, p, d1);
				DynRankView ConstructWithLabel(in_f_p_d_d, f, p, d1, d2);
				DynRankView ConstructWithLabel(in_c_f_p, c, f, p);
				DynRankView ConstructWithLabel(in_c_f_p_d, c, f, p, d1);
				DynRankView ConstructWithLabel(in_c_f_p_d_d, c, f, p, d1, d2);
				DynRankView ConstructWithLabel(data_c_p, c, p);
				DynRankView ConstructWithLabel(datainv_c_p, c, p);
				DynRankView ConstructWithLabel(data_c_1, c, 1);
				DynRankView ConstructWithLabel(datainv_c_1, c, 1);
				DynRankView ConstructWithLabel(data_c_p_one, c, p);
				DynRankView ConstructWithLabel(data_c_1_one, c, 1);
				DynRankView ConstructWithLabel(out_c_f_p, c, f, p);
				DynRankView ConstructWithLabel(outi_c_f_p, c, f, p);
				DynRankView ConstructWithLabel(out_c_f_p_d, c, f, p, d1);
				DynRankView ConstructWithLabel(outi_c_f_p_d, c, f, p, d1);
				DynRankView ConstructWithLabel(out_c_f_p_d_d, c, f, p, d1, d2);
				DynRankView ConstructWithLabel(outi_c_f_p_d_d, c, f, p, d1, d2);


        // fill with random numbers
        for (auto i=0;i<f;++i)
          for (auto j=0;j<p;++j) {
              in_f_p(i,j) = Teuchos::ScalarTraits<value_type>::random();
              for (auto m=0;m<d1;++m) {
                in_f_p_d(i,j,m) = Teuchos::ScalarTraits<value_type>::random();
                for (auto n=0;n<d2;++n)
                  in_f_p_d_d(i,j,m,n) = Teuchos::ScalarTraits<value_type>::random();
              }
            }

        for (auto i=0;i<c;++i) {
          for (auto j=0;j<f;++j)
            for (auto k=0;k<p;++k){
              in_c_f_p(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
              for (auto m=0;m<d1;++m) {
                in_c_f_p_d(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
                for (auto n=0;n<d2;++n)
                  in_c_f_p_d_d(i,j,k,m,n) = Teuchos::ScalarTraits<value_type>::random();
              }
            }
          for (auto j=0;j<p;++j) {
            data_c_p(i,j)  = Teuchos::ScalarTraits<value_type>::random();
            datainv_c_p(i,j) = 1./data_c_p(i,j);
            data_c_p_one(i,j) = 1.;
          }
          data_c_1(i,0)  = Teuchos::ScalarTraits<value_type>::random();
          datainv_c_1(i,0) = 1./data_c_1(i,0);
          data_c_1_one(i,0) = 1.;
        }

      art::scalarMultiplyDataField(out_c_f_p, data_c_p, in_f_p, true);
      art::scalarMultiplyDataField(outi_c_f_p, datainv_c_p, out_c_f_p, true);
      art::scalarMultiplyDataField(in_c_f_p, data_c_p_one, in_f_p);
      rst::subtract(outi_c_f_p, in_c_f_p);
      if (rst::Serial::vectorNorm(outi_c_f_p, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField(out_c_f_p_d, data_c_p, in_f_p_d, true);
      art::scalarMultiplyDataField(outi_c_f_p_d, datainv_c_p, out_c_f_p_d, true);
      art::scalarMultiplyDataField(in_c_f_p_d, data_c_p_one, in_f_p_d);
      rst::subtract(outi_c_f_p_d, in_c_f_p_d);
      if (rst::Serial::vectorNorm(outi_c_f_p_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_p, in_f_p_d_d, true);
      art::scalarMultiplyDataField(outi_c_f_p_d_d, datainv_c_p, out_c_f_p_d_d, true);
      art::scalarMultiplyDataField(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_1, in_f_p_d_d, true);
      art::scalarMultiplyDataField(outi_c_f_p_d_d, datainv_c_1, out_c_f_p_d_d, true);
      art::scalarMultiplyDataField(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
      rst::subtract(outi_c_f_p_d_d, in_c_f_p_d_d);
      if (rst::Serial::vectorNorm(outi_c_f_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
      deep_copy(in_f_p_d_d, 1.0); deep_copy(data_c_p, 5.0); deep_copy(data_c_1, 5.0);

      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_p, in_f_p_d_d, true);
      if (std::abs(rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) - (1.0/data_c_p(0,0))*in_f_p_d_d(0,0,0,0)*c*p*f*d1*d2)/rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (5): check result: "
                   << rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) << " != "
                   << (1.0/data_c_p(0,0))*in_f_p_d_d(0,0,0,0)*c*p*f*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataField(out_c_f_p_d_d, data_c_1, in_f_p_d_d, true);
      if (std::abs(rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) - (1.0/data_c_1(0,0))*in_f_p_d_d(0,0,0,0)*c*p*f*d1*d2)/rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataField (6): check result: "
                   << rst::Serial::vectorNorm(out_c_f_p_d_d, NORM_ONE) << " != "
                   << (1.0/data_c_1(0,0))*in_f_p_d_d(0,0,0,0)*c*p*f*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope




      { // start scope
      *outStream << "\n************ Checking scalarMultiplyDataData, reciprocal=false (branch 1) ************\n";

      int c=5, p=9, d1=7, d2=13;

				DynRankView ConstructWithLabel(in_c_p, c, p);
				DynRankView ConstructWithLabel(in_c_p_d, c, p, d1);
				DynRankView ConstructWithLabel(in_c_p_d_d, c, p, d1, d2);
				DynRankView ConstructWithLabel(data_c_p, c, p);
				DynRankView ConstructWithLabel(datainv_c_p, c, p);
				DynRankView ConstructWithLabel(data_c_1, c, 1);
				DynRankView ConstructWithLabel(datainv_c_1, c, 1);
				DynRankView ConstructWithLabel(out_c_p, c, p);
				DynRankView ConstructWithLabel(outi_c_p, c, p);
				DynRankView ConstructWithLabel(out_c_p_d, c, p, d1);
				DynRankView ConstructWithLabel(outi_c_p_d, c, p, d1);
				DynRankView ConstructWithLabel(out_c_p_d_d, c, p, d1, d2);
				DynRankView ConstructWithLabel(outi_c_p_d_d, c, p, d1, d2);


        for (auto i=0;i<c;++i) {
          for (auto j=0;j<p;++j) {
              in_c_p(i,j) = Teuchos::ScalarTraits<value_type>::random();
              data_c_p(i,j)  = Teuchos::ScalarTraits<value_type>::random();
              datainv_c_p(i,j) = 1./data_c_p(i,j);
              for (auto m=0;m<d1;++m) {
                in_c_p_d(i,j,m) = Teuchos::ScalarTraits<value_type>::random();
                for (auto n=0;n<d2;++n)
                  in_c_p_d_d(i,j,m,n) = Teuchos::ScalarTraits<value_type>::random();
              }
            }
          data_c_1(i,0)  = Teuchos::ScalarTraits<value_type>::random();
          datainv_c_1(i,0) = 1./data_c_1(i,0);
        }


      art::scalarMultiplyDataData(out_c_p, data_c_p, in_c_p);
      art::scalarMultiplyDataData(outi_c_p, datainv_c_p, out_c_p);
      rst::subtract(outi_c_p, in_c_p);
      if (rst::Serial::vectorNorm(outi_c_p, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData(out_c_p_d, data_c_p, in_c_p_d);
      art::scalarMultiplyDataData(outi_c_p_d, datainv_c_p, out_c_p_d);
      rst::subtract(outi_c_p_d, in_c_p_d);
      if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData(out_c_p_d_d, data_c_p, in_c_p_d_d);
      art::scalarMultiplyDataData(outi_c_p_d_d, datainv_c_p, out_c_p_d_d);
      rst::subtract(outi_c_p_d_d, in_c_p_d_d);
      if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData(out_c_p_d_d, data_c_1, in_c_p_d_d);
      art::scalarMultiplyDataData(outi_c_p_d_d, datainv_c_1, out_c_p_d_d);
      rst::subtract(outi_c_p_d_d, in_c_p_d_d);
      if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
      deep_copy(in_c_p_d_d, 1.0); deep_copy(data_c_p, 5.0); deep_copy(data_c_1, 5.0);

      art::scalarMultiplyDataData(out_c_p_d_d, data_c_p, in_c_p_d_d);
      if (std::abs(rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) - data_c_p(0,0)*in_c_p_d_d(0,0,0,0)*c*p*d1*d2) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (5): check result: "
                   << rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                   << data_c_p(0,0)*in_c_p_d_d(0,0,0,0)*c*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData(out_c_p_d_d, data_c_1, in_c_p_d_d);
      if (std::abs(rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) - data_c_1(0,0)*in_c_p_d_d(0,0,0,0)*c*p*d1*d2) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (6): check result: "
                   << rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                   << data_c_p(0,0)*in_c_p_d_d(0,0,0,0)*c*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking scalarMultiplyDataData, reciprocal=false (branch 2) ************\n";

      int c=5, p=9, d1=7, d2=13;

				DynRankView ConstructWithLabel(in_p, p);
				DynRankView ConstructWithLabel(in_p_d, p, d1);
				DynRankView ConstructWithLabel(in_p_d_d, p, d1, d2);
				DynRankView ConstructWithLabel(in_c_p, c, p);
				DynRankView ConstructWithLabel(in_c_p_d, c, p, d1);
				DynRankView ConstructWithLabel(in_c_p_d_d, c, p, d1, d2);
				DynRankView ConstructWithLabel(data_c_p, c, p);
				DynRankView ConstructWithLabel(datainv_c_p, c, p);
				DynRankView ConstructWithLabel(data_c_1, c, 1);
				DynRankView ConstructWithLabel(datainv_c_1, c, 1);
				DynRankView ConstructWithLabel(data_c_p_one, c, p);
				DynRankView ConstructWithLabel(data_c_1_one, c, 1);
				DynRankView ConstructWithLabel(out_c_p, c, p);
				DynRankView ConstructWithLabel(outi_c_p, c, p);
				DynRankView ConstructWithLabel(out_c_p_d, c, p, d1);
				DynRankView ConstructWithLabel(outi_c_p_d, c, p, d1);
				DynRankView ConstructWithLabel(out_c_p_d_d, c, p, d1, d2);
				DynRankView ConstructWithLabel(outi_c_p_d_d, c, p, d1, d2);

        // fill with random numbers
        for (auto i=0;i<c;++i) {
          in_p(i) = Teuchos::ScalarTraits<value_type>::random();
          for (auto m=0;m<d1;++m) {
            in_p_d(i,m) = Teuchos::ScalarTraits<value_type>::random();
            for (auto n=0;n<d2;++n)
              in_p_d_d(i,m,n) = Teuchos::ScalarTraits<value_type>::random();
          }
        }

        for (auto i=0;i<c;++i) {
          for (auto j=0;j<p;++j) {
              in_c_p(i,j) = Teuchos::ScalarTraits<value_type>::random();
              data_c_p(i,j)  = Teuchos::ScalarTraits<value_type>::random();
              datainv_c_p(i,j) = 1./data_c_p(i,j);
              data_c_p_one(i,j) = 1.0;
              for (auto m=0;m<d1;++m) {
                in_c_p_d(i,j,m) = Teuchos::ScalarTraits<value_type>::random();
                for (auto n=0;n<d2;++n)
                  in_c_p_d_d(i,j,m,n) = Teuchos::ScalarTraits<value_type>::random();
              }
            }
          data_c_1(i,0)  = Teuchos::ScalarTraits<value_type>::random();
          datainv_c_1(i,0) = 1./data_c_1(i,0);
          data_c_1_one(i,0) = 1.0;
        }

      art::scalarMultiplyDataData(out_c_p, data_c_p, in_p);
      art::scalarMultiplyDataData(outi_c_p, datainv_c_p, out_c_p);
      art::scalarMultiplyDataData(in_c_p, data_c_p_one, in_p);
      rst::subtract(outi_c_p, in_c_p);
      if (rst::Serial::vectorNorm(outi_c_p, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData(out_c_p_d, data_c_p, in_p_d);
      art::scalarMultiplyDataData(outi_c_p_d, datainv_c_p, out_c_p_d);
      art::scalarMultiplyDataData(in_c_p_d, data_c_p_one, in_p_d);
      rst::subtract(outi_c_p_d, in_c_p_d);
      if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData(out_c_p_d_d, data_c_p, in_p_d_d);
      art::scalarMultiplyDataData(outi_c_p_d_d, datainv_c_p, out_c_p_d_d);
      art::scalarMultiplyDataData(in_c_p_d_d, data_c_p_one, in_p_d_d);
      rst::subtract(outi_c_p_d_d, in_c_p_d_d);
      if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData(out_c_p_d_d, data_c_1, in_p_d_d);
      art::scalarMultiplyDataData(outi_c_p_d_d, datainv_c_1, out_c_p_d_d);
      art::scalarMultiplyDataData(in_c_p_d_d, data_c_p_one, in_p_d_d);
      rst::subtract(outi_c_p_d_d, in_c_p_d_d);
      if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
      deep_copy(in_p_d_d, 1.0); deep_copy(data_c_p,5.0); deep_copy(data_c_1, 5.0);

      art::scalarMultiplyDataData(out_c_p_d_d, data_c_p, in_p_d_d);
      if (std::abs(rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) - data_c_p(0,0)*in_p_d_d(0,0,0)*c*p*d1*d2) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (5): check result: "
                   << rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                   << data_c_p(0,0)*in_p_d_d(0,0,0)*c*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData(out_c_p_d_d, data_c_1, in_p_d_d);
      if (std::abs(rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) - data_c_1(0,0)*in_p_d_d(0,0,0)*c*p*d1*d2) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (6): check result: "
                   << rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                   << data_c_1(0,0)*in_p_d_d(0,0,0)*c*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking scalarMultiplyDataData, reciprocal=true, i.e. division (branch 1) ************\n";

      int c=5, p=9, d1=7, d2=13;

				DynRankView ConstructWithLabel(in_c_p, c, p);
				DynRankView ConstructWithLabel(in_c_p_d, c, p, d1);
				DynRankView ConstructWithLabel(in_c_p_d_d, c, p, d1, d2);
				DynRankView ConstructWithLabel(data_c_p, c, p);
				DynRankView ConstructWithLabel(datainv_c_p, c, p);
				DynRankView ConstructWithLabel(data_c_1, c, 1);
				DynRankView ConstructWithLabel(datainv_c_1, c, 1);
				DynRankView ConstructWithLabel(out_c_p, c, p);
				DynRankView ConstructWithLabel(outi_c_p, c, p);
				DynRankView ConstructWithLabel(out_c_p_d, c, p, d1);
				DynRankView ConstructWithLabel(outi_c_p_d, c, p, d1);
				DynRankView ConstructWithLabel(out_c_p_d_d, c, p, d1, d2);
				DynRankView ConstructWithLabel(outi_c_p_d_d, c, p, d1, d2);

        // fill with random numbers
        for (auto i=0;i<c;++i) {
          for (auto j=0;j<p;++j) {
              in_c_p(i,j) = Teuchos::ScalarTraits<value_type>::random();
              data_c_p(i,j)  = Teuchos::ScalarTraits<value_type>::random();
              datainv_c_p(i,j) = 1./data_c_p(i,j);
              for (auto m=0;m<d1;++m) {
                in_c_p_d(i,j,m) = Teuchos::ScalarTraits<value_type>::random();
                for (auto n=0;n<d2;++n)
                  in_c_p_d_d(i,j,m,n) = Teuchos::ScalarTraits<value_type>::random();
              }
            }
          data_c_1(i,0)  = Teuchos::ScalarTraits<value_type>::random();
          datainv_c_1(i,0) = 1./data_c_1(i,0);
        }

      art::scalarMultiplyDataData(out_c_p, data_c_p, in_c_p, true);
      art::scalarMultiplyDataData(outi_c_p, datainv_c_p, out_c_p, true);
      rst::subtract(outi_c_p, in_c_p);
      if (rst::Serial::vectorNorm(outi_c_p, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (1): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData(out_c_p_d, data_c_p, in_c_p_d, true);
      art::scalarMultiplyDataData(outi_c_p_d, datainv_c_p, out_c_p_d, true);
      rst::subtract(outi_c_p_d, in_c_p_d);
      if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (2): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData(out_c_p_d_d, data_c_p, in_c_p_d_d, true);
      art::scalarMultiplyDataData(outi_c_p_d_d, datainv_c_p, out_c_p_d_d, true);
      rst::subtract(outi_c_p_d_d, in_c_p_d_d);
      if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (3): check scalar inverse property\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData(out_c_p_d_d, data_c_1, in_c_p_d_d, true);
      art::scalarMultiplyDataData(outi_c_p_d_d, datainv_c_1, out_c_p_d_d, true);
      rst::subtract(outi_c_p_d_d, in_c_p_d_d);
      if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (4): check scalar inverse property\n\n";
        errorFlag = -1000;
      }

      // fill with constants
      deep_copy(in_c_p_d_d, 1.0); deep_copy(data_c_p,5.0); deep_copy(data_c_1, 5.0);

      art::scalarMultiplyDataData(out_c_p_d_d, data_c_p, in_c_p_d_d, true);
      if (std::abs(rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) - (1.0/data_c_p(0,0))*in_c_p_d_d(0,0,0,0)*c*p*d1*d2)/rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (5): check result: "
                   << rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                   << (1.0/data_c_p(0,0))*in_c_p_d_d(0,0,0,0)*c*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      art::scalarMultiplyDataData(out_c_p_d_d, data_c_1, in_c_p_d_d, true);
      if (std::abs(rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) - (1.0/data_c_1(0,0))*in_c_p_d_d(0,0,0,0)*c*p*d1*d2)/rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) > tol) {
        *outStream << "\n\nINCORRECT scalarMultiplyDataData (6): check result: "
                   << rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                   << (1.0/data_c_p(0,0))*in_c_p_d_d(0,0,0,0)*c*p*d1*d2 << "\n\n";
        errorFlag = -1000;
      }
      } // end scope

      { // start scope
      *outStream << "\n************ Checking scalarMultiplyDataData, reciprocal=true, i.e. division (branch 2) ************\n";

      int c=5, p=9, d1=7, d2=13;

				DynRankView ConstructWithLabel(in_p, p);
				DynRankView ConstructWithLabel(in_p_d, p, d1);
				DynRankView ConstructWithLabel(in_p_d_d, p, d1, d2);
				DynRankView ConstructWithLabel(in_c_p, c, p);
				DynRankView ConstructWithLabel(in_c_p_d, c, p, d1);
				DynRankView ConstructWithLabel(in_c_p_d_d, c, p, d1, d2);
				DynRankView ConstructWithLabel(data_c_p, c, p);
				DynRankView ConstructWithLabel(datainv_c_p, c, p);
				DynRankView ConstructWithLabel(data_c_1, c, 1);
				DynRankView ConstructWithLabel(datainv_c_1, c, 1);
				DynRankView ConstructWithLabel(data_c_p_one, c, p);
				DynRankView ConstructWithLabel(data_c_1_one, c, 1);
				DynRankView ConstructWithLabel(out_c_p, c, p);
				DynRankView ConstructWithLabel(outi_c_p, c, p);
				DynRankView ConstructWithLabel(out_c_p_d, c, p, d1);
				DynRankView ConstructWithLabel(outi_c_p_d, c, p, d1);
				DynRankView ConstructWithLabel(out_c_p_d_d, c, p, d1, d2);
				DynRankView ConstructWithLabel(outi_c_p_d_d, c, p, d1, d2);

        // fill with random numbers
        for (auto i=0;i<p;++i) {
          in_p(i) = Teuchos::ScalarTraits<value_type>::random();
          for (auto m=0;m<d1;++m) {
            in_p_d(i,m) = Teuchos::ScalarTraits<value_type>::random();
            for (auto n=0;n<d2;++n)
              in_p_d_d(i,m,n) = Teuchos::ScalarTraits<value_type>::random();
            }
          }
        for (auto i=0;i<c;++i) {
          for (auto j=0;j<p;++j) {
              in_c_p(i,j) = Teuchos::ScalarTraits<value_type>::random();
              data_c_p(i,j)  = Teuchos::ScalarTraits<value_type>::random();
              datainv_c_p(i,j) = 1./data_c_p(i,j);
              data_c_p_one(i,j) = 1.0;
              for (auto m=0;m<d1;++m) {
                in_c_p_d(i,j,m) = Teuchos::ScalarTraits<value_type>::random();
                for (auto n=0;n<d2;++n)
                  in_c_p_d_d(i,j,m,n) = Teuchos::ScalarTraits<value_type>::random();
              }
            }
          data_c_1(i,0)  = Teuchos::ScalarTraits<value_type>::random();
          datainv_c_1(i,0) = 1./data_c_1(i,0);
          data_c_1_one(i,0) = 1.0;
        }

        art::scalarMultiplyDataData(out_c_p, data_c_p, in_p, true);
        art::scalarMultiplyDataData(outi_c_p, datainv_c_p, out_c_p, true);
        art::scalarMultiplyDataData(in_c_p, data_c_p_one, in_p);
        rst::subtract(outi_c_p, in_c_p);
        if (rst::Serial::vectorNorm(outi_c_p, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT scalarMultiplyDataData (1): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        art::scalarMultiplyDataData(out_c_p_d, data_c_p, in_p_d, true);
        art::scalarMultiplyDataData(outi_c_p_d, datainv_c_p, out_c_p_d, true);
        art::scalarMultiplyDataData(in_c_p_d, data_c_p_one, in_p_d);
        rst::subtract(outi_c_p_d, in_c_p_d);
        if (rst::Serial::vectorNorm(outi_c_p_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT scalarMultiplyDataData (2): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        art::scalarMultiplyDataData(out_c_p_d_d, data_c_p, in_p_d_d, true);
        art::scalarMultiplyDataData(outi_c_p_d_d, datainv_c_p, out_c_p_d_d, true);
        art::scalarMultiplyDataData(in_c_p_d_d, data_c_p_one, in_p_d_d);
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT scalarMultiplyDataData (3): check scalar inverse property\n\n";
          errorFlag = -1000;
        }
        art::scalarMultiplyDataData(out_c_p_d_d, data_c_1, in_p_d_d, true);
        art::scalarMultiplyDataData(outi_c_p_d_d, datainv_c_1, out_c_p_d_d, true);
        art::scalarMultiplyDataData(in_c_p_d_d, data_c_p_one, in_p_d_d);
        rst::subtract(outi_c_p_d_d, in_c_p_d_d);
        if (rst::Serial::vectorNorm(outi_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT scalarMultiplyDataData (4): check scalar inverse property\n\n";
          errorFlag = -1000;
        }

        // fill with constants
        deep_copy(in_p_d_d, 1.0); deep_copy(data_c_p,5.0); deep_copy(data_c_1, 5.0);

        art::scalarMultiplyDataData(out_c_p_d_d, data_c_p, in_p_d_d, true);
        if (std::abs(rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) - (1.0/data_c_p(0,0))*in_p_d_d(0,0,0)*c*p*d1*d2)/rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT scalarMultiplyDataData (5): check result: "
                     << rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                     << (1.0/data_c_p(0,0))*in_p_d_d(0,0,0)*c*p*d1*d2 << "\n\n";
          errorFlag = -1000;
        }
        art::scalarMultiplyDataData(out_c_p_d_d, data_c_1, in_p_d_d, true);
        if (std::abs(rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) - (1.0/data_c_1(0,0))*in_p_d_d(0,0,0)*c*p*d1*d2)/rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) > tol) {
          *outStream << "\n\nINCORRECT scalarMultiplyDataData (6): check result: "
                     << rst::Serial::vectorNorm(out_c_p_d_d, NORM_ONE) << " != "
                     << (1.0/data_c_1(0,0))*in_p_d_d(0,0,0)*c*p*d1*d2 << "\n\n";
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
