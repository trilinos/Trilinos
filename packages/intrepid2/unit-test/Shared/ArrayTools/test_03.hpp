// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_03.hpp
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
      catch (std::logic_error &err) {                                   \
        *outStream << "Expected Error ----------------------------------------------------------------\n"; \
        *outStream << err.what() << '\n';                               \
        *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
      };                                                                \
    }
    
    template<typename ValueType, typename DeviceType>
    int ArrayTools_Test03(const bool verbose) {
      
      typedef ValueType value_type;

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose) 
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      using DeviceExecSpaceType = typename DeviceType::execution_space;
      using HostExecSpaceType = Kokkos::DefaultHostExecutionSpace;

      *outStream << "DeviceSpace::  "; DeviceExecSpaceType().print_configuration(std::cout, false);
      *outStream << "HostSpace::    ";   HostExecSpaceType().print_configuration(std::cout, false);

      *outStream                                                        \
        << "===============================================================================\n" \
        << "|                                                                             |\n" \
        << "|                       Unit Test (ArrayTools)                                |\n" \
        << "|                                                                             |\n" \
        << "|     1) Array operations: dot multiply                                       |\n" \
        << "|                                                                             |\n" \
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
        << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
        << "|                                                                             |\n" \
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
        << "|                                                                             |\n" \
        << "===============================================================================\n";      
      
      typedef RealSpaceTools<DeviceType> rst;
      typedef ArrayTools<DeviceType> art; 
      typedef Kokkos::DynRankView<value_type,DeviceType> DynRankView;

#if defined(INTREPID2_ENABLE_SACADO_ETI_TEST)
#define NumDerivative 10
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__, NumDerivative+1)
#else
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)
#endif 

      const value_type tol = tolerence()*10000.0;
      int errorFlag = 0;

      *outStream                                \
        << "\n"
        << "===============================================================================\n" \
        << "| TEST 1: exceptions                                                   |\n" \
        << "===============================================================================\n";

      try{

#ifdef HAVE_INTREPID2_DEBUG
        DynRankView ConstructWithLabel(a_2, 2);
        DynRankView ConstructWithLabel(a_2_2, 2, 2);
        DynRankView ConstructWithLabel(a_3_2, 3, 2);
        DynRankView ConstructWithLabel(a_2_2_2, 2, 2, 2);
        DynRankView ConstructWithLabel(a_2_2_2_2, 2, 2, 2, 2);
        DynRankView ConstructWithLabel(a_10_1, 10, 1);
        DynRankView ConstructWithLabel(a_10_2, 10, 2);
        DynRankView ConstructWithLabel(a_10_3, 10, 3);
        DynRankView ConstructWithLabel(a_10_1_2, 10, 1, 2);
        DynRankView ConstructWithLabel(a_10_2_2, 10, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_2_2, 10, 2, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_2_2_2, 10, 2, 2, 2, 2);
        DynRankView ConstructWithLabel(a_9_2_2, 9, 2, 2);
        DynRankView ConstructWithLabel(a_10_3_2, 10, 3, 2);
        DynRankView ConstructWithLabel(a_10_2_3, 10, 2, 3);
        DynRankView ConstructWithLabel(a_10_2_2_3, 10, 2, 2, 3);
        DynRankView ConstructWithLabel(a_10_1_2_2, 10, 1, 2, 2);

        *outStream << "-> dotMultiplyDataField:\n";
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_2_2, a_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_2_2, a_10_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2, a_10_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2_2, a_9_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2_2, a_10_3_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2_2, a_10_2_3, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2_2, a_10_2_2_3, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_2_2_2, a_10_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_3_2, a_10_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2_3, a_10_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2_2, a_10_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2_2, a_10_1_2_2, a_10_2_2_2_2) );
        //
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_2_2, a_2, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_2_2, a_10_2_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2_2, a_10_2_2, a_10_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2_2, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2_3, a_10_2_2, a_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2_2, a_9_2_2, a_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2_2, a_10_2_3, a_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2_2, a_10_2_2_3, a_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2_2, a_10_2_2_2, a_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataField(a_10_2_2, a_10_1, a_2_2) );

        *outStream << "-> dotMultiplyDataData:\n";
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_2_2, a_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_2_2, a_10_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2_2, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_9_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_3_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_2_3, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_2_2_2, a_10_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_9_2_2, a_9_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_3_2, a_10_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_2_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_1, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_1_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_1_2_2, a_10_2_2_2) );
        //
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_2_2, a_10_2_2_2_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_2_2, a_10_2_2, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2_2, a_10_2_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_2_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_3, a_10_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_9_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_2_3, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_2_2_3, a_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_2, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_2_2_2, a_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_1, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_1_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::dotMultiplyDataData(a_10_2, a_10_1_2_2, a_2_2_2) );
#endif

      }
      catch (std::logic_error &err) {
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
          *outStream << "\n************ Checking dotMultiplyDataField, (branch 1) ************\n";

          int c=5, p=9, f=7, d1=6, d2=14;

          DynRankView ConstructWithLabel(in_c_f_p, c, f, p);
          DynRankView ConstructWithLabel(in_c_f_p_d, c, f, p, d1);
          DynRankView ConstructWithLabel(in_c_f_p_d_d, c, f, p, d1, d2);
          DynRankView ConstructWithLabel(data_c_p, c, p);
          DynRankView ConstructWithLabel(data_c_p_d, c, p, d1);
          DynRankView ConstructWithLabel(data_c_p_d_d, c, p, d1, d2);
          DynRankView ConstructWithLabel(data_c_1, c, 1);
          DynRankView ConstructWithLabel(data_c_1_d, c, 1, d1);
          DynRankView ConstructWithLabel(data_c_1_d_d, c, 1, d1, d2);
          DynRankView ConstructWithLabel(out_c_f_p, c, f, p);
          DynRankView ConstructWithLabel(outSM_c_f_p, c, f, p);
          DynRankView ConstructWithLabel(outDM_c_f_p, c, f, p);

          const auto in_c_f_p_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_c_f_p);
          const auto in_c_f_p_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_c_f_p_d);

          const auto data_c_p_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), data_c_p);
          const auto data_c_1_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), data_c_1);

          // fill with random numbers
          for (auto i=0;i<c;++i) {
            for (auto j=0;j<f;++j)
              for (auto k=0;k<p;++k){
                in_c_f_p_host(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
                for (auto m=0;m<d1;++m)  //fill with alternating 1's and -1's
                  in_c_f_p_d_host(i,j,k,m) = 1.0 -2.0*((i*f*p*d1+j*p*d1+k*d1+m)%2);
              }
            for (auto j=0;j<p;++j) {
              data_c_p_host(i,j)  = Teuchos::ScalarTraits<value_type>::random();
            }
            data_c_1_host(i,0)  = Teuchos::ScalarTraits<value_type>::random();
          }

          Kokkos::deep_copy(in_c_f_p, in_c_f_p_host);
          Kokkos::deep_copy(in_c_f_p_d, in_c_f_p_d_host);
          Kokkos::deep_copy(data_c_p, data_c_p_host);
          Kokkos::deep_copy(data_c_1, data_c_1_host);

          // fill with 1's
          Kokkos::deep_copy(in_c_f_p_d_d, 1.0);
          Kokkos::deep_copy(data_c_p_d, 1.0); 
          Kokkos::deep_copy(data_c_1_d, 1.0);
          Kokkos::deep_copy(data_c_p_d_d, 1.0); 
          Kokkos::deep_copy(data_c_1_d_d, 1.0);

          art::scalarMultiplyDataField(outSM_c_f_p, data_c_p, in_c_f_p);
          art::dotMultiplyDataField(outDM_c_f_p, data_c_p, in_c_f_p);
          rst::subtract(out_c_f_p, outSM_c_f_p, outDM_c_f_p);
          auto out_c_f_p_h = Kokkos::create_mirror_view(out_c_f_p);
          Kokkos::deep_copy(out_c_f_p_h, out_c_f_p);
          if (rst::Serial::vectorNorm(out_c_f_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataField (1): check dot multiply for scalars vs. scalar multiply\n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataField(out_c_f_p, data_c_p_d, in_c_f_p_d);
          Kokkos::deep_copy(out_c_f_p_h, out_c_f_p);
          if (rst::Serial::vectorNorm(out_c_f_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataField (2): check dot multiply of orthogonal vectors\n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataField(out_c_f_p, data_c_p_d_d, in_c_f_p_d_d);
          Kokkos::deep_copy(out_c_f_p_h, out_c_f_p);
          if ((rst::Serial::vectorNorm(out_c_f_p_h, NORM_INF) - d1*d2) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataField (3): check dot multiply for tensors of 1s\n\n";
            errorFlag = -1000;
          }
          art::scalarMultiplyDataField(outSM_c_f_p, data_c_1, in_c_f_p);
          art::dotMultiplyDataField(outDM_c_f_p, data_c_1, in_c_f_p);
          rst::subtract(out_c_f_p, outSM_c_f_p, outDM_c_f_p);
          Kokkos::deep_copy(out_c_f_p_h, out_c_f_p);
          if (rst::Serial::vectorNorm(out_c_f_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataField (4): check dot multiply for scalars vs. scalar multiply\n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataField(out_c_f_p, data_c_1_d, in_c_f_p_d);
          Kokkos::deep_copy(out_c_f_p_h, out_c_f_p);
          if (rst::Serial::vectorNorm(out_c_f_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataField (5): check dot multiply of orthogonal vectors\n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataField(out_c_f_p, data_c_1_d_d, in_c_f_p_d_d);
          Kokkos::deep_copy(out_c_f_p_h, out_c_f_p);
          if ((rst::Serial::vectorNorm(out_c_f_p_h, NORM_INF) - d1*d2) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataField (6): check dot multiply for tensors of 1s\n\n";
            errorFlag = -1000;
          }
        } // end scope


        { // start scope
          *outStream << "\n************ Checking dotMultiplyDataField, (branch 2) ************\n";

          int c=5, p=9, f=7, d1=6, d2=14;

          DynRankView ConstructWithLabel(in_f_p, f, p);
          DynRankView ConstructWithLabel(in_f_p_d, f, p, d1);
          DynRankView ConstructWithLabel(in_f_p_d_d, f, p, d1, d2);
          DynRankView ConstructWithLabel(data_c_p, c, p);
          DynRankView ConstructWithLabel(data_c_p_d, c, p, d1);
          DynRankView ConstructWithLabel(data_c_p_d_d, c, p, d1, d2);
          DynRankView ConstructWithLabel(data_c_1, c, 1);
          DynRankView ConstructWithLabel(data_c_1_d, c, 1, d1);
          DynRankView ConstructWithLabel(data_c_1_d_d, c, 1, d1, d2);
          DynRankView ConstructWithLabel(out_c_f_p, c, f, p);
          DynRankView ConstructWithLabel(outSM_c_f_p, c, f, p);
          DynRankView ConstructWithLabel(outDM_c_f_p, c, f, p);

          const auto in_f_p_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_f_p);
          const auto in_f_p_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_f_p_d);

          const auto data_c_p_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), data_c_p);
          const auto data_c_1_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), data_c_1);

          // fill with random numbers
          for (auto i=0;i<f;++i)
            for (auto j=0;j<p;++j) {
              in_f_p_host(i,j) = Teuchos::ScalarTraits<value_type>::random();
              for (auto m=0;m<d1;++m)  //fill with alternating 1's and -1's
                in_f_p_d_host(i,j,m) = 1.0 -2.0*((i*p*d1+j*d1+m)%2);
            }

          for (auto i=0;i<c;++i) {
            for (auto j=0;j<p;++j) {
              data_c_p_host(i,j)  = Teuchos::ScalarTraits<value_type>::random();
            }
            data_c_1_host(i,0)  = Teuchos::ScalarTraits<value_type>::random();
          }

          Kokkos::deep_copy(in_f_p, in_f_p_host);
          Kokkos::deep_copy(in_f_p_d, in_f_p_d_host);
          Kokkos::deep_copy(data_c_p, data_c_p_host);
          Kokkos::deep_copy(data_c_1, data_c_1_host);

          // fill with 1's
          Kokkos::deep_copy(in_f_p_d_d, 1.0);
          Kokkos::deep_copy(data_c_p_d, 1.0); 
          Kokkos::deep_copy(data_c_1_d, 1.0);
          Kokkos::deep_copy(data_c_p_d_d, 1.0); 
          Kokkos::deep_copy(data_c_1_d_d, 1.0);

          art::scalarMultiplyDataField(outSM_c_f_p, data_c_p, in_f_p);
          art::dotMultiplyDataField(outDM_c_f_p, data_c_p, in_f_p);
          rst::subtract(out_c_f_p, outSM_c_f_p, outDM_c_f_p);
          auto out_c_f_p_h = Kokkos::create_mirror_view(out_c_f_p);
          Kokkos::deep_copy(out_c_f_p_h, out_c_f_p);
          if (rst::Serial::vectorNorm(out_c_f_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataField (7): check dot multiply for scalars vs. scalar multiply\n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataField(out_c_f_p, data_c_p_d, in_f_p_d);
          Kokkos::deep_copy(out_c_f_p_h, out_c_f_p);
          if (rst::Serial::vectorNorm(out_c_f_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataField (8): check dot multiply of orthogonal vectors\n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataField(out_c_f_p, data_c_p_d_d, in_f_p_d_d);
          Kokkos::deep_copy(out_c_f_p_h, out_c_f_p);
          if ((rst::Serial::vectorNorm(out_c_f_p_h, NORM_INF) - d1*d2) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataField (9): check dot multiply for tensors of 1s\n\n";
            errorFlag = -1000;
          }
          art::scalarMultiplyDataField(outSM_c_f_p, data_c_1, in_f_p);
          art::dotMultiplyDataField(outDM_c_f_p, data_c_1, in_f_p);
          rst::subtract(out_c_f_p, outSM_c_f_p, outDM_c_f_p);
          Kokkos::deep_copy(out_c_f_p_h, out_c_f_p);
          if (rst::Serial::vectorNorm(out_c_f_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataField (10): check dot multiply for scalars vs. scalar multiply\n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataField(out_c_f_p, data_c_1_d, in_f_p_d);
          Kokkos::deep_copy(out_c_f_p_h, out_c_f_p);
          if (rst::Serial::vectorNorm(out_c_f_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataField (11): check dot multiply of orthogonal vectors\n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataField(out_c_f_p, data_c_1_d_d, in_f_p_d_d);
          Kokkos::deep_copy(out_c_f_p_h, out_c_f_p);
          if ((rst::Serial::vectorNorm(out_c_f_p_h, NORM_INF) - d1*d2) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataField (12): check dot multiply for tensors of 1s\n\n";
            errorFlag = -1000;
          }
        } // end scope




        { // start scope
          *outStream << "\n************ Checking dotMultiplyDataData, (branch 1) ************\n";

          int c=5, p=9, d1=6, d2=14;

          DynRankView ConstructWithLabel(in_c_p, c, p);
          DynRankView ConstructWithLabel(in_c_p_d, c, p, d1);
          DynRankView ConstructWithLabel(in_c_p_d_d, c, p, d1, d2);
          DynRankView ConstructWithLabel(data_c_p, c, p);
          DynRankView ConstructWithLabel(data_c_p_d, c, p, d1);
          DynRankView ConstructWithLabel(data_c_p_d_d, c, p, d1, d2);
          DynRankView ConstructWithLabel(data_c_1, c, 1);
          DynRankView ConstructWithLabel(data_c_1_d, c, 1, d1);
          DynRankView ConstructWithLabel(data_c_1_d_d, c, 1, d1, d2);
          DynRankView ConstructWithLabel(out_c_p, c, p);
          DynRankView ConstructWithLabel(outSM_c_p, c, p);
          DynRankView ConstructWithLabel(outDM_c_p, c, p);

          const auto in_c_p_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_c_p);
          const auto in_c_p_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_c_p_d);

          const auto data_c_p_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), data_c_p);
          const auto data_c_1_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), data_c_1);

          for (auto i=0;i<c;++i) {
            for (auto j=0;j<p;++j) {
              in_c_p_host(i,j) = Teuchos::ScalarTraits<value_type>::random();
              data_c_p_host(i,j)  = Teuchos::ScalarTraits<value_type>::random();
              for (auto m=0;m<d1;++m)  //fill with alternating 1's and -1's
                in_c_p_d_host(i,j,m) = 1.0 -2.0*(m%2);
            }
            data_c_1_host(i,0)  = Teuchos::ScalarTraits<value_type>::random();
          }

          Kokkos::deep_copy(in_c_p, in_c_p_host);
          Kokkos::deep_copy(in_c_p_d, in_c_p_d_host);
          Kokkos::deep_copy(data_c_p, data_c_p_host);
          Kokkos::deep_copy(data_c_1, data_c_1_host);

          // fill with 1's
          Kokkos::deep_copy(in_c_p_d_d, 1.0);
          Kokkos::deep_copy(data_c_p_d, 1.0); 
          Kokkos::deep_copy(data_c_1_d, 1.0);
          Kokkos::deep_copy(data_c_p_d_d, 1.0); 
          Kokkos::deep_copy(data_c_1_d_d, 1.0);

          art::scalarMultiplyDataData(outSM_c_p, data_c_p, in_c_p);
          art::dotMultiplyDataData(outDM_c_p, data_c_p, in_c_p);
          rst::subtract(out_c_p, outSM_c_p, outDM_c_p);
          auto out_c_p_h = Kokkos::create_mirror_view(out_c_p);
          Kokkos::deep_copy(out_c_p_h, out_c_p);
          if (rst::Serial::vectorNorm(out_c_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataData (1): check dot multiply for scalars vs. scalar multiply\n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataData(out_c_p, data_c_p_d, in_c_p_d);
          Kokkos::deep_copy(out_c_p_h, out_c_p);
          if (rst::Serial::vectorNorm(out_c_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataData (2): check dot multiply of orthogonal vectors"<<rst::Serial::vectorNorm(out_c_p_h, NORM_ONE)<<" \n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataData(out_c_p, data_c_p_d_d, in_c_p_d_d);
          Kokkos::deep_copy(out_c_p_h, out_c_p);
          if ((rst::Serial::vectorNorm(out_c_p_h, NORM_INF) - d1*d2) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataData (3): check dot multiply for tensors of 1s\n\n";
            errorFlag = -1000;
          }
          art::scalarMultiplyDataData(outSM_c_p, data_c_1, in_c_p);
          art::dotMultiplyDataData(outDM_c_p, data_c_1, in_c_p);
          rst::subtract(out_c_p, outSM_c_p, outDM_c_p);
          Kokkos::deep_copy(out_c_p_h, out_c_p);
          if (rst::Serial::vectorNorm(out_c_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataData (4): check dot multiply for scalars vs. scalar multiply\n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataData(out_c_p, data_c_1_d, in_c_p_d);
          Kokkos::deep_copy(out_c_p_h, out_c_p);
          if (rst::Serial::vectorNorm(out_c_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataData (5): check dot multiply of orthogonal vectors\n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataData(out_c_p, data_c_1_d_d, in_c_p_d_d);
          Kokkos::deep_copy(out_c_p_h, out_c_p);
          if ((rst::Serial::vectorNorm(out_c_p_h, NORM_INF) - d1*d2) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataData (6): check dot multiply for tensors of 1s\n\n";
            errorFlag = -1000;
          }
        } // end scope


        { // start scope
          *outStream << "\n************ Checking dotMultiplyDataData, (branch 2) ************\n";

          int c=5, p=9, d1=6, d2=14;

          DynRankView ConstructWithLabel(in_p, p);
          DynRankView ConstructWithLabel(in_p_d, p, d1);
          DynRankView ConstructWithLabel(in_p_d_d, p, d1, d2);
          DynRankView ConstructWithLabel(data_c_p, c, p);
          DynRankView ConstructWithLabel(data_c_p_d, c, p, d1);
          DynRankView ConstructWithLabel(data_c_p_d_d, c, p, d1, d2);
          DynRankView ConstructWithLabel(data_c_1, c, 1);
          DynRankView ConstructWithLabel(data_c_1_d, c, 1, d1);
          DynRankView ConstructWithLabel(data_c_1_d_d, c, 1, d1, d2);
          DynRankView ConstructWithLabel(out_c_p, c, p);
          DynRankView ConstructWithLabel(outSM_c_p, c, p);
          DynRankView ConstructWithLabel(outDM_c_p, c, p);

          const auto in_p_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_p);
          const auto in_p_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_p_d);

          const auto data_c_p_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), data_c_p);
          const auto data_c_1_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), data_c_1);

          // fill with random numbers
          for (auto i=0;i<c;++i) {
            for (auto j=0;j<p;++j)
              data_c_p_host(i,j)  = 1;//Teuchos::ScalarTraits<value_type>::random();
            data_c_1_host(i,0)  = Teuchos::ScalarTraits<value_type>::random();
          }

          for (auto i=0;i<p;++i) {
            in_p_host(i) = 3; //Teuchos::ScalarTraits<value_type>::random();
            for (auto m=0;m<d1;++m)  //fill with alternating 1's and -1's
              in_p_d_host(i,m) = 1.0 -2.0*(m%2);
          }

          Kokkos::deep_copy(in_p, in_p_host);
          Kokkos::deep_copy(in_p_d, in_p_d_host);
          Kokkos::deep_copy(data_c_p, data_c_p_host);
          Kokkos::deep_copy(data_c_1, data_c_1_host);


          // fill with 1's
          Kokkos::deep_copy(in_p_d_d, 1.0);
          Kokkos::deep_copy(data_c_p_d, 1.0); 
          Kokkos::deep_copy(data_c_1_d, 1.0);
          Kokkos::deep_copy(data_c_p_d_d, 1.0); 
          Kokkos::deep_copy(data_c_1_d_d, 1.0);

          art::scalarMultiplyDataData(outSM_c_p, data_c_p, in_p);
          art::dotMultiplyDataData(outDM_c_p, data_c_p, in_p);
          
          rst::subtract(out_c_p, outSM_c_p, outDM_c_p);
          auto out_c_p_h = Kokkos::create_mirror_view(out_c_p);
          Kokkos::deep_copy(out_c_p_h, out_c_p);
          if (rst::Serial::vectorNorm(out_c_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataData (7): check dot multiply for scalars vs. scalar multiply\n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataData(out_c_p, data_c_p_d, in_p_d);
          Kokkos::deep_copy(out_c_p_h, out_c_p);
          if (rst::Serial::vectorNorm(out_c_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataData (8): check dot multiply of orthogonal vectors\n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataData(out_c_p, data_c_p_d_d, in_p_d_d);
          Kokkos::deep_copy(out_c_p_h, out_c_p);
          if ((rst::Serial::vectorNorm(out_c_p_h, NORM_INF) - d1*d2) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataData (9): check dot multiply for tensors of 1s\n\n";
            errorFlag = -1000;
          }
          art::scalarMultiplyDataData(outSM_c_p, data_c_1, in_p);
          art::dotMultiplyDataData(outDM_c_p, data_c_1, in_p);
          rst::subtract(out_c_p, outSM_c_p, outDM_c_p);
          Kokkos::deep_copy(out_c_p_h, out_c_p);
          if (rst::Serial::vectorNorm(out_c_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataData (10): check dot multiply for scalars vs. scalar multiply\n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataData(out_c_p, data_c_1_d, in_p_d);
          Kokkos::deep_copy(out_c_p_h, out_c_p);
          if (rst::Serial::vectorNorm(out_c_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataData (11): check dot multiply of orthogonal vectors\n\n";
            errorFlag = -1000;
          }
          art::dotMultiplyDataData(out_c_p, data_c_1_d_d, in_p_d_d);
          Kokkos::deep_copy(out_c_p_h, out_c_p);
          if ((rst::Serial::vectorNorm(out_c_p_h, NORM_INF) - d1*d2) > tol) {
            *outStream << "\n\nINCORRECT dotMultiplyDataData (12): check dot multiply for tensors of 1s\n\n";
            errorFlag = -1000;
          }
        } // end scope

          /******************************************/
        *outStream << "\n";
      }
      catch (std::logic_error &err) {
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
