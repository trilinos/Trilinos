// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
      catch (std::logic_error &err) {                                    \
        *outStream << "Expected Error ----------------------------------------------------------------\n"; \
        *outStream << err.what() << '\n';                               \
        *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
      };                                                                \
    }
    
    template<typename ValueType, typename DeviceType>
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

      using DeviceExecSpaceType = typename DeviceType::execution_space;
      using HostExecSpaceType = Kokkos::DefaultHostExecutionSpace;

      *outStream << "DeviceSpace::  "; DeviceExecSpaceType().print_configuration(std::cout, false);
      *outStream << "HostSpace::    ";   HostExecSpaceType().print_configuration(std::cout, false);

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
        DynRankView ConstructWithLabel(a_2_2, 2, 2);
        DynRankView ConstructWithLabel(a_10_2, 10, 2);
        DynRankView ConstructWithLabel(a_10_3, 10, 3);
        DynRankView ConstructWithLabel(a_10_2_2, 10, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_3, 10, 2, 3);        
        DynRankView ConstructWithLabel(a_10_3_2, 10, 3, 2);
        DynRankView ConstructWithLabel(a_9_2_2, 9, 2, 2);

        *outStream << "-> contractFieldFieldScalar:\n";
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldScalar(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldScalar(a_2_2, a_10_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldScalar(a_2_2, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldScalar(a_10_2_2, a_9_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldScalar(a_10_2_2, a_10_2_2, a_10_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldScalar(a_9_2_2, a_10_3_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldScalar(a_10_3_2, a_10_2_2, a_10_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldScalar(a_10_2_2, a_10_2_2, a_10_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldScalar(a_10_2_3, a_10_2_2, a_10_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldScalar(a_10_2_3, a_10_2_2, a_10_3_2) );

        DynRankView ConstructWithLabel(a_10_2_2_2, 10, 2, 2, 2);
        DynRankView ConstructWithLabel(a_9_2_2_2, 9, 2, 2, 2);
        DynRankView ConstructWithLabel(a_10_3_2_2, 10, 3, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_3_2, 10, 2, 3, 2);
        DynRankView ConstructWithLabel(a_10_2_2_3, 10, 2, 2, 3);

        *outStream << "-> contractFieldFieldVector:\n";
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldVector(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldVector(a_2_2, a_10_2_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldVector(a_2_2, a_10_2_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldVector(a_10_2_2, a_9_2_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldVector(a_10_2_2, a_10_2_2_2, a_10_2_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldVector(a_10_2_2, a_10_2_2_2, a_10_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldVector(a_9_2_2, a_10_2_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldVector(a_10_3_2, a_10_2_2_2, a_10_3_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldVector(a_10_2_2, a_10_2_2_2, a_10_3_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldVector(a_10_2_3, a_10_2_2_2, a_10_3_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldVector(a_10_2_3, a_10_2_2_2, a_10_3_2_2) );

        DynRankView ConstructWithLabel(a_10_2_2_2_2, 10, 2, 2, 2, 2);
        DynRankView ConstructWithLabel(a_9_2_2_2_2, 9, 2, 2, 2, 2);
        DynRankView ConstructWithLabel(a_10_3_2_2_2, 10, 3, 2, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_3_2_2, 10, 2, 3, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_2_3_2, 10, 2, 2, 3, 2);
        DynRankView ConstructWithLabel(a_10_2_2_2_3, 10, 2, 2, 2, 3);

        *outStream << "-> contractFieldFieldTensor:\n";
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldTensor(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldTensor(a_2_2, a_10_2_2_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldTensor(a_2_2, a_10_2_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldTensor(a_10_2_2, a_9_2_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldTensor(a_10_2_2, a_10_2_2_2_2, a_10_2_3_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldTensor(a_10_2_2, a_10_2_2_2_2, a_10_2_2_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldTensor(a_10_2_2, a_10_2_2_2_2, a_10_2_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldTensor(a_9_2_2, a_10_2_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldTensor(a_10_3_2, a_10_2_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldTensor(a_10_2_2, a_10_2_2_2_2, a_10_3_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldTensor(a_10_2_3, a_10_2_2_2_2, a_10_3_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractFieldFieldTensor(a_10_2_3, a_10_2_2_2_2, a_10_3_2_2_2) );

        DynRankView ConstructWithLabel(a_9_2, 9, 2);
        DynRankView ConstructWithLabel(a_10_1, 10, 1);

        *outStream << "-> contractDataFieldScalar:\n";
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldScalar(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldScalar(a_2_2, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldScalar(a_10_2_2, a_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldScalar(a_2_2, a_10_2, a_9_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldScalar(a_2_2, a_10_3, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldScalar(a_9_2, a_10_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldScalar(a_10_2, a_10_2, a_10_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldScalar(a_10_2, a_10_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldScalar(a_10_2, a_10_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldScalar(a_10_2, a_10_1, a_10_2_2) );

        DynRankView ConstructWithLabel(a_10_1_2, 10, 1, 2);
        DynRankView ConstructWithLabel(a_10_1_3, 10, 1, 3);

        *outStream << "-> contractDataFieldVector:\n";
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldVector(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldVector(a_2_2, a_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldVector(a_10_2_2, a_10_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldVector(a_2_2, a_10_2_2, a_9_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldVector(a_2_2, a_10_2_2, a_10_2_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldVector(a_2_2, a_10_2_2, a_10_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldVector(a_9_2, a_10_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldVector(a_10_2, a_10_2_2, a_10_3_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldVector(a_10_2, a_10_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldVector(a_10_2, a_10_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldVector(a_10_2, a_10_1_2, a_10_2_2_2) );

        DynRankView ConstructWithLabel(a_10_1_2_2, 10, 1, 2, 2);

        *outStream << "-> contractDataFieldTensor:\n";
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldTensor(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldTensor(a_2_2, a_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldTensor(a_10_2_2, a_10_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldTensor(a_2_2, a_10_2_2_2, a_9_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldTensor(a_2_2, a_10_2_2_2, a_10_2_3_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldTensor(a_2_2, a_10_2_2_2, a_10_2_2_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldTensor(a_2_2, a_10_2_2_2, a_10_2_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldTensor(a_9_2, a_10_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldTensor(a_10_2, a_10_2_2_2, a_10_3_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldTensor(a_10_2, a_10_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldTensor(a_10_2, a_10_2_2_2, a_10_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataFieldTensor(a_10_2, a_10_1_2_2, a_10_2_2_2_2) );

        DynRankView ConstructWithLabel(a_2, 2);
        DynRankView ConstructWithLabel(a_10, 10);

        *outStream << "-> contractDataDataScalar:\n";
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataScalar(a_2_2, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataScalar(a_2_2, a_10_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataScalar(a_2_2, a_10_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataScalar(a_2, a_9_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataScalar(a_2, a_10_2, a_10_3) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataScalar(a_2, a_10_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataScalar(a_10, a_10_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataScalar(a_10, a_10_2, a_10_2) );

        *outStream << "-> contractDataDataVector:\n";
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataVector(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataVector(a_2_2, a_10_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataVector(a_2_2, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataVector(a_10, a_9_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataVector(a_10, a_10_3_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataVector(a_10, a_10_2_3, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataVector(a_2, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataVector(a_10, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataVector(a_10, a_10_2_2, a_10_2_2) );

        *outStream << "-> contractDataDataTensor:\n";
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataTensor(a_2_2, a_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataTensor(a_2_2, a_10_2_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataTensor(a_2_2, a_10_2_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataTensor(a_10, a_9_2_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataTensor(a_10, a_10_2_2_2, a_10_3_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataTensor(a_10, a_10_2_2_2, a_10_2_3_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataTensor(a_10, a_10_2_2_2, a_10_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataTensor(a_2, a_10_2_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataTensor(a_10, a_10_2_2_2, a_10_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::contractDataDataTensor(a_10, a_10_2_2_2, a_10_2_2_2) );

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
          *outStream << "\n************ Checking contractFieldFieldScalar ************\n";

          int c=5, p=9, l=3, r=7;

          DynRankView ConstructWithLabel(in_c_l_p, c, l, p);
          DynRankView ConstructWithLabel(in_c_r_p, c, r, p);
          DynRankView ConstructWithLabel(out1_c_l_r, c, l, r);
          DynRankView ConstructWithLabel(out2_c_l_r, c, l, r);
          

          const auto in_c_l_p_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_c_l_p);
          const auto in_c_r_p_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_c_r_p);
          const auto out1_c_l_r_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out1_c_l_r);
          const auto out2_c_l_r_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out2_c_l_r);
          
          // fill with random numbers
          //Kokkos::Random_XorShift64_Pool<DeviceExecSpaceType> random(13718); 
          //Kokkos::fill_random(in_c_l_p, random, value_type(1));
          //Kokkos::fill_random(in_c_r_p, random, value_type(1));

          // fill with random numbers
          for (auto i=0;i<c;++i) {
            for (auto j=0;j<l;++j)
              for (auto k=0;k<p;++k)
                in_c_l_p_host(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
            for (auto j=0;j<r;++j)
              for (auto k=0;k<p;++k)
                in_c_r_p_host(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
          }
          
          Kokkos::deep_copy(in_c_l_p, in_c_l_p_host);
          Kokkos::deep_copy(in_c_r_p, in_c_r_p_host);
          
          art::contractFieldFieldScalar(out1_c_l_r, in_c_l_p, in_c_r_p);
          art::contractFieldFieldScalar(out2_c_l_r, in_c_l_p, in_c_r_p);
          rst::subtract(out1_c_l_r, out2_c_l_r);
          
          Kokkos::deep_copy(out1_c_l_r_host, out1_c_l_r);
          if (rst::Serial::vectorNorm(out1_c_l_r_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractFieldFieldScalar (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_r_host, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          // with sumInto:
          
          //fill with 2.0
          Kokkos::deep_copy(out1_c_l_r, value_type(2));
          Kokkos::deep_copy(out2_c_l_r, value_type(2));

          art::contractFieldFieldScalar(out1_c_l_r, in_c_l_p, in_c_r_p, true);
          art::contractFieldFieldScalar(out2_c_l_r, in_c_l_p, in_c_r_p, true);
          rst::subtract(out1_c_l_r, out2_c_l_r);
          Kokkos::deep_copy(out1_c_l_r_host, out1_c_l_r);
          if (rst::Serial::vectorNorm(out1_c_l_r_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractFieldFieldScalar (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_r_host, NORM_ONE) << "\n\n";
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

          const auto in_c_l_p_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_c_l_p_d);
          const auto in_c_r_p_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_c_r_p_d);
          const auto out1_c_l_r_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out1_c_l_r);
          const auto out2_c_l_r_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out2_c_l_r);

          // fill with random numbers
          for (auto i=0;i<c;++i) {
            for (auto j=0;j<l;++j)
              for (auto k=0;k<p;++k)
                for (auto m=0;m<d;++m)
                  in_c_l_p_d_host(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
            for (auto j=0;j<r;++j)
              for (auto k=0;k<p;++k)
                for (auto m=0;m<d;++m)
                  in_c_r_p_d_host(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
          }

          Kokkos::deep_copy(in_c_l_p_d, in_c_l_p_d_host);
          Kokkos::deep_copy(in_c_r_p_d, in_c_r_p_d_host);

          art::contractFieldFieldVector(out1_c_l_r, in_c_l_p_d, in_c_r_p_d);
          art::contractFieldFieldVector(out2_c_l_r, in_c_l_p_d, in_c_r_p_d);

          rst::subtract(out1_c_l_r, out2_c_l_r);
          Kokkos::deep_copy(out1_c_l_r_host, out1_c_l_r);

          if (rst::Serial::vectorNorm(out1_c_l_r_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractFieldFieldVector (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_r_host, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          
          // with sumInto:

          //fill with 2.0
          Kokkos::deep_copy(out1_c_l_r, value_type(2));
          Kokkos::deep_copy(out2_c_l_r, value_type(2));

          art::contractFieldFieldVector(out1_c_l_r, in_c_l_p_d, in_c_r_p_d, true);
          art::contractFieldFieldVector(out2_c_l_r, in_c_l_p_d, in_c_r_p_d, true);

          rst::subtract(out1_c_l_r, out2_c_l_r);
          Kokkos::deep_copy(out1_c_l_r_host, out1_c_l_r);

          if (rst::Serial::vectorNorm(out1_c_l_r_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractFieldFieldVector (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_r_host, NORM_ONE) << "\n\n";
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

          const auto in_c_l_p_d_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_c_l_p_d_d);
          const auto in_c_r_p_d_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_c_r_p_d_d);
          const auto out1_c_l_r_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out1_c_l_r);
          const auto out2_c_l_r_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out2_c_l_r);

          // fill with random numbers
          for (auto i=0;i<c;++i) {
            for (auto j=0;j<l;++j)
              for (auto k=0;k<p;++k)
                for (auto m=0;m<d1;++m)
                  for (auto n=0;n<d2;++n)
                    in_c_l_p_d_d_host(i,j,k,m,n) = Teuchos::ScalarTraits<value_type>::random();
            for (auto j=0;j<r;++j)
              for (auto k=0;k<p;++k)
                for (auto m=0;m<d1;++m)
                  for (auto n=0;n<d2;++n)
                    in_c_r_p_d_d_host(i,j,k,m,n) = Teuchos::ScalarTraits<value_type>::random();
          }

          Kokkos::deep_copy(in_c_l_p_d_d, in_c_l_p_d_d_host);
          Kokkos::deep_copy(in_c_r_p_d_d, in_c_r_p_d_d_host);

          art::contractFieldFieldTensor(out1_c_l_r, in_c_l_p_d_d, in_c_r_p_d_d);
          art::contractFieldFieldTensor(out2_c_l_r, in_c_l_p_d_d, in_c_r_p_d_d);

          rst::subtract(out1_c_l_r, out2_c_l_r);
          Kokkos::deep_copy(out1_c_l_r_host, out1_c_l_r);
          if (rst::Serial::vectorNorm(out1_c_l_r_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractFieldFieldTensor (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_r_host, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }

          // with sumInto:

          //fill with 2.0
          Kokkos::deep_copy(out1_c_l_r, value_type(2));
          Kokkos::deep_copy(out2_c_l_r, value_type(2));

          art::contractFieldFieldTensor(out1_c_l_r, in_c_l_p_d_d, in_c_r_p_d_d, true);
          art::contractFieldFieldTensor(out2_c_l_r, in_c_l_p_d_d, in_c_r_p_d_d, true);

          rst::subtract(out1_c_l_r, out2_c_l_r);
          Kokkos::deep_copy(out1_c_l_r_host, out1_c_l_r);
          if (rst::Serial::vectorNorm(out1_c_l_r_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractFieldFieldTensor (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_r_host, NORM_ONE) << "\n\n";
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

          const auto in_c_l_p_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_c_l_p);
          const auto data_c_p_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), data_c_p);
          const auto data_c_1_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), data_c_1);
          const auto out1_c_l_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out1_c_l);
          const auto out2_c_l_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out2_c_l);

          // fill with random numbers
          for (auto i=0;i<c;++i) {
            for (auto j=0;j<l;++j)
              for (auto k=0;k<p;++k)
                in_c_l_p_host(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
            for (auto j=0;j<p;++j)
              data_c_p_host(i,j) = Teuchos::ScalarTraits<value_type>::random();
            data_c_1_host(i,0) = Teuchos::ScalarTraits<value_type>::random();
          }

          // nonconstant data
          art::contractDataFieldScalar(out1_c_l, data_c_p, in_c_l_p);
          art::contractDataFieldScalar(out2_c_l, data_c_p, in_c_l_p);
          rst::subtract(out1_c_l, out2_c_l);

          Kokkos::deep_copy(out1_c_l_host, out1_c_l);
          if (rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractDataFieldScalar (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          // constant data
          art::contractDataFieldScalar(out1_c_l, data_c_1, in_c_l_p);
          art::contractDataFieldScalar(out2_c_l, data_c_1, in_c_l_p);
          rst::subtract(out1_c_l, out2_c_l);
          Kokkos::deep_copy(out1_c_l_host, out1_c_l);
          if (rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractDataFieldScalar (2): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }

          // nonconstant data with sumInto
          Kokkos::deep_copy(out1_c_l, value_type(2));
          Kokkos::deep_copy(out2_c_l, value_type(2));

          art::contractDataFieldScalar(out1_c_l, data_c_p, in_c_l_p, true);
          art::contractDataFieldScalar(out2_c_l, data_c_p, in_c_l_p, true);
          rst::subtract(out1_c_l, out2_c_l);

          Kokkos::deep_copy(out1_c_l_host, out1_c_l);
          if (rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractDataFieldScalar (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) << "\n\n";
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

          const auto in_c_l_p_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_c_l_p_d);
          const auto data_c_p_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), data_c_p_d);
          const auto data_c_1_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), data_c_1_d);
          const auto out1_c_l_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out1_c_l);
          const auto out2_c_l_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out2_c_l);

          // fill with random numbers
          for (auto i=0;i<c;++i) {
            for (auto j=0;j<l;++j)
              for (auto k=0;k<p;++k)
                for (auto m=0;m<d;++m)
                  in_c_l_p_d_host(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
            for (auto j=0;j<p;++j)
              for (auto k=0;k<d;++k)
                data_c_p_d_host(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
            for (auto k=0;k<d;++k)
              data_c_1_d_host(i,0,k) = Teuchos::ScalarTraits<value_type>::random();
          }

          Kokkos::deep_copy(in_c_l_p_d, in_c_l_p_d_host);
          Kokkos::deep_copy(data_c_p_d, data_c_p_d_host);
          Kokkos::deep_copy(data_c_1_d, data_c_1_d_host);

          // nonconstant data
          art::contractDataFieldVector(out1_c_l, data_c_p_d, in_c_l_p_d);
          art::contractDataFieldVector(out2_c_l, data_c_p_d, in_c_l_p_d);
          rst::subtract(out1_c_l, out2_c_l);

          Kokkos::deep_copy(out1_c_l_host, out1_c_l);
          if (rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractDataFieldVector (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          // constant data
          art::contractDataFieldVector(out1_c_l, data_c_1_d, in_c_l_p_d);
          art::contractDataFieldVector(out2_c_l, data_c_1_d, in_c_l_p_d);
          rst::subtract(out1_c_l, out2_c_l);
          Kokkos::deep_copy(out1_c_l_host, out1_c_l);
          if (rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractDataFieldVector (2): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          // nonconstant data with sumInto

          //fill with 2.0
          Kokkos::deep_copy(out1_c_l, value_type(2));
          Kokkos::deep_copy(out2_c_l, value_type(2));
          
          art::contractDataFieldVector(out1_c_l, data_c_p_d, in_c_l_p_d, true);
          art::contractDataFieldVector(out2_c_l, data_c_p_d, in_c_l_p_d, true);
          rst::subtract(out1_c_l, out2_c_l);
          Kokkos::deep_copy(out1_c_l_host, out1_c_l);
          if (rst::Serial::vectorNorm(out1_c_l_host,NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractDataFieldVector (3): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) << "\n\n";
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

          const auto in_c_l_p_d_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_c_l_p_d_d);          
          const auto data_c_p_d_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), data_c_p_d_d);          
          const auto data_c_1_d_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), data_c_1_d_d);          
          const auto out1_c_l_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out1_c_l);          
          const auto out2_c_l_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out2_c_l);          

          // fill with random numbers
          for (auto i=0;i<c;++i) {
            for (auto j=0;j<l;++j)
              for (auto k=0;k<p;++k)
                for (auto m=0;m<d1;++m)
                  for (auto n=0;n<d2;++n)
                    in_c_l_p_d_d_host(i,j,k,m,n) = Teuchos::ScalarTraits<value_type>::random();
            for (auto j=0;j<p;++j)
              for (auto k=0;k<d1;++k)
                for (auto m=0;m<d2;++m)
                  data_c_p_d_d_host(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
            for (auto k=0;k<d1;++k)
              for (auto m=0;m<d2;++m)
                data_c_1_d_d_host(i,0,k,m) = Teuchos::ScalarTraits<value_type>::random();
          }

          // nonconstant data
          art::contractDataFieldTensor(out1_c_l, data_c_p_d_d, in_c_l_p_d_d);
          art::contractDataFieldTensor(out2_c_l, data_c_p_d_d, in_c_l_p_d_d);
          rst::subtract(out1_c_l, out2_c_l);
          Kokkos::deep_copy(out1_c_l_host, out1_c_l);
          if (rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractDataFieldTensor (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }

          // constant data
          art::contractDataFieldTensor(out1_c_l, data_c_1_d_d, in_c_l_p_d_d);
          art::contractDataFieldTensor(out2_c_l, data_c_1_d_d, in_c_l_p_d_d);
          rst::subtract(out1_c_l, out2_c_l);
          Kokkos::deep_copy(out1_c_l_host, out1_c_l);
          if (rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractDataFieldTensor (2): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }

          // nonconstant data with sumInto

          //fill with 2.0
          Kokkos::deep_copy(out1_c_l, value_type(2));
          Kokkos::deep_copy(out2_c_l, value_type(2));

          art::contractDataFieldTensor(out1_c_l, data_c_p_d_d, in_c_l_p_d_d, true);
          art::contractDataFieldTensor(out2_c_l, data_c_p_d_d, in_c_l_p_d_d, true);
          rst::subtract(out1_c_l, out2_c_l);
          Kokkos::deep_copy(out1_c_l_host, out1_c_l);
          if (rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractDataFieldTensor (3): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_l_host, NORM_ONE) << "\n\n";
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

          const auto inl_c_p_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), inl_c_p);
          const auto inr_c_p_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), inr_c_p);
          const auto out1_c_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out1_c);
          const auto out2_c_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out2_c);

          // fill with random numbers
          for (auto i=0;i<c;++i)
            for (auto j=0;j<p;++j) {
              inl_c_p_host(i,j) = Teuchos::ScalarTraits<value_type>::random();
              inr_c_p_host(i,j) = Teuchos::ScalarTraits<value_type>::random();
            }

          Kokkos::deep_copy(inl_c_p, inl_c_p_host);
          Kokkos::deep_copy(inr_c_p, inr_c_p_host);

          art::contractDataDataScalar(out1_c, inl_c_p, inr_c_p);
          art::contractDataDataScalar(out2_c, inl_c_p, inr_c_p);
          rst::subtract(out1_c, out2_c);

          Kokkos::deep_copy(out1_c_host, out1_c);
          if (rst::Serial::vectorNorm(out1_c_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractDataDataScalar (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_host, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }
          // with sumInto:

          //fill with 2.0
          Kokkos::deep_copy(out1_c, value_type(2));
          Kokkos::deep_copy(out2_c, value_type(2));

          art::contractDataDataScalar(out1_c, inl_c_p, inr_c_p, true);
          art::contractDataDataScalar(out2_c, inl_c_p, inr_c_p, true);
          rst::subtract(out1_c, out2_c);
          Kokkos::deep_copy(out1_c_host, out1_c);
          if (rst::Serial::vectorNorm(out1_c_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractDataDataScalar (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_host, NORM_ONE) << "\n\n";
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

          const auto inl_c_p_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), inl_c_p_d);
          const auto inr_c_p_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), inr_c_p_d);
          const auto out1_c_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out1_c);
          const auto out2_c_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out2_c);

          // fill with random numbers
          for (auto i=0;i<c;++i)
            for (auto j=0;j<p;++j)
              for (auto k=0;k<d;++k) {
              inl_c_p_d_host(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
              inr_c_p_d_host(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
            }

          art::contractDataDataVector(out1_c, inl_c_p_d, inr_c_p_d);
          art::contractDataDataVector(out2_c, inl_c_p_d, inr_c_p_d);

          rst::subtract(out1_c, out2_c);

          Kokkos::deep_copy(out1_c_host, out1_c);
          if (rst::Serial::vectorNorm(out1_c_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractDataDataVector (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_host, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }

          // with sumInto:

          //fill with 2.0
          Kokkos::deep_copy(out1_c, value_type(2));
          Kokkos::deep_copy(out2_c, value_type(2));
          
          art::contractDataDataVector(out1_c, inl_c_p_d, inr_c_p_d, true);
          art::contractDataDataVector(out2_c, inl_c_p_d, inr_c_p_d, true);

          rst::subtract(out1_c, out2_c);
          Kokkos::deep_copy(out1_c_host, out1_c);
          if (rst::Serial::vectorNorm(out1_c_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractDataDataVector (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_host, NORM_ONE) << "\n\n";
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

          const auto inl_c_p_d_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), inl_c_p_d_d);
          const auto inr_c_p_d_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), inr_c_p_d_d);
          const auto out1_c_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out1_c);
          const auto out2_c_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), out2_c);

          // fill with random numbers
          for (auto i=0;i<c;++i)
            for (auto j=0;j<p;++j)
              for (auto k=0;k<d1;++k)
                for (auto m=0;m<d2;++m) {
                  inl_c_p_d_d_host(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
                  inr_c_p_d_d_host(i,j,k,m) = Teuchos::ScalarTraits<value_type>::random();
                }

          art::contractDataDataTensor(out1_c, inl_c_p_d_d, inr_c_p_d_d);
          art::contractDataDataTensor(out2_c, inl_c_p_d_d, inr_c_p_d_d);

          rst::subtract(out1_c, out2_c);

          Kokkos::deep_copy(out1_c_host, out1_c);
          if (rst::Serial::vectorNorm(out1_c_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractDataDataTensor (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_host, NORM_ONE) << "\n\n";
            errorFlag = -1000;
          }

          // with sumInto:

          //fill with 2.0
          Kokkos::deep_copy(out1_c, value_type(2));
          Kokkos::deep_copy(out2_c, value_type(2));

          art::contractDataDataTensor(out1_c, inl_c_p_d_d, inr_c_p_d_d, true);
          art::contractDataDataTensor(out2_c, inl_c_p_d_d, inr_c_p_d_d, true);

          rst::subtract(out1_c, out2_c);
          Kokkos::deep_copy(out1_c_host, out1_c);
          if (rst::Serial::vectorNorm(out1_c_host, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT contractDataDataTensor (1): check COMP_CPP vs. COMP_BLAS; "
                       << " diff-1norm = " << rst::Serial::vectorNorm(out1_c_host, NORM_ONE) << "\n\n";
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
