// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_05.cpp
    \brief  Unit test for the clone / scale operations of the ArrayTools class.
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
    int ArrayTools_Test05(const bool verbose) {
      
      typedef ValueType value_type;

      Teuchos::RCP<std::ostream> outStream;
      Teuchos::oblackholestream bhs; // outputs nothing

      if (verbose) 
        outStream = Teuchos::rcp(&std::cout, false);
      else
        outStream = Teuchos::rcp(&bhs,       false);

      Teuchos::oblackholestream oldFormatState;
      oldFormatState.copyfmt(std::cout);

      using DeviceSpaceType = typename DeviceType::execution_space;
      using HostSpaceType = Kokkos::DefaultHostExecutionSpace;
      
      *outStream << "DeviceSpace::  "; DeviceSpaceType().print_configuration(std::cout, false);
      *outStream << "HostSpace::    ";   HostSpaceType().print_configuration(std::cout, false);

      *outStream      \
        << "===============================================================================\n" \
        << "|                                                                             |\n" \
        << "|                       Unit Test (ArrayTools)                                |\n" \
        << "|                                                                             |\n" \
        << "|     1) Array operations: clone / scale                                       |\n" \
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

      const value_type tol = tolerence()*100.0;
      int errorFlag = 0;

      *outStream                                \
        << "\n"
        << "===============================================================================\n" \
        << "| TEST 1: exceptions                                                   |\n" \
        << "===============================================================================\n";




      try{

#ifdef HAVE_INTREPID2_DEBUG
        DynRankView ConstructWithLabel(a_2, 2);
        DynRankView ConstructWithLabel(a_9_2, 9, 2);
        DynRankView ConstructWithLabel(a_10_2, 10, 2);
        DynRankView ConstructWithLabel(a_10_3, 10, 3);
        DynRankView ConstructWithLabel(a_10_2_2, 10, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_2_2, 10, 2, 2, 2);
        DynRankView ConstructWithLabel(a_10_3_2_2, 10, 3, 2, 2);
        DynRankView ConstructWithLabel(a_10_3_2, 10, 3, 2);
        DynRankView ConstructWithLabel(a_2_2, 2, 2);
        DynRankView ConstructWithLabel(a_2_3_2_2, 2, 3, 2, 2);
        DynRankView ConstructWithLabel(a_2_2_2_2, 2, 2, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_2_2_2, 10, 2, 2, 2, 2);
        DynRankView ConstructWithLabel(a_10_3_2_2_2, 10, 3, 2, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_3_2_2, 10, 2, 3, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_2_3_2, 10, 2, 2, 3, 2);
        DynRankView ConstructWithLabel(a_10_2_2_2_3, 10, 2, 2, 2, 3);

        *outStream << "-> cloneFields:\n";
        INTREPID2_TEST_ERROR_EXPECTED( art::cloneFields(a_10_2_2_2, a_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::cloneFields(a_10_2_2_2, a_10_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::cloneFields(a_10_3_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::cloneFields(a_10_2_2_2_2, a_2_3_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::cloneFields(a_10_2_2_3_2, a_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::cloneFields(a_10_2_2_2_3, a_2_2_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::cloneFields(a_10_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( art::cloneFields(a_10_2_2_2_2, a_2_2_2_2) );
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
          *outStream << "\n************ Checking cloneFields ************\n";

          int c=5, p=9, f=7, d1=7, d2=13;

          DynRankView ConstructWithLabel(in_f_p, f, p);
          DynRankView ConstructWithLabel(in_f_p_d, f, p, d1);
          DynRankView ConstructWithLabel(in_f_p_d_d, f, p, d1, d2);
          DynRankView ConstructWithLabel(in_c_f_p, c, f, p);
          DynRankView ConstructWithLabel(in_c_f_p_d, c, f, p, d1);
          DynRankView ConstructWithLabel(in_c_f_p_d_d, c, f, p, d1, d2);
          DynRankView ConstructWithLabel(data_c_p_one, c, p);
          DynRankView ConstructWithLabel(out_c_f_p, c, f, p);
          DynRankView ConstructWithLabel(out_c_f_p_d, c, f, p, d1);
          DynRankView ConstructWithLabel(out_c_f_p_d_d, c, f, p, d1, d2);

          const auto in_f_p_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_f_p);
          const auto in_f_p_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_f_p_d);
          const auto in_f_p_d_d_host = Kokkos::create_mirror_view(Kokkos::HostSpace(), in_f_p_d_d);

          // fill with random numbers
          for (int i=0; i < f; i++)
            for (int j=0; j < p; j++) {
              in_f_p_host(i, j) = Teuchos::ScalarTraits<value_type>::random();
              for (int m=0; m<d1; m++) {
                in_f_p_d_host(i, j, m) = Teuchos::ScalarTraits<value_type>::random();
                for (int n=0; n<d2; n++)
                  in_f_p_d_d_host(i, j, m, n)  = Teuchos::ScalarTraits<value_type>::random();
              }
            }

          Kokkos::deep_copy(in_f_p, in_f_p_host);
          Kokkos::deep_copy(in_f_p_d, in_f_p_d_host);
          Kokkos::deep_copy(in_f_p_d_d, in_f_p_d_d_host);

          Kokkos::deep_copy(data_c_p_one, 1.0);

          art::cloneFields(out_c_f_p, in_f_p);
          art::scalarMultiplyDataField(in_c_f_p, data_c_p_one, in_f_p);
          rst::subtract(out_c_f_p, in_c_f_p);
          auto out_c_f_p_h = Kokkos::create_mirror_view(out_c_f_p);
          Kokkos::deep_copy(out_c_f_p_h, out_c_f_p);
          if (rst::Serial::vectorNorm(out_c_f_p_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT cloneFields (1): check multiplyScalarData vs. cloneFields\n\n";
            errorFlag = -1000;
          }
          art::cloneFields(out_c_f_p_d, in_f_p_d);
          art::scalarMultiplyDataField(in_c_f_p_d, data_c_p_one, in_f_p_d);
          rst::subtract(out_c_f_p_d, in_c_f_p_d);
          auto out_c_f_p_d_h = Kokkos::create_mirror_view(out_c_f_p_d);
          Kokkos::deep_copy(out_c_f_p_d_h, out_c_f_p_d);
          if (rst::Serial::vectorNorm(out_c_f_p_d_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT cloneFields (2): check multiplyScalarData vs. cloneFields\n\n";
            errorFlag = -1000;
          }
          art::cloneFields(out_c_f_p_d_d, in_f_p_d_d);
          art::scalarMultiplyDataField(in_c_f_p_d_d, data_c_p_one, in_f_p_d_d);
          rst::subtract(out_c_f_p_d_d, in_c_f_p_d_d);
          auto out_c_f_p_d_d_h = Kokkos::create_mirror_view(out_c_f_p_d_d);
          Kokkos::deep_copy(out_c_f_p_d_d_h, out_c_f_p_d_d);
          if (rst::Serial::vectorNorm(out_c_f_p_d_d_h, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT cloneFields (3): check multiplyScalarData vs. cloneFields\n\n";
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
