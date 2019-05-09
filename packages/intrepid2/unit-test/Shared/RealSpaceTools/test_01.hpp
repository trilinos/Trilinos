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
\brief  Unit test for the RealSpaceTools class.
\author Created by P. Bochev and D. Ridzal and Kyungjoo Kim.
*/

#include "Intrepid2_config.h"

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
#endif

#include "Intrepid2_RealSpaceTools.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

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
    int RealSpaceTools_Test01(const bool verbose) {
      
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

      *outStream << "DeviceSpace::  "; DeviceSpaceType::print_configuration(*outStream, false);
      *outStream << "HostSpace::    ";   HostSpaceType::print_configuration(*outStream, false);

      *outStream                                                        \
        << "===============================================================================\n" \
        << "|                                                                             |\n" \
        << "|                       Unit Test (RealSpaceTools)                            |\n" \
        << "|                                                                             |\n" \
        << "|     1) Vector operations in 1D, 2D, and 3D real space                       |\n" \
        << "|     2) Matrix / matrix-vector operations in 1D, 2D, and 3D real space       |\n" \
        << "|                                                                             |\n" \
        << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
        << "|                      Denis Ridzal (dridzal@sandia.gov) or                   |\n" \
        << "|                      Kyungjoo Kim (kyukim@sandia.gov).                      |\n" \
        << "|                                                                             |\n" \
        << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
        << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
        << "|                                                                             |\n" \
        << "===============================================================================\n";

      typedef RealSpaceTools<DeviceSpaceType> rst;
      typedef Kokkos::DynRankView<value_type,DeviceSpaceType> DynRankView;
#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

      int errorFlag = 0;

      *outStream                                \
        << "\n"
        << "===============================================================================\n" \
        << "| TEST 1: vector exceptions                                                   |\n" \
        << "===============================================================================\n";

      try {
        DynRankView ConstructWithLabel(a_2_2, 2, 2);
        DynRankView ConstructWithLabel(a_10_2, 10, 2);
        DynRankView ConstructWithLabel(a_10_3, 10, 3);
        DynRankView ConstructWithLabel(a_10_2_2, 10, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_3, 10, 2, 3);
        DynRankView ConstructWithLabel(a_10_2_2_3, 10, 2, 2, 3);
        
#ifdef HAVE_INTREPID2_DEBUG
        *outStream << "-> Errors are expected due to rank mismatches:\n";

        *outStream << "-> vector norm with multidimensional arrays:\n";

        rst::Serial::vectorNorm(a_2_2, NORM_TWO);
        INTREPID2_TEST_ERROR_EXPECTED( rst::vectorNorm(a_10_2_2, a_10_2_2, NORM_TWO) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::vectorNorm(a_10_2_2, a_10_2_2_3, NORM_TWO) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::vectorNorm(a_10_3, a_10_2_2, NORM_TWO) );

        *outStream << "-> add with multidimensional arrays:\n";

        INTREPID2_TEST_ERROR_EXPECTED( rst::add(a_10_2_2, a_10_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::add(a_10_2_3, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::add(a_10_2_2, a_10_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::add(a_10_2_3, a_10_2_2) );

        *outStream << "-> subtract with multidimensional arrays:\n";

        INTREPID2_TEST_ERROR_EXPECTED( rst::subtract(a_10_2_2, a_10_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::subtract(a_10_2_3, a_10_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::subtract(a_10_2_2, a_10_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::subtract(a_10_2_3, a_10_2_2) );
        
        *outStream << "-> dot product norm with multidimensional arrays:\n";
        
        INTREPID2_TEST_ERROR_EXPECTED( rst::dot(a_10_2, a_10_2_2_3, a_10_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::dot(a_10_2, a_10_2_2, a_10_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::dot(a_10_2_2, a_10_2_2_3, a_10_2_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::dot(a_10_2, a_10_2_2, a_10_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::dot(a_10_3, a_10_2_3, a_10_2_3) );

        *outStream << "-> absolute value with multidimensional arrays:\n";

        INTREPID2_TEST_ERROR_EXPECTED( rst::absval(a_10_3, a_10_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::absval(a_10_2_2, a_10_2_3) );
#endif
      } catch (std::logic_error err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      }
      
      *outStream                                \
        << "\n"
        << "===============================================================================\n" \
        << "| TEST 2: matrix / matrix-vector exceptions                                   |\n" \
        << "===============================================================================\n";

      try {
        DynRankView ConstructWithLabel(a_10_1_2_3_4, 10, 1, 2, 3, 4);
        DynRankView ConstructWithLabel(b_10_1_2_3_4, 10, 1, 2, 3, 4);
        DynRankView ConstructWithLabel(a_10, 10);
        DynRankView ConstructWithLabel(a_9, 9);
        DynRankView ConstructWithLabel(b_10, 10);
        DynRankView ConstructWithLabel(a_10_15_4_4, 10, 15, 4, 4);
        DynRankView ConstructWithLabel(b_10_15_4_4, 10, 15, 4, 4);
        DynRankView ConstructWithLabel(a_10_2_2, 10, 2, 2);
        DynRankView ConstructWithLabel(a_10_2_3, 10, 2, 3);
        DynRankView ConstructWithLabel(b_10_2_3, 10, 2, 3);
        
        DynRankView ConstructWithLabel(a_1_1, 1, 1);
        DynRankView ConstructWithLabel(b_1_1, 1, 1);
        DynRankView ConstructWithLabel(a_2_2, 2, 2);
        DynRankView ConstructWithLabel(b_2_2, 2, 2);
        DynRankView ConstructWithLabel(a_3_3, 3, 3);
        DynRankView ConstructWithLabel(b_3_3, 3, 3);
        DynRankView ConstructWithLabel(a_2_3, 2, 3);
        DynRankView ConstructWithLabel(a_4_4, 4, 4);
        
        DynRankView ConstructWithLabel(a_10_15_1_1, 10, 15, 1, 1);
        DynRankView ConstructWithLabel(b_10_15_1_1, 10, 15, 1, 1);
        DynRankView ConstructWithLabel(a_10_15_2_2, 10, 15, 2, 2);
        DynRankView ConstructWithLabel(b_10_15_2_2, 10, 15, 2, 2);
        DynRankView ConstructWithLabel(a_10_15_3_3, 10, 15, 3, 3);
        DynRankView ConstructWithLabel(a_10_15_3_2, 10, 15, 3, 2);
        DynRankView ConstructWithLabel(b_10_15_3_3, 10, 15, 3, 3);
        DynRankView ConstructWithLabel(b_10_14, 10, 14);
        DynRankView ConstructWithLabel(b_10_15, 10, 15);
        DynRankView ConstructWithLabel(b_10_14_3, 10, 14, 3);
        DynRankView ConstructWithLabel(b_10_15_3, 10, 15, 3);

#ifdef HAVE_INTREPID2_DEBUG
#define INTREPID2_TEST_FOR_ABORT_OVERRIDE_TO_CONTINUE
        *outStream << "-> Errors are expected due to rank mismatches:\n";
        
        *outStream << "-> inverse with multidimensional arrays:\n";
        
        INTREPID2_TEST_ERROR_EXPECTED( rst::inverse(a_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::inverse(b_10_1_2_3_4, a_10_1_2_3_4) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::inverse(b_10, a_10) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::inverse(a_10_2_2, a_10_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::inverse(b_10_2_3, a_10_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::inverse(b_10_15_4_4, a_10_15_4_4) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::inverse(b_1_1, a_1_1) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::inverse(b_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::inverse(b_3_3, a_3_3) );
        
        a_2_2(0, 0) = 1.0;
        a_3_3(0, 0) = 1.0;
        INTREPID2_TEST_ERROR_EXPECTED( rst::inverse(b_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::inverse(b_3_3, a_3_3) );
        
        *outStream << "-> transpose with multidimensional arrays:\n";

        INTREPID2_TEST_ERROR_EXPECTED( rst::transpose(a_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::transpose(b_10_1_2_3_4, a_10_1_2_3_4) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::transpose(b_10, a_10) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::transpose(a_10_2_2, a_10_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::transpose(b_10_2_3, a_10_2_3) );
        
        *outStream << "-> determinant with multidimensional arrays:\n";

        INTREPID2_TEST_ERROR_EXPECTED( rst::det(a_2_2, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::det(a_10_2_2, a_10_1_2_3_4) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::det(b_10_14, a_10_15_3_3) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::det(a_9, a_10_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::det(b_10, a_10_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::det(b_10_15, a_10_15_4_4) );
        rst::Serial::det(a_10_15_4_4);
        rst::Serial::det(a_2_3);
        rst::Serial::det(a_4_4);
        
        *outStream << "-> matrix-vector product with multidimensional arrays:\n";
        
        INTREPID2_TEST_ERROR_EXPECTED( rst::matvec(a_10_2_2, a_10_2_3, b_10_2_3) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::matvec(a_2_2, a_2_2, a_10) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::matvec(a_9, a_10_2_2, a_2_2) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::matvec(b_10_15_3, a_10_15_3_3, b_10_14_3) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::matvec(b_10_14_3, a_10_15_3_3, b_10_15_3) );
        INTREPID2_TEST_ERROR_EXPECTED( rst::matvec(b_10_15_3, a_10_15_3_2, b_10_15_3) );
#undef INTREPID2_TEST_FOR_ABORT_OVERRIDE_TO_CONTINUE 
#endif
      } catch (std::logic_error err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      }

      *outStream                                \
        << "\n"
        << "===============================================================================\n" \
        << "| TEST 2: correctness of math operations                                      |\n" \
        << "===============================================================================\n";
      
      outStream->precision(20);

      try {
        const value_type tol = tolerence()*100.0;
        
        // two-dimensional base containers
        for (auto dim=3; dim>0; dim--) {
          size_type i0=4, i1=5;
          DynRankView ConstructWithLabel(ma_x_x_d_d, i0, i1, dim, dim);
          DynRankView ConstructWithLabel(mb_x_x_d_d, i0, i1, dim, dim);
          DynRankView ConstructWithLabel(mc_x_x_d_d, i0, i1, dim, dim);

          DynRankView ConstructWithLabel(va_x_x_d, i0, i1, dim);
          DynRankView ConstructWithLabel(vb_x_x_d, i0, i1, dim);
          DynRankView ConstructWithLabel(vc_x_x_d, i0, i1, dim);
          DynRankView ConstructWithLabel(vdot_x_x, i0, i1);

          DynRankView ConstructWithLabel(vnorms_x_x, i0, i1);
          DynRankView ConstructWithLabel(vnorms_x, i0);

          
          // fill with random numbers
          {
            const auto iend = ma_x_x_d_d.extent(0);
            const auto jend = ma_x_x_d_d.extent(1);
            const auto kend = ma_x_x_d_d.extent(2);
            const auto lend = ma_x_x_d_d.extent(3);
            
            for (size_type i=0;i<iend;++i)
              for (size_type j=0;j<jend;++j)
                for (size_type k=0;k<kend;++k)
                  for (size_type l=0;l<lend;++l)
                    ma_x_x_d_d(i,j,k,l) = Teuchos::ScalarTraits<value_type>::random();
          }
          {
            const auto iend = va_x_x_d.extent(0);
            const auto jend = va_x_x_d.extent(1);
            const auto kend = va_x_x_d.extent(2);

            {
              for (size_type i=0;i<iend;++i)
                for (size_type j=0;j<jend;++j)
                  for (size_type k=0;k<kend;++k)
                    va_x_x_d(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
            }
          }
          
          *outStream << "\n-- Checking vectorNorm \n";

          rst::vectorNorm(vnorms_x_x, va_x_x_d, NORM_TWO);
          if ( std::abs(rst::Serial::vectorNorm(vnorms_x_x, NORM_TWO) - 
                        rst::Serial::vectorNorm(va_x_x_d, NORM_TWO)) > tol) {
            *outStream << "\n\nINCORRECT vectorNorm NORM_TWO\n\n";
            errorFlag = -1000;
          }
          
          rst::vectorNorm(vnorms_x_x, va_x_x_d, NORM_ONE);
          if ( std::abs(rst::Serial::vectorNorm(vnorms_x_x, NORM_ONE) - 
                        rst::Serial::vectorNorm(va_x_x_d, NORM_ONE)) > tol) {
            *outStream << "\n\nINCORRECT vectorNorm NORM_ONE\n\n";
            errorFlag = -1000;
          }
          
          rst::vectorNorm(vnorms_x_x, va_x_x_d, NORM_INF);
          if ( std::abs(rst::Serial::vectorNorm(vnorms_x_x, NORM_INF) - 
                        rst::Serial::vectorNorm(va_x_x_d, NORM_INF)) > tol) {
            *outStream << "\n\nINCORRECT vectorNorm NORM_INF\n\n";
            errorFlag = -1000;
          }
          
          *outStream << "\n-- Checking inverse, subtract, and vectorNorm \n";
          
          rst::inverse(mb_x_x_d_d, ma_x_x_d_d); // B = inv(A)
          rst::inverse(mc_x_x_d_d, mb_x_x_d_d); // C = inv(B) ~= A
          
          rst::subtract(mc_x_x_d_d, ma_x_x_d_d); // C = C - A ~= 0 
          
          if (rst::Serial::vectorNorm(mc_x_x_d_d, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT inverse OR subtract OR vectorNorm\n\n";
            errorFlag = -1000;
          }
          
          *outStream << "\n-- Checking determinant \n";
          
          DynRankView ConstructWithLabel(detA_x_x, i0, i1);
          DynRankView ConstructWithLabel(detB_x_x, i0, i1);
          
          rst::det(detA_x_x, ma_x_x_d_d);
          rst::det(detB_x_x, mb_x_x_d_d);
          
          { 
            value_type val = 0;
            for (size_type i=0;i<detA_x_x.extent(0);++i) 
              val += rst::Serial::dot( Kokkos::subdynrankview(detA_x_x, i, Kokkos::ALL()),
                                       Kokkos::subdynrankview(detB_x_x, i, Kokkos::ALL()) );

            const value_type diff = std::abs(val - value_type(i0*i1));
            if (diff  > tol) {
              *outStream << "\n\nINCORRECT det\n\n" ;
              errorFlag = -1000;
            }
          }
          
          *outStream << "\n-- Checking transpose and subtract \n";
          
          rst::transpose(mb_x_x_d_d, ma_x_x_d_d); // B = A^T
          rst::transpose(mc_x_x_d_d, mb_x_x_d_d); // C = B^T = A
          
          rst::subtract(mc_x_x_d_d, ma_x_x_d_d); // C = C - A = 0 
          
          if (rst::Serial::vectorNorm(mc_x_x_d_d, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT transpose OR subtract OR vectorNorm\n\n" ;
            errorFlag = -1000;
          }
          
          *outStream << "\n-- Checking matvec, vectorNorm, subtract, and inverse \n";
          
          rst::inverse(mb_x_x_d_d, ma_x_x_d_d); // B = inv(A)
          rst::inverse(mc_x_x_d_d, mb_x_x_d_d); // C = inv(B) ~= A
          rst::matvec(vb_x_x_d, ma_x_x_d_d, va_x_x_d); // b = A*a
          rst::matvec(vc_x_x_d, mb_x_x_d_d, vb_x_x_d); // c = inv(A)*(A*a) ~= a
          rst::subtract(vc_x_x_d, va_x_x_d); // c = c - a ~= 0
          
          rst::vectorNorm(vnorms_x_x, vc_x_x_d, NORM_ONE);
          rst::vectorNorm(vnorms_x, vnorms_x_x, NORM_INF);
          if (rst::Serial::vectorNorm(vnorms_x, NORM_TWO) > tol) {
            *outStream << "\n\nINCORRECT matvec OR inverse OR subtract OR vectorNorm\n\n";
            errorFlag = -1000;
          }
          
          *outStream << "\n-- Checking add, subtract, absval, and scale \n";
          
          const value_type x = 1.234;
          rst::add(vc_x_x_d, va_x_x_d, vb_x_x_d); // c = a + b
          rst::subtract(vc_x_x_d, vb_x_x_d); // c = c - b = a

          rst::scale(vb_x_x_d, vc_x_x_d, x); // b = c*x;
          rst::scale(vc_x_x_d, vb_x_x_d, (1.0/x)); // c = b*(1/x) = a;
          rst::subtract(vb_x_x_d, vc_x_x_d, va_x_x_d); // b = c - a ~= 0

          rst::absval(vc_x_x_d, vb_x_x_d); // c = |b|
          rst::scale(vb_x_x_d, vc_x_x_d, -1.0); // b = -c
          rst::absval(vc_x_x_d, vb_x_x_d); // c = |b|
          rst::add(vc_x_x_d, vb_x_x_d); // c = c + b === 0
          
          rst::vectorNorm(vnorms_x_x, vc_x_x_d, NORM_ONE);
          rst::vectorNorm(vnorms_x, vnorms_x_x, NORM_INF);
          if (rst::Serial::vectorNorm(vnorms_x, NORM_TWO) > tol) {
            *outStream << "\n\nSign flips combined with std::abs might not be invertible on this platform!\n"
                       << "Potential IEEE compliance issues!\n\n";
            if (rst::Serial::vectorNorm(vnorms_x, NORM_TWO) > tol) {
              *outStream << "\n\nINCORRECT add OR subtract OR scale OR absval OR vectorNorm\n\n";
              errorFlag = -1000;
            }
          }
          
          *outStream << "\n-- Checking dot and vectorNorm \n";
          
          // fill with random numbers
          {
            const auto iend = va_x_x_d.extent(0);
            const auto jend = va_x_x_d.extent(1);
            const auto kend = va_x_x_d.extent(2);
            
            for (size_type i=0;i<iend;++i)
              for (size_type j=0;j<jend;++j)
                for (size_type k=0;k<kend;++k)
                  va_x_x_d(i,j,k) = 2.0;
          }
          
          rst::dot(vdot_x_x, va_x_x_d, va_x_x_d); // dot = a'*a

          rst::vectorNorm(vnorms_x, vdot_x_x, NORM_ONE);
          if (rst::Serial::vectorNorm(vnorms_x, NORM_ONE) - (value_type)(4.0*dim*i0*i1) > tol) {
            *outStream << "\n\nINCORRECT dot OR vectorNorm\n\n";
            errorFlag = -1000;
          }
          
          *outStream << "\n";
        }
    
        // one-dimensional base containers
        for (auto dim=3;dim>0;--dim) {
          ordinal_type i0=7;
          DynRankView ConstructWithLabel(ma_x_d_d, i0, dim, dim);
          DynRankView ConstructWithLabel(mb_x_d_d, i0, dim, dim);
          DynRankView ConstructWithLabel(mc_x_d_d, i0, dim, dim);
          DynRankView ConstructWithLabel(va_x_d, i0, dim);
          DynRankView ConstructWithLabel(vb_x_d, i0, dim);
          DynRankView ConstructWithLabel(vc_x_d, i0, dim);
          DynRankView ConstructWithLabel(vdot_x, i0);
          DynRankView ConstructWithLabel(vnorms_x, i0);
          
          // fill with random numbers
          // fill with random numbers
          {
            const auto iend = ma_x_d_d.extent(0);
            const auto jend = ma_x_d_d.extent(1);
            const auto kend = ma_x_d_d.extent(2);
            const auto lend = ma_x_d_d.extent(3);
            
            for (size_type i=0;i<iend;++i)
              for (size_type j=0;j<jend;++j)
                for (size_type k=0;k<kend;++k)
                  for (size_type l=0;l<lend;++l)
                    ma_x_d_d(i,j,k,l) = Teuchos::ScalarTraits<value_type>::random();
          }
          {
            const auto iend = va_x_d.extent(0);
            const auto jend = va_x_d.extent(1);
            const auto kend = va_x_d.extent(2);
            
            for (size_type i=0;i<iend;++i)
              for (size_type j=0;j<jend;++j)
                for (size_type k=0;k<kend;++k)
                  va_x_d(i,j,k) = Teuchos::ScalarTraits<value_type>::random();
          }
          
          *outStream << "\n-- Checking vectorNorm \n";
          
          rst::vectorNorm(vnorms_x, va_x_d, NORM_TWO);
          if ( std::abs(rst::Serial::vectorNorm(vnorms_x, NORM_TWO) - 
                        rst::Serial::vectorNorm(va_x_d, NORM_TWO)) > tol) {
            *outStream << "\n\nINCORRECT vectorNorm NORM_TWO\n\n";
            errorFlag = -1000;
          }
          
          rst::vectorNorm(vnorms_x, va_x_d, NORM_ONE);
          if ( std::abs(rst::Serial::vectorNorm(vnorms_x, NORM_ONE) - 
                        rst::Serial::vectorNorm(va_x_d, NORM_ONE)) > tol) {
            *outStream << "\n\nINCORRECT vectorNorm NORM_ONE\n\n";
            errorFlag = -1000;
          }
          
          rst::vectorNorm(vnorms_x, va_x_d, NORM_INF);
          if ( std::abs(rst::Serial::vectorNorm(vnorms_x, NORM_INF) - 
                        rst::Serial::vectorNorm(va_x_d, NORM_INF)) > tol) {
            *outStream << "\n\nINCORRECT vectorNorm NORM_INF\n\n";
            errorFlag = -1000;
          }
          
          *outStream << "\n-- Checking inverse, subtract, and vectorNorm \n";
          
          rst::inverse(mb_x_d_d, ma_x_d_d); // B = inv(A)
          rst::inverse(mc_x_d_d, mb_x_d_d); // C = inv(B) ~= A
          rst::subtract(mc_x_d_d, ma_x_d_d); // C = C - A ~= 0 
          
          if (rst::Serial::vectorNorm(mc_x_d_d, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT inverse OR subtract OR vectorNorm\n\n";
            errorFlag = -1000;
          }
          
          *outStream << "\n-- Checking determinant\n";
          
          DynRankView ConstructWithLabel(detA_x, i0);
          DynRankView ConstructWithLabel(detB_x, i0);
          
          rst::det(detA_x, ma_x_d_d);
          rst::det(detB_x, mb_x_d_d);
          
          if ( (rst::Serial::dot(detA_x, detB_x) - (value_type)i0) > tol) {
            *outStream << "\n\nINCORRECT det\n\n" ;
            errorFlag = -1000;
          }
          
          *outStream << "\n-- Checking transpose and subtract \n";
          
          rst::transpose(mb_x_d_d, ma_x_d_d); // B = A^T
          rst::transpose(mc_x_d_d, mb_x_d_d); // C = B^T = A
          rst::subtract(mc_x_d_d, ma_x_d_d); // C = C - A = 0 
          
          if (rst::Serial::vectorNorm(mc_x_d_d, NORM_ONE) > tol) {
            *outStream << "\n\nINCORRECT transpose OR subtract OR vectorNorm\n\n" ;
            errorFlag = -1000;
          }
          
          *outStream << "\n-- Checking matvec, vectorNorm, subtract, and inverse \n";
          
          rst::inverse(mb_x_d_d, ma_x_d_d); // B = inv(A)
          rst::inverse(mc_x_d_d, mb_x_d_d); // C = inv(B) ~= A
          rst::matvec(vb_x_d, ma_x_d_d, va_x_d); // b = A*a
          rst::matvec(vc_x_d, mb_x_d_d, vb_x_d); // c = inv(A)*(A*a) ~= a
          rst::subtract(vc_x_d, va_x_d); // c = c - a ~= 0
          
          rst::vectorNorm(vnorms_x, vc_x_d, NORM_ONE);
          if (rst::Serial::vectorNorm(vnorms_x, NORM_TWO) > tol) {
            *outStream << "\n\nINCORRECT matvec OR inverse OR subtract OR vectorNorm\n\n";
            errorFlag = -1000;
          }
          
          *outStream << "\n-- Checking add, subtract, absval, and scale \n";
          
          const value_type x = 1.234;
          rst::add(vc_x_d, va_x_d, vb_x_d); // c = a + b
          rst::subtract(vc_x_d, vb_x_d); // c = c - b = a
          rst::scale(vb_x_d, vc_x_d, x); // b = c*x;
          rst::scale(vc_x_d, vb_x_d, (1.0/x)); // c = b*(1/x) = a;
          rst::subtract(vb_x_d, vc_x_d, va_x_d); // b = c - a ~= 0
          rst::absval(vc_x_d, vb_x_d); // c = |b|
          rst::scale(vb_x_d, vc_x_d, -1.0); // b = -c
          rst::absval(vc_x_d, vb_x_d); // c = |b|
          rst::add(vc_x_d, vb_x_d); // c = c + b === 0
          
          rst::vectorNorm(vnorms_x, vc_x_d, NORM_ONE);
          if (rst::Serial::vectorNorm(vnorms_x, NORM_TWO) > tol) {
            *outStream << "\n\nSign flips combined with std::abs might not be invertible on this platform!\n"
                       << "Potential IEEE compliance issues!\n\n";
            if (rst::Serial::vectorNorm(vnorms_x, NORM_TWO) > tol) {
              *outStream << "\n\nINCORRECT add OR subtract OR scale OR absval OR vectorNorm\n\n";
              errorFlag = -1000;
            }
          }
          
          *outStream << "\n-- Checking dot and vectorNorm \n";

          {
            const auto iend = va_x_d.extent(0);
            const auto jend = va_x_d.extent(1);
            
            for (size_type i=0;i<iend;++i)
              for (size_type j=0;j<jend;++j)
                va_x_d(i,j) = 2.0;
          }
          rst::dot(vdot_x, va_x_d, va_x_d); // dot = a'*a
          
          if (rst::Serial::vectorNorm(vdot_x, NORM_ONE) - (double)(4.0*dim*i0) > tol) {
            *outStream << "\n\nINCORRECT dot OR vectorNorm\n\n";
            errorFlag = -1000;
          }
          
          *outStream << "\n";
        }
      } catch (std::logic_error err) {
        *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
        *outStream << err.what() << '\n';
        *outStream << "-------------------------------------------------------------------------------" << "\n\n";
        errorFlag = -1000;
      }
      
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
