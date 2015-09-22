// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file test_01.cpp
\brief  Unit test for the RealSpaceTools class.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_FieldContainer.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include <Kokkos_Core.hpp>


using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
Kokkos::initialize();

  // This little trick lets us print to std::cout only if
  // a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  // Save the format state of the original std::cout.
  Teuchos::oblackholestream oldFormatState;
  oldFormatState.copyfmt(std::cout);

  *outStream \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                       Unit Test (RealSpaceTools)                            |\n" \
  << "|                                                                             |\n" \
  << "|     1) Vector operations in 1D, 2D, and 3D real space                       |\n" \
  << "|     2) Matrix / matrix-vector operations in 1D, 2D, and 3D real space       |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n";


  typedef RealSpaceTools<double> rst;


  int errorFlag = 0;
#ifdef HAVE_INTREPID_DEBUG
  int beginThrowNumber = Teuchos::TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 49;
#ifndef HAVE_INTREPID_DEBUG_INF_CHECK
  endThrowNumber = beginThrowNumber + 44;
#endif

#endif
  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 1: vector exceptions                                                   |\n"\
  << "===============================================================================\n";

  try{

    Kokkos::View<double**> a_2_2("a_2_2",2, 2);
    Kokkos::View<double**> a_10_2("a_10_2",10, 2);
    Kokkos::View<double**> a_10_3("a_10_3",10, 3);
    Kokkos::View<double***> a_10_2_2("a_10_2_2",10, 2, 2);
    Kokkos::View<double***> a_10_2_3("a_10_2_3",10, 2, 3);
    Kokkos::View<double****> a_10_2_2_3("a_10_2_2_3",10, 2, 2, 3);

#ifdef HAVE_INTREPID_DEBUG

    *outStream << "-> vector norm with multidimensional arrays:\n";

    try {
      rst::vectorNorm(a_2_2, NORM_TWO);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::vectorNorm(a_10_2_2, a_10_2_2, NORM_TWO);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::vectorNorm(a_10_2_2, a_10_2_2_3, NORM_TWO);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::vectorNorm(a_10_3, a_10_2_2, NORM_TWO);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

    *outStream << "-> add with multidimensional arrays:\n";

    try {
      rst::add(a_10_2_2, a_10_2, a_2_2);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::add(a_10_2_3, a_10_2_2, a_10_2_2);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::add(a_10_2_2, a_10_2_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::add(a_10_2_3, a_10_2_2);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

    *outStream << "-> subtract with multidimensional arrays:\n";

    try {
      rst::subtract(a_10_2_2, a_10_2, a_2_2);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::subtract(a_10_2_3, a_10_2_2, a_10_2_2);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::subtract(a_10_2_2, a_10_2_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::subtract(a_10_2_3, a_10_2_2);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

    *outStream << "-> dot product norm with multidimensional arrays:\n";

    try {
      rst::dot(a_10_2, a_10_2_2_3, a_10_2_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::dot(a_10_2, a_10_2_2, a_10_2_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::dot(a_10_2_2, a_10_2_2_3, a_10_2_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::dot(a_10_2, a_10_2_2, a_10_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::dot(a_10_3, a_10_2_3, a_10_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

    *outStream << "-> absolute value with multidimensional arrays:\n";

    try {
      rst::absval(a_10_3, a_10_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::absval(a_10_2_2, a_10_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
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
  << "| TEST 2: matrix / matrix-vector exceptions                                   |\n"\
  << "===============================================================================\n";

  try{

    Kokkos::View<double*****> a_10_1_2_3_4("a_10_1_2_3_4",10, 1, 2, 3, 4);
    Kokkos::View<double*****> b_10_1_2_3_4("b_10_1_2_3_4",10, 1, 2, 3, 4);
    Kokkos::View<double*> a_10("a_10",10);
    Kokkos::View<double*> a_9("a_9",9);
    Kokkos::View<double*> b_10("b_10",10);
    Kokkos::View<double****> a_10_15_4_4("a_10_15_4_4",10, 15, 4, 4);
    Kokkos::View<double****> b_10_15_4_4("b_10_15_4_4",10, 15, 4, 4);
    Kokkos::View<double***> a_10_2_2("a_10_2_2",10, 2, 2);
    Kokkos::View<double***> a_10_2_3("a_10_2_3",10, 2, 3);
    Kokkos::View<double***> b_10_2_3("b_10_2_3",10, 2, 3);

    Kokkos::View<double**> a_1_1("a_1_1",1, 1);
    Kokkos::View<double**> b_1_1("b_1_1",1, 1);
    Kokkos::View<double**> a_2_2("a_2_2",2, 2);
    Kokkos::View<double**> b_2_2("b_2_2",2, 2);
    Kokkos::View<double**> a_3_3("a_3_3",3, 3);
    Kokkos::View<double**> b_3_3("b_3_3",3, 3);
    Kokkos::View<double**> a_2_3("a_2_3",2, 3);
    Kokkos::View<double**> a_4_4("a_4_4",4, 4);

    Kokkos::View<double****> a_10_15_1_1("a_10_15_1_1",10, 15, 1, 1);
    Kokkos::View<double****> b_10_15_1_1("b_10_15_1_1",10, 15, 1, 1);
    Kokkos::View<double****> a_10_15_2_2("a_10_15_2_2",10, 15, 2, 2);
    Kokkos::View<double****> b_10_15_2_2("b_10_15_2_2",10, 15, 2, 2);
    Kokkos::View<double****> a_10_15_3_3("a_10_15_3_3",10, 15, 3, 3);
    Kokkos::View<double****> a_10_15_3_2("a_10_15_3_2",10, 15, 3, 2);
    Kokkos::View<double****> b_10_15_3_3("b_10_15_3_3",10, 15, 3, 3);
    Kokkos::View<double**> b_10_14("b_10_14",10, 14);
    Kokkos::View<double**> b_10_15("b_10_15",10, 15);
    Kokkos::View<double***> b_10_14_3("b_10_14_3",10, 14, 3);
    Kokkos::View<double***> b_10_15_3("b_10_15_3",10, 15, 3);
    
   
#ifdef HAVE_INTREPID_DEBUG

    *outStream << "-> inverse with multidimensional arrays:\n";

    try {
      rst::inverse(a_2_2, a_10_2_2);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::inverse(b_10_1_2_3_4, a_10_1_2_3_4);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::inverse(b_10, a_10);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::inverse(a_10_2_2, a_10_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::inverse(b_10_2_3, a_10_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::inverse(b_10_15_4_4, a_10_15_4_4);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::inverse(b_1_1, a_1_1);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::inverse(b_2_2, a_2_2);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::inverse(b_3_3, a_3_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    a_2_2(0,0) = 1.0;
    a_3_3(0,0) = 1.0;
    try {
      rst::inverse(b_2_2, a_2_2);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::inverse(b_3_3, a_3_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

    *outStream << "-> transpose with multidimensional arrays:\n";

    try {
      rst::transpose(a_2_2, a_10_2_2);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::transpose(b_10_1_2_3_4, a_10_1_2_3_4);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::transpose(b_10, a_10);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::transpose(a_10_2_2, a_10_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::transpose(b_10_2_3, a_10_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

    *outStream << "-> determinant with multidimensional arrays:\n";

    try {
      rst::det(a_2_2, a_10_2_2);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::det(a_10_2_2, a_10_1_2_3_4);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::det(b_10_14, a_10_15_3_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::det(a_9, a_10_2_2);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::det(b_10, a_10_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::det(b_10_15, a_10_15_4_4);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::det(a_10_15_4_4);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::det(a_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::det(a_4_4);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

    *outStream << "-> matrix-vector product with multidimensional arrays:\n";

    try {
      rst::matvec(a_10_2_2, a_10_2_3, b_10_2_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::matvec(a_2_2, a_2_2, a_10);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::matvec(a_9, a_10_2_2, a_2_2);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::matvec(b_10_15_3, a_10_15_3_3, b_10_14_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::matvec(b_10_14_3, a_10_15_3_3, b_10_15_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
    try {
      rst::matvec(b_10_15_3, a_10_15_3_2, b_10_15_3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };

#endif

  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1000;
  };
#ifdef HAVE_INTREPID_DEBUG
  if (Teuchos::TestForException_getThrowNumber() != endThrowNumber)
    errorFlag++;
#endif

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 2: correctness of math operations                                      |\n"\
  << "===============================================================================\n";

  outStream->precision(20);

  try{
 
     // two-dimensional base containers
     for (int dim=3; dim>0; dim--) {
      int i0=4, i1=5;
      Kokkos::View<double****> ma_x_x_d_d("ma_x_x_d_d",i0, i1, dim, dim);
      Kokkos::View<double****> mb_x_x_d_d("mb_x_x_d_d",i0, i1, dim, dim);
      Kokkos::View<double****> mc_x_x_d_d("mc_x_x_d_d",i0, i1, dim, dim);
      Kokkos::View<double***> va_x_x_d("va_x_x_d",i0, i1, dim);
      Kokkos::View<double***> vb_x_x_d("vb_x_x_d",i0, i1, dim);
      Kokkos::View<double***> vc_x_x_d("vc_x_x_d",i0, i1, dim);
      Kokkos::View<double**> vdot_x_x("vdot_x_x",i0, i1);
      Kokkos::View<double**> vnorms_x_x("vnorms_x_x",i0, i1);
      Kokkos::View<double*> vnorms_x("vnorms_x",i0);
      double zero = INTREPID_TOL*100.0;

      // fill with random numbers
       for (unsigned int i=0; i<ma_x_x_d_d.dimension(0); i++) {
		   for (unsigned int j=0; j<ma_x_x_d_d.dimension(1); j++) {
			   for (unsigned int k=0; k<ma_x_x_d_d.dimension(2); k++) {
				   for (unsigned int l=0; l<ma_x_x_d_d.dimension(3); l++) {
					 ma_x_x_d_d(i,j,k,l) = Teuchos::ScalarTraits<double>::random();
  
				   }
			   }
		   }
	   }
       for (unsigned int i=0; i<va_x_x_d.dimension(0); i++) {
		   for (unsigned int j=0; j<va_x_x_d.dimension(1); j++) {
			   for (unsigned int k=0; k<va_x_x_d.dimension(2); k++) {
               va_x_x_d(i,j,k) = Teuchos::ScalarTraits<double>::random();
				   }
			   }
		   }

      *outStream << "\n************ Checking vectorNorm ************\n";
     
      rst::vectorNorm(vnorms_x_x, va_x_x_d, NORM_TWO);
      *outStream << "va_x_x_d";
      *outStream << "vnorms_x_x";
      if ( std::abs(rst::vectorNorm(vnorms_x_x, NORM_TWO) - 
                    rst::vectorNorm(va_x_x_d, NORM_TWO)) > zero) {
        *outStream << "\n\nINCORRECT vectorNorm NORM_TWO\n\n";
        errorFlag = -1000;
      }
    
      rst::vectorNorm(vnorms_x_x, va_x_x_d, NORM_ONE);
      *outStream << "va_x_x_d";
      *outStream << "vnorms_x_x";
      if ( std::abs(rst::vectorNorm(vnorms_x_x, NORM_ONE) - 
                    rst::vectorNorm(va_x_x_d, NORM_ONE)) > zero) {
        *outStream << "\n\nINCORRECT vectorNorm NORM_ONE\n\n";
        errorFlag = -1000;
      }
      rst::vectorNorm(vnorms_x_x, va_x_x_d, NORM_INF);
      *outStream << "va_x_x_d";
      *outStream << "vnorms_x_x";
      if ( std::abs(rst::vectorNorm(vnorms_x_x, NORM_INF) - 
                    rst::vectorNorm(va_x_x_d, NORM_INF)) > zero) {
        *outStream << "\n\nINCORRECT vectorNorm NORM_INF\n\n";
        errorFlag = -1000;
      }

      /******************************************/


      *outStream << "\n************ Checking inverse, subtract, and vectorNorm ************\n";
     
      rst::inverse(mb_x_x_d_d, ma_x_x_d_d); // B = inv(A)
      rst::inverse(mc_x_x_d_d, mb_x_x_d_d); // C = inv(B) ~= A
      *outStream << "ma_x_x_d_d" << "mb_x_x_d_d" << "mc_x_x_d_d";

      rst::subtract(mc_x_x_d_d, ma_x_x_d_d); // C = C - A ~= 0 

      if (rst::vectorNorm(mc_x_x_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT inverse OR subtract OR vectorNorm\n\n";
        errorFlag = -1000;
      }

      /******************************************/


      *outStream << "\n********** Checking determinant **********\n";

      Kokkos::View<double**> detA_x_x("detA_x_x",i0, i1);
      Kokkos::View<double**> detB_x_x("detB_x_x",i0, i1);
      Kokkos::View<double*> detA_x_xlinear("detA_x_xlinear",i0*i1);
      Kokkos::View<double*> detB_x_xlinear("detB_x_xlinear",i0*i1);

      rst::det(detA_x_x, ma_x_x_d_d);
      rst::det(detB_x_x, mb_x_x_d_d);
      *outStream << "detA_x_x" << "detB_x_x";
       for(int i=0;i<i0;i++){
		  for(int j=0;j<i1;j++){
			detA_x_xlinear(i1*i+j)=detA_x_x(i,j);
			detB_x_xlinear(i1*i+j)=detB_x_x(i,j);
		  }
	  }
      if ( (rst::dot(detA_x_xlinear, detB_x_xlinear) - (double)(i0*i1)) > zero) {
        *outStream << "\n\nINCORRECT det\n\n" ;
        errorFlag = -1000;
      }
  /*    Kokkos::View<double*> ma_x_x_d_dlinear("ma_x_x_d_d",i0*i1*dim*dim);
      Kokkos::View<double*> mb_x_x_d_dlinear("mb_x_x_d_d",i0*i1*dim*dim);  
      for(int i=0;i<i0;i++){
		  for(int j=0;j<i1;j++){
			  for(int k=0;k<dim;k++){
		          for(int l=0;l<dim;l++){
		ma_x_x_d_dlinear(i*i1*dim*dim+j*dim*dim+k*dim+l)=ma_x_x_d_d(i,j,k,l);			  
		mb_x_x_d_dlinear(i*i1*dim*dim+j*dim*dim+k*dim+l)=mb_x_x_d_d(i,j,k,l);			  
				  }	  
			  }
		  }
	  }
      *outStream << "\n det(A)*det(inv(A)) = " <<
                    rst::det(ma_x_x_d_dlinear)*rst::det(mb_x_x_d_dlinear)
                 << "\n";


      if ( (rst::det(ma_x_x_d_dlinear)*
            rst::det(mb_x_x_d_dlinear) - (double)1) > zero) {
        *outStream << "\n\nINCORRECT det\n\n" ;
        errorFlag = -1000;
      }*/

      /******************************************/
   

      *outStream << "\n************ Checking transpose and subtract ************\n";
     
      rst::transpose(mb_x_x_d_d, ma_x_x_d_d); // B = A^T
      rst::transpose(mc_x_x_d_d, mb_x_x_d_d); // C = B^T = A
      *outStream << "ma_x_x_d_d" << "mb_x_x_d_d" << "mc_x_x_d_d";

      rst::subtract(mc_x_x_d_d, ma_x_x_d_d); // C = C - A = 0 

      if (rst::vectorNorm(mc_x_x_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT transpose OR subtract OR vectorNorm\n\n" ;
        errorFlag = -1000;
      }

      /******************************************/


      *outStream << "\n************ Checking matvec, vectorNorm, subtract, and inverse ************\n";

      rst::inverse(mb_x_x_d_d, ma_x_x_d_d); // B = inv(A)
      rst::inverse(mc_x_x_d_d, mb_x_x_d_d); // C = inv(B) ~= A
      rst::matvec(vb_x_x_d, ma_x_x_d_d, va_x_x_d); // b = A*a
      rst::matvec(vc_x_x_d, mb_x_x_d_d, vb_x_x_d); // c = inv(A)*(A*a) ~= a
      rst::subtract(vc_x_x_d, va_x_x_d); // c = c - a ~= 0
      *outStream << "vc_x_x_d";

      rst::vectorNorm(vnorms_x_x, vc_x_x_d, NORM_ONE);
      rst::vectorNorm(vnorms_x, vnorms_x_x, NORM_INF);
      if (rst::vectorNorm(vnorms_x, NORM_TWO) > zero) {
        *outStream << "\n\nINCORRECT matvec OR inverse OR subtract OR vectorNorm\n\n";
        errorFlag = -1000;
      }
     
      /******************************************/


      *outStream << "\n************ Checking add, subtract, absval, and scale ************\n";

      double x = 1.234;
      rst::add(vc_x_x_d, va_x_x_d, vb_x_x_d); // c = a + b
      rst::subtract(vc_x_x_d, vb_x_x_d); // c = c - b = a
      rst::scale(vb_x_x_d, vc_x_x_d, x); // b = c*x;
      rst::scale(vc_x_x_d, vb_x_x_d, (1.0/x)); // c = b*(1/x) = a;
      rst::subtract(vb_x_x_d, vc_x_x_d, va_x_x_d); // b = c - a ~= 0
      rst::absval(vc_x_x_d, vb_x_x_d); // c = |b|
      rst::scale(vb_x_x_d, vc_x_x_d, -1.0); // b = -c
      rst::absval(vc_x_x_d, vb_x_x_d); // c = |b|
      rst::add(vc_x_x_d, vb_x_x_d); // c = c + b === 0
      *outStream << "vc_x_x_d";

      rst::vectorNorm(vnorms_x_x, vc_x_x_d, NORM_ONE);
      rst::vectorNorm(vnorms_x, vnorms_x_x, NORM_INF);
      if (rst::vectorNorm(vnorms_x, NORM_TWO) > (double)0) {
        *outStream << "\n\nSign flips combined with std::abs might not be invertible on this platform!\n"
                   << "Potential IEEE compliance issues!\n\n";
        if (rst::vectorNorm(vnorms_x, NORM_TWO) > zero) {
          *outStream << "\n\nINCORRECT add OR subtract OR scale OR absval OR vectorNorm\n\n";
          errorFlag = -1000;
        }
      }

      /******************************************/


      *outStream << "\n************ Checking dot and vectorNorm ************\n";
      for (unsigned int i=0; i<va_x_x_d.dimension(0); i++) {
		   for (unsigned int j=0; j<va_x_x_d.dimension(1); j++) {
			   for (unsigned int k=0; k<va_x_x_d.dimension(2); k++) {
			va_x_x_d(i,j,k) = 2.0;
			   }
		   }
	   }

    
      rst::dot(vdot_x_x, va_x_x_d, va_x_x_d); // dot = a'*a
      *outStream << "vdot_x_x";
   
      rst::vectorNorm(vnorms_x, vdot_x_x, NORM_ONE);
      if (rst::vectorNorm(vnorms_x, NORM_ONE) - (double)(4.0*dim*i0*i1) > zero) {
        *outStream << "\n\nINCORRECT dot OR vectorNorm\n\n";
        errorFlag = -1000;
      }
    
      /******************************************/

      *outStream << "\n";
    }
    
    // one-dimensional base containers
    for (int dim=3; dim>0; dim--) {
      int i0=7;
      Kokkos::View<double***> ma_x_d_d("ma_x_d_d",i0, dim, dim);
      Kokkos::View<double***> mb_x_d_d("mb_x_d_d",i0, dim, dim);
      Kokkos::View<double***> mc_x_d_d("mc_x_d_d",i0, dim, dim);
      Kokkos::View<double**> va_x_d("va_x_d",i0, dim);
      Kokkos::View<double**> vb_x_d("vb_x_d",i0, dim);
      Kokkos::View<double**> vc_x_d("vc_x_d",i0, dim);
      Kokkos::View<double*> vdot_x("vdot_x",i0);
      Kokkos::View<double*> vnorms_x("vnorms_x",i0);
      double zero = INTREPID_TOL*100.0;

      // fill with random numbers
      
      for (unsigned int i=0; i<ma_x_d_d.dimension(0); i++) {
		   for (unsigned int j=0; j<ma_x_d_d.dimension(1); j++) {
			   for (unsigned int k=0; k<ma_x_d_d.dimension(2); k++) {
              ma_x_d_d(i,j,k) = Teuchos::ScalarTraits<double>::random();
				   
			   }
		   }
	   }      
      for (unsigned int i=0; i<va_x_d.dimension(0); i++) {
		   for (unsigned int j=0; j<va_x_d.dimension(1); j++) {
        va_x_d(i,j) = Teuchos::ScalarTraits<double>::random();
			 			   
		   }
	   }

    
      *outStream << "\n************ Checking vectorNorm ************\n";
     
      rst::vectorNorm(vnorms_x, va_x_d, NORM_TWO);
      *outStream << "va_x_d";
      *outStream << "vnorms_x";
      if ( std::abs(rst::vectorNorm(vnorms_x, NORM_TWO) - 
                    rst::vectorNorm(va_x_d, NORM_TWO)) > zero) {
        *outStream << "\n\nINCORRECT vectorNorm NORM_TWO\n\n";
        errorFlag = -1000;
      }

      rst::vectorNorm(vnorms_x, va_x_d, NORM_ONE);
      *outStream << "va_x_d";
      *outStream << "vnorms_x";
      if ( std::abs(rst::vectorNorm(vnorms_x, NORM_ONE) - 
                    rst::vectorNorm(va_x_d, NORM_ONE)) > zero) {
        *outStream << "\n\nINCORRECT vectorNorm NORM_ONE\n\n";
        errorFlag = -1000;
      }

      rst::vectorNorm(vnorms_x, va_x_d, NORM_INF);
      *outStream << "va_x_d";
      *outStream << "vnorms_x";
      if ( std::abs(rst::vectorNorm(vnorms_x, NORM_INF) - 
                    rst::vectorNorm(va_x_d, NORM_INF)) > zero) {
        *outStream << "\n\nINCORRECT vectorNorm NORM_INF\n\n";
        errorFlag = -1000;
      }

      /******************************************/

    
      *outStream << "\n************ Checking inverse, subtract, and vectorNorm ************\n";
     
      rst::inverse(mb_x_d_d, ma_x_d_d); // B = inv(A)
      rst::inverse(mc_x_d_d, mb_x_d_d); // C = inv(B) ~= A
      *outStream << "ma_x_d_d" << "mb_x_d_d" << "mc_x_d_d";

      rst::subtract(mc_x_d_d, ma_x_d_d); // C = C - A ~= 0 

      if (rst::vectorNorm(mc_x_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT inverse OR subtract OR vectorNorm\n\n";
        errorFlag = -1000;
      }

      /******************************************/


      *outStream << "\n********** Checking determinant **********\n";

      FieldContainer<double> detA_x(i0);
      FieldContainer<double> detB_x(i0);
      
      rst::det(detA_x, ma_x_d_d);
      rst::det(detB_x, mb_x_d_d);
      *outStream << "detA_x" << "detB_x";
      
      if ( (rst::dot(detA_x, detB_x) - (double)i0) > zero) {
        *outStream << "\n\nINCORRECT det\n\n" ;
        errorFlag = -1000;
      }
/*
      *outStream << "\n det(A)*det(inv(A)) = " <<
                    rst::det(ma_x_d_d)*rst::det(mb_x_d_d)
                 << "\n";

      if ( (rst::det(ma_x_d_d)*
            rst::det(mb_x_d_d) - (double)1) > zero) {
        *outStream << "\n\nINCORRECT det\n\n" ;
        errorFlag = -1000;
      }
*/
      /******************************************/

   
      *outStream << "\n************ Checking transpose and subtract ************\n";
     
      rst::transpose(mb_x_d_d, ma_x_d_d); // B = A^T
      rst::transpose(mc_x_d_d, mb_x_d_d); // C = B^T = A
      *outStream << "ma_x_d_d" << "mb_x_d_d" << "mc_x_d_d";

      rst::subtract(mc_x_d_d, ma_x_d_d); // C = C - A = 0 

      if (rst::vectorNorm(mc_x_d_d, NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT transpose OR subtract OR vectorNorm\n\n" ;
        errorFlag = -1000;
      }

      /******************************************/


      *outStream << "\n************ Checking matvec, vectorNorm, subtract, and inverse ************\n";

      rst::inverse(mb_x_d_d, ma_x_d_d); // B = inv(A)
      rst::inverse(mc_x_d_d, mb_x_d_d); // C = inv(B) ~= A
      rst::matvec(vb_x_d, ma_x_d_d, va_x_d); // b = A*a
      rst::matvec(vc_x_d, mb_x_d_d, vb_x_d); // c = inv(A)*(A*a) ~= a
      rst::subtract(vc_x_d, va_x_d); // c = c - a ~= 0
      *outStream << "vc_x_d";

      rst::vectorNorm(vnorms_x, vc_x_d, NORM_ONE);
      if (rst::vectorNorm(vnorms_x, NORM_TWO) > zero) {
        *outStream << "\n\nINCORRECT matvec OR inverse OR subtract OR vectorNorm\n\n";
        errorFlag = -1000;
      }

      /******************************************/


      *outStream << "\n************ Checking add, subtract, absval, and scale ************\n";

      double x = 1.234;
      rst::add(vc_x_d, va_x_d, vb_x_d); // c = a + b
      rst::subtract(vc_x_d, vb_x_d); // c = c - b = a
      rst::scale(vb_x_d, vc_x_d, x); // b = c*x;
      rst::scale(vc_x_d, vb_x_d, (1.0/x)); // c = b*(1/x) = a;
      rst::subtract(vb_x_d, vc_x_d, va_x_d); // b = c - a ~= 0
      rst::absval(vc_x_d, vb_x_d); // c = |b|
      rst::scale(vb_x_d, vc_x_d, -1.0); // b = -c
      rst::absval(vc_x_d, vb_x_d); // c = |b|
      rst::add(vc_x_d, vb_x_d); // c = c + b === 0
      *outStream << "vc_x_d";
    
      rst::vectorNorm(vnorms_x, vc_x_d, NORM_ONE);
      if (rst::vectorNorm(vnorms_x, NORM_TWO) > (double)0) {
        *outStream << "\n\nSign flips combined with std::abs might not be invertible on this platform!\n"
                   << "Potential IEEE compliance issues!\n\n";
        if (rst::vectorNorm(vnorms_x, NORM_TWO) > zero) {
          *outStream << "\n\nINCORRECT add OR subtract OR scale OR absval OR vectorNorm\n\n";
          errorFlag = -1000;
        }
      }
 
      /******************************************/
    

      *outStream << "\n************ Checking dot and vectorNorm ************\n";
      for (unsigned int i=0; i<va_x_d.dimension(0); i++) {
		   for (unsigned int j=0; j<va_x_d.dimension(1); j++) {
		   va_x_d(i,j) = 2.0;
		   }
	   }
      rst::dot(vdot_x, va_x_d, va_x_d); // dot = a'*a
      *outStream << "vdot_x";

      if (rst::vectorNorm(vdot_x, NORM_ONE) - (double)(4.0*dim*i0) > zero) {
        *outStream << "\n\nINCORRECT dot OR vectorNorm\n\n";
        errorFlag = -1000;
      }
      
      /******************************************/

      *outStream << "\n";
    }
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
  Kokkos::finalize();

  return errorFlag;
}
