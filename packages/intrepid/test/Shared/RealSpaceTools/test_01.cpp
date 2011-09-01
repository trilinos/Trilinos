// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
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


using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

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
  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 50;
#ifndef HAVE_INTREPID_DEBUG_INF_CHECK
  endThrowNumber = beginThrowNumber + 45;
#endif

#endif

  *outStream \
  << "\n"
  << "===============================================================================\n"\
  << "| TEST 1: vector exceptions                                                   |\n"\
  << "===============================================================================\n";

  try{

    FieldContainer<double> a_2_2(2, 2);
    FieldContainer<double> a_10_2(10, 2);
    FieldContainer<double> a_10_3(10, 3);
    FieldContainer<double> a_10_2_2(10, 2, 2);
    FieldContainer<double> a_10_2_3(10, 2, 3);
    FieldContainer<double> a_10_2_2_3(10, 2, 2, 3);

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

    FieldContainer<double> a_10_1_2_3_4(10, 1, 2, 3, 4);
    FieldContainer<double> b_10_1_2_3_4(10, 1, 2, 3, 4);
    FieldContainer<double> a_10(10);
    FieldContainer<double> a_9(9);
    FieldContainer<double> b_10(10);
    FieldContainer<double> a_10_15_4_4(10, 15, 4, 4);
    FieldContainer<double> b_10_15_4_4(10, 15, 4, 4);
    FieldContainer<double> a_10_2_2(10, 2, 2);
    FieldContainer<double> a_10_2_3(10, 2, 3);
    FieldContainer<double> b_10_2_3(10, 2, 3);

    FieldContainer<double> a_1_1(1, 1);
    FieldContainer<double> b_1_1(1, 1);
    FieldContainer<double> a_2_2(2, 2);
    FieldContainer<double> b_2_2(2, 2);
    FieldContainer<double> a_3_3(3, 3);
    FieldContainer<double> b_3_3(3, 3);
    FieldContainer<double> a_2_3(2, 3);
    FieldContainer<double> a_4_4(4, 4);

    FieldContainer<double> a_10_15_1_1(10, 15, 1, 1);
    FieldContainer<double> b_10_15_1_1(10, 15, 1, 1);
    FieldContainer<double> a_10_15_2_2(10, 15, 2, 2);
    FieldContainer<double> b_10_15_2_2(10, 15, 2, 2);
    FieldContainer<double> a_10_15_3_3(10, 15, 3, 3);
    FieldContainer<double> a_10_15_3_2(10, 15, 3, 2);
    FieldContainer<double> b_10_15_3_3(10, 15, 3, 3);
    FieldContainer<double> b_10_14(10, 14);
    FieldContainer<double> b_10_15(10, 15);
    FieldContainer<double> b_10_14_3(10, 14, 3);
    FieldContainer<double> b_10_15_3(10, 15, 3);
    

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
    a_2_2[0] = 1.0;
    a_3_3[0] = 1.0;
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
  if (TestForException_getThrowNumber() != endThrowNumber)
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
      FieldContainer<double> ma_x_x_d_d(i0, i1, dim, dim);
      FieldContainer<double> mb_x_x_d_d(i0, i1, dim, dim);
      FieldContainer<double> mc_x_x_d_d(i0, i1, dim, dim);
      FieldContainer<double> va_x_x_d(i0, i1, dim);
      FieldContainer<double> vb_x_x_d(i0, i1, dim);
      FieldContainer<double> vc_x_x_d(i0, i1, dim);
      FieldContainer<double> vdot_x_x(i0, i1);
      FieldContainer<double> vnorms_x_x(i0, i1);
      FieldContainer<double> vnorms_x(i0);
      double zero = INTREPID_TOL*100.0;

      // fill with random numbers
      for (int i=0; i<ma_x_x_d_d.size(); i++) {
        ma_x_x_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<va_x_x_d.size(); i++) {
        va_x_x_d[i] = Teuchos::ScalarTraits<double>::random();
      }

      *outStream << "\n************ Checking vectorNorm ************\n";
     
      rst::vectorNorm(vnorms_x_x, va_x_x_d, NORM_TWO);
      *outStream << va_x_x_d;
      *outStream << vnorms_x_x;
      if ( std::abs(rst::vectorNorm(&vnorms_x_x[0], vnorms_x_x.size(), NORM_TWO) - 
                    rst::vectorNorm(&va_x_x_d[0], va_x_x_d.size(), NORM_TWO)) > zero) {
        *outStream << "\n\nINCORRECT vectorNorm NORM_TWO\n\n";
        errorFlag = -1000;
      }

      rst::vectorNorm(vnorms_x_x, va_x_x_d, NORM_ONE);
      *outStream << va_x_x_d;
      *outStream << vnorms_x_x;
      if ( std::abs(rst::vectorNorm(&vnorms_x_x[0], vnorms_x_x.size(), NORM_ONE) - 
                    rst::vectorNorm(&va_x_x_d[0], va_x_x_d.size(), NORM_ONE)) > zero) {
        *outStream << "\n\nINCORRECT vectorNorm NORM_ONE\n\n";
        errorFlag = -1000;
      }

      rst::vectorNorm(vnorms_x_x, va_x_x_d, NORM_INF);
      *outStream << va_x_x_d;
      *outStream << vnorms_x_x;
      if ( std::abs(rst::vectorNorm(&vnorms_x_x[0], vnorms_x_x.size(), NORM_INF) - 
                    rst::vectorNorm(&va_x_x_d[0], va_x_x_d.size(), NORM_INF)) > zero) {
        *outStream << "\n\nINCORRECT vectorNorm NORM_INF\n\n";
        errorFlag = -1000;
      }

      /******************************************/


      *outStream << "\n************ Checking inverse, subtract, and vectorNorm ************\n";
     
      rst::inverse(mb_x_x_d_d, ma_x_x_d_d); // B = inv(A)
      rst::inverse(mc_x_x_d_d, mb_x_x_d_d); // C = inv(B) ~= A
      *outStream << ma_x_x_d_d << mb_x_x_d_d << mc_x_x_d_d;

      rst::subtract(&mc_x_x_d_d[0], &ma_x_x_d_d[0], ma_x_x_d_d.size()); // C = C - A ~= 0 

      if (rst::vectorNorm(&mc_x_x_d_d[0], mc_x_x_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT inverse OR subtract OR vectorNorm\n\n";
        errorFlag = -1000;
      }

      /******************************************/


      *outStream << "\n********** Checking determinant **********\n";

      FieldContainer<double> detA_x_x(i0, i1);
      FieldContainer<double> detB_x_x(i0, i1);
      
      rst::det(detA_x_x, ma_x_x_d_d);
      rst::det(detB_x_x, mb_x_x_d_d);
      *outStream << detA_x_x << detB_x_x;
      
      if ( (rst::dot(&detA_x_x[0], &detB_x_x[0], detA_x_x.size()) - (double)(i0*i1)) > zero) {
        *outStream << "\n\nINCORRECT det\n\n" ;
        errorFlag = -1000;
      }

      *outStream << "\n det(A)*det(inv(A)) = " <<
                    rst::det(&ma_x_x_d_d[0], ma_x_x_d_d.dimension(3))*rst::det(&mb_x_x_d_d[0], mb_x_x_d_d.dimension(3))
                 << "\n";

      if ( (rst::det(&ma_x_x_d_d[0], ma_x_x_d_d.dimension(3))*
            rst::det(&mb_x_x_d_d[0], mb_x_x_d_d.dimension(3)) - (double)1) > zero) {
        *outStream << "\n\nINCORRECT det\n\n" ;
        errorFlag = -1000;
      }

      /******************************************/


      *outStream << "\n************ Checking transpose and subtract ************\n";
     
      rst::transpose(mb_x_x_d_d, ma_x_x_d_d); // B = A^T
      rst::transpose(mc_x_x_d_d, mb_x_x_d_d); // C = B^T = A
      *outStream << ma_x_x_d_d << mb_x_x_d_d << mc_x_x_d_d;

      rst::subtract(&mc_x_x_d_d[0], &ma_x_x_d_d[0], ma_x_x_d_d.size()); // C = C - A = 0 

      if (rst::vectorNorm(&mc_x_x_d_d[0], mc_x_x_d_d.size(), NORM_ONE) > zero) {
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
      *outStream << vc_x_x_d;

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
      *outStream << vc_x_x_d;

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

      for (int i=0; i<va_x_x_d.size(); i++) {
        va_x_x_d[i] = 2.0;
      }
      rst::dot(vdot_x_x, va_x_x_d, va_x_x_d); // dot = a'*a
      *outStream << vdot_x_x;

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
      FieldContainer<double> ma_x_d_d(i0, dim, dim);
      FieldContainer<double> mb_x_d_d(i0, dim, dim);
      FieldContainer<double> mc_x_d_d(i0, dim, dim);
      FieldContainer<double> va_x_d(i0, dim);
      FieldContainer<double> vb_x_d(i0, dim);
      FieldContainer<double> vc_x_d(i0, dim);
      FieldContainer<double> vdot_x(i0);
      FieldContainer<double> vnorms_x(i0);
      double zero = INTREPID_TOL*100.0;

      // fill with random numbers
      for (int i=0; i<ma_x_d_d.size(); i++) {
        ma_x_d_d[i] = Teuchos::ScalarTraits<double>::random();
      }
      for (int i=0; i<va_x_d.size(); i++) {
        va_x_d[i] = Teuchos::ScalarTraits<double>::random();
      }

      *outStream << "\n************ Checking vectorNorm ************\n";
     
      rst::vectorNorm(vnorms_x, va_x_d, NORM_TWO);
      *outStream << va_x_d;
      *outStream << vnorms_x;
      if ( std::abs(rst::vectorNorm(&vnorms_x[0], vnorms_x.size(), NORM_TWO) - 
                    rst::vectorNorm(&va_x_d[0], va_x_d.size(), NORM_TWO)) > zero) {
        *outStream << "\n\nINCORRECT vectorNorm NORM_TWO\n\n";
        errorFlag = -1000;
      }

      rst::vectorNorm(vnorms_x, va_x_d, NORM_ONE);
      *outStream << va_x_d;
      *outStream << vnorms_x;
      if ( std::abs(rst::vectorNorm(&vnorms_x[0], vnorms_x.size(), NORM_ONE) - 
                    rst::vectorNorm(&va_x_d[0], va_x_d.size(), NORM_ONE)) > zero) {
        *outStream << "\n\nINCORRECT vectorNorm NORM_ONE\n\n";
        errorFlag = -1000;
      }

      rst::vectorNorm(vnorms_x, va_x_d, NORM_INF);
      *outStream << va_x_d;
      *outStream << vnorms_x;
      if ( std::abs(rst::vectorNorm(&vnorms_x[0], vnorms_x.size(), NORM_INF) - 
                    rst::vectorNorm(&va_x_d[0], va_x_d.size(), NORM_INF)) > zero) {
        *outStream << "\n\nINCORRECT vectorNorm NORM_INF\n\n";
        errorFlag = -1000;
      }

      /******************************************/


      *outStream << "\n************ Checking inverse, subtract, and vectorNorm ************\n";
     
      rst::inverse(mb_x_d_d, ma_x_d_d); // B = inv(A)
      rst::inverse(mc_x_d_d, mb_x_d_d); // C = inv(B) ~= A
      *outStream << ma_x_d_d << mb_x_d_d << mc_x_d_d;

      rst::subtract(&mc_x_d_d[0], &ma_x_d_d[0], ma_x_d_d.size()); // C = C - A ~= 0 

      if (rst::vectorNorm(&mc_x_d_d[0], mc_x_d_d.size(), NORM_ONE) > zero) {
        *outStream << "\n\nINCORRECT inverse OR subtract OR vectorNorm\n\n";
        errorFlag = -1000;
      }

      /******************************************/


      *outStream << "\n********** Checking determinant **********\n";

      FieldContainer<double> detA_x(i0);
      FieldContainer<double> detB_x(i0);
      
      rst::det(detA_x, ma_x_d_d);
      rst::det(detB_x, mb_x_d_d);
      *outStream << detA_x << detB_x;
      
      if ( (rst::dot(&detA_x[0], &detB_x[0], detA_x.size()) - (double)i0) > zero) {
        *outStream << "\n\nINCORRECT det\n\n" ;
        errorFlag = -1000;
      }

      *outStream << "\n det(A)*det(inv(A)) = " <<
                    rst::det(&ma_x_d_d[0], ma_x_d_d.dimension(2))*rst::det(&mb_x_d_d[0], mb_x_d_d.dimension(2))
                 << "\n";

      if ( (rst::det(&ma_x_d_d[0], ma_x_d_d.dimension(2))*
            rst::det(&mb_x_d_d[0], mb_x_d_d.dimension(2)) - (double)1) > zero) {
        *outStream << "\n\nINCORRECT det\n\n" ;
        errorFlag = -1000;
      }

      /******************************************/


      *outStream << "\n************ Checking transpose and subtract ************\n";
     
      rst::transpose(mb_x_d_d, ma_x_d_d); // B = A^T
      rst::transpose(mc_x_d_d, mb_x_d_d); // C = B^T = A
      *outStream << ma_x_d_d << mb_x_d_d << mc_x_d_d;

      rst::subtract(&mc_x_d_d[0], &ma_x_d_d[0], ma_x_d_d.size()); // C = C - A = 0 

      if (rst::vectorNorm(&mc_x_d_d[0], mc_x_d_d.size(), NORM_ONE) > zero) {
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
      *outStream << vc_x_d;

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
      *outStream << vc_x_d;

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

      for (int i=0; i<va_x_d.size(); i++) {
        va_x_d[i] = 2.0;
      }
      rst::dot(vdot_x, va_x_d, va_x_d); // dot = a'*a
      *outStream << vdot_x;

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

  return errorFlag;
}
