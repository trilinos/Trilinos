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

/** \file test_02.cpp
\brief  Unit test for the Matrix class.
\author Created by P. Bochev and D. Ridzal.
*/

#include "Intrepid_RealSpace.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {

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
  << "|                              Unit Test (Matrix)                             |\n" \
  << "|                                                                             |\n" \
  << "|     1) Matrix creation, math operations                                     |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 1: matrix creation, exception testing                                  |\n"\
  << "===============================================================================\n";

  int errorFlag = 0;
#ifdef HAVE_INTREPID_DEBUG
  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 18;
#endif
  
  // Create arrays of coefficients
  double vec[]  = {1.0, 2.0, 3.0, 4.0, 8.0, 7.0, 6.0, 5.0, 9.0};
  double vec2[] = {-1.0, 2.0, -3.0, 4.0, -8.0, 7.0, -6.0, 5.0, -9.0};

  try{

#ifdef HAVE_INTREPID_DEBUG
    try {
      Matrix<double> m00(0);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    Matrix<double> m01(1);
    Matrix<double> m02(2);
    Matrix<double> m03(3);

#ifdef HAVE_INTREPID_DEBUG
    try {
      Matrix<double> m04(4);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

#ifdef HAVE_INTREPID_DEBUG
    try {
      Matrix<double> m05(vec, 0);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    Matrix<double> m06(vec, 1);
    Matrix<double> m07(vec, 2);
    Matrix<double> m08(vec, 3);

#ifdef HAVE_INTREPID_DEBUG
    try {
      Matrix<double> m09(vec, 4);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

#ifdef HAVE_INTREPID_DEBUG
    try {
      m08.setElements(vec2, 2);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    m08.setElements(vec2, 3);

#ifdef HAVE_INTREPID_DEBUG
    try {
      m08.getElement(0,3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    m08.getElement(0,2);

#ifdef HAVE_INTREPID_DEBUG
    try {
      m08.getColumn(3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    m08.getColumn(2);

#ifdef HAVE_INTREPID_DEBUG
    try {
      m08.getRow(3);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    m08.getRow(1);

#ifdef HAVE_INTREPID_DEBUG
    try {
      m02 = m03;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif
    
    m03 = m08;

#ifdef HAVE_INTREPID_DEBUG
    try {
      m03 += m03;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

#ifdef HAVE_INTREPID_DEBUG
    try {
      m02 += m03;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    m03 += m08;

#ifdef HAVE_INTREPID_DEBUG
    try {
      m03 -= m03;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

#ifdef HAVE_INTREPID_DEBUG
    try {
      m02 -= m03;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    m03 -= m08;

#ifdef HAVE_INTREPID_DEBUG
    try {
      Point<double> v01(1);
      v01 = m03 * v01;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    Point<double> v03(3);
    v03 = m03 * v03;

#ifdef HAVE_INTREPID_DEBUG
    try {
      Point<double> v01(1);
      v01 = v01 * m03;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    v03 = v03 * m03;

#ifdef HAVE_INTREPID_DEBUG
    try {
      m03 = m02 * m03;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    m03 = m03 * m03;

    Matrix<double> mZero(3);

#ifdef HAVE_INTREPID_DEBUG
    try {
      m03 = mZero.getInverse();
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    double singvec[] = {1.0, 2.0, 2.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
    Matrix<double> mSing(singvec, 3);

#ifdef HAVE_INTREPID_DEBUG
    try {
      m03 = mSing.getInverse();
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    outStream->precision(4);
    *outStream << m01 << std::endl;
    *outStream << m02 << std::endl;
    *outStream << m03 << std::endl;

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
  << "| TEST 2: correctness of math operations in 1-3D                              |\n"\
  << "===============================================================================\n";

  outStream->precision(20);

  try{
    for (int dim=3; dim>0; dim--) {
      Matrix<double> m01(dim);
      Matrix<double> m02(dim);
      Matrix<double> m03(dim);
      Point<double>  v01(dim);
      Point<double>  v02(dim);
      m01.setElements(vec, dim);
      m02.setElements(vec2, dim);
      v01.setCoordinates(vec2, dim);
      double zero = INTREPID_TOL, mult = 123.123;

      *outStream << "m01 = \n" << m01 << "\n";
      *outStream << "m02 = \n" << m02 << "\n";
      *outStream << "v01 = "   << v01 << "\n\n";

      m03 = (mult*m01 + m02) - m02 - mult*m01;
      *outStream << "m03 =   (mult*m01 + m02) - m02 - mult*m01 = " << "\n" << m03;
      if (m03.norm(NORM_FRO) > zero || m03.norm(NORM_ONE) > zero || m03.norm(NORM_INF) > zero) {
        *outStream << " ^^^^----FAILURE!" << "\n";
        errorFlag++;
      }

      m03 = m01 - (m01.getInverse()).getInverse();
      *outStream << "\nm03 =   m01 - inv(inv(m01)) = " << "\n" << m03;
      if (m03.norm(NORM_FRO) > zero || m03.norm(NORM_ONE) > zero || m03.norm(NORM_INF) > zero) {
        *outStream << " ^^^^----FAILURE!" << "\n";
        errorFlag++;
      }

      m03 = m01 - (m01.getTranspose()).getTranspose();
      *outStream << "\nm03 =   m01 - (m01^T)^T = " << "\n" << m03;
      if (m03.norm(NORM_FRO) > zero || m03.norm(NORM_ONE) > zero || m03.norm(NORM_INF) > zero) {
        *outStream << " ^^^^----FAILURE!" << "\n";
        errorFlag++;
      }

      double res = 1.0 - m01.det() * (m01.getInverse()).det();
      *outStream << "\n1.0 - det(m01) * det(inv(m01)) = " << res << '\n';
      if (std::abs(res) > zero) {
        *outStream << " ^^^^----FAILURE!" << "\n";
        errorFlag++;
      }

      v02 = v01*m01.getTranspose() - m01*v01;
      *outStream << "\nv02 =   (v01^T * m01^T)^T - m01*v01 = " << v02 << "\n";
      if (v02.norm(NORM_TWO) > zero || v02.norm(NORM_ONE) > zero || v02.norm(NORM_INF) > zero) {
        *outStream << " ^^^^----FAILURE!" << "\n";
        errorFlag++;
      }

      m03 = (m01*m02).getInverse() - (m02.getInverse())*(m01.getInverse());
      *outStream << "\nm03 =   inv(m01*m02) - inv(m02)*inv(m01) = " << "\n" << m03;
      if (m03.norm(NORM_FRO) > zero || m03.norm(NORM_ONE) > zero || m03.norm(NORM_INF) > zero) {
        *outStream << " ^^^^----FAILURE!" << "\n";
        errorFlag++;
      }
      
      *outStream << "\n";
    }
  }
  catch (std::logic_error err) {
    *outStream << "UNEXPECTED ERROR !!! ----------------------------------------------------------\n";
    *outStream << err.what() << '\n';
    *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    errorFlag = -1;
  };

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  // reset format state of std::cout
  std::cout.copyfmt(oldFormatState);

  return errorFlag;
}
