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
\brief  Unit tests for the Point class.
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
  << "|                              Unit Test (Point)                              |\n" \
  << "|                                                                             |\n" \
  << "|     1) Point creation, math operations                                      |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n"\
  << "| TEST 1: point creation, exception testing                                   |\n"\
  << "===============================================================================\n";

  int errorFlag = 0;
#ifdef HAVE_INTREPID_DEBUG
  int beginThrowNumber = TestForException_getThrowNumber();
  int endThrowNumber = beginThrowNumber + 15;
#endif
  
  // Create arrays of coefficients
  double vec[]  = {1.0, 2.0, 3.0};
  double vec2[] = {-3.0, -1.0, 15.0};

  try{

#ifdef HAVE_INTREPID_DEBUG
    try {
      Point<double> p00(0);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    Point<double> p01(1);
    Point<double> p02(2);
    Point<double> p03(3);

#ifdef HAVE_INTREPID_DEBUG
    try {
      Point<double> p00(4);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

#ifdef HAVE_INTREPID_DEBUG
    try {
      Point<double> p00(vec, 0);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    Point<double> p04(vec, 1);
    Point<double> p05(vec, 2);
    Point<double> p06(vec, 3);

#ifdef HAVE_INTREPID_DEBUG
    try {
      Point<double> p00(vec, 4);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

#ifdef HAVE_INTREPID_DEBUG
    try {
      p06.setCoordinates(vec2, 4);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    p06.setCoordinates(vec2, 3);

#ifdef HAVE_INTREPID_DEBUG
    try {
      p06.distance(p01);
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    p06.distance(p03);

#ifdef HAVE_INTREPID_DEBUG
    try {
      p06[3];
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    p06[2];

#ifdef HAVE_INTREPID_DEBUG
    try {
      p05 = p06;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    p06 = p03;

#ifdef HAVE_INTREPID_DEBUG
    try {
      p03 += p03;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

#ifdef HAVE_INTREPID_DEBUG
    try {
      p03 += p01;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    p06 += p03;

#ifdef HAVE_INTREPID_DEBUG
    try {
      p03 ^= p03;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

#ifdef HAVE_INTREPID_DEBUG
    try {
      p03 ^= p01;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    p06 ^= p03;

#ifdef HAVE_INTREPID_DEBUG
    try {
      p03 -= p03;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

#ifdef HAVE_INTREPID_DEBUG
    try {
      p03 -= p01;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    p06 -= p03;

    p06 ^ p06;
    p06 - p06;
    p06 + p06;

#ifdef HAVE_INTREPID_DEBUG
    try {
      p06 * p01;
    }
    catch (std::logic_error err) {
      *outStream << "Expected Error ----------------------------------------------------------------\n";
      *outStream << err.what() << '\n';
      *outStream << "-------------------------------------------------------------------------------" << "\n\n";
    };
#endif

    p06 * p06;

    p06[0] * p06;
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
      Point<double> p01(dim);
      Point<double> p02(dim);
      Point<double> p03(dim);
      Point<double> p_cross(33.0, -24.0, 5.0);
      p01.setCoordinates(vec, dim);
      p02.setCoordinates(vec2, dim);
      double zero = INTREPID_TOL, mult = 123.123, p01p02dist[] = {4.0, 5.0, 13.0};

      *outStream << "p01 = " << p01 << "\n";
      *outStream << "p02 = " << p02 << "\n";

      p03 = (p01 + p02) - p02 - p01;
      *outStream << "p03 =   (p01 + p02) - p02 - p01 = " << p03 << "\n";
      if (p03.norm(NORM_TWO) > zero || p03.norm(NORM_ONE) > zero || p03.norm(NORM_INF) > zero) {
        *outStream << std::setw(52) << "^^^^----FAILURE!" << "\n";
        errorFlag++;
      }

      p03 = (1.0/mult) * (mult * p01) - p01;
      *outStream << "p03 =   (1.0/mult) * (mult * p01) - p01 = " << p03 << "\n";
      if (p03.norm(NORM_TWO) > zero || p03.norm(NORM_ONE) > zero || p03.norm(NORM_INF) > zero) {
        *outStream << std::setw(60) << "^^^^----FAILURE!" << "\n";
        errorFlag++;
      }

      *outStream << "d   =   distance(p01, p02) = " << p01.distance(p02) << "\n";
      if (std::abs(p01.distance(p02)-p01p02dist[dim-1]) > zero) {
        *outStream << std::setw(45) << "^^^^----FAILURE!" << "\n";
        errorFlag++;
      }

      if (dim==3) {
        p03 = p01^p02;
        *outStream << "p03 =   p01 x p02  = " << p03 << "\n";
        p03 -= p_cross;
        if (p03.norm(NORM_TWO) > zero || p03.norm(NORM_ONE) > zero || p03.norm(NORM_INF) > zero) {
          *outStream << std::setw(39) << "^^^^----FAILURE!";
          *outStream << "   ||difference||_2 = " << p03.norm(NORM_TWO) << "\n";
          errorFlag++;
        }
        p03 = p01^p01;
        *outStream << "p03 =   p01 x p01  = " << p03 << "\n";
        if (p03.norm(NORM_TWO) > zero || p03.norm(NORM_ONE) > zero || p03.norm(NORM_INF) > zero) {
          *outStream << std::setw(39) << "^^^^----FAILURE!" << "\n";
          errorFlag++;
        }
      }
      
      *outStream << "\n";
    }
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n\n";
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
