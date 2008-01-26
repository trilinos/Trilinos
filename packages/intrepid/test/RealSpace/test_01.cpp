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

/** \file
\brief  Illustrates use of the Point class.
\author Created by P. Bochev and D. Ridzal
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
  << "| TEST 1: point creation, math operations                                     |\n"\
  << "===============================================================================\n";

  int errorFlag = 0;
  
  // Create arrays of coefficients
  double vec[]  = {1.0, 2.0, 3.0};
  double vec2[] = {-4.0, -1.0, 10.0};

  try{
    Point<double> p01(1);
    Point<double> p02(2);
    Point<double> p03(3);
    Point<double> p04(vec, 1);
    Point<double> p05(vec, 2);
    Point<double> p06(vec, 3);
    p06.setCoordinates(vec2, 3);
    p06.distance(p03);
    p06[2];
    p06 = p03;
    p06 += p03;
    p06 ^= p03;
    p06 -= p03;
    p06 ^ p06;
    p06 - p06;
    p06 + p06;
    p06 * p06;
    p06[0] * p06;
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n\n";
    errorFlag = -1;
  };
  
  
  // Using constructor that takes pointer, dimension and uses default for the frame kind.
  Point<double> p1(vec, 3);
  cout << "\tCreated point:  p1 = " << p1 << endl;
  Point<double> p2(vec2, 3);
  cout << "\tCreated point:  p2 = " << p2 << "\n\n";
  
  // Using overloaded [] to access Point coordinates
  cout << "Using overloaded [] to access Point coordinates:\n ";
  cout << "\t x-coordinate: p2[0] = " << p2[0] <<"\n";
  cout << "\t y-coordinate: p2[1] = " << p2[1] <<"\n";
  cout << "\t z-coordinate: p2[2] = " << p2[2] <<"\n\n";
  
  cout << "Examples of operations between points in 3D:\n";
  cout << "\t p1 + p2           = " << p1+p2 << endl;
  cout << "\t p1 - p2           = " << p1-p2 << endl;
  cout << "\t p1 . p2           = " << p1*p2 << endl;
  cout << "\t p1 x p2           = " << (p1^p2) << endl;
  cout << "\t 5.0 * p1          = " << 5.0*p1 << endl;
  cout << "\t 1/(p1.p2)*(p1xp2) = " << (1.0 / (p1*p2)) * (p1 ^ p2) << endl;
  cout << "Computing distance between points:"<< endl;
  cout << "\t |p1-p2| = " << p1.distance(p2) << endl;
  cout << "\t |p2-p1| = " << p2.distance(p1) << endl;
  cout << "\t |p1-p1| = " << p1.distance(p1) << "\n";
  cout << "Computing norms of points:"<< endl;
  cout << "\t || p1 ||_2  = " << p1.norm(NORM_TWO) \
    << "\t || p2 ||_2  = " << p2.norm(NORM_TWO)<< endl;
  cout << "\t || p1 ||_1  = " << p1.norm(NORM_ONE) \
    << "\t || p2 ||_1  = " << p2.norm(NORM_ONE)<< endl;
  cout << "\t || p1 ||_I  = " << p1.norm(NORM_INF) \
    << "\t || p2 ||_I  = " << p2.norm(NORM_INF)<< endl;
  //
  cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 2: class Point in 2D                                                |\n"\
    << "===============================================================================\n\n";
  
  // Using constructor that takes pointer, dimension and uses default for the frame kind.
  Point<double> p4(vec, 2);
  cout << "\tCreated point:  p4 = " << p4 << endl;
  Point<double> p5(vec2, 2);
  cout << "\tCreated point:  p5 = " << p5 << endl;
  
  cout << "Examples of operations between points in 2D:\n";
  cout << "\t p4 + p5           = " << p4+p5 << endl;
  cout << "\t p4 - p5           = " << p4-p5 << endl;
  cout << "\t p4 . p5           = " << p4*p5 << endl;
  cout << "\t 5.0 * p4          = " << 5.0*p4 << endl;
  cout << "\t 1/(p4.p5)*(p4+p5) = " << (1.0 / (p4*p5)) * (p4 + p5) << endl;
  cout << "Computing distance between points:"<< endl;
  cout << "\t |p4-p5| = " << p4.distance(p5) << endl;
  cout << "\t |p5-p4| = " << p5.distance(p4) << endl;
  cout << "\t |p4-p4| = " << p4.distance(p4) << endl;
  cout << "Computing norms of 2D points:"<< endl;
  cout << "\t || p4 ||_2  = " << p4.norm(NORM_TWO) \
    << "\t || p5 ||_2  = " << p5.norm(NORM_TWO)<< endl;
  cout << "\t || p4 ||_1  = " << p4.norm(NORM_ONE) \
    << "\t || p5 ||_1  = " << p5.norm(NORM_ONE)<< endl;
  cout << "\t || p4 ||_I  = " << p4.norm(NORM_INF) \
    << "\t || p5 ||_I  = " << p5.norm(NORM_INF)<< endl;
  //
  cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 3: class Point in 2D                                                |\n"\
    << "===============================================================================\n\n";
  
  // Using constructor that takes pointer, dimension and uses default for the frame kind.
  Point<double> p7(vec, 1);
  cout << "\tCreated point:  p7 = " << p7 << endl;
  Point<double> p8(vec2, 1);
  cout << "\tCreated point:  p8 = " << p8 << endl;
  cout << "Examples of operations between points in 2D:\n";
  cout << "\t p7 + p8           = " << p7+p8 << endl;
  cout << "\t p7 - p8           = " << p7-p8 << endl;
  cout << "\t p7 . p8           = " << p7*p8 << endl;
  cout << "\t 5.0 * p7          = " << 5.0*p7 << endl;
  cout << "\t 1/(p7.p8)*(p7+p8) = " << (1.0 / (p7*p8)) * (p7 + p8) << endl;
  cout << "Computing distance between points in 1D:"<< endl;
  cout << "\t |p7-p8| = " << p7.distance(p8) << endl;
  cout << "\t |p8-p7| = " << p8.distance(p7) << endl;
  cout << "\t |p7-p7| = " << p7.distance(p7) << endl;
  cout << "Computing norms of 1D points:"<< endl;
  cout << "\t || p7 ||_2  = " << p7.norm(NORM_TWO) \
    << "\t || p8 ||_2  = " << p8.norm(NORM_TWO)<< endl;
  cout << "\t || p7 ||_1  = " << p7.norm(NORM_ONE) \
    << "\t || p8 ||_1  = " << p8.norm(NORM_ONE)<< endl;
  cout << "\t || p7 ||_I  = " << p7.norm(NORM_INF) \
    << "\t || p8 ||_I  = " << p8.norm(NORM_INF)<< endl;  
  //
  cout << "\n" \
    << "===============================================================================\n"\
    << "| EXAMPLE 4: other operations on Point objects                                |\n"\
    << "===============================================================================\n\n";
  //
  cout << "Changing the frame kind of a point: " << endl;
  cout << "\t Original Point p1 is : " << p1 << endl;
  p1.setFrameKind(FRAME_REFERENCE);
  cout << "\t      Now Point p1 is : " << p1 << endl;
  
  return 0;
}
