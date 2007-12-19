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
\brief  Test of the Point class.
\author Created by P. Bochev and D. Ridzal
*/
#include "Intrepid_RealSpace.hpp"

using namespace std;
using namespace Intrepid;


int main(int argc, char *argv[]) {

  Teuchos::SerialDenseMatrix<int,double> My_Matrix( 3, 4 );
  My_Matrix.random();             // random numbers
  std::cout<< My_Matrix << std::endl;

  std::cout << "Machine epsilon from Teuchos =" << DBL_EPSILON <<std::endl;
  std::ostream.precision(10);
  std::cout << "THRESHOLD =" << THRESHOLD <<std::endl;
  std::cout << "MINUS_ONE =" << MINUS_ONE <<std::endl;
  std::cout << "PLUS_ONE  =" << PLUS_ONE <<std::endl;

  
  cout << "\nTEST 1: class Point in 3D\n\n";

  double vec[]  = {1.0, 2.0, 3.0};
  double vec2[] = {-4.0, -1.0, 10.0};

  Point<double> v1(vec, 3);
  cout << "Created vector v1:\n" << v1 << endl;
  Point<double> v2(vec2, 3);
  cout << "Created vector v2:\n" << v2 << endl;
  cout << "Computing: v1 + v2\n";
  cout << v1+v2 << endl;
  cout << "Computing: v1 - v2\n";
  cout << v1-v2 << endl;
  cout << "Computing: v1 dot v2\n    ";
  cout << v1*v2 << endl;
  cout << "Computing: (v1 cross v2)\n";
  cout << (v1^v2) << endl;
  cout << "Computing: 5.0 * v1\n";
  cout << 5.0*v1 << endl;
  cout << "Computing: v3 = (1.0 / (v1 dot v2)) * v1 cross v2\n";
  Point<double> v3 = (1.0 / (v1*v2)) * v1 ^ v2;
  cout << v3 << "\n";
  cout << "Computing distance:"<< endl;
  cout << "\t |v1-v2| = " << v1.distance(v2) << endl;
  cout << "\t |v2-v1| = " << v2.distance(v1) << endl;
  cout << "\t |v1-v1| = " << v1.distance(v1) << endl;
  cout << "Computing norm:"<< endl;
  cout << "\t || v1 ||_2   = " << v1.norm(NORM_TWO) << "|| v1 ||_INF = " << v1.norm(NORM_INF);

  cout << "\nEND TEST 1: class Point in 3D\n\n";


  cout << "\nTEST 2: class Point in 2D\n\n";

  Point<double> v4(vec, 2);
  cout << "Created vector v4:\n" << v4 << endl;
  Point<double> v5(vec2, 2);
  cout << "Created vector v5:\n" << v5 << endl;
  cout << "Computing: v4 + v5\n";
  cout << v4+v5 << endl;
  cout << "Computing: v4 - v5\n";
  cout << v4-v5 << endl;
  cout << "Computing: v4 dot v5\n    ";
  cout << v4*v5 << endl;
  cout << "Computing: 5.0 * v4\n";
  cout << 5.0*v4 << endl;
  cout << "Computing: v6 = (1.0 / (v4 dot v5)) * (v4 + v5)\n";
  Point<double> v6 = (1.0 / (v4*v5)) * (v4 + v5);
  cout << v6 << "\n\n";
  cout << "Computing distance:"<< endl;
  cout << "\t |v4-v5| = " << v4.distance(v5) << endl;
  cout << "\t |v5-v4| = " << v5.distance(v4) << endl;
  cout << "\t |v4-v4| = " << v4.distance(v4) << endl;
  
  cout << "\nEND TEST 2: class Point in 2D\n\n";


  cout << "\nTEST 3: class Point in 1D\n\n";

  Point<double> v7(vec, 1);
  cout << "Created vector v7:\n" << v7 << endl;
  Point<double> v8(vec2, 1);
  cout << "Created vector v8:\n" << v8 << endl;
  cout << "Computing: v7 + v8\n";
  cout << v7+v8 << endl;
  cout << "Computing: v7 - v8\n";
  cout << v7-v8 << endl;
  cout << "Computing: v7 dot v8\n    ";
  cout << v7*v8 << endl;
  cout << "Computing: 5.0 * v7\n";
  cout << 5.0*v7 << endl;
  cout << "Computing: v9 = (1.0 / (v7 dot v8)) * (v7 + v8)\n";
  Point<double> v9 = (1.0 / (v7*v8)) * (v7 + v8);
  cout << v9 << "\n";
  cout << "Computing distance:"<< endl;
  cout << "\t |v7-v8| = " << v7.distance(v8) << endl;
  cout << "\t |v8-v7| = " << v8.distance(v7) << endl;
  cout << "\t |v7-v7| = " << v7.distance(v7) << endl;
  
  cout << "\nEND TEST 3: class Point in 1D\n\n";



  return 0;
}
