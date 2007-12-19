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
\brief  Test of the Matrix class.
\author Created by P. Bochev and D. Ridzal
*/

#include "Intrepid_RealSpace.hpp"

using namespace std;
using namespace Intrepid;


int main(int argc, char *argv[]) {

  cout << "\nTEST 1: class Matrix in 3D\n\n";

  double mat3[3][3] = {{1,2,3},{4,5,7},{7,8,10}};

  cout << "Created linear map lmap3\n";
  Matrix<double> lmap3(&mat3[0][0],3);
  cout << lmap3 << endl;

  cout << "Compute the transpose of lmap3\n";
  cout << lmap3.getTranspose() << endl;

  cout << "Transpose lmap3 (in place -- mutator)\n";
  lmap3.Transpose();
  cout << lmap3 << endl;

  cout << "Computing: lmap3 * lmap3\n";
  Matrix<double> prod3 = lmap3 * lmap3;
  cout << prod3 << endl;

  cout << "Compute the inverse of lmap3\n";
  cout << lmap3.getInverse() << endl;

  cout << "Invert lmap3 (in place -- mutator)\n";
  lmap3.Invert();
  cout << lmap3 << endl;

  cout << "Invert lmap3 again (in place -- mutator)\n";
  lmap3.Invert();
  cout << lmap3 << endl;

  cout << "\nEND TEST 1: class Matrix in 3D\n\n";


  cout << "\nTEST 2: class Matrix in 2D\n\n";

  double mat2[] = {1,2,4,6};

  cout << "Created linear map lmap2\n";
  Matrix<double> lmap2(mat2,2);
  cout << lmap2 << endl;

  cout << "Compute the transpose of lmap2\n";
  cout << lmap2.getTranspose() << endl;

  cout << "Transpose lmap2 (in place -- mutator)\n";
  lmap2.Transpose();
  cout << lmap2 << endl;

  cout << "Computing: lmap2 * lmap2\n";
  Matrix<double> prod2 = lmap2 * lmap2;
  cout << prod2 << endl;

  cout << "Compute the inverse of lmap2\n";
  cout << lmap2.getInverse() << endl;

  cout << "Invert lmap2 (in place -- mutator)\n";
  lmap2.Invert();
  cout << lmap2 << endl;

  cout << "Invert lmap2 again (in place -- mutator)\n";
  lmap2.Invert();
  cout << lmap2 << endl;

  cout << "\nEND TEST 2: class Matrix in 2D\n\n";


  cout << "\nTEST 3: class Matrix in 1D\n\n";

  double mat1 = 4;

  cout << "Created linear map lmap1\n";
  Matrix<double> lmap1(&mat1,1);
  cout << lmap1 << endl;

  cout << "Compute the transpose of lmap1\n";
  cout << lmap1.getTranspose() << endl;

  cout << "Transpose lmap1 (in place -- mutator)\n";
  lmap1.Transpose();
  cout << lmap1 << endl;

  cout << "Computing: lmap1 * lmap1\n";
  Matrix<double> prod1 = lmap1 * lmap1;
  cout << prod1 << endl;

  cout << "Compute the inverse of lmap1\n";
  cout << lmap1.getInverse() << endl;

  cout << "Invert lmap1 (in place -- mutator)\n";
  lmap1.Invert();
  cout << lmap1 << endl;

  cout << "Invert lmap1 again (in place -- mutator)\n";
  lmap1.Invert();
  cout << lmap1 << endl;

  cout << "\nEND TEST 3: class Matrix in 1D\n\n";

  return 0;
}
