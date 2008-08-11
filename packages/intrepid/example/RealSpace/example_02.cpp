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
\brief  Example of the Matrix class.
\author Created by P. Bochev and D. Ridzal
*/
#include "Intrepid_RealSpace.hpp"

using namespace std;
using namespace Intrepid;

int main(int argc, char *argv[]) {
  cout \
  << "===============================================================================\n" \
  << "|                                                                             |\n" \
  << "|                      Example use of the Matrix class                        |\n" \
  << "|                                                                             |\n" \
  << "|   1) Creating Matrix objects and accessing their elements                   |\n" \
  << "|   2) Matrix operations, inverse, transpose, norms and determinant           |\n" \
  << "|   3) Matrix - Point operations                                              |\n" \
  << "|                                                                             |\n" \
  << "|  Questions? Contact  Pavel Bochev (pbboche@sandia.gov) or                   |\n" \
  << "|                      Denis Ridzal (dridzal@sandia.gov).                     |\n" \
  << "|                                                                             |\n" \
  << "|  Intrepid's website: http://trilinos.sandia.gov/packages/intrepid           |\n" \
  << "|  Trilinos website:   http://trilinos.sandia.gov                             |\n" \
  << "|                                                                             |\n" \
  << "===============================================================================\n" \
  << "| EXAMPLE 1: class Matrix in 3D                                               |\n"\
  << "===============================================================================\n\n";
  
  // Matrix objects can be created from square or flat lists of data, arranged by row. If
  double mat3[3][3] = {{1,2,3},{4,5,7},{7,8,10}};
  double arr3[]     = { 1,2,3,  4,5,7,  7,8,10};
  // then Matrix objects created from these lists correspond to the same 3-by-3 matrix:
  //          | 1 2 3 |
  //          | 4 5 6 |
  //          | 7 8 9 |
  //
  cout << "Created Matrix my_3D_Matrix from a 3-by-3 array of numbers:";
  Matrix<double> my_3D_Matrix(&mat3[0][0],3);
  cout << my_3D_Matrix << endl;

  cout << "Created Matrix my_other_3D_Matrix from a 1-by-9 array with the same numbers:";
  Matrix<double> my_other_3D_Matrix(arr3,3);
  cout << my_3D_Matrix << endl;
  
  cout << "Compute the transpose of my_3D_Matrix:";
  cout << my_3D_Matrix.getTranspose() << endl;

  cout << "Transpose my_3D_Matrix (in place -- mutator):";
  my_3D_Matrix.transpose();
  cout << my_3D_Matrix << endl;

  cout << "Computing: my_3D_Matrix * my_other_3D_Matrix:";
  Matrix<double> prod3 = my_3D_Matrix * my_other_3D_Matrix;
  cout << prod3 << endl;

  cout << "Compute the inverse of my_3D_Matrix:";
  cout << my_3D_Matrix.getInverse() << endl;

  cout << "Invert my_3D_Matrix (in place -- mutator):";
  my_3D_Matrix.invert();
  cout << my_3D_Matrix << endl;
  cout << "\t Determinant of the inverted matrix = " << my_3D_Matrix.det() << "\n\n";

  cout << "Invert my_3D_Matrix again (in place -- mutator):";
  my_3D_Matrix.invert();
  cout << my_3D_Matrix << endl;
  cout << "\t Determinant of the inverted matrix = " << my_3D_Matrix.det() << "\n\n";

  // Norms of matrices
  cout << "Computing norms of my_3D_Matrix: " << endl;
  cout << "\t First     matrix norm (max column sum) = " << my_3D_Matrix.norm(NORM_ONE) << endl;
  cout << "\t Infinity  matrix norm (max row    sum) = " << my_3D_Matrix.norm(NORM_INF) << endl;
  cout << "\t Frobenius matrix norm                  = " << my_3D_Matrix.norm(NORM_FRO) << "\n\n";
  
  
  // Extraction of rows and columns as Point objects (they are indexed from 0!)
  cout << "Get first row of my_3D_Matrix as a Point object:\n";
  Point<double> my_1st_Row = my_3D_Matrix.getRow(0);
  cout << my_1st_Row << "\n\n";

  cout << "Get third column of my_3D_Matrix as a Point object:\n";
  Point<double> my_3rd_Column  = my_3D_Matrix.getColumn(2);
  cout << my_3rd_Column << "\n\n";
    
  // Changing data in an existing matrix object
  double my_New_3D_Data[] = {1,1,1,2,2,2,3,3,3};
  cout << "Changing data in an existing matrix (my_3D_Matrix):";
  my_3D_Matrix.setElements(my_New_3D_Data,3);
  cout << my_3D_Matrix << endl;

  cout << "\n" \
  << "===============================================================================\n" \
  << "| EXAMPLE 1a: overloaded operators for Matrix objects                         |\n"\
  << "===============================================================================\n\n";
  
  Matrix<double> first_3D_Matrix(&mat3[0][0],3);
  Matrix<double> second_3D_Matrix(arr3,3);
  Matrix<double> result(3);

  // Subtraction
  cout << "Subtract two matrices holding the same coefficients:";
  result = first_3D_Matrix - second_3D_Matrix;
  cout << result << "\n\n";
  
  // Addition
  cout << "Add two matrices holding the same coefficients:"; 
    result = first_3D_Matrix + second_3D_Matrix;
  cout << result << "\n\n";
  
  // Scalar multiplication
  cout << "Multiply my_3D_Matrix by 10: ";
  cout << 10.0*my_3D_Matrix << endl;

  // Matrix times Point
  cout << "Multiply my_3D_Matrix by a Point that is its first row:\n";
  cout << my_3D_Matrix*my_3D_Matrix.getRow(0) << "\n\n";

  // Point times Matrix
  cout << "Multiply Point that is first column of my_3D_Matrix by the matrix:\n";
  cout << my_3D_Matrix.getColumn(0)*my_3D_Matrix << endl;
  
  
  cout << "\n" \
  << "===============================================================================\n" \
  << "| EXAMPLE 2: class Matrix in 2D                                                |\n"\
  << "===============================================================================\n\n";
  
  // A flat list of matrix data arranged by row
  double mat2[] = {1,2,4,6};

  // Uisng constructor that takes flat list of elements arranged by row and matrix dimension
  cout << "Created Matrix my_2D_Matrix:";
  Matrix<double> my_2D_Matrix(mat2,2);
  cout << my_2D_Matrix << endl;

  cout << "Compute the transpose of my_2D_Matrix:";
  cout << my_2D_Matrix.getTranspose() << endl;

  cout << "Transpose my_2D_Matrix (in place -- mutator):";
  my_2D_Matrix.transpose();
  cout << my_2D_Matrix << endl;

  cout << "Computing: my_2D_Matrix * my_2D_Matrix:";
  Matrix<double> prod2 = my_2D_Matrix * my_2D_Matrix;
  cout << prod2 << endl;

  cout << "Compute the inverse of my_2D_Matrix:";
  cout << my_2D_Matrix.getInverse() << endl;

  cout << "Invert my_2D_Matrix (in place -- mutator):";
  my_2D_Matrix.invert();
  cout << my_2D_Matrix << endl;

  cout << "Invert my_2D_Matrix again (in place -- mutator):";
  my_2D_Matrix.invert();
  cout << my_2D_Matrix << endl;

  cout << "\n" \
  << "===============================================================================\n" \
  << "| EXAMPLE 3: class Matrix in 1D                                                |\n"\
  << "===============================================================================\n\n";
  
  double mat1 = 4;

  cout << "Created Matrix my_1D_Matrix:";
  Matrix<double> my_1D_Matrix(&mat1,1);
  cout << my_1D_Matrix << endl;

  cout << "Compute the transpose of my_1D_Matrix:";
  cout << my_1D_Matrix.getTranspose() << endl;

  cout << "Transpose my_1D_Matrix (in place -- mutator):";
  my_1D_Matrix.transpose();
  cout << my_1D_Matrix << endl;

  cout << "Computing: my_1D_Matrix * my_1D_Matrix:";
  Matrix<double> prod1 = my_1D_Matrix * my_1D_Matrix;
  cout << prod1 << endl;

  cout << "Compute the inverse of my_1D_Matrix:";
  cout << my_1D_Matrix.getInverse() << endl;

  cout << "Invert my_1D_Matrix (in place -- mutator):";
  my_1D_Matrix.invert();
  cout << my_1D_Matrix << endl;

  cout << "Invert my_1D_Matrix again (in place -- mutator):";
  my_1D_Matrix.invert();
  cout << my_1D_Matrix << endl;

  
  
  
  return 0;
}
