// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "NOX_LAPACK_Matrix.H"

using namespace NOX;
using namespace NOX::LAPACK;

Matrix::Matrix() :
  p(0), q(0),             	
  entries()			
{
}

Matrix::Matrix(int m, int n) :
  p(m), q(n),             	// set the size of the matrix
  entries(m*n)			// create an array of m*n entries
{
}

Matrix::Matrix(const Matrix& a, CopyType type) :
  p(a.p), q(a.q),                      // set the size of the square matrix
  entries(a.entries)		      // create a copy of a
{
}

Matrix::~Matrix()
{
}

double& Matrix::operator()(int i, int j)
{
  return entries(i + (p*j));
}

const double& Matrix::operator()(int i, int j) const
{
  return entries(i + (p*j));
}

void Matrix::scale(double value)
{
  entries.scale(value);
  return;
}

ostream& Matrix::leftshift(ostream& stream) const
{
  for (int i = 0; i < p; i++) {
    stream << "[ ";
    for (int j = 0; j < q; j++)
      stream << operator()(i,j) << " ";
    stream << "]" << endl;
  }
  return stream;
}

ostream& operator<<(ostream& stream, const Matrix& m)
{
  return m.leftshift(stream);
}

bool Matrix::print() const
{
  cout << *this << endl;
  return true;
}

int Matrix::numRows() const
{
  return p;
}

int Matrix::numCols() const
{
  return q;
}

