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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#include "NOX_BLAS_Matrix.H"

using namespace NOX;
using namespace NOX::BLAS;

Matrix::Matrix(int m, int n) :
  p(m), q(n),			// set the size of the square matrix
  entries(m*n)			// create an array of n*n entries
{
}

Matrix::Matrix(const Matrix& a, CopyType type) :
  p(a.p), q(a.q),		// set the size of the square matrix
  entries(a.entries)		// create a copy of a
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
