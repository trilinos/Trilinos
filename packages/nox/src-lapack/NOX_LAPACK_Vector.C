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

#include "NOX_Common.H"
#include "NOX_LAPACK_Vector.H"
#include "NOX_BLAS_Wrappers.H"

using namespace NOX;
using namespace NOX::LAPACK;

Vector::Vector(int n) : 
  x(n, 0.0)			// initialize a zero vector of length n
{
}

Vector::Vector(const Vector& source, CopyType type) :
  x(source.x)
{
}

Vector::~Vector()
{
}

Abstract::Vector& Vector::operator=(const vector<double>& source)
{
  x = source;
  return *this;
}

Abstract::Vector& Vector::operator=(const Abstract::Vector& source)
{
  return operator=(dynamic_cast<const Vector&>(source));
}

Abstract::Vector& Vector::operator=(const Vector& source)
{
  x = source.x;
  return *this;
}

Abstract::Vector& Vector::init(double value)
{
  for (int i = 0; i < x.size(); i ++)
    x[i] = value;
  return *this;
}

Abstract::Vector& Vector::abs(const Abstract::Vector& base)
{
  return abs(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& Vector::abs(const Vector& base)
{
  for (int i = 0; i < x.size(); i ++)
    x[i] = fabs(base[i]);
  return *this;
}

Abstract::Vector& Vector::reciprocal(const Abstract::Vector& base)
{
  return reciprocal(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& Vector::reciprocal(const Vector& base)
{
  for (int i = 0; i < x.size(); i ++)
    x[i] = 1.0 / base[i];
  return *this;
}

Abstract::Vector& Vector::scale(double alpha)
{
  for (int i = 0; i < x.size(); i ++)
    x[i] *= alpha;
  return *this;
}

Abstract::Vector& Vector::update(double alpha, const Abstract::Vector& a, 
				 double gamma)
{
  return update(alpha, dynamic_cast<const Vector&>(a), gamma);
}

Abstract::Vector& Vector::update(double alpha, const Vector& a, 
				 double gamma)
{
  for (int i = 0; i < x.size(); i ++)
    x[i] = alpha * a[i] + gamma * x[i];
  return *this;
}

Abstract::Vector& Vector::update(double alpha, const Abstract::Vector& a, 
				 double beta, const Abstract::Vector& b,
				 double gamma)
{
  return update(alpha, dynamic_cast<const Vector&>(a), 
		beta, dynamic_cast<const Vector&>(b), gamma);
}

Abstract::Vector& Vector::update(double alpha, const Vector& a, 
				 double beta, const Vector& b,
				 double gamma)
{
  for (int i = 0; i < x.size(); i ++)
    x[i] = alpha * a[i] + beta * b[i] + gamma * x[i];
  return *this;
}

Abstract::Vector& Vector::scale(const Abstract::Vector& a)
{  
  return scale(dynamic_cast<const Vector&>(a));
}

Abstract::Vector& Vector::scale(const Vector& a)
{  
  for (int i = 0; i < x.size(); i ++)
    x[i] = a[i] * x[i];
  return *this;
}

Abstract::Vector* Vector::clone(CopyType type) const
{
  NOX::LAPACK::Vector* ptr = new Vector(*this, type);
  return ptr;
}

double Vector::norm(Abstract::Vector::NormType type) const
{
  
  if (x.empty())
    return 0.0;

  int i;			// counter
  int n = x.size();		// size of x
  double value;			// final answer

  switch (type) {
  case MaxNorm:
    value = fabs(x[0]);
    for (i = 1; i < n; i ++)
      if (value < fabs(x[i]))
	value = fabs(x[i]);
    break;
  case OneNorm:
    value = DASUM_F77(&n, &x[0], &i_one);
    break;
  case TwoNorm:
  default:
    value = DNRM2_F77(&n, &x[0], &i_one);
   break;
  }

  return value;
}

double Vector::norm(const Abstract::Vector& weights) const
{
  return norm(dynamic_cast<const Vector&>(weights));
}

double Vector::norm(const Vector& weights) const
{
  if (weights.length() != x.size()) {
    cerr << "NOX::LAPACK::Vector::norm - size mismatch for weights vector" << endl;
    throw "NOX::LAPACK Error";
  }

  int n = x.size();		// size of x
  double value = 0;		// final answer

  for (int i = 0; i < n; i ++)
    value += weights[i] * x[i] * x[i];

  return value;
}

double Vector::dot(const Abstract::Vector& y) const
{
  return dot(dynamic_cast<const Vector&>(y));
}

double Vector::dot(const Vector& y) const
{
  if (y.length() != x.size()) {
    cerr << "NOX::LAPACK::Vector::norm - size mismatch for weights vector" << endl;
    throw "NOX::LAPACK Error";
  }

  int n = x.size();		// size of x
  
  return DDOT_F77(&n, &x[0], &i_one, &y[0], &i_one);
}

int Vector::length() const
{
  return x.size();
}

double& Vector::operator[] (int i)
{
  return x[i];
}

const double& Vector::operator[] (int i) const
{
  return x[i];
}

double& Vector::operator() (int i)
{
  return x[i];
}

const double& Vector::operator() (int i) const
{
  return x[i];
}

ostream& Vector::leftshift(ostream& stream) const
{
  stream << "[ ";
  for (int i = 0; i < x.size(); i ++)
    stream << x[i] << " ";
  stream << "]";
  return stream;
}

ostream& operator<<(ostream& stream, const Vector& v)
{
  return v.leftshift(stream);
}

void Vector::print() const
{
  cout << *this << endl;
}
