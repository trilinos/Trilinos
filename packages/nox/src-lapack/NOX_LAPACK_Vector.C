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

#include "NOX_Common.H"
#include "NOX_LAPACK_Vector.H"
#include "NOX_BLAS_Wrappers.H"
#include "NOX_Random.H" // for Random class

NOX::LAPACK::Vector::Vector() :
  n(0),
  x()
{
}

NOX::LAPACK::Vector::Vector(int N) :
  n(N),
  x(N, 0.0)
{
}

NOX::LAPACK::Vector::Vector(int N, double *v) : 
  n(N),
  x(v,v+N)
{
}

NOX::LAPACK::Vector::Vector(const NOX::LAPACK::Vector& source, 
			    NOX::CopyType type) :
  n(source.n),
  x(source.x)
{
}

NOX::LAPACK::Vector::~Vector()
{
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::operator=(
						const vector<double>& source)
{
  x = source;
  return *this;
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::operator=(
					   const NOX::Abstract::Vector& source)
{
  return operator=(dynamic_cast<const NOX::LAPACK::Vector&>(source));
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::operator=(
					   const NOX::LAPACK::Vector& source)
{
  x = source.x;
  return *this;
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::init(double value)
{
  for (int i = 0; i < n; i ++)
    x[i] = value;
  return *this;
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::random(bool useSeed, int seed) 
{
  if (useSeed)
    NOX::Random::setSeed(seed);

  for (int i = 0; i < n; i ++) 
    x[i] = NOX::Random::number();

  return *this;
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::abs(
					     const NOX::Abstract::Vector& base)
{
  return abs(dynamic_cast<const NOX::LAPACK::Vector&>(base));
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::abs(
					     const NOX::LAPACK::Vector& base)
{
  for (int i = 0; i < n; i ++)
    x[i] = fabs(base[i]);
  return *this;
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::reciprocal(
					    const NOX::Abstract::Vector& base)
{
  return reciprocal(dynamic_cast<const NOX::LAPACK::Vector&>(base));
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::reciprocal(
					    const NOX::LAPACK::Vector& base)
{
  for (int i = 0; i < n; i ++)
    x[i] = 1.0 / base[i];
  return *this;
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::scale(double alpha)
{
  for (int i = 0; i <n; i ++)
    x[i] *= alpha;
  return *this;
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::update(
					       double alpha, 
					       const NOX::Abstract::Vector& a, 
					       double gamma)
{
  return update(alpha, dynamic_cast<const NOX::LAPACK::Vector&>(a), gamma);
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::update(
						 double alpha, 
						 const NOX::LAPACK::Vector& a, 
						 double gamma)
{
  for (int i = 0; i < n; i ++)
    x[i] = alpha * a[i] + gamma * x[i];
  return *this;
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::update(
					      double alpha, 
					      const NOX::Abstract::Vector& a, 
					      double beta, 
					      const NOX::Abstract::Vector& b,
					      double gamma)
{
  return update(alpha, dynamic_cast<const NOX::LAPACK::Vector&>(a), 
		beta, dynamic_cast<const NOX::LAPACK::Vector&>(b), gamma);
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::update(
					       double alpha, 
					       const NOX::LAPACK::Vector& a, 
					       double beta, 
					       const NOX::LAPACK::Vector& b,
					       double gamma)
{
  for (int i = 0; i < n; i ++)
    x[i] = alpha * a[i] + beta * b[i] + gamma * x[i];
  return *this;
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::scale(
					      const NOX::Abstract::Vector& a)
{  
  return scale(dynamic_cast<const NOX::LAPACK::Vector&>(a));
}

NOX::Abstract::Vector& NOX::LAPACK::Vector::scale(const NOX::LAPACK::Vector& a)
{  
  for (int i = 0; i < n; i ++)
    x[i] = a[i] * x[i];
  return *this;
}

NOX::Abstract::Vector* NOX::LAPACK::Vector::clone(CopyType type) const
{
  return new NOX::LAPACK::Vector(*this, type);
}

double NOX::LAPACK::Vector::norm(NOX::Abstract::Vector::NormType type) const
{
  
  if (n == 0)
    return 0.0;

  int i;			// counter
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

double NOX::LAPACK::Vector::norm(const NOX::Abstract::Vector& weights) const
{
  return norm(dynamic_cast<const NOX::LAPACK::Vector&>(weights));
}

double NOX::LAPACK::Vector::norm(const NOX::LAPACK::Vector& weights) const
{
  if (weights.length() != n) {
    cerr << "NOX::LAPACK::Vector::norm - size mismatch for weights vector" << endl;
    throw "NOX::LAPACK Error";
  }

  double value = 0;		// final answer

  for (int i = 0; i < n; i ++)
    value += weights[i] * x[i] * x[i];

  value = sqrt(value);

  return value;
}

double NOX::LAPACK::Vector::dot(const NOX::Abstract::Vector& y) const
{
  return dot(dynamic_cast<const NOX::LAPACK::Vector&>(y));
}

double NOX::LAPACK::Vector::dot(const NOX::LAPACK::Vector& y) const
{
  if (y.length() != n) {
    cerr << "NOX::LAPACK::Vector::dot - size mismatch for y vector" << endl;
    throw "NOX::LAPACK Error";
  }

  return DDOT_F77(&n, &x[0], &i_one, &y[0], &i_one);
}

int NOX::LAPACK::Vector::length() const
{
  return n;
}

double& NOX::LAPACK::Vector::operator[] (int i)
{
  return x[i];
}

const double& NOX::LAPACK::Vector::operator[] (int i) const
{
  return x[i];
}

double& NOX::LAPACK::Vector::operator() (int i)
{
  return x[i];
}

const double& NOX::LAPACK::Vector::operator() (int i) const
{
  return x[i];
}

ostream& NOX::LAPACK::Vector::leftshift(ostream& stream) const
{
  stream << "[ ";
  for (int i = 0; i < n; i ++)
    stream << x[i] << " ";
  stream << "]";
  return stream;
}

ostream& operator<<(ostream& stream, const NOX::LAPACK::Vector& v)
{
  return v.leftshift(stream);
}

void NOX::LAPACK::Vector::print() const
{
  cout << *this << endl;
}

