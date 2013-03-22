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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "NOX_Common.H"
#include "NOX_LAPACK_Vector.H"
#include "Teuchos_BLAS_wrappers.hpp"  // for blas access
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

NOX::Abstract::Vector& 
NOX::LAPACK::Vector::operator=(const std::vector<double>& source)
{
  x = source;
  return *this;
}

NOX::Abstract::Vector& 
NOX::LAPACK::Vector::operator=(const NOX::Abstract::Vector& source)
{
  return operator=(dynamic_cast<const NOX::LAPACK::Vector&>(source));
}

NOX::Abstract::Vector&
NOX::LAPACK::Vector::operator=(const NOX::LAPACK::Vector& source)
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

Teuchos::RCP<NOX::Abstract::Vector> NOX::LAPACK::Vector::
clone(CopyType type) const
{
  Teuchos::RCP<NOX::Abstract::Vector> tmp;
  tmp = Teuchos::rcp(new NOX::LAPACK::Vector(*this, type));
  return tmp;
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
    std::cerr << "NOX::LAPACK::Vector::norm - size mismatch for weights vector" << std::endl;
    throw "NOX::LAPACK Error";
  }

  double value = 0;		// final answer

  for (int i = 0; i < n; i ++)
    value += weights[i] * x[i] * x[i];

  value = sqrt(value);

  return value;
}

double NOX::LAPACK::Vector::innerProduct(const NOX::Abstract::Vector& y) const
{
  return innerProduct(dynamic_cast<const NOX::LAPACK::Vector&>(y));
}

double NOX::LAPACK::Vector::innerProduct(const NOX::LAPACK::Vector& y) const
{
  if (y.length() != n) {
    std::cerr << "NOX::LAPACK::Vector::innerProduct - size mismatch for y vector" 
	 << std::endl;
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

std::ostream& NOX::LAPACK::Vector::leftshift(std::ostream& stream) const
{
  stream << "[ ";
  for (int i = 0; i < n; i ++)
    stream << x[i] << " ";
  stream << "]";
  return stream;
}

std::ostream& NOX::LAPACK::operator<<(std::ostream& stream, const NOX::LAPACK::Vector& v)
{
  return v.leftshift(stream);
}

void NOX::LAPACK::Vector::print(std::ostream& stream) const
{
  stream << *this << std::endl;
}

