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

#include "NOX_Random.H" //for NOX::Random
#include "LOCA_ExtendedVector.H"  // Class definition

LOCA::ExtendedVector::ExtendedVector(int nvecs, int nscalars) :
  vectorPtrs(nvecs), scalars(nscalars)
{
}

LOCA::ExtendedVector::ExtendedVector(const LOCA::ExtendedVector& source, 
			       NOX::CopyType type) :
  vectorPtrs(source.vectorPtrs.size()), 
  scalars(source.scalars)
{
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i] = source.vectorPtrs[i]->clone(type);

  if (type != NOX::DeepCopy)
    init(0.0);
}


LOCA::ExtendedVector::~ExtendedVector()
{
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    delete vectorPtrs[i];
}

NOX::Abstract::Vector& 
LOCA::ExtendedVector::operator=(const NOX::Abstract::Vector& y)
{
  return operator=(dynamic_cast<const LOCA::ExtendedVector&>(y));
}

LOCA::ExtendedVector& 
LOCA::ExtendedVector::operator=(const LOCA::ExtendedVector& y)
{ 
  if (this != &y) {
    for (unsigned int i=0; i<vectorPtrs.size(); i++)
      delete vectorPtrs[i];

    vectorPtrs = y.vectorPtrs;
    scalars = y.scalars;

    for (unsigned int i=0; i<vectorPtrs.size(); i++)
      vectorPtrs[i] = y.vectorPtrs[i]->clone(NOX::DeepCopy);
  }
  return *this;
}

NOX::Abstract::Vector* 
LOCA::ExtendedVector::clone(NOX::CopyType type) const
{
  return new LOCA::ExtendedVector(*this, type);
}

NOX::Abstract::Vector& 
LOCA::ExtendedVector::init(double gamma)
{
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->init(gamma);
  for (unsigned int i=0; i<scalars.size(); i++)
    scalars[i] = gamma;
  return *this;
}

NOX::Abstract::Vector& 
LOCA::ExtendedVector::random(bool useSeed, int seed) {
  if (useSeed)
    NOX::Random::setSeed(seed);
  for (unsigned int i=0; i<scalars.size(); i++)
    scalars[i] = NOX::Random::number();
  if (vectorPtrs.size() > 0)
    vectorPtrs[0]->random(useSeed, seed);
  for (unsigned int i=1; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->random();
  return *this;
}

NOX::Abstract::Vector& 
LOCA::ExtendedVector::abs(const NOX::Abstract::Vector& y)
{
  const LOCA::ExtendedVector& Y = dynamic_cast<const LOCA::ExtendedVector&>(y);
  
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->abs(*(Y.vectorPtrs[i]));
  for (unsigned int i=0; i<scalars.size(); i++)
    scalars[i] = fabs(Y.scalars[i]);

  return *this;
}

NOX::Abstract::Vector& 
LOCA::ExtendedVector::reciprocal(const NOX::Abstract::Vector& y)
{
  const LOCA::ExtendedVector& Y = dynamic_cast<const LOCA::ExtendedVector&>(y);
  
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->reciprocal(*(Y.vectorPtrs[i]));
  for (unsigned int i=0; i<scalars.size(); i++)
    scalars[i] = 1.0 / Y.scalars[i];
  
  return *this;
}

NOX::Abstract::Vector& 
LOCA::ExtendedVector::scale(double gamma)
{
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->scale(gamma);
  for (unsigned int i=0; i<scalars.size(); i++)
    scalars[i] *= gamma;

  return *this;
}

NOX::Abstract::Vector& 
LOCA::ExtendedVector::scale(const NOX::Abstract::Vector& y)
{
  const LOCA::ExtendedVector& Y = dynamic_cast<const LOCA::ExtendedVector&>(y);
  
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->scale(*(Y.vectorPtrs[i]));
  for (unsigned int i=0; i<scalars.size(); i++)
    scalars[i] *= Y.scalars[i];

  return *this;
}

NOX::Abstract::Vector& 
LOCA::ExtendedVector::update(double alpha, const NOX::Abstract::Vector& y, 
			  double gamma)
{
  const LOCA::ExtendedVector& Y = dynamic_cast<const LOCA::ExtendedVector&>(y);
  
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->update(alpha, *(Y.vectorPtrs[i]), gamma);
  for (unsigned int i=0; i<scalars.size(); i++)
    scalars[i] = alpha*(Y.scalars[i]) + gamma*scalars[i];

  return *this;
}

NOX::Abstract::Vector& 
LOCA::ExtendedVector::update(double alpha, 
			  const NOX::Abstract::Vector& y, 
			  double beta, const NOX::Abstract::Vector& b,
			  double gamma)
{
  const LOCA::ExtendedVector& Y = dynamic_cast<const LOCA::ExtendedVector&>(y);
  const LOCA::ExtendedVector& B = dynamic_cast<const LOCA::ExtendedVector&>(b);
  
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->update(alpha, *(Y.vectorPtrs[i]), beta, *(B.vectorPtrs[i]),
			  gamma);
  for (unsigned int i=0; i<scalars.size(); i++)
    scalars[i] = alpha*(Y.scalars[i]) + beta*(B.scalars[i]) + gamma*scalars[i];

  return *this;
}

double 
LOCA::ExtendedVector::norm(NormType type) const
{
  double n = 0.0;
  double nv;

  switch (type) {

  case NOX::Abstract::Vector::MaxNorm:
    for (unsigned int i=0; i<vectorPtrs.size(); i++) 
      n = max(n, vectorPtrs[i]->norm(type));
    for (unsigned int i=0; i<scalars.size(); i++) 
      n = max(n, fabs(scalars[i]));
    break;

  case NOX::Abstract::Vector::OneNorm:
    for (unsigned int i=0; i<vectorPtrs.size(); i++) 
      n += vectorPtrs[i]->norm(type);
    for (unsigned int i=0; i<scalars.size(); i++) 
      n += fabs(scalars[i]);
    break;

  case NOX::Abstract::Vector::TwoNorm:
  default:
    for (unsigned int i=0; i<vectorPtrs.size(); i++) {
      nv = vectorPtrs[i]->norm(type);
      n += nv*nv;
    }
    for (unsigned int i=0; i<scalars.size(); i++) 
      n += scalars[i]*scalars[i];
    n = sqrt(n);
   break;

  }

  return n;
}

double 
LOCA::ExtendedVector::norm(const NOX::Abstract::Vector& weights) const
{
  const LOCA::ExtendedVector& W = dynamic_cast<const LOCA::ExtendedVector&>(weights);
  double n = 0.0;
  double nv;

  for (unsigned int i=0; i<vectorPtrs.size(); i++) {
    nv = vectorPtrs[i]->norm(*(W.vectorPtrs[i]));
    n += nv*nv;
  }
  for (unsigned int i=0; i<scalars.size(); i++) 
    n += W.scalars[i]*scalars[i]*scalars[i];
  n = sqrt(n);

  return n;
}

double 
LOCA::ExtendedVector::dot(const NOX::Abstract::Vector& y) const
{
  const LOCA::ExtendedVector& Y = dynamic_cast<const LOCA::ExtendedVector&>(y);
  double d = 0.0;
  
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    d += vectorPtrs[i]->dot(*(Y.vectorPtrs[i]));
  for (unsigned int i=0; i<scalars.size(); i++)
    d += scalars[i]*(Y.scalars[i]);

  return d;
}

int 
LOCA::ExtendedVector::length() const
{
  int l = 0;
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    l += vectorPtrs[i]->length();
  l += scalars.size();

  return l;
}

void 
LOCA::ExtendedVector::print() const
{
  for (unsigned int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->print();
  for (unsigned int i=0; i<scalars.size(); i++)
    cout << scalars[i] << " ";
  cout << endl;

}

void 
LOCA::ExtendedVector::setVector(int i, const NOX::Abstract::Vector& v)
{
  if (vectorPtrs[i] != NULL)
    *vectorPtrs[i] = v;
  else
    vectorPtrs[i] = v.clone();
}

void 
LOCA::ExtendedVector::setScalar(int i, double s)
{
  scalars[i] = s;
}

const NOX::Abstract::Vector& 
LOCA::ExtendedVector::getVector(int i) const
{
  return *vectorPtrs[i];
}

NOX::Abstract::Vector& 
LOCA::ExtendedVector::getVector(int i)
{
  return *vectorPtrs[i];
}

double 
LOCA::ExtendedVector::getScalar(int i) const
{
  return scalars[i];
}

double& 
LOCA::ExtendedVector::getScalar(int i)
{
  return scalars[i];
}
