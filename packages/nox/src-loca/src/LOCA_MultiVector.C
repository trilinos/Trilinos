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

#include "NOX_Random.H" //for NOX::random()
#include "LOCA_MultiVector.H"  // Class definition

LOCA::MultiVector::MultiVector(int nvecs, int nscalars) :
  vectorPtrs(nvecs), scalars(nscalars)
{
}

LOCA::MultiVector::MultiVector(const LOCA::MultiVector& source, 
			       NOX::CopyType type) :
  vectorPtrs(source.vectorPtrs.size()), 
  scalars(source.scalars)
{
  for (int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i] = source.vectorPtrs[i]->clone(type);

  if (type != NOX::DeepCopy)
    init(0.0);
}


LOCA::MultiVector::~MultiVector()
{
  for (int i=0; i<vectorPtrs.size(); i++)
    delete vectorPtrs[i];
}

NOX::Abstract::Vector& 
LOCA::MultiVector::operator=(const NOX::Abstract::Vector& y)
{
  return operator=(dynamic_cast<const LOCA::MultiVector&>(y));
}

LOCA::MultiVector& 
LOCA::MultiVector::operator=(const LOCA::MultiVector& y)
{ 
  if (this != &y) {
    for (int i=0; i<vectorPtrs.size(); i++)
      delete vectorPtrs[i];

    vectorPtrs = y.vectorPtrs;
    scalars = y.scalars;

    for (int i=0; i<vectorPtrs.size(); i++)
      vectorPtrs[i] = y.vectorPtrs[i]->clone(NOX::DeepCopy);
  }
  return *this;
}

NOX::Abstract::Vector* 
LOCA::MultiVector::clone(NOX::CopyType type) const
{
  return new LOCA::MultiVector(*this, type);
}

NOX::Abstract::Vector& 
LOCA::MultiVector::init(double gamma)
{
  for (int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->init(gamma);
  for (int i=0; i<scalars.size(); i++)
    scalars[i] = gamma;
  return *this;
}

NOX::Abstract::Vector& 
LOCA::MultiVector::random(bool useSeed, double seed) {
  for (int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->random(useSeed, seed);
  if (useSeed)
    urand.setSeed(seed);
  for (int i=0; i<scalars.size(); i++)
    scalars[i] = urand();
  return *this;
}

NOX::Abstract::Vector& 
LOCA::MultiVector::abs(const NOX::Abstract::Vector& y)
{
  const LOCA::MultiVector& Y = dynamic_cast<const LOCA::MultiVector&>(y);
  
  for (int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->abs(*(Y.vectorPtrs[i]));
  for (int i=0; i<scalars.size(); i++)
    scalars[i] = fabs(Y.scalars[i]);

  return *this;
}

NOX::Abstract::Vector& 
LOCA::MultiVector::reciprocal(const NOX::Abstract::Vector& y)
{
  const LOCA::MultiVector& Y = dynamic_cast<const LOCA::MultiVector&>(y);
  
  for (int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->reciprocal(*(Y.vectorPtrs[i]));
  for (int i=0; i<scalars.size(); i++)
    scalars[i] = 1.0 / Y.scalars[i];
  
  return *this;
}

NOX::Abstract::Vector& 
LOCA::MultiVector::scale(double gamma)
{
  for (int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->scale(gamma);
  for (int i=0; i<scalars.size(); i++)
    scalars[i] *= gamma;

  return *this;
}

NOX::Abstract::Vector& 
LOCA::MultiVector::scale(const NOX::Abstract::Vector& y)
{
  const LOCA::MultiVector& Y = dynamic_cast<const LOCA::MultiVector&>(y);
  
  for (int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->scale(*(Y.vectorPtrs[i]));
  for (int i=0; i<scalars.size(); i++)
    scalars[i] *= Y.scalars[i];

  return *this;
}

NOX::Abstract::Vector& 
LOCA::MultiVector::update(double alpha, const NOX::Abstract::Vector& y, 
			  double gamma)
{
  const LOCA::MultiVector& Y = dynamic_cast<const LOCA::MultiVector&>(y);
  
  for (int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->update(alpha, *(Y.vectorPtrs[i]), gamma);
  for (int i=0; i<scalars.size(); i++)
    scalars[i] = alpha*(Y.scalars[i]) + gamma*scalars[i];

  return *this;
}

NOX::Abstract::Vector& 
LOCA::MultiVector::update(double alpha, 
			  const NOX::Abstract::Vector& y, 
			  double beta, const NOX::Abstract::Vector& b,
			  double gamma)
{
  const LOCA::MultiVector& Y = dynamic_cast<const LOCA::MultiVector&>(y);
  const LOCA::MultiVector& B = dynamic_cast<const LOCA::MultiVector&>(b);
  
  for (int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->update(alpha, *(Y.vectorPtrs[i]), beta, *(B.vectorPtrs[i]),
			  gamma);
  for (int i=0; i<scalars.size(); i++)
    scalars[i] = alpha*(Y.scalars[i]) + beta*(B.scalars[i]) + gamma*scalars[i];

  return *this;
}

double 
LOCA::MultiVector::norm(NormType type) const
{
  double n = 0.0;
  double nv;

  switch (type) {

  case NOX::Abstract::Vector::MaxNorm:
    for (int i=0; i<vectorPtrs.size(); i++) 
      n = max(n, vectorPtrs[i]->norm(type));
    for (int i=0; i<scalars.size(); i++) 
      n = max(n, fabs(scalars[i]));
    break;

  case NOX::Abstract::Vector::OneNorm:
    for (int i=0; i<vectorPtrs.size(); i++) 
      n += vectorPtrs[i]->norm(type);
    for (int i=0; i<scalars.size(); i++) 
      n += fabs(scalars[i]);
    break;

  case NOX::Abstract::Vector::TwoNorm:
  default:
    for (int i=0; i<vectorPtrs.size(); i++) {
      nv = vectorPtrs[i]->norm(type);
      n += nv*nv;
    }
    for (int i=0; i<scalars.size(); i++) 
      n += scalars[i]*scalars[i];
    n = sqrt(n);
   break;

  }

  return n;
}

double 
LOCA::MultiVector::norm(const NOX::Abstract::Vector& weights) const
{
  const LOCA::MultiVector& W = dynamic_cast<const LOCA::MultiVector&>(weights);
  double n = 0.0;
  double nv;

  for (int i=0; i<vectorPtrs.size(); i++) {
    nv = vectorPtrs[i]->norm(*(W.vectorPtrs[i]));
    n += nv*nv;
  }
  for (int i=0; i<scalars.size(); i++) 
    n += W.scalars[i]*scalars[i]*scalars[i];
  n = sqrt(n);

  return n;
}

double 
LOCA::MultiVector::dot(const NOX::Abstract::Vector& y) const
{
  const LOCA::MultiVector& Y = dynamic_cast<const LOCA::MultiVector&>(y);
  double d = 0.0;
  
  for (int i=0; i<vectorPtrs.size(); i++)
    d += vectorPtrs[i]->dot(*(Y.vectorPtrs[i]));
  for (int i=0; i<scalars.size(); i++)
    d += scalars[i]*(Y.scalars[i]);

  return d;
}

int 
LOCA::MultiVector::length() const
{
  int l = 0;
  for (int i=0; i<vectorPtrs.size(); i++)
    l += vectorPtrs[i]->length();
  l += scalars.size();

  return l;
}

void 
LOCA::MultiVector::print() const
{
  for (int i=0; i<vectorPtrs.size(); i++)
    vectorPtrs[i]->print();
  for (int i=0; i<scalars.size(); i++)
    cout << scalars[i] << " ";
  cout << endl;

}

void 
LOCA::MultiVector::setVector(int i, const NOX::Abstract::Vector& v)
{
  if (vectorPtrs[i] != NULL)
    *vectorPtrs[i] = v;
  else
    vectorPtrs[i] = v.clone();
}

void 
LOCA::MultiVector::setScalar(int i, double s)
{
  scalars[i] = s;
}

const NOX::Abstract::Vector& 
LOCA::MultiVector::getVector(int i) const
{
  return *vectorPtrs[i];
}

NOX::Abstract::Vector& 
LOCA::MultiVector::getVector(int i)
{
  return *vectorPtrs[i];
}

double 
LOCA::MultiVector::getScalar(int i) const
{
  return scalars[i];
}

double& 
LOCA::MultiVector::getScalar(int i)
{
  return scalars[i];
}
