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

#include "LOCA_Bifurcation_ArcLengthVector.H"  // Class definition

#include "NOX_Common.H"    // for abs()

using namespace LOCA;
using namespace LOCA::Bifurcation; 

ArcLengthVector::ArcLengthVector(const NOX::Abstract::Vector& xVec,
			   double arcParam) :
  xVectorPtr(xVec.clone()),
  arcParam(arcParam) 
{

}

ArcLengthVector::ArcLengthVector(const ArcLengthVector& source) :
  xVectorPtr(source.xVectorPtr->clone(NOX::DeepCopy)),
  arcParam(source.arcParam)  
{
  
}

ArcLengthVector::~ArcLengthVector()
{
  delete xVectorPtr;
}

NOX::Abstract::Vector& ArcLengthVector::init(double gamma)
{
  xVectorPtr->init(gamma);
  arcParam = gamma;
  return *this;
}

NOX::Abstract::Vector& ArcLengthVector::abs(const NOX::Abstract::Vector& y)
{
  const ArcLengthVector* tpY = 0;
  tpY = dynamic_cast<const ArcLengthVector*>(&y);
  
  xVectorPtr->abs(*(tpY->xVectorPtr));
  arcParam = fabs(tpY->getArcParam());
  return *this;
}

NOX::Abstract::Vector& ArcLengthVector::operator=(const NOX::Abstract::Vector& y)
{
 const LOCA::Bifurcation::ArcLengthVector* aly = dynamic_cast<const LOCA::Bifurcation::ArcLengthVector*>(&y);

  return (*this = *aly);
}

LOCA::Bifurcation::ArcLengthVector& ArcLengthVector::operator=(const LOCA::Bifurcation::ArcLengthVector& y)
{
  const ArcLengthVector* tpY = 0;
  tpY = dynamic_cast<const ArcLengthVector*>(&y);
  
  *xVectorPtr = *(tpY->xVectorPtr);
  arcParam = tpY->arcParam;
  return *this;
}

NOX::Abstract::Vector& ArcLengthVector::reciprocal(const NOX::Abstract::Vector& y)
{
  const ArcLengthVector* tpY = 0;
  tpY = dynamic_cast<const ArcLengthVector*>(&y);
  
  xVectorPtr->reciprocal(*(tpY->xVectorPtr));
  arcParam = 1.0 / (tpY->arcParam);
  return *this;
}

NOX::Abstract::Vector& ArcLengthVector::scale(double gamma)
{
  xVectorPtr->scale(gamma);
  arcParam *= gamma;
  return *this;
}

NOX::Abstract::Vector& ArcLengthVector::scale(const NOX::Abstract::Vector& y)
{
  const ArcLengthVector* tpY = 0;
  tpY = dynamic_cast<const ArcLengthVector*>(&y);
  
  xVectorPtr->scale(*(tpY->xVectorPtr));
  arcParam *= tpY->arcParam;
  return *this;
}

NOX::Abstract::Vector& ArcLengthVector::update(double alpha, const NOX::Abstract::Vector& y, double gamma)
{
  const ArcLengthVector* tpY = 0;
  tpY = dynamic_cast<const ArcLengthVector*>(&y);
  
  xVectorPtr->update(alpha, *(tpY->xVectorPtr), gamma);
  arcParam = alpha * tpY->arcParam + gamma;
  return *this;
}

NOX::Abstract::Vector& ArcLengthVector::update(double alpha, 
			 const NOX::Abstract::Vector& y, 
                         double beta, const NOX::Abstract::Vector& b,
                         double gamma)
{
  const ArcLengthVector* tpY = 0;
  tpY = dynamic_cast<const ArcLengthVector*>(&y);

  const ArcLengthVector* tpB = 0;
  tpB = dynamic_cast<const ArcLengthVector*>(&b);  

  xVectorPtr->update(alpha, *(tpY->xVectorPtr), beta, 
		     *(tpB->xVectorPtr), gamma);

  arcParam = alpha * tpY->arcParam + beta * tpB->arcParam + gamma;
  return *this;
}

NOX::Abstract::Vector* ArcLengthVector::clone(NOX::CopyType type) const
{
  ArcLengthVector* newCopy = 0;

  newCopy = new ArcLengthVector(*this);
  
  if (type != NOX::DeepCopy)
    newCopy->init(0.0);
      
  return newCopy;
}

double ArcLengthVector::norm(NormType type) const
{
  double n, nx;
  switch (type) {
  case NOX::Abstract::Vector::MaxNorm:
    nx = xVectorPtr->norm(NOX::Abstract::Vector::MaxNorm);
    n = max(nx, arcParam);
    break;
  case NOX::Abstract::Vector::OneNorm:
    nx = xVectorPtr->norm(NOX::Abstract::Vector::OneNorm);
    n = nx + fabs(arcParam);
    break;
  case NOX::Abstract::Vector::TwoNorm:
  default:
   nx = xVectorPtr->norm(NOX::Abstract::Vector::TwoNorm);
   n = sqrt(nx*nx + arcParam*arcParam);
   break;
  }
  return n;
}

double ArcLengthVector::norm(const NOX::Abstract::Vector& weights) const
{
  cout << "ERROR: LOCA::Bifurcation::ArcLengthVector::norm(weights) - "
       << "Not implemented yet!!!" << endl;
  throw "LOCA Error";
}

double ArcLengthVector::dot(const NOX::Abstract::Vector& y) const
{
  const ArcLengthVector* tpY = 0;
  tpY = dynamic_cast<const ArcLengthVector*>(&y);

  return( xVectorPtr->dot(*tpY->xVectorPtr) + arcParam*arcParam);

}

int ArcLengthVector::length() const
{
  return (xVectorPtr->length() + 1);
}

bool ArcLengthVector::print() const
{
  xVectorPtr->print();
  cout << "ArcLength Parameter = " << arcParam << endl;
  return true;
}

void ArcLengthVector::setVec(const NOX::Abstract::Vector& xVec,
                             double arcPar)
{
  *xVectorPtr = xVec;
  arcParam = arcPar;
}

  //! Returns the solution vector component of extended vector
const NOX::Abstract::Vector& ArcLengthVector::getXVec() const
{
  return *xVectorPtr;
}

  //! Get Arclength parameter
double ArcLengthVector::getArcParam() const
{
  return arcParam;
}

  //! Returns the solution vector component of extended vector
NOX::Abstract::Vector& ArcLengthVector::getXVec()
{
  return *xVectorPtr;
}

  //! Get Arclength parameter
double& ArcLengthVector::getArcParam()
{
  return arcParam;
}
