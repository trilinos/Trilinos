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

#include "LOCA_Bifurcation_TPBordVector.H"  // Class definition

#include "NOX_Common.H"    // for abs()

using namespace LOCA;
using namespace LOCA::Bifurcation; 

TPBordVector::TPBordVector(const NOX::Abstract::Vector& xVec,
			   const NOX::Abstract::Vector& nullVec,
			   double bifParam) :
  xVectorPtr(xVec.clone()),
  nullVectorPtr(nullVec.clone()),
  bifParam(bifParam) 
{

}

TPBordVector::TPBordVector(const TPBordVector& source) :
  xVectorPtr(source.xVectorPtr->clone(NOX::DeepCopy)),
  nullVectorPtr(source.nullVectorPtr->clone(NOX::DeepCopy)),
  bifParam(source.bifParam)  
{
  
}


TPBordVector::~TPBordVector()
{
  delete xVectorPtr;
  delete nullVectorPtr;
}


NOX::Abstract::Vector& TPBordVector::init(double gamma)
{
  xVectorPtr->init(gamma);
  nullVectorPtr->init(gamma);
  bifParam = gamma;
  return *this;
}

NOX::Abstract::Vector& TPBordVector::abs(const NOX::Abstract::Vector& y)
{
  const TPBordVector* tpY = 0;
  tpY = dynamic_cast<const TPBordVector*>(&y);
  
  xVectorPtr->abs(*(tpY->xVectorPtr));
  nullVectorPtr->abs(*(tpY->nullVectorPtr));
  bifParam = fabs(tpY->getBifParam());
  return *this;
}

NOX::Abstract::Vector& TPBordVector::operator=(const NOX::Abstract::Vector& y)
{
  return operator=(dynamic_cast<const TPBordVector&>(y));
}

TPBordVector& TPBordVector::operator=(const TPBordVector& y)
{ 
  *xVectorPtr = *(y.xVectorPtr);
  *nullVectorPtr = *(y.nullVectorPtr);
  bifParam = y.bifParam;
  return *this;
}

NOX::Abstract::Vector& TPBordVector::reciprocal(const NOX::Abstract::Vector& y)
{
  const TPBordVector* tpY = 0;
  tpY = dynamic_cast<const TPBordVector*>(&y);
  
  xVectorPtr->reciprocal(*(tpY->xVectorPtr));
  nullVectorPtr->reciprocal(*(tpY->nullVectorPtr));
  bifParam = 1.0 / (tpY->bifParam);
  return *this;
}

NOX::Abstract::Vector& TPBordVector::scale(double gamma)
{
  xVectorPtr->scale(gamma);
  nullVectorPtr->scale(gamma);
  bifParam *= gamma;
  return *this;
}

NOX::Abstract::Vector& TPBordVector::scale(const NOX::Abstract::Vector& y)
{
  const TPBordVector* tpY = 0;
  tpY = dynamic_cast<const TPBordVector*>(&y);
  
  xVectorPtr->scale(*(tpY->xVectorPtr));
  nullVectorPtr->scale(*(tpY->nullVectorPtr));
  bifParam *= tpY->bifParam;
  return *this;
}

NOX::Abstract::Vector& TPBordVector::update(double alpha, const NOX::Abstract::Vector& y, double gamma)
{
  const TPBordVector* tpY = 0;
  tpY = dynamic_cast<const TPBordVector*>(&y);
  
  xVectorPtr->update(alpha, *(tpY->xVectorPtr), gamma);
  nullVectorPtr->update(alpha, *(tpY->nullVectorPtr), gamma);
  bifParam = alpha * tpY->bifParam + gamma*bifParam;
  return *this;
}

NOX::Abstract::Vector& TPBordVector::update(double alpha, 
			 const NOX::Abstract::Vector& y, 
                         double beta, const NOX::Abstract::Vector& b,
                         double gamma)
{
  const TPBordVector* tpY = 0;
  tpY = dynamic_cast<const TPBordVector*>(&y);

  const TPBordVector* tpB = 0;
  tpB = dynamic_cast<const TPBordVector*>(&b);  

  xVectorPtr->update(alpha, *(tpY->xVectorPtr), beta, 
		     *(tpB->xVectorPtr), gamma);

  nullVectorPtr->update(alpha, *(tpY->nullVectorPtr), beta, 
			*(tpB->nullVectorPtr), gamma);

  bifParam = alpha * tpY->bifParam + beta * tpB->bifParam + gamma * bifParam;
  return *this;
}

NOX::Abstract::Vector* TPBordVector::clone(NOX::CopyType type) const
{
  TPBordVector* newCopy = 0;

  newCopy = new TPBordVector(*this);
  
  if (type != NOX::DeepCopy)
    newCopy->init(0.0);
      
  return newCopy;
}

double TPBordVector::norm(NormType type) const
{
  double n, nx, nn;
  switch (type) {
  case NOX::Abstract::Vector::MaxNorm:
    nx = xVectorPtr->norm(NOX::Abstract::Vector::MaxNorm);
    nn = nullVectorPtr->norm(NOX::Abstract::Vector::MaxNorm);
    n = max(nx, nn);
    n = max(n, bifParam);
    break;
  case NOX::Abstract::Vector::OneNorm:
    nx = xVectorPtr->norm(NOX::Abstract::Vector::OneNorm);
    nn = nullVectorPtr->norm(NOX::Abstract::Vector::OneNorm);
    n = nx + nn + fabs(bifParam);
    break;
  case NOX::Abstract::Vector::TwoNorm:
  default:
   nx = xVectorPtr->norm(NOX::Abstract::Vector::TwoNorm);
   nn = nullVectorPtr->norm(NOX::Abstract::Vector::TwoNorm);
   n = sqrt(nx*nx + nn*nn + bifParam*bifParam);
   break;
  }
  return n;
}

double TPBordVector::norm(const NOX::Abstract::Vector& weights) const
{
  cout << "ERROR: LOCA::Bifurcation::TPBordVector::norm(weights) - "
       << "Not implemented yet!!!" << endl;
  throw "LOCA Error";
}

double TPBordVector::dot(const NOX::Abstract::Vector& y) const
{
  const TPBordVector* tpY = 0;
  tpY = dynamic_cast<const TPBordVector*>(&y);

  return( xVectorPtr->dot(*tpY->xVectorPtr) 
	  + nullVectorPtr->dot(*tpY->nullVectorPtr) 
	  + bifParam*bifParam);

}

int TPBordVector::length() const
{
  return (xVectorPtr->length() + nullVectorPtr->length() + 1);
}

bool TPBordVector::print() const
{
  xVectorPtr->print();
  nullVectorPtr->print();
  cout << "Bifurcation Parameter = " << bifParam << endl;
  return true;
}

void TPBordVector::setVec(const NOX::Abstract::Vector& xVec,
                    const NOX::Abstract::Vector& nullVec,
                    double bifPar)
{
  *xVectorPtr = xVec;
  *nullVectorPtr = nullVec;
  bifParam = bifPar;
}

  //! Returns the solution vector component of extended vector
const NOX::Abstract::Vector& TPBordVector::getXVec() const
{
  return *xVectorPtr;
}

  //! Returns the null vector component of extended vector
const NOX::Abstract::Vector& TPBordVector::getNullVec() const
{
  return *nullVectorPtr;
}

  //! Get Bifurcation parameter
double TPBordVector::getBifParam() const
{
  return bifParam;
}

  //! Returns the solution vector component of extended vector
NOX::Abstract::Vector& TPBordVector::getXVec()
{
  return *xVectorPtr;
}

  //! Returns the null vector component of extended vector
NOX::Abstract::Vector& TPBordVector::getNullVec()
{
  return *nullVectorPtr;
}

  //! Get Bifurcation parameter
double& TPBordVector::getBifParam()
{
  return bifParam;
}
