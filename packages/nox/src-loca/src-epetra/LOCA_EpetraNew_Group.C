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

#include "LOCA_EpetraNew_Group.H"	          // class definition

#include "LOCA_EpetraNew_Interface_Required.H"        // class data members
#include "NOX_Parameter_List.H"
#include "Epetra_Vector.h"
#include "LOCA_Utils.H"
#include "LOCA_ErrorCheck.H"

LOCA::EpetraNew::Group::Group(NOX::Parameter::List& printingParams, 
			      LOCA::EpetraNew::Interface::Required& i, 
			      NOX::Epetra::Vector& initialGuess,
			      const LOCA::ParameterVector& p) :
  NOX::EpetraNew::Group(printingParams, i, initialGuess),
  LOCA::Abstract::Group(),
  params(p),
  userInterface(i),
  tmpVectorPtr2(0),
  scaleVecPtr(NULL)
{
}

LOCA::EpetraNew::Group::Group(NOX::Parameter::List& printingParams, 
			      LOCA::EpetraNew::Interface::Required& i, 
			      NOX::Epetra::Vector& initialGuess, 
			      NOX::EpetraNew::LinearSystem& linSys,
			      const LOCA::ParameterVector& p) :
  NOX::EpetraNew::Group(printingParams, i, initialGuess, linSys),
  LOCA::Abstract::Group(),
  params(p),
  userInterface(i),
  tmpVectorPtr2(0),
  scaleVecPtr(NULL)
{
}

LOCA::EpetraNew::Group::Group(const LOCA::EpetraNew::Group& source, 
			   NOX::CopyType type) :
  NOX::EpetraNew::Group(source, type),
  LOCA::Abstract::Group(source, type),
  params(source.params),
  userInterface(source.userInterface),
  tmpVectorPtr2(0),
  scaleVecPtr(NULL)
{
  if (source.scaleVecPtr != NULL)
    scaleVecPtr = source.scaleVecPtr->clone(NOX::DeepCopy);
}

LOCA::EpetraNew::Group::~Group() 
{
  delete tmpVectorPtr2;
  if (scaleVecPtr != NULL)
    delete scaleVecPtr;
}

NOX::Abstract::Group* 
LOCA::EpetraNew::Group::clone(NOX::CopyType type) const 
{
  return new Group(*this, type);
}

NOX::Abstract::Group& 
LOCA::EpetraNew::Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

LOCA::Abstract::Group& 
LOCA::EpetraNew::Group::operator=(const LOCA::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

NOX::Abstract::Group& 
LOCA::EpetraNew::Group::operator=(const NOX::EpetraNew::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

LOCA::EpetraNew::Group& 
LOCA::EpetraNew::Group::operator=(const LOCA::EpetraNew::Group& source)
{
  params = source.params;
  NOX::EpetraNew::Group::operator=(source);
  LOCA::Abstract::Group::operator=(source);
  if (source.scaleVecPtr != NULL)
    scaleVecPtr = source.scaleVecPtr->clone(NOX::DeepCopy);
  return *this;
}

void 
LOCA::EpetraNew::Group::setParams(const LOCA::ParameterVector& p)
{
  resetIsValid();
  params = p;
}

void
LOCA::EpetraNew::Group::setParam(int paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::EpetraNew::Group::getParam(int paramID) const
{
  return params.getValue(paramID);
}

void
LOCA::EpetraNew::Group::setParam(string paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::EpetraNew::Group::getParam(string paramID) const
{
  return params.getValue(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::EpetraNew::Group::computeF() 
{

  if (isF())
    return Abstract::Group::Ok;
  
  // Set the parameters prior to computing F
  userInterface.setParameters(params);
  
  return NOX::EpetraNew::Group::computeF();
}

NOX::Abstract::Group::ReturnType 
LOCA::EpetraNew::Group::computeJacobian() 
{

  if (isJacobian())
    return Abstract::Group::Ok;
  
  // Set the parameters prior to computing F
  userInterface.setParameters(params);

  return NOX::EpetraNew::Group::computeJacobian();
}

const LOCA::ParameterVector& 
LOCA::EpetraNew::Group::getParams() const 
{
  return params;
}

NOX::EpetraNew::Interface::Required& 
LOCA::EpetraNew::Group::getUserInterface()
{
  return userInterface;
}

void
LOCA::EpetraNew::Group::printSolution(const double conParam) const
{
  printSolution(xVector, conParam);
}

void
LOCA::EpetraNew::Group::printSolution(const NOX::Epetra::Vector& x_,
				      const double conParam) const
{
  userInterface.printSolution(x_.getEpetraVector(), conParam);
}

void
LOCA::EpetraNew::Group::printSolution(const NOX::Abstract::Vector& x_,
				      const double conParam) const
{
  printSolution(dynamic_cast<const NOX::Epetra::Vector&>(x_), conParam);
}

double
LOCA::EpetraNew::Group::computeScaledDotProduct(
				       const NOX::Abstract::Vector& a,
				       const NOX::Abstract::Vector& b) const
{
  if (scaleVecPtr == NULL)
    return a.dot(b) / a.length();
  else {
    NOX::Abstract::Vector* as = a.clone(NOX::DeepCopy);
    NOX::Abstract::Vector* bs = b.clone(NOX::DeepCopy);
    double d;

    as->scale(*scaleVecPtr);
    bs->scale(*scaleVecPtr);
    d = as->dot(*bs);

    delete as;
    delete bs;

    return d;
  }
}
   
void
LOCA::EpetraNew::Group::setScaleVector(const NOX::Abstract::Vector& s)
{
  if (scaleVecPtr != NULL)
    delete scaleVecPtr;

  scaleVecPtr = s.clone(NOX::DeepCopy);

  return;
}
