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

#include "LOCA_Epetra_Group.H"	          // class definition

#include "LOCA_Epetra_Interface.H"        // class data members
#include "NOX_Parameter_List.H"

LOCA::Epetra::Group::Group(NOX::Parameter::List& printParams,
			   NOX::Parameter::List& par, 
			   LOCA::Epetra::Interface& i, 
			   const LOCA::ParameterVector& p, 
			   NOX::Epetra::Vector& x, 
			   Epetra_Operator& J) :
  NOX::Epetra::Group(printParams, par, i, x, J),
  params(p),
  userInterface(i),
  scaleVec(x)
{
  computeScaleVec();  // use default method for computing scale vector
}

LOCA::Epetra::Group::Group(NOX::Parameter::List& printParams,
			   NOX::Parameter::List& par, 
			   LOCA::Epetra::Interface& i, 
			   const LOCA::ParameterVector& p, 
			   NOX::Epetra::Vector& x, 
			   Epetra_Operator& J, 
			   Epetra_Operator& M) :
  NOX::Epetra::Group(printParams, par, i, x, J, M),
  params(p),
  userInterface(i),
  scaleVec(x)
{
  computeScaleVec();  // use default method for computing scale vector
}

LOCA::Epetra::Group::Group(NOX::Parameter::List& printParams,
			   NOX::Parameter::List& par, 
			   LOCA::Epetra::Interface& i, 
			   const LOCA::ParameterVector& p, 
			   NOX::Epetra::Vector& x, 
			   Epetra_Operator& J,
			   const NOX::Epetra::Vector& s) :
  NOX::Epetra::Group(printParams, par, i, x, J),
  params(p),
  userInterface(i),
  scaleVec(s)
{
}

LOCA::Epetra::Group::Group(NOX::Parameter::List& printParams,
			   NOX::Parameter::List& par, 
			   LOCA::Epetra::Interface& i, 
			   const LOCA::ParameterVector& p, 
			   NOX::Epetra::Vector& x, 
			   Epetra_Operator& J, 
			   Epetra_Operator& M,
			   const NOX::Epetra::Vector& s) :
  NOX::Epetra::Group(printParams, par, i, x, J, M),
  params(p),
  userInterface(i),
  scaleVec(s)
{
}

LOCA::Epetra::Group::Group(const LOCA::Epetra::Group& source, 
			   NOX::CopyType type) :
  NOX::Epetra::Group(source, type),
  params(source.params),
  userInterface(source.userInterface),
  scaleVec(source.scaleVec)
{
}

LOCA::Epetra::Group::~Group() 
{
}

NOX::Abstract::Group* 
LOCA::Epetra::Group::clone(NOX::CopyType type) const 
{
  return new Group(*this, type);
}

NOX::Abstract::Group& 
LOCA::Epetra::Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

LOCA::Abstract::Group& 
LOCA::Epetra::Group::operator=(const LOCA::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

LOCA::Epetra::Group& 
LOCA::Epetra::Group::operator=(const LOCA::Epetra::Group& source)
{
  params = source.params;
  scaleVec = source.scaleVec;
  NOX::Epetra::Group::operator=(source);
  return *this;
}

void 
LOCA::Epetra::Group::setParams(const LOCA::ParameterVector& p)
{
  resetIsValid();
  params = p;
}

void
LOCA::Epetra::Group::setParam(int paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::Epetra::Group::getParam(int paramID) const
{
  return params.getValue(paramID);
}

void
LOCA::Epetra::Group::setParam(string paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::Epetra::Group::getParam(string paramID) const
{
  return params.getValue(paramID);
}

NOX::Abstract::Group::ReturnType
LOCA::Epetra::Group::computeF() 
{
  
  // Set the parameters prior to computing F
  userInterface.setParameters(params);
  
  return NOX::Epetra::Group::computeF();
}

NOX::Abstract::Group::ReturnType 
LOCA::Epetra::Group::computeJacobian() 
{
  
  // Set the parameters prior to computing F
  userInterface.setParameters(params);

  return NOX::Epetra::Group::computeJacobian();
}

const LOCA::ParameterVector& 
LOCA::Epetra::Group::getParams() const 
{
  return params;
}

NOX::Epetra::Interface& 
LOCA::Epetra::Group::getUserInterface()
{
  return userInterface;
}

void
LOCA::Epetra::Group::printSolution(const double conParam) const
{
  printSolution(xVector, conParam);
}

void
LOCA::Epetra::Group::printSolution(const NOX::Epetra::Vector& x_,
                                   const double conParam) const
{
  userInterface.printSolution(x_.getEpetraVector(), conParam);
}

void
LOCA::Epetra::Group::printSolution(const NOX::Abstract::Vector& x_,
                                   const double conParam) const
{
  printSolution(dynamic_cast<const NOX::Epetra::Vector&>(x_), conParam);
}

void
LOCA::Epetra::Group::setScaleVec(const NOX::Abstract::Vector& s) {
  setScaleVec( dynamic_cast<const NOX::Epetra::Vector&>(s) );
}

void
LOCA::Epetra::Group::setScaleVec(const NOX::Epetra::Vector& s) {
  scaleVec = s;
}

const NOX::Abstract::Vector&
LOCA::Epetra::Group::getScaleVec() const {
  return scaleVec;
}

void
LOCA::Epetra::Group::computeScaleVec() {
  scaleVec.init(1.0);
}
