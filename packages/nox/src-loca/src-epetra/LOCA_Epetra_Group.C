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
//#include "NOX_Epetra_SharedOperator.H"
#include "NOX_Parameter_List.H"

using namespace LOCA;
using namespace LOCA::Epetra;

Group::Group(const NOX::Parameter::List& params, Interface& i, 
	     const ParameterVector& p, NOX::Epetra::Vector& x, 
	     Epetra_Operator& J) :
  NOX::Epetra::Group(params, i, x, J),
  pVectorPtr(new ParameterVector(p)),
  userInterface(i)
{
}

Group::Group(const NOX::Parameter::List& params, Interface& i, 
	     const ParameterVector& p, NOX::Epetra::Vector& x, 
	     Epetra_Operator& J, 
	     Epetra_Operator& M) :
  NOX::Epetra::Group(params, i, x, J, M),
  pVectorPtr(new ParameterVector(p)),
  userInterface(i)
{
}

Group::Group(const Group& source, NOX::CopyType type) :
  NOX::Epetra::Group(source, type),
  pVectorPtr(new ParameterVector(*(source.pVectorPtr))),
  userInterface(source.userInterface)
{
}

Group::~Group() 
{
  delete pVectorPtr;
}

NOX::Abstract::Group* Group::clone(NOX::CopyType type) const 
{
  return new Group(*this, type);
}

NOX::Abstract::Group& Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

Abstract::Group& Group::operator=(const Abstract::Group& source)
{
  return operator=(dynamic_cast<const Group&> (source));
}

Group& Group::operator=(const Group& source)
{
  *pVectorPtr = *source.pVectorPtr;
  NOX::Epetra::Group::operator=(source);
  return *this;
}

bool Group::setParams(const ParameterVector& p)
{
  *pVectorPtr = p;
  resetIsValid();
  return true;
}

bool Group::computeParams(const ParameterVector& oldParams,
			  const ParameterVector& direction, 
			  double step)
{
  *pVectorPtr = oldParams;
  pVectorPtr->update(step, direction, 1.0);
  resetIsValid();
  return true;
}


bool Group::computeF() 
{

  bool status = false;
  
  // Set the parameters prior to computing F
  status = userInterface.setParameters(*pVectorPtr);

  if (status == false)
    return false;
  
  return NOX::Epetra::Group::computeF();;
}

bool Group::computeJacobian() 
{
  bool status = false;
  
  // Set the parameters prior to computing F
  status = userInterface.setParameters(*pVectorPtr);

  if (status == false)
    return false;
  
  return NOX::Epetra::Group::computeJacobian();
}

const ParameterVector& Group::getParams() const 
{
  return *pVectorPtr;
}

NOX::Epetra::Interface& Group::getUserInterface()
{
  return userInterface;
}
