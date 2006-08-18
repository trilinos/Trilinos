//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#include "NOX_Multiphysics_Group.H"

NOX::Multiphysics::Group::Group(
          const Teuchos::RefCountPtr< vector<NOX::Solver::Manager*> >& solvers, 
          const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& t, 
          const Teuchos::RefCountPtr<Teuchos::ParameterList>& p) :
  solversVecPtr(solvers),
  normRHS(0.0)
{
  // Create our version of the composite solution vector
  vector<const NOX::Abstract::Vector*> vecPtrs;

  for( unsigned int i = 0; i < solvers->size(); ++i )
  {
    cout << " .. .. .. received solver # " << i << endl;
    vecPtrs.push_back( &((*solvers)[i]->getSolutionGroup().getX()) );
  }

  resetIsValid();
}

NOX::Multiphysics::Group::Group( const Group & source, NOX::CopyType type )
{
  switch (type) 
  {
    case DeepCopy:
      
      isValidRHS = source.isValidRHS;
      normRHS = source.normRHS;
      
      break;

    case ShapeCopy:
      resetIsValid();
      break;

    default:
      cerr << "ERROR: Invalid ConstructorType for group copy constructor." << endl;
      throw "NOX Error";
  }
}

NOX::Multiphysics::Group::~Group()
{
}

void 
NOX::Multiphysics::Group::resetIsValid()
{
  isValidRHS = false;
  return;
}

NOX::Abstract::Group & 
NOX::Multiphysics::Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const NOX::Multiphysics::Group&> (source));
}

NOX::Abstract::Group & 
NOX::Multiphysics::Group::operator=(const Group& source)
{
  return *this;
}

void
NOX::Multiphysics::Group::setX(const NOX::Abstract::Vector& y)
{
  resetIsValid();
}

void
NOX::Multiphysics::Group::computeX(const NOX::Abstract::Group& grp, 
			 const NOX::Abstract::Vector& d, double step)
{
  resetIsValid();
}

NOX::Abstract::Group::ReturnType 
NOX::Multiphysics::Group::computeF()
{
  NOX::Abstract::Group::ReturnType status;

  for( unsigned int i = 0; i < (*solversVecPtr).size(); ++i )
  {
    status = const_cast<NOX::Abstract::Group&>((*solversVecPtr)[i]->getSolutionGroup()).computeF();
    if( NOX::Abstract::Group::Ok != status )
      return status;
  }

  isValidRHS = true;

  // Compute the composite normF
  normRHS = 0;

  for( unsigned int i = 0; i < (*solversVecPtr).size(); ++i )
    normRHS += (*solversVecPtr)[i]->getSolutionGroup().getNormF() *
             (*solversVecPtr)[i]->getSolutionGroup().getNormF()  ;
  normRHS = sqrt(normRHS);

  return NOX::Abstract::Group::Ok;
}

bool
NOX::Multiphysics::Group::isF() const
{
  return isValidRHS;
}

const NOX::Abstract::Vector& 
NOX::Multiphysics::Group::getX() const
{
  return (*solversVecPtr)[0]->getSolutionGroup().getX();
}

const NOX::Abstract::Vector& 
NOX::Multiphysics::Group::getF() const
{
  return (*solversVecPtr)[0]->getSolutionGroup().getX();
}

const NOX::Abstract::Vector& 
NOX::Multiphysics::Group::getGradient() const
{
  return (*solversVecPtr)[0]->getSolutionGroup().getX();
}

const NOX::Abstract::Vector& 
NOX::Multiphysics::Group::getNewton() const
{
  return (*solversVecPtr)[0]->getSolutionGroup().getX();
}

double
NOX::Multiphysics::Group::getNormF() const
{
  if (!isF()) {
    cerr << "ERROR: NOX::Epetra::Group::getNormF() - invalid RHS" << endl;
    throw "NOX Error";
  }
    
  return normRHS;

}

Teuchos::RefCountPtr<NOX::Abstract::Group> 
NOX::Multiphysics::Group::clone(NOX::CopyType type) const 
{
  Teuchos::RefCountPtr<NOX::Abstract::Group> newgrp = 
    Teuchos::rcp(new NOX::Multiphysics::Group(*this, type));

  return newgrp;
}
