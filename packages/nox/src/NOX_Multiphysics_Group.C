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

#include "NOX_Multiphysics_Group.H"
#include <iostream>

NOX::Multiphysics::Group::Group(
          const Teuchos::RCP<std::vector<Teuchos::RCP<NOX::Solver::Generic> > >& solvers, 
          const Teuchos::RCP<NOX::StatusTest::Generic>& t, 
          const Teuchos::RCP<Teuchos::ParameterList>& p) :
  solversVecPtr(solvers),
  normRHS(0.0)
{
  // Create our version of the composite solution vector
  std::vector<const NOX::Abstract::Vector*> vecPtrs;

  for( unsigned int i = 0; i < solvers->size(); ++i )
  {
    std::cout << " .. .. .. received solver # " << i << std::endl;
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
      std::cerr << "ERROR: Invalid ConstructorType for group copy constructor." << std::endl;
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
  return (*solversVecPtr)[0]->getSolutionGroup().getF();
}

const NOX::Abstract::Vector& 
NOX::Multiphysics::Group::getGradient() const
{
  return (*solversVecPtr)[0]->getSolutionGroup().getGradient();
}

const NOX::Abstract::Vector& 
NOX::Multiphysics::Group::getNewton() const
{
  return (*solversVecPtr)[0]->getSolutionGroup().getNewton();
}

Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Multiphysics::Group::getXPtr() const
{
  return (*solversVecPtr)[0]->getSolutionGroup().getXPtr();
}

Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Multiphysics::Group::getFPtr() const
{
  return (*solversVecPtr)[0]->getSolutionGroup().getFPtr();
}

Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Multiphysics::Group::getGradientPtr() const
{
  return (*solversVecPtr)[0]->getSolutionGroup().getGradientPtr();
}

Teuchos::RCP< const NOX::Abstract::Vector >
NOX::Multiphysics::Group::getNewtonPtr() const
{
  return (*solversVecPtr)[0]->getSolutionGroup().getNewtonPtr();
}

double
NOX::Multiphysics::Group::getNormF() const
{
  if (!isF()) {
    std::cerr << "ERROR: NOX::Epetra::Group::getNormF() - invalid RHS" << std::endl;
    throw "NOX Error";
  }
    
  return normRHS;

}

Teuchos::RCP<NOX::Abstract::Group> 
NOX::Multiphysics::Group::clone(NOX::CopyType type) const 
{
  Teuchos::RCP<NOX::Abstract::Group> newgrp = 
    Teuchos::rcp(new NOX::Multiphysics::Group(*this, type));

  return newgrp;
}
