// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Multiphysics_Group.H"
#include <iostream>

NOX::Multiphysics::Group::Group(
          const Teuchos::RCP<std::vector<Teuchos::RCP<NOX::Solver::Generic> > >& solvers,
          const Teuchos::RCP<NOX::StatusTest::Generic>& /* t */,
          const Teuchos::RCP<Teuchos::ParameterList>& /* p */) :
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

NOX::Multiphysics::Group::Group( const Group & source, NOX::CopyType type ):
  normRHS(0.0)
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
      throw std::runtime_error("NOX Error");
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
NOX::Multiphysics::Group::operator=(const Group& /* source */)
{
  return *this;
}

void
NOX::Multiphysics::Group::setX(const NOX::Abstract::Vector& /* y */)
{
  resetIsValid();
}

void
NOX::Multiphysics::Group::computeX(const NOX::Abstract::Group& /* grp */,
             const NOX::Abstract::Vector& /* d */, double /* step */)
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
    throw std::runtime_error("NOX Error");
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
