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

#include "NOX_Parameter_List.H"
#include "LOCA_BorderedSystem_Manager.H"
#include "LOCA_BorderedSystem_Bordering.H"
#include "LOCA_Utils.H"

LOCA::BorderedSystem::Manager::Manager(NOX::Parameter::List& params): 
  method(),
  borderedPtr(NULL)
{
  reset(params);
}

LOCA::BorderedSystem::Manager::Manager(
			       const LOCA::BorderedSystem::Manager& source) :
  method(source.method),
  borderedPtr(source.borderedPtr->clone())
{
}

LOCA::BorderedSystem::Manager::~Manager()
{
  delete borderedPtr;
}

LOCA::BorderedSystem::Manager&
LOCA::BorderedSystem::Manager::operator=(
				const LOCA::BorderedSystem::Manager& source)
{
  if (this != &source) {
    method = source.method;
    *borderedPtr = *source.borderedPtr;
  }
  
  return *this;
}

LOCA::BorderedSystem::Generic*
LOCA::BorderedSystem::Manager::clone() const
{
  return new Manager(*this);
}

LOCA::BorderedSystem::Generic& 
LOCA::BorderedSystem::Manager::operator=(
				  const LOCA::BorderedSystem::Generic& source)
{
  return 
    operator=(dynamic_cast<const LOCA::BorderedSystem::Manager&>(source));
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::Manager::reset(NOX::Parameter::List& params)
{
  string newmethod = params.getParameter("Bordered Solver Method", 
					 "Bordering");

  if (method != newmethod) {
    delete borderedPtr;

    method = newmethod;

    if (method == "Bordering")
      borderedPtr = new LOCA::BorderedSystem::Bordering(params);
    else {
      if (LOCA::Utils::doPrint(LOCA::Utils::Error)) {
	cout << "LOCA::BorderedSystem::Manager::reset() - invalid choice (" 
	     << method << ") for bordered solver method " << endl;
      }
      return NOX::Abstract::Group::Failed;
    }
  }
  return NOX::Abstract::Group::Ok;
}

void
LOCA::BorderedSystem::Manager::setIsZero(bool flagA, bool flagB, bool flagC,
					 bool flagF, bool flagG)
{
  borderedPtr->setIsZero(flagA, flagB, flagC, flagF, flagG);
}

void
LOCA::BorderedSystem::Manager::setIsContiguous(bool flag)
{
  borderedPtr->setIsContiguous(flag);
}

void
LOCA::BorderedSystem::Manager::setMatrixBlocks(
			const NOX::Abstract::Group* group,
			const NOX::Abstract::MultiVector* blockA,
			const NOX::Abstract::MultiVector* blockB,
			const NOX::Abstract::MultiVector::DenseMatrix* blockC)
{
  borderedPtr->setMatrixBlocks(group, blockA, blockB, blockC);
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::Manager::apply(
			  const NOX::Abstract::MultiVector& X,
			  const NOX::Abstract::MultiVector::DenseMatrix& Y,
			  NOX::Abstract::MultiVector& U,
			  NOX::Abstract::MultiVector::DenseMatrix& V) const
{
  return borderedPtr->apply(X, Y, U, V);
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::Manager::applyTranspose(
			  const NOX::Abstract::MultiVector& X,
			  const NOX::Abstract::MultiVector::DenseMatrix& Y,
			  NOX::Abstract::MultiVector& U,
			  NOX::Abstract::MultiVector::DenseMatrix& V) const
{
  return borderedPtr->applyTranspose(X, Y, U, V);
}

NOX::Abstract::Group::ReturnType 
LOCA::BorderedSystem::Manager::applyInverse(
			      NOX::Parameter::List& params,
			      const NOX::Abstract::MultiVector* F,
			      const NOX::Abstract::MultiVector::DenseMatrix* G,
			      NOX::Abstract::MultiVector& X,
			      NOX::Abstract::MultiVector::DenseMatrix& Y) const
{
  return borderedPtr->applyInverse(params, F, G, X, Y);
}
