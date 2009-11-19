// $Id: LOCA_PhaseTransition_ExtendedVector.C,v 1.7 2007/06/21 16:22:52 rhoope Exp $ 
// $Source: /space/CVS/Trilinos/packages/nox/src-loca/src/LOCA_PhaseTransition_ExtendedVector.C,v $ 

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
//  $Source: /space/CVS/Trilinos/packages/nox/src-loca/src/LOCA_PhaseTransition_ExtendedVector.C,v $
//  $Author: rhoope $
//  $Date: 2007/06/21 16:22:52 $
//  $Revision: 1.7 $
// ************************************************************************
//@HEADER

#include "LOCA_PhaseTransition_ExtendedVector.H"  // Class definition
#include "LOCA_PhaseTransition_ExtendedMultiVector.H"

LOCA::PhaseTransition::ExtendedVector::ExtendedVector(
		    const Teuchos::RCP<LOCA::GlobalData>& global_data,
		    const NOX::Abstract::Vector& x1Vec,
		    const NOX::Abstract::Vector& x2vec,
		    double ptp) :
  LOCA::Extended::Vector(global_data,2,1)
{
  setVector(0, x1Vec);
  setVector(1, x2vec);
  setScalar(0, ptp);
}

LOCA::PhaseTransition::ExtendedVector::ExtendedVector(
                const LOCA::PhaseTransition::ExtendedVector& source,
		NOX::CopyType type) :
  LOCA::Extended::Vector(source, type)
{
}


LOCA::PhaseTransition::ExtendedVector::~ExtendedVector()
{
}

NOX::Abstract::Vector& 
LOCA::PhaseTransition::ExtendedVector::operator=(
					      const NOX::Abstract::Vector& y)
{
  operator=(dynamic_cast<const LOCA::PhaseTransition::ExtendedVector&>(y));
  return *this;
}

LOCA::Extended::Vector& 
LOCA::PhaseTransition::ExtendedVector::operator=(
					     const LOCA::Extended::Vector& y)
{
  operator=(dynamic_cast<const LOCA::PhaseTransition::ExtendedVector&>(y));
  return *this;
}

LOCA::PhaseTransition::ExtendedVector& 
LOCA::PhaseTransition::ExtendedVector::operator=(
                     const LOCA::PhaseTransition::ExtendedVector& y)
{ 
  LOCA::Extended::Vector::operator=(y);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::PhaseTransition::ExtendedVector::clone(
						    NOX::CopyType type) const
{
  return 
    Teuchos::rcp(new LOCA::PhaseTransition::ExtendedVector(*this, 
								     type));
}

void 
LOCA::PhaseTransition::ExtendedVector::setVec(
					const NOX::Abstract::Vector& x1Vec,
					const NOX::Abstract::Vector& x2vec,
					double ptp)
{
  setVector(0, x1Vec);
  setVector(1, x2vec);
  setScalar(0, ptp);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::PhaseTransition::ExtendedVector::X1() const
{
  return getVector(0);
}

Teuchos::RCP<const NOX::Abstract::Vector>
LOCA::PhaseTransition::ExtendedVector::X2() const
{
  return getVector(1);
}

double 
LOCA::PhaseTransition::ExtendedVector::PTP() const
{
  return getScalar(0);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::PhaseTransition::ExtendedVector::X1()
{
  return getVector(0);
}

Teuchos::RCP<NOX::Abstract::Vector>
LOCA::PhaseTransition::ExtendedVector::X2()
{
  return getVector(1);
}

double& 
LOCA::PhaseTransition::ExtendedVector::PTP()
{
  return getScalar(0);
}

LOCA::PhaseTransition::ExtendedVector::ExtendedVector(
		  const Teuchos::RCP<LOCA::GlobalData>& global_data) :
  LOCA::Extended::Vector(global_data,2,1)
{
}

Teuchos::RCP<LOCA::Extended::MultiVector>
LOCA::PhaseTransition::ExtendedVector::generateMultiVector(
							int nColumns, 
							int nVectorRows, 
							int nScalarRows) const
{
  return 
    Teuchos::rcp(new LOCA::PhaseTransition::ExtendedMultiVector(
								    globalData,
								    nColumns));
}
