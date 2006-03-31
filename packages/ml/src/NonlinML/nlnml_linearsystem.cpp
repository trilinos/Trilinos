/*
#@HEADER
# ************************************************************************
#
#               ML: A Multilevel Preconditioner Package
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
// ML-headers
#include "ml_common.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

// ----------   Includes   ----------
#include <ctime>
#include <iostream>


// ----------   User Defined Includes   ----------
#include "nlnml_linearsystem.H"
#include "nlnml_preconditioner.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_LinearSystem::NLNML_LinearSystem()
{
}

/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::applyJacobian(
                                         const NOX::Epetra::Vector& input, 
	                                 NOX::Epetra::Vector& result) const
{

  return true;
}                     


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::applyJacobianTranspose(
                                         const NOX::Epetra::Vector& input, 
		                         NOX::Epetra::Vector& result) const
{

  return true;
}                     



/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::applyJacobianInverse(
                                         NOX::Parameter::List &params, 
		                         const NOX::Epetra::Vector &input, 
		                         NOX::Epetra::Vector &result)
{

  return true;
}                     



/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::applyRightPreconditioning(
                                 bool useTranspose,
				 NOX::Parameter::List& params, 
				 const NOX::Epetra::Vector& input, 
				 NOX::Epetra::Vector& result) const
{

  return true;
}                     


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
Teuchos::RefCountPtr<NOX::Epetra::Scaling>  NLNML::NLNML_LinearSystem::getScaling()
{

  return Teuchos::null;
}                     


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_LinearSystem::resetScaling(
                     const Teuchos::RefCountPtr<NOX::Epetra::Scaling>& s)
{

  return;
}                     


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::computeJacobian(const NOX::Epetra::Vector& x)
{

  return true;
}                     



/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::createPreconditioner(
                            const NOX::Epetra::Vector& x, 
			    NOX::Parameter::List& p,
			    bool recomputeGraph)
{

  return true;
}                     



/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::destroyPreconditioner()
{

  return true;
}                     


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::checkPreconditionerReuse()
{

  return true;
}                     


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::hasPreconditioner() const
{

  return false;
}                     


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_LinearSystem::isPreconditionerConstructed() const
{

  return false;
}                     



/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
Teuchos::RefCountPtr<const Epetra_Operator> NLNML::NLNML_LinearSystem::getJacobianOperator() const
{

  return Teuchos::null;
}                     



/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
Teuchos::RefCountPtr<const Epetra_Operator> NLNML::NLNML_LinearSystem::getGeneratedPrecOperator() const
{

  return Teuchos::null;
}                     



/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
Teuchos::RefCountPtr<Epetra_Operator> NLNML::NLNML_LinearSystem::getGeneratedPrecOperator()
{

  return Teuchos::null;
}                     


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_LinearSystem::setJacobianOperatorForSolve(
                const Teuchos::RefCountPtr<const Epetra_Operator>& solveJacOp)
{

  return;
}                     


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_LinearSystem::setPrecOperatorForSolve(
                 const Teuchos::RefCountPtr<const Epetra_Operator>& solvePrecOp)
{

  return;
}                     




#endif









