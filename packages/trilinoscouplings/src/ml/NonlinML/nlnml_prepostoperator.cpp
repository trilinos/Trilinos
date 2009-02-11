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
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
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
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

// ----------   Includes   ----------
#include <ctime>
#include <iostream>
#include "nlnml_prepostoperator.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_PrePostOperator::NLNML_PrePostOperator(
              RefCountPtr<NLNML::NLNML_CoarseLevelNoxInterface> coarseinterface,
              RefCountPtr<NLNML::NLNML_FineLevelNoxInterface> finterface) :
fineinterface_(finterface),
coarseinterface_(coarseinterface)
{
  type_ = "NLNML::NLNML_PrePostOperator";
  return;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_PrePostOperator::NLNML_PrePostOperator(
              const NLNML::NLNML_PrePostOperator& old) :
type_(old.type_),
fineinterface_(old.fineinterface_),
coarseinterface_(old.coarseinterface_)
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_PrePostOperator::~NLNML_PrePostOperator()
{
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_PrePostOperator::runPreSolve(const NOX::Solver::Generic& solver)
{
  const NOX::Epetra::Group& solgroup = 
             dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
  const NOX::Epetra::Vector* noxsolution = 
             dynamic_cast<const NOX::Epetra::Vector*>(&solgroup.getX());
  Epetra_Vector& epetrasolution = 
             const_cast<Epetra_Vector&>(noxsolution->getEpetraVector());
  
  coarseinterface_->ApplyAllConstraints(epetrasolution);
  
  // put the vector back into place
  const NOX::Epetra::Group& solgroup2 = 
             dynamic_cast<const NOX::Epetra::Group&>(solver.getSolutionGroup());
  const NOX::Epetra::Vector* noxsolution2 = 
             dynamic_cast<const NOX::Epetra::Vector*>(&solgroup2.getX());
  Epetra_Vector& epetrasolution2 = 
             const_cast<Epetra_Vector&>(noxsolution2->getEpetraVector());
  epetrasolution2.Scale(1.0,epetrasolution); 
  return;
}










#endif
