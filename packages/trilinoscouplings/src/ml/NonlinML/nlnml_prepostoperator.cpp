// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
