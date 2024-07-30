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
#include "ml_nox_coarse_prepostoperator.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            m.gee 10/05|
 *----------------------------------------------------------------------*/
ML_NOX::Ml_Nox_CoarsePrePostOperator::Ml_Nox_CoarsePrePostOperator(
                               ML_NOX::Nox_CoarseProblem_Interface& coarseinterface,
                               ML_NOX::Ml_Nox_Fineinterface& fineinterface) :
coarseinterface_(coarseinterface),
fineinterface_(fineinterface)                               
{
  type_ = "Ml_Nox_CoarsePrePostOperator";
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       m.gee 10/05|
 *----------------------------------------------------------------------*/
ML_NOX::Ml_Nox_CoarsePrePostOperator::Ml_Nox_CoarsePrePostOperator(const ML_NOX::Ml_Nox_CoarsePrePostOperator& source) :
coarseinterface_(source.coarseinterface_),
fineinterface_(source.fineinterface_)                               
{
  type_ = source.type_;
  return;
}

/*----------------------------------------------------------------------*
 |  clone (public)                                           m.gee 10/05|
 *----------------------------------------------------------------------*/
NOX::Parameter::Arbitrary* ML_NOX::Ml_Nox_CoarsePrePostOperator::clone() const
{
  ML_NOX::Ml_Nox_CoarsePrePostOperator* tmp 
    = new ML_NOX::Ml_Nox_CoarsePrePostOperator(*this);
  return tmp;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            m.gee 10/05|
 *----------------------------------------------------------------------*/
ML_NOX::Ml_Nox_CoarsePrePostOperator::~Ml_Nox_CoarsePrePostOperator()
{ 
   return; 
}

/*----------------------------------------------------------------------*
 |       (public)                                            m.gee 10/05|
 *----------------------------------------------------------------------*/
void ML_NOX::Ml_Nox_CoarsePrePostOperator::runPreSolve(const NOX::Solver::Generic& solver)
{
#if 0
  cout << "Ml_Nox_CoarsePrePostOperator::runPreSolve called\n";
  
  // get the current solution out of the solver
  NOX::Abstract::Group& solgroup = const_cast<NOX::Abstract::Group&>(solver.getSolutionGroup());
  const NOX::Epetra::Vector* noxsolution = 
    dynamic_cast<const NOX::Epetra::Vector*>(&solgroup.getX());
  Epetra_Vector& epetrasolution = const_cast<Epetra_Vector&>(noxsolution->getEpetraVector());

#if 1 // this is expensive and might be removed
  // test whether we habe a level 0 (fine level) vector here
  if (!epetrasolution.Map().PointSameAs(fineinterface_.getMap()))
  {
    cout << "**ERR**: Ml_Nox_CoarsePrePostOperator::runPreSolve:\n";
    cout << "**ERR**: incoming vector map does not match fine interface.Map()\n";
    cout << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; throw -1;
  }
#endif  

  // make the application apply its constraints to this vector
  fineinterface_.ApplyAllConstraints(epetrasolution);
  
  // put the vector back into the solver
  solgroup = const_cast<NOX::Abstract::Group&>(solver.getSolutionGroup());
  noxsolution = dynamic_cast<const NOX::Epetra::Vector*>(&solgroup.getX());
  Epetra_Vector& epetrasolutionnew = const_cast<Epetra_Vector&>(noxsolution->getEpetraVector());

  epetrasolutionnew.Scale(1.0,epetrasolution);  
      
#endif
  return;
}
/*----------------------------------------------------------------------*
 |       (public)                                            m.gee 10/05|
 *----------------------------------------------------------------------*/
void ML_NOX::Ml_Nox_CoarsePrePostOperator::runPostSolve(const NOX::Solver::Generic& solver)
{
#if 0
  cout << "Ml_Nox_CoarsePrePostOperator::runPostSolve called\n";
#endif
  return;
}

/*----------------------------------------------------------------------*
 |       (public)                                            m.gee 10/05|
 *----------------------------------------------------------------------*/
void ML_NOX::Ml_Nox_CoarsePrePostOperator::runPreIterate(const NOX::Solver::Generic& solver)
{
#if 0
  cout << "Ml_Nox_CoarsePrePostOperator::runPreIterate called\n";
#endif  
  return;
}

/*----------------------------------------------------------------------*
 |       (public)                                            m.gee 10/05|
 *----------------------------------------------------------------------*/
void ML_NOX::Ml_Nox_CoarsePrePostOperator::runPostIterate(const NOX::Solver::Generic& solver)
{
#if 0
  cout << "Ml_Nox_CoarsePrePostOperator::runPostIterate called\n";
#endif
  return;
}


#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)
