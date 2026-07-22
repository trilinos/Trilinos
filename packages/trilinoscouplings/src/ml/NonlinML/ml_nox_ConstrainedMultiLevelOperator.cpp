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


// ----------   User Defined Includes   ----------
#include "ml_nox_ConstrainedMultiLevelOperator.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            m.gee 11/05|
 *----------------------------------------------------------------------*/
ML_NOX::ML_Nox_ConstrainedMultiLevelOperator::ML_Nox_ConstrainedMultiLevelOperator(
ML_Epetra::MultiLevelOperator* ml_operator, 
ML_NOX::Nox_CoarseProblem_Interface& coarseinterface) :
comm_(ml_operator->Comm()),
coarseinterface_(coarseinterface),
ml_operator_(ml_operator)
{
  label_  = "ML_Nox_ConstrainedMultiLevelOperator";
  usetranspose_ = false;
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            m.gee 11/05|
 *----------------------------------------------------------------------*/
ML_NOX::ML_Nox_ConstrainedMultiLevelOperator::~ML_Nox_ConstrainedMultiLevelOperator() 
{
  if (ml_operator_) delete ml_operator_; ml_operator_ = NULL;
}

/*----------------------------------------------------------------------*
 |  apply multigrid linear preconditioner (public)           m.gee 11/05|
 *----------------------------------------------------------------------*/
int ML_NOX::ML_Nox_ConstrainedMultiLevelOperator::ApplyInverse(
                     const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
{
  int err = ml_operator_->ApplyInverse(X,Y);
#if 0
  Epetra_Vector tmp(View,Y,0);
  coarseinterface_.ApplyAllConstraints(tmp);
#endif
  return err;
}



#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)
