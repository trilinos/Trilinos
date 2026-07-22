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
/*!
 * \file nlnml_preconditioner1.cpp
 *
 * \class ML_Nox_Preconditioner
 *
 * \brief ML nonlinear preconditioner and solver
 *
 * \date Last update do Doxygen: 31-Mar-05
 *
 */
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

// ----------   Includes   ----------
#include <ctime>
#include <iostream>

// this class
#include "nlnml_ConstrainedMultiLevelOperator.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_ConstrainedMultiLevelOperator::NLNML_ConstrainedMultiLevelOperator(
              RefCountPtr<ML_Epetra::MultiLevelOperator> ml_operator,
              RefCountPtr<NLNML::NLNML_CoarseLevelNoxInterface> coarseinterface,
              bool applyconstraints) :
comm_(ml_operator->Comm()),
coarseinterface_(coarseinterface),
ml_operator_(ml_operator),
applyconstraints_(applyconstraints)
{
  label_ = "NLNML_ConstrainedMultiLevelOperator";
  return;
}


/*----------------------------------------------------------------------*
 |                                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
int NLNML::NLNML_ConstrainedMultiLevelOperator::ApplyInverse(
                     const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int err = ml_operator_->ApplyInverse(X,Y);
  if (applyconstraints_)
  {
    Epetra_Vector tmp(View,Y,0);
    coarseinterface_->ApplyAllConstraints(tmp);
  }
  return err;
}















#endif
