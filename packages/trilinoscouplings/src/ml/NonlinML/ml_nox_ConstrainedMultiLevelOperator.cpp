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
