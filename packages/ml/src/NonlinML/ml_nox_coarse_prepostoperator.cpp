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
