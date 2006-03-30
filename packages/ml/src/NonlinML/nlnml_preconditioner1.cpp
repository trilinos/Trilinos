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
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

// ----------   Includes   ----------
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "ml_common.h"
#ifdef HAVE_ML_MLAPI
#include "MLAPI_Space.h"
#include "MLAPI_Operator.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Gallery.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_MultiLevelAdaptiveSA.h" 
#include "MLAPI_DistributedMatrix.h"
#endif

// this class
#include "nlnml_preconditioner.H" 

/*----------------------------------------------------------------------*
 |  ctor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_Preconditioner::NLNML_Preconditioner(
                                NLNML::NLNML_FineLevelNoxInterface& interface,
                                Teuchos::ParameterList& mlparams,
                                const Epetra_Comm& comm) : 
Epetra_Operator(),
NOX::Epetra::Interface::Preconditioner(),                                
interface_(interface),
comm_(comm)                                
{
  params_ = Teuchos::rcp(new Teuchos::ParameterList(mlparams));
  return;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_Preconditioner::NLNML_Preconditioner(int dummy)
{
  cout << "Dummy constructor called\n"; fflush(stdout);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_Preconditioner::~NLNML_Preconditioner()
{
  return;
}


/*----------------------------------------------------------------------*
 |  register outer nox solver                                 m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_Preconditioner::SetNoxSolver()
{
  return true;
}


/*----------------------------------------------------------------------*
 |  compute this preconditioner (public, derived)             m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_Preconditioner::computePreconditioner(
                                     const Epetra_Vector& x, 
				     Epetra_Operator& M,
				     NOX::Parameter::List* precParams)
{
  return;
}



/*----------------------------------------------------------------------*
 |  compute this preconditioner                                 m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_Preconditioner::compPrec()
{
  return ;
}



/*----------------------------------------------------------------------*
 |  compute linear preconditioner                             m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_Preconditioner::ML_Nox_compute_Jacobian_Linearpreconditioner()
{
  return ;
}



/*----------------------------------------------------------------------*
 |  compute nonlinear preconditioner                          m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_Preconditioner::ML_Nox_compute_Jacobian_Nonlinearpreconditioner()
{
  return ;
}



/*----------------------------------------------------------------------*
 |  apply this preconditioner (public, derived)               m.gee 3/06|
 *----------------------------------------------------------------------*/
int NLNML::NLNML_Preconditioner::ApplyInverse(
                                     const Epetra_MultiVector& X, 
                                     Epetra_MultiVector& Y) const
{
  return 0;
}



/*----------------------------------------------------------------------*
 |  apply inverse for linear preconditioner                   m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_Preconditioner::ML_Nox_ApplyInverse_Linear() const
{
  return ;
}



/*----------------------------------------------------------------------*
 |  apply inverse for nonlinear preconditioner                   m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_Preconditioner::ML_Nox_ApplyInverse_NonLinear() const
{
  return ;
}




/*----------------------------------------------------------------------*
 |  run FAS preconditioner as a solver                        m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_Preconditioner::solve() const
{
  return ;
}




/*----------------------------------------------------------------------*
 |  choose smoothers                                          m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_Preconditioner::Set_Smoothers()
{
  return ;
}



/*----------------------------------------------------------------------*
 |  fix main diagonal of Jacobian                             m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_Preconditioner::fix_MainDiagonal()
{
  return ;
}




/*----------------------------------------------------------------------*
 |  compute Jacobian on fine grid                             m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_Preconditioner::ML_Nox_computeFineLevelJacobian()
{
  return ;
}



/*----------------------------------------------------------------------*
 |  FAS-preconditioner                                        m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_Preconditioner::ML_Nox_FAS_cycle() const
{
  return ;
}



/*----------------------------------------------------------------------*
 |  FAS-solver                                                m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_Preconditioner::ML_Nox_FAS_cycle1() const
{
  return ;
}



/*----------------------------------------------------------------------*
 |  adaptive setup                             m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_Preconditioner::Ml_Nox_adaptivesetup()
{
  return ;
}
















#endif
