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
#include "ml_common.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

// ----------   Includes   ----------
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

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
                  RefCountPtr<NLNML::NLNML_FineLevelNoxInterface> interface,
                  ParameterList& mlparams,
                  const Epetra_Comm& comm) : 
isinit_(false),
comm_(comm),                                
interface_(interface)
{
  label_  = "nlnML_Preconditioner";
  params_ = rcp(new Teuchos::ParameterList(mlparams));
  return;
}


/*----------------------------------------------------------------------*
 |  compute this preconditioner (public, derived)             m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_Preconditioner::computePreconditioner(
                                     const Epetra_Vector& x, 
				     Epetra_Operator& M,
				     NOX::Parameter::List* precParams)
{
  if (&M != this)
  {
    cout << "**ERR**: NLNML::NLNML_Preconditioner::computePreconditioner:\n"
         << "**ERR**: supplied preconditioner is not this\n"  
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  
  bool flag = true;
  int offset = getoffset();
  
  if (offset && ncalls_NewPrec_)
    if (ncalls_NewPrec_ % offset == 0)
     setinit(false);
  
  else if ( params_->get("nlnML adaptive recompute",0.0) &&
            ncalls_NewPrec_                              &&
            !params_->get("nlnML is linear preconditioner",true))
  {
    if (noxsolver_ == null)
    {
      cout << "**ERR**: NLNML::NLNML_Preconditioner::computePreconditioner:\n"
           << "**ERR**: outer nox solver not registered\n"  
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
    }
    const NOX::Epetra::Group& finalGroup =
      dynamic_cast<const NOX::Epetra::Group&>(noxsolver_->getSolutionGroup());
    const Epetra_Vector& currentF = 
      (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getF())).getEpetraVector();
    double norm;
    currentF.Norm2(&norm);
    if (norm > params_->get("nlnML adaptive recompute",0.0))
      setinit(false);    
  }  
  
  if (!isinit())
  {
    if (Comm().MyPID() && OutLevel())
      cout << "ML: NLNML_Preconditioner::computePreconditioner: (re)computing ML-Preconditioner\n";
      
    // save number of calls to computeF
    int ncalls = interface_->getnumcallscomputeF();
    
    interface_->setnumcallscomputeF(0);
    
    double t0 = GetClock();
    flag = compPrec(x);
    double t1 = GetClock();
    
    if (Comm().MyPID() && OutLevel())
    {
      cout << "ML: Setup time for preconditioner: " << (t1-t0) << " sec\n";
      cout << "ML: Number of calls to fineinterface.computeF() in setup: " 
           << interface_->getnumcallscomputeF() << endl;
    }
    
    // reset the number of calls to computeF
    interface_->setnumcallscomputeF(ncalls);
    
    if (flag) setinit(true);
    else
    {
      cout << "**ERR**: NLNML::NLNML_Preconditioner::computePreconditioner:\n"
           << "**ERR**: setup of Preconditioner failed\n"  
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
    }
    setinit(true);
  }
  ++ncalls_NewPrec_;
  return flag;  
}



/*----------------------------------------------------------------------*
 |  compute this preconditioner                                 m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_Preconditioner::compPrec(const Epetra_Vector& x)
{
  return true;
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
