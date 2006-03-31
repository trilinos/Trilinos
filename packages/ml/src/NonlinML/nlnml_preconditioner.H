/*
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
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file nlnml_preconditioner.H
 *
 * \class NLNML_Preconditioner
 *
 * \brief ML nonlinear preconditioner and solver
 *
 * \date Last update do Doxygen: 31-Mar-05
 *
 */
#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef NLNML_PRECONDITIONER_H
#define NLNML_PRECONDITIONER_H

// ML-headers
#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelOperator.h"
#include "ml_agg_VBMETIS.h"

// epetra-headers
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialComm.h"

// Teuchos
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"


// NOX-headers
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "NOX_Epetra_Group.H"
#include "NOX_Solver_Manager.H"

// own headers
#include "nlnml_preconditioner_utils.H"
#include "nlnml_finelevelnoxinterface.H"

// ---------- Forward Declarations ----------
class NLNML::NLNML_FineLevelNoxInterface;

using namespace std;
using namespace Teuchos;
namespace NLNML
{

class NLNML_Preconditioner : public virtual Epetra_Operator,
                             public virtual NOX::Epetra::Interface::Preconditioner
{
public:

  explicit NLNML_Preconditioner(
                     RefCountPtr<NLNML::NLNML_FineLevelNoxInterface> interface,
                     ParameterList& mlparams,
                     const Epetra_Comm& comm);
  
  virtual ~NLNML_Preconditioner() {}
  
  void SetNoxSolver(RefCountPtr<NOX::Solver::Manager> s) { noxsolver_ = s; } 
  
  const char* Label() const {return(&(label_[0]));}
  
  const Epetra_Comm & Comm() const {return(comm_);}
  
  const Epetra_Map & OperatorDomainMap() const {return interface_->getMap();}
  
  const Epetra_Map & OperatorRangeMap() const {return interface_->getMap();}
  
  int OutLevel() { return params_->get("nlnML output",0); } 
  
  bool isinit() const {return(isinit_);}
  
  bool setinit(bool flag) 
  {isinit_=flag; ncalls_NewPrec_ = 0; return(true);}
  
  int getoffset() const {return params_->get("nlnML offset recompute",0);} 
  
  virtual bool computePreconditioner(const Epetra_Vector& x, 
				     Epetra_Operator& M,
				     NOX::Parameter::List* precParams = 0);
  
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const ;
  
  void solve() const ;
  
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
  { cout << "ML_Nox_Preconditioner::Apply does not make sense for an ML-Multilevel Operator\n"
         << "Use ApplyInverse\n";return(-1);}
  
  int SetUseTranspose(bool UseTranspose) { return -1; }

  double NormInf() const {return(-1.0);}
  
  bool UseTranspose() const {return false;}
  
  bool HasNormInf() const {return(false);}

private:

  explicit NLNML_Preconditioner(NLNML::NLNML_Preconditioner& old);
  NLNML_Preconditioner operator = (const NLNML_Preconditioner& old);
    
  void Initialize();
  
  bool compPrec(const Epetra_Vector& x);
  void ML_Nox_compute_Jacobian_Linearpreconditioner();      
  void ML_Nox_compute_Jacobian_Nonlinearpreconditioner();   

  void Set_Smoothers();
  void fix_MainDiagonal();
  void ML_Nox_computeFineLevelJacobian();

  void ML_Nox_ApplyInverse_Linear() const;
  void ML_Nox_ApplyInverse_NonLinear() const;

  //the FAS-preconditioner
  void ML_Nox_FAS_cycle() const;
  // the FAS-solver
  void ML_Nox_FAS_cycle1() const;

  void Ml_Nox_adaptivesetup();

private:

  char* label_;
  bool  isinit_; 
  int   ncalls_NewPrec_;

  const Epetra_Comm&                              comm_;
  RefCountPtr<NLNML::NLNML_FineLevelNoxInterface> interface_;
  RefCountPtr<Teuchos::ParameterList>             params_;
  RefCountPtr<NOX::Solver::Manager>               noxsolver_;

}; // class NLNML_Preconditioner
} // namespace NLNML

#endif
#endif
