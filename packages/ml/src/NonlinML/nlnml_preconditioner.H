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
 * \date Last update do Doxygen: 6-April-06
 *
 */
#include "ml_common.h"
#include "ml_include.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef NLNML_PRECONDITIONER_H
#define NLNML_PRECONDITIONER_H

// ML-headers
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
#include "Epetra_IntVector.h"
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
#include "NOX_Epetra_FiniteDifferenceColoring.H"
#include "NOX_Parameter_List.H"
#include "NOX_Epetra_Vector.H"

// own headers
#include "nlnml_preconditioner_utils.H"
#include "nlnml_finelevelnoxinterface.H"
#include "nlnml_nonlinearlevel.H"

// ---------- Forward Declarations ----------
class NLNML::NLNML_FineLevelNoxInterface;

using namespace std;
using namespace Teuchos;
//! NLNML: default namespace for all ML/NOX nonlinear methods
namespace NLNML
{

/*!
 
   \brief NLNML_Preconditioner: a class to define a multilevel nonlinear preconditioner and solver

   Class NLNML::NLNML_Preconditioner defines an algebraic nonlinear preconditioner
   or solver that is based on a full approximation scheme (FAS) approach and
   is implemented making use of NOX' nonlinear solver capabilities.
   It operates on Epetra objects and implements an Epetra_Operator to be easily
   used as a preconditioner to some outer nonlinear iteration.
   It also implements the NOX::EpetraNew::Interface::Preconditioner from NOX
   so it can serve as it's own interface (see NOX documentation for details).
   In this class several methods and options can be chosen, which are described in the
   parameter set functions
   
  
    To make use of the full functionality of this class it 
    requires ML to be configured with the following options:
    - \c --with-ml_nox
    - \c --enable-epetra
    - \c --enable-epetraext
    - \c --enable-nox
    - \c --enable-nox-epetra
    - \c --enable-prerelease
    - \c --enable-aztecoo
    - \c --enable-amesos
    
    \note
    an example is given in the directory ml/examples/NonlinML

    \author Michael Gee, SNL 9214
*/  
class NLNML_Preconditioner : public virtual Epetra_Operator,
                             public virtual NOX::Epetra::Interface::Preconditioner
{
public:

  explicit NLNML_Preconditioner(
                     RefCountPtr<NLNML::NLNML_FineLevelNoxInterface> interface,
                     ParameterList& mlparams,
                     const Epetra_Comm& comm);
  
  virtual ~NLNML_Preconditioner();
  
  void SetNoxSolver(RefCountPtr<NOX::Solver::Manager> s) { noxsolver_ = s; } 
  
  const char* Label() const {return(&(label_[0]));}
  
  const Epetra_Comm & Comm() const {return(comm_);}
  
  const Epetra_Map & OperatorDomainMap() const {return interface_->getMap();}
  
  const Epetra_Map & OperatorRangeMap() const {return interface_->getMap();}
  
  inline int OutLevel() const { return params_->get("nlnML output",0); } 
  
  inline bool isinit() const {return(isinit_);}
  
  bool setinit(bool flag) 
  {isinit_=flag; ncalls_NewPrec_ = 0; return(true);}
  
  int getoffset() const {return params_->get("nlnML offset recompute",0);} 
  
  virtual bool computePreconditioner(const Epetra_Vector& x, 
				     Epetra_Operator& M,
				     NOX::Parameter::List* precParams = 0);
  
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const ;
  
  int solve() const ;
  
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
  bool Compute_Linearpreconditioner(const Epetra_Vector& x);      
  bool Compute_Nonlinearpreconditioner(const Epetra_Vector& x);   

  void Set_Smoothers();
  void fix_MainDiagonal(RefCountPtr<Epetra_CrsMatrix> A, int level);
  RefCountPtr<Epetra_CrsMatrix> ComputeFineLevelJacobian(const Epetra_Vector& x);

  int ApplyInverse_Linear(const Epetra_MultiVector& X, 
                          Epetra_MultiVector& Y) const;
  int ApplyInverse_NonLinear(const Epetra_MultiVector& X, 
                             Epetra_MultiVector& Y) const;

  //the FAS-preconditioner
  bool FAS_cycle(Epetra_Vector* f, Epetra_Vector* x, int level, bool* converged, double* finenorm) const;
  // the FAS-solver
  bool FAS_cycle_solver() const;

  bool Adaptivesetup(double** oldns,Epetra_CrsMatrix* Jac,int oldnumpde,int olddimns);
  
  template<typename kind>
  kind getParameter(string name,kind dvalue) const
  { return params_->get(name,dvalue); }

  template<typename kind>
  void setParameter(string name,kind dvalue) const
  { params_->set(name,dvalue); return; }

  bool CheckInputParameters(ParameterList& params);
  
private:

  char* label_;
  bool  isinit_; 
  int   ncalls_NewPrec_;
  
  const Epetra_Comm&                                              comm_;
  RefCountPtr<NLNML::NLNML_FineLevelNoxInterface>                 interface_;
  RefCountPtr<Teuchos::ParameterList>                             params_;
  RefCountPtr<NOX::Solver::Manager>                               noxsolver_;
  RefCountPtr<Epetra_CrsMatrix>                                   fineJac_;
  RefCountPtr<Epetra_Map>                                         rowmap_;
  RefCountPtr<ML_Epetra::MultiLevelOperator>                      linPrec_;
  RefCountPtr<vector<RefCountPtr<NLNML::NLNML_NonlinearLevel> > > nlnlevel_;
  ML*                                                             ml_;
  ML_Aggregate*                                                   ag_;


}; // class NLNML_Preconditioner
} // namespace NLNML

#endif
#endif
