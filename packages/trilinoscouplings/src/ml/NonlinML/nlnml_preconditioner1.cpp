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
interface_(interface),
ml_(NULL),
ag_(NULL),
label_("nlnML_Preconditioner")
{
  CheckInputParameters(mlparams);
  params_   = rcp(new Teuchos::ParameterList(mlparams));
  // we make a backup of the nullspace dimension as it might be
  // overwritten by the adaptive ns procedure
  int dimns = getParameter("nlnML null space: dimension",1);
              setParameter("nlnML null space: dimension2",dimns);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_Preconditioner::~NLNML_Preconditioner()
{
  if (ag_) ML_Aggregate_Destroy(&ag_);
  if (ml_) ML_Destroy(&ml_);
  return;
}


/*----------------------------------------------------------------------*
 |   (private)                                                m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_Preconditioner::CheckInputParameters(ParameterList& params)
{
  int missing = 0;
  vector<string> missingparams(38);
  
  if (!params.isParameter("nlnML output")) 
    { missingparams[missing] = "\"nlnML output\" [int]"; ++missing; }
  if (!params.isParameter("nlnML max levels")) 
    { missingparams[missing] = "\"nlnML max levels\" [int]"; ++missing; }
  if (!params.isParameter("nlnML coarse: max size")) 
    { missingparams[missing] = "\"nlnML coarse: max size\" [int]"; ++missing; }
  if (!params.isParameter("nlnML is linear preconditioner")) 
    { missingparams[missing] = "\"nlnML is linear preconditioner\" [bool]"; ++missing; }
  if (!params.isParameter("nlnML apply constraints")) 
    { missingparams[missing] = "\"nlnML apply constraints\" [bool]"; ++missing; }
  if (!params.isParameter("nlnML is matrixfree")) 
    { missingparams[missing] = "\"nlnML is matrixfree\" [bool]"; ++missing; }
  if (!params.isParameter("nlnML finite difference fine level")) 
    { missingparams[missing] = "\"nlnML finite difference fine level\" [bool]"; ++missing; }
  if (!params.isParameter("nlnML finite difference alpha")) 
    { missingparams[missing] = "\"nlnML finite difference alpha\" [double]"; ++missing; }
  if (!params.isParameter("nlnML finite difference beta")) 
    { missingparams[missing] = "\"nlnML finite difference beta\" [double]"; ++missing; }
  if (!params.isParameter("nlnML finite difference centered")) 
    { missingparams[missing] = "\"nlnML finite difference centered\" [bool]"; ++missing; }
  if (!params.isParameter("nlnML Jacobian fix diagonal")) 
    { missingparams[missing] = "\"nlnML Jacobian fix diagonal\" [bool]"; ++missing; }
  if (!params.isParameter("nlnML absolute residual tolerance")) 
    { missingparams[missing] = "\"nlnML absolute residual tolerance\" [double]"; ++missing; }
  if (!params.isParameter("nlnML max cycles")) 
    { missingparams[missing] = "\"nlnML max cycles\" [int]"; ++missing; }
  if (!params.isParameter("nlnML adaptive recompute")) 
    { missingparams[missing] = "\"nlnML adaptive recompute\" [double]"; ++missing; }
  if (!params.isParameter("nlnML offset recompute")) 
    { missingparams[missing] = "\"nlnML offset recompute\" [int]"; ++missing; }
  if (!params.isParameter("nlnML additional adaptive nullspace")) 
    { missingparams[missing] = "\"nlnML additional adaptive nullspace\" [int]"; ++missing; }
  if (!params.isParameter("nlnML PDE equations")) 
    { missingparams[missing] = "\"nlnML PDE equations\" [int]"; ++missing; }
  if (!params.isParameter("nlnML null space: dimension")) 
    { missingparams[missing] = "\"nlnML null space: dimension\" [int]"; ++missing; }
  if (!params.isParameter("nlnML spatial dimension")) 
    { missingparams[missing] = "\"nlnML spatial dimension\" [int]"; ++missing; }
  if (!params.isParameter("nlnML coarse: type")) 
    { missingparams[missing] = "\"nlnML coarse: type\" [string]"; ++missing; }
  if (!params.isParameter("nlnML nodes per aggregate")) 
    { missingparams[missing] = "\"nlnML nodes per aggregate\" [int]"; ++missing; }
  if (!params.isParameter("nlnML use nlncg on fine level")) 
    { missingparams[missing] = "\"nlnML use nlncg on fine level\" [bool]"; ++missing; }
  if (!params.isParameter("nlnML use nlncg on medium level")) 
    { missingparams[missing] = "\"nlnML use nlncg on medium level\" [bool]"; ++missing; }
  if (!params.isParameter("nlnML use nlncg on coarsest level")) 
    { missingparams[missing] = "\"nlnML use nlncg on coarsest level\" [bool]"; ++missing; }
  if (!params.isParameter("nlnML max iterations newton-krylov fine level")) 
    { missingparams[missing] = "\"nlnML max iterations newton-krylov fine level\" [int]"; ++missing; }
  if (!params.isParameter("nlnML max iterations newton-krylov medium level")) 
    { missingparams[missing] = "\"nlnML max iterations newton-krylov medium level\" [int]"; ++missing; }
  if (!params.isParameter("nlnML max iterations newton-krylov coarsest level")) 
    { missingparams[missing] = "\"nlnML max iterations newton-krylov coarsest level\" [int]"; ++missing; }
  if (!params.isParameter("nlnML linear smoother type fine level")) 
    { missingparams[missing] = "\"nlnML linear smoother type fine level\" [string]"; ++missing; }
  if (!params.isParameter("nlnML linear smoother type medium level")) 
    { missingparams[missing] = "\"nlnML linear smoother type medium level\" [string]"; ++missing; }
  if (!params.isParameter("nlnML linear smoother type coarsest level")) 
    { missingparams[missing] = "\"nlnML linear smoother type coarsest level\" [string]"; ++missing; }
  if (!params.isParameter("nlnML linear smoother sweeps fine level")) 
    { missingparams[missing] = "\"nlnML linear smoother sweeps fine level\" [int]"; ++missing; }
  if (!params.isParameter("nlnML linear smoother sweeps medium level")) 
    { missingparams[missing] = "\"nlnML linear smoother sweeps medium level\" [int]"; ++missing; }
  if (!params.isParameter("nlnML linear smoother sweeps coarsest level")) 
    { missingparams[missing] = "\"nlnML linear smoother sweeps coarsest level\" [int]"; ++missing; }
  if (!params.isParameter("nlnML nonlinear presmoothing sweeps fine level")) 
    { missingparams[missing] = "\"nlnML nonlinear presmoothing sweeps fine level\" [int]"; ++missing; }
  if (!params.isParameter("nlnML nonlinear presmoothing sweeps medium level")) 
    { missingparams[missing] = "\"nlnML nonlinear presmoothing sweeps medium level\" [int]"; ++missing; }
  if (!params.isParameter("nlnML nonlinear smoothing sweeps coarse level")) 
    { missingparams[missing] = "\"nlnML nonlinear smoothing sweeps coarse level\" [int]"; ++missing; }
  if (!params.isParameter("nlnML nonlinear postsmoothing sweeps medium level")) 
    { missingparams[missing] = "\"nlnML nonlinear postsmoothing sweeps medium level\" [int]"; ++missing; }
  if (!params.isParameter("nlnML nonlinear postsmoothing sweeps fine level")) 
    { missingparams[missing] = "\"nlnML nonlinear postsmoothing sweeps fine level\" [int]"; ++missing; }
  // 38 parameters right now
  
  
  if (missing && Comm().MyPID()==0)
  {
    cout << "======================================================\n";
    cout << "nlnML (level 0): missing parameters in parameter list:\n";
    cout << "------------------------------------------------------\n";
    for (int i=0; i<missing; ++i)
      cout << missingparams[i] << endl;
    cout << "------------------------------------------------------\n";
    cout << "nlnML (level 0): Note that depending on what method you would like\n"
         << "                 to run it's ok to miss some while for others\n"
         << "                 suboptimal default values will be used.\n";
    cout << "======================================================\n";
    fflush(stdout);     
  }
  
  return true;
}


/*----------------------------------------------------------------------*
 |  compute this preconditioner (public, derived)             m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_Preconditioner::computePreconditioner(
                                     const Epetra_Vector& x, 
				     Epetra_Operator& M,
				     Teuchos::ParameterList* precParams)
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
  if ( getParameter("nlnML finite difference fine level",false)==true)
    setParameter("nlnML is matrixfree",true);
    
  // Jacobian is supposed to be supplied by application
  if (getParameter("nlnML is matrixfree",false)==false)
  {
    // recompute the Jacobian
    fineJac_ = rcp(interface_->getJacobian());
    fineJac_.release();
    interface_->computeJacobian(x,*fineJac_);
    fineJac_ = rcp(interface_->getJacobian());
    fineJac_.release();
    // check for zero diagonal entries
    if (getParameter("nlnML Jacobian fix diagonal",true)==true)
      fix_MainDiagonal(fineJac_,0);
      
    fineJac_ = NLNML::StripZeros(fineJac_,1.0e-11);
  }
  else if (getParameter("nlnML is matrixfree",false)==true &&
           getParameter("nlnML finite difference fine level",false)==true)
  {
    fineJac_ = ComputeFineLevelJacobian(x);
    setParameter("nlnML is matrixfree",false);
    if (getParameter("nlnML Jacobian fix diagonal",true)==true)
      fix_MainDiagonal(fineJac_,0);
  }
  else
  {
    cout << "**ERR**: NLNML::NLNML_Preconditioner::compPrec:\n"
         << "**ERR**: unknown parameter combination for Jacobian\n"  
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  
  // create ml
  if (ag_) ML_Aggregate_Destroy(&ag_);
  if (ml_) ML_Destroy(&ml_);
  ML_Aggregate_Create(&ag_);
  int maxnlevel = getParameter("nlnML max levels",10);
  ML_Create(&ml_,maxnlevel);
  
  // set fine grid matrix
  EpetraMatrix2MLMatrix(ml_,0,(dynamic_cast<Epetra_RowMatrix*>(fineJac_.get())));
  
  // set coarsening type
  string coarsentype = getParameter("nlnML coarse: type",(string)"Uncoupled");
  if (coarsentype=="Uncoupled")
    ML_Aggregate_Set_CoarsenScheme_Uncoupled(ag_);
  else if (coarsentype=="MIS")
    ML_Aggregate_Set_CoarsenScheme_MIS(ag_);
  else if (coarsentype=="METIS")
  {
    ML_Aggregate_Set_CoarsenScheme_METIS(ag_);
    int nnodeperagg = getParameter("nlnML nodes per aggregate",27);
    for (int i=0; i<maxnlevel; ++i)
      ML_Aggregate_Set_NodesPerAggr(ml_,ag_,i,nnodeperagg);
  }
  else if (coarsentype=="VBMETIS")
  {
    vector<int> blocks(0);
    vector<int> block_pde(0);
    int         nblocks = 0;
    bool ok = interface_->getBlockInfo(&nblocks,blocks,block_pde);
    if (!ok)
    {
      cout << "**WRN**: NLNML::NLNML_Preconditioner::compPrec:\n"
           << "**WRN**: interface returned no blocks,\n"
           << "**WRN**: using aggregation scheme METIS instead of VBMETIS\n";
      setParameter("nlnML coarse: type",(string)"METIS");
      ML_Aggregate_Set_CoarsenScheme_METIS(ag_);
    }
    else
    {
      ML_Aggregate_Set_CoarsenScheme_VBMETIS(ag_);
      int nupdate = interface_->getMap().NumMyElements();
      ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS(ag_,0,maxnlevel,
                                                     nblocks,&blocks[0],
                                                     &block_pde[0],nupdate);
      blocks.clear();
      block_pde.clear();
    }
    int nnodeperagg = getParameter("nlnML nodes per aggregate",27);
    for (int i=0; i<maxnlevel; ++i)
      ML_Aggregate_Set_NodesPerAggr(ml_,ag_,i,nnodeperagg);
  }
  
  // set damping factor and threshold for smoothed aggregation
  ML_Aggregate_Set_DampingFactor(ag_, 1.3333);
  ML_Aggregate_Set_Threshold(ag_, 0.0);
  
  // set max coarse grid size
  ML_Aggregate_Set_MaxCoarseSize(ag_,getParameter("nlnML coarse: max size",128));
  
  // Calculate spectral norm
  ML_Set_SpectralNormScheme_Calc(ml_);
  
  // Set output level
  ML_Set_PrintLevel(OutLevel());
  
  // get nullspace
  int nummyrows = fineJac_->NumMyRows();
  int dimnullsp = getParameter("nlnML null space: dimension",1);
  int blocksize = getParameter("nlnML PDE equations",1);
  double* nullsp = interface_->Get_Nullspace(nummyrows,blocksize,dimnullsp);
  
  // run adaptive nullspace procedure
  int adaptns = getParameter("nlnML additional adaptive nullspace",0);
  if (adaptns)
  {
    Adaptivesetup(&nullsp,fineJac_.get(),blocksize,dimnullsp);
    setParameter("nlnML null space: dimension",(dimnullsp+adaptns));
  }
  
  // Pass nullspace to ml
  if (nullsp)
  {
    dimnullsp = getParameter("nlnML null space: dimension",1);
    ML_Aggregate_Set_NullSpace(ag_,blocksize,dimnullsp,nullsp,nummyrows);
    delete [] nullsp; nullsp = NULL;
  }
  else
  {
     dimnullsp = getParameter("nlnML null space: dimension",1);
     cout << "**WRN**: NLNML_Preconditioner::compPrec:\n"
          << "**WRN**: interface returned no nullspace,\n"
          << "**WRN**: using ML's default nullspace\n";
     if (blocksize != dimnullsp)
     {
       cout << "**WRN**: ML_Nox_Preconditioner::compPrec:\n"
            << "**WRN**: with default nullspace, nullspace-dimension must match number of PDEs per node\n"
            << "**WRN**: numPDE = " << blocksize << ", dimnullsp = " << dimnullsp << "\n"
            << "**WRN**: continue with setting dimnullsp = numPDE\n";
       setParameter("nlnML null space: dimension",blocksize);
       setParameter("nlnML null space: dimension2",blocksize);
     }
     ML_Aggregate_Set_NullSpace(ag_,blocksize,dimnullsp,NULL,nummyrows);
  }
  
  // keep the aggregation info
  ag_->keep_agg_information = 1;
  
  // build amg hierarchy
  int nlevel = ML_Gen_MGHierarchy_UsingAggregation(ml_,0,ML_INCREASING,ag_);
  setParameter("nlnML max levels",nlevel);
  
  if (nlevel<2)
  {
     cout << "**ERR**: NLNML_Preconditioner::compPrec:\n"
          << "**ERR**: number of levels generated is " << nlevel << "\n"
          << "**ERR**: this algorithm relies on at least nlevel >=2 !\n"
          << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; throw -1;
  }
  
  // do rest of setup
  bool islinearPrec = getParameter("nlnML is linear preconditioner",true);
  if (islinearPrec)
  {
    isinit_ = Compute_Linearpreconditioner(x);
    if (isinit_) return true;
    else         return false;
  }
  else
  {
    isinit_ = Compute_Nonlinearpreconditioner(x);
    if (isinit_) return true;
    else         return false;
  }
}



/*----------------------------------------------------------------------*
 |  compute linear preconditioner                             m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_Preconditioner::Compute_Linearpreconditioner(
                                                  const Epetra_Vector& x)
{
  // choose some smoothers
  Set_Smoothers();
  rowmap_  = rcp(new Epetra_Map(interface_->getMap()));
  linPrec_ = rcp(new ML_Epetra::MultiLevelOperator(ml_,Comm(),*rowmap_,*rowmap_));
  return true;
}



/*----------------------------------------------------------------------*
 |  compute nonlinear preconditioner                          m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_Preconditioner::Compute_Nonlinearpreconditioner(
                                                  const Epetra_Vector& x)
{
  // extract the prolongations from the hierachy
  int maxlevel = getParameter("nlnML max levels",10);
  
  RefCountPtr< vector< RefCountPtr<Epetra_CrsMatrix> > > P = 
    rcp( new vector< RefCountPtr<Epetra_CrsMatrix> >(maxlevel) );
    
  vector< RefCountPtr<Epetra_CrsMatrix> >& Pref = *P;
  for (int i=0; i<maxlevel; ++i) Pref[i] = null;
  for (int i=1; i<maxlevel; ++i) // there is not Pmat on level 0
  {
    double t1     = GetClock();
    int    maxnnz = 0;
    double cputime;
    Epetra_CrsMatrix* tmp;
    ML_Operator2EpetraCrsMatrix(&(ml_->Pmat[i]), tmp, maxnnz, 
                                false, cputime);
    (*P)[i] = rcp(tmp);
    (*P)[i]->OptimizeStorage();
    double t2 = GetClock() - t1;
    if (OutLevel()>5 && Comm().MyPID()==0)
      cout << "nlnML (level " << i << "): extraction of P in " << t2 << " sec\n";
  }
  
  // create the vector of nonlinear levels
  nlnlevel_ = rcp( new vector<RefCountPtr<NLNML::NLNML_NonlinearLevel> >(maxlevel));
  
  // loop all levels and allocate the nonlinear level class  
  for (int i=0; i<maxlevel; ++i)
  {
    int spatial = getParameter("nlnML spatial dimension",1);
    int dimns   = getParameter("nlnML null space: dimension",1);
    int numpde  = 0;
    if (i==0) numpde = getParameter("nlnML PDE equations",1);
    else if (spatial==1)
      numpde=1;
    else if (spatial==2)
      numpde=3;
    else if (spatial==3)
      numpde=6;
    else
    {
      cout << "**ERR**: NLNML::NLNML_Preconditioner::Compute_Nonlinearpreconditioner:\n"
           << "**ERR**: nlnML spatial dimension = " << spatial << " unknown\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
    }
    
    // choose the nonlinear solver
    bool isnlncg  = true;
    int  niterscg = 0;
    if (i==0) // fine
    {
      isnlncg  = getParameter("nlnML use nlncg on fine level",true);
      niterscg = getParameter("nlnML max iterations newton-krylov fine level",40);
    }
    else if (i==maxlevel-1) // coarse
    {
      isnlncg  = getParameter("nlnML use nlncg on coarsest level",true);
      niterscg = getParameter("nlnML max iterations newton-krylov coarsest level",40);
    }
    else // intermediate
    {
      isnlncg  = getParameter("nlnML use nlncg on medium level",true);
      niterscg = getParameter("nlnML max iterations newton-krylov medium level",40);
    }
    
    // Allocate the nonlinear level class
    (*nlnlevel_)[i] = rcp( new NLNML::NLNML_NonlinearLevel(i,params_,ml_,ag_,
                                                           P,interface_,Comm(),
                                                           fineJac_,
                                                           x,isnlncg,niterscg,
                                                           numpde,dimns));
  } // for (int i=0; i<maxlevel; ++i)
  
  // don't need the ag_ and ml_ objects anymore
  if (ag_)
  {
     ML_Aggregate_Destroy(&ag_);
     ag_ = 0;
  }
  if (ml_)
  {
     ML_Destroy(&ml_);
     ml_ = 0;
  }
  return true;
}



/*----------------------------------------------------------------------*
 |  apply this preconditioner (public, derived)               m.gee 3/06|
 *----------------------------------------------------------------------*/
int NLNML::NLNML_Preconditioner::ApplyInverse(
                                     const Epetra_MultiVector& X, 
                                     Epetra_MultiVector& Y) const
{
  if (!isinit())
  {
    Epetra_Vector x(View,X,0);
    NLNML::NLNML_Preconditioner* tmp = 
                          const_cast<NLNML::NLNML_Preconditioner*>(this);
    tmp->computePreconditioner(x,*tmp);
  } 
  
  double t0 = GetClock();
  int ncalls0 = interface_->getnumcallscomputeF();
  
  int err = 0;
  if (getParameter("nlnML is linear preconditioner",true)==true)
    err = ApplyInverse_Linear(X,Y);
  else
    err = ApplyInverse_NonLinear(X,Y);
  
  // let application enforce constraints on the gradient
  Epetra_Vector epetragradient(View,Y,0);
  interface_->ApplyAllConstraints(epetragradient,0);
  
  int ncalls1 = interface_->getnumcallscomputeF();
  double t1 = GetClock();
  
  if (OutLevel()>7 && Comm().MyPID()==0)
    cout << "nlnML (level 0): Preconditioner time " << (t1-t0) 
         << " sec, # evaluateF total/this " 
         << ncalls1 << "/" << (ncalls1-ncalls0) << endl;
  
  if (!err) return 0;
  else      return -1;
}



/*----------------------------------------------------------------------*
 |  apply inverse for linear preconditioner                   m.gee 3/06|
 *----------------------------------------------------------------------*/
int NLNML::NLNML_Preconditioner::ApplyInverse_Linear(const Epetra_MultiVector& X, 
                                                     Epetra_MultiVector& Y) const
{
  if (linPrec_==null)
  {
    cout << "**ERR**: NLNML_Preconditioner::ApplyInverse_Linear:\n";
    cout << "**ERR**: ptr to ML_Epetra::MultiLevelOperator is NULL\n";
    cout << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; throw -1;
  }
  
  // apply ML
  return linPrec_->ApplyInverse(X,Y);   
}



/*----------------------------------------------------------------------*
 |  apply inverse for nonlinear preconditioner                   m.gee 3/06|
 *----------------------------------------------------------------------*/
int NLNML::NLNML_Preconditioner::ApplyInverse_NonLinear(const Epetra_MultiVector& X, 
                                                        Epetra_MultiVector& Y) const
{
  if (noxsolver_==null)
  {
    cout << "**ERR**: NLNML::NLNML_Preconditioner::ApplyInverse_NonLinear:\n"
         << "**ERR**: noxsolver not registered, use SetNoxSolver(solver)!\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  
  const NOX::Epetra::Group& finalGroup = 
    dynamic_cast<const NOX::Epetra::Group&>(noxsolver_->getSolutionGroup());
  const Epetra_Vector& currentSolution = 
   (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();      
  const Epetra_Vector& currentF = 
   (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getF())).getEpetraVector();
  double norm2;
   currentF.Norm2(&norm2);
  
  // make a copy of currentSolution and currentF
  RefCountPtr<Epetra_Vector> f = rcp(new Epetra_Vector(View,X,0));
  RefCountPtr<Epetra_Vector> x = rcp(new Epetra_Vector(Copy,currentSolution,0));

  // call the cycle
  if (OutLevel() && !Comm().MyPID())
    cout << "\n\nnlnML :============Entering Nonlinear V-cycle============\n";
  bool converged = false;
  double t3 = GetClock();
  FAS_cycle(f.get(),x.get(),0,&converged,&norm2); 
  double t4 = GetClock();
  if (OutLevel() && !Comm().MyPID())
  {
    cout << "nlnML :============V-cycle time is : " << (t4-t3) << " sec\n";
    if (converged)
      cout << "nlnML :============Nonlinear preconditioner converged====\n";
  }
  f = null;  
  Y.Scale(1.0,*x);
  return 0;
}




/*----------------------------------------------------------------------*
 |  run FAS preconditioner as a solver                        m.gee 3/06|
 *----------------------------------------------------------------------*/
int NLNML::NLNML_Preconditioner::solve() const
{
  // get starting solution form the interface
  const Epetra_Vector* xfine = interface_->getSolution();
  
  // compute preconditioner
  double t5 = GetClock();
  NLNML::NLNML_Preconditioner* tmp = 
                          const_cast<NLNML::NLNML_Preconditioner*>(this);
  tmp->computePreconditioner(*xfine,*tmp);
  double t6 = GetClock();
  
  // check sanity
  if (!isinit())
  {
    cout << "**ERR**: NLNML::NLNML_Preconditioner::solve:\n"
         << "**ERR**: initflag is false\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  if (getParameter("nlnML is linear preconditioner",true) != false)
  {
    cout << "**ERR**: NLNML::NLNML_Preconditioner::solve:\n"
         << "**ERR**: Preconditioner has to be nonlinear preconditioner to run as solver\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  
  // compute current residual
  RefCountPtr<Epetra_Vector> f = rcp(new Epetra_Vector(xfine->Map(),false));
  RefCountPtr<Epetra_Vector> x = rcp(new Epetra_Vector(Copy,*xfine,0));
  (*nlnlevel_)[0]->setModifiedSystem(false,NULL,NULL);
  (*nlnlevel_)[0]->computeF(*x,*f,NOX::Epetra::Interface::Required::Residual);
  
  // max cycles
  int maxcycle = getParameter("nlnML max cycles",10);
  
  {
    double t1 = GetClock();
    bool   converged = false;
    for (int i=1; i<=maxcycle; ++i)
    {
      if (OutLevel() && Comm().MyPID()==0)
      {
        cout << "\n\nnlnML (level 0):============ nonlinear V-cycle # " << i << " ================\n";
        fflush(stdout);
      }
      double t3 = GetClock();
      double norm = getParameter("nlnML absolute residual tolerance",1.0e-06);
      
      //===== presmooth level 0=================================================
      int nsmooth = getParameter("nlnML nonlinear presmoothing sweeps fine level",2);
      double prenorm = norm;
      if (nsmooth)
        converged = (*nlnlevel_)[0]->Iterate(f.get(),x.get(),nsmooth,&prenorm);
      if (converged) break;
      
      //===== restrict to level 0===============================================
      RefCountPtr<Epetra_Vector> xcoarse = rcp((*nlnlevel_)[0]->restrict_to_next_coarser_level(x.get(),0,1));
      RefCountPtr<Epetra_Vector> fcoarse = rcp((*nlnlevel_)[0]->restrict_to_next_coarser_level(f.get(),0,1));
      
      //===== call FAS-cycle on level 1=========================================
      bool coarseconverged=false;
      FAS_cycle(fcoarse.get(),xcoarse.get(),1,&coarseconverged,&prenorm);
      fcoarse = null;
      
      //===== prolongate correction=============================================
      RefCountPtr<Epetra_Vector> xcorrect = rcp((*nlnlevel_)[0]->prolong_to_this_level(xcoarse.get(),0,1));
      xcoarse = null;
      
      //=====apply correction===================================================
      x->Update(1.0,*xcorrect,1.0);
      xcorrect = null;
      
      //=====postsmoothing level 0==============================================
      nsmooth = getParameter("nlnML nonlinear postsmoothing sweeps fine level",2);
      double postnorm = norm;
      if (nsmooth)
        converged = (*nlnlevel_)[0]->Iterate(f.get(),x.get(),nsmooth,&postnorm);
      if (converged) break;
      if (OutLevel() && Comm().MyPID()==0)
      {
        double t4 = GetClock();
        cout << "nlnML (level 0):============ time this cycle: " << t4-t3 << " sec\n";
        fflush(stdout);
      }
    }
    double t2 = GetClock();
    if (OutLevel() && Comm().MyPID()==0)
    {
      cout << "nlnML (level 0):============solve time incl. setup " 
           << t2-t1+t6-t5 << " sec\n";
      fflush(stdout);
    }
    if (converged) return 0;
    else           return -1;
  }
}




/*----------------------------------------------------------------------*
 |  choose smoothers                                          m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_Preconditioner::Set_Smoothers()
{
  // choose fine grid smoother
  int    maxnlevel   = getParameter("nlnML max levels",2);
  string fsmoother   = getParameter("nlnML linear smoother type fine level",(string)"SGS");
  int    nsmoothfine = getParameter("nlnML linear smoother sweeps fine level",1);
  string smoother    = getParameter("nlnML linear smoother type medium level",(string)"SGS");
  int    nsmooth     = getParameter("nlnML linear smoother sweeps medium level",1);
  string csmoother   = getParameter("nlnML linear smoother type coarsest level",(string)"AmesosKLU");
  int    ncsmooth    = getParameter("nlnML linear smoother sweeps coarsest level",1);
  int    coarsegrid  = maxnlevel - 1;  
  
  // choose fine level smoother
  if (fsmoother=="SGS")
    ML_Gen_Smoother_SymGaussSeidel(ml_,0,ML_BOTH,nsmoothfine,0.67);
  else if (fsmoother=="Jacobi")
  {
    ML_Gen_Smoother_Jacobi(ml_,0,ML_PRESMOOTHER, nsmoothfine,0.2);
    ML_Gen_Smoother_Jacobi(ml_,0,ML_POSTSMOOTHER,nsmoothfine,0.2);
  }
  else if (fsmoother=="BSGS")
  {
    int  nblocks  = 0;
    int* blocks   = NULL;
    int* blockpde = NULL;
    // try to get nodal blocks from the VBMETIS aggregation scheme
    ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS(ag_,0,maxnlevel,
                                                   &nblocks,&blocks,&blockpde);
    if (!nblocks && !blocks)
       ML_Gen_Blocks_Aggregates(ag_,0,&nblocks,&blocks);
       
    ML_Gen_Smoother_VBlockSymGaussSeidel(ml_,0,ML_BOTH,nsmoothfine,0.67,
                                         nblocks,blocks);
    if (nblocks && blocks)
    {
       ML_free(blocks); 
       ML_free(blockpde);
    }
  }
  else if (fsmoother=="Bcheby")
  {
    int  nblocks  = 0;
    int* blocks   = NULL;
    int* blockpde = NULL;
    // try to get nodal blocks from the VBMETIS aggregation scheme
    ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS(ag_,0,maxnlevel,
                                                   &nblocks,&blocks,&blockpde);
    if (!nblocks && !blocks)
       ML_Gen_Blocks_Aggregates(ag_,0,&nblocks,&blocks);
       
     ML_Gen_Smoother_BlockDiagScaledCheby(ml_,0,ML_BOTH,30.,nsmoothfine,
                                          nblocks,blocks);
    if (nblocks && blocks)
    {
       ML_free(blocks); 
       ML_free(blockpde);
    }
  }
  else if ((fsmoother=="MLS")||(fsmoother=="Cheby"))
    ML_Gen_Smoother_Cheby(ml_,0,ML_BOTH,30.,nsmoothfine);
  else if (fsmoother=="AmesosKLU")
    ML_Gen_Smoother_Amesos(ml_,0,ML_AMESOS_KLU,-1,0.0);
  else
  {
    cout << "**ERR**: NLNML_Preconditioner::Set_Smoothers:\n"
         << "**ERR**: smoother " << fsmoother << " not recognized\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  
  // for intermediate levels)
  for (int i=1; i<coarsegrid; ++i)
  {
    if (smoother=="SGS")
      ML_Gen_Smoother_SymGaussSeidel(ml_,i,ML_BOTH,nsmooth,0.67);
    else if (smoother=="Jacobi")
    {
      ML_Gen_Smoother_Jacobi(ml_,i,ML_PRESMOOTHER, nsmooth,.2);
      ML_Gen_Smoother_Jacobi(ml_,i,ML_POSTSMOOTHER,nsmooth,.2);
    }
    else if (smoother=="BSGS")
    {
      int  nblocks  = 0;
      int* blocks   = NULL;
      int* blockpde = NULL;
      // try to get nodal blocks from the VBMETIS aggregation scheme
      ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS(ag_,i,maxnlevel,
                                                  &nblocks,&blocks,&blockpde);
      if (!nblocks && !blocks)
         ML_Gen_Blocks_Aggregates(ag_,i,&nblocks,&blocks);

      ML_Gen_Smoother_VBlockSymGaussSeidel(ml_,i,ML_BOTH,nsmooth,0.67,nblocks,blocks);
      if (nblocks && blocks)
      {
         ML_free(blocks); 
         ML_free(blockpde);
      }
    }
    else if (smoother=="Bcheby")
    {
      int  nblocks  = 0;
      int* blocks   = NULL;
      int* blockpde = NULL;
      // try to get nodal blocks from the VBMETIS aggregation scheme
      ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS(ag_,i,maxnlevel,
                                                  &nblocks,&blocks,&blockpde);
      if (!nblocks && !blocks)
         ML_Gen_Blocks_Aggregates(ag_,i,&nblocks,&blocks);

        ML_Gen_Smoother_BlockDiagScaledCheby(ml_,i,ML_BOTH,30.,nsmooth,
                                             nblocks,blocks);
      if (nblocks && blocks)
      {
         ML_free(blocks); 
         ML_free(blockpde);
      }
    }
    else if ((smoother=="MLS")||(smoother=="Cheby"))
      ML_Gen_Smoother_Cheby(ml_,i,ML_BOTH,30.,nsmooth);
    else if (smoother=="AmesosKLU")
      ML_Gen_Smoother_Amesos(ml_,i,ML_AMESOS_KLU,-1,0.0);
    else
    {
      cout << "**ERR**: NLNML_Preconditioner::Set_Smoothers:\n"
           << "**ERR**: smoother " << smoother << " not recognized\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
    }
  } // for (int i=1; i<coarsegrid; ++i)
  
  // choose coarse grid smoother
  if (csmoother=="SGS")
    ML_Gen_Smoother_SymGaussSeidel(ml_,coarsegrid,ML_BOTH,ncsmooth,0.67);
  else if (csmoother=="Jacobi")
  {
    ML_Gen_Smoother_Jacobi(ml_,coarsegrid,ML_PRESMOOTHER, ncsmooth,0.2);
    ML_Gen_Smoother_Jacobi(ml_,coarsegrid,ML_POSTSMOOTHER,ncsmooth,0.2);
  }
  else if (csmoother=="BSGS")
  {
    int  nblocks  = 0;
    int* blocks   = NULL;
    int* blockpde = NULL;
    // try to get nodal blocks from the VBMETIS aggregation scheme
    ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS(ag_,coarsegrid,maxnlevel,
                                                   &nblocks,&blocks,&blockpde);
    if (!nblocks && !blocks)
       ML_Gen_Blocks_Aggregates(ag_,0,&nblocks,&blocks);
       
    ML_Gen_Smoother_VBlockSymGaussSeidel(ml_,coarsegrid,ML_BOTH,ncsmooth,0.67,
                                         nblocks,blocks);
    if (nblocks && blocks)
    {
       ML_free(blocks); 
       ML_free(blockpde);
    }
  }
  else if (csmoother=="Bcheby")
  {
    int  nblocks  = 0;
    int* blocks   = NULL;
    int* blockpde = NULL;
    // try to get nodal blocks from the VBMETIS aggregation scheme
    ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS(ag_,coarsegrid,maxnlevel,
                                                   &nblocks,&blocks,&blockpde);
    if (!nblocks && !blocks)
       ML_Gen_Blocks_Aggregates(ag_,0,&nblocks,&blocks);
       
     ML_Gen_Smoother_BlockDiagScaledCheby(ml_,coarsegrid,ML_BOTH,30.,ncsmooth,
                                          nblocks,blocks);
    if (nblocks && blocks)
    {
       ML_free(blocks); 
       ML_free(blockpde);
    }
  }
  else if ((csmoother=="MLS")||(csmoother=="Cheby"))
    ML_Gen_Smoother_Cheby(ml_,coarsegrid,ML_BOTH,30.,ncsmooth);
  else if (csmoother=="AmesosKLU")
    ML_Gen_Smoother_Amesos(ml_,coarsegrid,ML_AMESOS_KLU,-1,0.0);
  else
  {
    cout << "**ERR**: NLNML_Preconditioner::Set_Smoothers:\n"
         << "**ERR**: smoother " << csmoother << " not recognized\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  return;
}



/*----------------------------------------------------------------------*
 |  fix main diagonal of Jacobian                             m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_Preconditioner::fix_MainDiagonal(
                                RefCountPtr<Epetra_CrsMatrix> A, int level)
{
  if (!A->Filled()) A->FillComplete();
  
  Epetra_Vector diag(A->RowMap(),false);
  int err = A->ExtractDiagonalCopy(diag);
  
  double average=0.0;
  for (int i=0; i<diag.MyLength(); i++) average += diag[i];
  average /= diag.MyLength();

  for (int i=0; i<diag.MyLength(); i++)
  {
     // fix zero entries to be real dirichlet boundary constraints because they come from mesh tying
     if (abs(diag[i])<1.0e-10)
     {
        if (OutLevel()>8)
          printf("found zero diagonal entry %20.12e in row %d, fixing to be %e\n",diag[i],i,average);
        //check whether there are nonzero off-diagonal entries in that row
        int numentries;
        double* values;
        err = A->ExtractMyRowView(i,numentries,values);
        if (err)
        {
           cout << "**ERR**: NLNML::NLNML_Preconditioner::fix_MainDiagonal:\n"
                << "**ERR**: A->ExtractMyRowView returned " << err << endl 
                << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
        }
        double sum=0.0;
        for (int j=0; j<numentries; j++)
           sum += abs(values[j]);
        
        if (sum>1.0e-9) 
          continue;
        
        //double small = 10.0;
        err = A->ReplaceMyValues(i,1,&average,&i);
        if (err)
        {
           cout << "**ERR**: NLNML::NLNML_Preconditioner::fix_MainDiagonal:\n"
                << "**ERR**: A->ReplaceMyValues returned " << err << endl 
                << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
        }
     }
     // fix small values to be resonably sized because they come from frictionless contact
     else if (abs(diag[i])<1.0)
     {
       double small=10.0;
       if (abs(diag[i])>small) small = abs(diag[i]);
       if (OutLevel()>8)
         printf("found tiny diagonal value %20.12e in row %d, fixing to be %e\n",diag[i],i,small);

       err = A->ReplaceMyValues(i,1,&small,&i);
       if (err)
       {
         cout << "**ERR**: NLNML::NLNML_Preconditioner::fix_MainDiagonal:\n"
              << "**ERR**: A->ReplaceMyValues returned " << err << endl 
              << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
       }
     }
     if (diag[i]<0.0)
     {
       cout << "ML: ***WRN*** Found negative main diag entry! Fixing...\n"; fflush(stdout);
       double small=10.0;
       err = A->ReplaceMyValues(i,1,&small,&i);
       if (err)
       {
         cout << "**ERR**: NLNML::NLNML_Preconditioner::fix_MainDiagonal:\n"
              << "**ERR**: A->ReplaceMyValues returned " << err << endl 
              << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
       }
     }
  }
  return;
}




/*----------------------------------------------------------------------*
 |  compute Jacobian on fine grid                             m.gee 3/06|
 *----------------------------------------------------------------------*/
RefCountPtr<Epetra_CrsMatrix> NLNML::NLNML_Preconditioner::ComputeFineLevelJacobian(const Epetra_Vector& x)
{
  // get graph and constrained modified graph
  Epetra_CrsGraph* graph  = 
                const_cast<Epetra_CrsGraph*>(interface_->getGraph());
  Epetra_CrsGraph* cgraph = 
                const_cast<Epetra_CrsGraph*>(interface_->getModifiedGraph());
  
  // get the block size of the problem
  int bsize = getParameter("nlnML PDE equations",1);
  
  double t0 = GetClock();

  if (OutLevel()>0 && Comm().MyPID()==0)
  {
     cout << "nlnML (level 0): Entering Coloring on level 0\n";
     fflush(stdout);
  }  
  
  RefCountPtr<Epetra_MapColoring> colorMap = 
                   NLNML::Collapsedcoloring(cgraph,bsize,false,OutLevel());
  if (colorMap==null)
    colorMap = NLNML::Standardcoloring(cgraph,false);
    
  RefCountPtr<EpetraExt::CrsGraph_MapColoringIndex> colorMapIndex = 
                  rcp(new EpetraExt::CrsGraph_MapColoringIndex(*colorMap));
  
  RefCountPtr< vector<Epetra_IntVector> > colorcolumns = 
                                           rcp(&(*colorMapIndex)(*graph));
  
  double t1 = GetClock();
  if (OutLevel()>0 && Comm().MyPID()==0)
  {
     cout << "nlnML (level 0): Proc " << comm_.MyPID() <<" Coloring time is " << (t1-t0) << " sec\n";
     cout << "nlnML (level 0): Entering Construction FD-Operator on level 0\n";
     fflush(stdout);
  }
  
  t0 = GetClock();
  int ncalls = interface_->getnumcallscomputeF();
  interface_->setnumcallscomputeF(0);

  RefCountPtr<Teuchos::ParameterList> dummylist = 
                                             rcp(new Teuchos::ParameterList());
  
  NOX::Epetra::Vector nx(x); 
  
  RefCountPtr<Epetra_CrsGraph> rcpgraph = rcp(graph); 
  rcpgraph.release();
  
  double alpha = getParameter("nlnML finite difference alpha",1.0e-07);
  double beta  = getParameter("nlnML finite difference beta",1.0e-06);
  
  RefCountPtr<NOX::Epetra::FiniteDifferenceColoring> FD = 
              rcp(new NOX::Epetra::FiniteDifferenceColoring(*dummylist,
                                                            interface_,
                                                            nx,
                                                            rcpgraph,
                                                            colorMap,
                                                            colorcolumns,
                                                            true,false,
                                                            beta,alpha));
  
  if (getParameter("nlnML finite difference centered",false)==true)
    FD->setDifferenceMethod(NOX::Epetra::FiniteDifferenceColoring::Centered);                                                            
  else                                                     
    FD->setDifferenceMethod(NOX::Epetra::FiniteDifferenceColoring::Forward);                                                            
  
  bool err = FD->computeJacobian(x); 
  if (err==false)
  {
    cout << "**ERR**: NLNML::NLNML_Preconditioner::ComputeFineLevelJacobian:\n"
         << "**ERR**: NOX::Epetra::FiniteDifferenceColoring returned an error on level 0" 
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  
  t1 = GetClock();
  if (OutLevel()>0 && Comm().MyPID()==0)
  {
     cout << "nlnML (level 0): colored Finite Differencing time :" << (t1-t0) << " sec\n";
     cout << "nlnML (level 0): colored Finite Differencing number of calls to computeF : " 
          << interface_->getnumcallscomputeF() << endl;
     fflush(stdout);
  }
  interface_->setnumcallscomputeF(ncalls);
  Epetra_CrsMatrix& tmp = FD->getUnderlyingMatrix();
  RefCountPtr<Epetra_CrsMatrix> rcptmp = rcp(&tmp); rcptmp.release();
  RefCountPtr<Epetra_CrsMatrix> A = NLNML::StripZeros(rcptmp,1.0e-11);
  return A;
}



/*----------------------------------------------------------------------*
 |  FAS-preconditioner                                        m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_Preconditioner::FAS_cycle(Epetra_Vector* f, 
                                            Epetra_Vector* x, 
                                            int level, 
                                            bool* converged, 
                                            double* finenorm) const
{
  int coarsegrid   = getParameter("nlnML max levels",2)-1;
  double thisnorm  = *finenorm/10000;
  double tolerance = getParameter("nlnML absolute residual tolerance",1.0e-06);
  if (thisnorm < tolerance) thisnorm = tolerance;
  
  RefCountPtr<Epetra_Vector> fbar  = null;
  RefCountPtr<Epetra_Vector> xbar  = null;
  RefCountPtr<Epetra_Vector> fxbar = null;
  
  //=====coarsest grid==================================================
  if (level==coarsegrid)
  {
    fbar  = rcp(new Epetra_Vector(Copy,*f,0));
    xbar  = rcp(new Epetra_Vector(Copy,*x,0));
    fxbar = rcp(new Epetra_Vector(xbar->Map(),false));
    (*nlnlevel_)[level]->setModifiedSystem(false,NULL,NULL);
    (*nlnlevel_)[level]->computeF(*xbar,*fxbar,NOX::Epetra::Interface::Required::Residual);
    (*nlnlevel_)[level]->setModifiedSystem(true,fbar.get(),fxbar.get());
    int nsmooth = getParameter("nlnML nonlinear smoothing sweeps coarse level",2);
    if (nsmooth)
      *converged = (*nlnlevel_)[level]->Iterate(f,x,nsmooth,&thisnorm);
    x->Update(-1.0,*xbar,1.0);
    (*nlnlevel_)[level]->setModifiedSystem(false,NULL,NULL);
    return true;
  }
  
  //=====setup FAS on this level========================================
  fbar  = rcp(new Epetra_Vector(Copy,*f,0));
  xbar  = rcp(new Epetra_Vector(Copy,*x,0));
  fxbar = rcp(new Epetra_Vector(xbar->Map(),false));
  (*nlnlevel_)[level]->setModifiedSystem(false,NULL,NULL);
  (*nlnlevel_)[level]->computeF(*xbar,*fxbar,NOX::Epetra::Interface::Required::Residual);
  (*nlnlevel_)[level]->setModifiedSystem(true,fbar.get(),fxbar.get());
  
  //=====presmoothing on the FAS system=================================
  double prenorm = thisnorm;
  int    nsmooth;
  if (!level) nsmooth = getParameter("nlnML nonlinear presmoothing sweeps fine level",2);
  else        nsmooth = getParameter("nlnML nonlinear presmoothing sweeps medium level",2);
  if (nsmooth)
    *converged = (*nlnlevel_)[level]->Iterate(f,x,nsmooth,&prenorm);
  if (*converged)
  {
    (*nlnlevel_)[level]->setModifiedSystem(false,NULL,NULL);
    x->Update(-1.0,*xbar,1.0);
    return true;
  }
  
  //=====restrict to next coarser level=================================
  RefCountPtr<Epetra_Vector> xcoarse = rcp((*nlnlevel_)[level]->restrict_to_next_coarser_level(x,level,level+1));
  RefCountPtr<Epetra_Vector> fcoarse = rcp((*nlnlevel_)[level]->restrict_to_next_coarser_level(f,level,level+1));
  
  
  //======call cycle on coarser level===================================
  bool coarseconverged=false;
  FAS_cycle(fcoarse.get(),xcoarse.get(),level+1,&coarseconverged,&prenorm);
  fcoarse = null;
  
  //=====prolongate correction to this level============================
  RefCountPtr<Epetra_Vector> xcorrect = rcp((*nlnlevel_)[level]->prolong_to_this_level(xcoarse.get(),level,level+1));
  xcoarse = null;
  
  //=====apply correction===============================================  
  x->Update(1.0,*xcorrect,1.0);
  xcorrect = null;
  
  //=====postsmoothing on the FAS system================================
  double postnorm = thisnorm;
  if (!level) nsmooth = getParameter("nlnML nonlinear postsmoothing sweeps fine level",2);
  else        nsmooth = getParameter("nlnML nonlinear postsmoothing sweeps medium level",2);
  if (nsmooth)
    *converged = (*nlnlevel_)[level]->Iterate(f,x,nsmooth,&postnorm);
  (*nlnlevel_)[level]->setModifiedSystem(false,NULL,NULL);
  x->Update(-1.0,*xbar,1.0);
  
  return true;
}



/*----------------------------------------------------------------------*
 |  FAS-solver                                                m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_Preconditioner::FAS_cycle_solver() const
{
  return true;
}



/*----------------------------------------------------------------------*
 |  adaptive setup                                            m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_Preconditioner::Adaptivesetup(double** oldns,
                                                Epetra_CrsMatrix* Jac,
                                                int oldnumpde,
                                                int olddimns)
{
#ifdef HAVE_ML_MLAPI

  // init the mlapi
  MLAPI::Init();
  MLAPI::Space FineSpace(Jac->RowMap());

  // fine grid operator
  MLAPI::Operator A(FineSpace,FineSpace,Jac,false);

  // create the cycle options
  Teuchos::ParameterList List;
  List.set("PDE equations", oldnumpde);
  List.set("use default null space", false);
  List.set("smoother: type", "symmetric Gauss-Seidel"); 
  List.set("smoother: sweeps", 1);
  List.set("smoother: damping factor", 1.0);
  List.set("coarse: type", "Amesos-KLU");
  List.set("coarse: max size", getParameter("nlnML coarse: max size",128));
  List.set("max levels", getParameter("nlnML max levels",10));
  List.set("adapt: max reduction", 0.05);
  List.set("adapt: iters fine", 30); // 35
  List.set("adapt: iters coarse", 25); // 20
  List.set("aggregation: damping", 1.33);
  List.set("aggregation: type", "Uncoupled");  // or "METIS", not "VBMETIS"
  
  
  // create the adaptive class
  MLAPI::MultiLevelAdaptiveSA Prec(A,List,oldnumpde,getParameter("nlnML max levels",10));

  // if there's no nullspace vector, create the default one
  if (!(*oldns))
  {
    if (olddimns<oldnumpde) olddimns = oldnumpde;
    int nmyrows = Jac->NumMyRows();
    (*oldns) = new double[olddimns*nmyrows];
    for (int v=0; v<olddimns; ++v)
    {
      for (int i=0; i<nmyrows; ++i) (*oldns)[i] = 0.0;
      int i=v;
      for (; i<nmyrows; )
      {
        (*oldns)[i] = 1.0;
        i +=  olddimns;
      }
    }
  }

  // set the input nullspace
  MLAPI::MultiVector NSfine(FineSpace,olddimns);
  for (int v=0; v<olddimns; ++v)
    for (int i=0; i<NSfine.GetMyLength(); ++i)
      NSfine(i,v) = (*oldns)[v*(NSfine.GetMyLength())+i];
  Prec.SetNullSpace(NSfine);

  // set input numpde
  Prec.SetInputNumPDEEqns(oldnumpde);
  
  // run adatpive setup
  Prec.AdaptCompute(true,getParameter("nlnML additional adaptive nullspace",0));
  
  // get the new nullspace
  MLAPI::MultiVector NSnew = Prec.GetNullSpace();
  int newdimns = NSnew.GetNumVectors();
  
  // delete the old one and copy the new one
  delete [] (*oldns);
  (*oldns) = new double[newdimns*NSnew.GetMyLength()];
  
  for (int v=0; v<newdimns; ++v)
    for (int i=0; i<NSnew.GetMyLength(); ++i)
      (*oldns)[v*(NSnew.GetMyLength())+i] = NSnew(i,v);
  
#if 1
  // scale new candidates to reasonible size
  int dimns = getParameter("nlnML null space: dimension2",1);
  for (int v=dimns; v<newdimns; ++v)
  {
    double length = 0.;
    for (int i=0; i<NSnew.GetMyLength(); ++i)
      length += (*oldns)[v*(NSnew.GetMyLength())+i]*(*oldns)[v*(NSnew.GetMyLength())+i];
    length = sqrt(length);
    length /= 10.;
    for (int i=0; i<NSnew.GetMyLength(); ++i)
      (*oldns)[v*(NSnew.GetMyLength())+i] /= length;
  }
#endif  

  // change the dimension of the nullspace
  dimns = newdimns;
  setParameter("nlnML null space: dimension",dimns);
  
  
#endif
  return true;
}
















#endif
