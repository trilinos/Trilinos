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
 * \file ml_nox_preconditioner1.cpp
 *
 * \class ML_Nox_Preconditioner
 *
 * \brief ML nonlinear preconditioner and solver
 *
 * \date Last update do Doxygen: 31-Mar-05
 *
 */
// ML-headers
#include "ml_common.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) 

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#ifdef ML_MPI // FIXME: do we need this here?
#include <mpi.h>
#endif

// this class
#include "ml_nox_preconditioner.H"

/*----------------------------------------------------------------------*
 |  the class defining the nln-ml-preconditioner             m.gee 11/04|
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     m.gee 11/04|
 *----------------------------------------------------------------------*/
ML_NOX::ML_Nox_Preconditioner::ML_Nox_Preconditioner(
                              ML_NOX::Ml_Nox_Fineinterface&    interface,
                              Epetra_Map&                      dm, 
                              Epetra_Map&                      rm, 
                              Epetra_Comm&                     comm):
interface_(interface),
DomainMap_(dm),
RangeMap_(rm),
comm_(comm)                                             
{
  // label 
  label_         = "ML_Nox_Preconditioner";

  // some flags
  ismatrixfree_ = true;
  matfreelev0_  = true;
  islinearPrec_ = false;
  isnlnCG_      = true;

  // set some default values
  // default values (some of them derived and not supported)
  usetranspose_        = false;
  isinit_              = false;
  fineJac_             = 0;
  destroyfineJac_      = false;
  fineGraph_           = 0;
  ml_graphwrap_        = 0;
  ml_linPrec_          = 0;
  nmatfreelevel_       = 0;
  ml_matfreelevel_     = 0;
  ncalls_NewPrec_      = 0;
  n_nlnlevel_          = 0;
  nlnLevel_            = 0;
  fd_alpha_            = 1.0e-07;
  fd_beta_             = 1.0e-06;
  fd_centered_         = false;
  FAS_normF_           = 1.0e-05;
  FAS_nupdate_         = 1.0e-05;
  FAS_prefinesmooth_   = 2;
  FAS_presmooth_       = 2;
  FAS_postsmooth_      = 2;
  FAS_postfinesmooth_  = 2;
  FAS_maxcycle_        = 250;
  noxsolver_           = 0;
  nitersCG_            = 80;
  offset_newPrec_      = 100;

  // ML stuff
  ml_               = 0;
  ag_               = 0;
  ml_coarsestlev_   = 0;
  ml_nlevel_        = 0;
  ml_nblocks_       = 0;
  ml_blocks_        = 0;
  ml_block_pde_     = 0;
  ml_N_levels_      = 3;
  ml_numPDE_        = 3;
  ml_dim_nullsp_    = 3;
  ml_coarsentype_   = "Uncoupled";
  ml_printlevel_    = 8;
  ml_nnodeperagg_   = 9;
  ml_smoothertype_  = "SGS";
  ml_fsmoothertype_ = "SGS";
  ml_coarsesolve_   = "AmesosKLU";
  ml_maxcoarsesize_ = 50;
  ml_nsmooth_       = new int[ml_N_levels_];
  for (int i=0; i<ml_N_levels_; i++) ml_nsmooth_[i] = 3;

  return;
}

/*----------------------------------------------------------------------*
 |  Set methods for flags/data (public)                      m.gee 03/05|
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::SetCoarsenType(string coarsentype, 
                                                   int maxlevel,
                                                   int maxcoarsesize,
                                                   int nnodeperagg)
{
  if (coarsentype != "Uncoupled" && coarsentype != "MIS" && 
      coarsentype != "METIS" && coarsentype != "VBMETIS")
  {
    cout << "**ERR**: ML_Nox_Preconditioner::SetCoarsenType:\n"
         << "**ERR**: coarsen type: " << coarsentype << " not recognized!\n"
         << "**ERR**: Using previously set coarsen method: " << ml_coarsentype_ << endl
         << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; 
    coarsentype = ml_coarsentype_;
  }
  
  if (maxlevel<2)
  {
    cout << "**ERR**: ML_Nox_Preconditioner::SetCoarsenType:\n"
         << "**ERR**: maxlevel < 2 out of range, minimum is 2 level!\n"
         << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; 
    maxlevel=2;
  }
  
  ml_N_levels_      = maxlevel;
  ml_coarsentype_   = coarsentype;
  ml_nnodeperagg_   = nnodeperagg;
  ml_maxcoarsesize_ = maxcoarsesize;

  return true;
}                              
/*----------------------------------------------------------------------*
 |  Set methods for flags/data (public)                      m.gee 03/05|
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::SetSmoothers(string finesmoothertype, 
                                                 string smoothertype, 
                                                 string coarsesolve)
{
  if (smoothertype != "SGS" && 
      smoothertype != "Jacobi" &&
      smoothertype != "AmesosKLU")
  {
    cout << "**ERR**: ML_Nox_Preconditioner::SetSmoothers:\n"
         << "**ERR**: smoother type: " << smoothertype << " not recognized!\n"
         << "**ERR**: using previous smoother: " << ml_smoothertype_ << endl
         << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n";   
  }
  else
     ml_smoothertype_ = smoothertype;

  if (finesmoothertype != "SGS" && 
      finesmoothertype != "Jacobi" &&
      finesmoothertype != "AmesosKLU")
  {
    cout << "**ERR**: ML_Nox_Preconditioner::SetSmoothers:\n"
         << "**ERR**: fine level smoother type: " << finesmoothertype << " not recognized!\n"
         << "**ERR**: using previous smoother: " << ml_fsmoothertype_ << endl
         << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; 
  }
  else
     ml_fsmoothertype_ = finesmoothertype;
  
  if (coarsesolve != "SGS" && 
      coarsesolve != "Jacobi" && 
      coarsesolve != "AmesosKLU")
  {
    cout << "**ERR**: ML_Nox_Preconditioner::ML_Nox_Preconditioner:\n"
         << "**ERR**: coarse solver: " << coarsesolve << " not recognized!\n"
         << "**ERR**: using previous coarse solver: " << ml_coarsesolve_ << endl
         << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n";
  }
  else
     ml_coarsesolve_ = coarsesolve;
  return true;
}                              
/*----------------------------------------------------------------------*
 |  Set methods for flags/data (public)                      m.gee 03/05|
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::SetPrintLevel(int printlevel)
{
  if (printlevel<0 || printlevel>10)
  {
    cout << "**ERR**: ML_Nox_Preconditioner::SetPrintLevel:\n"
         << "**ERR**: printlevel out of range, using previous one\n"
         << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; 
    return false;     
  }
  ml_printlevel_ = printlevel;
  return true;
}                              

/*----------------------------------------------------------------------*
 |  Set methods for flags/data (public)                      m.gee 03/05|
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::SetRecomputeOffset(int offset)
{ 
  if (offset<0)
  {
    cout << "**ERR**: ML_Nox_Preconditioner::SetRecomputeOffset:\n"
         << "**ERR**: offset out of range, using previous one\n"
         << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; 
    return false;     
  }
  offset_newPrec_ = offset;
  return true;
}                              

/*----------------------------------------------------------------------*
 |  Set methods for flags/data (public)                      m.gee 03/05|
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::SetDimensions(int numPDE, int dimNS)
{ 
  if (numPDE<1)
  {
    cout << "**ERR**: ML_Nox_Preconditioner::SetDimensions:\n"
         << "**ERR**: numPDE out of range, using previous one\n"
         << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; 
  }
  else
     ml_numPDE_ = numPDE;
     
  if (dimNS<0)
  {
    cout << "**ERR**: ML_Nox_Preconditioner::SetDimensions:\n"
         << "**ERR**: dimNS<0 out of range using dimNS = numPDE\n"
         << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n";
    dimNS = numPDE;
  }
  if (dimNS<numPDE)
  {
    cout << "**ERR**: ML_Nox_Preconditioner::SetDimensions:\n"
         << "**ERR**: dimNS<numPDE bad choice\n"
         << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n";
    
  }
  ml_dim_nullsp_ = dimNS;
  return true;
}                              

/*----------------------------------------------------------------------*
 |  Set methods for flags/data (public)                      m.gee 03/05|
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::SetSmootherSweeps(int* nsmooth)
{ 
  for (int i=0; i<ml_N_levels_; i++)
  if (nsmooth[i]<0)
  {
    cout << "**ERR**: ML_Nox_Preconditioner::SetSmootherSweeps:\n"
         << "**ERR**: nsmooth[i]<0 out of range, using nsmooth[i]=0\n"
         << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; 
    nsmooth[i] = 0;
  }
  if (ml_nsmooth_) delete [] ml_nsmooth_;
  ml_nsmooth_ = new int[ml_N_levels_];
  for (int i=0; i<ml_N_levels_; i++)
     ml_nsmooth_[i] = nsmooth[i];
     
  return true;
}                              

/*----------------------------------------------------------------------*
 |  Set methods for flags/data (public)                      m.gee 03/05|
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::SetConvergenceCriteria(double FAS_normF, 
                                                           double FAS_nupdate)
{ 
  FAS_normF_   = FAS_normF;
  FAS_nupdate_ = FAS_nupdate;
  return true;
}                              

/*----------------------------------------------------------------------*
 |  Set methods for flags/data (public)                      m.gee 03/05|
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::SetFiniteDifferencing(bool centered,
                                                          double alpha,
                                                          double beta)
{ 
  fd_alpha_    = alpha;
  fd_beta_     = beta;
  fd_centered_ = centered;
  return true;
}                              

/*----------------------------------------------------------------------*
 |  Set methods for flags/data (public)                      m.gee 03/05|
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::SetFAScycle(int prefsmooth,int presmooth,
                                                int postsmooth,int postfsmooth,
                                                int maxcycle)
{ 
  FAS_prefinesmooth_    = prefsmooth;
  FAS_presmooth_        = presmooth;
  FAS_postsmooth_       = postsmooth;
  FAS_postfinesmooth_   = postfsmooth;
  FAS_maxcycle_         = maxcycle;
  return true;
}                              

/*----------------------------------------------------------------------*
 |  Set methods for flags/data (public)                      m.gee 03/05|
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::SetNonlinearMethod(bool islinPrec, 
                                                       bool isnlnCG, 
                                                       bool ismatrixfree, 
                                                       bool ismatfreelev0)
{ 
  islinearPrec_ = islinPrec;
  ismatrixfree_ = ismatrixfree;
  matfreelev0_  = ismatfreelev0;
  isnlnCG_      = isnlnCG;
  if (islinPrec && ismatrixfree && !ismatfreelev0 &&
      (ml_fsmoothertype_== "Jacobi" || 
       ml_smoothertype_ == "Jacobi" ||
       ml_coarsesolve_  == "Jacobi"))
  {
    cout << "**ERR**: ML_Nox_Preconditioner::SetNonlinearMethod:\n"
         << "**ERR**: linearPrec && matrixfree & Jacobi smoother doesn't work!\n"
         << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; throw -1;
  }
  if (!islinPrec && !isnlnCG && ismatrixfree && !ismatfreelev0 &&
      (ml_fsmoothertype_== "Jacobi" || 
       ml_smoothertype_ == "Jacobi" ||
       ml_coarsesolve_  == "Jacobi"))
  {
    cout << "**ERR**: ML_Nox_Preconditioner::SetNonlinearMethod:\n"
         << "**ERR**: NonlinPrec && matrixfree & Newton smoother && Jacobi doesn't work!\n"
         << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; throw -1;
  }
     
  return true;
}                              

/*----------------------------------------------------------------------*
 |  Destructor (public)                                      m.gee 11/04|
 *----------------------------------------------------------------------*/
ML_NOX::ML_Nox_Preconditioner::~ML_Nox_Preconditioner()
{
  if (ml_linPrec_)
     delete ml_linPrec_;
  ml_linPrec_ = 0;
  
  if (ml_nsmooth_)
    delete [] ml_nsmooth_;
  ml_nsmooth_ = 0;

  if (ml_graphwrap_)
    delete ml_graphwrap_;
  ml_graphwrap_ = 0;

  if (ml_matfreelevel_)
    destroy_matfreelevels(&ml_matfreelevel_,nmatfreelevel_);
  
  if (ml_blocks_)
     delete [] ml_blocks_;
  ml_blocks_ = 0;
  
  if (ag_)
     ML_Aggregate_Destroy(&ag_);
  ag_ = 0;
  
  if (ml_)
     ML_Destroy(&ml_);
  ml_ = 0;
  
  if (nlnLevel_)
     destroy_nonlinearlevels(&nlnLevel_,n_nlnlevel_);
     
  if (ml_blocks_)
     delete [] ml_blocks_;
  ml_blocks_ = 0;
  
  if (ml_block_pde_)
     delete [] ml_block_pde_;
  ml_block_pde_ = 0;
  
  if (destroyfineJac_ && fineJac_)
     delete fineJac_;
  fineJac_ = 0;

  return;
}

/*----------------------------------------------------------------------*
 |  destroy vector of nonlinearlevels                        m.gee 1/05 |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::destroy_nonlinearlevels(
                                  ML_NOX::ML_Nox_NonlinearLevel*** nln, 
                                  int nlevel)
{
  int i;
  if (*nln)
  {
     (*nln)[0]->destroyP();
     for (i=0; i<nlevel; i++)
     {
        if ((*nln)[i])
        {
           (*nln)[i]->setP(NULL);
           delete (*nln)[i];
           (*nln)[i] = 0;
        }
     }
     delete [] (*nln);
     *nln        = 0;
     n_nlnlevel_ = 0;
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  destroy vector of matrixfreelevels                       m.gee 1/05 |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::destroy_matfreelevels(
                                  ML_NOX::ML_Nox_MatrixfreeLevel*** m, 
                                  int nlevel)
{
  int i;
  if (*m)
  {
     (*m)[0]->destroyP();
     for (i=0; i<nlevel; i++)
     {
        if ((*m)[i])
        {
           (*m)[i]->setP(NULL);
           delete (*m)[i];
           (*m)[i] = 0;
        }
     }
     delete [] (*m);
     *m             = 0;
     nmatfreelevel_ = 0;
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  compute the multilevel preconditioner (public)           m.gee 11/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::computePreconditioner(
                                              const Epetra_Vector& x,
                                              NOX::Parameter::List* precParams)
{
   bool flag;
   int offset = getoffset();
   if (ncalls_NewPrec_ % offset == 0)
         setinit(false);
   if (isinit() == false)
   {
      if (comm_.MyPID()==0 && ml_printlevel_ > 0 )
         cout << "ML: ML_Nox_Preconditioner::computePreconditioner: (re)computing ML-Preconditioner\n";
      
      // save number of calls to computeF up to now
      int ncalls = interface_.getnumcallscomputeF();
      
      // reset number of calls to computeF to zero to measure number of calls in compPrec()
      interface_.setnumcallscomputeF(0);
      
      double t0 = GetClock();
      flag = compPrec(x);
      double t1 = GetClock();
      
      if (comm_.MyPID()==0 && ml_printlevel_ > 0 )
      {
         cout << "ML: Setup time for preconditioner: " << (t1-t0) << " sec\n";
         cout << "ML: Number of calls to fineinterface.computeF() in setup: " 
              << interface_.getnumcallscomputeF() << endl;
      }
      
      // reset the number of calls to computeF to what it was before
      interface_.setnumcallscomputeF(ncalls);
      
      if (flag==true)
         setinit(true);
      else
      {
         cout << "**ERR**: ML_Nox_Preconditioner::computePreconditioner:\n"
              << "**ERR**: compPrec returned false\n"
              << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; throw -1;
      }
      setinit(true);
   }
   ++ncalls_NewPrec_;
   return flag;
}                                                  
/*----------------------------------------------------------------------*
 |  compute the multilevel preconditioner (private)          m.gee 11/04|
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::compPrec(const Epetra_Vector& x)
{
  int i;

  // build hierarchy with given Jacobian
  if (ismatrixfree_ == false)
  {
    // check for valid Jacobian, if not, recompute
    if (interface_.isnewJacobian()==false)
       interface_.computeJacobian(x);

    // get the fine grid Jacobian
    fineJac_ = interface_.getJacobian();
    if (fineJac_ == NULL)
    {
      cout << "**ERR**: ML_Nox_Preconditioner::compPrec:\n"
           << "**ERR**: interface_.getJacobian() returned NULL\n"
           << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; throw -1;
    }
  }

  // build matrixfree hierachy and probe for all operators
  else if (ismatrixfree_ == true && matfreelev0_ == false ) 
  {
    // get the graph from the user-interface
    fineGraph_ = interface_.getGraph();
    if (fineGraph_ == NULL)
    {
      cout << "**ERR**: ML_Nox_Preconditioner::compPrec:\n"
           << "**ERR**: interface_.getGraph() returned NULL\n"
           << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; throw -1;
    } 

     // Maybe there exists a ml_graphwrap_ from a previous call?
     if (ml_graphwrap_ != NULL)
       delete ml_graphwrap_;

     // this guy wraps an Epetra_CrsGraph as an Epetra_RowMatrix
     // (there exists a ML_Operator - wrapper for an Epetra_RowMatrix
     //  we can then use)
     ml_graphwrap_ = new ML_Epetra::CrsGraphWrapper(*fineGraph_,DomainMap_,RangeMap_,comm_);
  }

  // probe for operator on level 0 and switch to ismatrixfree_==false
  else if (ismatrixfree_ == true && matfreelev0_ == true)
  {
    fineJac_ = ML_Nox_computeFineLevelJacobian(x);
    // this class is in charge of destroying the operator fineJac_
    destroyfineJac_ = true;
    // turn to ismatrixfree_==false
    ismatrixfree_ = false;
  }
  
  // check whether ML and ML_Aggregate handles exist and (re)create them
  if (ag_ != NULL)
  {
    ML_Aggregate_Destroy(&ag_);
    ML_Aggregate_Create(&ag_);
  }
  else
     ML_Aggregate_Create(&ag_);
     
  if (ml_ != NULL)
  {
    ML_Destroy(&ml_);
    ML_Create(&ml_,ml_N_levels_);
  }
  else
    ML_Create(&ml_,ml_N_levels_);
  
  // set matrix for fine grid
  // NOTE: this is either a wrapped graph or the real fine grid Jacobian
  if (ismatrixfree_ == false)
    EpetraMatrix2MLMatrix(ml_,0,(dynamic_cast<Epetra_RowMatrix*>(fineJac_)));
  else
    EpetraMatrix2MLMatrix(ml_,0,(dynamic_cast<Epetra_RowMatrix*>(ml_graphwrap_)));
  
  // set coarsening scheme
  if (ml_coarsentype_ == "Uncoupled")
    ML_Aggregate_Set_CoarsenScheme_Uncoupled(ag_);
  if (ml_coarsentype_ == "MIS")
    ML_Aggregate_Set_CoarsenScheme_MIS(ag_);
  if (ml_coarsentype_ == "METIS")
  {
    ML_Aggregate_Set_CoarsenScheme_METIS(ag_);
    for (i=0; i<ml_N_levels_; i++)
       ML_Aggregate_Set_NodesPerAggr(ml_,ag_,i,ml_nnodeperagg_);
  } 
  if (ml_coarsentype_ == "VBMETIS")
  {
    if (ml_blocks_)    delete [] ml_blocks_;    ml_blocks_ = 0;
    if (ml_block_pde_) delete [] ml_block_pde_; ml_block_pde_ = 0;

    // get the block information from the application
    bool ok = interface_.getBlockInfo(&ml_nblocks_,&ml_blocks_,&ml_block_pde_);
    if (!ok) // we failed to get the blocks, do plain METIS instead
    {
       cout << "**WRN**: ML_Nox_Preconditioner::compPrec:\n"
            << "**WRN**: interface returned no blocks,\n"
            << "**WRN**: using aggregation scheme METIS instead of VBMETIS\n";
       ml_coarsentype_ = "METIS";
       ML_Aggregate_Set_CoarsenScheme_METIS(ag_);
    }
    else
    {
       ML_Aggregate_Set_CoarsenScheme_VBMETIS(ag_);
       int N_update = interface_.getMap().NumMyElements();
       ML_Aggregate_Set_Vblocks_CoarsenScheme_VBMETIS(ag_,0,ml_N_levels_,
                                                      ml_nblocks_,ml_blocks_,
                                                      ml_block_pde_,N_update);
    if (ml_blocks_)    delete [] ml_blocks_;    ml_blocks_ = 0;
    if (ml_block_pde_) delete [] ml_block_pde_; ml_block_pde_ = 0;
    }
    for (i=0; i<ml_N_levels_; i++)
       ML_Aggregate_Set_NodesPerAggr(ml_,ag_,i,ml_nnodeperagg_);
  }
  
  // set default damping factor
  if (ismatrixfree_ == true)
    ML_Aggregate_Set_DampingFactor(ag_, 0.0); // cannot do smoothed aggregation on a graph
  else
    ML_Aggregate_Set_DampingFactor(ag_, 1.5); // do smoothed aggregation on a Jacobian

  // set default threshold to zero
  ML_Aggregate_Set_Threshold(ag_, 0.0);
  
  // set max coarsest grid size
  ML_Aggregate_Set_MaxCoarseSize(ag_,ml_maxcoarsesize_);
  
  // Calculate spectral norm
  ML_Aggregate_Set_SpectralNormScheme_Calc(ag_);
  
  // set ML printlevel
  ML_Set_PrintLevel(ml_printlevel_);
  
  // set the nullspace (currently default nullspace)
  if (ismatrixfree_ == true)
     i = ml_graphwrap_->NumMyRows();
  else
     i = fineJac_->NumMyRows();
  // get the nullspace
  double* nullsp = interface_.Get_Nullspace(i,ml_numPDE_,ml_dim_nullsp_);
  if (nullsp)
  {
     ML_Aggregate_Set_NullSpace(ag_,ml_numPDE_,ml_dim_nullsp_,nullsp,i);
     // delete [] nullsp; nullsp = 0; // FIXME, who deletes this?
  }
  else
  {
     cout << "**WRN**: ML_Nox_Preconditioner::compPrec:\n"
          << "**WRN**: interface returned no nullspace,\n"
          << "**WRN**: using ML's default nullspace\n";
     if (ml_numPDE_ != ml_dim_nullsp_)
     {
        cout << "**WRN**: ML_Nox_Preconditioner::compPrec:\n"
             << "**WRN**: with default nullspace, nullspace-dimension must match number of PDEs per node\n"
             << "**WRN**: numPDE = " << ml_numPDE_ << ", dim_nullsp = " << ml_dim_nullsp_ << "\n"
             << "**WRN**: continue with setting dim_nullsp = ml_numPDE_\n";
        ml_dim_nullsp_ = ml_numPDE_;
     }
     ML_Aggregate_Set_NullSpace(ag_,ml_numPDE_,ml_dim_nullsp_,NULL,i);
  }
  // keep the aggregation information
  ag_->keep_agg_information = 0;
  
  // build hierarchy
  ml_nlevel_ = ML_Gen_MGHierarchy_UsingAggregation(ml_,0,ML_INCREASING,ag_);

  if (ml_nlevel_<2)
  {
     cout << "**ERR**: ML_Nox_Preconditioner::compPrec:\n"
          << "**ERR**: number of levels generated is " << ml_nlevel_ << "\n"
          << "**ERR**: this algorithm relies on at least nlevel >=2 !\n"
          << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; throw -1;
  }

  ml_coarsestlev_ = ml_nlevel_-1;
  
  if (ismatrixfree_ == false && islinearPrec_ == true)
  {
     isinit_ = ML_Nox_compute_Jacobian_Linearpreconditioner(x);
     if (isinit_ == true)
       return(true);
     else
       return(false);
  }
  else if (ismatrixfree_ == true && islinearPrec_ == true)
  {
     isinit_ = ML_Nox_compute_Matrixfree_Linearpreconditioner(x);
     if (isinit_ == true)
       return(true);
     else
       return(false);
  }
  else if (ismatrixfree_ == false && islinearPrec_ == false)
  {
     isinit_ = ML_Nox_compute_Jacobian_Nonlinearpreconditioner(x);
     if (isinit_ == true)
       return(true);
     else
       return(false);
  }
  else if (ismatrixfree_ == true && islinearPrec_ == false)
  {
     isinit_ = ML_Nox_compute_Matrixfree_Nonlinearpreconditioner(x);
     if (isinit_ == true)
       return(true);
     else
       return(false);
  }
  
  cout << "**ERR**: ML_Nox_Preconditioner::compPrec:\n";
  cout << "**ERR**: something wrong with the flags ismatrixfree_ and islinearPrec_\n";
  cout << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; throw -1;
  isinit_ = false;
  return(false);
}

/*----------------------------------------------------------------------*
 |  compute the linear multilevel preconditioner (private)   m.gee 11/04|
 |  with a given Jacobian                                               |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::ML_Nox_compute_Jacobian_Linearpreconditioner(
                                                   const Epetra_Vector& x)
{
   int i;

   // choose some smoothers on level ==0
   if (ml_fsmoothertype_ == "SGS")
      ML_Gen_Smoother_SymGaussSeidel(ml_,0,ML_BOTH,ml_nsmooth_[0],1.);
   if (ml_fsmoothertype_ == "Jacobi")
   {
      ML_Gen_Smoother_Jacobi(ml_,0,ML_PRESMOOTHER, ml_nsmooth_[0],.6);
      ML_Gen_Smoother_Jacobi(ml_,0,ML_POSTSMOOTHER,ml_nsmooth_[0],.6);
   }
   if (ml_fsmoothertype_ == "BSGS")
   {
     ML_Gen_Blocks_Aggregates(ag_,0,&ml_nblocks_,&ml_blocks_);
     ML_Gen_Smoother_VBlockSymGaussSeidel(ml_,0,ML_BOTH,ml_nsmooth_[0],1.,ml_nblocks_,ml_blocks_);
   }
   if (ml_fsmoothertype_ == "AmesosKLU")
     ML_Gen_Smoother_Amesos(ml_,0,ML_AMESOS_KLU,-1,0.0);
   
   // choose some smoothers on level > 0
   for (i=1; i<ml_coarsestlev_; i++)
   {
      if (ml_smoothertype_ == "SGS")
      {
         ML_Gen_Smoother_SymGaussSeidel(ml_,i,ML_BOTH,ml_nsmooth_[i],1.);
         continue;
      }
      if (ml_smoothertype_ == "Jacobi")
      {
         ML_Gen_Smoother_Jacobi(ml_,i,ML_PRESMOOTHER, ml_nsmooth_[i],.4);
         ML_Gen_Smoother_Jacobi(ml_,i,ML_POSTSMOOTHER,ml_nsmooth_[i],.4);
         continue;
      }
      if (ml_smoothertype_ == "BSGS")
      {
        ML_Gen_Blocks_Aggregates(ag_,i,&ml_nblocks_,&ml_blocks_);
        ML_Gen_Smoother_VBlockSymGaussSeidel(ml_,i,ML_BOTH,ml_nsmooth_[i],1.,ml_nblocks_,ml_blocks_);
         continue;
      }
      if (ml_smoothertype_ == "AmesosKLU")
        ML_Gen_Smoother_Amesos(ml_,i,ML_AMESOS_KLU,-1,0.0);
   }
   // choose a coarse grid solver
   if (ml_coarsesolve_ == "AmesosKLU")
      ML_Gen_Smoother_Amesos(ml_,ml_coarsestlev_,ML_AMESOS_KLU,-1,0.0);
   if (ml_coarsesolve_ == "SGS")
      ML_Gen_Smoother_SymGaussSeidel(ml_,ml_coarsestlev_,ML_BOTH,ml_nsmooth_[ml_coarsestlev_],1.);
   if (ml_coarsesolve_ == "Jacobi")
   {
      ML_Gen_Smoother_Jacobi(ml_,ml_coarsestlev_,ML_PRESMOOTHER, ml_nsmooth_[ml_coarsestlev_],.4);
      ML_Gen_Smoother_Jacobi(ml_,ml_coarsestlev_,ML_POSTSMOOTHER,ml_nsmooth_[ml_coarsestlev_],.4);
   }   
   if (ml_coarsesolve_ == "BSGS")
   {
      ML_Gen_Blocks_Aggregates(ag_,ml_coarsestlev_,&ml_nblocks_,&ml_blocks_);
      ML_Gen_Smoother_VBlockSymGaussSeidel(ml_,ml_coarsestlev_,ML_BOTH,ml_nsmooth_[ml_coarsestlev_],1.,ml_nblocks_,ml_blocks_);
   }

   // build the ML_MultilevelOperator
   ml_linPrec_ = new ML_Epetra::MultiLevelOperator(ml_,comm_,DomainMap_,RangeMap_);


#if 0  // test the MultiLevelOperator and exit
   {
   Epetra_Vector xthis(x);
   cout << "Test of preconditioner " << endl;
   Epetra_Vector *out = new Epetra_Vector(xthis);
   out->PutScalar(0.0);
   cout << "Input\n";
   xthis.PutScalar(1.0);
   interface_.getJacobian()->Multiply(false,xthis,*out);
   double norm;
   out->Norm1(&norm);
   cout << "Norm = " << norm << endl;
   xthis.PutScalar(3.0);
   cout << "out between\n";
   cout << *out;
   ml_linPrec_->ApplyInverse(*out,xthis);
   cout << "Output\n";
   cout << xthis;
   delete out; out = 0;
   exit(0);
   }
#endif   


#if 0 // print out all level's Jacobians
   for (i=0; i<=ml_coarsestlev_; i++)
   {
      cout << "Matrix on level " << i << "\n";
      Epetra_CrsMatrix* tmpMat  = 0;
      int               maxnnz  = 0;
      double            cputime = 0.0;
      ML_Operator2EpetraCrsMatrix(&(ml_->Amat[i]), tmpMat, maxnnz, false, cputime);
      cout << *tmpMat;
      delete tmpMat;
      tmpMat = 0;
   }
   exit(0);
#endif
   
return (true);
}

/*----------------------------------------------------------------------*
 |  compute the linear multilevel preconditioner (private)   m.gee 11/04|
 |  without Jacobian                                                    |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::ML_Nox_compute_Matrixfree_Linearpreconditioner(
                                                 const Epetra_Vector& xfine)
{
   int i;
   int j;
   
   // extract the Prolongators from the ml-hierarchy and create Epetra_CrsMatrices from them
   // NOTE: as we are doing the linear preconditioner here, we keep the P's in the hierarchy
   //       for further use in the MG-cycle
   Epetra_CrsMatrix** P  = new (Epetra_CrsMatrix*)[ml_nlevel_];
   for (i=0; i<ml_nlevel_; i++)
      P[i] = 0;
   for (i=1; i<ml_nlevel_; i++) // there is no Pmat on level 0
   {
      int    maxnnz  = 0;
      double cputime = 0.0;
      ML_Operator2EpetraCrsMatrix(&(ml_->Pmat[i]), P[i], maxnnz, false, cputime);
      if (ml_printlevel_>0 && 0 == comm_.MyPID())
            cout << "matrixfreeML (level " << i << "): extraction of P in " << cputime << " sec\n";
   }
   
   // construct the vector of coarse level problems
   // the coarse level problems also create the coarse level interface to the user code
   if (ml_matfreelevel_) // the matfreelevels already exist, clean them for reuse
   {
      ml_matfreelevel_[0]->destroyP();
      for (i=0; i<nmatfreelevel_; i++)
         ml_matfreelevel_[i]->setP(NULL);
   }
   else // create new vector of matfreelevels
   {
      nmatfreelevel_   = ml_nlevel_;
      ml_matfreelevel_ = new (ML_NOX::ML_Nox_MatrixfreeLevel*)[nmatfreelevel_];
      for (i=0; i<nmatfreelevel_; i++)
         ml_matfreelevel_[i] = 0;
   }
   
   // loop levels and create the Operators and coarse interfaces
   for (i=0; i<ml_nlevel_; i++)
   {
      if (comm_.MyPID()==0 && ml_printlevel_ > 5 )
         cout << "\nmatrixfreeML (level " << i << "): Entering FD-coarselevel (re)construction\n";
      fflush(stdout);

      bool isJacobismoother=false;
      if (i==0 && ml_fsmoothertype_ == "Jacobi")
         isJacobismoother=true;
      else if (i<ml_coarsestlev_ && ml_smoothertype_ == "Jacobi")
         isJacobismoother=true;
      else
         isJacobismoother=false;
      if (i==ml_coarsestlev_ && ml_coarsesolve_ == "Jacobi")
         isJacobismoother=true;
      
      if (!(ml_matfreelevel_[i])) // create a new level   
         ml_matfreelevel_[i] = new ML_NOX::ML_Nox_MatrixfreeLevel(i,ml_nlevel_,ml_printlevel_,ml_,
                                                                  ag_,P,interface_,comm_,xfine,
                                                                  fd_alpha_,fd_beta_,fd_centered_,
                                                                  isJacobismoother);
      else // redo an existing level
         ml_matfreelevel_[i]->recreateLevel(i,ml_nlevel_,ml_printlevel_,ml_,
                                            ag_,P,interface_,comm_,xfine);
      // get this level's operator 
      Epetra_CrsMatrix* tmpMat = ml_matfreelevel_[i]->getunderlyingMatrix();

      if (!tmpMat)
      {
         cout << "**ERR**: ML_Nox_Preconditioner::ML_Nox_compute_Matrixfree_Linearpreconditioner:\n"
              << "**ERR**: ML_Epetra::ML_Nox_MatrixfreeLevel::getunderlyingMatrix() on level " << i << "\n"
              << "**ERR**: returned NULL-ptr\n"
              << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      // create a ML_Operator
      // NOTE:
      // As I understand this will wrap the Epetra_CrsMatrix from ML_Nox_MatrixfreeLevel,
      // so we must not destroy the ML_Nox_MatrixfreeLevel class unless the ML-hierarchy 
      // is no longer needed
      ML_Operator* tmpOperator = ML_Operator_Create(ml_->comm);
      Epetra2MLMatrix(tmpMat,tmpOperator);
      // move it to the hierarchy, this destroys tmpOperator
      ML_Operator_Move2HierarchyAndDestroy(&tmpOperator,&(ml_->Amat[i]));
   }

   // the hierarchy is done, so set smoothers 

   // choose some smoothers on level ==0
   if (ml_fsmoothertype_ == "SGS")
      ML_Gen_Smoother_SymGaussSeidel(ml_,0,ML_BOTH,ml_nsmooth_[0],1.);
   if (ml_fsmoothertype_ == "Jacobi")
   {
      ML_Gen_Smoother_Jacobi(ml_,0,ML_PRESMOOTHER, ml_nsmooth_[0],.4);
      ML_Gen_Smoother_Jacobi(ml_,0,ML_POSTSMOOTHER,ml_nsmooth_[0],.4);
   }
   if (ml_fsmoothertype_ == "BSGS")
   {
     ML_Gen_Blocks_Aggregates(ag_,0,&ml_nblocks_,&ml_blocks_);
     ML_Gen_Smoother_VBlockSymGaussSeidel(ml_,0,ML_BOTH,ml_nsmooth_[0],1.,ml_nblocks_,ml_blocks_);
   }
   if (ml_fsmoothertype_ == "AmesosKLU")
     ML_Gen_Smoother_Amesos(ml_,0,ML_AMESOS_KLU,-1,0.0);
   
   // choose some smoothers on level > 0
   for (i=1; i<ml_coarsestlev_; i++)
   {
      if (ml_smoothertype_ == "SGS")
      {
         ML_Gen_Smoother_SymGaussSeidel(ml_,i,ML_BOTH,ml_nsmooth_[i],1.);
         continue;
      }
      if (ml_smoothertype_ == "Jacobi")
      {
         ML_Gen_Smoother_Jacobi(ml_,i,ML_PRESMOOTHER, ml_nsmooth_[i],.4);
         ML_Gen_Smoother_Jacobi(ml_,i,ML_POSTSMOOTHER,ml_nsmooth_[i],.4);
         continue;
      }
      if (ml_smoothertype_ == "BSGS")
      {
        ML_Gen_Blocks_Aggregates(ag_,i,&ml_nblocks_,&ml_blocks_);
        ML_Gen_Smoother_VBlockSymGaussSeidel(ml_,i,ML_BOTH,ml_nsmooth_[i],1.,ml_nblocks_,ml_blocks_);
         continue;
      }
      if (ml_smoothertype_ == "AmesosKLU")
        ML_Gen_Smoother_Amesos(ml_,i,ML_AMESOS_KLU,-1,0.0);
   }
   
   // choose a coarse solver
   if (ml_coarsesolve_ == "AmesosKLU")
      ML_Gen_Smoother_Amesos(ml_,ml_coarsestlev_,ML_AMESOS_KLU,-1,0.0);
   if (ml_coarsesolve_ == "SGS")
      ML_Gen_Smoother_SymGaussSeidel(ml_,ml_coarsestlev_,ML_BOTH,ml_nsmooth_[ml_coarsestlev_],1.);
   if (ml_coarsesolve_ == "Jacobi")
   {
      ML_Gen_Smoother_Jacobi(ml_,ml_coarsestlev_,ML_PRESMOOTHER, ml_nsmooth_[ml_coarsestlev_],.4);
      ML_Gen_Smoother_Jacobi(ml_,ml_coarsestlev_,ML_POSTSMOOTHER,ml_nsmooth_[ml_coarsestlev_],.4);
   }   
   if (ml_coarsesolve_ == "BSGS")
   {
      ML_Gen_Blocks_Aggregates(ag_,ml_coarsestlev_,&ml_nblocks_,&ml_blocks_);
      ML_Gen_Smoother_VBlockSymGaussSeidel(ml_,ml_coarsestlev_,ML_BOTH,ml_nsmooth_[ml_coarsestlev_],1.,ml_nblocks_,ml_blocks_);
   }
   
   // build the ML_MultilevelOperator
   if (ml_linPrec_) 
      delete ml_linPrec_;
   ml_linPrec_ = new ML_Epetra::MultiLevelOperator(ml_,comm_,DomainMap_,RangeMap_);
   
   return(true);
}

/*----------------------------------------------------------------------*
 |  apply inverse of operator (public)                       m.gee 11/04|
 *----------------------------------------------------------------------*/
Epetra_CrsMatrix* ML_NOX::ML_Nox_Preconditioner::ML_Nox_computeFineLevelJacobian(
                                                  const Epetra_Vector& x)
{
  // make a copy of the graph
  Epetra_CrsGraph* graph = deepcopy_graph(interface_.getGraph());

  // create coloring of the graph
  if (ml_printlevel_>0 && comm_.MyPID()==0)
     cout << "matrixfreeML (level 0): Entering Coloring on level 0\n";
  double t0 = GetClock();
  EpetraExt::CrsGraph_MapColoring::ColoringAlgorithm algType = 
                                  EpetraExt::CrsGraph_MapColoring::GREEDY;

  EpetraExt::CrsGraph_MapColoring* MapColoring = 
                   new EpetraExt::CrsGraph_MapColoring(algType,0,false,0);

  Epetra_MapColoring* colorMap = &(*MapColoring)(*graph);

  EpetraExt::CrsGraph_MapColoringIndex* colorMapIndex = 
                      new EpetraExt::CrsGraph_MapColoringIndex(*colorMap);
                      
  vector<Epetra_IntVector>* colorcolumns = &(*colorMapIndex)(*graph);
  double t1 = GetClock();
  if (ml_printlevel_>0 && comm_.MyPID()==0)
     cout << "matrixfreeML (level 0): Proc " << comm_.MyPID() <<" Coloring time is " << (t1-t0) << " sec\n";
  
  // construct the FiniteDifferenceColoring-Matrix
  if (ml_printlevel_>0 && comm_.MyPID()==0)
  {
     cout << "matrixfreeML (level 0): Entering Construction FD-Operator on level 0\n";
     fflush(stdout);
  }
  
  t0 = GetClock();
  int ncalls = interface_.getnumcallscomputeF();
  interface_.setnumcallscomputeF(0);
  
  NOX::EpetraNew::FiniteDifferenceColoring* FD = 
           new NOX::EpetraNew::FiniteDifferenceColoring(interface_,x,*graph,
                                                        *colorMap,*colorcolumns,
                                                        true,false,
                                                        fd_beta_,fd_alpha_);
  if (fd_centered_)
    FD->setDifferenceMethod(NOX::EpetraNew::FiniteDifferenceColoring::Centered);

  bool err = FD->computeJacobian(x); 
  if (err==false)
  {
    cout << "**ERR**: ML_NOX::ML_Nox_Preconditioner::ML_Nox_computeFineLevelJacobian:\n"
         << "**ERR**: NOX::Epetra::FiniteDifferenceColoring returned an error on level 0" 
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }

  t1 = GetClock();
  if (ml_printlevel_>0 && comm_.MyPID()==0)
  {
     cout << "matrixfreeML (level 0): colored Finite Differencing time :" << (t1-t0) << " sec\n";
     cout << "matrixfreeML (level 0): colored Finite Differencing number of calls to computeF : " 
          << interface_.getnumcallscomputeF() << endl;
     fflush(stdout);
  }
  
  interface_.setnumcallscomputeF(ncalls);

  Epetra_CrsMatrix* B = dynamic_cast<Epetra_CrsMatrix*>(&(FD->getUnderlyingMatrix()));                       
  Epetra_CrsMatrix* A = new Epetra_CrsMatrix(*B);
  
  // tidy up
  delete FD;            FD = 0;
  delete MapColoring;   MapColoring = 0;
  delete colorMap;      colorMap = 0;
  delete colorMapIndex; colorMapIndex = 0;
  delete colorcolumns;  colorcolumns = 0;
  
  return A;
}

/*----------------------------------------------------------------------*
 |  apply inverse of operator (public)                       m.gee 11/04|
 *----------------------------------------------------------------------*/
int ML_NOX::ML_Nox_Preconditioner::ApplyInverse(
                                        const Epetra_MultiVector& X, 
                                        Epetra_MultiVector& Y) const
{
  int err=1;
  if (isinit_==false)
  {
    cout << "**ERR**: ML_Nox_Preconditioner::ApplyInverse:\n";
    cout << "**ERR**: tried to call ApplyInverse with isinit_ == false !\n";
    cout << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; throw -1;
  }
  if (islinearPrec_ == true)
  {
     err = ML_Nox_ApplyInverse_Linear(X,Y);
     if (err==0)
       return(0);
     else
       return(-1);
  }
  else if (islinearPrec_ == false)
  {
     err = ML_Nox_ApplyInverse_NonLinear(X,Y);
     if (err==0)
       return(0);
     else
       return(-1);
  }
  cout << "**ERR**: ML_Nox_Preconditioner::ApplyInverse:\n";
  cout << "**ERR**: something wrong with the flags ismatrixfree_ and islinearPrec_\n";
  cout << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; throw -1;
  return(-1);
}

/*----------------------------------------------------------------------*
 |  apply preconditioner                         (private)   m.gee 11/04|
 |  linear, with Jacobian supplied                                      |
 *----------------------------------------------------------------------*/
int ML_NOX::ML_Nox_Preconditioner::ML_Nox_ApplyInverse_Linear(
                                              const Epetra_MultiVector& X, 
                                              Epetra_MultiVector& Y) const
{
   int i,err;
   if (ml_linPrec_==NULL)
   {
     cout << "**ERR**: ML_Nox_Preconditioner::ML_Nox_ApplyInverse_Linear:\n";
     cout << "**ERR**: ptr to ML_Epetra::MultiLevelOperator is NULL\n";
     cout << "**ERR**: file/line: " << __FILE__ << "(" << __LINE__ << ")\n"; throw -1;
   }
   // apply ML
   err = ml_linPrec_->ApplyInverse(X,Y);   
   return(err);
}

/*----------------------------------------------------------------------*
 |  use transpose (public)                                   m.gee 11/04|
 *----------------------------------------------------------------------*/
int ML_NOX::ML_Nox_Preconditioner::SetUseTranspose(bool UseTranspose)
{
  if (UseTranspose==true) usetranspose_=true;
  else                    usetranspose_=false;
  return(-1); //  the implementation does not support use of transpose
}

/*----------------------------------------------------------------------*
 |  make a deep copy of a graph                              m.gee 01/05|
 |  allocate the new graph                                              |
 *----------------------------------------------------------------------*/
Epetra_CrsGraph* ML_NOX::ML_Nox_Preconditioner::deepcopy_graph(const Epetra_CrsGraph* oldgraph)
{
   int  i,j,ierr;
   int  nrows = oldgraph->NumMyRows();
   int* nIndicesperRow = new int[nrows];

   for (i=0; i<nrows; i++)
      nIndicesperRow[i] = oldgraph->NumMyIndices(i);
   Epetra_CrsGraph* graph = new Epetra_CrsGraph(Copy,oldgraph->RowMap(),oldgraph->ColMap(),
                                                &(nIndicesperRow[0]));
   delete [] nIndicesperRow;
   nIndicesperRow = 0;
   
   for (i=0; i<nrows; i++)
   {
      int  numIndices;
      int* Indices=0;
      ierr = oldgraph->ExtractMyRowView(i,numIndices,Indices);
      ierr = graph->InsertMyIndices(i,numIndices,Indices);
   }

   graph->TransformToLocal();
   return graph;
}                                                     


#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) 
