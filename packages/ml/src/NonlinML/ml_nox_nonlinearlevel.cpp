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

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) 

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#ifdef ML_MPI // FIXME: do I need this?
#include <mpi.h>
#endif

// this class
#include "ml_nox_nonlinearlevel.H"

/*----------------------------------------------------------------------*
 |  the class defining a nonlinear coarse level              m.gee 01/05|
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  Constructor (public)                                     m.gee 01/05|
 |  IMPORTANT:                                                          |
 |  No matter on which level we are here, the vector xfine is ALWAYS    |
 |  a fine grid vector here!                                            |
 |  this is the constructor for the ismatrixfree==false case
 *----------------------------------------------------------------------*/
ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel(
                          int level, int nlevel, int printlevel, ML* ml, 
                          ML_Aggregate* ag,Epetra_CrsMatrix** P, 
                          ML_NOX::Ml_Nox_Fineinterface& interface,
                          Epetra_Comm& comm,  const Epetra_Vector& xfine, 
                          bool ismatrixfree, bool matfreelev0, Epetra_CrsMatrix* Jac,
                          string fsmoothertype, string smoothertype, string coarsesolvetype, 
                          int *nsmooth, double conv_normF, 
                          double conv_nupdate, int conv_maxiter,
                          int numPDE, int nullspdim) 
: fineinterface_(interface),
  comm_(comm)
{
   level_            = level;        // this level
   nlevel_           = nlevel;       // number of total levels
   ml_printlevel_    = printlevel;   // printlevel
   ml_               = ml;           // the global ML object
   ag_               = ag;           // the global ML_Aggregate object
   thislevel_prec_   = 0;            // this level's linear preconditioner
   thislevel_ml_     = 0;            // this level's local ML object 
   thislevel_ag_     = 0;            // this level's local ML_Aggregate object
   coarseinterface_  = 0;            // this level's coarse interface
   xthis_            = 0;            // this level's current solution matching this level's map!!!!
   thislevel_A_      = 0;            // this level's NOX Matrixfree operator
   SmootherA_        = 0;            // this level's Epetra_CrsMatrix for thislevel_prec_
   ismatrixfree_     = ismatrixfree; // matrixfree flag
   conv_normF_       = conv_normF;   // NOX convergence test stuff
   conv_nupdate_     = conv_nupdate;
   conv_maxiter_     = conv_maxiter;
   absresid_         = 0;            
   nupdate_          = 0;
   fv_               = 0;
   maxiters_         = 0;
   combo1_           = 0;
   combo2_           = 0;
   thislevel_linSys_ = 0;            // this level's NOX linear system
   nlParams_         = 0;            // NOX parameters
   initialGuess_     = 0;            // NOX initial guess
   group_            = 0;            // NOX group
   solver_           = 0;            // NOX solver

   if (ismatrixfree_==true)
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
           << "**ERR**: ismatrixfree_==true on level " << level_ << "\n"
           << "**ERR**: in constructor for ismatrixfree_==false - case\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   
   // get the Jacobian of this level
   const Epetra_CrsGraph* graph = 0;
   if (level_==0)
   {
      graph = fineinterface_.getGraph();
      // On fine level this is the fineinterface's Jacobian
      if (matfreelev0==false)
         SmootherA_ = fineinterface_.getJacobian();
      else if (matfreelev0==true && Jac)
         SmootherA_ = Jac;
      else
      {
        cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
             << "**ERR**: something weired happened\n"
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
   }
   else
   {
      // On coarse levels get Jacobian from hierarchy
      // Note: On levels>0 SmootherA_ is a real copy of the Jacobian
      int maxnnz=0;
      double cputime=0.0;
      ML_Operator2EpetraCrsMatrix(&(ml_->Amat[level_]), SmootherA_, maxnnz, 
                                  false, cputime);
      graph = &(SmootherA_->Graph());
   }
   
   // just to be save
   if (!SmootherA_ || !graph)
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
           << "**ERR**: Smoother==NULL on level " << level_ << "\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }

   // generate this level's coarse interface
   coarseinterface_ = new ML_NOX::Nox_CoarseProblem_Interface(
                                      fineinterface_,level_,ml_printlevel_,
                                      P,&(graph->RowMap()),nlevel_);  

   // get the current solution to this level
   xthis_ = coarseinterface_->restrict_fine_to_this(xfine);
   

   // create this level's preconditioner
   // We use a 1-level ML-hierarchy for that
   ML_Aggregate_Create(&thislevel_ag_);
   ML_Create(&thislevel_ml_,1);
   
   // set the Jacobian on level 0 of the local ml
   EpetraMatrix2MLMatrix(thislevel_ml_,0,
                         (dynamic_cast<Epetra_RowMatrix*>(SmootherA_)));
   
   // construct a 1-level ML-hierarchy on this level as a smoother   
   ML_Set_PrintLevel(ml_printlevel_);  
   ML_Aggregate_Set_CoarsenScheme_Uncoupled(thislevel_ag_); 
   ML_Aggregate_Set_DampingFactor(thislevel_ag_, 0.0);
   ML_Aggregate_Set_Threshold(thislevel_ag_, 0.0);
   ML_Aggregate_Set_MaxCoarseSize(thislevel_ag_,1);
   ML_Aggregate_Set_NullSpace(thislevel_ag_,numPDE,nullspdim,NULL,
                              SmootherA_->NumMyRows());
   int thislevel_nlevel = ML_Gen_MGHierarchy_UsingAggregation(thislevel_ml_,0,
                                                   ML_INCREASING,thislevel_ag_);
   if (thislevel_nlevel != 1)
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
           << "**ERR**: ML generated a local hierarchy of " <<  thislevel_nlevel << " on level " << level_ << "\n" 
           << "**ERR**: this is supposed to be 1 Level only!\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   
   // set the preconditioner
   if (level_==0)
   {
      if (fsmoothertype == "SGS")
         ML_Gen_Smoother_SymGaussSeidel(thislevel_ml_,0,ML_BOTH,
                                        nsmooth[level_],0.6);
      else if (fsmoothertype == "Jacobi")
         ML_Gen_Smoother_Jacobi(thislevel_ml_,0,ML_BOTH, 
                                nsmooth[level_],0.25);
      else if (fsmoothertype == "AmesosKLU")
         ML_Gen_Smoother_Amesos(thislevel_ml_,0,ML_AMESOS_KLU,-1);
      else
      {
         cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
              << "**ERR**: unknown type of fine smoother: " <<  fsmoothertype << " on level " << level_ << "\n" 
              << "**ERR**: implemented are 'SGS' 'Jacobi' 'AmesosKLU'\n"
              << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
   }
   else if (level_ != nlevel_-1) // set the smoother from the input
   {
      if (smoothertype == "SGS")
         ML_Gen_Smoother_SymGaussSeidel(thislevel_ml_,0,ML_BOTH,
                                        nsmooth[level_],0.6);
      else if (smoothertype == "Jacobi")
         ML_Gen_Smoother_Jacobi(thislevel_ml_,0,ML_BOTH, 
                                nsmooth[level_],0.25);
      else if (smoothertype == "AmesosKLU")
         ML_Gen_Smoother_Amesos(thislevel_ml_,0,ML_AMESOS_KLU,-1);
      else
      {
         cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
              << "**ERR**: unknown type of smoother: " <<  smoothertype << " on level " << level_ << "\n" 
              << "**ERR**: implemented are 'SGS' 'Jacobi' 'AmesosKLU'\n"
              << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
   }
   else // set the coarse solver from the input
   {
      if (coarsesolvetype == "SGS")
         ML_Gen_Smoother_SymGaussSeidel(thislevel_ml_,0,ML_BOTH,
                                        nsmooth[level_],0.6);
      else if (coarsesolvetype == "Jacobi")
         ML_Gen_Smoother_Jacobi(thislevel_ml_,0,ML_BOTH, 
                                nsmooth[level_],.25);
      else if (coarsesolvetype == "AmesosKLU")
         ML_Gen_Smoother_Amesos(thislevel_ml_,0,ML_AMESOS_KLU,-1);
      else
      {
         cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
              << "**ERR**: unknown type of coarsesolve: " <<  coarsesolvetype << " on level " << level_ << "\n" 
              << "**ERR**: implemented are 'SGS' 'Jacobi' 'AmesosKLU'\n"
              << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
   }
   
  
   // create this level's preconditioner class
   thislevel_prec_ = new ML_Epetra::MultiLevelOperator(
                                                thislevel_ml_,comm_,
                                                SmootherA_->OperatorDomainMap(),
                                                SmootherA_->OperatorRangeMap());
   if (!thislevel_prec_)
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
           << "**ERR**: thislevel_prec_==NULL on level " << level_ << "\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
                                                                 
   // intensive test of this level's ML-smoother
#if 0
   {
   cout << "Test of smoother on level " << level_ << endl;
   Epetra_Vector *out = new Epetra_Vector(Copy,*xthis_,0);
   out->PutScalar(0.0);
   cout << "Input\n";
   xthis_->PutScalar(1.0);
   SmootherA_->Multiply(false,*xthis_,*out);
   xthis_->PutScalar(3.0);
   cout << "rhs\n";
   cout << *out;
   double norm = 0.0;
   out->Norm1(&norm);
   cout << "Norm of rhs = " << norm << endl;
   thislevel_prec_->ApplyInverse(*out,*xthis_);
   cout << "result after smoother\n";
   cout << *xthis_;
   delete out; out = 0;
   }
   if (level_==2) exit(0);
#endif   

   // set up NOX's nlnCG on this level   
   nlParams_ = new NOX::Parameter::List();
   NOX::Parameter::List& printParams = nlParams_->sublist("Printing");        
   printParams.setParameter("MyPID", comm_.MyPID()); 
   printParams.setParameter("Output Precision", 9);
   printParams.setParameter("Output Processor", 0);
   if (ml_printlevel_>9)
      printParams.setParameter("Output Information",
   	                       NOX::Utils::OuterIteration + 
			       //NOX::Utils::OuterIterationStatusTest + 
			       //NOX::Utils::InnerIteration +
			       //NOX::Utils::Parameters + 
			       //NOX::Utils::Details + 
			       NOX::Utils::Warning);
  else if (ml_printlevel_>8)
      printParams.setParameter("Output Information",
			       NOX::Utils::Warning);
  else
      printParams.setParameter("Output Information",0);

  nlParams_->setParameter("Nonlinear Solver", "Line Search Based");         
  NOX::Parameter::List& searchParams = nlParams_->sublist("Line Search");
  searchParams.setParameter("Method", "NonlinearCG");
  NOX::Parameter::List& dirParams = nlParams_->sublist("Direction"); 
  dirParams.setParameter("Method", "NonlinearCG");
  NOX::Parameter::List& nlcgParams = dirParams.sublist("Nonlinear CG");
  nlcgParams.setParameter("Restart Frequency", 10);                         
  nlcgParams.setParameter("Precondition", "On");
  nlcgParams.setParameter("Orthogonalize", "Polak-Ribiere");
  //nlcgParams.setParameter("Orthogonalize", "Fletcher-Reeves");

  NOX::Parameter::List& lsParams = nlcgParams.sublist("Linear Solver");     
  lsParams.setParameter("Aztec Solver", "CG"); 
  lsParams.setParameter("Max Iterations", 1);  
  lsParams.setParameter("Tolerance", 1e-11);
  lsParams.setParameter("Output Frequency", 50);   
  lsParams.setParameter("Preconditioning", "User Supplied Preconditioner");   
  lsParams.setParameter("Preconditioner","User Defined");
   
  // create the matrixfree operator used in the nlnCG
  thislevel_A_ = new NOX::EpetraNew::MatrixFree(*coarseinterface_,*xthis_,false);
  
  // create the necessary interfaces
  NOX::EpetraNew::Interface::Preconditioner* iPrec = 0; 
  NOX::EpetraNew::Interface::Required*       iReq  = coarseinterface_;
  NOX::EpetraNew::Interface::Jacobian*       iJac  = thislevel_A_;
  
  // create the linear system 
  thislevel_linSys_ = new ML_NOX::Ml_Nox_LinearSystem(
                                 *iJac,*thislevel_A_,*iPrec,*thislevel_prec_,
                                 *xthis_,ismatrixfree_,level_,ml_printlevel_);

  // create the initial guess     
  initialGuess_ = new NOX::Epetra::Vector(*xthis_, NOX::DeepCopy, true);
  // NOTE: do not delete xthis_, it's used and destroyed by initialGuess_

  // create the group
  group_ = new NOX::EpetraNew::Group(printParams,*iReq,*initialGuess_,*thislevel_linSys_);

  // create convergence test
  create_Nox_Convergencetest(conv_normF_,conv_nupdate_,conv_maxiter_);

  // create the solver
  solver_ = new NOX::Solver::Manager(*group_,*combo2_,*nlParams_);

  return;
}

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     m.gee 01/05|
 |  IMPORTANT:                                                          |
 |  No matter on which level we are here, the vector xfine is ALWAYS    |
 |  a fine grid vector here!                                            |
 |  this is the constructor for the ismatrixfree==true case
 *----------------------------------------------------------------------*/
ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel(
                          int level, int nlevel, int printlevel, ML* ml, 
                          ML_Aggregate* ag,Epetra_CrsMatrix** P, 
                          ML_NOX::Ml_Nox_Fineinterface& interface,
                          Epetra_Comm& comm,  const Epetra_Vector& xfine, 
                          bool ismatrixfree, 
                          string fsmoothertype, string smoothertype, string coarsesolvetype, 
                          int *nsmooth, double conv_normF, 
                          double conv_nupdate, int conv_maxiter,
                          int numPDE, int nullspdim, Epetra_CrsMatrix* Mat,
                          ML_NOX::Nox_CoarseProblem_Interface* coarseinterface) 
: fineinterface_(interface),
  comm_(comm)
{
   level_            = level;           // this level
   nlevel_           = nlevel;          // number of total levels
   ml_printlevel_    = printlevel;      // printlevel
   ml_               = ml;              // the global ML object
   ag_               = ag;              // the global ML_Aggregate object
   thislevel_prec_   = 0;               // this level's linear preconditioner
   thislevel_ml_     = 0;               // this level's local ML object 
   thislevel_ag_     = 0;               // this level's local ML_Aggregate object
   coarseinterface_  = coarseinterface; // this level's coarse interface
   xthis_            = 0;               // this level's current solution matching this level's map!!!!
   thislevel_A_      = 0;               // this level's NOX Matrixfree operator
   SmootherA_        = 0;               // this level's Epetra_CrsMatrix for thislevel_prec_
   ismatrixfree_     = ismatrixfree;    // matrixfree flag
   conv_normF_       = conv_normF;      // NOX convergence test stuff
   conv_nupdate_     = conv_nupdate;
   conv_maxiter_     = conv_maxiter;
   absresid_         = 0;            
   nupdate_          = 0;
   fv_               = 0;
   maxiters_         = 0;
   combo1_           = 0;
   combo2_           = 0;
   thislevel_linSys_ = 0;               // this level's NOX linear system
   nlParams_         = 0;               // NOX parameters
   initialGuess_     = 0;               // NOX initial guess
   group_            = 0;               // NOX group
   solver_           = 0;               // NOX solver
   SmootherA_        = Mat;
   
   if (ismatrixfree_==false)
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
           << "**ERR**: ismatrixfree_==false on level " << level_ << "\n"
           << "**ERR**: in constructor for ismatrixfree_==true - case\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   if (!coarseinterface_)
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
           << "**ERR**: ptr to coarseinterface=NULL on level " << level_ << "\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   if (!Mat)
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
           << "**ERR**: ptr to Matrix Mat=NULL on level " << level_ << "\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   
   // get the current solution to this level
   xthis_ = coarseinterface_->restrict_fine_to_this(xfine);

   // create this level's preconditioner
   // We use a 1-level ML-hierarchy for that
   ML_Aggregate_Create(&thislevel_ag_);
   ML_Create(&thislevel_ml_,1);
   
   // set the Jacobian on level 0 of the local ml
   EpetraMatrix2MLMatrix(thislevel_ml_,0,
                         (dynamic_cast<Epetra_RowMatrix*>(Mat)));
   
   // construct a 1-level ML-hierarchy on this level as a smoother   
   ML_Set_PrintLevel(ml_printlevel_);  
   ML_Aggregate_Set_CoarsenScheme_Uncoupled(thislevel_ag_); 
   ML_Aggregate_Set_DampingFactor(thislevel_ag_, 0.0);
   ML_Aggregate_Set_Threshold(thislevel_ag_, 0.0);
   ML_Aggregate_Set_MaxCoarseSize(thislevel_ag_,1);
   ML_Aggregate_Set_NullSpace(thislevel_ag_,numPDE,nullspdim,NULL,Mat->NumMyRows());
   int thislevel_nlevel = ML_Gen_MGHierarchy_UsingAggregation(thislevel_ml_,0,
                                                   ML_INCREASING,thislevel_ag_);
   if (thislevel_nlevel != 1)
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
           << "**ERR**: ML generated a local hierarchy of " <<  thislevel_nlevel << " on level " << level_ << "\n" 
           << "**ERR**: this is supposed to be 1 Level only!\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   
   // set the preconditioner
   if (level_==0)
   {
      if (fsmoothertype == "SGS")
         ML_Gen_Smoother_SymGaussSeidel(thislevel_ml_,0,ML_BOTH,
                                        nsmooth[level_],1.);
      else if (fsmoothertype == "Jacobi")
         ML_Gen_Smoother_Jacobi(thislevel_ml_,0,ML_BOTH, 
                                nsmooth[level_],0.6);
      else if (fsmoothertype == "AmesosKLU")
         ML_Gen_Smoother_Amesos(thislevel_ml_,0,ML_AMESOS_KLU,-1);
      else
      {
         cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
              << "**ERR**: unknown type of fine smoother: " <<  fsmoothertype << " on level " << level_ << "\n" 
              << "**ERR**: implemented are 'SGS' 'Jacobi' 'AmesosKLU'\n"
              << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
   }
   else if (level_ != nlevel_-1) // set the smoother from the input
   {
      if (smoothertype == "SGS")
         ML_Gen_Smoother_SymGaussSeidel(thislevel_ml_,0,ML_BOTH,
                                        nsmooth[level_],1.);
      else if (smoothertype == "Jacobi")
         ML_Gen_Smoother_Jacobi(thislevel_ml_,0,ML_BOTH, 
                                nsmooth[level_],0.6);
      else if (smoothertype == "AmesosKLU")
         ML_Gen_Smoother_Amesos(thislevel_ml_,0,ML_AMESOS_KLU,-1);
      else
      {
         cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
              << "**ERR**: unknown type of smoother: " <<  smoothertype << " on level " << level_ << "\n" 
              << "**ERR**: implemented are 'SGS' 'Jacobi' 'AmesosKLU'\n"
              << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      
   }
   else // set the coarse solver from the input
   {
      if (coarsesolvetype == "SGS")
         ML_Gen_Smoother_SymGaussSeidel(thislevel_ml_,0,ML_BOTH,
                                        nsmooth[level_],1.);
      else if (coarsesolvetype == "Jacobi")
         ML_Gen_Smoother_Jacobi(thislevel_ml_,0,ML_BOTH, 
                                nsmooth[level_],.6);
      else if (coarsesolvetype == "AmesosKLU")
         ML_Gen_Smoother_Amesos(thislevel_ml_,0,ML_AMESOS_KLU,-1);
      else
      {
         cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
              << "**ERR**: unknown type of coarsesolve: " <<  coarsesolvetype << " on level " << level_ << "\n" 
              << "**ERR**: implemented are 'SGS' 'Jacobi' 'AmesosKLU'\n"
              << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
   }
   
  
   // create this level's preconditioner class
   thislevel_prec_ = new ML_Epetra::MultiLevelOperator(thislevel_ml_,comm_,
                                                       Mat->OperatorDomainMap(),
                                                       Mat->OperatorRangeMap());
   if (!thislevel_prec_)
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
           << "**ERR**: thislevel_prec_==NULL on level " << level_ << "\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
                                                                 
   // intensive test of this level's ML-smoother
#if 0
   {
   cout << "Test of smoother on level " << level_ << endl;
   Epetra_Vector *out = new Epetra_Vector(Copy,*xthis_,0);
   out->PutScalar(0.0);
   cout << "Input\n";
   xthis_->PutScalar(1.0);
   Mat->Multiply(false,*xthis_,*out);
   xthis_->PutScalar(3.0);
   cout << "rhs\n";
   cout << *out;
   double norm = 0.0;
   out->Norm1(&norm);
   cout << "Norm of rhs = " << norm << endl;
   thislevel_prec_->ApplyInverse(*out,*xthis_);
   cout << "result after smoother\n";
   cout << *xthis_;
   delete out; out = 0;
   }
   if (level_==2) exit(0);
#endif   

   // set up NOX's nlnCG on this level   
   nlParams_ = new NOX::Parameter::List();
   NOX::Parameter::List& printParams = nlParams_->sublist("Printing");        
   printParams.setParameter("MyPID", comm_.MyPID()); 
   printParams.setParameter("Output Precision", 9);
   printParams.setParameter("Output Processor", 0);
   if (ml_printlevel_>9)
      printParams.setParameter("Output Information",
   	                       NOX::Utils::OuterIteration + 
			       //NOX::Utils::OuterIterationStatusTest + 
			       //NOX::Utils::InnerIteration +
			       //NOX::Utils::Parameters + 
			       //NOX::Utils::Details + 
			       NOX::Utils::Warning);
  else if (ml_printlevel_>8)
      printParams.setParameter("Output Information",
			       NOX::Utils::Warning);
  else
      printParams.setParameter("Output Information",0);

  nlParams_->setParameter("Nonlinear Solver", "Line Search Based");         
  NOX::Parameter::List& searchParams = nlParams_->sublist("Line Search");
  searchParams.setParameter("Method", "NonlinearCG");
  NOX::Parameter::List& dirParams = nlParams_->sublist("Direction"); 
  dirParams.setParameter("Method", "NonlinearCG");
  NOX::Parameter::List& nlcgParams = dirParams.sublist("Nonlinear CG");
  nlcgParams.setParameter("Restart Frequency", 50);                         
  nlcgParams.setParameter("Precondition", "On");
  nlcgParams.setParameter("Orthogonalize", "Polak-Ribiere");
  //nlcgParams.setParameter("Orthogonalize", "Fletcher-Reeves");

  NOX::Parameter::List& lsParams = nlcgParams.sublist("Linear Solver");     
  lsParams.setParameter("Aztec Solver", "CG"); 
  lsParams.setParameter("Max Iterations", 1);  
  lsParams.setParameter("Tolerance", 1e-11);
  lsParams.setParameter("Output Frequency", 50);   
  lsParams.setParameter("Preconditioning", "User Supplied Preconditioner");   
  lsParams.setParameter("Preconditioner","User Defined");
   
  // create the matrixfree operator used in the nlnCG
  thislevel_A_ = new NOX::EpetraNew::MatrixFree(*coarseinterface_,*xthis_,false);
  
  // create the necessary interfaces
  NOX::EpetraNew::Interface::Preconditioner* iPrec = 0; 
  NOX::EpetraNew::Interface::Required*       iReq  = coarseinterface_;
  NOX::EpetraNew::Interface::Jacobian*       iJac  = thislevel_A_;
  
  // create the linear system 
  thislevel_linSys_ = new ML_NOX::Ml_Nox_LinearSystem(
                                 *iJac,*thislevel_A_,*iPrec,*thislevel_prec_,
                                 *xthis_,ismatrixfree_,level_,ml_printlevel_);

  // create the initial guess   
  initialGuess_ = new NOX::Epetra::Vector(*xthis_, NOX::DeepCopy, true);
  // NOTE: do not delete xthis_, it's used and destroyed by initialGuess_

  // create the group
  group_ = new NOX::EpetraNew::Group(printParams,*iReq,*initialGuess_,*thislevel_linSys_);

  // create convergence test
  create_Nox_Convergencetest(conv_normF_,conv_nupdate_,conv_maxiter_);

  // create the solver
  solver_ = new NOX::Solver::Manager(*group_,*combo2_,*nlParams_);
  
  return;
}

/*----------------------------------------------------------------------*
 |  Destructor (public)                                      m.gee 01/05|
 *----------------------------------------------------------------------*/
ML_NOX::ML_Nox_NonlinearLevel::~ML_Nox_NonlinearLevel()
{

   // destroying this level's operator depends on ismatrixfree_ and on the level
   if (ismatrixfree_==false)
   {
      if (level_ == 0)
      {
         SmootherA_ = 0;
      }
      else if (SmootherA_)
      {
         delete SmootherA_;
         SmootherA_ = 0;
      }
      else
      {
        cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::~ML_Nox_NonlinearLevel:\n"
             << "**ERR**: something weird happened while destroying SmootherA_ on level " << level_ << "\n" 
             << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
   }
   else  // ismatrixfree==true;
   {
      // this is just a ptr to the finite difference constructed matrix
      // the matrixfree level is charge of destroying it
      SmootherA_ = 0;
   }

   // in the matrixfree case, the coarseinterface is owned by the matrixfree
   // class of this level and is destroyed there
   if (ismatrixfree_==false)
   {
      if (coarseinterface_)
         delete coarseinterface_;
      coarseinterface_ = 0;
   }
   else
      coarseinterface_ = 0;
   
   if (thislevel_ag_)
      ML_Aggregate_Destroy(&thislevel_ag_);
   thislevel_ag_ = 0;
   
   if (thislevel_ml_)
      ML_Destroy(&thislevel_ml_);
   thislevel_ml_ = 0;
 
   if (thislevel_prec_)
      delete thislevel_prec_;
   thislevel_prec_ = 0;
   
   if (xthis_)
      delete xthis_;
   xthis_ = 0;
   
   if (thislevel_A_)
      delete thislevel_A_;
   thislevel_A_ = 0;
   
   if (thislevel_linSys_)
      delete thislevel_linSys_;
   thislevel_linSys_ = 0;
   
   if (nlParams_)
      delete nlParams_;
   nlParams_ = 0;
   
   if (absresid_)
      delete absresid_;
   absresid_ = 0;

   if (nupdate_)
      delete nupdate_;
   nupdate_ = 0;

   if (fv_)
      delete fv_;
   fv_ = 0;

   if (maxiters_)
      delete maxiters_;
   maxiters_ = 0;

   if (combo1_)
      delete combo1_;
   combo1_ = 0;

   if (combo2_)
      delete combo2_;
   combo2_ = 0;
   
   if (group_)
      delete group_;
   group_ = 0;
   
   if (initialGuess_)
      delete initialGuess_;
   initialGuess_ = 0;
   
   if (solver_)
      delete solver_;
   solver_ = 0;
 
   return;
}

/*----------------------------------------------------------------------*
 |  (public)                                                 m.gee 01/05|
 |  iterate on f with initial guess x numiter times                     |
 |  returns absolut value of result in x                                |
 |  returns f matching x                                                |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_NonlinearLevel::iterate(Epetra_Vector* f, Epetra_Vector* x, 
                                            int numiter)
{
   if (!group_ || !solver_ || !combo2_)
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::iterate:\n"
           << "**ERR**: group_ || solver_ || combo2_ or even some of them are NULL\n" 
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   
   // put x into group_ 
   NOX::Epetra::Vector X(*x,NOX::DeepCopy,true);
   group_->setX(X);
   
   // create the convergence test
   create_Nox_Convergencetest(conv_normF_,conv_nupdate_,numiter);
  
   // make a soft reset of the NOX solver (see NOX manual)
   solver_->reset(*group_,*combo2_);
   
   // iterate
   if (ml_printlevel_ > 5 && comm_.MyPID() == 0 && coarseinterface_->isFAS()==false)
      cout << "ML (level " << level_ << "): Entering NOX iteration\n";
   else if (ml_printlevel_ > 5 && comm_.MyPID() == 0 && coarseinterface_->isFAS()==true)
      cout << "ML (level " << level_ << "): Entering FAS-NOX iteration\n";
   NOX::StatusTest::StatusType status = solver_->solve();
   
   // get the solution and the new F
   const NOX::EpetraNew::Group& finalGroup = 
   dynamic_cast<const NOX::EpetraNew::Group&>(solver_->getSolutionGroup());
   const Epetra_Vector& finalSolution = 
   (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();      
   const Epetra_Vector& finalF = 
   (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getF())).getEpetraVector();

   bool returnstatus=false;
   double norm2=0.0;
   if (ml_printlevel_ > 5)
      finalF.Norm2(&norm2);
      
   if (status == NOX::StatusTest::Converged)
   {
      returnstatus = true;
      if (ml_printlevel_ > 5 && comm_.MyPID() == 0)
         cout << "ML (level " << level_ << "): NOX Converged, Norm(F)=" 
              << norm2 << "\n";
   }
   else if (status == NOX::StatusTest::Unconverged)
   {
      returnstatus = false;
      if (ml_printlevel_ > 5 && comm_.MyPID() == 0)
         cout << "ML (level " << level_ << "): NOX Unconverged, Norm(F)=" 
              << norm2 << "\n";
   }
   else if (status == NOX::StatusTest::Failed)
   {
      returnstatus = false;
      if (ml_printlevel_ > 5 && comm_.MyPID() == 0)
         cout << "ML (level " << level_ << "): NOX Unconverged after maxiter=" 
              << numiter << ", Norm(F)=" << norm2 << "\n";
   }
   else
   {
      returnstatus = false;
      if (comm_.MyPID() == 0)
         cout << "ML (level " << level_ << "): ***WRN*** NOX returned unknown status, Norm(F)=" 
              << norm2 << "\n";
   }
   
   // reset number of calls to coarseinterface->computeF
   coarseinterface_->resetnumcallscomputeF();
   
   // update the solution
   x->Update(1.0,finalSolution,0.0);
   f->Update(1.0,finalF,0.0);
   
   return returnstatus;
}

/*----------------------------------------------------------------------*
 |  (public)                                                 m.gee 01/05|
 |  iterate on f with initial guess x numiter times                     |
 |  returns absolut value of result in x                                |
 |  returns f matching x                                                |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_NonlinearLevel::applySmoother(Epetra_Vector* f, Epetra_Vector* x, 
                                                  int numiter)
{
   if (!thislevel_prec_)
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::applySmoother:\n"
           << "**ERR**: thislevel_prec_ is NULL\n" 
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   
   int i;
   
   Epetra_Vector out(Copy,*x,0);
   
   for (i=0; i<numiter; i++)
   {
      thislevel_prec_->ApplyInverse(*f,out);
      x->Update(1.0,out,1.0);
      computeF(*x,*f,NOX::EpetraNew::Interface::Required::Residual);
   }
      
   return true;
}

/*----------------------------------------------------------------------*
 |  (private)                                                m.gee 01/05|
 | (re)create a NOX convergence test                                    |
 | deletes combo1_ and combo2_ and recreates them                       |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_NonlinearLevel::create_Nox_Convergencetest(double conv_normF, 
                                                               double conv_nupdate, 
                                                               int    conv_maxiter)
{
   if (absresid_)
      delete absresid_;
   if (nupdate_)
      delete nupdate_;
   if (fv_)
      delete fv_;
   if (maxiters_)
      delete maxiters_;
   if (combo1_)
      delete combo1_;
   if (combo2_)
      delete combo2_;
   
   absresid_ = new NOX::StatusTest::NormF(conv_normF); 
   
   nupdate_ = new NOX::StatusTest::NormUpdate(conv_nupdate);
   
   combo1_ = new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND);
   combo1_->addStatusTest(*absresid_);
   combo1_->addStatusTest(*nupdate_);
   
   fv_ = new NOX::StatusTest::FiniteValue();
   
   maxiters_ = new NOX::StatusTest::MaxIters(conv_maxiter);

   combo2_ = new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR);
   combo2_->addStatusTest(*maxiters_);
   combo2_->addStatusTest(*combo1_);
   combo2_->addStatusTest(*fv_);
    
   return true;
}

/*----------------------------------------------------------------------*
 |                      (public,derived)                     m.gee 01/05|
 | compute this level's prec to this level's nlnCG                      |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_NonlinearLevel::computePreconditioner(const Epetra_Vector& x,
			                                  NOX::Parameter::List* precParams)
{
   cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::computePreconditioner:\n"
        << "**ERR**: not impl. on level " << level_ << "\n"
        << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   
   return true;
}

#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) 
