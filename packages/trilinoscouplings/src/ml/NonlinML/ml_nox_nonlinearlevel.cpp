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
// ML-headers
#include "ml_common.h"
#include "ml_include.h"
#include "ml_agg_VBMETIS.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

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
                          const Epetra_Comm& comm,  const Epetra_Vector& xfine, 
                          bool ismatrixfree, bool matfreelev0, bool isnlnCG,
                          int nitersCG, bool broyden, Epetra_CrsMatrix* Jac, 
                          string fsmoothertype, string smoothertype, 
                          string coarsesolvetype, 
                          int nsmooth_fine, int nsmooth, int nsmooth_coarse,  
                          double conv_normF, double conv_nupdate, 
                          int conv_maxiter,int numPDE, int nullspdim) 
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
   coarseprepost_    = 0;
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
   isnlnCG_          = isnlnCG;
   azlinSys_         = 0;
   clone_            = 0;
   nitersCG_         = nitersCG;
   broyden_          = broyden;
   Broyd_            = 0;

   if (ismatrixfree_==true)
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
           << "**ERR**: ismatrixfree_==true on level " << level_ << "\n"
           << "**ERR**: in constructor for ismatrixfree_==false - case\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   
   // ------------------------------------------------------------------------
   // get the Jacobian of this level
   const Epetra_CrsGraph* graph = 0;
   // ------------------------------------------------------------------------
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
   // ------------------------------------------------------------------------
   else
   {
      // On coarse levels get Jacobian from hierarchy
      // Note: On levels>0 SmootherA_ is a real copy of the Jacobian
      int maxnnz=0;
      double cputime=0.0;
      ML_Operator2EpetraCrsMatrix(&(ml_->Amat[level_]), SmootherA_, maxnnz, 
                                  false, cputime);
      SmootherA_->OptimizeStorage();
      graph = &(SmootherA_->Graph());
   }
   
   // just to be save
   if (!SmootherA_ || !graph)
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
           << "**ERR**: Smoother==NULL on level " << level_ << "\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }

   // ------------------------------------------------------------------------
   // generate this level's coarse interface
   coarseinterface_ = new ML_NOX::Nox_CoarseProblem_Interface(
                                      fineinterface_,level_,ml_printlevel_,
                                      P,&(graph->RowMap()),nlevel_);  

   // ------------------------------------------------------------------------
   // generate this level's coarse prepostoperator
   if (level_==0)
     coarseprepost_ = new ML_NOX::Ml_Nox_CoarsePrePostOperator(*coarseinterface_,
                                                               fineinterface_);    
   
   // ------------------------------------------------------------------------
   // get the current solution to this level
   xthis_ = coarseinterface_->restrict_fine_to_this(xfine);
   

   // ------------------------------------------------------------------------
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
   
   // set the smoother
   if (level_==0)
      Set_Smoother(ml,ag,level_,nlevel,thislevel_ml_,thislevel_ag_,fsmoothertype,nsmooth_fine);

   else if (level_ != nlevel_-1) // set the smoother from the input
      Set_Smoother(ml,ag,level_,nlevel,thislevel_ml_,thislevel_ag_,smoothertype,nsmooth);

   else // set the coarse solver from the input
      Set_Smoother(ml,ag,level_,nlevel,thislevel_ml_,thislevel_ag_,coarsesolvetype,nsmooth_coarse);
  
   // create this level's preconditioner class
   ML_Epetra::MultiLevelOperator* ml_tmp = new ML_Epetra::MultiLevelOperator(
                                                thislevel_ml_,comm_,
                                                SmootherA_->OperatorDomainMap(),
                                                SmootherA_->OperatorRangeMap());
   thislevel_prec_ = new ML_NOX::ML_Nox_ConstrainedMultiLevelOperator(ml_tmp,*coarseinterface_);
   
   if (!thislevel_prec_)
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
           << "**ERR**: thislevel_prec_==NULL on level " << level_ << "\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
                                                                 
   // ------------------------------------------------------------------------
   // set up NOX on this level   
   // ------------------------------------------------------------------------
   nlParams_ = new Teuchos::ParameterList();
   Teuchos::ParameterList& printParams = nlParams_->sublist("Printing");        
   printParams.setParameter("MyPID", comm_.MyPID()); 
   printParams.setParameter("Output Precision", 14);
   printParams.setParameter("Output Processor", 0);
   if (ml_printlevel_>9)
      printParams.setParameter("Output Information",
   	                       NOX::Utils::OuterIteration + 
			       NOX::Utils::Warning);
  else if (ml_printlevel_>8)
      printParams.setParameter("Output Information",
			       NOX::Utils::Warning);
  else
      printParams.setParameter("Output Information",0);

  if (level_==0)
    nlParams_->sublist("Solver Options").setParameter("User Defined Pre/Post Operator", *coarseprepost_);
  nlParams_->setParameter("Nonlinear Solver", "Line Search Based");         
  Teuchos::ParameterList& searchParams = nlParams_->sublist("Line Search");
  Teuchos::ParameterList* lsParamsptr  = 0;
  if (isnlnCG_)
  {
     searchParams.setParameter("Method", "NonlinearCG");
     Teuchos::ParameterList& dirParams = nlParams_->sublist("Direction"); 
     dirParams.setParameter("Method", "NonlinearCG");
     Teuchos::ParameterList& nlcgParams = dirParams.sublist("Nonlinear CG");
     nlcgParams.setParameter("Restart Frequency", 10);                         
     nlcgParams.setParameter("Precondition", "On");
     nlcgParams.setParameter("Orthogonalize", "Polak-Ribiere");
     //nlcgParams.setParameter("Orthogonalize", "Fletcher-Reeves");

     Teuchos::ParameterList& lsParams = nlcgParams.sublist("Linear Solver");     
     lsParams.setParameter("Aztec Solver", "CG"); 
     lsParams.setParameter("Max Iterations", 1);  
     lsParams.setParameter("Tolerance", 1e-11);
     lsParams.setParameter("Output Frequency", 0);   
     lsParams.setParameter("Preconditioning", "User Supplied Preconditioner");   
     lsParams.setParameter("Preconditioner","User Defined");
  }
  else // Newton's method using ML-preconditioned Aztec as linear solver
  {
     searchParams.setParameter("Method", "Full Step");
     // Sublist for direction
     Teuchos::ParameterList& dirParams = nlParams_->sublist("Direction");
     dirParams.setParameter("Method", "Newton");
     Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
     newtonParams.setParameter("Forcing Term Method", "Constant");
     //newtonParams.setParameter("Forcing Term Method", "Type 1");
     //newtonParams.setParameter("Forcing Term Method", "Type 2");
     newtonParams.setParameter("Forcing Term Minimum Tolerance", 1.0e-6);
     newtonParams.setParameter("Forcing Term Maximum Tolerance", 0.1);

     Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
     lsParamsptr = &lsParams;
     lsParams.setParameter("Size of Krylov Subspace", 100);
     lsParams.setParameter("Aztec Solver", "GMRES"); 
     lsParams.setParameter("Max Iterations", nitersCG_);  
     lsParams.setParameter("Tolerance", conv_normF_); // FIXME? is this correct?
     if (ml_printlevel_>8)
        lsParams.setParameter("Output Frequency", 50);   
     else
        lsParams.setParameter("Output Frequency", 0);   
     lsParams.setParameter("Preconditioning", "User Supplied Preconditioner");
     lsParams.setParameter("Preconditioner","User Defined");
  }

  // create the initial guess     
  initialGuess_ = new NOX::Epetra::Vector(*xthis_, NOX::DeepCopy, true);
  // NOTE: do not delete xthis_, it's used and destroyed by initialGuess_

   
  // create the necessary interfaces
  NOX::EpetraNew::Interface::Preconditioner* iPrec = 0; 
  NOX::EpetraNew::Interface::Required*       iReq  = 0;
  NOX::EpetraNew::Interface::Jacobian*       iJac  = 0;

  if (isnlnCG_)
  {
     // create the matrixfree operator used in the nlnCG
     thislevel_A_ = new NOX::EpetraNew::MatrixFree(*coarseinterface_,*xthis_,false);
  
     // create the necessary interfaces
     iPrec = 0; 
     iReq  = coarseinterface_;
     iJac  = thislevel_A_;
  
     // create the linear system 
     thislevel_linSys_ = new ML_NOX::Ml_Nox_LinearSystem(
                                    *iJac,*thislevel_A_,*iPrec,
                                    coarseinterface_,*thislevel_prec_,
                                    *xthis_,ismatrixfree_,level_,ml_printlevel_);
     // create the group
     group_ = new NOX::EpetraNew::Group(printParams,*iReq,*initialGuess_,*thislevel_linSys_);
  }
  else // Modified Newton's method
  {
     if (!broyden_)
     {
       // create the necessary interfaces   
       iPrec = this; 
       iReq  = coarseinterface_;
       //iJac  = this;
       thislevel_A_ = new NOX::EpetraNew::MatrixFree(*coarseinterface_,*xthis_,false);
     
       // create the initial guess vector
       //clone_  = new Epetra_Vector(*xthis_);

       // create the linear system 
       //azlinSys_ = new NOX::EpetraNew::LinearSystemAztecOO(
       //                                               printParams,*lsParamsptr,
       //                                               *iJac,*SmootherA_,*iPrec,
       //                                               *thislevel_prec_,*clone_);
       azlinSys_ = new NOX::EpetraNew::LinearSystemAztecOO(
                                                      printParams,*lsParamsptr,
                                                      *thislevel_A_,*thislevel_A_,*iPrec,
                                                      *thislevel_prec_,*xthis_);
     }
     else // use a Broyden update for the Jacobian
     {
       // create the initial guess vector
       //clone_  = new Epetra_Vector(*xthis_);

       // create the necessary interfaces   
       iPrec  = this; 
       iReq   = coarseinterface_;
       Broyd_ = new NOX::EpetraNew::BroydenOperator(*nlParams_,*xthis_,
                                                    *SmootherA_,false);
       
       // create the linear system 
       azlinSys_ = new NOX::EpetraNew::LinearSystemAztecOO(
                                                   printParams,*lsParamsptr,
                                                   *Broyd_,*SmootherA_,*iPrec,
                                                   *thislevel_prec_,*xthis_);
       
     }
     // create the group
     group_ = new NOX::EpetraNew::Group(printParams,*iReq,*initialGuess_,
                                        *azlinSys_);
  }

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
                          const Epetra_Comm& comm,  const Epetra_Vector& xfine, 
                          bool ismatrixfree, bool isnlnCG, int nitersCG, bool broyden,
                          string fsmoothertype, string smoothertype, string coarsesolvetype, 
                          int nsmooth_fine, int nsmooth, int nsmooth_coarse,  
                          double conv_normF, 
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
   coarseprepost_    = 0;
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
   isnlnCG_          = isnlnCG;
   azlinSys_         = 0;
   clone_            = 0;
   nitersCG_         = nitersCG;
   broyden_          = broyden;
   Broyd_            = 0;

#if 0   
   if (isnlnCG_==false && 
      (fsmoothertype   == "Jacobi" || 
       smoothertype    == "Jacobi" || 
       coarsesolvetype == "Jacobi" ))
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::ML_Nox_NonlinearLevel:\n"
           << "**ERR**: Modified Newton's method not supported for \n"
           << "**ERR**: ismatrixfree_==true && smoothertype == Jacobi-Smoother\n"
           << "**ERR**: because no full Jacobian exists!\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
#endif   
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
   
   // ------------------------------------------------------------------------
   Mat->OptimizeStorage();

   // ------------------------------------------------------------------------
   // get the current solution to this level
   xthis_ = coarseinterface_->restrict_fine_to_this(xfine);

   // ------------------------------------------------------------------------
   // create this level's preconditioner
   // We use a 1-level ML-hierarchy for that
   ML_Aggregate_Create(&thislevel_ag_);
   ML_Create(&thislevel_ml_,1);
   
   // ------------------------------------------------------------------------
   // set the Jacobian on level 0 of the local ml
   EpetraMatrix2MLMatrix(thislevel_ml_,0,
                         (dynamic_cast<Epetra_RowMatrix*>(Mat)));
   
   // ------------------------------------------------------------------------
   // construct a 1-level ML-hierarchy on this level as a smoother   
   // ------------------------------------------------------------------------
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
   
   
   // set the smoother
   if (level_==0)
      Set_Smoother(ml,ag,level_,nlevel,thislevel_ml_,thislevel_ag_,fsmoothertype,nsmooth_fine);

   else if (level_ != nlevel_-1) // set the smoother from the input
      Set_Smoother(ml,ag,level_,nlevel,thislevel_ml_,thislevel_ag_,smoothertype,nsmooth);

   else // set the coarse solver from the input
      Set_Smoother(ml,ag,level_,nlevel,thislevel_ml_,thislevel_ag_,coarsesolvetype,nsmooth_coarse);
  
   // create this level's preconditioner class
   ML_Epetra::MultiLevelOperator* ml_tmp = new ML_Epetra::MultiLevelOperator(
                                                       thislevel_ml_,comm_,
                                                       Mat->OperatorDomainMap(),
                                                       Mat->OperatorRangeMap());
   thislevel_prec_ = new ML_NOX::ML_Nox_ConstrainedMultiLevelOperator(ml_tmp,*coarseinterface_);

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

   // ------------------------------------------------------------------------
   // generate this level's coarse prepostoperator
   if (level_==0)
      coarseprepost_ = new ML_NOX::Ml_Nox_CoarsePrePostOperator(*coarseinterface_,
                                                                fineinterface_);    

   // ------------------------------------------------------------------------
   // set up NOX on this level   
   // ------------------------------------------------------------------------
   nlParams_ = new Teuchos::ParameterList();
   Teuchos::ParameterList& printParams = nlParams_->sublist("Printing");        
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

  if (level_==0)
    nlParams_->sublist("Solver Options").setParameter("User Defined Pre/Post Operator", *coarseprepost_);
  nlParams_->setParameter("Nonlinear Solver", "Line Search Based");         
  Teuchos::ParameterList& searchParams = nlParams_->sublist("Line Search");
  Teuchos::ParameterList* lsParamsptr  = 0;
  if (isnlnCG_)
  {
     searchParams.setParameter("Method", "NonlinearCG");
     Teuchos::ParameterList& dirParams = nlParams_->sublist("Direction"); 
     dirParams.setParameter("Method", "NonlinearCG");
     Teuchos::ParameterList& nlcgParams = dirParams.sublist("Nonlinear CG");
     nlcgParams.setParameter("Restart Frequency", 10);                         
     nlcgParams.setParameter("Precondition", "On");
     nlcgParams.setParameter("Orthogonalize", "Polak-Ribiere");
     //nlcgParams.setParameter("Orthogonalize", "Fletcher-Reeves");

     Teuchos::ParameterList& lsParams = nlcgParams.sublist("Linear Solver");     
     lsParams.setParameter("Aztec Solver", "CG"); 
     lsParams.setParameter("Max Iterations", 1);  
     lsParams.setParameter("Tolerance", 1e-11);
     lsParams.setParameter("Output Frequency", 0);   
     lsParams.setParameter("Preconditioning", "User Supplied Preconditioner");   
     lsParams.setParameter("Preconditioner","User Defined");
  }
  else // Newton's method using ML-preconditioned Aztec as linear solver
  {
     searchParams.setParameter("Method", "Full Step");
     // Sublist for direction
     Teuchos::ParameterList& dirParams = nlParams_->sublist("Direction");
     dirParams.setParameter("Method", "Newton");
     Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
     newtonParams.setParameter("Forcing Term Method", "Constant");
     //newtonParams.setParameter("Forcing Term Method", "Type 1");
     //newtonParams.setParameter("Forcing Term Method", "Type 2");
     newtonParams.setParameter("Forcing Term Minimum Tolerance", 1.0e-6);
     newtonParams.setParameter("Forcing Term Maximum Tolerance", 0.1);

     Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
     lsParamsptr = &lsParams;
     lsParams.setParameter("Aztec Solver", "CG"); 
     lsParams.setParameter("Max Iterations", nitersCG_);  
     lsParams.setParameter("Tolerance", conv_normF_); // FIXME? is this correct?
     if (ml_printlevel_>8)
        lsParams.setParameter("Output Frequency", 50);   
     else
        lsParams.setParameter("Output Frequency", 0);   
     lsParams.setParameter("Preconditioning", "User Supplied Preconditioner");
     lsParams.setParameter("Preconditioner","User Defined");
  }
   
  // create the initial guess   
  initialGuess_ = new NOX::Epetra::Vector(*xthis_, NOX::DeepCopy, true);
  // NOTE: do not delete xthis_, it's used and destroyed by initialGuess_

  // create the necessary interfaces
  NOX::EpetraNew::Interface::Preconditioner* iPrec = 0; 
  NOX::EpetraNew::Interface::Required*       iReq  = 0;
  NOX::EpetraNew::Interface::Jacobian*       iJac  = 0;

  if (isnlnCG_)
  {
     // create the matrixfree operator used in the nlnCG
     thislevel_A_ = new NOX::EpetraNew::MatrixFree(*coarseinterface_,*xthis_,false);
  
     // create the necessary interfaces
     iPrec = 0; 
     iReq  = coarseinterface_;
     iJac  = thislevel_A_;
  
     // create the linear system 
     thislevel_linSys_ = new ML_NOX::Ml_Nox_LinearSystem(
                                    *iJac,*thislevel_A_,*iPrec,
                                    coarseinterface_,*thislevel_prec_,
                                    *xthis_,ismatrixfree_,level_,ml_printlevel_);

     // create the group
     group_ = new NOX::EpetraNew::Group(printParams,*iReq,*initialGuess_,*thislevel_linSys_);
  }
  else // Modified Newton's method
  {
     if (!broyden_)
     {
       // create the necessary interfaces   
       iPrec = this; 
       iReq  = coarseinterface_;
       iJac  = this;
     
       // create the initial guess vector
       clone_  = new Epetra_Vector(*xthis_);

       // create the linear system 
       azlinSys_ = new NOX::EpetraNew::LinearSystemAztecOO(
                                                      printParams,*lsParamsptr,
                                                      *iJac,*SmootherA_,*iPrec,
                                                      *thislevel_prec_,*clone_);
     }
     else
     {
       // create the initial guess vector
       clone_  = new Epetra_Vector(*xthis_);

       // create the necessary interfaces   
       iPrec = this; 
       iReq  = coarseinterface_;
       Broyd_ = new NOX::EpetraNew::BroydenOperator(*nlParams_,*clone_,
                                                    *SmootherA_,false);
     
       // create the linear system 
       azlinSys_ = new NOX::EpetraNew::LinearSystemAztecOO(
                                                   printParams,*lsParamsptr,
                                                   *Broyd_,*SmootherA_,*iPrec,
                                                   *thislevel_prec_,*clone_);
     }
     // create the group
     group_ = new NOX::EpetraNew::Group(printParams,*iReq,*initialGuess_,*azlinSys_);
  }


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
      
   if (coarseprepost_) 
     delete coarseprepost_;
   coarseprepost_ = 0;
   
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
   
   if (azlinSys_)
      delete azlinSys_;
   azlinSys_ = 0;
   
   if (clone_)
      delete clone_;
   clone_ = 0;
   
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
   
   if (Broyd_)
      delete Broyd_;
   Broyd_ = 0;
 
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
  double norm = conv_normF_;
  return iterate(f,x,numiter,&norm);
}
/*----------------------------------------------------------------------*
 |  (public)                                                 m.gee 01/05|
 |  iterate on f with initial guess x numiter times                     |
 |  returns absolut value of result in x                                |
 |  returns f matching x                                                |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_NonlinearLevel::iterate(Epetra_Vector* f, Epetra_Vector* x, 
                                            int numiter, double* norm)
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
   create_Nox_Convergencetest(*norm,*norm,numiter);
  
   // make a soft reset of the NOX solver (see NOX manual)
   solver_->reset(*group_,*combo2_);
   
   // iterate
   if (ml_printlevel_ > 0 && comm_.MyPID() == 0 && coarseinterface_->isFAS()==false)
   {
      printf("ML (level %d): Entering Nonlinear Smoother, Goal: %12.8e\n",level_,*norm); fflush(stdout);
      //cout << "ML (level " << level_ << "): Entering Nonlinear Smoother, Goal: " << *norm << "\n"; fflush(stdout);
   }
   else if (ml_printlevel_ > 0 && comm_.MyPID() == 0 && coarseinterface_->isFAS()==true)
   {
      printf("ML (level %d): Entering FAS-Nonlinear Smoother, Goal: %12.8e\n",level_,*norm); fflush(stdout);
      //cout << "ML (level " << level_ << "): Entering FAS-Nonlinear Smoother, Goal: " << *norm << "\n"; fflush(stdout);
   }

   NOX::StatusTest::StatusType status = solver_->solve();
   
   // get number of iterations done
   int niter = solver_->getNumIterations();
   // get the solution and the new F
   const NOX::EpetraNew::Group& finalGroup = 
   dynamic_cast<const NOX::EpetraNew::Group&>(solver_->getSolutionGroup());
   const Epetra_Vector& finalSolution = 
   (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();      
   const Epetra_Vector& finalF = 
   (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getF())).getEpetraVector();

   bool returnstatus=false;
   double norm2;
   finalF.Norm2(&norm2);
   *norm = norm2;
      
   if (status == NOX::StatusTest::Converged)
   {
      returnstatus = true;
      if (ml_printlevel_ > 0 && comm_.MyPID() == 0) {
         //cout << "ML (level " << level_ << "): NOX: " 
         //     << niter << " iterations, Norm(F)=" 
         //     << norm2 << " , Converged\n"; fflush(stdout);
         printf("ML (level %d): NOX: %d iterations, Norm(F)=%12.8e , Converged\n",level_,niter,norm2);
         fflush(stdout);
      }
   }
   else if (status == NOX::StatusTest::Unconverged)
   {
      returnstatus = false;
      if (ml_printlevel_ > 0 && comm_.MyPID() == 0) {
         //cout << "ML (level " << level_ << "): NOX: "
         //     << niter << " iterations, Norm(F)=" 
         //     << norm2 << ", Unconverged\n"; fflush(stdout);
         printf("ML (level %d): NOX: %d iterations, Norm(F)=%12.8e , Unonverged\n",level_,niter,norm2);
         fflush(stdout);
      }
   }
   else if (status == NOX::StatusTest::Failed)
   {
      returnstatus = false;
      if (ml_printlevel_ > 0 && comm_.MyPID() == 0) {
         //cout << "ML (level " << level_ << "): NOX: " 
         //     << niter << " iterations, Norm(F)=" << norm2 << ", Failed\n"; fflush(stdout);
         printf("ML (level %d): NOX: %d iterations, Norm(F)=%12.8e , Failed\n",level_,niter,norm2);
         fflush(stdout);
      }
   }
   else
   {
      returnstatus = false;
      if (comm_.MyPID() == 0)
         //cout << "ML (level " << level_ << "): ***WRN*** NOX returned unknown status, Norm(F)=" 
         //     << norm2 << "\n"; fflush(stdout);
         printf("ML (level %d): ***WRN*** NOX: return status unknown, Norm(F)=%12.8e , Failed\n",level_,norm2);
         fflush(stdout);
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
 |  (private)                                                m.gee 04/05|
 | set the smoother on this nonlinear level                             |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_NonlinearLevel::Set_Smoother(ML*           ml, 
                                                 ML_Aggregate* ag, 
                                                 int           level,
                                                 int           nlevel,
                                                 ML*           thislevel_ml,
                                                 ML_Aggregate* thislevel_ag, 
                                                 string        smoothertype, 
                                                 int           nsmooth)
{
   if (smoothertype == "SGS")
      ML_Gen_Smoother_SymGaussSeidel(thislevel_ml,0,ML_BOTH,nsmooth,1.0);
   else if (smoothertype == "Jacobi")
      ML_Gen_Smoother_Jacobi(thislevel_ml,0,ML_BOTH,nsmooth,0.25);
   else if (smoothertype == "AmesosKLU")
      ML_Gen_Smoother_Amesos(thislevel_ml,0,ML_AMESOS_KLU,-1,0.0,1);
   else if ( (smoothertype == "MLS") || (smoothertype == "Cheby") )
      ML_Gen_Smoother_Cheby(thislevel_ml,0,ML_BOTH,30.,nsmooth);
   else if (smoothertype == "BSGS")
   {
      int  nblocks  = 0;
      int* blocks   = NULL;
      int* blockpde = NULL;
      bool needfree = false;

      // try to get nodal blocks from the VBMETIS aggregation scheme
      ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS(ag,level,nlevel,
                                                  &nblocks,&blocks,&blockpde);

      if (nblocks && blocks)
         needfree=true;
      else
         ML_Gen_Blocks_Aggregates(ag,level,&nblocks,&blocks);

      ML_Gen_Smoother_VBlockSymGaussSeidel(thislevel_ml,0,ML_BOTH,nsmooth,1.,
                                           nblocks,blocks);
      if (needfree)
      {
         ML_free(blocks); 
         ML_free(blockpde);
      }
   }
   else if (smoothertype == "Bcheby")
   {
      int  nblocks  = 0;
      int* blocks   = NULL;
      int* blockpde = NULL;
      bool needfree = false;

      // try to get nodal blocks from the VBMETIS aggregation scheme
      ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS(ag,level,nlevel,
                                                  &nblocks,&blocks,&blockpde);

      if (nblocks && blocks)
         needfree=true;
      else
         ML_Gen_Blocks_Aggregates(ag,level,&nblocks,&blocks);

      ML_Gen_Smoother_BlockDiagScaledCheby(thislevel_ml,0,ML_BOTH,30.,nsmooth,
                                           nblocks,blocks);
      if (needfree)
      {
         ML_free(blocks); 
         ML_free(blockpde);
      }
   }
   else
   {
      cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::Setsmoother:\n"
           << "**ERR**: unknown type of smoother: " <<  smoothertype << "\n" 
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
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
   
   absresid_ = new NOX::StatusTest::NormF(conv_normF,NOX::StatusTest::NormF::Unscaled); 
   
   nupdate_ = new NOX::StatusTest::NormUpdate(conv_nupdate,NOX::StatusTest::NormUpdate::Unscaled);
   
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
			                                  Teuchos::ParameterList* precParams)
{
#if 0
   cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::computePreconditioner:\n"
        << "**ERR**: not impl. on level " << level_ << "\n"
        << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
#endif   
   return true;
}

/*----------------------------------------------------------------------*
 |                      (public,derived)                     m.gee 03/05|
 | compute this level's Jacobian                                        |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_NonlinearLevel::computeJacobian(const Epetra_Vector& x)
{
#if 0
   cout << "**ERR**: ML_NOX::ML_Nox_NonlinearLevel::computeJacobian:\n"
        << "**ERR**: not impl. on level " << level_ << "\n"
        << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
#endif
   return true;
}
#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) 
