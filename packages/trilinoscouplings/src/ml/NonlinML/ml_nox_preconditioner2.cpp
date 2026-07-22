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
/*!
 * \file ml_nox_preconditioner2.cpp
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
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#ifdef ML_MPI // FIXME: Do we need this here?
#include <mpi.h>
#endif

// this class
#include "ml_nox_preconditioner.H"

/*----------------------------------------------------------------------*
 |  compute the nonlin multilevel preconditioner (private)   m.gee 01/05|
 |  with Jacobian supplied                                              |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::ML_Nox_compute_Jacobian_Nonlinearpreconditioner(
                                                const Epetra_Vector& xfine)
{
   int i;

   // extract the prolongators from the hierarchy 
   Epetra_CrsMatrix** P  = new Epetra_CrsMatrix*[ml_nlevel_];   
   for (i=0; i<ml_nlevel_; i++)
      P[i] = 0;
   for (i=1; i<ml_nlevel_; i++) // there is no Pmat on level 0
   {
      double t1 = GetClock();
      int    maxnnz  = 0;
      double cputime;
      ML_Operator2EpetraCrsMatrix(&(ml_->Pmat[i]), P[i], maxnnz, false, cputime);
      P[i]->OptimizeStorage();
      double t2 = GetClock() - t1;
      if (ml_printlevel_>5 && 0 == comm_.MyPID())
            cout << "ML (level " << i << "): extraction of P in " << t2 << " sec\n";
   }

   // create the vector of nonlinear levels if it does not exist yet
   n_nlnlevel_ = ml_nlevel_;
   if (!nlnLevel_)
   {
      nlnLevel_   = new ML_NOX::ML_Nox_NonlinearLevel*[n_nlnlevel_];
      for (i=0; i<n_nlnlevel_; i++)
         nlnLevel_[i] = 0;
   }
   else // the vector of nonlinear levels already exists
   {
      nlnLevel_[0]->destroyP();
      for (i=0; i<n_nlnlevel_; i++)
      {
         nlnLevel_[i]->setP(NULL);
         delete nlnLevel_[i];
         nlnLevel_[i] = 0;
      }
   }
   
   // loop all levels and create the nonlinear level class
   for (i=0; i<n_nlnlevel_; i++)
   {
      // choose a blocks size
      int nullspdim = ml_dim_nullsp_;
      int numPDE    = 0;
      if (i==0) 
         numPDE = ml_numPDE_;
      else if (ml_spatialDimension_==1)
         numPDE = 1;
      else if (ml_spatialDimension_==2)
         numPDE = 3;
      else if (ml_spatialDimension_==3)
         numPDE = 6;
      else
      {
          cout << "**ERR**: ML_NOX::ML_Nox_Preconditioner::ML_Nox_compute_Jacobian_Nonlinearpreconditioner:\n"
               << "**ERR**: spatialDimension=" << ml_spatialDimension_ << " unknown\n"
               << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      
      // choose the nonlinear solver
      bool isnlnCG  = true;
      int  nitersCG = 100;
      if (i==0)
      {
         isnlnCG  = usenlnCG_fine_;
         nitersCG = nitersCG_fine_;
      }
      else if (i == n_nlnlevel_-1)
      {
         isnlnCG  = usenlnCG_coarse_;
         nitersCG = nitersCG_coarse_;
      }
      else
      {
         isnlnCG  = usenlnCG_;
         nitersCG = nitersCG_;
      }
      
      // create new nonlinear level  
      if (!nlnLevel_[i])  
         nlnLevel_[i] = new ML_NOX::ML_Nox_NonlinearLevel(
                                     i,ml_nlevel_,ml_printlevel_,ml_,ag_,P,
                                     interface_,comm_,xfine,ismatrixfree_,
                                     matfreelev0_,isnlnCG,nitersCG,
                                     useBroyden_,fineJac_,ml_fsmoothertype_,
                                     ml_smoothertype_,ml_coarsesolve_,
                                     nsmooth_fine_,nsmooth_,nsmooth_coarse_,
                                     FAS_normF_,FAS_nupdate_,
                                     3000,numPDE,nullspdim);
      else
      {
          cout << "**ERR**: ML_NOX::ML_Nox_Preconditioner::ML_Nox_compute_Jacobian_Nonlinearpreconditioner:\n"
               << "**ERR**: nlnLevel_[i] != 0 ????????\n"
               << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      } 
      
   }

   // in the nonlinear preconditioner, we don't need the ml hierarchy anymore
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
   
   return(true);
}

/*----------------------------------------------------------------------*
 |  compute the nonlin multilevel preconditioner (private)   m.gee 02/05|
 |  without Jacobian supplied                                           |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::ML_Nox_compute_Matrixfree_Nonlinearpreconditioner(
                                                const Epetra_Vector& xfine)
{
   int i;
   
   // extract the prolongators from the hierarchy 
   Epetra_CrsMatrix** P  = new Epetra_CrsMatrix*[ml_nlevel_];   
   for (i=0; i<ml_nlevel_; i++)
      P[i] = 0;
   for (i=1; i<ml_nlevel_; i++) // there is no Pmat on level 0
   {
      double t1 = GetClock();
      int    maxnnz  = 0;
      double cputime;
      ML_Operator2EpetraCrsMatrix(&(ml_->Pmat[i]), P[i], maxnnz, false, cputime);
      P[i]->OptimizeStorage();
      double t2 = GetClock() - t1;
      if (ml_printlevel_>5 && 0 == comm_.MyPID())
            cout << "ML (level " << i << "): extraction of P in " << t2 << " sec\n";
   }

   // construct the vector of coarse level matrixfree problems
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
      ml_matfreelevel_ = new ML_NOX::ML_Nox_MatrixfreeLevel*[nmatfreelevel_];
      for (i=0; i<nmatfreelevel_; i++)
         ml_matfreelevel_[i] = 0;
   }

   // create the vector of nonlinear levels if it does not exist yet
   n_nlnlevel_ = ml_nlevel_;
   if (!nlnLevel_)
   {
      nlnLevel_   = new ML_NOX::ML_Nox_NonlinearLevel*[n_nlnlevel_];
      for (i=0; i<n_nlnlevel_; i++)
         nlnLevel_[i] = 0;
   }
   else // the vector of nonlinear levels already exists
   {
      // NOTE: the matfree level and the nonlinear level share the
      //       coarseinterface, so they also share the vector of P's
      //       -> no need to destroyP again here
      for (i=0; i<n_nlnlevel_; i++)
      {
         nlnLevel_[i]->setP(NULL);
         delete nlnLevel_[i];
         nlnLevel_[i] = 0;
      }
      delete [] nlnLevel_;
      nlnLevel_   = new ML_NOX::ML_Nox_NonlinearLevel*[n_nlnlevel_];
      for (i=0; i<n_nlnlevel_; i++)
         nlnLevel_[i] = 0;
   }

   // the number of matrixfree levels has to be equal to the number of nonlinear levels
   if (nmatfreelevel_ !=  n_nlnlevel_)
   {
       cout << "**ERR**: ML_NOX::ML_Nox_Preconditioner::ML_Nox_compute_Matrixfree_Nonlinearpreconditioner:\n"
            << "**ERR**: nmatfreelevel_=" << nmatfreelevel_ << " != n_nlnlevel_=" << n_nlnlevel_ << endl
            << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }  
   
   // loop all levels, create the matrixfree and then the nonlinear level class
   for (i=0; i<n_nlnlevel_; i++)
   {
      if (comm_.MyPID()==0 && ml_printlevel_ > 5 )
         cout << "ML (level " << i << "): Entering FD-coarselevel (re)construction\n"; fflush(stdout);
      bool isJacobismoother=false;
      
      // for Jacobi smoothing, compute main diagonal only
      if (i==0 && ml_fsmoothertype_ == "Jacobi")
         isJacobismoother=true;
      else if (i<ml_coarsestlev_ && ml_smoothertype_ == "Jacobi")
         isJacobismoother=true;
      else
         isJacobismoother=false;
      if (i==ml_coarsestlev_ && ml_coarsesolve_ == "Jacobi")
         isJacobismoother=true;
      
      // when using Newton's method, compute complete Jacobian
      if (i==0 && usenlnCG_fine_==false)
         isJacobismoother=false;
      if (i==n_nlnlevel_-1 && usenlnCG_coarse_==false)
         isJacobismoother=false;
      if (i!=0 && i!=n_nlnlevel_-1 && usenlnCG_==false)
         isJacobismoother=false;
      
      int bsize;
      if (i==0) bsize = ml_numPDE_;
      else      bsize = ml_dim_nullsp_;
      
      if (!(ml_matfreelevel_[i])) // create a new matrixfree level   
         ml_matfreelevel_[i] = new ML_NOX::ML_Nox_MatrixfreeLevel(
                                           i,ml_nlevel_,ml_printlevel_,ml_,
                                           ag_,P,interface_,comm_,xfine,
                                           fd_alpha_,fd_beta_,fd_centered_,
                                           isJacobismoother,bsize);
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
      
      // check the matrix for zero rows and fix the main diagonal
      if (fixdiagonal_)
         fix_MainDiagonal(&tmpMat,i);
      tmpMat->OptimizeStorage();
      
      // get the coarse interface from the matfreelevel
      ML_NOX::Nox_CoarseProblem_Interface* coarseinterface = 
                                     ml_matfreelevel_[i]->getCoarseinterface();
      if (!coarseinterface)
      {
         cout << "**ERR**: ML_Nox_Preconditioner::ML_Nox_compute_Matrixfree_Linearpreconditioner:\n"
              << "**ERR**: ML_Epetra::ML_Nox_MatrixfreeLevel::getCoarseinterface() on level " << i << "\n"
              << "**ERR**: returned NULL-ptr\n"
              << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      
      // choose a blocks size
      int nullspdim = ml_dim_nullsp_;
      int numPDE    = 0;
      if (i==0) 
         numPDE = ml_numPDE_;
      else if (ml_spatialDimension_==1)
         numPDE = 1;
      else if (ml_spatialDimension_==2)
         numPDE = 3;
      else if (ml_spatialDimension_==3)
         numPDE = 6;
      else
      {
          cout << "**ERR**: ML_NOX::ML_Nox_Preconditioner::ML_Nox_compute_Jacobian_Nonlinearpreconditioner:\n"
               << "**ERR**: spatialDimension=" << ml_spatialDimension_ << " unknown\n"
               << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      }
      
      // choose the nonlinear solver
      bool isnlnCG  = true;
      int  nitersCG = 100;
      if (i==0)
      {
         isnlnCG  = usenlnCG_fine_;
         nitersCG = nitersCG_fine_;
      }
      else if (i == n_nlnlevel_-1)
      {
         isnlnCG  = usenlnCG_coarse_;
         nitersCG = nitersCG_coarse_;
      }
      else
      {
         isnlnCG  = usenlnCG_;
         nitersCG = nitersCG_;
      }
      
      // create new nonlinear level  
      if (!nlnLevel_[i]) // create new nonlinear level   
         nlnLevel_[i] = new ML_NOX::ML_Nox_NonlinearLevel(
                                     i,ml_nlevel_,ml_printlevel_,ml_,ag_,P,
                                     interface_,comm_,xfine,ismatrixfree_, 
                                     isnlnCG,nitersCG,useBroyden_,
                                     ml_fsmoothertype_,ml_smoothertype_,
                                     ml_coarsesolve_,
                                     nsmooth_fine_,nsmooth_,nsmooth_coarse_,
                                     FAS_normF_,
                                     FAS_nupdate_,3000,numPDE,nullspdim,tmpMat,
                                     coarseinterface);
      else
      {
          cout << "**ERR**: ML_NOX::ML_Nox_Preconditioner::ML_Nox_compute_Jacobian_Nonlinearpreconditioner:\n"
               << "**ERR**: nlnLevel_[i] != 0 ????????\n"
               << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
      } 
   } // loop over all levels

   // in the nonlinear preocnditioner, we don't need the ml hierarchy anymore
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

   return(true);
}

/*----------------------------------------------------------------------*
 |  apply preconditioner                         (private)   m.gee 02/05|
 |  nonlinear,                                                          |
 *----------------------------------------------------------------------*/
int ML_NOX::ML_Nox_Preconditioner::ML_Nox_ApplyInverse_NonLinear(
                                              const Epetra_MultiVector& X, 
                                              Epetra_MultiVector& Y) const
{
   if (!noxsolver_)
   {
      cout << "**ERR**: ML_NOX::ML_Nox_Preconditioner::ML_Nox_ApplyInverse_NonLinear:\n"
           << "**ERR**: noxsolver not registered, use set_nox_solver(solver)!\n"
           << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   
   // get the solution group from the nox solver
   const NOX::EpetraNew::Group& finalGroup = 
   dynamic_cast<const NOX::EpetraNew::Group&>(noxsolver_->getSolutionGroup());
   // get current guess and F
   const Epetra_Vector& currentSolution = 
   (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();      
   const Epetra_Vector& currentF = 
   (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getF())).getEpetraVector();
   double norm2;
   currentF.Norm2(&norm2);
   
   // make a copy of currentSolution and currentF
   Epetra_Vector* f = new Epetra_Vector(View,X,0);
   Epetra_Vector* x = new Epetra_Vector(Copy,currentSolution,0);
   
   // call the cycle
   if (ml_printlevel_>0 && comm_.MyPID()==0)
      cout << "\n\nML :============Entering Nonlinear V-cycle============\n";
   bool converged = false;
   double t3 = GetClock();
   ML_Nox_FAS_cycle(f,x,0,&converged,&norm2); 
   double t4 = GetClock();
   if (ml_printlevel_>0 && comm_.MyPID()==0)
      cout << "ML :============V-cycle time is : " << (t4-t3) << " sec\n";
   if (converged && ml_printlevel_>0 && comm_.MyPID()==0)
      cout << "ML :============Nonlinear preconditioner converged====\n";
   
   // copy correction to Y
   Y.Scale(1.0,*x);
   
   // tidy up
   if (f) delete f;
   if (x) delete x;
   
   return(0);
}

/*----------------------------------------------------------------------*
 |                                               (public)    m.gee 01/05|
 | run a nonlinear FAS cycle solver with nlnCG on each level            |
 *----------------------------------------------------------------------*/
int ML_NOX::ML_Nox_Preconditioner::solve() 
{
   int i;
   bool status=false;

   // get the initial solution from the interface
   const Epetra_Vector* xfine = interface_.getSolution();

   // compute the preconditioner
   double t5 = GetClock();
   status = computePreconditioner(*xfine);
   double t6 = GetClock();
   if (ml_printlevel_>0 && comm_.MyPID()==0)
      cout << "ML :============setup time : " << (t6-t5) << " sec\n\n\n";
   
   
   
   // make sure everything is sane
   if (isinit() != true)
   {
       cout << "**ERR**: ML_NOX::ML_Nox_Preconditioner::solve:\n"
            << "**ERR**: initflag is false\n"
            << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   if (islinearPrec_ != false)
   {
       cout << "**ERR**: ML_NOX::ML_Nox_Preconditioner::solve:\n"
            << "**ERR**: 'this' has not been initialized as nonlinear multigrid\n"
            << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   // make sure that nlnLevel_[level] exists for all levels we are going to use
   if (!nlnLevel_)
   {
       cout << "**ERR**: ML_NOX::ML_Nox_Preconditioner::solve:\n"
            << "**ERR**: nlnLevel_ is NULL\n"
            << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
   }
   
   // get the residual matching this solution
   Epetra_Vector* f = new Epetra_Vector(xfine->Map(),false);
   Epetra_Vector* x = new Epetra_Vector(Copy,*xfine,0);
 
   nlnLevel_[0]->setModifiedSystem(false, NULL, NULL);
   nlnLevel_[0]->computeF(*x,*f,NOX::EpetraNew::Interface::Required::Residual);

   double t1 = GetClock();
   {   
   for (i=1; i<=FAS_maxcycle_; i++)
   {
      if (ml_printlevel_>0 && comm_.MyPID()==0)
         cout << "\n\n\nML :============Entering FAS-V-cycle # " << i << "============\n";
      bool converged = false;
      double t3 = GetClock();
      status = ML_Nox_FAS_cycle1(f,x,0,&converged); 
      double t4 = GetClock();
      if (ml_printlevel_>0 && comm_.MyPID()==0)
         cout << "ML :============V-cycle time is : " << (t4-t3) << " sec\n";
      if (converged)
      { 
         double t2 = GetClock();
         if (ml_printlevel_>0 && comm_.MyPID()==0)
            cout << "ML :============FAS converged after " << i << " V-cycles============\n"
                 << "ML :============solution time : " << (t2-t1) << " sec\n\n\n";
         if (f) delete f;
         if (x) delete x;
         return(0);
      }
      computePreconditioner(*x);
   }
   if (f) delete f;
   if (x) delete x;
   return(-1);
   }
}


/*----------------------------------------------------------------------*
 |                                               (private)   m.gee 02/05|
 | the FAS cycle as a preconditioner                                    |
 |                                                                      |
 | f : on input the residual                                            |
 | x : on input the solution vector                                     |    
 | both vectors match the level level                                   |
 |                                                                      |
 | f : on output nothing at the moment                                  |
 | x : on output a correction to the solution                           |
 |     on output the correction vector (level==0)                       |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::ML_Nox_FAS_cycle(Epetra_Vector* f, Epetra_Vector* x,
                                                     int level, bool* converged, double* finenorm) const
{
   // we want to converge this level some digits better than the finer one
   double thisnorm = *finenorm/10000;
   if (thisnorm<FAS_normF_) thisnorm = FAS_normF_;
   
   Epetra_Vector* fbar    = 0;
   Epetra_Vector* xbar    = 0;
   Epetra_Vector* fxbar   = 0;

   //======reached coarsest level========================================
   if (level==ml_coarsestlev_)
   {
      // create the FAS-problem
      fbar  = new Epetra_Vector(Copy,*f,0);
      xbar  = new Epetra_Vector(Copy,*x,0);
      fxbar = new Epetra_Vector(xbar->Map(),false);
      nlnLevel_[level]->setModifiedSystem(false,NULL,NULL);
      nlnLevel_[level]->computeF(*xbar,*fxbar,NOX::EpetraNew::Interface::Required::Residual);
      nlnLevel_[level]->setModifiedSystem(true,fbar,fxbar);
      // iterate on the FAS-problem
      if (FAS_coarsesmooth_>0)
        *converged = nlnLevel_[level]->iterate(f,x,FAS_coarsesmooth_,&thisnorm);
      // calculate the correction
      x->Update(-1.0,*xbar,1.0);
      // reset the system
      nlnLevel_[level]->setModifiedSystem(false,NULL,NULL);
      // tidy up
      delete fbar;  fbar  = 0;
      delete xbar;  xbar  = 0;
      delete fxbar; fxbar = 0;
      return true;
   }

   //======set up FAS on this level=================
   // copy input residual f for use in FAS-postsmoothing
   fbar  = new Epetra_Vector(Copy,*f,0);
   xbar  = new Epetra_Vector(Copy,*x,0);
   fxbar = new Epetra_Vector(xbar->Map(),false);
   nlnLevel_[level]->setModifiedSystem(false,NULL,NULL);
   nlnLevel_[level]->computeF(*xbar,*fxbar,NOX::EpetraNew::Interface::Required::Residual);
   nlnLevel_[level]->setModifiedSystem(true,fbar,fxbar);

   //======presmoothing on the FAS system===========================
   double prenorm = thisnorm;
   if (level > 0 && FAS_presmooth_>0)
      *converged = nlnLevel_[level]->iterate(f,x,FAS_presmooth_,&prenorm);
   else if (level==0 && FAS_prefinesmooth_>0)
      *converged = nlnLevel_[level]->iterate(f,x,FAS_prefinesmooth_,&prenorm);
   // a converged level is enough to be happy, don't go any coarser
   if (*converged) 
   {
      nlnLevel_[level]->setModifiedSystem(false,NULL,NULL);
      x->Update(-1.0,*xbar,1.0);
      if (fbar)  delete fbar;  fbar  = 0;
      if (xbar)  delete xbar;  xbar  = 0;
      if (fxbar) delete fxbar; fxbar = 0;
      return true;
   }
   
   //======restrict to next coarse level=================================
   Epetra_Vector* xcoarse = nlnLevel_[level]->restrict_to_next_coarser_level(x,level,level+1);
   Epetra_Vector* fcoarse = nlnLevel_[level]->restrict_to_next_coarser_level(f,level,level+1);
   
   //======call this cycle on next coarser level=========================
   bool coarseconverged=false;
   ML_Nox_FAS_cycle(fcoarse,xcoarse,level+1,&coarseconverged,&prenorm);
   delete fcoarse; fcoarse = 0;
   
   //======prolong correction to this level==============================
   Epetra_Vector *xcorrect = nlnLevel_[level]->prolong_to_this_level(xcoarse,level,level+1);
   delete xcoarse; xcoarse = 0;
   
   //===== apply constraints (this may or may not be in here, it converges better without,
   // but is more stable with)
   // nlnLevel_[level]->ApplyAllConstraints(*xcorrect);

   //======update this level's solution==================================
   x->Update(1.0,*xcorrect,1.0);
   delete xcorrect; xcorrect = 0;
   
   //======do this level's FAS-iteration or postsmoothing================
   double postnorm = thisnorm;
   // iterate
   if (level>0 && FAS_postsmooth_>0)
      *converged = nlnLevel_[level]->iterate(f,x,FAS_postsmooth_,&postnorm);
   else if (level==0 && FAS_postfinesmooth_>0)
      *converged = nlnLevel_[level]->iterate(f,x,FAS_postfinesmooth_,&postnorm);
   // reset the system
   nlnLevel_[level]->setModifiedSystem(false,NULL,NULL);
   // for FAS calulate the correction and tidy up
   x->Update(-1.0,*xbar,1.0);
   if (fbar)  delete fbar;  fbar  = 0;
   if (xbar)  delete xbar;  xbar  = 0;
   if (fxbar) delete fxbar; fxbar = 0;

   return true;      
}

/*----------------------------------------------------------------------*
 |                                               (private)   m.gee 01/05|
 | the FAS cycle as a solver                                            |
 |                                                                      |
 | f : on input the residual                                            |
 | x : on input the solution vector                                     |    
 | both vectors match the level level                                   |
 |                                                                      |
 | f : on output nothing at the moment                                  |
 | x : on output a correction to the solution (level>0)                 |
 |     on output the solution vector (level==0)                         |
 *----------------------------------------------------------------------*/
bool ML_NOX::ML_Nox_Preconditioner::ML_Nox_FAS_cycle1(Epetra_Vector* f, Epetra_Vector* x,
                                                     int level, bool* converged)
{
   Epetra_Vector* fbar    = 0;
   Epetra_Vector* xbar    = 0;
   Epetra_Vector* fxbar   = 0;
   
   //======reached coarsest level========================================
   if (level==ml_coarsestlev_)
   {
      // create the FAS-problem
      fbar  = new Epetra_Vector(Copy,*f,0);
      xbar  = new Epetra_Vector(Copy,*x,0);
      fxbar = new Epetra_Vector(xbar->Map(),false);
      nlnLevel_[level]->setModifiedSystem(false,NULL,NULL);
      nlnLevel_[level]->computeF(*xbar,*fxbar,NOX::EpetraNew::Interface::Required::Residual);
      nlnLevel_[level]->setModifiedSystem(true,fbar,fxbar);
      // iterate on the FAS-problem
      if (FAS_coarsesmooth_>0)
        *converged = nlnLevel_[level]->iterate(f,x,FAS_coarsesmooth_);
      // calculate the correction
      x->Update(-1.0,*xbar,1.0);
      // reset the system
      nlnLevel_[level]->setModifiedSystem(false,NULL,NULL);
      // tidy up
      delete fbar;  fbar  = 0;
      delete xbar;  xbar  = 0;
      delete fxbar; fxbar = 0;
      return true;
   }
   
   //======if not finest level, set up FAS on this level=================
   if (level>0)
   {
      // copy input residual f for use in FAS-postsmoothing
      fbar  = new Epetra_Vector(Copy,*f,0);
      xbar  = new Epetra_Vector(Copy,*x,0);
      fxbar = new Epetra_Vector(xbar->Map(),false);
      nlnLevel_[level]->setModifiedSystem(false,NULL,NULL);
      nlnLevel_[level]->computeF(*xbar,*fxbar,NOX::EpetraNew::Interface::Required::Residual);
      nlnLevel_[level]->setModifiedSystem(true,fbar,fxbar);
   }

   //======presmoothing on the original system===========================
   if (level > 0 && FAS_presmooth_>0)
   {
      *converged = nlnLevel_[level]->iterate(f,x,FAS_presmooth_);
      if (*converged)
      {
         // for FAS calulate the correction and tidy up
         x->Update(-1.0,*xbar,1.0);
         nlnLevel_[level]->setModifiedSystem(false,NULL,NULL);
         delete fbar;  fbar  = 0;
         delete xbar;  xbar  = 0;
         delete fxbar; fxbar = 0;
         return true;
      }
   }
   else if (level==0 && FAS_prefinesmooth_>0)
   {
      *converged = nlnLevel_[level]->iterate(f,x,FAS_prefinesmooth_);
      // a converged fine level is enough to be happy....
      if (*converged)
         return true;
   }
   
   //======restrict to next coarse level=================================
   Epetra_Vector* xcoarse = nlnLevel_[level]->restrict_to_next_coarser_level(x,level,level+1);
   Epetra_Vector* fcoarse = nlnLevel_[level]->restrict_to_next_coarser_level(f,level,level+1);
   
   //======call this cycle on next coarser level=========================
   bool coarseconverged=false;
   ML_Nox_FAS_cycle1(fcoarse,xcoarse,level+1,&coarseconverged);
   delete fcoarse; fcoarse = 0;
   
   //======prolong correction to this level==============================
   Epetra_Vector *xcorrect = nlnLevel_[level]->prolong_to_this_level(xcoarse,level,level+1);
   delete xcoarse; xcoarse = 0;
   
   //======update this level's solution==================================
   x->Update(1.0,*xcorrect,1.0);
   delete xcorrect; xcorrect = 0;
   
   //======do this level's FAS-iteration or postsmoothing================
   if (level>0 && FAS_postsmooth_>0)
   {
      // iterate
      *converged = nlnLevel_[level]->iterate(f,x,FAS_postsmooth_);
      // reset the system
      nlnLevel_[level]->setModifiedSystem(false,NULL,NULL);
      // for FAS calulate the correction and tidy up
      x->Update(-1.0,*xbar,1.0);
      delete fbar;  fbar  = 0;
      delete xbar;  xbar  = 0;
      delete fxbar; fxbar = 0;
   }
   else if (FAS_postfinesmooth_>0) // level==0
   {
      // iterate
      *converged = nlnLevel_[level]->iterate(f,x,FAS_postfinesmooth_);
   }
   return true;      
}

#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) 
