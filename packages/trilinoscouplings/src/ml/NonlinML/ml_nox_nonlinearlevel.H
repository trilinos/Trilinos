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
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef ML_NOX_NONLINEARLEVEL_H
#define ML_NOX_NONLINEARLEVEL_H

#include <vector>

// epetra-headers
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_MapColoring.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

// NOX-headers
#include <NOX.H>
#include <NOX_Epetra.H> 
#include "NOX_EpetraNew_Interface_Required.H"
#include "NOX_EpetraNew_Interface_Jacobian.H"
#include "NOX_EpetraNew_Interface_Preconditioner.H"
#include "NOX_EpetraNew_FiniteDifference.H"
#include "NOX_EpetraNew_FiniteDifferenceColoring.H"
#include "NOX_EpetraNew_LinearSystemAztecOO.H"
#include "NOX_EpetraNew_BroydenOperator.H"

// ML-headers
#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelOperator.h"

// the NOX::EpetraNew::LinearSystem - derived linear System  
#include "ml_nox_linearsystem.H"

// interface to user code on fine level
#include "ml_nox_fineinterface.H"

// interface to user code on coarse levels
#include "ml_nox_coarseinterface.H"
#include "ml_nox_coarse_prepostoperator.H"

#include "ml_nox_ConstrainedMultiLevelOperator.H"

// the graph coloring class - EpetraExt headers
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"

// ---------- Forward Declarations ----------
class Epetra_MultiVector;
class Epetra_BlockMap;
class Epetra_Comm;
class Nox_Problem_Interface;
class Nox_CoarseProblem_Interface;
class Epetra_MapColoring;
class Epetra_IntVector;

namespace ML_NOX
{
class  ML_Nox_NonlinearLevel : public NOX::EpetraNew::Interface::Preconditioner,
                               public NOX::EpetraNew::Interface::Jacobian
{

public:

  // constructor for case ismatrixfree==false
  explicit ML_Nox_NonlinearLevel(
       int level, int nlevel, int printlevel, ML* ml, ML_Aggregate* ag,
       Epetra_CrsMatrix** P, ML_NOX::Ml_Nox_Fineinterface& interface,
       const Epetra_Comm& comm,  const Epetra_Vector& xfine, bool ismatrixfree,
       bool matfreelev0, bool isnlnCG, int nitersCG, bool broyden, 
       Epetra_CrsMatrix* Jac, string fsmoothertype, string smoothertype, 
       string coarsesolvetype, 
       int nsmooth_fine, int nsmooth, int nsmooth_coarse, 
       double conv_normF, 
       double conv_nupdate, int conv_maxiter, int numPDE, int nullspdim);

  // constructor for case ismatrixfree==true
  explicit ML_Nox_NonlinearLevel(
       int level, int nlevel, int printlevel, ML* ml, ML_Aggregate* ag,
       Epetra_CrsMatrix** P, ML_NOX::Ml_Nox_Fineinterface& interface,
       const Epetra_Comm& comm,  const Epetra_Vector& xfine, bool ismatrixfree, 
       bool isnlnCG, int nitersCG, bool broyden, 
       string fsmoothertype, string smoothertype, string coarsesolvetype, 
       int nsmooth_fine, int nsmooth, int nsmooth_coarse,  
       double conv_normF, double conv_nupdate, int conv_maxiter,
       int numPDE, int nullspdim, Epetra_CrsMatrix* Mat, 
       ML_NOX::Nox_CoarseProblem_Interface* coarseinterface);

  // destructor
  ~ML_Nox_NonlinearLevel();
  
  //! Computes a user defined preconditioner. (derived)
  bool computePreconditioner(const Epetra_Vector& x,
			     NOX::Parameter::List* precParams = 0);

  // destroy the P-hierarchy in the coarse interface
  bool destroyP() { if (coarseinterface_) return coarseinterface_->destroyP(); 
                    else                  return false; }
  // set the ptr to the P-hierarchy in the coarseinterface_
  bool setP(Epetra_CrsMatrix** P) { if (coarseinterface_) return coarseinterface_->setP(P); 
                                    else return false; }
                                    
  // set the FAS-modified system flag
  bool setModifiedSystem(bool flag, Epetra_Vector* fbar, Epetra_Vector* fxbar)
  { return coarseinterface_->setModifiedSystem(flag,fbar,fxbar);}
  
  // iterate on this level's nln problem numiter times or until NOX converged
  bool iterate(Epetra_Vector* f, Epetra_Vector* x, int numiter); 
  bool iterate(Epetra_Vector* f, Epetra_Vector* x, int numiter, double* norm); 
  
  // restrict vector to next coarser level
  Epetra_Vector* restrict_to_next_coarser_level(Epetra_Vector* finevec, int current, int next) 
  { return coarseinterface_->restrict_to_next_coarser_level(finevec,current,next);}                                   
  
  // prolong vector from next coarser level to current level
  Epetra_Vector* prolong_to_this_level(Epetra_Vector* coarsevec, int current, int next) 
  { return coarseinterface_->prolong_to_this_level(coarsevec,current,next);}                                   
  

  // evaluate the nonlinear function on this level
  bool computeF(const Epetra_Vector& x, Epetra_Vector& F, 
                NOX::EpetraNew::Interface::Required::FillType fillFlag)
  { return coarseinterface_->computeF(x,F,fillFlag);}
  
  // Compute Jacobian (derived from NOX::EpetraNew::Interface::Jacobian)
  bool computeJacobian(const Epetra_Vector& x);
  
  // apply constrints using the coarse interface
  void ApplyAllConstraints(Epetra_Vector& gradient) 
  { coarseinterface_->ApplyAllConstraints(gradient); return; }

private:
  // don't want copy-ctor and = operator
  ML_Nox_NonlinearLevel(ML_NOX::ML_Nox_NonlinearLevel& old);
  ML_Nox_NonlinearLevel operator = (const ML_NOX::ML_Nox_NonlinearLevel& old);

  // create a NOX convergence test
  bool create_Nox_Convergencetest(double conv_normF, double conv_nupdate, int conv_maxiter);

  // set the smoother for one nonlinear level
  bool Set_Smoother(ML* ml, ML_Aggregate* ag, int level, int nlevel,
                    ML* thislevel_ml, ML_Aggregate* thislevel_ag, 
                    string smoothertype, int nsmooth);

  // iterate on this level's nln problem numiter times or until NOX converged
  bool applySmoother(Epetra_Vector* f, Epetra_Vector* x, int numiter); 
  
  int                                          level_;            // the level we are on here
  int                                          nlevel_;           // total number of levels
  int                                          ml_printlevel_;    // ML's output level
  ML*                                          ml_;               // the overall ML-hierarchy
  ML_Aggregate*                                ag_;               // the overall ML-hierarchy
  Epetra_CrsMatrix*                            SmootherA_;        // this level's preconditioner to this level's nlnCG
  ML_NOX::Ml_Nox_Fineinterface&                fineinterface_;    // the fineinterface
  const Epetra_Comm&                           comm_;             // Epetra_Comm
  ML_NOX::Nox_CoarseProblem_Interface*         coarseinterface_;  // this level's coarse interface
  ML_NOX::Ml_Nox_CoarsePrePostOperator*        coarseprepost_;
  ML*                                          thislevel_ml_;     // this level's 1-Level ML as a smoother to this level's nlnCG
  ML_Aggregate*                                thislevel_ag_;     // this level's 1-Level ML as a smoother to this level's nlnCG
  ML_NOX::ML_Nox_ConstrainedMultiLevelOperator* thislevel_prec_;   // this level's 1-Level ML-smoother
  Epetra_Vector*                               xthis_;            // this level's current solution

  NOX::EpetraNew::MatrixFree*                  thislevel_A_;      // this level's Matrixfree operator
  ML_NOX::Ml_Nox_LinearSystem*                 thislevel_linSys_; // this level's NOX-linearSystem
  NOX::Parameter::List*                        nlParams_;         // this level's NOX parameter list

  NOX::StatusTest::NormF*                      absresid_;         // a combination of NOX-convergence criteria
  NOX::StatusTest::NormUpdate*                 nupdate_;
  NOX::StatusTest::FiniteValue*                fv_;
  NOX::StatusTest::MaxIters*                   maxiters_;
  NOX::StatusTest::Combo*                      combo1_;           
  NOX::StatusTest::Combo*                      combo2_;           

  NOX::Epetra::Vector*                         initialGuess_;     // the initial guess for NOX
  NOX::EpetraNew::Group*                       group_;            // whatever that is, one needs it
  NOX::Solver::Manager*                        solver_;           // this level's NOX-solver (nlnCG)
  double                                       conv_normF_;       // this level's NOX's convergence test
  double                                       conv_nupdate_;
  int                                          conv_maxiter_;
  
  // --- for the matrixfree method ---
  bool                                         ismatrixfree_;    // flag to indicate whether this level is matrixfree or not
  
  // -- for the Newton method
  bool                                         isnlnCG_;          // flag to switch nlnCG/Neweton's method
  int                                          nitersCG_;         // number of Aztec-CG iterations as linear solve in mod. Newton's method
  bool                                         broyden_;          // use Broyden update for jacobian in case of Newton's method
  NOX::EpetraNew::LinearSystemAztecOO*         azlinSys_;         // the aztec linear system to be used with mod. Newton's method  
  Epetra_Vector*                               clone_;
  NOX::EpetraNew::BroydenOperator*             Broyd_;            // Broyden operator
}; // class  ML_Nox_NonlinearLevel 
}  // namespace ML_NOX


#endif // ML_NOX_NONLINEARLEVEL_H
#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) 
