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
 * \file nlnml_preconditioner1.cpp
 *
 * \class NLNML_NonlinearLevel
 *
 * \brief a nonlinear coarse grid class
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

// this class
#include "nlnml_nonlinearlevel.H"

using namespace std;

/*----------------------------------------------------------------------*
 |  ctor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_NonlinearLevel::NLNML_NonlinearLevel(
                int level, RefCountPtr<Teuchos::ParameterList> params,
                ML* ml, ML_Aggregate* ag,
                RefCountPtr< vector< RefCountPtr<Epetra_CrsMatrix> > > P,
                RefCountPtr<NLNML::NLNML_FineLevelNoxInterface> finterface,
                const Epetra_Comm& comm,
                RefCountPtr<Epetra_CrsMatrix> fineJac,
                const Epetra_Vector& xfine,
                bool isnlncg, int niterscg,
                int numpde, int dimns) :
level_(level),
isnlncg_(isnlncg),
ml_(ml),
ag_(ag),
thislevel_ml_(NULL),
thislevel_ag_(NULL),
params_(params),
P_(P),
fineinterface_(finterface),
comm_(comm),
fineJac_(fineJac)
{
  // get Jacobian of this level
  const Epetra_CrsGraph* graph = 0;
  if (!Level()) // fine grid
  {
    graph      = &(fineJac->Graph());
    SmootherA_ = fineJac;
  }
  else
  {
    // on coarse grids get Jacobian from hierarchy
    int    maxnnz  = 0;
    double cputime = 0.0;
    Epetra_CrsMatrix* tmp;
    ML_Operator2EpetraCrsMatrix(&(ml_->Amat[level_]),tmp,maxnnz,
                                false, cputime);
    SmootherA_ = rcp(tmp);
    SmootherA_->OptimizeStorage();
    graph = &(SmootherA_->Graph());
  }

  // generate the coarse level interface
  int outlevel = getParameter("nlnML output",0);
  int maxlevel = getParameter("nlnML max levels",2);
  coarseinterface_ =
    rcp(new NLNML::NLNML_CoarseLevelNoxInterface(fineinterface_,Level(),outlevel,
                                                 P_,graph->RowMap(),maxlevel));

  // on the fine grid, generate a pre/postoperator to apply constraints
  bool applyconstraints = false;
  if (!Level() && getParameter("nlnML apply constraints",false)==true)
    applyconstraints = true;

  if (!Level() && applyconstraints)
    prepost_ = rcp(new NLNML::NLNML_PrePostOperator(coarseinterface_,
                                                    fineinterface_));

  // get the current solution to this level
  xthis_ = rcp(coarseinterface_->restrict_fine_to_this(xfine));

  // create this level's preconditioner
  ML_Aggregate_Create(&thislevel_ag_);
  ML_Create(&thislevel_ml_,1);

  // set the Jacobian on level 0 of the local ml
  EpetraMatrix2MLMatrix(thislevel_ml_,0,SmootherA_.get());
  ML_Set_PrintLevel(OutLevel());
  ML_Aggregate_Set_CoarsenScheme_Uncoupled(thislevel_ag_);
  ML_Aggregate_Set_DampingFactor(thislevel_ag_, 0.0);
  ML_Aggregate_Set_Threshold(thislevel_ag_, 0.0);
  ML_Aggregate_Set_MaxCoarseSize(thislevel_ag_,1);
  ML_Aggregate_Set_NullSpace(thislevel_ag_,numpde,dimns,NULL,
                             SmootherA_->NumMyRows());
  int thislevel_nlevel = ML_Gen_MGHierarchy_UsingAggregation(thislevel_ml_,0,
                                                   ML_INCREASING,thislevel_ag_);
  if (thislevel_nlevel != 1)
  {
     cout << "**ERR**: NLNML::NLNML_NonlinearLevel::NLNML_NonlinearLevel:\n"
          << "**ERR**: ML generated a local hierarchy of " <<  thislevel_nlevel << " on level " << level_ << "\n"
          << "**ERR**: this is supposed to be 1 Level only!\n"
          << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }

  // set the smoother
  string smoothertype;
  int    nsmooth;
  if (!Level())
  {
    smoothertype = getParameter("nlnML linear smoother type fine level",(string)"SGS");
    nsmooth      = getParameter("nlnML linear smoother sweeps fine level",2);
  }
  else if (Level() != maxlevel-1)
  {
    smoothertype = getParameter("nlnML linear smoother type coarsest level",(string)"SGS");
    nsmooth      = getParameter("nlnML linear smoother sweeps coarsest level",2);
  }
  else
  {
    smoothertype = getParameter("nlnML linear smoother type medium level",(string)"SGS");
    nsmooth      = getParameter("nlnML linear smoother sweeps medium level",2);
  }
  Set_Smoother(ml_,ag_,Level(),maxlevel,thislevel_ml_,thislevel_ag_,smoothertype,nsmooth);

  // create this level's preconditioner class
  ML_Epetra::MultiLevelOperator* ml_tmp = new ML_Epetra::MultiLevelOperator(
                                               thislevel_ml_,comm_,
                                               SmootherA_->OperatorDomainMap(),
                                               SmootherA_->OperatorRangeMap());
  RefCountPtr<ML_Epetra::MultiLevelOperator> rcpml_tmp = rcp(ml_tmp);

  thislevel_prec_ = rcp(new NLNML::NLNML_ConstrainedMultiLevelOperator(
                                                             rcpml_tmp,
                                                             coarseinterface_,
                                                             applyconstraints));

  // ------------------------------------------------------------------------
  // set up NOX on this level
  // ------------------------------------------------------------------------
  nlparams_ = rcp(new Teuchos::ParameterList());
  Teuchos::ParameterList& printParams = nlparams_->sublist("Printing");
  printParams.set("MyPID", comm_.MyPID());
  printParams.set("Output Precision", 14);
  printParams.set("Output Processor", 0);
  if (OutLevel()>9)
    printParams.set("Output Information",
                             NOX::Utils::OuterIteration +
           NOX::Utils::Warning);
  else if (OutLevel()>8)
    printParams.set("Output Information",
                       NOX::Utils::Warning);
  else
    printParams.set("Output Information",0);

  if (!Level() && applyconstraints)
    nlparams_->sublist("Solver Options").set("User Defined Pre/Post Operator",prepost_);

  nlparams_->set("Nonlinear Solver", "Line Search Based");
  Teuchos::ParameterList& searchParams = nlparams_->sublist("Line Search");
  Teuchos::ParameterList* lsParamsptr  = 0;
  if (isnlncg)
  {
    searchParams.set("Method", "NonlinearCG");
    Teuchos::ParameterList& dirParams = nlparams_->sublist("Direction");
    dirParams.set("Method", "NonlinearCG");
    Teuchos::ParameterList& nlcgParams = dirParams.sublist("Nonlinear CG");
    nlcgParams.set("Restart Frequency", 10);
    nlcgParams.set("Precondition", "On");
    nlcgParams.set("Orthogonalize", "Polak-Ribiere");
    //nlcgParams.set("Orthogonalize", "Fletcher-Reeves");

    Teuchos::ParameterList& lsParams = nlcgParams.sublist("Linear Solver");
    lsParams.set("Aztec Solver", "CG");
    lsParams.set("Max Iterations", 1);
    lsParams.set("Tolerance", 1e-11);
    lsParams.set("Output Frequency", 0);
    lsParams.set("Preconditioning", "User Supplied Preconditioner");
    lsParams.set("Preconditioner","User Defined");
  }
  else
  {
    searchParams.set("Method", "Full Step");
    Teuchos::ParameterList& dirParams = nlparams_->sublist("Direction");
    dirParams.set("Method", "Newton");
    Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
    newtonParams.set("Forcing Term Method", "Constant");
    //newtonParams.set("Forcing Term Method", "Type 1");
    //newtonParams.set("Forcing Term Method", "Type 2");
    newtonParams.set("Forcing Term Minimum Tolerance", 1.0e-6);
    newtonParams.set("Forcing Term Maximum Tolerance", 0.1);

    Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
    lsParamsptr = &lsParams;
    lsParams.set("Size of Krylov Subspace", 100);
    lsParams.set("Aztec Solver", "GMRES");
    lsParams.set("Max Iterations", niterscg);
    double norm = getParameter("nlnML absolute residual tolerance",1.0e-06);
    lsParams.set("Tolerance",norm);
    if (OutLevel()>9)
      lsParams.set("Output Frequency",10);
    else
      lsParams.set("Output Frequency",0);
    lsParams.set("Preconditioning", "User Supplied Preconditioner");
    lsParams.set("Preconditioner","User Defined");
  }

  // create the initial guess
  NOX::Epetra::Vector initialguess(xthis_);

  // create interfaces
  thislevel_A_ = rcp(new NOX::Epetra::MatrixFree(printParams,coarseinterface_,*xthis_,false));
  if (isnlncg)
  {
    thislevel_linSys_ =
      rcp( new NLNML::NLNML_LinearSystem(thislevel_A_,thislevel_A_,null,
                                         coarseinterface_,thislevel_prec_,
                                         true,Level(),OutLevel()));
    group_ = rcp(new NOX::Epetra::Group(printParams,coarseinterface_,initialguess,
                                        thislevel_linSys_));
  }
  else // Newton's method
  {
    RefCountPtr<NOX::Epetra::Interface::Preconditioner> iPrec = rcp(this);
    iPrec.release();
    azlinSys_ =
      rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,*lsParamsptr,
                                               thislevel_A_,thislevel_A_,
                                               iPrec,thislevel_prec_,
                                               initialguess));
    group_ = rcp(new NOX::Epetra::Group(printParams,coarseinterface_,initialguess,
                                        azlinSys_));
  }

  // create some convergence test
  double norm = getParameter("nlnML absolute residual tolerance",1.0e-06);
  create_Nox_Convergencetest(norm,norm,1);

  // create the solver
  //JJH
  //solver_ = rcp(new NOX::Solver::Manager(group_,combo2_,nlparams_));
  solver_ = NOX::Solver::buildSolver(group_,combo2_,nlparams_);


  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             m.gee 3/06|
 *----------------------------------------------------------------------*/
NLNML::NLNML_NonlinearLevel::~NLNML_NonlinearLevel()
{
  if (thislevel_ag_) {
    ML_Aggregate_Destroy(&thislevel_ag_);
    thislevel_ag_ = NULL;
  }
  if (thislevel_ml_) {
    ML_Destroy(&thislevel_ml_);
    thislevel_ml_ = NULL;
  }
  return;
}


/*----------------------------------------------------------------------*
 |  compute this preconditioner (public, derived)             m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_NonlinearLevel::computePreconditioner(
                                     const Epetra_Vector& x,
             Epetra_Operator& M,
             Teuchos::ParameterList* precParams)
{
  return true;
}


/*----------------------------------------------------------------------*
 |  compute this Jacobian (public, derived)                   m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_NonlinearLevel::computeJacobian(const Epetra_Vector& x,
                                                  Epetra_Operator& Jac)
{
  return true;
}

/*----------------------------------------------------------------------*
 |  (private)                                                 m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_NonlinearLevel::Set_Smoother(ML* ml, ML_Aggregate* ag,
                                               int level, int nlevel,
                                               ML* thislevel_ml,
                                               ML_Aggregate* thislevel_ag,
                                               string smoothertype,
                                               int nsmooth)
{
  if (smoothertype=="SGS")
    ML_Gen_Smoother_SymGaussSeidel(thislevel_ml,0,ML_BOTH,nsmooth,0.67);
  else if (smoothertype == "Jacobi")
    ML_Gen_Smoother_Jacobi(thislevel_ml,0,ML_BOTH,nsmooth,0.2);
  else if (smoothertype == "AmesosKLU")
    ML_Gen_Smoother_Amesos(thislevel_ml,0,ML_AMESOS_KLU,-1,0.0,1);
  else if ((smoothertype == "MLS")||(smoothertype == "Cheby"))
    ML_Gen_Smoother_Cheby(thislevel_ml,0,ML_BOTH,30.,nsmooth);
  else if (smoothertype == "BSGS")
  {
    int  nblocks  = 0;
    int* blocks   = NULL;
    int* blockpde = NULL;

    // try to get nodal blocks from the VBMETIS aggregation scheme
    ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS(ag,level,nlevel,
                                                &nblocks,&blocks,&blockpde);

    if (nblocks && blocks);
    else
      ML_Gen_Blocks_Aggregates(ag,level,&nblocks,&blocks);

    ML_Gen_Smoother_VBlockSymGaussSeidel(thislevel_ml,0,ML_BOTH,nsmooth,1.,
                                         nblocks,blocks);
    if (nblocks && blocks)
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

    // try to get nodal blocks from the VBMETIS aggregation scheme
    ML_Aggregate_Get_Vblocks_CoarsenScheme_VBMETIS(ag,level,nlevel,
                                                &nblocks,&blocks,&blockpde);

    if (nblocks && blocks);
    else
      ML_Gen_Blocks_Aggregates(ag,level,&nblocks,&blocks);

    ML_Gen_Smoother_BlockDiagScaledCheby(thislevel_ml,0,ML_BOTH,30.,nsmooth,
                                         nblocks,blocks);
    if (nblocks && blocks)
    {
      ML_free(blocks);
      ML_free(blockpde);
    }
  }
  else
  {
    cout << "**ERR**: NLNML::NLNML_NonlinearLevel::Set_Smoother:\n"
         << "**ERR**: unknown type of smoother: " <<  smoothertype << " on level " << Level() << "\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }
  return true;
}



/*----------------------------------------------------------------------*
 |  (private)                                                 m.gee 3/06|
 *----------------------------------------------------------------------*/
void NLNML::NLNML_NonlinearLevel::create_Nox_Convergencetest(double normf,
                                                             double norm_update,
                                                             int maxiter)
{
  absresid_ = rcp(new NOX::StatusTest::NormF(normf,NOX::StatusTest::NormF::Unscaled));
  nupdate_  = rcp(new NOX::StatusTest::NormUpdate(norm_update,NOX::StatusTest::NormUpdate::Unscaled));
  combo1_   = rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
  combo1_->addStatusTest(absresid_);
  combo1_->addStatusTest(nupdate_);

  fv_ = rcp(new NOX::StatusTest::FiniteValue());
  maxiters_ = rcp(new NOX::StatusTest::MaxIters(maxiter));
  combo2_   = rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  combo2_->addStatusTest(fv_);
  combo2_->addStatusTest(maxiters_);
  combo2_->addStatusTest(combo1_);
  return;
}

/*----------------------------------------------------------------------*
 |  (private)                                                 m.gee 3/06|
 *----------------------------------------------------------------------*/
bool NLNML::NLNML_NonlinearLevel::Iterate(Epetra_Vector* f,
                                          Epetra_Vector* x,
                                          int numiter,
                                          double* norm)
{
  if (solver_==null || group_==null || combo2_==null)
  {
    cout << "**ERR**: NLNML::NLNML_NonlinearLevel::Iterate:\n"
         << "**ERR**: group_ || solver_ || combo2_ or even some of them are NULL\n"
         << "**ERR**: file/line: " << __FILE__ << "/" << __LINE__ << "\n"; throw -1;
  }

  // put x into group_
  NOX::Epetra::Vector X(*x,NOX::DeepCopy);
  group_->setX(X);

  // recreate the convergence test
  create_Nox_Convergencetest(*norm,*norm,numiter);

  // make a soft reset of the NOX solver
  //JJH
  //solver_->reset(group_,combo2_);
  solver_->reset(group_->getX(),combo2_);

  if (OutLevel() && !Comm().MyPID() && !coarseinterface_->isFAS())
  {
    printf("nlnML (level %d): Entering Nonlinear Smoother, Goal: %12.8e\n",level_,*norm);
    fflush(stdout);
  }
  else if (OutLevel() && !Comm().MyPID() && coarseinterface_->isFAS())
  {
    printf("nlnML (level %d): Entering FAS-Nonlinear Smoother, Goal: %12.8e\n",level_,*norm);
    fflush(stdout);
  }

  // iterate
  NOX::StatusTest::StatusType status = solver_->solve();

  // get # iterations
  int niter = solver_->getNumIterations();

  // get solution and F
  const NOX::Epetra::Group& finalGroup =
    dynamic_cast<const NOX::Epetra::Group&>(solver_->getSolutionGroup());
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
     if (OutLevel() > 0 && Comm().MyPID() == 0) {
        //cout << "ML (level " << level_ << "): NOX: "
        //     << niter << " iterations, Norm(F)="
        //     << norm2 << " , Converged\n"; fflush(stdout);
        printf("nlnML (level %d): NOX: %d iterations, Norm(F)=%12.8e , Converged\n",level_,niter,norm2);
        fflush(stdout);
     }
  }
  else if (status == NOX::StatusTest::Unconverged)
  {
     returnstatus = false;
     if (OutLevel() > 0 && Comm().MyPID() == 0) {
        //cout << "ML (level " << level_ << "): NOX: "
        //     << niter << " iterations, Norm(F)="
        //     << norm2 << ", Unconverged\n"; fflush(stdout);
        printf("nlnML (level %d): NOX: %d iterations, Norm(F)=%12.8e , Unonverged\n",level_,niter,norm2);
        fflush(stdout);
     }
  }
  else if (status == NOX::StatusTest::Failed)
  {
     returnstatus = false;
     if (OutLevel() > 0 && Comm().MyPID() == 0) {
        //cout << "ML (level " << level_ << "): NOX: "
        //     << niter << " iterations, Norm(F)=" << norm2 << ", Failed\n"; fflush(stdout);
        printf("nlnML (level %d): NOX: %d iterations, Norm(F)=%12.8e , Failed\n",level_,niter,norm2);
        fflush(stdout);
     }
  }
  else
  {
     returnstatus = false;
     if (Comm().MyPID() == 0) {
        //cout << "ML (level " << level_ << "): ***WRN*** NOX returned unknown status, Norm(F)="
        //     << norm2 << "\n"; fflush(stdout);
        printf("nlnML (level %d): ***WRN*** NOX: return status unknown, Norm(F)=%12.8e , Failed\n",level_,norm2);
        fflush(stdout);
     }
  }

  // reset number of calls to coarseinterface->computeF
  coarseinterface_->resetnumcallscomputeF();

  // update the solution
  x->Update(1.0,finalSolution,0.0);
  f->Update(1.0,finalF,0.0);

  return returnstatus;
}
#endif
