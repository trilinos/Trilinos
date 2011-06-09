/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
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
# Questions? Contact Glen Hansen (Glen.Hansen@inl.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#include "mrtr_solver.H"
#include "mrtr_manager.H"
#include "mrtr_utils.H"
#include "Epetra_Time.h"
#include <EpetraExt_Transpose_RowMatrix.h>

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 12/05|
 *----------------------------------------------------------------------*/
MOERTEL::Solver::Solver(Epetra_Comm& comm, int outlevel) :
outlevel_(outlevel),
comm_(comm),
matrix_(Teuchos::null),
matrixisnew_(true),
x_(Teuchos::null),
b_(Teuchos::null),
linearproblem_(Teuchos::null),
amesossolver_(Teuchos::null),
mlprec_(Teuchos::null),
aztecsolver_(Teuchos::null),
origmatrix_(Teuchos::null),
WT_(Teuchos::null),
B_(Teuchos::null),
I_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 12/05|
 *----------------------------------------------------------------------*/
MOERTEL::Solver::~Solver()
{
}

/*----------------------------------------------------------------------*
 |  set a linear system (public)                             mwgee 12/05|
 *----------------------------------------------------------------------*/
void MOERTEL::Solver::SetSystem(Teuchos::RCP<Epetra_CrsMatrix> matrix,
                                Teuchos::RCP<Epetra_Vector> x,
                                Teuchos::RCP<Epetra_Vector> b)
{
  matrix_ = matrix;
  x_      = x;
  b_      = b;
  return;
}

/*----------------------------------------------------------------------*
 |  solve a linear system (public)                           mwgee 12/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Solver::Solve(Teuchos::RCP<Teuchos::ParameterList> params,
                            Teuchos::RCP<Epetra_CrsMatrix> matrix,
                            Teuchos::RCP<Epetra_Vector> x,
                            Teuchos::RCP<Epetra_Vector> b,
                            MOERTEL::Manager& manager)
{
  SetParameters(params.get());
  SetSystem(matrix,x,b);
  origmatrix_ = manager.inputmatrix_;
  WT_         = manager.WT_;     
  B_          = manager.B_;      
  I_          = manager.I_;
  Annmap_     = manager.annmap_; 
  lm_to_dof_  = manager.lm_to_dof_;     
  return Solve();
}

/*----------------------------------------------------------------------*
 |  solve a linear system (public)                           mwgee 12/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Solver::Solve()
{
  bool ok = false;
  
  //---------------------------------------------------------------------------
  // check the linear system  
  if (x_==Teuchos::null || b_==Teuchos::null || matrix_==Teuchos::null)
  {
    cout << "***ERR*** MOERTEL::Solver::Solve:\n"
         << "***ERR*** matrix and/or rhs and/or solution vector are Teuchos::null\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  //---------------------------------------------------------------------------
  // check for parameters
  if (params_==NULL)
  {
    cout << "***ERR*** MOERTEL::Solver::Solve:\n"
         << "***ERR*** solver parameters are Teuchos::null\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }  

  //---------------------------------------------------------------------------
  // (re)create a linear problem
  if (linearproblem_==Teuchos::null)
    linearproblem_ = Teuchos::rcp(new Epetra_LinearProblem());

  linearproblem_->SetLHS(x_.get());
  linearproblem_->SetRHS(b_.get());
  if (matrixisnew_)
    linearproblem_->SetOperator(matrix_.get());
    
  linearproblem_->CheckInput();
  //---------------------------------------------------------------------------
  // get type of solver to be used
  string solver = params_->get("Solver","None");
  
  //---------------------------------------------------------------------------
  // time the solution process
  Epetra_Time time(Comm());
  time.ResetStartTime();

  //---------------------------------------------------------------------------
  // use Amesos
  if (solver=="Amesos" || solver=="amesos" || solver=="AMESOS")
  {
	Teuchos::ParameterList& amesosparams = params_->sublist("Amesos");
    ok = Solve_Amesos(amesosparams);
    if (!ok)
    {
      cout << "***WRN*** MOERTEL::Solver::Solve:\n"
           << "***WRN*** Solve_Amesos returned false\n"
           << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    }
  }

  //---------------------------------------------------------------------------
  // use ML/Aztec
  else if (solver=="ML/Aztec" || solver=="ml/aztec" || solver=="ML" ||
           solver=="ml"       || solver=="Ml"       || solver=="Aztec" ||
           solver=="AZTEC"    || solver=="aztec")
  {
    // see whether we have a spd system
    string system = params_->get("System","None");
    if (system!="SPDSystem"  && system!="spdsystem" && system!="spd_system" && 
        system!="SPD_System" && system!="SPDSYSTEM" && system!="SPD_SYSTEM")
    {
      cout << "***ERR*** MOERTEL::Solver::Solve:\n"
           << "***ERR*** To use ML?Aztec for solution, parameter \"System\" hast to be \"SPDSystem\" \n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
	Teuchos::ParameterList& mlparams    = params_->sublist("ML");
	Teuchos::ParameterList& aztecparams = params_->sublist("Aztec");
    ok = Solve_MLAztec(mlparams,aztecparams);
    if (!ok)
    {
      cout << "***WRN*** MOERTEL::Solver::Solve:\n"
           << "***WRN*** Solve_MLAztec returned false\n"
           << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    }
  }


  //---------------------------------------------------------------------------
  // unknown solver
  else
  {
    cout << "***ERR*** MOERTEL::Solver::Solve:\n"
         << "***ERR*** solver type is unknown\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    ok = false;
  }


  //---------------------------------------------------------------------------
  // time the solution process
  double t = time.ElapsedTime();
  if (OutLevel()>5 && Comm().MyPID()==0)
    cout << "MOERTEL (Proc 0): Solution of system of equations in " << t << " sec\n";

  return ok;
}

/*----------------------------------------------------------------------*
 |  solve a linear system (private)                          mwgee 12/05|
 |  using Amesos                                                        |
 *----------------------------------------------------------------------*/
bool MOERTEL::Solver::Solve_Amesos(Teuchos::ParameterList& amesosparams)
{
  int ok = 0;

  //---------------------------------------------------------------------------
  // which amesos solver
  string solver       = amesosparams.get("Solver","None");
  bool   usetranspose = amesosparams.get("UseTranspose",false);
  if (solver=="None")
  {
    cout << "***ERR*** MOERTEL::Solver::Solve_Amesos:\n"
         << "***ERR*** No Amesos solver chosen\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //---------------------------------------------------------------------------
  // new amesos solver
  if (matrixisnew_ || amesossolver_==Teuchos::null)  
  {
    amesossolver_ = Teuchos::null;
    Amesos Factory;
    if (!Factory.Query(solver))
    {
      cout << "***ERR*** MOERTEL::Solver::Solve_Amesos:\n"
           << "***ERR*** Amesos solver '" << solver << "' not supported\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    amesossolver_ = Teuchos::rcp(Factory.Create(solver,*linearproblem_));
    if (amesossolver_.get()==0)
    {
      cout << "***ERR*** MOERTEL::Solver::Solve_Amesos:\n"
           << "***ERR*** Could not create Amesos solver '" << solver << "'\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    amesossolver_->SetUseTranspose(usetranspose);
    ok += amesossolver_->SetParameters(amesosparams);
    ok += amesossolver_->SymbolicFactorization();
    ok += amesossolver_->NumericFactorization();
    matrixisnew_ = false;
    ok += amesossolver_->Solve();
  }
  // neither the matrix is new nor the solver
  else
  {
    ok = amesossolver_->Solve();
  }
  
  if (ok==0) return true;
  else
  {
    cout << "***ERR*** MOERTEL::Solver::Solve_Amesos:\n"
         << "***ERR*** Amesos returned " << ok << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  return false;
}

/*----------------------------------------------------------------------*
 |  solve a linear system (private)                          mwgee 01/06|
 |  using ML and AztecOO                                                |
 *----------------------------------------------------------------------*/
bool MOERTEL::Solver::Solve_MLAztec(Teuchos::ParameterList& mlparams, 
                                    Teuchos::ParameterList& aztecparams)
{
  
  // create ML preconditioner if aztec parameter indicates user preconditioner
  string preconditioner = aztecparams.get("AZ_precond","none");
  if (preconditioner=="none")
  {
    cout << "***ERR*** MOERTEL::Solver::Solve_MLAztec:\n"
         << "***ERR*** Aztec parameter \"AZ_precond\" is not set\n"
         << "***ERR*** set to \"AZ_user_precond\" to use ML or to some Aztec internal method\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  // make initial guess x satisfy constraints: x = (I - W * B^T) x
  Epetra_Vector xtmp(B_->DomainMap(),false);
  Epetra_Vector xtmp2(x_->Map(),false);
 
  B_->Multiply(true,*x_,xtmp);
  WT_->Multiply(true,xtmp,xtmp2);
  x_->Update(-1.0,xtmp2,1.0);
  
  // make rhs satisfy constraints: b = (I - B * W^T) b
  WT_->Multiply(false,*b_,xtmp);
  B_->Multiply(false,xtmp,xtmp2);
  b_->Update(-1.0,xtmp2,1.0);
  
  if (preconditioner=="AZ_user_precond")
    if (mlprec_==Teuchos::null || matrixisnew_);
    {
#if 1
      mlprec_ = Teuchos::rcp(new MOERTEL::Mortar_ML_Preconditioner(matrix_,
                                                          origmatrix_,
                                                          WT_,B_,
                                                          Annmap_,
                                                          mlparams));
#else // change mlprec_ in mrtr_solver.H as well to test black box ML
      mlprec_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*matrix_,mlparams,true));
#endif
    }
    
  // create the Aztec solver
  aztecsolver_ = Teuchos::rcp(new AztecOO());  
  aztecsolver_->SetAztecDefaults();
  aztecsolver_->SetProblem(*linearproblem_);
  aztecsolver_->SetParameters(aztecparams,true);
  if (mlprec_ != Teuchos::null)
    aztecsolver_->SetPrecOperator(mlprec_.get());
  
  // solve it
  double tol  = aztecparams.get("AZ_tol",1.0e-05);
  int maxiter = aztecparams.get("AZ_max_iter",1000);
  aztecsolver_->Iterate(maxiter,tol);
  matrixisnew_ = false;
  const double* azstatus = aztecsolver_->GetAztecStatus();
  if (azstatus[AZ_why] == AZ_normal)
    return true;
  else if (azstatus[AZ_why] == AZ_breakdown)
  {
    if (Comm().MyPID() == 0)
    cout << "MOERTEL: ***WRN*** MOERTEL::Solver::Solve_MLAztec:\n"
         << "MOERTEL: ***WRN*** Aztec returned status AZ_breakdown\n"
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  else if (azstatus[AZ_why] == AZ_loss)
  {
    if (Comm().MyPID() == 0)
    cout << "MOERTEL: ***WRN*** MOERTEL::Solver::Solve_MLAztec:\n"
         << "MOERTEL: ***WRN*** Aztec returned status AZ_loss\n"
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  else if (azstatus[AZ_why] == AZ_ill_cond)
  {
    if (Comm().MyPID() == 0)
    cout << "MOERTEL: ***WRN*** MOERTEL::Solver::Solve_MLAztec:\n"
         << "MOERTEL: ***WRN*** Aztec returned status AZ_ill_cond\n"
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  else if (azstatus[AZ_why] == AZ_maxits)
  {
    if (Comm().MyPID() == 0)
    cout << "MOERTEL: ***WRN*** MOERTEL::Solver::Solve_MLAztec:\n"
         << "MOERTEL: ***WRN*** Aztec returned status AZ_maxits\n"
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  else
  {
    if (Comm().MyPID() == 0)
    cout << "MOERTEL: ***WRN*** MOERTEL::Solver::Solve_MLAztec:\n"
         << "MOERTEL: ***WRN*** Aztec returned unknown status: " << azstatus[AZ_why] << "\n"
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

#if 0
  // test whether the solution e (or x as initial guess was zero) satisfies
  // B^T e = 0 ok this is true
  Epetra_Vector* BTe = new Epetra_Vector(x_->Map(),true);
  B_->Multiply(true,*x_,*BTe);
  cout << *BTe;

#endif

  return true;
}





















