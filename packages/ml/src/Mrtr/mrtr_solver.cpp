/*
#@HEADER
# ************************************************************************
#
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
#include "mrtr_solver.H"
#include "Epetra_Time.h"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 12/05|
 *----------------------------------------------------------------------*/
MOERTEL::Solver::Solver(Epetra_Comm& comm, int outlevel) :
outlevel_(outlevel),
comm_(comm),
matrix_(null),
matrixisnew_(true),
x_(null),
b_(null),
amesossolver_(null)
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
void MOERTEL::Solver::SetSystem(RefCountPtr<Epetra_CrsMatrix> matrix,
                                RefCountPtr<Epetra_Vector> x,
                                RefCountPtr<Epetra_Vector> b)
{
  matrix_ = matrix;
  x_      = x;
  b_      = b;
  return;
}

/*----------------------------------------------------------------------*
 |  solve a linear system (public)                           mwgee 12/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Solver::Solve(RefCountPtr<Teuchos::ParameterList> params,
                            RefCountPtr<Epetra_CrsMatrix> matrix,
                            RefCountPtr<Epetra_Vector> x,
                            RefCountPtr<Epetra_Vector> b)
{
  SetParameters(params.get());
  SetSystem(matrix,x,b);
  return Solve();
}

/*----------------------------------------------------------------------*
 |  solve a linear system (public)                           mwgee 12/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Solver::Solve()
{
  //---------------------------------------------------------------------------
  // check the linear system  
  if (x_==null || b_==null || matrix_==null)
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
  if (linearproblem_==null)
    linearproblem_ = rcp(new Epetra_LinearProblem());

  //b_.get()->PutScalar(0.0);
  //x_.get()->PutScalar(0.0);
  linearproblem_->SetLHS(x_.get());
  linearproblem_->SetRHS(b_.get());
    
  if (matrixisnew_)
    linearproblem_->SetOperator(matrix_.get());
    
  linearproblem_->CheckInput();
  
  //---------------------------------------------------------------------------
  // get type of solver to be used
  string solver = params_->get("Solver","None");
  
  //---------------------------------------------------------------------------
  // use Amesos
  if (solver=="Amesos" || solver=="amesos" || solver=="AMESOS")
  {
    ParameterList amesosparams = params_->sublist("Amesos");
    bool ok = Solve_Amesos(amesosparams);
    if (!ok)
    {
      cout << "***ERR*** MOERTEL::Solver::Solve:\n"
           << "***ERR*** Solve_Amesos returned false\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }

  //---------------------------------------------------------------------------
  // unknown solver
  else
  {
    cout << "***ERR*** MOERTEL::Solver::Solve:\n"
         << "***ERR*** solver type is unknown\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  solve a linear system (public)                           mwgee 12/05|
 |  using Amesos                                                        |
 *----------------------------------------------------------------------*/
bool MOERTEL::Solver::Solve_Amesos(ParameterList& amesosparams)
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
  if (matrixisnew_ || amesossolver_==null)  
  {
    amesossolver_ = null;
    Amesos Factory;
    if (!Factory.Query(solver))
    {
      cout << "***ERR*** MOERTEL::Solver::Solve_Amesos:\n"
           << "***ERR*** Amesos solver '" << solver << "' not supported\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    amesossolver_ = rcp(Factory.Create(solver,*linearproblem_));
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





















