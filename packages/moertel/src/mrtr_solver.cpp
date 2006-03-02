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
matrix_(null),
matrixisnew_(true),
x_(null),
b_(null),
linearproblem_(null),
amesossolver_(null),
mlprec_(null),
moertelprec_(null),
aztecsolver_(null),
origmatrix_(null),
WT_(null),
B_(null),
I_(null)
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
                            RefCountPtr<Epetra_Vector> b,
                            MOERTEL::Manager& manager)
{
  SetParameters(params.get());
  SetSystem(matrix,x,b);
  origmatrix_ = manager.inputmatrix_;
  WT_         = manager.WT_;     
  B_          = manager.B_;      
  I_          = manager.I_;      
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
    ParameterList& amesosparams = params_->sublist("Amesos");
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
    ParameterList& mlparams    = params_->sublist("ML");
    ParameterList& aztecparams = params_->sublist("Aztec");
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

/*----------------------------------------------------------------------*
 |  solve a linear system (private)                          mwgee 01/06|
 |  using ML and AztecOO                                                |
 *----------------------------------------------------------------------*/
bool MOERTEL::Solver::Solve_MLAztec(ParameterList& mlparams, 
                                    ParameterList& aztecparams)
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
  // make initial guess x satisfy constraints: x = (I-WBT)x
  Epetra_Vector xtmp(B_->DomainMap(),false);
  B_->Multiply(true,*x_,xtmp);
  Epetra_Vector xtmp2(x_->Map(),false);
  WT_->Multiply(true,xtmp,xtmp2);
  x_->Update(-1.0,xtmp2,1.0);
  
  // make rhs satisfy constraints
  WT_->Multiply(false,*b_,xtmp);
  B_->Multiply(false,xtmp,xtmp2);
  b_->Update(-1.0,xtmp2,1.0);
  
  if (preconditioner=="AZ_user_precond")
  {
    if (mlprec_==null || matrixisnew_);
    {
      mlprec_ = rcp(new ML_Epetra::MultiLevelPreconditioner(*matrix_,mlparams),true);
      // create a constrained mortar-ml-preconditioner
      moertelprec_ = rcp(new MOERTEL::ConstrainedPreconditioner(mlprec_,I_,WT_,B_));
    }  
  }
  
#if 0
  
  // serial and on 1 level only
  // in parallel the gids of range and domain map of P are not ok as
  // ML creates a linear map
  const ML* ml = mlprec_->GetML();
  int nlevel = ml->ML_num_actual_levels; 
  Epetra_CrsMatrix* P;
  int maxnnz=0;
  double cputime;
  ML_Operator2EpetraCrsMatrix(&(ml->Pmat[1]),P,maxnnz,false,cputime);
  Epetra_CrsMatrix* P2;
  ML_Operator2EpetraCrsMatrix(&(ml->Pmat[2]),P2,maxnnz,false,cputime);
  //MOERTEL::Print_Matrix("zzz_Ifine",*I_,1);
  //MOERTEL::Print_Matrix("zzz_Psmooth",*P,1);
  //MOERTEL::Print_Matrix("zzz_PPsmooth",*P2,1);
  //MOERTEL::Print_Matrix("zzz_B",*B_,1);
  //MOERTEL::Print_Matrix("zzz_WT",*WT_,1);
  //MOERTEL::Print_Matrix("zzz_Atilde",*matrix_,1);
  //MOERTEL::Print_Matrix("zzz_A",*origmatrix_,1);
  //MOERTEL::Print_Vector("zzz_x",*x_,1);
  //MOERTEL::Print_Vector("zzz_b",*b_,1);
  //exit(0);
  // test BT x = 0 <-- this is correct
  //Epetra_Vector* BTx = new Epetra_Vector(B_->DomainMap(),true);
  //B_->Multiply(true,*x_,*BTx);
  //cout << *BTx;
  
  //------------------------------------------------------------------- form BWT
  Epetra_CrsMatrix* BWT = MOERTEL::MatMatMult(*B_,false,*WT_,false,OutLevel());
  //cout << *BWT;
  
  //------------------------------------------------------- form fine grid ImBWT
  Epetra_CrsMatrix* ImBWT = new Epetra_CrsMatrix(Copy,I_->RowMap(),10,false);
  MOERTEL::MatrixMatrixAdd(*I_,false,1.0,*ImBWT,0.0);
  MOERTEL::MatrixMatrixAdd(*BWT,false,-1.0,*ImBWT,1.0);
  ImBWT->FillComplete();
  //cout << *ImBWT;
  
  //-------------------------------------------------------------- form coarse I
  /*
  Epetra_CrsMatrix* RI      = MOERTEL::MatMatMult(*P,true,*I_,false,OutLevel());  
  Epetra_CrsMatrix* tmp     = MOERTEL::MatMatMult(*RI,false,*P,false,OutLevel());  
  delete RI; RI = NULL;
  Epetra_CrsMatrix* Icoarse = MOERTEL::StripZeros(*tmp,1.0e-12);
  delete tmp; tmp = NULL;
  */
  Epetra_CrsMatrix* Icoarse = new Epetra_CrsMatrix(Copy,P->DomainMap(),10,false);
  for (int i=0; i<Icoarse->NumMyRows(); ++i)
  {
    double one = 1.0;
    Icoarse->InsertGlobalValues(i,1,&one,&i);
  }
  Icoarse->FillComplete();
  //cout << *Icoarse; // correct

#if 1 //-------------------------------------------------------- to get WTcoarse
  // padd WT to be of full size and square
  Epetra_CrsMatrix* WTsquare = new Epetra_CrsMatrix(Copy,I_->RowMap(),1,false);
  // padd in a zero in every row
  for (int i=0; i<WTsquare->NumMyRows(); ++i)
  { 
    double zero = 0.0;
    int grid = WTsquare->GRID(i);
    WTsquare->InsertGlobalValues(grid,1,&zero,&grid);
  }
  MOERTEL::MatrixMatrixAdd(*WT_,false,1.0,*WTsquare,0.0);
  WTsquare->FillComplete(I_->OperatorDomainMap(),I_->OperatorRangeMap());
  //cout << *WT_;
  //cout << *WTsquare;
  
  // form WTcoarse = R WTsquare P
  Epetra_CrsMatrix* RWTsquare = MOERTEL::MatMatMult(*P,true,*WTsquare,false,OutLevel());
  Epetra_CrsMatrix* WTcoarse  = MOERTEL::MatMatMult(*RWTsquare,false,*P,false,OutLevel());
  delete RWTsquare;
  //cout << *WTsquare;
  //cout << *WTcoarse;
#endif

#if 1 //--------------------------------------------------------- to get Bsquare
  // padd B to be full size and square
  Epetra_CrsMatrix* Bsquare = new Epetra_CrsMatrix(Copy,I_->RowMap(),1,false);
  // padd in a zero in ervery row
  for (int i=0; i<Bsquare->NumMyRows(); ++i)
  { 
    double zero = 0.0;
    int grid = Bsquare->GRID(i);
    Bsquare->InsertGlobalValues(grid,1,&zero,&grid);
  }
  MOERTEL::MatrixMatrixAdd(*B_,false,1.0,*Bsquare,0.0);
  Bsquare->FillComplete(I_->OperatorDomainMap(),I_->OperatorRangeMap());
#endif

  //------------------------------------------ there are 2 ways to form BWTcoarse
  
  //=== form BWTcoarse = R * BWT * P
  Epetra_CrsMatrix* tmp5      = MOERTEL::MatMatMult(*P,true,*BWT,false,OutLevel());  
  Epetra_CrsMatrix* BWTcoarse = MOERTEL::MatMatMult(*tmp5,false,*P,false,OutLevel());  
  delete tmp5; tmp5 = NULL;
  //cout << *BWTcoarse;
  
  //=== form BWTCoarse = R * B * P * R * WT * P
  
  
  
  
  //------------------------------------------------- form coarse grid ImBWTcoarse
  Epetra_CrsMatrix* ImBWTcoarse = new Epetra_CrsMatrix(Copy,Icoarse->RowMap(),10,false);
  MOERTEL::MatrixMatrixAdd(*Icoarse,false,1.0,*ImBWTcoarse,0.0);
  MOERTEL::MatrixMatrixAdd(*BWTcoarse,false,-1.0,*ImBWTcoarse,1.0);
  ImBWTcoarse->FillComplete();
  //cout << *ImBWTcoarse;

  //------------------- make r, note that b_ and x_ satisfy constraints, so does r
  Epetra_Vector* r = new Epetra_Vector(b_->Map(),true);
  matrix_->Multiply(false,*x_,*r);
  r->Update(1.0,*b_,-1.0);
  //cout << *r;
  //MOERTEL::Print_Vector("zzz_r",*r,0);

  //---------------------------------------- make sure WT * r is zero <-- correct!
  //Epetra_Vector* WTr = new Epetra_Vector(WT_->RangeMap(),true);
  //WT_->Multiply(false,*r,*WTr);
  //cout << *r;
  //cout << *WTr;

  //------------------------------------------------------------- allocate rcoarse
  Epetra_Vector* rcoarse = new Epetra_Vector(Icoarse->RowMap(),true);

  //------------------------------------------------------------------- build Pmod
  Epetra_CrsMatrix* tmp1 = MOERTEL::MatMatMult(*P,true,*ImBWT,false,OutLevel()); 
  Epetra_CrsMatrix* tmp2 = MOERTEL::MatMatMult(*ImBWTcoarse,false,*tmp1,false,OutLevel());
  delete tmp1; tmp1 = NULL;
  
  //Epetra_CrsMatrix* tmp3 = MOERTEL::MatMatMult(*P,true,*BWT,false,OutLevel());
  //Epetra_CrsMatrix* tmp4 = MOERTEL::MatMatMult(*BWTcoarse,false,*tmp3,false,OutLevel());
  //delete tmp3; tmp3 = NULL;

  Epetra_CrsMatrix* Pmod = new Epetra_CrsMatrix(Copy,P->OperatorRangeMap(),10,false);  
  MOERTEL::MatrixMatrixAdd(*tmp2,true,1.0,*Pmod,0.0);
  //MOERTEL::MatrixMatrixAdd(*tmp4,true,1.0,*Pmod,1.0);
  delete tmp2; tmp2 = NULL;
  //delete tmp4; tmp4 = NULL;
  Pmod->FillComplete(P->OperatorDomainMap(),P->OperatorRangeMap());    
  Pmod->Multiply(true,*r,*rcoarse);
  //cout << *rcoarse;  

  // fine grid matrix is matrix_, fine grid lhs, rhs are b_, x_
  //Constraints are satisfied on the fine grid by the modified system
  //   Atilde x = btilde
  //          r = btilde - Atilde x
  //   Atilde e = r  
  // to satisfy constraints, we have to have WT r = 0
  // then we also have 
  //   B^T e = 0 (which is checked true after the solve)
  
  // evaluate RP = I and PR != I <-- correct
  //Epetra_CrsMatrix* tmp6 = MOERTEL::MatMatMult(*P,true,*P,false,OutLevel());  
  //Epetra_CrsMatrix* RP   = MOERTEL::StripZeros(*tmp6,1.0e-12); delete tmp6;
  //cout << *RP;
  //Epetra_CrsMatrix* tmp7 = MOERTEL::MatMatMult(*P,false,*P,true,OutLevel());
  //Epetra_CrsMatrix* PR   = MOERTEL::StripZeros(*tmp7,1.0e-12); delete tmp7;
  //cout << *PR;
  
  
  // restrict rcoarse = R r
  //P->Multiply(true,*r,*rcoarse);
  //cout << *rcoarse;
  
  Epetra_Vector* WTcoarse_rcoarse = new Epetra_Vector(Icoarse->RowMap(),true);
  
  // see that WTcoarse * rcoarse != 0 ( with standard prolongator ) as expected
  //WTcoarse->Multiply(false,*rcoarse,*WTcoarse_rcoarse);
  //cout << *WTcoarse_rcoarse;
  
  // see that WTcoarse rcoarse = 0 ( with modified prolongator ) <-- not true
  Pmod->Multiply(true,*r,*rcoarse);
  //cout << *rcoarse;
  WTcoarse->Multiply(false,*rcoarse,*WTcoarse_rcoarse);
  cout << *WTcoarse_rcoarse;

  exit(0);

#endif


#if 1
  // create the Aztec solver
  aztecsolver_ = rcp(new AztecOO());  
  aztecsolver_->SetAztecDefaults();
  aztecsolver_->SetProblem(*linearproblem_);
  aztecsolver_->SetParameters(aztecparams,true);
  if (mlprec_ != null && moertelprec_ != null)
    aztecsolver_->SetPrecOperator(moertelprec_.get());
  
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

#endif  


#if 0
  // test whether the solution e (or x as initial guess was zero) satisfies
  // B^T e = 0 ok this is true
  Epetra_Vector* BTe = new Epetra_Vector(x_->Map(),true);
  B_->Multiply(true,*x_,*BTe);
  cout << *BTe;

#endif





  return true;
}





















