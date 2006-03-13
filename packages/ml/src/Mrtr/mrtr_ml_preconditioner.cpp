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
#include "mrtr_ml_preconditioner.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            m.gee 03/06|
 *----------------------------------------------------------------------*/
MOERTEL::Mortar_ML_Preconditioner::Mortar_ML_Preconditioner(
                                    RefCountPtr<Epetra_CrsMatrix> Atilde,
                                    RefCountPtr<Epetra_CrsMatrix> A,
                                    RefCountPtr<Epetra_CrsMatrix> WT,
                                    RefCountPtr<Epetra_CrsMatrix> B,
                                    ParameterList& mlparams,
                                    bool constructnow) :
iscomputed_(false),
mlparams_(mlparams),
Atilde_(Atilde),
A_(A),
WT_(WT),
B_(B)                                    
{
  label_  = "MOERTEL::Mortar_ML_Preconditioner";

  mlapiImBWT_.resize(0);                       
  mlapiImWBT_.resize(0);                       
  mlapiRmod_.resize(0);                       
  mlapiPmod_.resize(0);                       
  mlapiAtilde_.resize(0);
  mlapiS_.resize(0);

  if (constructnow) Compute();
    
  return;
}

/*----------------------------------------------------------------------*
 |  apply multigrid linear preconditioner (public)           m.gee 03/06|
 *----------------------------------------------------------------------*/
int MOERTEL::Mortar_ML_Preconditioner::ApplyInverse(
                     const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
{
  if (!iscomputed_) ML_THROW("Method Compute() must be called before Apply()", -1);

  Epetra_Vector x(View,X,0);
  Epetra_Vector xtmp(B_->DomainMap(),false);
  Epetra_Vector xtmp2(x.Map(),false);
#if 0
  // make X (residual) satisfy constraints 
  // do X = (I-BW^T)X
  WT_->Multiply(false,x,xtmp);
  B_->Multiply(false,xtmp,xtmp2);
  x.Update(-1.0,xtmp2,1.0);
#endif

  // create a Space
  const Epetra_BlockMap& bmap = X.Map();
  Space space;
  space.Reshape(bmap.NumGlobalElements(),bmap.NumMyElements(),bmap.MyGlobalElements());
  
  // create input/output mlapi multivectors
  MultiVector b_f(space,1,false);
  MultiVector x_f(space,1,false);
  const int nele = X.Map().NumMyElements();
  for (int i=0; i<nele; ++i)
  {
    x_f(i) = Y[0][i];
    b_f(i) = X[0][i];
  }
  
  // call AMG
  MultiLevelSA(b_f,x_f,0);
  
  // copy solution back
  for (int i=0; i<nele; ++i)
    Y[0][i] = x_f(i);
  
#if 0
  // make Y (search direction) satisfy constraints
  // do Y = (I-WB^T)Y
  Epetra_Vector y(View,Y,0);
  B_->Multiply(true,y,xtmp);
  WT_->Multiply(true,xtmp,xtmp2);
  y.Update(-1.0,xtmp2,1.0);
#endif

#if 0
  Epetra_MultiVector Ytmp(B_->DomainMap(),1,false);
  B_->Multiply(true,Y,Ytmp);
  cout << Ytmp; exit(0);
#endif
  return 0;
}

/*----------------------------------------------------------------------*
 |  apply multigrid linear preconditioner (private)          m.gee 03/06|
 *----------------------------------------------------------------------*/
int MOERTEL::Mortar_ML_Preconditioner::MultiLevelSA(
                                            const MultiVector& b_f,
                                                  MultiVector& x_f, int level) const 
{
  if (level == maxlevels_-1)
  {
    x_f = S(level) * b_f;
    return 0;
  }
  
  MultiVector r_f(P(level).GetRangeSpace());
  MultiVector r_c(P(level).GetDomainSpace());
  MultiVector z_c(P(level).GetDomainSpace());
  
  // pre-smoothing
  x_f = S(level) * b_f;
  //x_f = ImWBT(level) * x_f;
  
  // new residual
  r_f = b_f - A(level) * x_f;
  
  // restrict
  r_c = R(level) * r_f;
  
  // solve coarse problem
  MultiLevelSA(r_c,z_c,level+1);
  
  // prolongate
  x_f = x_f + P(level) * z_c;
  
  // post-smooth
  S(level).Apply(b_f,x_f);
  //x_f = ImWBT(level) * x_f;

  
  return 0;
}

/*----------------------------------------------------------------------*
 |  compute the preconditioner (public)                      m.gee 03/06|
 *----------------------------------------------------------------------*/
bool MOERTEL::Mortar_ML_Preconditioner::Compute()
{

  iscomputed_ = false;
  
  MLAPI::Init();
  
  // get parameters
  int     maxlevels     = mlparams_.get("max levels",10);
  int     maxcoarsesize = mlparams_.get("coarse: max size",10);
  double* nullspace     = mlparams_.get("null space: vectors",(double*)NULL);
  int     nsdim         = mlparams_.get("null space: dimension",1);
  int     numpde        = mlparams_.get("PDE equations",1);
  double  damping       = mlparams_.get("aggregation: damping factor",1.33);
  string  eigenanalysis = mlparams_.get("eigen-analysis: type", "Anorm");
  string  smoothertype  = mlparams_.get("smoother: type","symmetric Gauss-Seidel");
  string  coarsetype    = mlparams_.get("coarse: type","Amesos-KLU");
  
  Space space(A_->RowMatrixRowMap());
  Operator mlapiA(space,space,A_.get(),false);
  Operator mlapiAtilde(space,space,Atilde_.get(),false);
  
  // make the multiplication of BWT
  RefCountPtr<Epetra_CrsMatrix> BWT = rcp(MOERTEL::MatMatMult(*B_,false,*WT_,false,0));
  RefCountPtr<Epetra_CrsMatrix> tmp = rcp(new Epetra_CrsMatrix(Copy,BWT->RowMap(),5,false));
  for (int i=0; i<tmp->NumMyRows(); ++i)
  {
    int grid = tmp->GRID(i);
    double zero = 0.0;
    int err = tmp->InsertGlobalValues(grid,1,&zero,&grid);
    if (err<0) ML_THROW("Padding of BWT matrix failed", -1);
  }
  MOERTEL::MatrixMatrixAdd(*BWT,false,1.0,*tmp,1.0);
  tmp->FillComplete(BWT->DomainMap(),BWT->RangeMap());
  BWT = tmp;
  tmp = null;

  Operator mlapiBWT(space,space,BWT.get(),false);
  
  mlapiImBWT_.resize(maxlevels);                       
  mlapiImWBT_.resize(maxlevels);                       
  mlapiRmod_.resize(maxlevels);                       
  mlapiPmod_.resize(maxlevels);                       
  mlapiAtilde_.resize(maxlevels);
  mlapiS_.resize(maxlevels);
  
  // build nullspace;
  MultiVector NS;
  MultiVector NextNS;
  
  NS.Reshape(mlapiA.GetRangeSpace(),nsdim);
  if (nullspace)
  {
    for (int i=0; i<nsdim; ++i)
      for (int j=0; j<NS.GetMyLength(); ++j)
        NS(j,i) = nullspace[i*NS.GetMyLength()+j];
  }
  else
  {
    if (numpde==1) NS = 1.0;
    else
    {
      NS = 0.0;
      for (int i=0; i<NS.GetMyLength(); ++i)
        for (int j=0; j<numpde; ++j)
          if ( i % numpde == j)
            NS(i,j) = 1.0;
    }
  }

  double lambdamax;
  Operator Ptent;
  Operator P;
  Operator Rtent;
  Operator R;
  Operator IminusA;
  Operator C;

  Operator Pmod;
  Operator Rmod;
  Operator ImBWTfine;
  Operator ImBWTcoarse;
  Operator mlapiBWTcoarse;
  InverseOperator S;

  mlapiAtilde_[0] = mlapiAtilde;

  int level;
  for (level=0; level<maxlevels-1; ++level)
  {
    // this level's operator
    mlapiAtilde = mlapiAtilde_[level];

    // build smoother
    S.Reshape(mlapiAtilde,smoothertype,mlparams_);
    
    if (level) mlparams_.set("PDE equations", NS.GetNumVectors());
    
  
    if (Comm().MyPID()==0)
    {
      ML_print_line("-", 80);
      cout << "MOERTEL/ML : creating level " << level+1 << endl;
      ML_print_line("-", 80);
      fflush(stdout);
    }

    mlparams_.set("workspace: current level",level);
    GetPtent(mlapiA,mlparams_,NS,Ptent,NextNS);
    NS = NextNS;
    fflush(stdout);
    
    if (damping)
    {
      if (eigenanalysis == "Anorm")
        lambdamax = MaxEigAnorm(mlapiA,true);
      else if (eigenanalysis == "cg")
        lambdamax = MaxEigCG(mlapiA,true);
      else if (eigenanalysis == "power-method")
        lambdamax = MaxEigPowerMethod(mlapiA,true);
      else ML_THROW("incorrect parameter (" + eigenanalysis + ")", -1);
    
      IminusA = GetJacobiIterationOperator(mlapiA,damping/lambdamax);
      P = IminusA * Ptent;
    }
    else
    {
      P = Ptent;
      lambdamax = -1.0;
    }
    
    R = GetTranspose(P);
    if (damping)
      Rtent = GetTranspose(Ptent);
    else
      Rtent = R;
      
    // variational coarse grid
    C = GetRAP(R,mlapiA,P);
    
    // compute fine mortar projection operator
    ImBWTfine = GetIdentity(mlapiA.GetDomainSpace(),mlapiA.GetRangeSpace());
    ImBWTfine = ImBWTfine - mlapiBWT;
    
    // compute fine mortar projection operator
    mlapiBWTcoarse = GetRAP(Rtent,mlapiBWT,Ptent);
    ImBWTcoarse = GetIdentity(C.GetDomainSpace(),C.GetRangeSpace());
    ImBWTcoarse = ImBWTcoarse - mlapiBWTcoarse;
    
    // make modified restriction/prolongation
    Rmod = ImBWTcoarse * R * ImBWTfine;
    //Rmod = R * ImBWTfine; 
    //Rmod = R; 
    Pmod = GetTranspose(Rmod);
    
    // store original matrix for construction of next level
    mlapiA = C;
    
    // make final coarse grid operator
    C = GetRAP(Rmod,mlapiAtilde,Pmod);
    
    // store values
    mlapiImBWT_[level]    = ImBWTfine;
    mlapiImBWT_[level+1]  = ImBWTcoarse;
    mlapiImWBT_[level]    = GetTranspose(ImBWTfine);
    mlapiImWBT_[level+1]  = GetTranspose(ImBWTcoarse);
    mlapiRmod_[level]     = Rmod;  
    mlapiPmod_[level]     = Pmod;
    mlapiAtilde_[level+1] = C;
    mlapiS_[level]        = S;

    // prepare for next level
    mlapiBWT = mlapiBWTcoarse;
    
    // break if coarsest level is below specified size
    if (C.GetNumGlobalRows() <= maxcoarsesize)
    {
      ++level;
      break;
    }
  
  } // for (level=0; level<maxlevels-1; ++level)
  
  // set coarse solver
  S.Reshape(mlapiAtilde_[level],coarsetype,mlparams_);
  mlapiS_[level] = S;
  
  // store number of levels
  maxlevels_ = level+1;

  iscomputed_ = true;
  return true;
}
