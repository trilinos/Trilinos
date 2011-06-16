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
#include "mrtr_ml_preconditioner.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            m.gee 03/06|
 *----------------------------------------------------------------------*/
MOERTEL::Mortar_ML_Preconditioner::Mortar_ML_Preconditioner(
                                    Teuchos::RCP<Epetra_CrsMatrix> Atilde,
                                    Teuchos::RCP<Epetra_CrsMatrix> A,
                                    Teuchos::RCP<Epetra_CrsMatrix> WT,
                                    Teuchos::RCP<Epetra_CrsMatrix> B,
                                    Teuchos::RCP<Epetra_Map>       Annmap,
                                    Teuchos::ParameterList& mlparams,
                                    bool constructnow) :
iscomputed_(false),
mlparams_(mlparams),
Atilde_(Atilde),
A_(A),
WT_(WT),
B_(B),
Annmap_(Annmap)
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

#if 1 // working version
/*----------------------------------------------------------------------*
 |  apply multigrid linear preconditioner (public)           m.gee 03/06|
 *----------------------------------------------------------------------*/
int MOERTEL::Mortar_ML_Preconditioner::ApplyInverse(
                     const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
{
  if (!iscomputed_)
  {
    MOERTEL::Mortar_ML_Preconditioner& tmp = 
                       const_cast<MOERTEL::Mortar_ML_Preconditioner&>(*this);
    tmp.Compute();
  }

#if 0
  Epetra_Vector x(View,X,0);
  Epetra_Vector xtmp(B_->DomainMap(),false);
  Epetra_Vector xtmp2(x.Map(),false);
  // make X (residual) satisfy constraints 
  // do X = (I-BW^T)X
  WT_->Multiply(false,x,xtmp);
  B_->Multiply(false,xtmp,xtmp2);
  x.Update(-1.0,xtmp2,1.0);
#endif

  // create a Space
  const Epetra_BlockMap& bmap = X.Map();
  MLAPI::Space space;
  space.Reshape(bmap.NumGlobalElements(),bmap.NumMyElements(),bmap.MyGlobalElements());
  
  // create input/output mlapi multivectors
  MLAPI::MultiVector b_f(space,1,false);
  MLAPI::MultiVector x_f(space,1,false);
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
  Epetra_MultiVector Ytmp(B_->DomainMap(),1,true);
  B_->Multiply(true,Y,Ytmp);
  std::stringstream oss;
  oss << Ytmp; throw ReportError(oss);
#endif
  return 0;
}
#endif

#if 1 // working version
/*----------------------------------------------------------------------*
 |  apply multigrid linear preconditioner (private)          m.gee 03/06|
 *----------------------------------------------------------------------*/
int MOERTEL::Mortar_ML_Preconditioner::MultiLevelSA(
                                            const MLAPI::MultiVector& b_f,
                                                  MLAPI::MultiVector& x_f, int level) const 
{
  if (level == maxlevels_-1)
  {
    x_f = S(level) * b_f;
    //x_f = ImWBT(level) * x_f;
    return 0;
  }
  
  MLAPI::MultiVector r_f(P(level).GetRangeSpace(),1,false);
  MLAPI::MultiVector r_c(P(level).GetDomainSpace(),1,false);
  MLAPI::MultiVector z_c(P(level).GetDomainSpace(),1,true);
  
  // pre-smoothing
  x_f = 0;
  S(level).Apply(b_f,x_f);
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
#endif

#if 1 // working version with standard smoothers
/*----------------------------------------------------------------------*
 |  compute the preconditioner (public)                      m.gee 03/06|
 *----------------------------------------------------------------------*/
bool MOERTEL::Mortar_ML_Preconditioner::Compute()
{

  using namespace MLAPI;

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
  string  ptype         = mlparams_.get("prolongator: type","mod_full");
  
  MLAPI::Space space(A_->RowMatrixRowMap());
  MLAPI::Operator mlapiA(space,space,A_.get(),false);
  MLAPI::Operator mlapiAtilde(space,space,Atilde_.get(),false);
  
  // make the multiplication of BWT
  Teuchos::RCP<Epetra_CrsMatrix> BWT = Teuchos::rcp(MOERTEL::MatMatMult(*B_,false,*WT_,false,0));
  Teuchos::RCP<Epetra_CrsMatrix> tmp = Teuchos::rcp(MOERTEL::PaddedMatrix(BWT->RowMap(),0.0,25));
  MOERTEL::MatrixMatrixAdd(*BWT,false,1.0,*tmp,1.0);
  tmp->FillComplete(BWT->DomainMap(),BWT->RangeMap());
  BWT = tmp;
  tmp = Teuchos::null;

  MLAPI::Operator mlapiBWT(space,space,BWT.get(),false);
  
  mlapiImBWT_.resize(maxlevels);                       
  mlapiImWBT_.resize(maxlevels);                       
  mlapiRmod_.resize(maxlevels);                       
  mlapiPmod_.resize(maxlevels);                       
  mlapiAtilde_.resize(maxlevels);
  mlapiS_.resize(maxlevels);
  
  // build nullspace;
  MLAPI::MultiVector NS;
  MLAPI::MultiVector NextNS;
  
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
  MLAPI::Operator Ptent;
  MLAPI::Operator P;
  MLAPI::Operator Rtent;
  MLAPI::Operator R;
  MLAPI::Operator IminusA;
  MLAPI::Operator C;

  MLAPI::Operator Pmod;
  MLAPI::Operator Rmod;
  MLAPI::Operator ImBWTfine;
  MLAPI::Operator ImBWTcoarse;
  MLAPI::Operator mlapiBWTcoarse;
  MLAPI::InverseOperator S;

  mlapiAtilde_[0] = mlapiAtilde;

  int level;
  for (level=0; level<maxlevels-1; ++level)
  {
    // this level's operator
    mlapiAtilde = mlapiAtilde_[level];

    // build smoother
    if (Comm().MyPID()==0)
    {
      ML_print_line("-", 78);
      cout << "MOERTEL/ML : creating smoother level " << level << endl; 
      fflush(stdout);
    }
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
    if (ptype=="mod_full")
      Rmod = ImBWTcoarse * ( R * ImBWTfine ) + mlapiBWTcoarse * ( R * mlapiBWT ); 
    else if (ptype=="mod_middle")
      Rmod = ImBWTcoarse * ( R * ImBWTfine ); 
    else if (ptype=="mod_simple")
      Rmod = R * ImBWTfine; 
    else if (ptype=="original")
      Rmod = R; 
    else
      ML_THROW("incorrect parameter ( " + ptype + " )", -1);
    
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
  if (Comm().MyPID()==0)
  {
    ML_print_line("-", 78);
    cout << "MOERTEL/ML : creating coarse solver level " << level << endl; 
    fflush(stdout);
  }
  S.Reshape(mlapiAtilde_[level],coarsetype,mlparams_);
  mlapiS_[level] = S;
  
  // store number of levels
  maxlevels_ = level+1;

  iscomputed_ = true;
  return true;
}
#endif

#if 0 // experimental version I
/*----------------------------------------------------------------------*
 |  apply multigrid linear preconditioner (public)           m.gee 03/06|
 *----------------------------------------------------------------------*/
int MOERTEL::Mortar_ML_Preconditioner::ApplyInverse(
                     const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
{
  if (!iscomputed_)
  {
    MOERTEL::Mortar_ML_Preconditioner& tmp = 
                       const_cast<MOERTEL::Mortar_ML_Preconditioner&>(*this);
    tmp.Compute();
  }
  
  Epetra_Vector x(View,X,0);
  Epetra_Vector* x1;
  Epetra_Vector* x2;
  MOERTEL::SplitVector(x,*Arrmap_,x1,*Annmap_,x2);

  Epetra_Vector y(View,Y,0);
  Epetra_Vector* y1;
  Epetra_Vector* y2;
  MOERTEL::SplitVector(y,*Arrmap_,y1,*Annmap_,y2);
  
  // create a mlapi space for x1 and x2
  Space space1;
  Space space2;
  const Epetra_BlockMap& map1 = x1->Map();
  const Epetra_BlockMap& map2 = x2->Map();
  space1.Reshape(map1.NumGlobalElements(),map1.NumMyElements(),map1.MyGlobalElements());
  space2.Reshape(map2.NumGlobalElements(),map2.NumMyElements(),map2.MyGlobalElements());
  
  // create input/output vectors in mlapi
  MultiVector b1_f(space1,1,false);
  MultiVector x1_f(space1,1,false);
  MultiVector b2_f(space2,1,false);
  MultiVector x2_f(space2,1,false);
  const int nele1 = map1.NumMyElements();
  const int nele2 = map2.NumMyElements();
  for (int i=0; i<nele1; ++i)
  {
    x1_f(i) = (*y1)[i];
    b1_f(i) = (*x1)[i];
  }
  for (int i=0; i<nele2; ++i)
  {
    x2_f(i) = (*y2)[i];
    b2_f(i) = (*x2)[i];
  }
  
  // call AMG
  MultiLevelSA(b1_f,b2_f,x1_f,x2_f,0);
  
  // copy solution back
  for (int i=0; i<nele1; ++i)
    (*y1)[i] = x1_f(i); 
  for (int i=0; i<nele2; ++i)
    (*y2)[i] = x2_f(i); 
  MOERTEL::MergeVector(*y1,*y2,y);
  
  return 0;
}
#endif

#if 0 // new experimental version I
/*----------------------------------------------------------------------*
 |  apply multigrid linear preconditioner (private)          m.gee 03/06|
 *----------------------------------------------------------------------*/
int MOERTEL::Mortar_ML_Preconditioner::MultiLevelSA(
                                            const MultiVector& b1_f,
                                            const MultiVector& b2_f,
                                                  MultiVector& x1_f, 
                                                  MultiVector& x2_f, 
                                                  int level) const 
{
  MultiVector r1_f(b1_f.GetVectorSpace(),1,false);
  MultiVector z1_f(b1_f.GetVectorSpace(),1,false);
  
  // presmoothing
  x1_f = 0;
  G_.Apply(b1_f,x1_f);
  x2_f = mlapiMT_ * x1_f;
  x2_f.Scale(-1.0);
  
  // compute residual
  r1_f = b1_f - mlapiAhat11_ * x1_f;

  // postsmoothing
  z1_f = 0;
  G_.Apply(r1_f,z1_f);
  x1_f = x1_f + z1_f;
  x2_f = mlapiMT_ * x1_f;
  x2_f.Scale(-1.0);



  return 0;
}
#endif



#if 0 // new experimental version I
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
  string  ptype         = mlparams_.get("prolongator: type","mod_full");
  
  // create the 2 rowmaps
  Arrmap_ = Teuchos::rcp(MOERTEL::SplitMap(A_->RowMap(),*Annmap_));
  Teuchos::RCP<Epetra_Map> map1 = Arrmap_;
  Teuchos::RCP<Epetra_Map> map2 = Annmap_;

  // split Atilde
  //
  //  Atilde11 Atilde12
  //  Atilde21 Atilde22
  //
  Teuchos::RCP<Epetra_CrsMatrix> Atilde11;
  Teuchos::RCP<Epetra_CrsMatrix> Atilde12;
  Teuchos::RCP<Epetra_CrsMatrix> Atilde21;
  Teuchos::RCP<Epetra_CrsMatrix> Atilde22;
  MOERTEL::SplitMatrix2x2(Atilde_,map1,map2,Atilde11,Atilde12,Atilde21,Atilde22);
  Atilde11_ = Atilde11;
  
  // build BWT (padded to full size)
  //
  //  0   Mr Dinv
  //  0    I
  //
  Teuchos::RCP<Epetra_CrsMatrix> BWT = Teuchos::rcp(MOERTEL::MatMatMult(*B_,false,*WT_,false,0));
  Teuchos::RCP<Epetra_CrsMatrix> tmp = Teuchos::rcp(MOERTEL::PaddedMatrix(BWT->RowMap(),0.0,25));
  MOERTEL::MatrixMatrixAdd(*BWT,false,1.0,*tmp,0.0);
  tmp->FillComplete(BWT->DomainMap(),BWT->RangeMap());
  BWT = tmp;
  tmp = Teuchos::null;
  
  // split BWT to obtain M = Mr Dinv
  Teuchos::RCP<Epetra_CrsMatrix> Zero11;
  Teuchos::RCP<Epetra_CrsMatrix> M;
  Teuchos::RCP<Epetra_CrsMatrix> Zero21;
  Teuchos::RCP<Epetra_CrsMatrix> I22;
  MOERTEL::SplitMatrix2x2(BWT,map1,map2,Zero11,M,Zero21,I22);
  M_ = M;
  
  // transpose BWT to get WBT and split again
  tmp = Teuchos::rcp(MOERTEL::PaddedMatrix(BWT->RowMap(),0.0,25));
  MOERTEL::MatrixMatrixAdd(*BWT,true,1.0,*tmp,0.0);
  tmp->FillComplete();
  Teuchos::RCP<Epetra_CrsMatrix> Zero12;
  MOERTEL::SplitMatrix2x2(tmp,map1,map2,Zero11,Zero12,MT_,I22);
  
  // build matrix Ahat11 = Atilde11 + M Atilde22 M^T    
  Teuchos::RCP<Epetra_CrsMatrix> Ahat11 = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*map1,50,false));
  MOERTEL::MatrixMatrixAdd(*Atilde11,false,1.0,*Ahat11,0.0);
  Teuchos::RCP<Epetra_CrsMatrix> tmp1 = Teuchos::rcp(MOERTEL::MatMatMult(*Atilde22,false,*M,true,0));
  Teuchos::RCP<Epetra_CrsMatrix> tmp2 = Teuchos::rcp(MOERTEL::MatMatMult(*M,false,*tmp1,false,0));
  MOERTEL::MatrixMatrixAdd(*tmp2,false,-1.0,*Ahat11,1.0);
  Ahat11->FillComplete();
  Ahat11->OptimizeStorage();
  Ahat11_ = Ahat11;
  
  // build mlapi objects
  Space space1(*map1);
  Space space2(*map2);
  mlapiAtilde11_.Reshape(space1,space1,Atilde11_.get(),false);
  mlapiAhat11_.Reshape(space1,space1,Ahat11_.get(),false);
  mlapiM_.Reshape(space2,space1,M_.get(),false);
  mlapiMT_.Reshape(space1,space2,MT_.get(),false);
  
  // build the smoother G(Atilde11)
  G_.Reshape(mlapiAtilde11_,smoothertype,mlparams_);

  iscomputed_ = true;
  return true;
}
#endif


#if 0 // experimental version II
/*----------------------------------------------------------------------*
 |  apply multigrid linear preconditioner (public)           m.gee 03/06|
 *----------------------------------------------------------------------*/
int MOERTEL::Mortar_ML_Preconditioner::ApplyInverse(
                     const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
{
  if (!iscomputed_)
  {
    MOERTEL::Mortar_ML_Preconditioner& tmp = 
                       const_cast<MOERTEL::Mortar_ML_Preconditioner&>(*this);
    tmp.Compute();
  }
  
#if 0
  Epetra_Vector x(View,X,0);
  Epetra_Vector xtmp(B_->DomainMap(),false);
  Epetra_Vector xtmp2(x.Map(),false);
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
  Epetra_MultiVector Ytmp(B_->DomainMap(),1,true);
  B_->Multiply(true,Y,Ytmp);
  std::stringstream oss;
  oss << Ytmp; throw ReportError(oss);
#endif
  
  return 0;
}
#endif

#if 0 // new experimental version II
/*----------------------------------------------------------------------*
 |  apply multigrid linear preconditioner (private)          m.gee 03/06|
 *----------------------------------------------------------------------*/
int MOERTEL::Mortar_ML_Preconditioner::MultiLevelSA(const MultiVector& b_f,
                                                          MultiVector& x_f, 
                                                          int level) const 
{
  if (level == maxlevels_-1)
  {
    x_f = 0;
    S(level).Apply(b_f,x_f);
    x_f = ImWBT(level) * x_f;
    return 0;
  }

  MultiVector r_f(b_f.GetVectorSpace(),1,false);
  MultiVector z_f(b_f.GetVectorSpace(),1,false);
  MultiVector r_c(P(level).GetDomainSpace(),1,false);
  MultiVector z_c(P(level).GetDomainSpace(),1,false);

  // presmoothing
  x_f = 0;
  S(level).Apply(b_f,x_f);
  x_f = ImWBT(level) * x_f; 


  // compute residual (different operator)
  r_f = b_f - Ahat(level) * x_f;
  
  // restrict
  r_c = R(level) * r_f;
  
  // solve coarser problem
  MultiLevelSA(r_c,z_c,level+1);
  
  // prolongate
  x_f = x_f + P(level) * z_c;
  x_f = ImWBT(level) * x_f; 
  
  // recompute residual using a different operator
  r_f = b_f - Ahat(level) * x_f;
  
  // postsmoothing
  z_f = 0;
  S(level).Apply(r_f,z_f);
  z_f = ImWBT(level) * z_f; 
  x_f = x_f + z_f;

  return 0;
}
#endif



#if 0 // new experimental version II
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
  string  ptype         = mlparams_.get("prolongator: type","mod_full");
  
  // create the missing rowmap Arrmap_ 
  Arrmap_ = Teuchos::rcp(MOERTEL::SplitMap(A_->RowMap(),*Annmap_));
  Teuchos::RCP<Epetra_Map> map1 = Arrmap_;
  Teuchos::RCP<Epetra_Map> map2 = Annmap_;

  // split Atilde
  //
  //  Atilde11 Atilde12
  //  Atilde21 Atilde22
  //
  Teuchos::RCP<Epetra_CrsMatrix> Atilde11;
  Teuchos::RCP<Epetra_CrsMatrix> Atilde12;
  Teuchos::RCP<Epetra_CrsMatrix> Atilde21;
  Teuchos::RCP<Epetra_CrsMatrix> Atilde22;
  MOERTEL::SplitMatrix2x2(Atilde_,map1,map2,Atilde11,Atilde12,Atilde21,Atilde22);
  
  // build Atildesplit
  //
  //  Atilde11  0
  //  0         I
  //
  Atildesplit_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A_->RowMap(),50,false));
  MOERTEL::MatrixMatrixAdd(*Atilde11,false,1.0,*Atildesplit_,0.0);
  Teuchos::RCP<Epetra_CrsMatrix> tmp = Teuchos::rcp(MOERTEL::PaddedMatrix(*map2,1.0,1));
  tmp->FillComplete();
  MOERTEL::MatrixMatrixAdd(*tmp,false,1.0,*Atildesplit_,1.0);
  Atildesplit_->FillComplete();
  Atildesplit_->OptimizeStorage();
  
  // split A
  //
  //  A11 A12
  //  A21 A22
  //
  Teuchos::RCP<Epetra_CrsMatrix> A11;
  Teuchos::RCP<Epetra_CrsMatrix> A12;
  Teuchos::RCP<Epetra_CrsMatrix> A21;
  Teuchos::RCP<Epetra_CrsMatrix> A22;
  MOERTEL::SplitMatrix2x2(A_,map1,map2,A11,A12,A21,A22);
  
  // build Asplit_
  //
  //  A11  0
  //  0    A22
  //
  Asplit_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,A_->RowMap(),50,false));
  MOERTEL::MatrixMatrixAdd(*A11,false,1.0,*Asplit_,0.0);
  MOERTEL::MatrixMatrixAdd(*A22,false,1.0,*Asplit_,1.0);
  Asplit_->FillComplete();
  Asplit_->OptimizeStorage();
  
  // build BWT (padded to full size)
  //
  //  0   Mr Dinv
  //  0    I
  //
  Teuchos::RCP<Epetra_CrsMatrix> BWT = Teuchos::rcp(MOERTEL::MatMatMult(*B_,false,*WT_,false,10));
                                tmp = Teuchos::rcp(MOERTEL::PaddedMatrix(BWT->RowMap(),0.0,25));
  MOERTEL::MatrixMatrixAdd(*BWT,false,1.0,*tmp,0.0);
  tmp->FillComplete(BWT->DomainMap(),BWT->RangeMap());
  BWT = tmp;
  tmp = Teuchos::null;
  
  // split BWT to obtain M = Mr Dinv
  Teuchos::RCP<Epetra_CrsMatrix> Zero11;
  Teuchos::RCP<Epetra_CrsMatrix> M;
  Teuchos::RCP<Epetra_CrsMatrix> Zero21;
  Teuchos::RCP<Epetra_CrsMatrix> I22;
  MOERTEL::SplitMatrix2x2(BWT,map1,map2,Zero11,M,Zero21,I22);
  
  
  // build matrix Ahat11 = Atilde11 - M Atilde22 M^T    
  Teuchos::RCP<Epetra_CrsMatrix> Ahat11 = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*map1,50,false));
  MOERTEL::MatrixMatrixAdd(*Atilde11,false,1.0,*Ahat11,0.0);
  Teuchos::RCP<Epetra_CrsMatrix> tmp1 = Teuchos::rcp(MOERTEL::MatMatMult(*Atilde22,false,*M,true,10));
  Teuchos::RCP<Epetra_CrsMatrix> tmp2 = Teuchos::rcp(MOERTEL::MatMatMult(*M,false,*tmp1,false,10));
  MOERTEL::MatrixMatrixAdd(*tmp2,false,-1.0,*Ahat11,1.0);
  Ahat11->FillComplete();
  Ahat11->OptimizeStorage();
  
  // build matrix Ahat
  //
  //  Ahat11   0   =   Atilde11 - M Atilde22 M^T   0
  //  0        0       0                           0
  //
  Ahat_ = Teuchos::rcp(MOERTEL::PaddedMatrix(A_->RowMap(),0.0,25));
  MOERTEL::MatrixMatrixAdd(*Ahat11,false,1.0,*Ahat_,0.0);
  Ahat_->FillComplete();
  Ahat_->OptimizeStorage();

  
  // build mlapi objects
  Space space(A_->RowMatrixRowMap());
  Operator mlapiAsplit(space,space,Asplit_.get(),false);
  Operator mlapiAtildesplit(space,space,Atildesplit_.get(),false);
  Operator mlapiAhat(space,space,Ahat_.get(),false);
  Operator mlapiBWT(space,space,BWT.get(),false);
  Operator mlapiBWTcoarse;
  Operator ImBWTfine = GetIdentity(space,space) - mlapiBWT;
  Operator ImBWTcoarse;
  Operator Ptent;
  Operator P;
  Operator Pmod;
  Operator Rtent;
  Operator R;
  Operator Rmod;
  Operator IminusA;
  Operator C;
  InverseOperator S;
  
  mlapiAtildesplit_.resize(maxlevels);
  mlapiAhat_.resize(maxlevels);
  mlapiImBWT_.resize(maxlevels);                       
  mlapiImWBT_.resize(maxlevels);                       
  mlapiRmod_.resize(maxlevels);                       
  mlapiPmod_.resize(maxlevels);                       
  mlapiS_.resize(maxlevels);

  mlapiAtildesplit_[0] = mlapiAtildesplit;
  mlapiAhat_[0]        = mlapiAhat;
  mlapiImBWT_[0]       = ImBWTfine;
  mlapiImWBT_[0]       = GetTranspose(ImBWTfine);
  
  
  // build nullspace
  MultiVector NS;
  MultiVector NextNS;
  NS.Reshape(mlapiAsplit.GetRangeSpace(),nsdim);
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
  
  // construct the hierarchy
  int level=0;
  for (level=0; level<maxlevels-1; ++level)
  {
    // this level's smoothing operator
    mlapiAtildesplit = mlapiAtildesplit_[level];

    // build smoother
    if (Comm().MyPID()==0)
    {
      ML_print_line("-", 78);
      cout << "MOERTEL/ML : creating smoother level " << level << endl; 
      fflush(stdout);
    }
    S.Reshape(mlapiAtildesplit,smoothertype,mlparams_);
    
    if (level) mlparams_.set("PDE equations", NS.GetNumVectors());
    
    if (Comm().MyPID()==0)
    {
      ML_print_line("-", 80);
      cout << "MOERTEL/ML : creating level " << level+1 << endl;
      ML_print_line("-", 80);
      fflush(stdout);
    }
    mlparams_.set("workspace: current level",level);
    
    // get tentative prolongator based on decoupled original system
    GetPtent(mlapiAsplit,mlparams_,NS,Ptent,NextNS);
    NS = NextNS;
    
    // do prolongator smoothing
    if (damping)
    {
      if (eigenanalysis == "Anorm")
        lambdamax = MaxEigAnorm(mlapiAsplit,true);
      else if (eigenanalysis == "cg")
        lambdamax = MaxEigCG(mlapiAsplit,true);
      else if (eigenanalysis == "power-method")
        lambdamax = MaxEigPowerMethod(mlapiAsplit,true);
      else ML_THROW("incorrect parameter (" + eigenanalysis + ")", -1);
      
      IminusA = GetJacobiIterationOperator(mlapiAsplit,damping/lambdamax);
      P       = IminusA * Ptent;
      R       = GetTranspose(P);
      Rtent   = GetTranspose(Ptent);
    }
    else
    {
      P     = Ptent;
      Rtent = GetTranspose(Ptent);
      R     = Rtent;
      lambdamax = -1.0;
    }
    
    // do variational coarse grid of split original matrix Asplit
    C = GetRAP(R,mlapiAsplit,P);
    
    // compute the mortar projections on coarse grid
    mlapiBWTcoarse = GetRAP(Rtent,mlapiBWT,Ptent); 
    ImBWTcoarse    = GetIdentity(C.GetDomainSpace(),C.GetRangeSpace());    
    ImBWTcoarse    = ImBWTcoarse - mlapiBWTcoarse;
    
    // do modified prolongation and restriction
    if (ptype=="mod_full")
      Rmod = ImBWTcoarse * ( R * ImBWTfine ) + mlapiBWTcoarse * ( R * mlapiBWT ); 
    else if (ptype=="mod_middle")
      Rmod = ImBWTcoarse * ( R * ImBWTfine ); 
    else if (ptype=="mod_simple")
      Rmod = R * ImBWTfine; 
    else if (ptype=="original")
      Rmod = R; 
    else
      ML_THROW("incorrect parameter ( " + ptype + " )", -1);
    Pmod = GetTranspose(Rmod);
    
    // store matrix for construction of next level
    mlapiAsplit = C;
    
    // make coarse smoothing operator
    // make coarse residual operator
    mlapiAtildesplit_[level+1] = GetRAP(Rmod,mlapiAtildesplit,Pmod);
    mlapiAhat_[level+1]        = GetRAP(Rmod,mlapiAhat_[level],Pmod);
    mlapiImBWT_[level]         = ImBWTfine;
    mlapiImBWT_[level+1]       = ImBWTcoarse;
    mlapiImWBT_[level]         = GetTranspose(ImBWTfine);
    mlapiImWBT_[level+1]       = GetTranspose(ImBWTcoarse);
    mlapiRmod_[level]          = Rmod;
    mlapiPmod_[level]          = Pmod;
    mlapiS_[level]             = S;
    
    // prepare for next level
    mlapiBWT  = mlapiBWTcoarse;
    ImBWTfine = ImBWTcoarse;
    
    // break if coarsest level is below specified size
    if (mlapiAsplit.GetNumGlobalRows() <= maxcoarsesize)
    {
      ++level;
      break;
    }
    
  } // for (level=0; level<maxlevels-1; ++level)
  
  // do coarse smoother
  S.Reshape(mlapiAtildesplit_[level],coarsetype,mlparams_);
  mlapiS_[level] = S;
  
  // store max number of levels
  maxlevels_ = level+1;

  iscomputed_ = true;
  return true;
}
#endif
