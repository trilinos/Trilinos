// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//                  LOCA Continuation Algorithm Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "Anasazi_LOCA_MultiVec.H"
#include "LOCA_ErrorCheck.H"

Anasazi::LOCAMultiVec::LOCAMultiVec(const NOX::Abstract::Vector& N_vec, 
				    int NumVecs) :
  mvPtrs(NumVecs), CV(Anasazi::Copy)
{
  for (int i=0; i<NumVecs; i++) {
    mvPtrs[i] = N_vec.clone(NOX::ShapeCopy);
    mvPtrs[i]->init(0.0);
  }
}

Anasazi::LOCAMultiVec::LOCAMultiVec(
			      const vector<NOX::Abstract::Vector*> N_vecPtrs,
			      Anasazi::DataAccess type) : 
  mvPtrs(N_vecPtrs.size()), CV(type)
{
  if (type == Anasazi::Copy)
    for (unsigned int i=0; i<mvPtrs.size(); i++)
      mvPtrs[i] = N_vecPtrs[i]->clone(NOX::DeepCopy);

  else
    for (unsigned int i=0; i<mvPtrs.size(); i++)
      mvPtrs[i] = N_vecPtrs[i];
}

Anasazi::LOCAMultiVec::LOCAMultiVec(const Anasazi::LOCAMultiVec& source, 
				    Anasazi::DataAccess type ) : 
  mvPtrs(source.mvPtrs.size()),
  CV(type)
{
  if (type == Anasazi::Copy)
    for (unsigned int i=0; i<mvPtrs.size(); i++)
      mvPtrs[i] = source.mvPtrs[i]->clone(NOX::DeepCopy);

  else
    for (unsigned int i=0; i<mvPtrs.size(); i++)
      mvPtrs[i] = source.mvPtrs[i];
}

Anasazi::LOCAMultiVec::LOCAMultiVec(Anasazi::DataAccess type, 
				    const Anasazi::LOCAMultiVec& source, 
				    const std::vector<int>& index): 
  mvPtrs(index.size()), CV(type)
{
  int i;

  if (type == Anasazi::Copy)
    for (i=0; i<index.size(); i++)
      mvPtrs[i] = source.mvPtrs[ index[i] ]->clone(NOX::DeepCopy);
  
  else
    for (i=0; i<index.size(); i++)
      mvPtrs[i] = source.mvPtrs[ index[i] ];
}

Anasazi::LOCAMultiVec::~LOCAMultiVec()
{
  if (CV == Anasazi::Copy)
    for (unsigned int i=0; i<mvPtrs.size(); i++)
      delete mvPtrs[i];
}

Anasazi::MultiVec<double>* 
Anasazi::LOCAMultiVec::Clone(const int NumVecs) const
{
  return new Anasazi::LOCAMultiVec(*(mvPtrs[0]),NumVecs);
}

Anasazi::MultiVec<double>* 
Anasazi::LOCAMultiVec::CloneCopy() const
{
  return new Anasazi::LOCAMultiVec(*this);
}

Anasazi::MultiVec<double>* 
Anasazi::LOCAMultiVec::CloneCopy( const std::vector<int>& index ) const
{
  return new Anasazi::LOCAMultiVec( Copy, *this, index );
}

Anasazi::MultiVec<double>* 
Anasazi::LOCAMultiVec::CloneView( const std::vector<int>& index ) 
{
  return  new Anasazi::LOCAMultiVec( View, *this, index );
}

int 
Anasazi::LOCAMultiVec::GetNumberVecs() const 
{
  return mvPtrs.size();
}

int 
Anasazi::LOCAMultiVec::GetVecLength() const 
{
  return mvPtrs[0]->length();
}

void 
Anasazi::LOCAMultiVec::SetBlock(const Anasazi::MultiVec<double>& A, 
				const std::vector<int>& index)
{
  int i, ind;
  Anasazi::LOCAMultiVec *A_vec = 
    dynamic_cast<Anasazi::LOCAMultiVec *>(&const_cast<Anasazi::MultiVec<double> &>(A)); assert(A_vec!=NULL);
  int MyNumVecs = mvPtrs.size();
  for (i=0; i<index.size(); i++) {
    ind = index[i];
    if (ind < MyNumVecs) {
      delete mvPtrs[ind];
      mvPtrs[ind] = A_vec->mvPtrs[i]->clone(NOX::DeepCopy);
    }
  }
}

void 
Anasazi::LOCAMultiVec::MvTimesMatAddMv(
			     const double alpha, 
			     const Anasazi::MultiVec<double>& A, 
			     const Teuchos::SerialDenseMatrix<int,double>& B, 
			     const double beta) 
{
  int i,j;
  Anasazi::LOCAMultiVec *A_vec = 
    dynamic_cast<Anasazi::LOCAMultiVec *>(&const_cast<Anasazi::MultiVec<double> &>(A)); assert(A_vec!=NULL);
  int m = B.numRows();
  int n = B.numCols();
  int ldb = B.stride();
  double *Bvals = B.values();  	
  Anasazi::LOCAMultiVec *temp_vec = 
    new Anasazi::LOCAMultiVec(*(mvPtrs[0]),n);
  temp_vec->MvInit(0.0);
  double one = 1.0;
  //
  //	*this <- alpha * A * B + beta *(*this)
  //
  for (j=0; j<n; j++) {
    for (i=0; i<m; i++) {
      temp_vec->mvPtrs[j]->update(Bvals[j*ldb+i], *(A_vec->mvPtrs[i]),one);
    }				
    mvPtrs[j]->update(alpha,*(temp_vec->mvPtrs[j]),beta);
  }
  delete temp_vec;
}

void 
Anasazi::LOCAMultiVec::MvAddMv(const double alpha, 
			       const Anasazi::MultiVec<double>& A, 
			       const double beta, 
			       const Anasazi::MultiVec<double>& B) 
{
  const double zero = 0.0;
  Anasazi::LOCAMultiVec *A_vec = 
    dynamic_cast<Anasazi::LOCAMultiVec *>(&const_cast<Anasazi::MultiVec<double> &>(A)); assert(A_vec!=NULL);
  Anasazi::LOCAMultiVec *B_vec = 
    dynamic_cast<Anasazi::LOCAMultiVec *>(&const_cast<Anasazi::MultiVec<double> &>(B)); assert(B_vec!=NULL);

  for (unsigned int i=0; i<mvPtrs.size(); i++)
    mvPtrs[i]->update(alpha, *(A_vec->mvPtrs[i]), beta, *(B_vec->mvPtrs[i]), 
		      zero);		
}

void 
Anasazi::LOCAMultiVec::MvTransMv(
			     const double alpha, 
			     const Anasazi::MultiVec<double>& A,
			     Teuchos::SerialDenseMatrix<int,double>& B) const
{
  int i,j;
  Anasazi::LOCAMultiVec *A_vec = 
    dynamic_cast<Anasazi::LOCAMultiVec *>(&const_cast<Anasazi::MultiVec<double> &>(A)); assert(A_vec!=NULL);
  int m = B.numRows();
  int n = B.numCols();
  int ldb = B.stride();
  double *Bvals = B.values();  	
  
  for (j=0; j<n; j++)
    for (i=0; i<m; i++)
      Bvals[j*ldb + i] = alpha * mvPtrs[j]->dot(*(A_vec->mvPtrs[i]));
}


void
Anasazi::LOCAMultiVec::MvDot(const Anasazi::MultiVec<double>& A, 
			     std::vector<double>* b ) const
{
  int j;
  int n = mvPtrs.size();
  Anasazi::LOCAMultiVec *A_vec = 
    dynamic_cast<Anasazi::LOCAMultiVec *>(&const_cast<Anasazi::MultiVec<double> &>(A)); assert(A_vec!=NULL);
  if (A_vec && b && ( b->size() >= n )) {
    for (j=0; j<n; j++)
      (*b)[j] = mvPtrs[j]->dot(*(A_vec->mvPtrs[j]));
  }
}


void 
Anasazi::LOCAMultiVec::MvNorm(std::vector<double>* normvec) const 
{
  if (normvec)
    for (unsigned int i=0; i<mvPtrs.size(); i++)
      (*normvec)[i] = mvPtrs[i]->norm();
}

void 
Anasazi::LOCAMultiVec::MvRandom() 
{
  for (unsigned int i=0; i<mvPtrs.size(); i++)
    mvPtrs[i]->random();	
}

void 
Anasazi::LOCAMultiVec::MvInit(double alpha) 
{
  for (unsigned int i=0; i<mvPtrs.size(); i++)
    mvPtrs[i]->init(alpha);	
}

void 
Anasazi::LOCAMultiVec::MvPrint() const
{
}

NOX::Abstract::Vector&
Anasazi::LOCAMultiVec::GetNOXVector(int index) 
{
  if (index < static_cast<int>(mvPtrs.size())) 
    return *(mvPtrs[index]);
  else {
    ::LOCA::ErrorCheck::throwError("Anasazi::LOCAMultiVec::GetNOXVector()",
				   "Invalid index");
    return *(mvPtrs[0]);
  }
}

const NOX::Abstract::Vector&
Anasazi::LOCAMultiVec::GetNOXVector(int index) const
{
  if (index < static_cast<int>(mvPtrs.size())) 
    return *(mvPtrs[index]);
  else {
    ::LOCA::ErrorCheck::throwError("Anasazi::LOCAMultiVec::GetNOXVector()",
				   "Invalid index");
    return *(mvPtrs[0]);
  }
}
