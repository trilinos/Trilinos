// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef ANASAZI_EPETRA_ADAPTER_HPP
#define ANASAZI_EPETRA_ADAPTER_HPP

#include "AnasaziMultiVec.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziReturnType.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"

namespace Anasazi {

//--------template class AnasaziEpetraVec-------------------------------------
class EpetraVec : public MultiVec<double>, public Epetra_MultiVector {
public:
// constructors
	EpetraVec(const Epetra_BlockMap&, double *, const int, const int stride=0);
	EpetraVec(const Epetra_BlockMap&, const int);
	EpetraVec(Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, int index[], int NumVecs );
	EpetraVec(const Epetra_MultiVector & P_vec);
	~EpetraVec();
	//
	//  member functions inherited from Anasazi::MultiVec
	//
	//  the following is a virtual copy constructor returning
	//  a pointer to the pure virtual class. vector values are
	//  not copied; instead a new MultiVec is created containing
	//  a non-zero amount of columns.
	//
	MultiVec<double> * Clone ( const int numvecs );
	//
	//  the following is a virtual copy constructor returning
	//  a pointer to the pure virtual class. vector values are
	//  copied and a new stand-alone MultiVector is created.
	//  (deep copy).
	//
	MultiVec<double> * CloneCopy ();
	//
	//  the following is a virtual copy constructor returning
	//  a pointer to the pure virtual class. vector values are
	//  copied and a new stand-alone MultiVector is created
	//  where only selected columns are chosen.  (deep copy).
	//
	MultiVec<double> * CloneCopy ( int index[], int numvecs );
	//
	//  the following is a virtual view constructor returning
	//  a pointer to the pure virtual class. vector values are 
	//  shared and hence no memory is allocated for the columns.
	//
	MultiVec<double> * CloneView ( int index[], int numvecs);
	//
	//  this routine sets a subblock of the multivector, which
	//  need not be contiguous, and is given by the indices.
	//
	void SetBlock ( MultiVec<double>& A, int index[], int numvecs );
	//
	int GetNumberVecs () const { return NumVectors(); }
	int GetVecLength () const { return MyLength(); }
	//
	// *this <- alpha * A * B + beta * (*this)
	//
	void MvTimesMatAddMv ( double alpha, MultiVec<double>& A, 
		Teuchos::SerialDenseMatrix<int,double>& B, double beta );
	//
	// *this <- alpha * A + beta * B
	//
	void MvAddMv ( double alpha , MultiVec<double>& A, double beta,
		MultiVec<double>& B);
	//
	// B <- alpha * A^T * (*this)
	//
	void MvTransMv ( double alpha, MultiVec<double>& A, Teuchos::SerialDenseMatrix<int,double>& B );
        //
        // b[i] = A[i]^T * this[i]
        // 
        void MvDot ( MultiVec<double>& A, double b[] );
	//
	// alpha[i] = norm of i-th column of (*this)
	//	
        void MvNorm ( double * normvec ) {
	  if (normvec)
	    assert( Norm2(normvec) == 0 );
	};
  	//
	// random vectors in i-th column of (*this)
	//
	void MvRandom() { assert( Random() == 0 ); };
        //
        // initializes each element of (*this) with alpha
        //
        void MvInit ( double alpha );
	//
	// print (*this)
	//
	void MvPrint();
private:
};
//-------------------------------------------------------------

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


EpetraVec::EpetraVec(const Epetra_BlockMap& Map, double * array, 
		   				const int numvec, const int stride)
	: Epetra_MultiVector(Copy, Map, array, stride, numvec) 
{
}


EpetraVec::EpetraVec(const Epetra_BlockMap& Map, const int numvec)
	: Epetra_MultiVector(Map, numvec) 
{
}


EpetraVec::EpetraVec(Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, 
						int index[], int NumVecs )
	: Epetra_MultiVector(CV, P_vec, index, NumVecs) 
{
}


EpetraVec::EpetraVec(const Epetra_MultiVector& P_vec)
	: Epetra_MultiVector(P_vec) 
{
}


EpetraVec::~EpetraVec() 
{
}
//
//  member functions inherited from Anasazi::MultiVec
//
//
//  Simulating a virtual copy constructor. If we could rely on the co-variance
//  of virtual functions, we could return a pointer to EpetraVec
//  (the derived type) instead of a pointer to the pure virtual base class.
//

MultiVec<double>* EpetraVec::Clone ( const int NumVecs ) 
{
	EpetraVec * ptr_apv = new EpetraVec(Map(),NumVecs);
	return ptr_apv; // safe upcast.
}
//
//  the following is a virtual copy constructor returning
//  a pointer to the pure virtual class. vector values are
//  copied.
//

MultiVec<double>* EpetraVec::CloneCopy() 
{
	EpetraVec *ptr_apv = new EpetraVec(*this);
	return ptr_apv; // safe upcast
}


MultiVec<double>* EpetraVec::CloneCopy ( int index[], int numvecs ) 
{
	EpetraVec * ptr_apv = new EpetraVec(Copy, *this, index, numvecs );
	return ptr_apv; // safe upcast.
}


MultiVec<double>* EpetraVec::CloneView ( int index[], int numvecs ) 
{
	EpetraVec * ptr_apv = new EpetraVec(View, *this, index, numvecs );
	return ptr_apv; // safe upcast.
}


void EpetraVec::SetBlock(MultiVec<double>& A, int index[], int numvecs ) 
{	
	int i,j,ind;
	EpetraVec *A_vec = dynamic_cast<EpetraVec *>(&A); assert(A_vec!=NULL);
	int MyNumVecs = (*this).GetNumberVecs();
	int VecLength = A.GetVecLength();

	// Set the vector values in the right order, careful that the index
	// doesn't go beyond the bounds of the multivector
	for ( j=0; j< numvecs; j++) {
		ind = index[j];
		if (ind < MyNumVecs) {
			for ( i=0; i<VecLength; i++) {
				(*this)[ind][i] = (*A_vec)[j][i];	
			}
		}
	}
}								
//
// *this <- alpha * A * B + beta * (*this)
//

void EpetraVec::MvTimesMatAddMv ( double alpha, MultiVec<double>& A, 
						   Teuchos::SerialDenseMatrix<int,double>& B, double beta ) 
{
	int info=0;
	const int izero=0;
	char* trans="N";
	Epetra_LocalMap LocalMap(B.numRows(), izero, Map().Comm());
	Epetra_MultiVector B_Pvec(Copy, LocalMap, B.values(), B.stride(), B.numCols());

	EpetraVec *A_vec = dynamic_cast<EpetraVec *>(&A); assert(A_vec!=NULL);

	info = Multiply( *trans, *trans, alpha, *A_vec, B_Pvec, beta );

	assert(info==0);
}
//
// *this <- alpha * A + beta * B
//

void EpetraVec::MvAddMv ( double alpha , MultiVec<double>& A, 
						   double beta, MultiVec<double>& B) 
{
	int info=0;
	const double zero = 0.0;

	EpetraVec *A_vec = dynamic_cast<EpetraVec *>(&A); assert(A_vec!=NULL);
	EpetraVec *B_vec = dynamic_cast<EpetraVec *>(&B); assert(B_vec!=NULL);

	info = Update( alpha, *A_vec, beta, *B_vec, zero ); assert(info==0);
}
//
// dense B <- alpha * A^T * (*this)
//

void EpetraVec::MvTransMv ( double alpha, MultiVec<double>& A,
						   Teuchos::SerialDenseMatrix<int,double>& B) 
{
	int info=0;
	const int izero=0;
	const double zero=0.0;
	//const double one=1.0;
	char* trans1="T";
	char* trans2="N";

	EpetraVec *A_vec = dynamic_cast<EpetraVec *>(&A);

	if (A_vec) {

		Epetra_LocalMap LocalMap(B.numRows(), izero, Map().Comm());
		Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());

		info = B_Pvec.Multiply( *trans1, *trans2, alpha, *A_vec, *this, zero ); 
		assert(info==0);
	}
}

//
// b[i] = A[i]^T * this[i]
// 

void EpetraVec::MvDot ( MultiVec<double>& A, double b[] )
{
  EpetraVec *A_vec = dynamic_cast<EpetraVec *>(&A); assert(A_vec!=NULL);
  if (b) {
    assert( this->Dot( *A_vec, b ) == 0 );
  }
}
//
// initializes each element of (*this) with alpha
//


void EpetraVec::MvInit( double alpha )
{	
	int i,j;
	int MyNumVecs = (*this).GetNumberVecs();
	int MyVecLength = (*this).GetVecLength();

	// Set the vector values in the right order, careful that the index
	// doesn't go beyond the bounds of the multivector
	for ( j=0; j< MyNumVecs; j++) {
		for ( i=0; i<MyVecLength; i++) {
			(*this)[j][i] = alpha;	
		}
	}
}								
//
//  print multivectors
//

void EpetraVec::MvPrint() 
{
	cout << *this << endl;
}

///////////////////////////////////////////////////////////////
//--------template class AnasaziEpetraOp---------------------

class EpetraOp : public virtual Operator<double> {
public:
  EpetraOp(const Epetra_Operator& );
  ~EpetraOp();
  ReturnType Apply ( const MultiVec<double>& x, 
		     MultiVec<double>& y ) const;
private:
  const Epetra_Operator & Epetra_Op;
};
//-------------------------------------------------------------
//
// implementation of the AnasaziEpetraOp class.
//
////////////////////////////////////////////////////////////////////
//
// AnasaziOperator constructors
//

EpetraOp::EpetraOp(const Epetra_Operator& Op) 
  : Epetra_Op(Op)
{
}

EpetraOp::~EpetraOp() 
{
}
//
// AnasaziOperator applications
//
ReturnType EpetraOp::Apply ( const MultiVec<double>& x, 
			     MultiVec<double>& y ) const 
{
  //
  // This standard operator computes y = A*x
  //
  int info=0;
  MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
  Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
  
  assert( vec_x!=NULL && vec_y!=NULL );
  //
  // Need to cast away constness because the member function Apply is not declared const.  
  // Change the transpose setting for the operator if necessary and change it back when done.
  //
  info = const_cast<Epetra_Operator&>(Epetra_Op).Apply( *vec_x, *vec_y );
  
  if (info==0) { 
    return Ok; 
  } else { 
    return Failed; 
  }	
}

///////////////////////////////////////////////////////////////
//--------template class AnasaziEpetraGenOp---------------------

class EpetraGenOp : public virtual Operator<double> {
public:
  EpetraGenOp(const Epetra_Operator&, const Epetra_Operator& );
  ~EpetraGenOp();
  ReturnType Apply ( const MultiVec<double>& x, 
		     MultiVec<double>& y ) const;
private:
  const Epetra_Operator & Epetra_AOp;
  const Epetra_Operator & Epetra_BOp;
};
//-------------------------------------------------------------
//
// implementation of the AnasaziEpetraGenOp class.
//
////////////////////////////////////////////////////////////////////
//
// AnasaziOperator constructors
//

EpetraGenOp::EpetraGenOp(const Epetra_Operator& AOp,
			 const Epetra_Operator& BOp) 
  : Epetra_AOp(AOp), Epetra_BOp(BOp) 
{
}


EpetraGenOp::~EpetraGenOp() 
{
}
//
// AnasaziOperator applications
//
ReturnType EpetraGenOp::Apply ( const MultiVec<double>& x, 
				MultiVec<double>& y ) const 
{
  //
  // This generalized operator computes y = A*B*x of y = (A*B)^T*x
  //
  int info=0;
  MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
  Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
  Epetra_MultiVector temp_y(*vec_y); 
  
  assert( vec_x!=NULL && vec_y!=NULL );
  //
  // Need to cast away constness because the member function Apply is not declared const.  
  // Change the transpose setting for the operator if necessary and change it back when done.
  //
  // Apply B
  info=const_cast<Epetra_Operator&>(Epetra_BOp).Apply( *vec_x, temp_y );
  assert(info==0);
  // Apply A
  info=const_cast<Epetra_Operator&>(Epetra_AOp).Apply( temp_y, *vec_y );
  if (info==0) { 
    return Ok; 
  } else { 
    return Failed; 
  }	
}

///////////////////////////////////////////////////////////////
//--------template class AnasaziEpetraSymOp---------------------
class EpetraSymOp : public virtual Operator<double> {
public:
  EpetraSymOp(const Epetra_Operator& Op );
  ~EpetraSymOp();
  ReturnType Apply ( const MultiVec<double>& x, 
		     MultiVec<double>& y ) const;
private:
  const Epetra_Operator& Epetra_Op;
};
//-------------------------------------------------------------
//
// implementation of the AnasaziEpetraSymOp class.
//
////////////////////////////////////////////////////////////////////
//
// AnasaziOperator constructors
//
EpetraSymOp::EpetraSymOp(const Epetra_Operator& Op) 
  : Epetra_Op(Op)
{
}

EpetraSymOp::~EpetraSymOp() 
{
}
//
// AnasaziOperator applications
//
ReturnType EpetraSymOp::Apply ( const MultiVec<double>& x, 
				MultiVec<double>& y ) const 
{
  int info=0;
  MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
  Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
  Epetra_MultiVector* temp_vec = new Epetra_MultiVector( Epetra_Op.OperatorRangeMap(), vec_x->NumVectors() );
  
  assert( vec_x!=NULL && vec_y!=NULL && temp_vec!=NULL );
  //
  // Need to cast away constness because the member function Apply
  // is not declared const.
  //
  // Compute A*x
  info=const_cast<Epetra_Operator&>(Epetra_Op).Apply( *vec_x, *temp_vec );
  if (info!=0) { delete temp_vec; return Failed; }
  
  // Transpose the operator
  info=const_cast<Epetra_Operator&>(Epetra_Op).SetUseTranspose( true );
  if (info!=0) { delete temp_vec; return Failed; }
  
  // Compute A^T*(A*x)
  info=const_cast<Epetra_Operator&>(Epetra_Op).Apply( *temp_vec, *vec_y );
  if (info!=0) { delete temp_vec; return Failed; }
  
  // Un-transpose the operator
  info=const_cast<Epetra_Operator&>(Epetra_Op).SetUseTranspose( false );
  delete temp_vec;
  
  if (info==0)
    return Ok; 
  else
    return Failed; 
}

} // end of Anasazi namespace 

#endif 
// end of file ANASAZI_EPETRA_ADAPTER_HPP
