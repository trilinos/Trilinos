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

#ifndef ANASAZI_PETRA_INTERFACE_HPP
#define ANASAZI_PETRA_INTERFACE_HPP

#include "AnasaziMultiVec.hpp"
#include "AnasaziMatrix.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziReturnType.hpp"

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"

namespace Anasazi {

//--------template class AnasaziPetraVec-------------------------------------
class PetraVec : public MultiVec<double>, public Epetra_MultiVector {
public:
// constructors
	PetraVec(const Epetra_BlockMap&, double *, const int, const int stride=0);
	PetraVec(const Epetra_BlockMap&, const int);
	PetraVec(Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, int index[], int NumVecs );
	PetraVec(const Epetra_MultiVector & P_vec);
	~PetraVec();
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
	// alpha[i] = norm of i-th column of (*this)
	//
	void MvNorm ( double* normvec);
	//
	// random vectors in i-th column of (*this)
	//
	void MvRandom();
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


PetraVec::PetraVec(const Epetra_BlockMap& Map, double * array, 
		   				const int numvec, const int stride)
	: Epetra_MultiVector(Copy, Map, array, stride, numvec) 
{
}


PetraVec::PetraVec(const Epetra_BlockMap& Map, const int numvec)
	: Epetra_MultiVector(Map, numvec) 
{
}


PetraVec::PetraVec(Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, 
						int index[], int NumVecs )
	: Epetra_MultiVector(CV, P_vec, index, NumVecs) 
{
}


PetraVec::PetraVec(const Epetra_MultiVector& P_vec)
	: Epetra_MultiVector(P_vec) 
{
}


PetraVec::~PetraVec() 
{
}
//
//  member functions inherited from Anasazi::MultiVec
//
//
//  Simulating a virtual copy constructor. If we could rely on the co-variance
//  of virtual functions, we could return a pointer to PetraVec
//  (the derived type) instead of a pointer to the pure virtual base class.
//

MultiVec<double>* PetraVec::Clone ( const int NumVecs ) 
{
	PetraVec * ptr_apv = new PetraVec(Map(),NumVecs);
	return ptr_apv; // safe upcast.
}
//
//  the following is a virtual copy constructor returning
//  a pointer to the pure virtual class. vector values are
//  copied.
//

MultiVec<double>* PetraVec::CloneCopy() 
{
	PetraVec *ptr_apv = new PetraVec(*this);
	return ptr_apv; // safe upcast
}


MultiVec<double>* PetraVec::CloneCopy ( int index[], int numvecs ) 
{
	PetraVec * ptr_apv = new PetraVec(Copy, *this, index, numvecs );
	return ptr_apv; // safe upcast.
}


MultiVec<double>* PetraVec::CloneView ( int index[], int numvecs ) 
{
	PetraVec * ptr_apv = new PetraVec(View, *this, index, numvecs );
	return ptr_apv; // safe upcast.
}


void PetraVec::SetBlock(MultiVec<double>& A, int index[], int numvecs ) 
{	
	int i,j,ind;
	PetraVec *A_vec = dynamic_cast<PetraVec *>(&A); assert(A_vec!=NULL);
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

void PetraVec::MvTimesMatAddMv ( double alpha, MultiVec<double>& A, 
						   Teuchos::SerialDenseMatrix<int,double>& B, double beta ) 
{
	int info=0;
	const int izero=0;
	char* trans="N";
	Epetra_LocalMap LocalMap(B.numRows(), izero, Map().Comm());
	Epetra_MultiVector B_Pvec(Copy, LocalMap, B.values(), B.stride(), B.numCols());

	PetraVec *A_vec = dynamic_cast<PetraVec *>(&A); assert(A_vec!=NULL);

	info = Multiply( *trans, *trans, alpha, *A_vec, B_Pvec, beta );

	assert(info==0);
}
//
// *this <- alpha * A + beta * B
//

void PetraVec::MvAddMv ( double alpha , MultiVec<double>& A, 
						   double beta, MultiVec<double>& B) 
{
	int info=0;
	const double zero = 0.0;

	PetraVec *A_vec = dynamic_cast<PetraVec *>(&A); assert(A_vec!=NULL);
	PetraVec *B_vec = dynamic_cast<PetraVec *>(&B); assert(B_vec!=NULL);

	info = Update( alpha, *A_vec, beta, *B_vec, zero ); assert(info==0);
}
//
// dense B <- alpha * A^T * (*this)
//

void PetraVec::MvTransMv ( double alpha, MultiVec<double>& A,
						   Teuchos::SerialDenseMatrix<int,double>& B) 
{
	int info=0;
	const int izero=0;
	const double zero=0.0;
	//const double one=1.0;
	char* trans1="T";
	char* trans2="N";

	PetraVec *A_vec = dynamic_cast<PetraVec *>(&A);

	if (A_vec) {

		Epetra_LocalMap LocalMap(B.numRows(), izero, Map().Comm());
		Epetra_MultiVector B_Pvec(View, LocalMap, B.values(), B.stride(), B.numCols());

		info = B_Pvec.Multiply( *trans1, *trans2, alpha, *A_vec, *this, zero ); 
		assert(info==0);
	}
}
//
// alpha[i] = norm of i-th column of (*this)
//

void PetraVec::MvNorm ( double * normvec ) 
{
	int info=0;
	if (normvec) {
		info = Norm2(normvec);
		assert(info==0);
	}
}
//
// random vectors in i-th column of (*this)
//

void PetraVec::MvRandom () 
{
	int info=0;
	info = Random();
	assert(info==0);
}
//
// initializes each element of (*this) with alpha
//


void PetraVec::MvInit( double alpha )
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

void PetraVec::MvPrint() 
{
	cout << *this << endl;
}

///////////////////////////////////////////////////////////////
//--------template class AnasaziPetraMat-----------------------
class PetraMat : public virtual Matrix<double> {
public:
	PetraMat(const Epetra_Operator& );
	~PetraMat();
	ReturnType ApplyMatrix ( const MultiVec<double>& x, 
					MultiVec<double>& y ) const;
private:
	const Epetra_Operator & Epetra_Mat;
};
//-------------------------------------------------------------
//
// implementation of the AnasaziPetraMat class.
//
////////////////////////////////////////////////////////////////////
//
// Anasazi::Matrix constructors
//
PetraMat::PetraMat(const Epetra_Operator& Matrix) 
	: Epetra_Mat(Matrix) 
{
}

PetraMat::~PetraMat() 
{
}

//
// Anasazi::Matrix matrix multiply
//
ReturnType PetraMat::ApplyMatrix ( const MultiVec<double>& x, 
						  MultiVec<double>& y ) const 
{
	int info=0;
	MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
	Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
	Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);

	assert( vec_x!=NULL && vec_y!=NULL );
	//
	// Need to cast away constness because the member function Apply
	// is not declared const.
	//
	info=const_cast<Epetra_Operator&>(Epetra_Mat).Apply( *vec_x, *vec_y );
	if (info==0) { 
		return Ok; 
	} else { 
		return Failed; 
	}	
}

///////////////////////////////////////////////////////////////
//--------template class AnasaziPetraStdOp---------------------

class PetraStdOp : public virtual Operator<double> {
public:
	PetraStdOp(const Epetra_Operator& );
	~PetraStdOp();
	ReturnType Apply ( const MultiVec<double>& x, 
					MultiVec<double>& y ) const;
private:
	const Epetra_Operator & Epetra_Op;
};
//-------------------------------------------------------------
//
// implementation of the AnasaziPetraStdOp class.
//
////////////////////////////////////////////////////////////////////
//
// AnasaziOperator constructors
//

PetraStdOp::PetraStdOp(const Epetra_Operator& Op) 
	: Epetra_Op(Op)
{
}

PetraStdOp::~PetraStdOp() 
{
}
//
// AnasaziOperator applications
//
ReturnType PetraStdOp::Apply ( const MultiVec<double>& x, 
						  MultiVec<double>& y ) const 
{
	int info=0;
	MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
	Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
	Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);

	assert( vec_x!=NULL && vec_y!=NULL );
	//
	// Need to cast away constness because the member function Apply
	// is not declared const.
	//
	info=const_cast<Epetra_Operator&>(Epetra_Op).Apply( *vec_x, *vec_y );
	if (info==0) { 
		return Ok; 
	} else { 
		return Failed; 
	}	
}

///////////////////////////////////////////////////////////////
//--------template class AnasaziPetraGenOp---------------------

class PetraGenOp : public virtual Operator<double> {
public:
	PetraGenOp(const Epetra_Operator&, const Epetra_Operator& );
	~PetraGenOp();
	ReturnType Apply ( const MultiVec<double>& x, 
			   MultiVec<double>& y ) const;
private:
	const Epetra_Operator & Epetra_AOp;
	const Epetra_Operator & Epetra_BOp;
};
//-------------------------------------------------------------
//
// implementation of the AnasaziPetraGenOp class.
//
////////////////////////////////////////////////////////////////////
//
// AnasaziOperator constructors
//

PetraGenOp::PetraGenOp(const Epetra_Operator& AOp,
		       const Epetra_Operator& BOp) 
	: Epetra_AOp(AOp), Epetra_BOp(BOp) 
{
}


PetraGenOp::~PetraGenOp() 
{
}
//
// AnasaziOperator applications
//

ReturnType PetraGenOp::Apply ( const MultiVec<double>& x, 
			       MultiVec<double>& y ) const 
{
	int info=0;
	MultiVec<double> & temp_x = const_cast<MultiVec<double> &>(x);
	Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
	Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
	Epetra_MultiVector temp_y(*vec_y); 

	assert( vec_x!=NULL && vec_y!=NULL );
	//
	// Need to cast away constness because the member function Apply
	// is not declared const.
	//
	info=const_cast<Epetra_Operator&>(Epetra_BOp).Apply( *vec_x, temp_y );
	assert(info==0);
	info=const_cast<Epetra_Operator&>(Epetra_AOp).Apply( temp_y, *vec_y );
	if (info==0) { 
		return Ok; 
	} else { 
		return Failed; 
	}	
}

} // end of Anasazi namespace 

#endif 
 // end of file ANASAZI_PETRA_INTERFACE_HPP
