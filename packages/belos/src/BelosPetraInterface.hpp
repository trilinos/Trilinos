// File BelosPetraInterface.hpp: interface for the BelosPetra class.
//
#ifndef BELOS_PETRA_HPP
#define BELOS_PETRA_HPP

#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"

#include "BelosConfigDefs.hpp"
#include "AnasaziMatrix.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziPrecondition.hpp"
#include "AnasaziReturnType.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Belos {

//--------template class BelosPetraVec-------------------------------------
template <class TYPE>
class PetraVec : public Anasazi::MultiVec<TYPE>, public Epetra_MultiVector {
public:
// constructors
	PetraVec(const Epetra_BlockMap&, TYPE *, const int, const int stride=0);
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
	Anasazi::MultiVec<TYPE> * Clone ( const int numvecs );
	//
	//  the following is a virtual copy constructor returning
	//  a pointer to the pure virtual class. vector values are
	//  copied and a new stand-alone MultiVector is created.
	//  (deep copy).
	//
	Anasazi::MultiVec<TYPE> * CloneCopy ();
	//
	//  the following is a virtual copy constructor returning
	//  a pointer to the pure virtual class. vector values are
	//  copied and a new stand-alone MultiVector is created
	//  where only selected columns are chosen.  (deep copy).
	//
	Anasazi::MultiVec<TYPE> * CloneCopy ( int index[], int numvecs );
	//
	//  the following is a virtual view constructor returning
	//  a pointer to the pure virtual class. vector values are 
	//  shared and hence no memory is allocated for the columns.
	//
	Anasazi::MultiVec<TYPE> * CloneView ( int index[], int numvecs );
	//
	//  this routine sets a subblock of the multivector, which
	//  need not be contiguous, and is given by the indices.
	//
	void SetBlock ( Anasazi::MultiVec<TYPE>& A, int index[], int numvecs );
	//
	int GetNumberVecs () const;
	int GetVecLength () const;
	//
	// *this <- alpha * A * B + beta * (*this)
	//
	void MvTimesMatAddMv ( TYPE alpha, Anasazi::MultiVec<TYPE>& A, 
		Teuchos::SerialDenseMatrix<int,TYPE>& B, TYPE beta );
	//
	// *this <- alpha * A + beta * B
	//
	void MvAddMv ( TYPE alpha , Anasazi::MultiVec<TYPE>& A, TYPE beta,
		Anasazi::MultiVec<TYPE>& B);
	//
	// B <- alpha * A^T * (*this)
	//
	void MvTransMv ( TYPE alpha, Anasazi::MultiVec<TYPE>& A, Teuchos::SerialDenseMatrix<int,TYPE>& B );
	//
	// alpha[i] = norm of i-th column of (*this)
	//
	void MvNorm ( TYPE* normvec);
	//
	// random vectors in i-th column of (*this)
	//
	void MvRandom();
        //
        // initializes each element of (*this) with alpha
        //
        void MvInit ( TYPE alpha );
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

template<class TYPE>
PetraVec<TYPE>::PetraVec(const Epetra_BlockMap& Map, TYPE * array, 
					   	const int numvec, const int stride)
	: Epetra_MultiVector(Copy, Map, array, stride, numvec) 
{
}

template<class TYPE>
PetraVec<TYPE>::PetraVec(const Epetra_BlockMap& Map, const int numvec)
	: Epetra_MultiVector(Map, numvec) 
{
}

template<class TYPE>
PetraVec<TYPE>::PetraVec(Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, 
						int index[], int NumVecs ) 
	: Epetra_MultiVector(CV, P_vec, index, NumVecs) 
{
}

template<class TYPE>
PetraVec<TYPE>::PetraVec(const Epetra_MultiVector& P_vec) 
	: Epetra_MultiVector(P_vec) 
{
}

template<class TYPE>
PetraVec<TYPE>::~PetraVec() 
{
}
//
//  member functions inherited from Anasazi::MultiVec
//
//
//  Simulating a virtual copy constructor. If we could rely on the co-variance
//  of virtual functions, we could return a pointer to BelosPetraVec<TYPE>
//  (the derived type) instead of a pointer to the pure virtual base class.
//
template<class TYPE>
Anasazi::MultiVec<TYPE>* PetraVec<TYPE>::Clone ( const int NumVecs ) {
	PetraVec * ptr_apv = new PetraVec(Map(),NumVecs);
	return ptr_apv; // safe upcast.
}
	//
	//  the following is a virtual copy constructor returning
	//  a pointer to the pure virtual class. vector values are
	//  copied.
	//
template<class TYPE>
Anasazi::MultiVec<TYPE>* PetraVec<TYPE>::CloneCopy() {
	PetraVec *ptr_apv = new PetraVec(*this);
	return ptr_apv; // safe upcast
}

template<class TYPE>
Anasazi::MultiVec<TYPE>* PetraVec<TYPE>::CloneCopy ( int index[], int numvecs ) {
	PetraVec * ptr_apv = new PetraVec(Copy, *this, index, numvecs );
	return ptr_apv; // safe upcast.
}

template<class TYPE>
Anasazi::MultiVec<TYPE>* PetraVec<TYPE>::CloneView ( int index[], int numvecs ) {
	PetraVec * ptr_apv = new PetraVec(View, *this, index, numvecs );
	return ptr_apv; // safe upcast.
}

template<class TYPE>
void PetraVec<TYPE>::SetBlock(Anasazi::MultiVec<TYPE>& A, int index[], int numvecs ) 
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
		
template<class TYPE>
int PetraVec<TYPE>::GetNumberVecs () const {
	return NumVectors();
}

template<class TYPE>
int PetraVec<TYPE>::GetVecLength () const {
	return MyLength();
}
//
// *this <- alpha * A * B + beta * (*this)
//
template<class TYPE>
void PetraVec<TYPE>::MvTimesMatAddMv ( TYPE alpha, Anasazi::MultiVec<TYPE>& A, 
						   Teuchos::SerialDenseMatrix<int,TYPE>& B, TYPE beta ) {
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
template<class TYPE>
void PetraVec<TYPE>::MvAddMv ( TYPE alpha , Anasazi::MultiVec<TYPE>& A, 
						   TYPE beta, Anasazi::MultiVec<TYPE>& B) {
	int info=0;
	const TYPE zero = 0.0;

	PetraVec *A_vec = dynamic_cast<PetraVec *>(&A); assert(A_vec!=NULL);
	PetraVec *B_vec = dynamic_cast<PetraVec *>(&B); assert(B_vec!=NULL);

	info = Update( alpha, *A_vec, beta, *B_vec, zero ); assert(info==0);
}
//
// dense B <- alpha * A^T * (*this)
//
template<class TYPE>
void PetraVec<TYPE>::MvTransMv ( TYPE alpha, Anasazi::MultiVec<TYPE>& A,
						   Teuchos::SerialDenseMatrix<int,TYPE>& B) {
	int info=0;
	const int izero=0;
	const TYPE zero=0.0;
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
template<class TYPE>
void PetraVec<TYPE>::MvNorm ( TYPE * normvec ) {
	int info=0;
	if (normvec) {
		info = Norm2(normvec);
		assert(info==0);
	}
}
//
// random vectors in i-th column of (*this)
//
template<class TYPE>
void PetraVec<TYPE>::MvRandom () {
	int info=0;
	info = Random();
	assert(info==0);
}
//
// initializes each element of (*this) with alpha
//

template<class TYPE>
void PetraVec<TYPE>::MvInit( TYPE alpha )
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
template<class TYPE>
void PetraVec<TYPE>::MvPrint() {
	cout << *this << endl;
}

///////////////////////////////////////////////////////////////
//--------template class PetraMat-----------------------
template <class TYPE> 
class PetraMat : public Anasazi::Matrix<TYPE> {
public:
	PetraMat(const Epetra_Operator& Matrix);
	~PetraMat();
	Anasazi::ReturnType ApplyMatrix ( const Anasazi::MultiVec<TYPE>& x, Anasazi::MultiVec<TYPE>& y ) const;
	const Epetra_Operator & GetMatrix () { return(Epetra_Mat); };
private:
	const Epetra_Operator & Epetra_Mat;
};
//-------------------------------------------------------------
//
// implementation of the Belos::PetraMat class.
//
////////////////////////////////////////////////////////////////////
//
// Anasazi::Matrix constructors
//
template <class TYPE>
PetraMat<TYPE>::PetraMat(const Epetra_Operator& Matrix) 
	: Epetra_Mat(Matrix) 
{
}

template <class TYPE>
PetraMat<TYPE>::~PetraMat() 
{
}
//
// Anasazi::Matrix matrix multiply
//
template <class TYPE>
Anasazi::ReturnType PetraMat<TYPE>::ApplyMatrix ( const Anasazi::MultiVec<TYPE>& x, 
					  Anasazi::MultiVec<TYPE>& y ) const {
	int info=0;
	Anasazi::MultiVec<TYPE> & temp_x = const_cast<Anasazi::MultiVec<TYPE> &>(x);
	Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
	Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);

	assert( vec_x!=NULL && vec_y!=NULL );
	//
	// Need to cast away constness because the member function Multiply
	// is not declared const.
	//
	info=const_cast<Epetra_Operator&>(Epetra_Mat).Apply( *vec_x, *vec_y );
	assert(info==0);
	return Anasazi::Ok;	
}
 
///////////////////////////////////////////////////////////////
//--------template class Belos::PetraPrec--------------------

template <class TYPE>
class PetraPrec : public Anasazi::Precondition<TYPE> {
public:
        PetraPrec(const Epetra_Operator& Prec);
        ~PetraPrec();
        void ApplyPrecondition ( const Anasazi::MultiVec<TYPE>& x, Anasazi::MultiVec<TYPE>& y ) const;
private:
   	const Epetra_Operator& Epetra_Prec;
};
//--------------------------------------------------------------
//
// implementation of the Belos::PetraPrec class.
//
// Constructor.
//
template <class TYPE>
PetraPrec<TYPE>::PetraPrec(const Epetra_Operator& Prec) : Epetra_Prec(Prec) 
{
}
//
// Destructor.
//
template <class TYPE>
PetraPrec<TYPE>::~PetraPrec() 
{
}
/////////////////////////////////////////////////////////////
//
// Anasazi::Precondition Operator application
//
template <class TYPE>
void PetraPrec<TYPE>::ApplyPrecondition ( const Anasazi::MultiVec<TYPE>& x,
                                                    Anasazi::MultiVec<TYPE>& y) const {
	int info=0;
	Anasazi::MultiVec<TYPE> & temp_x = const_cast<Anasazi::MultiVec<TYPE> &>(x);
     	Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
     	Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);

     	assert( vec_x!=NULL && vec_y!=NULL );
	//
	// Need to cast away constness because the member function Multiply
	// is not declared const.
	//
	info=const_cast<Epetra_Operator&>(Epetra_Prec).ApplyInverse( *vec_x, *vec_y );
	assert(info==0);
}

} // end Belos namespace

// end of file BELOS_PETRA_HPP
#endif 

