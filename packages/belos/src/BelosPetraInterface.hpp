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
#include "BelosBlockGmres.hpp"
#include "BelosBlockCG.hpp"

//--------template class BelosPetraVec-------------------------------------
template <class TYPE>
class BelosPetraVec : public AnasaziMultiVec<TYPE>, public Epetra_MultiVector {
public:
// constructors
	BelosPetraVec(const Epetra_BlockMap&, TYPE *, const int, const int stride=0);
	BelosPetraVec(const Epetra_BlockMap&, const int);
	BelosPetraVec(Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, int index[], int NumVecs );
	BelosPetraVec(const Epetra_MultiVector & P_vec);
	~BelosPetraVec();
	//
	//  member functions inherited from AnasaziMultiVec
	//
	//  the following is a virtual copy constructor returning
	//  a pointer to the pure virtual class. vector values are
	//  not copied; instead a new MultiVec is created containing
	//  a non-zero amount of columns.
	//
	AnasaziMultiVec<TYPE> * Clone ( const int numvecs );
	//
	//  the following is a virtual copy constructor returning
	//  a pointer to the pure virtual class. vector values are
	//  copied and a new stand-alone MultiVector is created.
	//  (deep copy).
	//
	AnasaziMultiVec<TYPE> * CloneCopy ();
	//
	//  the following is a virtual copy constructor returning
	//  a pointer to the pure virtual class. vector values are
	//  copied and a new stand-alone MultiVector is created
	//  where only selected columns are chosen.  (deep copy).
	//
	AnasaziMultiVec<TYPE> * CloneCopy ( int index[], int numvecs );
	//
	//  the following is a virtual view constructor returning
	//  a pointer to the pure virtual class. vector values are 
	//  shared and hence no memory is allocated for the columns.
	//
	AnasaziMultiVec<TYPE> * CloneView ( int index[], int numvecs );
	//
	//  this routine sets a subblock of the multivector, which
	//  need not be contiguous, and is given by the indices.
	//
	void SetBlock ( AnasaziMultiVec<TYPE>& A, int index[], int numvecs );
	//
	int GetNumberVecs () const;
	int GetVecLength () const;
	//
	// *this <- alpha * A * B + beta * (*this)
	//
	void MvTimesMatAddMv ( TYPE alpha, AnasaziMultiVec<TYPE>& A, 
		AnasaziDenseMatrix<TYPE>& B, TYPE beta );
	//
	// *this <- alpha * A + beta * B
	//
	void MvAddMv ( TYPE alpha , AnasaziMultiVec<TYPE>& A, TYPE beta,
		AnasaziMultiVec<TYPE>& B);
	//
	// B <- alpha * A^T * (*this)
	//
	void MvTransMv ( TYPE alpha, AnasaziMultiVec<TYPE>& A, AnasaziDenseMatrix<TYPE>& B );
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
BelosPetraVec<TYPE>::BelosPetraVec(const Epetra_BlockMap& Map, TYPE * array, 
					   	const int numvec, const int stride)
	: Epetra_MultiVector(Copy, Map, array, stride, numvec) 
{
}

template<class TYPE>
BelosPetraVec<TYPE>::BelosPetraVec(const Epetra_BlockMap& Map, const int numvec)
	: Epetra_MultiVector(Map, numvec) 
{
}

template<class TYPE>
BelosPetraVec<TYPE>::BelosPetraVec(Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, 
						int index[], int NumVecs ) 
	: Epetra_MultiVector(CV, P_vec, index, NumVecs) 
{
}

template<class TYPE>
BelosPetraVec<TYPE>::BelosPetraVec(const Epetra_MultiVector& P_vec) 
	: Epetra_MultiVector(P_vec) 
{
}

template<class TYPE>
BelosPetraVec<TYPE>::~BelosPetraVec() 
{
}
//
//  member functions inherited from AnasaziMultiVec
//
//
//  Simulating a virtual copy constructor. If we could rely on the co-variance
//  of virtual functions, we could return a pointer to BelosPetraVec<TYPE>
//  (the derived type) instead of a pointer to the pure virtual base class.
//
template<class TYPE>
AnasaziMultiVec<TYPE>* BelosPetraVec<TYPE>::Clone ( const int NumVecs ) {
	BelosPetraVec * ptr_apv = new BelosPetraVec(Map(),NumVecs);
	return ptr_apv; // safe upcast.
}
	//
	//  the following is a virtual copy constructor returning
	//  a pointer to the pure virtual class. vector values are
	//  copied.
	//
template<class TYPE>
AnasaziMultiVec<TYPE>* BelosPetraVec<TYPE>::CloneCopy() {
	BelosPetraVec *ptr_apv = new BelosPetraVec(*this);
	return ptr_apv; // safe upcast
}

template<class TYPE>
AnasaziMultiVec<TYPE>* BelosPetraVec<TYPE>::CloneCopy ( int index[], int numvecs ) {
	BelosPetraVec * ptr_apv = new BelosPetraVec(Copy, *this, index, numvecs );
	return ptr_apv; // safe upcast.
}

template<class TYPE>
AnasaziMultiVec<TYPE>* BelosPetraVec<TYPE>::CloneView ( int index[], int numvecs ) {
	BelosPetraVec * ptr_apv = new BelosPetraVec(View, *this, index, numvecs );
	return ptr_apv; // safe upcast.
}

template<class TYPE>
void BelosPetraVec<TYPE>::SetBlock(AnasaziMultiVec<TYPE>& A, int index[], int numvecs ) 
{	
	int i,j,ind;
	BelosPetraVec *A_vec = dynamic_cast<BelosPetraVec *>(&A); assert(A_vec);
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
int BelosPetraVec<TYPE>::GetNumberVecs () const {
	return NumVectors();
}

template<class TYPE>
int BelosPetraVec<TYPE>::GetVecLength () const {
	return MyLength();
}
//
// *this <- alpha * A * B + beta * (*this)
//
template<class TYPE>
void BelosPetraVec<TYPE>::MvTimesMatAddMv ( TYPE alpha, AnasaziMultiVec<TYPE>& A, 
								   AnasaziDenseMatrix<TYPE>& B, TYPE beta ) {
	int info=0;
	const int izero=0;
	char* trans="N";
	Epetra_LocalMap LocalMap(B.getrows(), izero, Map().Comm());
	Epetra_MultiVector B_Pvec(Copy, LocalMap, B.getarray(), B.getld(), B.getcols());

	BelosPetraVec *A_vec = dynamic_cast<BelosPetraVec *>(&A); assert(A_vec);

	info = Multiply( *trans, *trans, alpha, *A_vec, B_Pvec, beta );

	assert(info==0);
}
//
// *this <- alpha * A + beta * B
//
template<class TYPE>
void BelosPetraVec<TYPE>::MvAddMv ( TYPE alpha , AnasaziMultiVec<TYPE>& A, 
						   TYPE beta, AnasaziMultiVec<TYPE>& B) {
	int info=0;
	const TYPE one =1.0;
	const TYPE zero = 0.0;

	BelosPetraVec *A_vec = dynamic_cast<BelosPetraVec *>(&A); assert(A_vec);
	BelosPetraVec *B_vec = dynamic_cast<BelosPetraVec *>(&B); assert(B_vec);

	info = Update( alpha, *A_vec, beta, *B_vec, zero ); assert(info==0);
}
//
// dense B <- alpha * A^T * (*this)
//
template<class TYPE>
void BelosPetraVec<TYPE>::MvTransMv ( TYPE alpha, AnasaziMultiVec<TYPE>& A,
									   AnasaziDenseMatrix<TYPE>& B) {
	int info=0;
	const int izero=0;
	const TYPE zero=0.0;
	//const TYPE one=1.0;
	char* trans1="T";
	char* trans2="N";

	BelosPetraVec *A_vec = dynamic_cast<BelosPetraVec *>(&A);

	if (A_vec) {

		Epetra_LocalMap LocalMap(B.getrows(), izero, Map().Comm());
		Epetra_MultiVector B_Pvec(View, LocalMap, B.getarray(), B.getld(), B.getcols());

		info = B_Pvec.Multiply( *trans1, *trans2, alpha, *A_vec, *this, zero ); 
		assert(info==0);
	}
}
//
// alpha[i] = norm of i-th column of (*this)
//
template<class TYPE>
void BelosPetraVec<TYPE>::MvNorm ( TYPE * normvec ) {
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
void BelosPetraVec<TYPE>::MvRandom () {
	int info=0;
	info = Random();
	assert(info==0);
}
//
// initializes each element of (*this) with alpha
//

template<class TYPE>
void BelosPetraVec<TYPE>::MvInit( TYPE alpha )
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
void BelosPetraVec<TYPE>::MvPrint() {
	cout << *this << endl;
}

///////////////////////////////////////////////////////////////
//--------template class BelosPetraMat-----------------------
template <class TYPE> 
class BelosPetraMat : public AnasaziMatrix<TYPE> {
public:
	BelosPetraMat(const Epetra_Operator& Matrix);
	~BelosPetraMat();
	Anasazi_ReturnType ApplyMatrix ( const AnasaziMultiVec<TYPE>& x, AnasaziMultiVec<TYPE>& y ) const;
	const Epetra_Operator & GetMatrix () { return(Epetra_Mat); };
private:
	const Epetra_Operator & Epetra_Mat;
};
//-------------------------------------------------------------
//
// implementation of the BelosPetraMat class.
//
////////////////////////////////////////////////////////////////////
//
// AnasaziMatrix constructors
//
template <class TYPE>
BelosPetraMat<TYPE>::BelosPetraMat(const Epetra_Operator& Matrix) 
	: Epetra_Mat(Matrix) 
{
}

template <class TYPE>
BelosPetraMat<TYPE>::~BelosPetraMat() 
{
}
//
// AnasaziMatrix matrix multiply
//
template <class TYPE>
Anasazi_ReturnType BelosPetraMat<TYPE>::ApplyMatrix ( const AnasaziMultiVec<TYPE>& x, 
					  AnasaziMultiVec<TYPE>& y ) const {
	int info=0;
	AnasaziMultiVec<TYPE> & temp_x = const_cast<AnasaziMultiVec<TYPE> &>(x);
	Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
	Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);

	assert( vec_x && vec_y );
	//
	// Need to cast away constness because the member function Multiply
	// is not declared const.
	//
	info=const_cast<Epetra_Operator&>(Epetra_Mat).Apply( *vec_x, *vec_y );
	assert(info==0);
	return Ok;	
}
 
///////////////////////////////////////////////////////////////
//--------template class BelosPetraPrec--------------------

template <class TYPE>
class BelosPetraPrec : public AnasaziPrecondition<TYPE> {
public:
        BelosPetraPrec(const Epetra_Operator& Prec);
        ~BelosPetraPrec();
        void ApplyPrecondition ( const AnasaziMultiVec<TYPE>& x, AnasaziMultiVec<TYPE>& y ) const;
private:
   	const Epetra_Operator& Epetra_Prec;
};
//--------------------------------------------------------------
//
// implementation of the BelosPetraPrec class.
//
// Constructor.
//
template <class TYPE>
BelosPetraPrec<TYPE>::BelosPetraPrec(const Epetra_Operator& Prec) : Epetra_Prec(Prec) 
{
}
//
// Destructor.
//
template <class TYPE>
BelosPetraPrec<TYPE>::~BelosPetraPrec() 
{
}
/////////////////////////////////////////////////////////////
//
// AnasaziPrecondition Operator application
//
template <class TYPE>
void BelosPetraPrec<TYPE>::ApplyPrecondition ( const AnasaziMultiVec<TYPE>& x,
                                                    AnasaziMultiVec<TYPE>& y) const {
	int info=0;
	AnasaziMultiVec<TYPE> & temp_x = const_cast<AnasaziMultiVec<TYPE> &>(x);
     	Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
     	Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);

     	assert( vec_x && vec_y );
	//
	// Need to cast away constness because the member function Multiply
	// is not declared const.
	//
	info=const_cast<Epetra_Operator&>(Epetra_Prec).ApplyInverse( *vec_x, *vec_y );
	assert(info==0);
}

// end of file BELOS_PETRA_HPP
#endif 

