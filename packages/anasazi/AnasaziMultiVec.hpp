// file AnasaziMultiVec.hpp
#ifndef ANASAZI_MULTI_VEC_HPP
#define ANASAZI_MULTI_VEC_HPP

#include "AnasaziDenseMatrix.hpp"
#include <iostream>

template <class TYPE>
class AnasaziMultiVec {
public:
	AnasaziMultiVec() {
	//		std::cout << "ctor:AnasaziMultiVec " << this << std::endl; 
	};
	virtual ~AnasaziMultiVec () {
	//		std::cout << "dtor:AnasaziMultiVec " << this << std::endl;
	};
	//
	// The following are virtual constructors. Note that
	// until the co-variance of virtual functions is widely
	// supported, the implementation of this must return the
	// same type. 
	// 
	// Clone creates a new MultiVec containing numvecs columns.
	//
	virtual AnasaziMultiVec<TYPE> * Clone ( const int numvecs ) = 0;
	//
	// Deep copy (or copy) constructor.
	//
	virtual AnasaziMultiVec<TYPE> * CloneCopy () = 0;
	//
	// Selective deep copy (or copy) constructor.
	//
	virtual AnasaziMultiVec<TYPE> * CloneCopy ( int [], int ) = 0;
	//
	// CloneView is a view constructor. The MultiVec created shares
	// the columns of *this.
	// 
	virtual AnasaziMultiVec<TYPE> * CloneView ( int [], int ) = 0;
	//
	// Set a block of vectors of the multivec specified by the index
	//
	virtual void SetBlock ( AnasaziMultiVec<TYPE>& A, int index[], int NumVecs ) = 0;
	//
	// Obtain attributes from the multivec
	//
	virtual int GetVecLength () const = 0;
	virtual int GetNumberVecs () const = 0;
	//
	// *this <- alpha * A * B + beta * (*this)
	//
	virtual void MvTimesMatAddMv ( TYPE alpha, AnasaziMultiVec<TYPE>& A, 
		AnasaziDenseMatrix<TYPE>& B, TYPE beta ) = 0;
	//
	// *this <- alpha * A + beta * B
	//
	virtual void MvAddMv ( TYPE alpha, AnasaziMultiVec<TYPE>& A, 
		TYPE beta, AnasaziMultiVec<TYPE>& B ) = 0;
	//
	// B <- alpha * A^T * (*this)
	//
	virtual void MvTransMv ( TYPE alpha, AnasaziMultiVec<TYPE>& A,
							AnasaziDenseMatrix<TYPE>& B) = 0;
	//
	// alpha[i] = norm of i-th column of (*this)
	//
	virtual void MvNorm ( TYPE* ) = 0;
	//
	// random vectors in i-th column of (*this)
	//
	virtual void MvRandom () = 0;
	//
	// print the multi vector
	//
	virtual void MvPrint () = 0;
};
#endif
// end of file AnasaziMultiVec.hpp
