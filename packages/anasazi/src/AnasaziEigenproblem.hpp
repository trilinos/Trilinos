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

#ifndef ANASAZI_EIGENPROBLEM_H
#define ANASAZI_EIGENPROBLEM_H

#include "AnasaziMatrix.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziReturnType.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

/*! \class Anasazi::Eigenproblem
    \brief This class defines the standard or generalized eigenvalue problem
	to be solved.

	The type of eigenproblem is determined by the constructor used.  This
	class is passed to the solver and is responsible for providing the solver
	with the essential operations.
*/
 	

namespace Anasazi {

template<class TYPE>
class Eigenproblem {

    public:

	//@{ \name Constructors/Destructor.

	//! Empty constructor - allows Anasazi::Eigenproblem to be described at a later time through "Set Methods".
	Eigenproblem(void);

	//! Standard Eigenvalue Problem Constructor
	Eigenproblem( Matrix<TYPE>* A, MultiVec<TYPE>* Ivec );

	//! Standard Eigenvalue Problem Constructor w/ Operator.
	Eigenproblem( Matrix<TYPE>* A, Operator<TYPE>* Op,
			MultiVec<TYPE>* Ivec );

	//! Generalized Eigenvalue Problem Constructor
	Eigenproblem( Operator<TYPE>* Op , Matrix<TYPE>* B,
			MultiVec<TYPE>* Ivec );

	//! Generalized Eigenvalue Problem Constructor with Matrix A
	Eigenproblem( Matrix<TYPE>* A, Matrix<TYPE>* B,
			Operator<TYPE>* Op, MultiVec<TYPE>* Ivec );

	//! Eigenvalue Problem Constructor w/ Operator Only.
	Eigenproblem( Operator<TYPE>* Op, MultiVec<TYPE>* Ivec );

	//! Copy Constructor.
	Eigenproblem( const Eigenproblem<TYPE>* Problem );	

	//! Destructor.
	virtual ~Eigenproblem(void);
	//@}

	//@{ \name Set Methods.

	/*! \brief Set the initial guess.  This vector is required to create all the space needed 
	by Anasazi to solve the eigenvalue problem.  Even if an initial guess is not known
	by the user, an initial vector must be passed in.  Sets the pointer to the input
	%Anasazi::MultiVec, so no copy is made.  
	*/
	void SetInitVec( MultiVec<TYPE>* Ivec ) { _InitVec = Ivec; };

	/*! \brief Set the operator for which eigenvalues will be computed.  This may be different
	from the matrix \c A if a spectral transformation is employed, for example.  Sets
	the operator pointer to the input %Anasazi::MultiVec, so no copy is made.
	*/
	void SetOperator( Operator<TYPE>* Op ) { _Op = Op; };

	/*! \brief Set the matrix A (stiffness matrix) of the eigenvalue problem AX = BX\lambda.
	Sets the pointer to the input %Anasazi::Matrix, so no copy is made.
	*/
	void SetMatrixA( Matrix<TYPE>* A ) { _Amat = A; };

	/*! \brief Set the matrix B (mass matrix) of the eigenvalue problem AX = BX\lambda.
	Sets the pointer to the input %Anasazi::Matrix, so no copy is made.
	*/
	void SetMatrixB( Matrix<TYPE>* B ) { _Bmat = B; }

	//@}

	//@{ \name Accessor Methods.

	//! Get a pointer to the initial vector
	MultiVec<TYPE>* GetInitVec() const { return(_InitVec); };

	//! Get a pointer to the Operator.
	Operator<TYPE>* GetOperator() const { return(_Op); };

	//! Get a pointer to the stiffness matrix.
	Matrix<TYPE>* GetMatrixA() const { return(_Amat); };

	//! Get a pointer to the mass matrix.
	Matrix<TYPE>* GetMatrixB() const { return(_Bmat); };

	//! Get a pointer to the real part of the eigenvalues of the operator.
	TYPE* GetREvals() { return(_REvals); };
	
	//! Get a pointer to the imaginary part of the eigenvalues of the operator.
	TYPE* GetIEvals() { return(_IEvals); };

	//! Get a pointer to the real part of the eigenvectors of the operator.
	MultiVec<TYPE>* GetREvecs() { return(_REvecs); };

	//! Get a pointer to the imaginary part of the eigenvectors of the operator.
	MultiVec<TYPE>* GetIEvecs() { return(_IEvecs); };

	//@}	

	//@{ \name Matrix/Operator Application Method.

	/*! \brief This routine will apply the matrix/operator for this eigenproblem
	to the %Anasazi::MultiVec \c x and the resulting %Anasazi::MultiVec is \c y.  If
	the operator is not defined for this eigenproblem, it will apply the matrix A,
	otherwise it will return undefined.
	*/
	ReturnType ApplyOp (const MultiVec<TYPE>& X, 
						MultiVec<TYPE>& Y);	

	/*! \brief This routine will apply the stiffness matrix (A) to the 
	%Anasazi::MultiVec \c X and the resulting %Anasazi::MultiVec is \c Y.
	*/
	ReturnType ApplyMatrixA (const MultiVec<TYPE>& X, 
						MultiVec<TYPE>& Y);

	/*! \brief This routine will apply the mass matrix (B) to the %Anasazi::MultiVec
	\c X and the resulting %Anasazi::MultiVec is \c Y.
	*/
	ReturnType ApplyMatrixB (const MultiVec<TYPE>& X, 
						MultiVec<TYPE>& Y);

	//@}
	//@{ \name Inner Product Methods.

	/*! \brief Computes A inner product using definition of ApplyMatrixA.  

	\returns  Z = alpha*X^T*A*Y
	*/
	ReturnType AInProd( TYPE alpha, const MultiVec<TYPE>& X, 
						const MultiVec<TYPE>& Y,
						Teuchos::SerialDenseMatrix<int,TYPE>& Z );

	/*! \brief Computes B inner product using definition of ApplyMatrixB.

	\returns  Z = alpha*X^T*B*Y
	*/
	ReturnType BInProd( TYPE alpha, const MultiVec<TYPE>& X, 
						const MultiVec<TYPE>& Y,
						Teuchos::SerialDenseMatrix<int,TYPE>& Z );
	//@}

	//@{ \name Norm Methods.
	/*! \brief Computes the B norm of MultiVecs.
	*/
	ReturnType BMvNorm( MultiVec<TYPE>& X, TYPE* normvec );

	//@}	

    protected:

        Operator<TYPE> *_Op;
        Matrix<TYPE> *_Amat, *_Bmat; 
        MultiVec<TYPE> *_InitVec;
	TYPE *_REvals, *_IEvals;
	MultiVec<TYPE> *_REvecs, *_IEvecs; 	
};		

//=============================================================================
//	Implementations (Constructors / Destructors)
//=============================================================================

template<class TYPE>
Eigenproblem<TYPE>::Eigenproblem(void) : 
	_Op(0), _Amat(0), _Bmat(0), _InitVec(0),_REvals(0), 
	_IEvals(0), _REvecs(0), _IEvecs(0)
{
}

//=============================================================================

template<class TYPE>
Eigenproblem<TYPE>::Eigenproblem( Matrix<TYPE>* A, MultiVec<TYPE>* Ivec ) :
	_Op(0), _Amat(A), _Bmat(0), _InitVec(Ivec), _REvals(0), 
	_IEvals(0), _REvecs(0), _IEvecs(0)
{
}

//=============================================================================

template<class TYPE>
Eigenproblem<TYPE>::Eigenproblem( Matrix<TYPE>* A, Operator<TYPE>* Op,
			MultiVec<TYPE>* Ivec ) :
	_Op(Op), _Amat(A), _Bmat(0), _InitVec(Ivec), _REvals(0),
	_IEvals(0), _REvecs(0), _IEvecs(0)
{
}

//=============================================================================

template<class TYPE>
Eigenproblem<TYPE>::Eigenproblem( Operator<TYPE>* Op , Matrix<TYPE>* B,
			MultiVec<TYPE>* Ivec ) :
	_Op(Op), _Amat(0), _Bmat(B), _InitVec(Ivec), _REvals(0),
	_IEvals(0), _REvecs(0), _IEvecs(0)
{
}

//=============================================================================

template<class TYPE>
Eigenproblem<TYPE>::Eigenproblem( Matrix<TYPE>* A, Matrix<TYPE>* B,
			Operator<TYPE>* Op, MultiVec<TYPE>* Ivec ) :
	_Op(Op), _Amat(A), _Bmat(B), _InitVec(Ivec), _REvals(0),
	_IEvals(0), _REvecs(0), _IEvecs(0)
{
}

//=============================================================================

template<class TYPE>
Eigenproblem<TYPE>::Eigenproblem( Operator<TYPE>* Op, MultiVec<TYPE>* Ivec ) :
	_Op(Op), _Amat(0), _Bmat(0), _InitVec(Ivec), _REvals(0),
	_IEvals(0), _REvecs(0), _IEvecs(0)
{
}

//=============================================================================

template<class TYPE>
Eigenproblem<TYPE>::Eigenproblem( const Eigenproblem<TYPE>* Problem ) :
	_Op(Problem._Op), _Amat(Problem._Amat), _Bmat(Problem._Bmat),
	_InitVec(Problem._InitVec), _REvals(Problem._REvals),
	_IEvals(Problem._IEvals), _REvecs(Problem._REvecs),
	_IEvecs(Problem._IEvecs)
{
}
	
//=============================================================================

template<class TYPE>
Eigenproblem<TYPE>::~Eigenproblem(void)
{
}

//=============================================================================
//	Implementations (Matrix/Operator Application Methods)
//=============================================================================

template<class TYPE>
ReturnType Eigenproblem<TYPE>::ApplyOp (const MultiVec<TYPE>& X, 
							MultiVec<TYPE>& Y )
{ 
	if (_Op) {
		return(_Op->ApplyOp( X, Y )); 
	} else if (_Amat) {
		return(_Amat->ApplyMatrix( X, Y ));
	} else {
		return Undefined;
	}
}

template<class TYPE>
ReturnType Eigenproblem<TYPE>::ApplyMatrixB (const MultiVec<TYPE>& X, 
							MultiVec<TYPE>& Y)
{ 
	if (_Bmat) {
		return(_Bmat->ApplyMatrix( X, Y ));
	} else {
		//
                // First cast away the const on x.
                //
                MultiVec<TYPE>& tempX = const_cast<MultiVec<TYPE>& >(X);
                //
                // Now create the indexing for copying X into Y.
                //
                int i, numvecs = X.GetNumberVecs();
                int *index = new int[numvecs]; assert(index!=NULL);
                for (i=0; i<numvecs; i++) { index[i] = i; }
                Y.SetBlock(tempX, index, numvecs);

                delete [] index;
		return Ok;
	}
}

template<class TYPE>
ReturnType Eigenproblem<TYPE>::ApplyMatrixA (const MultiVec<TYPE>& X, 
							MultiVec<TYPE>& Y)
{ 
	if (_Amat) {
		return(_Amat->ApplyMatrix( X, Y ));
	} else {
		return Undefined;
	}
}

//=============================================================================
//	Implementations (Inner Product Methods)
//=============================================================================

template<class TYPE>
ReturnType Eigenproblem<TYPE>::AInProd( TYPE alpha, const MultiVec<TYPE>& X, 
						const MultiVec<TYPE>& Y,
						Teuchos::SerialDenseMatrix<int,TYPE>& Z )
{	
	if (_Amat) {
		//
                // First cast away the const on Y.
                //
                MultiVec<TYPE>& tempX = const_cast<MultiVec<TYPE>& >(X);
                MultiVec<TYPE>& tempY = const_cast<MultiVec<TYPE>& >(Y);
		MultiVec<TYPE>* AY = tempY.CloneCopy();

		// Apply A and check that it returned Ok.
		ReturnType ret = _Amat->ApplyMatrix( Y, *AY );
		if ( ret != Ok ) { return ret; }

		// Now perform inner product.  Result is stored in Z.
		AY->MvTransMv ( alpha, tempX, Z );

		delete AY;
		return Ok;		
	} else {
		return Undefined;
	}
}

template<class TYPE>
ReturnType Eigenproblem<TYPE>::BInProd( TYPE alpha, const MultiVec<TYPE>& X, 
						const MultiVec<TYPE>& Y,
						Teuchos::SerialDenseMatrix<int,TYPE>& Z )
{
	if (_Bmat) {
		//
                // First cast away the const on MultiVecs.
                //
                MultiVec<TYPE>& tempX = const_cast<MultiVec<TYPE>& >(X);
                MultiVec<TYPE>& tempY = const_cast<MultiVec<TYPE>& >(Y);
		MultiVec<TYPE>* BY = tempY.CloneCopy();

		// Apply B and check that it returned Ok.
		ReturnType ret = _Bmat->ApplyMatrix( Y, *BY );
		if ( ret != Ok ) { return ret; }

		// Now perform inner product.  Result is stored in Z.
		BY->MvTransMv( alpha, tempX, Z );

		delete BY;
		return Ok;				
	} else {
		// Perform the inner product, assume B=I.
                MultiVec<TYPE>& tempY = const_cast<MultiVec<TYPE>& >(Y);
                MultiVec<TYPE>& tempX = const_cast<MultiVec<TYPE>& >(X);
		tempY.MvTransMv ( alpha, tempX, Z );
		return Ok;				
	}
}

//=============================================================================
//	Implementations (Norm Methods)
//=============================================================================

template<class TYPE>
ReturnType Eigenproblem<TYPE>::BMvNorm( MultiVec<TYPE>& X, TYPE* normvec )
{
	int IntOne = 1;
	int numvecs = X.GetNumberVecs();
	Teuchos::SerialDenseVector<int,TYPE> DenseOne(IntOne);
	MultiVec<TYPE>* Xj = 0;
	int *index = new int[IntOne];	
	
	for (int i=0; i<numvecs; i++) {
		index[0] = i;
		Xj = X.CloneView( index, IntOne );
		BInProd( 1.0, *Xj, *Xj, DenseOne );
		normvec[i] = sqrt(DenseOne(0));
	}

	delete Xj;
	delete [] index;

	return Ok;
}

} // end Anasazi namespace
#endif

// end AnasaziEigenproblem.hpp
