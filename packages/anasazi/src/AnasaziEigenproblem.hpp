// file Anasazi_Eigenproblem.hpp
#ifndef ANASAZI_EIGENPROBLEM_H
#define ANASAZI_EIGENPROBLEM_H

#include "AnasaziMultiVec.hpp"
#include "AnasaziDenseMatrix.hpp"
#include "AnasaziOperator.hpp"
#include "AnasaziReturnType.hpp"

template<class TYPE>
class AnasaziEigenproblem {

    public:

	//@{ \name Constructors/Destructor.

	/*! Empty constructor - allows AnasaziEigenproblem to be described
		at a later time through "Set Methods".
	*/
	AnasaziEigenproblem(void);

	//! Standard Eigenvalue Problem Constructor
	AnasaziEigenproblem( AnasaziMatrix<TYPE>* A, AnasaziMultiVec<TYPE>* Ivec );

	//! Standard Eigenvalue Problem Constructor w/ Operator.
	AnasaziEigenproblem( AnasaziMatrix<TYPE>* A, AnasaziOperator<TYPE>* Op,
			AnasaziMultiVec<TYPE>* Ivec );

	//! Generalized Eigenvalue Problem Constructor
	AnasaziEigenproblem( AnasaziOperator<TYPE>* Op , AnasaziMatrix<TYPE>* B,
			AnasaziMultiVec<TYPE>* Ivec );

	//! Generalized Eigenvalue Problem Constructor with Matrix A
	AnasaziEigenproblem( AnasaziMatrix<TYPE>* A, AnasaziMatrix<TYPE>* B,
			AnasaziOperator<TYPE>* Op, AnasaziMultiVec<TYPE>* Ivec );

	//! Eigenvalue Problem Constructor w/ Operator Only.
	AnasaziEigenproblem( AnasaziOperator<TYPE>* Op, AnasaziMultiVec<TYPE>* Ivec );

	//! Copy Constructor.
	AnasaziEigenproblem( const AnasaziEigenproblem<TYPE>* Problem );	

	//! Destructor.
	virtual ~AnasaziEigenproblem(void);
	//@}

	//@{ \name Set Methods.

	/*! Set the initial guess.  This vector is required to create all the space needed 
	by Anasazi to solve the eigenvalue problem.  Even if an initial guess is not known
	by the user, an initial vector must be passed in.  Sets the pointer to the input
	%AnasaziMultiVec, so no copy is made.  
	*/
	void SetInitVec( AnasaziMultiVec<TYPE>* Ivec ) { _InitVec = Ivec; };

	/*! Set the operator for which eigenvalues will be computed.  This may be different
	from the matrix \c A if a spectral transformation is employed, for example.  Sets
	the operator pointer to the input %AnasaziMultiVec, so no copy is made.
	*/
	void SetOperator( AnasaziOperator<TYPE>* Op ) { _Op = Op; };

	//! Set the matrix A (stiffness matrix) of the eigenvalue problem AX = BX\lambda.
	//! Sets the pointer to the input %AnasaziMatrix, so no copy is made.
	void SetMatrixA( AnasaziMatrix<TYPE>* A ) { _A = A; };

	//! Set the matrix B (mass matrix) of the eigenvalue problem AX = BX\lambda.
	//! Sets the pointer to the input %AnasaziMatrix, so no copy is made.
	void SetMatrixB( AnasaziMatrix<TYPE>* B ) { _B = B; }

	//@}

	//@{ \name Accessor Methods.

	//! Get a pointer to the initial vector
	AnasaziMultiVec<TYPE>* GetInitVec() const { return(_InitVec); };

	//! Get a pointer to the Operator.
	AnasaziOperator<TYPE>* GetOperator() const { return(_Op); };

	//! Get a pointer to the stiffness matrix.
	AnasaziMatrix<TYPE>* GetMatrixA() const { return(_A); };

	//! Get a pointer to the mass matrix.
	AnasaziMatrix<TYPE>* GetMatrixB() const { return(_B); };

	//! Get a pointer to the real part of the eigenvalues of the operator.
	TYPE* GetREvals() { return(_REvals); };
	
	//! Get a pointer to the imaginary part of the eigenvalues of the operator.
	TYPE* GetIEvals() { return(_IEvals); };

	//! Get a pointer to the real part of the eigenvectors of the operator.
	AnasaziMultiVec<TYPE>* GetREvecs() { return(_REvecs); };

	//! Get a pointer to the imaginary part of the eigenvectors of the operator.
	AnasaziMultiVec<TYPE>* GetIEvecs() { return(_IEvecs); };

	//@}	

	//@{ \name Matrix/Operator Application Method.

	/*! \brief This routine will apply the matrix/operator for this eigenproblem
	to the %AnasaziMultiVec \c x and the resulting %AnasaziMultiVec is \c y.  If
	the operator is not defined for this eigenproblem, it will apply the matrix A,
	otherwise it will return undefined.
	*/
	Anasazi_ReturnType ApplyOp (const AnasaziMultiVec<TYPE>& X, 
						AnasaziMultiVec<TYPE>& Y);	

	/*! \brief This routine will apply the stiffness matrix (A) to the 
	%AnasaziMultiVec \c X and the resulting %AnasaziMultiVec is \c Y.
	*/
	Anasazi_ReturnType ApplyMatrixA (const AnasaziMultiVec<TYPE>& X, 
						AnasaziMultiVec<TYPE>& Y);

	/*! \brief This routine will apply the mass matrix (B) to the %AnasaziMultiVec
	\c X and the resulting %AnasaziMultiVec is \c Y.
	*/
	Anasazi_ReturnType ApplyMatrixB (const AnasaziMultiVec<TYPE>& X, 
						AnasaziMultiVec<TYPE>& Y);

	//@}
	//@{ \name Inner Product Method.

	/*! \brief Computes A inner product using definition of ApplyMatrixA.  

	\returns  Z = alpha*X^T*A*Y
	*/
	Anasazi_ReturnType AInProd( TYPE alpha, const AnasaziMultiVec<TYPE>& X, 
						const AnasaziMultiVec<TYPE>& Y,
						AnasaziDenseMatrix<TYPE>& Z );

	/*! \brief Computes B inner product using definition of ApplyMatrixB.

	\returns  Z = alpha*X^T*B*Y
	*/
	Anasazi_ReturnType BInProd( TYPE alpha, const AnasaziMultiVec<TYPE>& X, 
						const AnasaziMultiVec<TYPE>& Y,
						AnasaziDenseMatrix<TYPE>& Z );
	//@}

    protected:

	TYPE *_REvals, *_IEvals;
	AnasaziMultiVec<TYPE> *_REvecs, *_IEvecs; 	
	AnasaziMultiVec<TYPE> *_InitVec;
	AnasaziMatrix<TYPE> *_A, *_B; 
	AnasaziOperator<TYPE> *_Op;
};		

//=============================================================================
//	Implementations (Constructors / Destructors)
//=============================================================================

template<class TYPE>
AnasaziEigenproblem<TYPE>::AnasaziEigenproblem(void) : 
	_Op(0), _A(0), _B(0), _InitVec(0),_REvals(0), 
	_IEvals(0), _REvecs(0), _IEvecs(0)
{
}

//=============================================================================

template<class TYPE>
AnasaziEigenproblem<TYPE>::AnasaziEigenproblem( AnasaziMatrix<TYPE>* A, AnasaziMultiVec<TYPE>* Ivec ) :
	_Op(0), _A(A), _B(0), _InitVec(Ivec), _REvals(0), 
	_IEvals(0), _REvecs(0), _IEvecs(0)
{
}

//=============================================================================

template<class TYPE>
AnasaziEigenproblem<TYPE>::AnasaziEigenproblem( AnasaziMatrix<TYPE>* A, AnasaziOperator<TYPE>* Op,
			AnasaziMultiVec<TYPE>* Ivec ) :
	_Op(Op), _A(A), _B(0), _InitVec(Ivec), _REvals(0),
	_IEvals(0), _REvecs(0), _IEvecs(0)
{
}

//=============================================================================

template<class TYPE>
AnasaziEigenproblem<TYPE>::AnasaziEigenproblem( AnasaziOperator<TYPE>* Op , AnasaziMatrix<TYPE>* B,
			AnasaziMultiVec<TYPE>* Ivec ) :
	_Op(Op), _A(0), _B(B), _InitVec(Ivec), _REvals(0),
	_IEvals(0), _REvecs(0), _IEvecs(0)
{
}

//=============================================================================

template<class TYPE>
AnasaziEigenproblem<TYPE>::AnasaziEigenproblem( AnasaziMatrix<TYPE>* A, AnasaziMatrix<TYPE>* B,
			AnasaziOperator<TYPE>* Op, AnasaziMultiVec<TYPE>* Ivec ) :
	_Op(Op), _A(A), _B(B), _InitVec(Ivec), _REvals(0),
	_IEvals(0), _REvecs(0), _IEvecs(0)
{
}

//=============================================================================

template<class TYPE>
AnasaziEigenproblem<TYPE>::AnasaziEigenproblem( AnasaziOperator<TYPE>* Op, AnasaziMultiVec<TYPE>* Ivec ) :
	_Op(Op), _A(0), _B(0), _InitVec(Ivec), _REvals(0),
	_IEvals(0), _REvecs(0), _IEvecs(0)
{
}

//=============================================================================

template<class TYPE>
AnasaziEigenproblem<TYPE>::AnasaziEigenproblem( const AnasaziEigenproblem<TYPE>* Problem ) :
	_Op(Problem._Op), _A(Problem._A), _B(Problem._B),
	_InitVec(Problem._InitVec), _REvals(Problem._REvals),
	_IEvals(Problem._IEvals), _REvecs(Problem._REvecs),
	_IEvecs(Problem._IEvecs)
{
}
	
//=============================================================================

template<class TYPE>
AnasaziEigenproblem<TYPE>::~AnasaziEigenproblem(void)
{
}

//=============================================================================
//	Implementations (Matrix/Operator Application Methods)
//=============================================================================

template<class TYPE>
Anasazi_ReturnType AnasaziEigenproblem<TYPE>::ApplyOp (const AnasaziMultiVec<TYPE>& X, 
							AnasaziMultiVec<TYPE>& Y )
{ 
	if (_Op) {
		return(_Op->ApplyOp( X, Y )); 
	} else if (_A) {
		return(_A->ApplyMatrix( X, Y ));
	} else {
		return Undefined;
	}
}

template<class TYPE>
Anasazi_ReturnType AnasaziEigenproblem<TYPE>::ApplyMatrixB (const AnasaziMultiVec<TYPE>& X, 
							AnasaziMultiVec<TYPE>& Y)
{ 
	if (_B) {
		return(_B->ApplyMatrix( X, Y ));
	} else {
		//
                // First cast away the const on x.
                //
                AnasaziMultiVec<TYPE>& tempX = const_cast<AnasaziMultiVec<TYPE>& >(X);
                //
                // Now create the indexing for copying X into Y.
                //
                int i, numvecs = X.GetNumberVecs();
                int *index = new int[numvecs]; assert(index);
                for (i=0; i<numvecs; i++) { index[i] = i; }
                Y.SetBlock(tempX, index, numvecs);

                delete [] index;
		return Ok;
	}
}

template<class TYPE>
Anasazi_ReturnType AnasaziEigenproblem<TYPE>::ApplyMatrixA (const AnasaziMultiVec<TYPE>& X, 
							AnasaziMultiVec<TYPE>& Y)
{ 
	if (_A) {
		return(_A->ApplyMatrix( X, Y ));
	} else {
		return Undefined;
	}
}

//=============================================================================
//	Implementations (Inner Product Methods)
//=============================================================================

template<class TYPE>
Anasazi_ReturnType AnasaziEigenproblem<TYPE>::AInProd( TYPE alpha, const AnasaziMultiVec<TYPE>& X, 
						const AnasaziMultiVec<TYPE>& Y,
						AnasaziDenseMatrix<TYPE>& Z )
{	
	if (_A) {
		//
                // First cast away the const on Y.
                //
                AnasaziMultiVec<TYPE>& tempX = const_cast<AnasaziMultiVec<TYPE>& >(X);
                AnasaziMultiVec<TYPE>& tempY = const_cast<AnasaziMultiVec<TYPE>& >(Y);
		AnasaziMultiVec<TYPE>* AY = tempY.CloneCopy();

		// Apply A and check that it returned Ok.
		Anasazi_ReturnType ret = _A->ApplyMatrix( Y, *AY );
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
Anasazi_ReturnType AnasaziEigenproblem<TYPE>::BInProd( TYPE alpha, const AnasaziMultiVec<TYPE>& X, 
						const AnasaziMultiVec<TYPE>& Y,
						AnasaziDenseMatrix<TYPE>& Z )
{
	if (_B) {
		//
                // First cast away the const on AnasaziMultiVecs.
                //
                AnasaziMultiVec<TYPE>& tempX = const_cast<AnasaziMultiVec<TYPE>& >(X);
                AnasaziMultiVec<TYPE>& tempY = const_cast<AnasaziMultiVec<TYPE>& >(Y);
		AnasaziMultiVec<TYPE>* BY = tempY.CloneCopy();

		// Apply B and check that it returned Ok.
		Anasazi_ReturnType ret = _B->ApplyMatrix( Y, *BY );
		if ( ret != Ok ) { return ret; }

		// Now perform inner product.  Result is stored in Z.
		BY->MvTransMv( alpha, tempX, Z );

		delete BY;
		return Ok;				
	} else {
		// Perform the inner product, assume B=I.
                AnasaziMultiVec<TYPE>& tempY = const_cast<AnasaziMultiVec<TYPE>& >(Y);
                AnasaziMultiVec<TYPE>& tempX = const_cast<AnasaziMultiVec<TYPE>& >(X);
		tempY.MvTransMv ( alpha, tempX, Z );
		return Ok;				
	}
}

#endif

// end AnasaziEigenproblem.hpp
