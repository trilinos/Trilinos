// file AnasaziOperator.hpp
#ifndef ANASAZI_OPERATOR_HPP
#define ANASAZI_OPERATOR_HPP

#include "AnasaziMultiVec.hpp"
#include "AnasaziMatrix.hpp"
#include "AnasaziReturnType.hpp"
#include <iostream>

/*!	\class AnasaziOperator

	\brief Anasazi's templated virtual class for constructing the operator that is
	used by the eigensolver.
	
	A concrete implementation of this class is necessary.  The user can create their own implementation
	if those supplied are not suitable for their needs.

	\author Rich Lehoucq, Heidi Thornquist
*/

template <class TYPE>
class AnasaziOperator {
public:
	//@{ \name Constructor/Destructor.
	//! %AnasaziOperator constructor.
	AnasaziOperator() {
//		cout << "ctor:Anasazi::Operator " << this << endl; 
	}

	//! %AnasaziOperator destructor.
	virtual ~AnasaziOperator(void) {
//		cout << "dtor:Anasazi::Operator " << this << endl; 
	};
	//@}

	//@{ \name Matrix/Operator application method.

	/*! \brief This routine takes the %AnasaziMultiVec \c x and applies the matrix/operator
	to it resulting in the %AnasaziMultiVec \c y, which is returned.
	*/
	virtual Anasazi_ReturnType ApplyOp (const AnasaziMultiVec<TYPE>& X, 
						      AnasaziMultiVec<TYPE>& Y ) = 0;
	//@}
};

#endif
// end of file AnasaziOperator.hpp
