// file AnasaziOperator.hpp
#ifndef ANASAZI_OPERATOR_HPP
#define ANASAZI_OPERATOR_HPP

#include "AnasaziMultiVec.hpp"
#include "AnasaziReturnType.hpp"
#include "AnasaziConfigDefs.hpp"

namespace Anasazi {

/*!	\class Anasazi::Operator

	\brief Anasazi's templated virtual class for constructing the operator that is
	used by the eigensolver.
	
	A concrete implementation of this class is necessary.  The user can create their own implementation
	if those supplied are not suitable for their needs.

	\author Rich Lehoucq, Heidi Thornquist
*/

template <class TYPE>
class Operator {
public:
	//@{ \name Constructor/Destructor.
	//! %Anasazi::Operator constructor.
	Operator() {
//		cout << "ctor:Anasazi::Operator " << this << endl; 
	}

	//! %Anasazi::Operator destructor.
	virtual ~Operator(void) {
//		cout << "dtor:Anasazi::Operator " << this << endl; 
	};
	//@}

	//@{ \name Matrix/Operator application method.

	/*! \brief This routine takes the %Anasazi::MultiVec \c x and applies the matrix/operator
	to it resulting in the %Anasazi::MultiVec \c y, which is returned.
	*/
	virtual ReturnType ApplyOp (const MultiVec<TYPE>& x, MultiVec<TYPE>& y ) const = 0;
	//@}
};

} // end of Anasazi namespace
#endif
// end of file AnasaziOperator.hpp
