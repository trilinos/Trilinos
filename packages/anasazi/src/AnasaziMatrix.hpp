// file AnasaziMatrix.hpp
#ifndef ANASAZI_MATRIX_HPP
#define ANASAZI_MATRIX_HPP

#include "AnasaziMultiVec.hpp"
#include <iostream>

/*!	\class AnasaziMatrix

	\brief Anasazi's templated pure virtual class for constructing the matrix/operator that is
	used by the eigensolver.
	
	A concrete implementation of this class is necessary.  The user can create their own implementation
	if those supplied are not suitable for their needs.

	\author Rich Lehoucq, Heidi Thornquist
*/

template <class TYPE>
class AnasaziMatrix {
public:
	//@{ \name Constructor/Destructor.
	//! %AnasaziMatrix constructor.
	AnasaziMatrix() {
//		cout << "ctor:AnasaziMatrix " << this << endl; 
	}
	//! %AnasaziMatrix destructor.
	virtual ~AnasaziMatrix() {
//		cout << "dtor:AnasaziMatrix " << this << endl; 
	};
	//@}
	
	//@{ \name Matrix/Operator application method.

	/*! \brief This routine takes the %AnasaziMultiVec \c x and applies the matrix/operator
	to it resulting in the %AnasaziMultiVec \c y, which is returned.
	*/
	virtual void ApplyMatrix (const AnasaziMultiVec<TYPE>& x, 
						      AnasaziMultiVec<TYPE>& y ) const = 0;
	//@}
};

#endif
// end of file AnasaziMatrix.hpp
