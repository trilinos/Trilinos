// file AnasaziMatrix.hpp
#ifndef ANASAZI_MATRIX_HPP
#define ANASAZI_MATRIX_HPP

#include "AnasaziMultiVec.hpp"
#include "AnasaziReturnType.hpp"
#include "AnasaziCommon.hpp"

namespace Anasazi {

/*!	\class Anasazi::Matrix

	\brief Anasazi's templated pure virtual class for constructing the matrix/operator that is
	used by the eigensolver.
	
	A concrete implementation of this class is necessary.  The user can create their own implementation
	if those supplied are not suitable for their needs.

	\author Rich Lehoucq, Heidi Thornquist
*/

template <class TYPE>
class Matrix {
public:

	//@{ \name Constructor/Destructor.
	//! %Anasazi::Matrix constructor.
	Matrix() {
//		cout << "ctor:Anasazi::Matrix " << this << endl; 
	}
	//! %Anasazi::Matrix destructor.
	virtual ~Matrix() {
//		cout << "dtor:Anasazi::Matrix " << this << endl; 
	};
	//@}
	
	//@{ \name Matrix/Operator application method.

	/*! \brief This routine takes the %Anasazi::MultiVec \c x and applies the matrix/operator
	to it resulting in the %Anasazi::MultiVec \c y, which is returned.
	*/
	virtual ReturnType ApplyMatrix (const MultiVec<TYPE>& x, 
						      MultiVec<TYPE>& y ) const = 0;
	//@}
};

} // end Anasazi namespace
#endif
// end of file AnasaziMatrix.hpp
