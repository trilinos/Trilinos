// file AnasaziPrecondition.hpp
#ifndef ANASAZI_PRECONDITION_HPP
#define ANASAZI_PRECONDITION_HPP

#include "AnasaziMultiVec.hpp"
#include "BelosConfigDefs.hpp"

/*!	\class AnasaziPrecondition

	\brief Belos's templated pure virtual class for constructing the 
	preconditioner to the AnasaziMatrix that is used by the linear 
	solver.

	A concrete implementation of this class is necessary.  The user 
	can create their own implementation if those supplied are not
	suitable for their needs.

	\author Rich Lehoucq, Teri Barth
*/

template <class TYPE>
class AnasaziPrecondition {
public:
	//@{ \name Constructor/Destructor.
	//! %AnasaziPrecondition constructor.
	AnasaziPrecondition() {
//		std::cout << "ctor:AnasaziPrecondition " << this << endl; 
	}
	//! %AnasaziPrecondition destructor.
	virtual ~AnasaziPrecondition() {
//		std::cout << "dtor:AnasaziPrecondition " << this << endl; 
	};
	//@}

	//@{ \name Preconditioner application method.
	
	/*! \brief This routine takes the %AnasaziMultiVec \c x and 
	applies the preconditioner to it resulting in the %AnasaziMultiVec \c y, 
	which is returned.
	*/
	virtual void ApplyPrecondition (const AnasaziMultiVec<TYPE>& x, 
						      AnasaziMultiVec<TYPE>& y ) const = 0;
	//@}
};

#endif
// end of file AnasaziPrecondition.hpp
