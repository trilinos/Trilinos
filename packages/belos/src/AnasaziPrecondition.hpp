// file AnasaziPrecondition.hpp
#ifndef ANASAZI_PRECONDITION_HPP
#define ANASAZI_PRECONDITION_HPP

#include "AnasaziMultiVec.hpp"
#include "BelosConfigDefs.hpp"

/*!	\class AnasaziPrecondition

	\brief Belos's templated virtual class for constructing the 
	preconditioner to the AnasaziMatrix that is used by the linear 
	solver.

	A concrete implementation of this class is not necessary.  The user 
	can create their own implementation if those supplied are not
	suitable for their needs.  By default, this class will provide
	the identity as the preconditioner.  Thus, if the ApplyPrecondition
	function is not overridden, then the block linear solver will be
	unpreconditioned.

	\author Rich Lehoucq, Teri Barth
*/

template <class TYPE>
class AnasaziPrecondition {
public:
	//@{ \name Constructor/Destructor.
	//! %AnasaziPrecondition constructor.
	AnasaziPrecondition() {
//		cout << "ctor:AnasaziPrecondition " << this << endl; 
	}
	//! %AnasaziPrecondition destructor.
	virtual ~AnasaziPrecondition() {
//		cout << "dtor:AnasaziPrecondition " << this << endl; 
	};
	//@}

	//@{ \name Preconditioner application method.
	
	/*! \brief This routine takes the %AnasaziMultiVec \c x and 
	applies the preconditioner to it resulting in the %AnasaziMultiVec \c y, 
	which is returned.  If this routine is not overridden, then the %AnasaziMultiVec
	\c x will be passed directly to \c y, resulting in an unpreconditioned linear
	solver.
	*/
	virtual void ApplyPrecondition (const AnasaziMultiVec<TYPE>& x, 
						      AnasaziMultiVec<TYPE>& y ) const {
		//
		// First cast away the const on x.
		//
		AnasaziMultiVec<TYPE>& temp_x = const_cast<AnasaziMultiVec<TYPE>& >(x);
		//
		// Now create the indexing for copying x into y.
		//
		int i, numvecs = x.GetNumberVecs();
		int *index = new int[numvecs]; assert(index);
		for (i=0; i<numvecs; i++) { index[i] = i; }
		y.SetBlock(temp_x, index, numvecs);
		delete [] index;
	}
	//@}
};

#endif
// end of file AnasaziPrecondition.hpp
