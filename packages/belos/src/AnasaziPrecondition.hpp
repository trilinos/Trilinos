// file AnasaziPrecondition.hpp
#ifndef ANASAZI_PRECONDITION_HPP
#define ANASAZI_PRECONDITION_HPP

#include "AnasaziMultiVec.hpp"
#include <iostream>

template <class TYPE>
class AnasaziPrecondition {
public:
	AnasaziPrecondition() {
//		std::cout << "ctor:AnasaziPrecondition " << this << endl; 
	}
	virtual ~AnasaziPrecondition() {
//		std::cout << "dtor:AnasaziPrecondition " << this << endl; 
	};
	virtual void ApplyPrecondition (const AnasaziMultiVec<TYPE>& x, 
						      AnasaziMultiVec<TYPE>& y ) const = 0;
};

#endif
// end of file AnasaziPrecondition.hpp
