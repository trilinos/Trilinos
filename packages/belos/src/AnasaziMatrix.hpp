// file AnasaziMatrix.hpp
#ifndef ANASAZI_MATRIX_HPP
#define ANASAZI_MATRIX_HPP

#include "AnasaziMultiVec.hpp"
#include <iostream>

template <class TYPE>
class AnasaziMatrix {
public:
	AnasaziMatrix() {
//		std::cout << "ctor:AnasaziMatrix " << this << std::endl; 
	}
	virtual ~AnasaziMatrix() {
//		std::cout << "dtor:AnasaziMatrix " << this << std::endl; 
	};
	virtual void ApplyMatrix (const AnasaziMultiVec<TYPE>& x, 
						      AnasaziMultiVec<TYPE>& y ) const = 0;
};

#endif
// end of file AnasaziMatrix.hpp
