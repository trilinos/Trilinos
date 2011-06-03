// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
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

#ifndef TEUCHOS_2DARRAY_HPP
#define TEUCHOS_2DARRAY_HPP


/*! \file Teuchos_2DArray.hpp
  \brief A thin wrapper around the Teuchos Array class that allows for 
  2 dimensional arrays.
*/

#include "Teuchos_Array.hpp"


namespace Teuchos{


template<class T> 
class 2DArray{
public:
  2DArray(size_type width, size_type height):_width(width),_height(height),_data(Array(width*height, 0)){}
  2DArray(size_type width, size_type height, Array<T> data):
    _width(width),_height(height),_data(data){}
  virtual ~2DArray(){}


  inline ArrayView<T> 2DArray::operator[](size_type i);
  
  inline ArrayView<const T> 2DArray::operator[](size_type i) const;

  inline size_type getWidth() const{
    return _width;
  }

  inline size_type getHeight() const{
    return _height;
  }

  inline Array<const T> get1DArray() const{
    return _data;
  }

private:
  size_type _width, _height;
  Array<T> _data;
};

inline Array<const T> 2DArray::getDataArray() const{
  return _data;
}

template<class T> inline
ArrayView<T> 2DArray::operator[](size_type i){
  return data.view(width*i, width);  
}

template<class T> inline
ArrayView<const T> 2DArray::operator[](size_type i) const{
  return data.view(width*i, width);  
}

template<class T>
std::istringstream& operator>> (std::istringstream& in, 2DArray<T>& array){
  std::string input = in.str();
  size_type dimCharPos = input.find_first_of('x');
  size_type delimCharPos = input.find_first_of(':');
  std::string valueData = input.substr(delimCharPos+1);
  istringstream widthStream(input.substr(0,dimCharPos))
  istringstream heightStream(input.substr(dimCharPos+1, delimCharPos-dimCharPos-1));
  size_t width, height;
  widthStream >> width;
  heightStream >> height;
  array = fromStringToArray<T>(valueData);
  array = 2DArray(width, height, array);
  return in;
}

template<class T> inline
std::ostream& Teuchos::operator<<( std::ostream& os, const 2DArray<T>& array){
  return os << array.getWidth() << "x" << array.getHeight() << 
    ":" << Teuchos::toString(array.get1DArray());
}


} //namespace Teuchos


#endif // TEUCHOS_ARRAY_H
