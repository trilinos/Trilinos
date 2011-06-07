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

#ifndef TEUCHOS_TWODARRAY_HPP
#define TEUCHOS_TWODARRAY_HPP


/*! \file Teuchos_TwoDArray.hpp
  \brief A thin wrapper around the Teuchos Array class that allows for 
  2 dimensional arrays.
*/

#include "Teuchos_Array.hpp"


namespace Teuchos{


template<class T> 
class TwoDArray{
public:
  typedef Ordinal size_type;
  TwoDArray(size_type width, size_type height):
    _width(width),_height(height),_data(Array<T>(width*height, 0)){}
  TwoDArray(size_type width, size_type height, Array<T> data):
    _width(width),_height(height),_data(data){}
  TwoDArray():
    _width(0),_height(0),_data(Array<T>()){}
  virtual ~TwoDArray(){}


  inline ArrayView<T> operator[](size_type i);
  
  inline ArrayView<const T> operator[](size_type i) const;

  inline size_type getWidth() const{
    return _width;
  }

  inline size_type getHeight() const{
    return _height;
  }

  inline const Array<T>& getDataArray() const{
    return _data;
  }

  static const std::string& getDimensionsSeperator(){
    static const std::string dimensionSeperator = ":";
    return dimensionSeperator;
  }

  static const std::string& getDimensionsDelimiter(){
    static const std::string dimensionsDelimiter = "x";
    return dimensionsDelimiter;
  }

  static std::string toString(const TwoDArray<T> array);

  static TwoDArray<T> fromString(const std::string& string);

private:
  size_type _width, _height;
  Array<T> _data;
};


template<class T> inline
ArrayView<T> TwoDArray<T>::operator[](size_type i){
  return _data.view(_width*i, _width);  
}

template<class T> inline
ArrayView<const T> TwoDArray<T>::operator[](size_type i) const{
  return _data.view(_width*i, _width);  
}

template<class T>
std::string TwoDArray<T>::toString(const TwoDArray<T> array){
  std::stringstream widthStream;
  std::stringstream heightStream;
  widthStream << array.getWidth();
  heightStream << array.getHeight();
  return widthStream.str() +
    TwoDArray<T>::getDimensionsDelimiter() +
    heightStream.str() +
    TwoDArray<T>::getDimensionsSeperator() +
    array.getDataArray().toString();
}

template<class T>
TwoDArray<T> TwoDArray<T>::fromString(const std::string& string){
  size_t dimCharPos = 
    string.find_first_of(TwoDArray<T>::getDimensionsDelimiter());
  size_t seperatorPos = 
    string.find_first_of(TwoDArray<T>::getDimensionsSeperator());
  std::string valueData = string.substr(seperatorPos+1);
  std::istringstream widthStream(string.substr(0,dimCharPos));
  std::istringstream heightStream(string.substr(dimCharPos+1, seperatorPos-dimCharPos-1));
  size_t width, height;
  widthStream >> width;
  heightStream >> height;
  Array<T> array = fromStringToArray<T>(valueData);
  return TwoDArray<T>(width, height, array);
}

template<class T>
std::istringstream& operator>> (std::istringstream& in, TwoDArray<T>& array){
  array = TwoDArray<T>::fromString(in.str());
  return in;
}

template<class T> inline
std::ostream& operator<<(std::ostream& os, const TwoDArray<T>& array){
  return os << TwoDArray<T>::toString(array);
}

template<typename T> inline
bool operator==( const TwoDArray<T> &a1, const TwoDArray<T> &a2 ){
  return a1.getDataArray() == a2.getDataArray() && 
    a1.getWidth() == a2.getWidth() &&
    a1.getHeight() == a2.getHeight();
}

inline
std::string getTwoDArrayTypeNameTraitsFormat(){
  return "Array(*)";
}

template<typename T>
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<TwoDArray<T> > {
public:
  static std::string name(){ 
    std::string formatString = getTwoDArrayTypeNameTraitsFormat();
    size_t starPos = formatString.find("*");
    std::string prefix = formatString.substr(0,starPos);
    std::string postFix = formatString.substr(starPos+1);
    return prefix+TypeNameTraits<T>::name()+postFix;
  }
  static std::string concreteName(const TwoDArray<T>&)
    { return name(); }
};


} //namespace Teuchos


#endif // TEUCHOS_TWODARRAY_H
