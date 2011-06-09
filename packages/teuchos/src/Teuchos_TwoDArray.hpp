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

/**
 * A thin wrapper around the Array class which causes it to be interpurted
 * as a 2D Array.
 */
template<class T> 
class TwoDArray{
public:
  /**
   * \brief .
   */
  typedef Ordinal size_type;

  /** \name Constructors and Destructors */
  //@{

  /**
   * \brief Constructs a TwoDArray with the given number of rows and colums with each
   * entry being populated with the specified value.
   *
   * @param numCols The number of columns in the TwoDArray.
   * @param numRows The number of rows in the TwoDArray.
   * @param value The value with which to populate the TwoDArray.
   */
  TwoDArray(size_type numCols, size_type numRows, T value=T()):
    _numCols(numCols),_numRows(numRows),_data(Array<T>(numCols*numRows, value))
    {}
  /**
   * \brief Constructs an empty TwoDArray.
   */
  TwoDArray():
    _numCols(0),_numRows(0),_data(Array<T>()){}

  /** \brief . */
  virtual ~TwoDArray(){}

  //@}
  
  /** \name Getters and Setters */
  //@{
  
  /** \brief Returns an ArrayView containing the contents of row i */
  inline ArrayView<T> operator[](size_type i);
  
  /** \brief Returns a const ArrayView containing the contents of row i */
  inline const ArrayView<T> operator[](size_type i) const;

  /** \brief returns the number of columns in the TwoDArray. */
  inline size_type getNumCols() const{
    return _numCols;
  }

  /** \brief returns the number of rows in the TwoDArray. */
  inline size_type getNumRows() const{
    return _numRows;
  }

  /** \brief Returns the 1D array that is backing this TwoDArray. */
  inline const Array<T>& getDataArray() const{
    return _data;
  }

  //@}

  /** \name String conversion functions */
  //@{

  /** \brief returns the string used as the dimension dilimeter when convering
   * the TwoDArray to a string.
   */
  static const std::string& getDimensionsSeperator(){
    static const std::string dimensionSeperator = ":";
    return dimensionSeperator;
  }

  /** \brief returns the string used to seperate dimension information from
   * actual data information when converting a TwoDArray to a string.
   */
  static const std::string& getDimensionsDelimiter(){
    static const std::string dimensionsDelimiter = "x";
    return dimensionsDelimiter;
  }

  /** \brief Converts a given TwoDArray to a valid string representation. */
  static std::string toString(const TwoDArray<T> array);

  /** \brief Converts a valid string to it's corresponding TwoDArray. */
  static TwoDArray<T> fromString(const std::string& string);

  //@}
  
private:
  size_type _numCols, _numRows;
  Array<T> _data;
  TwoDArray(size_type numCols, size_type numRows, Array<T> data):
    _numCols(numCols),_numRows(numRows),_data(data){}
};

template<class T> inline
ArrayView<T> TwoDArray<T>::operator[](size_type i){
  return _data.view(_numCols*i, _numCols);  
}

template<class T> inline
const ArrayView<T> TwoDArray<T>::operator[](size_type i) const{
  return _data.view(_numCols*i, _numCols);  
}

template<class T>
std::string TwoDArray<T>::toString(const TwoDArray<T> array){
  std::stringstream numColsStream;
  std::stringstream numRowsStream;
  numColsStream << array.getNumCols();
  numRowsStream << array.getNumRows();
  return numColsStream.str() +
    TwoDArray<T>::getDimensionsDelimiter() +
    numRowsStream.str() +
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
  std::istringstream numColsStream(string.substr(0,dimCharPos));
  std::istringstream numRowsStream(string.substr(dimCharPos+1, seperatorPos-dimCharPos-1));
  size_t numCols, numRows;
  numColsStream >> numCols;
  numRowsStream >> numRows;
  Array<T> array = fromStringToArray<T>(valueData);
  return TwoDArray<T>(numCols, numRows, array);
}

/* \brief .
 * \relates TwoDArray
 */
template<class T>
std::istringstream& operator>> (std::istringstream& in, TwoDArray<T>& array){
  array = TwoDArray<T>::fromString(in.str());
  return in;
}

/* \brief .
 * \relates TwoDArray
 */
template<class T> inline
std::ostream& operator<<(std::ostream& os, const TwoDArray<T>& array){
  return os << TwoDArray<T>::toString(array);
}

/* \brief Returns true of the two TwoDArrays have the same contents and
 * their dimensions are the same.
 * \relates TwoDArray
 */
template<typename T> inline
bool operator==( const TwoDArray<T> &a1, const TwoDArray<T> &a2 ){
  return a1.getDataArray() == a2.getDataArray() && 
    a1.getNumCols() == a2.getNumCols() &&
    a1.getNumRows() == a2.getNumRows();
}

/**
 * \brief Get the format that is used for the specialization of the TypeName
 * traits class for TwoDArray. 
 *
 * The string returned will contain only one
 * "*" character. The "*" character should then be replaced with the actual 
 * template type of the array. 
 * \relates TwoDArray.
 */
inline
std::string getTwoDArrayTypeNameTraitsFormat(){
  return "TwoDArray(*)";
}

/** \brief TypeNameTraits specialization for Array.  */
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
