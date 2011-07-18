// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
 * \brief A thin wrapper around the Array class which causes it to be interpurted
 * as a 2D Array.
 *
 * 2D Array's can also be "symetric". This means that anyone viewing the
 * Array should only consider the lower half of the array as valid. The
 * 
 * \warning The TwoDArray will not enforce symetry. However, when two 
 * symetrical TwoDArrays are compared, only the the lower half of the
 * TwoDArray's will be compared.
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
  TwoDArray(size_type numRows, size_type numCols, T value=T()):
    _numRows(numRows),
    _numCols(numCols),
    _data(Array<T>(numCols*numRows, value)),
    _symetrical(false)
    {}
  /**
   * \brief Constructs an empty TwoDArray.
   */
  TwoDArray():
    _numRows(0),_numCols(0),_data(Array<T>()),_symetrical(false){}

  /** \brief . */
  virtual ~TwoDArray(){}

  //@}
  
  /** \name Getters and Setters */
  //@{
  
  /** \brief Returns an ArrayView containing the contents of row i */
  inline ArrayView<T> operator[](size_type i);
  
  /** \brief Returns a const ArrayView containing the contents of row i */
  inline const ArrayView<T> operator[](size_type i) const;

  /** \brief returns the number of rows in the TwoDArray. */
  inline size_type getNumRows() const{
    return _numRows;
  }

  /** \brief returns the number of columns in the TwoDArray. */
  inline size_type getNumCols() const{
    return _numCols;
  }

  /** \brief Returns the 1D array that is backing this TwoDArray. */
  inline const Array<T>& getDataArray() const{
    return _data;
  }

  /** \brief Returns the element located at i,j */
  inline T& operator()(size_type i, size_type j){
    return _data[(i*_numCols)+j];
  }

  /** \brief Returns the element located at i,j */
  inline const T& operator()(size_type i, size_type j) const{
    return _data[(i*_numCols)+j];
  }

  /** \brief delets all the entries from the TwoDArray */
  inline void clear(){
    _data.clear();
    _numRows =0;
    _numCols =0;
  }

  inline bool isEmpty(){
    return _numRows == 0 &&
      _numCols == 0 &&
      _data.size() == 0;
  }

  /** \brief A simple flag indicating whether or not
   * this TwoDArray should be interpurted as symetrical.
   *
   * \note A symetrical TwoDArray is defined as an TwoDArray where
   * entry i,j is the same as entry j,i.
   *
   * \note This does not change any of the TwoDArrays behavior.
   * It merely serves as an indicator to any one using a TwoDArray
   * that the TwoDArray can be read as if it were symetrical. In other words,
   * the TwoDArray class does not enforce the symetry.
   *
   * @return True if the array is "symetrical", false otherwise.
   */
  inline bool isSymetrical() const{
    return _symetrical;
  }

  /**
   * \brief Sets whether or not the the TwoDArray should be 
   * interpurted as symetric.
   *
   * \note A symetrical TwoDArray is defined as an TwoDArray where
   * entry i,j is the same as entry j,i.
   *
   * \note This does not change any of the TwoDArrays behavior.
   * It merely serves as an indicator to any one using a TwoDArray
   * that the TwoDArray can be read as if it were symetrical. In other words,
   * the TwoDArray class does not enforce the symetry.
   *
   * @param symetry Whether or not the matrix should be interpurted
   * as symetric.
   */
  inline void setSymetrical(bool symetrical){
    _symetrical = symetrical;
  }

  //@}

  /** @name Resizing Functions */
  //@{

  /** \brief Changes the number of rows in the matrix.
   *
   * If the new number of rows is less than the current number, the
   * last rows in the array will be deleted (i.e. if an array has 10 rows and
   * it is resized to have only 5 rows, rows 5-9 are deleted). If the new number
   * of rows is greater than the current number of rows, the rows are appended on
   * to the end of the array and the new entries are initialized to T's
   * default value.
   *
   * @param numberOfRows The new number of rows the TwoDArray should have.
   */
  void resizeRows(size_type numberOfRows);

  /** \brief Changes the number of rows in the matrix.
   *
   * If the new number of columns is less than the current number, the
   * last columns in the array will be deleted (i.e. if an array has 10 columns and
   * it is resized to have only 5 columns, columns 5-9 are deleted). If the new number
   * of columns is greater than the current number of columns, the columns are appended on
   * to the end of the array and the new entries are initialized to T's default value.
   *
   * \warning This operation has the potential to be very expensive 
   * as it essentially creates an entirely new 2DArray.
   * Please take this into account when using this function.
   *
   * @param numberOfCols The new number of rows the TwoDArray should have.
   */
  void resizeCols(size_type numberOfCols);

  //@}

  
  /** \name String conversion functions */
  //@{

  /** \brief returns the string used to seperate meta information from
   * actual data information when converting a TwoDArray to a string.
   */
  static const std::string& getMetaSeperator(){
    static const std::string metaSeperator = ":";
    return metaSeperator;
  }

  /** \brief returns the string used as the dimension dilimeter when convering
   * the TwoDArray to a string.
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
  size_type _numRows,_numCols;
  Array<T> _data;
  TwoDArray(size_type numRows, size_type numCols, Array<T> data):
    _numRows(numRows),_numCols(numCols),_data(data){}
  bool _symetrical;
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
void TwoDArray<T>::resizeRows(size_type numberOfRows){
  _data.resize(_numCols*numberOfRows);
  _numRows = numberOfRows;
}


template<class T>
void TwoDArray<T>::resizeCols(size_type numberOfCols){
  Array<T> newData(numberOfCols*_numRows);
  size_type colLimit = (numberOfCols < _numCols ? numberOfCols : _numCols);
  for(size_type i = 0; i<_numRows; i++){
    for(size_type j = 0; j<colLimit; j++){
      newData[i*numberOfCols+j] = _data[i*_numCols+j];
    }
  }
  _data = newData;
  _numCols=numberOfCols;
}


template<class T>
std::string TwoDArray<T>::toString(const TwoDArray<T> array){
  std::stringstream numColsStream;
  std::stringstream numRowsStream;
  numColsStream << array.getNumCols();
  numRowsStream << array.getNumRows();
  std::string metaSeperator = TwoDArray<T>::getMetaSeperator();
  return 
    numRowsStream.str() +
    TwoDArray<T>::getDimensionsDelimiter() +
    numColsStream.str() +
    metaSeperator +
    (array.isSymetrical() ? "sym"+metaSeperator : "") +
    array.getDataArray().toString();
}

template<class T>
TwoDArray<T> TwoDArray<T>::fromString(const std::string& string){
  std::string curString = string;
  std::string metaSeperator = TwoDArray<T>::getMetaSeperator();
  size_t curPos = curString.find(metaSeperator);
  std::string dimString = curString.substr(0, curPos);
  curString = curString.substr(curPos+1);

  //process dimensions
  size_t dimCharPos = 
    dimString.find(TwoDArray<T>::getDimensionsDelimiter());
  std::istringstream numRowsStream(dimString.substr(0,dimCharPos));
  std::istringstream numColsStream(dimString.substr(dimCharPos+1));
  size_t numRows, numCols;
  numRowsStream >> numRows;
  numColsStream >> numCols;

  //determine symetry state
  bool symetrical = false;
  curPos = curString.find(metaSeperator);
  if(curPos != std::string::npos){
    symetrical = true;
    curString = curString.substr(curPos+1);
  }

  //Get actual data
  Array<T> array = fromStringToArray<T>(curString);

  TEST_FOR_EXCEPTION(array.size() != (typename Array<T>::size_type)(numRows*numCols),
    InvalidArrayStringRepresentation,
    "Error: You've specified an TwoDArray as having the dimensions of "
    << numRows << "x" << numCols << ". This means you should have " <<
    (numRows*numCols) << " entries specified in your array. However you "
    "only specified " << array.size() << " entries."
  )

  //Construct object to return
  TwoDArray<T> toReturn(numRows, numCols, array);
  toReturn.setSymetrical(symetrical);
  return toReturn;
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


namespace TwoDDetails {

/**
 * \brief A function for comparing symetrical arrarys.
 *
 * @param a1 The first array to compare.
 * @param a2 The second array to compare.
 * @return True if the two TwoDArrays are symetricaly the same, false 
 * otherwise.
 */
template<typename T>
bool symetricCompare(const TwoDArray<T> &a1, const TwoDArray<T> &a2 ){
  if(a1.getNumRows() != a2.getNumRows() || 
    a1.getNumRows() != a2.getNumRows())
  {
    return false;
  }
  else{
    typedef typename TwoDArray<T>::size_type ST;
    for(ST i=0;i<a1.getNumRows(); ++i){
      for(ST j=0;j<a1.getNumCols()-a1.getNumRows()+i; ++j){
        if(a1(i,j) != a2(i,j)){
          return false;
        }
      }
    }
    return true;
  }
}


}

/* \brief Returns true of the two TwoDArrays have the same contents and
 * their dimensions are the same.
 *
 * \note If the arrays are symetrical, only the values in the upper half
 * of the array are compared. For example: in a 4x4 array, only the values 
 * indicated with x's in the figure below would be compared.
 *
 *   o o o o
 *   x o o o 
 *   x x o o  
 *   x x x o
 *
 * \relates TwoDArray
 */
template<typename T> 
bool operator==( const TwoDArray<T> &a1, const TwoDArray<T> &a2 ){
  if(a1.isSymetrical() != a2.isSymetrical()){
    return false;
  }
  if(a1.isSymetrical()){
    return TwoDDetails::symetricCompare(a1,a2);
  }
  else{
    return a1.getDataArray() == a2.getDataArray() && 
    a1.getNumRows() == a2.getNumRows() &&
    a1.getNumCols() == a2.getNumCols();
  }
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
