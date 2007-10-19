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

#ifndef TEUCHOS_ARRAY_H
#define TEUCHOS_ARRAY_H

/*! \file Teuchos_Array.hpp
  \brief Templated array class derived from the STL std::vector
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Utils.hpp"


namespace Teuchos {


/** \brief .
 * \relates Array
 */
class InvalidArrayStringRepresentation : public std::logic_error
{public:InvalidArrayStringRepresentation(const std::string& what_arg) : std::logic_error(what_arg) {}};


/** \brief .
 * \relates Array
 */
class RangeError : public std::range_error
{public:RangeError(const std::string& what_arg) : std::range_error(what_arg) {}};


/** \brief Memory-safe tempalted array class that encapsulates std::vector.
 *
 * ToDo: Finish documentation!
 */
template<typename T>
class Array
{
public:

  /** \name std::vector typedefs */
  //@{

  /** \brief . */
  typedef typename std::vector<T>::value_type value_type;
  /** \brief . */
  typedef typename std::vector<T>::pointer pointer;
  /** \brief . */
  typedef typename std::vector<T>::const_pointer const_pointer;
  /** \brief . */
  typedef typename std::vector<T>::reference reference;
  /** \brief . */
  typedef typename std::vector<T>::const_reference const_reference;

//#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
//  /** \brief . */
//  typedef ArrayRCP<T> iterator;
//  /** \brief . */
//  typedef typename ArrayRCP<const T> const_iterator;
//  /** \brief . */
//  typedef typename ArrayRCP<T> reverse_iterator;
//  /** \brief . */
//  typedef typename ArrayRCP<const T> const_reverse_iterator;
//#else
  /** \brief . */
  typedef typename std::vector<T>::iterator iterator;
  /** \brief . */
  typedef typename std::vector<T>::const_iterator const_iterator;
  /** \brief . */
  typedef typename std::vector<T>::reverse_iterator reverse_iterator;
  /** \brief . */
  typedef typename std::vector<T>::const_reverse_iterator const_reverse_iterator;
//#endif

  /** \brief . */
  typedef typename std::vector<T>::size_type size_type;
  /** \brief . */
  typedef typename std::vector<T>::difference_type difference_type;
  /** \brief . */
  typedef typename std::vector<T>::allocator_type allocator_type;

  //@}

  /** \name std::vector functions */
  //@{

  /** \brief . */
  Array();
  /** \brief . */
  explicit Array(size_type n, const value_type& value = value_type());
  /** \brief . */
  Array(const Array& x);
  /** \brief . */
  template<typename InputIterator>
  Array(InputIterator first, InputIterator last);
  /** \brief . */
  ~Array();
  /** \brief . */
  Array& operator=(const Array& x);
  /** \brief . */
  void assign(size_type n, const value_type& val);
  /** \brief . */
  template<typename InputIterator>
  void assign(InputIterator first, InputIterator last);
  /** \brief . */
  iterator begin();
  /** \brief . */
  iterator end();
  /** \brief . */
  const_iterator begin() const;
  /** \brief . */
  const_iterator end() const;
  /** \brief . */
  reverse_iterator rbegin();
  /** \brief . */
  reverse_iterator rend();
  /** \brief . */
  const_reverse_iterator rbegin() const;
  /** \brief . */
  const_reverse_iterator rend() const;
  /** \brief . */
  size_type size() const;
  /** \brief . */
  size_type max_size() const;
  /** \brief . */
  void resize(size_type new_size, const value_type& x = value_type());
  /** \brief . */
  size_type capacity() const;
  /** \brief . */
  bool empty() const;
  /** \brief . */
  void reserve(size_type n);
  /** \brief . */
  reference operator[](size_type i);
  /** \brief . */
  const_reference operator[](size_type i) const;
  /** \brief . */
  reference at(size_type i);
  /** \brief . */
  const_reference at(size_type i) const;
  /** \brief . */
  reference front();
  /** \brief . */
  const_reference front() const;
  /** \brief . */
  reference back();
  /** \brief . */
  const_reference back() const;
  /** \brief . */
  void push_back(const value_type& x);
  /** \brief . */
  void pop_back();
  /** \brief . */
  iterator insert(iterator position, const value_type& x);
  /** \brief . */
  void insert(iterator position, size_type n, const value_type& x);
  /** \brief . */
  template<typename InputIterator>
  void insert(iterator position, InputIterator first, InputIterator last);
  /** \brief . */
  iterator erase(iterator position);
  /** \brief . */
  iterator erase(iterator first, iterator last);
  /** \brief . */
  void swap(Array& x);
  /** \brief . */
  void clear();

  //@}

  /** \name Non-standard functions. */
  //@{

  /** \brief Add a new entry at the end of the array.
   *
   * Resize to allow space for the new entry.
   */
  Array<T>& append(const T& entry)
    {
      this->push_back(entry);
      return *this;
    }

  /** \brief Remove the i-th element from the array, with optional
   * boundschecking.
   */
  void remove(int i);

  /** \brief Return number of elements in the array.
   *
   * Equivalent to size(), but * included for backwards compatibility.
   */
  int length() const {return this->size();}

  /** \brief Convert an Array to an <tt>std::string</tt> */
  std::string toString() const;

  /** \brief Return true if Array has been compiled with boundschecking on. */ 
  static bool hasBoundsChecking();

  //@}

private:

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  RCP<std::vector<T> > vec_;
#else
  std::vector<T> vec_;
#endif

  std::vector<T>& vec()
    {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
      return *vec_;
#else
      return vec_;
#endif
    }

  const std::vector<T>& vec() const
    {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
      return *vec_;
#else
      return vec_;
#endif
    }
  
  void indexCheckCrash(int i) const;

};


/** \brief Write an Array to an ostream.
 *
 * This prints arrays in the form:
 
 \verbatim

 { 1.0, 2.0, 3.0 }

 \endverbatim

 * \relates Array 
 */
template<typename T>
std::ostream& operator<<(std::ostream& os, const Array<T>& array);

/** \relates Array */
template<typename T>
int hashCode(const Array<T>& array);


/** \brief Convert an array to a string representation.
 *
 * \relates Array
 */
template<typename T>
std::string toString(const Array<T>& array);


/** \brief Converts from std::string representation (as created by
 * <tt>toString()</tt>) back into the array object.
 *
 * \param arrayStr [in] The std::string representation of the array (see
 * below).
 *
 * <b>Exceptions:</b> If the std::string representation is not valid, then an
 * std::exception of type <tt>InvalidArrayStringRepresentation</tt> with be
 * thrown with a decent error message attached.
 *
 * The formating of the std::string <tt>arrayStr</tt> must look like:
   
 \verbatim

 {  val[0], val[1], val[2], val[3], ..., val[n-1] }

 \endverbatim

 * Currently <tt>operator>>()</tt> is used to convert the entries from their
 * std::string representation to objects of type <tt>T</tt>.  White space is
 * unimportant and the parser keys off of ',', '{' and '}' so even newlines
 * are allowed.  In the future, a traits class might be defined that will
 * allow for finer-grained control of how the conversion from strings to
 * values is performed in cases where <tt>operator>>()</tt> does not exist
 * for certain types.
 *
 * <b>Warning!</b> Currently this function only supports reading in flat
 * array objects for basic types like <tt>bool</tt>, <tt>int</tt>, and
 * <tt>double</tt> and does not yet support nested arrays (i.e. no
 * <tt>Array<Array<int> ></tt>) or other such fancy nested types.  Support
 * for nested arrays and other user defined types <tt>T</tt> can be added in
 * the future with no impact on user code.  Only the parser for the array
 * needs to be improved.  More specifically, the current implementation will
 * not work for any types <tt>T</tt> who's std::string representation contains
 * the characters <tt>','</tt> or <tt>'}'</tt>.  This implementation can be
 * modified to allow any such types by watching for the nesting of common
 * enclosing structures like <tt>[...]</tt>, <tt>{...}</tt> or
 * <tt>(...)</tt> within each entry of the std::string representation.  However,
 * this should all just work fine on most machines for the types
 * <tt>int</tt>, <tt>bool</tt>, <tt>float</tt>, <tt>double</tt> etc.
 *
 * <b>Warning!</b> Trying to read in an array in std::string format of doubles in
 * scientific notation such as <tt>{1e+2,3.53+6,...}</tt> into an array
 * object such as <tt>Array<int></tt> will not yield the correct results.
 * If one wants to allow a neutral std::string representation to be read in as an
 * <tt>Array<double></tt> object or an <tt>Array<int></tt> object, then
 * general formating such as <tt>{100,3530000,...}</tt> should be used.
 * This templated function is unable to deal std::complex type conversion issues.
 * 
 * \relates Array.
 */
template<typename T>
Array<T> fromStringToArray(const std::string& arrayStr);
                      
/** \brief Create an array with one entry.
 *
 * \relates Array
 */
template<typename T>
inline
Array<T> tuple(const T& a);

                      
/** \brief Create an array with two entries.
 *
 * \relates Array
 */
template<typename T>
inline
Array<T> tuple(const T& a, const T& b);


/** \brief Create an array with three entries.
 *
 * \relates Array
 */
template<typename T>
inline
Array<T> tuple(const T& a, const T& b, const T& c);


/** \brief Create an array with four entries.
 *
 * \relates Array
 */
template<typename T>
inline
Array<T> tuple(const T& a, const T& b, const T& c, const T& d);


/** \brief Create an array with five entries.
 *
 * \relates Array
 */
template<typename T>
inline
Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e);


/** \brief Create an array with six entries.
 *
 * \relates Array
 */
template<typename T>
inline
Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f);


/** \brief Create an array with seven entries.
 *
 * \relates Array
 */
template<typename T>
inline
Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g);


/** \brief Create an array with eight entries.
 *
 * \relates Array
 */
template<typename T>
inline
Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h);


/** \brief Create an array with nine entries.
 *
 * \relates Array
 */
template<typename T>
inline
Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i);


/** \brief Create an array with ten entries.
 *
 * \relates Array
 */
template<typename T> inline
Array<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i, const T& j);


} // namespace Teuchos


//
// Implementation
//


namespace Teuchos {


// std::vector functions


template<typename T>
inline
Array<T>::Array() :
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  vec_(rcp(new std::vector<T>()))
#endif
{}


template<typename T>
inline
Array<T>::Array(size_type n, const value_type& value) :
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  vec_(rcp(new std::vector<T>(n,value)))
#else
  vec_(n,value)
#endif
{}


template<typename T>
inline
Array<T>::Array(const Array& x) :
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  vec_(rcp(new std::vector<T>(x)))
#else
  vec_(x)
#endif
{}


template<typename T>
template<typename InputIterator>
inline
Array<T>::Array(InputIterator first, InputIterator last) :
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  vec_(rcp(new std::vector<T>(first,last)))
#else
  vec_(first,last)
#endif
{}


template<typename T>
inline
Array<T>::~Array()
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  TEST_FOR_EXCEPTION(
    vec_.count() > 1, std::runtime_error,
    "Error, there is some client with a dangling iterator when the object that "
    "created it is being deleted!"
    );
#endif
}


template<typename T>
inline
Array<T>& Array<T>::operator=(const Array& x)
{
  vec() = x.vec();
  return *this;
}


template<typename T>
inline
void Array<T>::assign(size_type n, const value_type& val)
{
  vec().assign(n,val);
}


template<typename T>
template<typename InputIterator>
inline
void Array<T>::assign(InputIterator first, InputIterator last)
{
  vec().assign(first,last);
}


template<typename T>
inline
typename Array<T>::iterator
Array<T>::begin()
{
  return vec().begin();
}


template<typename T>
inline
typename Array<T>::iterator
Array<T>::end()
{
  return vec().end();
}


template<typename T>
inline
typename Array<T>::const_iterator
Array<T>::begin() const
{
  return vec().begin();
}


template<typename T>
inline
typename Array<T>::const_iterator
Array<T>::end() const
{
  return vec().end();
}


template<typename T>
inline
typename Array<T>::reverse_iterator
Array<T>::rbegin()
{
  return vec().rbegin();
}


template<typename T>
inline
typename Array<T>::reverse_iterator
Array<T>::rend()
{
  return vec().rend();
}


template<typename T>
inline
typename Array<T>::const_reverse_iterator
Array<T>::rbegin() const
{
  return vec().rbegin();
}


template<typename T>
inline
typename Array<T>::const_reverse_iterator
Array<T>::rend() const
{
  return vec().rend();
}


template<typename T>
inline
typename Array<T>::size_type
Array<T>::size() const
{
  return vec().size();
}


template<typename T>
inline
typename Array<T>::size_type
Array<T>::max_size() const
{
  return vec().max_size();
}


template<typename T>
inline
void
Array<T>::resize(size_type new_size, const value_type& x)
{
  vec().resize(new_size,x);
}


template<typename T>
inline
typename Array<T>::size_type
Array<T>::capacity() const
{
  return vec().capacity();
}


template<typename T>
inline
bool Array<T>::empty() const
{
  return vec().empty();
}


template<typename T>
inline
void Array<T>::reserve(size_type n)
{
  vec().reserve(n);
}


template<typename T>
inline
typename Array<T>::reference
Array<T>::operator[](size_type i)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  indexCheckCrash(i);
#endif
  return vec()[i];
}


template<typename T>
inline
typename Array<T>::const_reference
Array<T>::operator[](size_type i) const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  indexCheckCrash(i);
#endif
  return vec()[i];
}


template<typename T>
inline
typename Array<T>::reference
Array<T>::at(size_type i)
{
  return vec().at(i);
}


template<typename T>
inline
typename Array<T>::const_reference
Array<T>::at(size_type i) const
{
  return vec().at(i);
}


template<typename T>
inline
typename Array<T>::reference
Array<T>::front()
{
  return vec().front();
}


template<typename T>
inline
typename Array<T>::const_reference
Array<T>::front() const
{
  return vec().front();
}


template<typename T>
inline
typename Array<T>::reference
Array<T>::back()
{
  return vec().back();
}


template<typename T>
inline
typename Array<T>::const_reference
Array<T>::back() const
{
  return vec().back();
}


template<typename T>
inline
void Array<T>::push_back(const value_type& x)
{
  vec().push_back(x);
}


template<typename T>
inline
void Array<T>::pop_back()
{
  vec().pop_back();
}


template<typename T>
inline
typename Array<T>::iterator
Array<T>::insert(iterator position, const value_type& x)
{
  return vec().insert(position,x);
}


template<typename T>
inline
void Array<T>::insert(iterator position, size_type n, const value_type& x)
{
  vec().insert(position,n,x);
}


template<typename T>
template<typename InputIterator>
inline
void Array<T>::insert(iterator position, InputIterator first, InputIterator last)
{
  vec().insert(position,first,last);
}


template<typename T>
inline
typename Array<T>::iterator
Array<T>::erase(iterator position)
{
  return vec().erase(position);
}


template<typename T>
inline
typename Array<T>::iterator
Array<T>::erase(iterator first, iterator last)
{
  return vec().erase(first,last);
}


template<typename T>
inline
void Array<T>::swap(Array& x)
{
  vec().swap(x.vec());
}


template<typename T>
inline
void Array<T>::clear()
{
  vec().clear();
}


// Non-standard functions


template<typename T>
void Array<T>::remove(int i)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  indexCheckCrash(i);
#endif
  // Erase the i-th element of this array.
  this->erase( this->begin() + i );
}


template<typename T>
inline
std::string Array<T>::toString() const
{
  std::ostringstream ss;
  ss << "{";

  for (int i=0; i<length(); i++)
  {
    ss << operator[](i);
    if (i<length()-1) ss << ", ";
  }
  ss << "}";

  return ss.str();
}


template<typename T> inline
bool Array<T>::hasBoundsChecking()
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK  
  return true;
#else
  return false;
#endif
}


template<typename T> inline
void Array<T>::indexCheckCrash(int i) const
{
  TEST_FOR_EXCEPTION(
    !( 0 <= i && i < length() ), RangeError,
    "Array<"<<TypeNameTraits<T>::name()<<">::indexCheckCrash: "
    "index " << i << " out of range [0, "<< length() << ")"
    );
}


} // namespace Teuchos


// Nonmember functions


template<typename T>
inline
std::ostream& Teuchos::operator<<(
  std::ostream& os, const Array<T>& array
  )
{
  return os << Teuchos::toString(array);
}


template<typename T>
inline
int Teuchos::hashCode(const Array<T>& array)
{
  int rtn = hashCode(array.length());
  for (int i=0; i<array.length(); i++)
  {
    rtn += hashCode(array[i]);
  }
  return rtn;
}


template<typename T>
inline
std::string Teuchos::toString(const Array<T>& array)
{
  return array.toString();
}


template<typename T>
Teuchos::Array<T>
Teuchos::fromStringToArray(const std::string& arrayStr)
{
  const std::string str = Utils::trimWhiteSpace(arrayStr);
  std::istringstream iss(str);
  TEST_FOR_EXCEPTION(
    ( str[0]!='{' || str[str.length()-1] != '}' )
    ,InvalidArrayStringRepresentation
    ,"Error, the std::string:\n"
    "----------\n"
    <<str<<
    "\n----------\n"
    "is not a valid array represntation!"
    );
  char c;
  c = iss.get(); // Read initial '{'
  TEST_FOR_EXCEPT(c!='{'); // Should not throw!
  // Now we are ready to begin reading the entries of the array!
  Array<T> a;
  while( !iss.eof() ) {
    // Get the basic entry std::string
    std::string entryStr;
    std::getline(iss,entryStr,','); // Get next entry up to ,!
    // ToDo: Above, we might have to be careful to look for the opening and
    // closing of parentheses in order not to pick up an internal ',' in the
    // middle of an entry (for a std::complex number for instance).  The above
    // implementation assumes that there will be no commas in the middle of
    // the std::string representation of an entry.  This is certainly true for
    // the types bool, int, float, and double.
    //
    // Trim whitespace from beginning and end
    entryStr = Utils::trimWhiteSpace(entryStr);
    // Remove the final '}' if this is the last entry and we did not
    // actually terminate the above getline(...) on ','
    bool found_end = false;
    if(entryStr[entryStr.length()-1]=='}') {
      entryStr = entryStr.substr(0,entryStr.length()-1);
      found_end = true;
      if( entryStr.length()==0 && a.size()==0 )
        return a; // This is the empty array "{}" (with any spaces in it!)
    }
    TEST_FOR_EXCEPTION(
      0 == entryStr.length()
      ,InvalidArrayStringRepresentation
      ,"Error, the std::string:\n"
      "----------\n"
      <<str<<
      "\n----------\n"
      "is not a valid array represntation!"
      );
    // Finally we can convert the entry and add it to the array!
    std::istringstream entryiss(entryStr);
    T entry;
    entryiss >> entry; // Assumes type has operator>>(...) defined!
    // ToDo: We may need to define a traits class to allow us to specialized
    // how conversion from a std::string to a object is done!
    a.push_back(entry);
    // At the end of the loop body here, if we have reached the last '}'
    // then the input stream iss should be empty and iss.eof() should be
    // true, so the loop should terminate.  We put an std::exception test here
    // just in case something has gone wrong.
    TEST_FOR_EXCEPTION(
      found_end && !iss.eof()
      ,InvalidArrayStringRepresentation
      ,"Error, the std::string:\n"
      "----------\n"
      <<str<<
      "\n----------\n"
      "is not a valid array represntation!"
      );
  }
  return a;
}

                      
template<typename T>
inline
Teuchos::Array<T>
Teuchos::tuple(const T& a)
{
  Array<T> rtn(1, a);
  return rtn;
}


template<typename T>
inline
Teuchos::Array<T>
Teuchos::tuple(const T& a, const T& b)
{
  Array<T> rtn(2);
  rtn[0] = a;
  rtn[1] = b;
  return rtn;
}


template<typename T>
inline
Teuchos::Array<T>
Teuchos::tuple(const T& a, const T& b, const T& c)
{
  Array<T> rtn(3);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  return rtn;
}


template<typename T>
inline
Teuchos::Array<T>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d)
{
  Array<T> rtn(4);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  return rtn;
}


template<typename T>
inline
Teuchos::Array<T>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e)
{
  Array<T> rtn(5);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  return rtn;
}


template<typename T>
inline
Teuchos::Array<T>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f)
{
  Array<T> rtn(6);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  return rtn;
}


template<typename T>
inline
Teuchos::Array<T>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g)
{
  Array<T> rtn(7);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  return rtn;
}


template<typename T>
inline
Teuchos::Array<T>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h)
{
  Array<T> rtn(8);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  rtn[7] = h;
  return rtn;
}


template<typename T>
inline
Teuchos::Array<T>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i)
{
  Array<T> rtn(9);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  rtn[7] = h;
  rtn[8] = i;
  return rtn;
}


template<typename T>
inline
Teuchos::Array<T>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i, const T& j)
{
  Array<T> rtn(10);
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  rtn[7] = h;
  rtn[8] = i;
  rtn[9] = j;
  return rtn;
}


#endif // TEUCHOS_ARRAY_H
