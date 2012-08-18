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

#ifndef TEUCHOS_ARRAY_H
#define TEUCHOS_ARRAY_H

/*! \file Teuchos_Array.hpp
  \brief Templated array class derived from the STL std::vector
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_Utils.hpp"
#include "Teuchos_Assert.hpp"


namespace Teuchos {


/** \brief .
 *
 * \ingroup teuchos_mem_mng_grp
 */
class InvalidArrayStringRepresentation : public std::logic_error
{public:InvalidArrayStringRepresentation(const std::string& what_arg) : std::logic_error(what_arg) {}};


template<typename T> class Array;


// 2007/11/30: rabartl: Below, I had to move the initial declaration of these
// non-member template functions outside of the Array class since the Sun
// compiler on sass9000 would not accept this.  However, this did work on a
// number of other compilers such a g++, Intel C++ etc.  The old in-class
// non-member friend definition is clearly ISO 98 C++ as shown in Item 46 of
// "Effective C++: Third Edition".  This is not the end of the world but this
// is something to remember for this platform.


/** \brief Equality operator.
 *
 * \relates Array
 */
template<typename T> inline
bool operator==( const Array<T> &a1, const Array<T> &a2 );


/** \brief Non-equality operator.
 *
 * \relates Array
 */
template<typename T> inline
bool operator!=( const Array<T> &a1, const Array<T> &a2 );


/** \brief Non-member swap (specializes default std version).
 *
 * \relates Array
 */
template<typename T> inline
void swap( Array<T> &a1, Array<T> &a2 );


/** \brief Less-than operator.
 *
 * \relates Array
 */
template<typename T> inline
bool operator<( const Array<T> &a1, const Array<T> &a2 );


/** \brief Less-than-or-equal operator.
 *
 * \relates Array
 */
template<typename T> inline
bool operator<=( const Array<T> &a1, const Array<T> &a2 );


/** \brief Greater-than operator.
 *
 * \relates Array
 */
template<typename T> inline
bool operator>( const Array<T> &a1, const Array<T> &a2 );


/** \brief Greater-than-or-equal operator.
 *
 * \relates Array
 */
template<typename T> inline
bool operator>=( const Array<T> &a1, const Array<T> &a2 );


/** \brief Memory-safe templated array class that encapsulates std::vector.
 *
 * ToDo: Finish documentation!
 *
 * \section Teuchos_Array_Tuple_sec Tuple Construction
 *
 * A user can create a Teuchos::Tuple object to initialize an Array object by
 * using one of the the convenient overloaded Teuchos::tuple() non-member
 * constructor functions.  For example, see Array_test.cpp for how this is
 * done.
 *
 * \section Teuchos_Array_DesignDiscussion_sec Design Discussion
 *
 * Currently, this class defines implicit conversions to ArrayView.  An
 * alternative design would be to have Array derive from ArrayView.  This is a
 * workable design but it would impart some extra storage and runtime
 * overhead.  Perhaps the most significant overhead would be having the reset
 * the base ArrayView pointer and size on each and every change in the
 * structure of the container.  This would import extra overhead beyond a
 * straight std::vector.
 *
 * The big advantage of deriving Array from ArrayView is that this would allow
 * Array to be used to call some functions taking ArrayView without requiring
 * an implicit conversion.  While the implicit shallow conversion from Array
 * to ArrayView is very cheap (just a pointer and int copy), it does cause
 * problems where the compiler will refuse to perform an implicit conversion
 * to call a templated function.  However, note that an implicit conversion to
 * an ArrayView<const T> would always have to be performed no matter what.
 *
 * In summary, having Array implicitly convert to ArrayView instead of having
 * Array derive from ArrayView results in faster and simpler code at the
 * expense of the compiler refusing the make implicit conversions in some
 * cases when calling template functions.  Such conversion problems can always
 * be dealt with by using explicit template arguments.
 *
 * \ingroup teuchos_mem_mng_grp
 */
template<typename T>
class Array
{
public:

  // 2007/11/30: rabartl: Below, note that the only reason that these
  // functions are declared as friends is so that the compiler will do
  // automatic type conversions as described in "Effective C++: Third Edition"
  // Item 46.

  /** \brief . */
  template<typename T2>
  friend bool Teuchos::operator==( const Array<T2> &a1, const Array<T2> &a2 );

  /** \brief . */
  template<typename T2>
  friend bool Teuchos::operator!=( const Array<T2> &a1, const Array<T2> &a2 );

  /** \brief . */
  template<typename T2>
  friend void swap( Array<T2> &a1, Array<T2> &a2 );

  /** \brief . */
  template<typename T2>
  friend bool Teuchos::operator<( const Array<T2> &a1, const Array<T2> &a2 );

  /** \brief . */
  template<typename T2>
  friend bool Teuchos::operator<=( const Array<T2> &a1, const Array<T2> &a2 );

  /** \brief . */
  template<typename T2>
  friend bool Teuchos::operator>( const Array<T2> &a1, const Array<T2> &a2 );

  /** \brief . */
  template<typename T2>
  friend bool Teuchos::operator>=( const Array<T2> &a1, const Array<T2> &a2 );

  /** \name std::vector typedefs */
  //@{

  /** \brief. */
  typedef Teuchos_Ordinal Ordinal;
  /** \brief . */
  typedef Ordinal size_type;
  /** \brief . */
  typedef Ordinal difference_type;
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
  /** \brief . */
  typedef typename std::vector<T>::allocator_type allocator_type;

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  /** \brief . */
  typedef ArrayRCP<T> iterator;
  /** \brief . */
  typedef ArrayRCP<const T> const_iterator;
  /** \brief . */
  typedef std::reverse_iterator<iterator> reverse_iterator;
  /** \brief . */
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
#else
  /** \brief . */
  typedef typename std::vector<T>::iterator iterator;
  /** \brief . */
  typedef typename std::vector<T>::const_iterator const_iterator;
  /** \brief . */
  typedef typename std::vector<T>::reverse_iterator reverse_iterator;
  /** \brief . */
  typedef typename std::vector<T>::const_reverse_iterator const_reverse_iterator;
#endif


  //@}

  /** \name All constructors */
  //@{

  /** \brief . */
  inline Array();
  /** \brief . */
  inline explicit Array(size_type n, const value_type& value = value_type());
  /** \brief . */
  inline Array(const Array<T>& x);
  /** \brief . */
  template<typename InputIterator>
  inline Array(InputIterator first, InputIterator last);
  /** \brief . */
  inline Array(const ArrayView<const T>& a);
  /** \brief . */
  template<int N>
  inline Array(const Tuple<T,N>& t);
  /** \brief . */
  inline ~Array();
  /** \brief . */
  inline Array& operator=(const Array<T>& a);

  //@}

  /** \name Other std::vector functions */
  //@{

  /** \brief . */
  inline void assign(size_type n, const value_type& val);
  /** \brief . */
  template<typename InputIterator>
  inline void assign(InputIterator first, InputIterator last);
  /** \brief . */
  inline iterator begin();
  /** \brief . */
  inline iterator end();
  /** \brief . */
  inline const_iterator begin() const;
  /** \brief . */
  inline const_iterator end() const;
  /** \brief . */
  inline reverse_iterator rbegin();
  /** \brief . */
  inline reverse_iterator rend();
  /** \brief . */
  inline const_reverse_iterator rbegin() const;
  /** \brief . */
  inline const_reverse_iterator rend() const;
  /** \brief . */
  inline size_type size() const;
  /** \brief . */
  inline size_type max_size() const;
  /** \brief . */
  inline void resize(size_type new_size, const value_type& x = value_type());
  /** \brief . */
  inline size_type capacity() const;
  /** \brief . */
  inline bool empty() const;
  /** \brief . */
  inline void reserve(size_type n);
  /** \brief . */
  inline reference operator[](size_type i);
  /** \brief . */
  inline const_reference operator[](size_type i) const;
  /** \brief . */
  inline reference at(size_type i);
  /** \brief . */
  inline const_reference at(size_type i) const;
  /** \brief . */
  inline reference front();
  /** \brief . */
  inline const_reference front() const;
  /** \brief . */
  inline reference back();
  /** \brief . */
  inline const_reference back() const;
  /** \brief . */
  inline void push_back(const value_type& x);
  /** \brief . */
  inline void pop_back();
  /** \brief . */
  inline iterator insert(iterator position, const value_type& x);
  /** \brief . */
  inline void insert(iterator position, size_type n, const value_type& x);
  /** \brief . */
  template<typename InputIterator>
  inline void insert(iterator position, InputIterator first, InputIterator last);
  /** \brief . */
  inline iterator erase(iterator position);
  /** \brief . */
  inline iterator erase(iterator first, iterator last);
  /** \brief . */
  inline void swap(Array& x);
  /** \brief . */
  inline void clear();

  //@}

  /** \name General non-standard functions. */
  //@{

  /** \brief Add a new entry at the end of the array.
   *
   * Resize to allow space for the new entry.
   */
  inline Array<T>& append(const T& x);

  /** \brief Remove the i-th element from the array, with optional
   * boundschecking.
   */
  inline void remove(int i);

  /** \brief Return number of elements in the array.
   *
   * Equivalent to size(), but * included for backwards compatibility.
   */
  inline int length() const;

  /** \brief Convert an Array to an <tt>std::string</tt> */
  inline std::string toString() const;

  /** \brief Return true if Array has been compiled with boundschecking on. */
  inline static bool hasBoundsChecking();

  /** \brief Return a raw pointer to beginning of array or NULL if unsized. */
  inline T* getRawPtr();

  /** \brief Return a const raw pointer to beginning of array or NULL if unsized. */
  inline const T* getRawPtr() const;

  //@}

  /** \name Conversions to and from std::vector. */
  //@{

  /** \brief Copy constructor from an std::vector. */
  inline Array( const std::vector<T> &v );

  /** \brief Explicit copy conversion to an std::vector. */
  inline std::vector<T> toVector() const;

  /** \brief Assignment operator for std::vector. */
  inline Array& operator=( const std::vector<T> &v );

  //@}

  //! @name Views
  //@{

        /** \brief Return non-const view of a contiguous range of elements.
         *
         * <b>Preconditions:</b><ul>
   * <li><tt>0 <= offset && offset + size <= this->size()</tt>
         * </ul>
         *
         * <b>Postconditions:</b><ul>
   * <li><tt>returnVal.size() == size</tt>
         * </ul>
   *
   * NOTE: A <tt>size==0</tt> view of even an empty Array is allowed and
   * returns a <tt>null</tt> view.
   */
        inline ArrayView<T> view( size_type offset, size_type size );

        /** \brief Return const view of a contiguous range of elements.
         *
         * <b>Preconditions:</b><ul>
   * <li><tt>0 <= offset && offset + size <= this->size()</tt>
         * </ul>
         *
         * <b>Postconditions:</b><ul>
   * <li><tt>returnVal.size() == size</tt>
         * </ul>
   *
   * NOTE: A <tt>size==0</tt> view of even an empty Array is allowed and
   * returns a <tt>null</tt> view.
   */
        inline ArrayView<const T> view( size_type offset, size_type size ) const;

        /** \brief Return a non-const view of a contiguous range of elements (calls
   * view(offset,size)).
   */
        inline ArrayView<T> operator()( size_type offset, size_type size );

        /** \brief Return a const view of a contiguous range of elements (calls
   * view(offset,size)).
   */
        inline ArrayView<const T> operator()( size_type offset, size_type size ) const;

        /** \brief Return an non-const ArrayView of *this.
   *
   * NOTE: This will return a null ArrayView if this->size() == 0.
   */
        inline ArrayView<T> operator()();

        /** \brief Return an const ArrayView of *this.
   *
   * NOTE: This will return a null ArrayView if this->size() == 0.
   */
        inline ArrayView<const T> operator()() const;

  /** \brief Perform an implicit conversion to a non-const ArrayView (calls
   * operator()()).
   */
        inline operator ArrayView<T>();

  /** \brief Perform an implicit conversion to a non-const ArrayView (calls
   * operator()()).
   */
        inline operator ArrayView<const T>() const;

  //@}

private:

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  RCP<std::vector<T> > vec_;
  mutable ArrayRCP<T> extern_arcp_;
  mutable ArrayRCP<const T> extern_carcp_;
#else
  std::vector<T> vec_;
#endif

  inline std::vector<T>& vec(
    bool isStructureBeingModified = false,
    bool activeIter = false
    );

  inline const std::vector<T>& vec() const;

  inline typename std::vector<T>::iterator
  raw_position( iterator position );

  inline void assertIndex(int i) const;

  inline void assertNotNull() const;

};


/** \brief Wrap an <tt>RCP<Array<T> ></tt> object as an <tt>ArrayRCP<T></tt>
 * object.
 *
 * \relates ArrayRCP
 */
template<class T>
ArrayRCP<T> arcp( const RCP<Array<T> > &v )
{
  if ( is_null(v) || !v->size() )
    return null;
  return arcpWithEmbeddedObjPostDestroy<T,RCP<Array<T> > >(
    &(*v)[0], 0, v->size(),
    v, false
    );
}


/** \brief Wrap a <tt>RCP<const Array<T> ></tt> object as an
 * <tt>ArrayRCP<const T></tt> object.
 *
 * \relates ArrayRCP
 */
template<class T>
ArrayRCP<const T> arcp( const RCP<const Array<T> > &v )
{
  if ( is_null(v) || !v->size() )
    return null;
  return arcpWithEmbeddedObjPostDestroy<const T,RCP<const Array<T> > >(
    &(*v)[0], 0, v->size(),
    v, false
    );
}


/** \brief Wrap an <tt>Array<T></tt> object as a non-owning
 * <tt>ArrayRCP<T></tt> object.
 *
 * \relates ArrayRCP
 */
template<class T>
ArrayRCP<T> arcpFromArray( Array<T> &a )
{
  if (a.size() == 0)
    return null;
#ifdef TEUCHOS_DEBUG
  return a.begin(); // Catch dangling reference!
#else
  return arcp(a.getRawPtr(), 0, a.size(), false);
#endif
}


/** \brief Wrap a <tt>const Array<T></tt> object as a non-owning
 * <tt>ArrayRCP<T></tt> object.
 *
 * \relates ArrayRCP
 */
template<class T>
ArrayRCP<const T> arcpFromArray( const Array<T> &a )
{
  if (a.size() == 0)
    return null;
#ifdef TEUCHOS_DEBUG
  return a.begin(); // Catch dangling reference!
#else
  return arcp(a.getRawPtr(), 0, a.size(), false);
#endif
}


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


/** \brief Return the hash code.
 *
 * \relates Array.
 */
template<typename T> inline
int hashCode(const Array<T>& array);


/** \brief Copy conversion to an std::vector.
 *
 * This function is included for consistency with ArrayView.
 *
 * \relates Array.
 */
template<typename T> inline
std::vector<T> createVector( const Array<T> &a );


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

/** \brief A wrapper around the \c fromStringToArray function
 * which allows the operator>> to be used on Arrays.
 *
 * \relates Array
 */
template<typename T>
std::istringstream& operator>> (std::istringstream& in, Array<T>& array){
  array = fromStringToArray<T>(in.str());
  return in;
}

/** \brief Extracts data from an istringstream object
 * \note This templated function is necessary for the proper extraction of
 *       data by the \c fromStringToArray function.
 * \relates Array.
 */
template<typename T> inline
void extractDataFromISS( std::istringstream& iss, T& data )
{
  iss >> data; // Assumes type has operator>>(...) defined!
}

/** \brief Extracts std::string data from an istringstream object
 * \note This function overloads the templated \c extractDataFromISS function
         and is necessary for the proper extraction of std::string objects
         by the \c fromStringToArray function.
 * \relates Array.
 */
inline
void extractDataFromISS( std::istringstream& iss, std::string& data )
{
  // grab unformatted string.
  data = iss.str();
  // remove white space from beginning and end of string.
  data = Utils::trimWhiteSpace(data);
}

/**
 * \brief Get the format that is used for the specialization of the TypeName
 * traits class for Array.
 *
 * The string returned will contain only one
 * "*" character. The "*" character should then be replaced with the actual
 * template type of the array.
 * \relates Array.
 */
inline
std::string getArrayTypeNameTraitsFormat(){
  return "Array(*)";
}



/** \brief TypeNameTraits specialization for Array.
 *
 * NOTE: Use of this class requires that either that the type T be defined or
 * that a TypeNameTraits<T> specialization exists.  In order to not restrict
 * the use of Array<T> for undefined pointer types (where T=U*), this
 * TypeNameTraits class specialization will not be used in core Array
 * functionality.  This might seem trivial except that some MPI
 * implementations use pointers to undefined structs and if you want to
 * portably story these undefined struct pointers in an Array, then you can't
 * use this traits class.  C++ is quite lacking in cases like this.
 *
 * \ingroup teuchos_mem_mng_grp
 */
template<typename T>
class TEUCHOS_LIB_DLL_EXPORT TypeNameTraits<Array<T> > {
public:
  static std::string name(){
    std::string formatString = getArrayTypeNameTraitsFormat();
    size_t starPos = formatString.find("*");
    std::string prefix = formatString.substr(0,starPos);
    std::string postFix = formatString.substr(starPos+1);
    return prefix+TypeNameTraits<T>::name()+postFix;
  }
  static std::string concreteName(const Array<T>&)
    { return name(); }
};


} // namespace Teuchos


//
// Implementation
//


namespace Teuchos {


// All constructors


template<typename T> inline
Array<T>::Array()
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  : vec_(rcp(new std::vector<T>()))
#endif
{}


template<typename T> inline
Array<T>::Array(size_type n, const value_type& value) :
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  vec_(rcp(new std::vector<T>(n,value)))
#else
  vec_(n, value)
#endif
{}


template<typename T> inline
Array<T>::Array(const Array<T>& x) :
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  vec_(rcp(new std::vector<T>(*x.vec_)))
#else
  vec_(x.vec_)
#endif
{}


template<typename T> template<typename InputIterator> inline
Array<T>::Array(InputIterator first, InputIterator last) :
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  vec_(rcp(new std::vector<T>(first, last)))
#else
  vec_(first, last)
#endif
{}


template<typename T> inline
Array<T>::~Array()
{}


template<typename T> inline
Array<T>::Array(const ArrayView<const T>& a)
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  : vec_(rcp(new std::vector<T>()))
#endif
{
  insert(begin(), a.begin(), a.end());
}


template<typename T>
template<int N>
inline
Array<T>::Array(const Tuple<T,N>& t)
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  : vec_(rcp(new std::vector<T>()))
#endif
{
  insert(begin(), t.begin(), t.end());
}


template<typename T> inline
Array<T>& Array<T>::operator=(const Array& a)
{
  vec(true) = a.vec();
  return *this;
}


// Other std::vector functions


template<typename T> inline
void Array<T>::assign(size_type n, const value_type& val)
{
  vec(true).assign(n,val);
}


template<typename T> template<typename InputIterator> inline
void Array<T>::assign(InputIterator first, InputIterator last)
{
  vec(true).assign(first,last);
}


template<typename T> inline
typename Array<T>::iterator
Array<T>::begin()
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  if (is_null(extern_arcp_)) {
    // Here we must use the same RCP to avoid creating two unrelated RCPNodes!
    extern_arcp_ = arcp(vec_); // Will be null if vec_ is sized!
  }
  // Returning a weak pointer will help to catch dangling references but still
  // keep the same behavior as optimized code.
  return extern_arcp_.create_weak();
#else
  return vec().begin();
#endif
}


template<typename T> inline
typename Array<T>::iterator
Array<T>::end()
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return begin() + size();
#else
  return vec().end();
#endif
}


template<typename T> inline
typename Array<T>::const_iterator
Array<T>::begin() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  if (is_null(extern_carcp_)) {
    extern_carcp_ = const_cast<Array<T>*>(this)->begin();
  }
  // Returning a weak pointer will help to catch dangling references but still
  // keep the same behavior as optimized code.
  return extern_carcp_.create_weak();
#else
  return vec().begin();
#endif
}


template<typename T> inline
typename Array<T>::const_iterator
Array<T>::end() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return begin() + size();
#else
  return vec().end();
#endif
}


template<typename T> inline
typename Array<T>::reverse_iterator
Array<T>::rbegin()
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return reverse_iterator(end());
#else
  return vec().rbegin();
#endif
}


template<typename T> inline
typename Array<T>::reverse_iterator
Array<T>::rend()
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return reverse_iterator(begin());
#else
  return vec().rend();
#endif
}


template<typename T> inline
typename Array<T>::const_reverse_iterator
Array<T>::rbegin() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return const_reverse_iterator(end());
#else
  return vec().rbegin();
#endif
}


template<typename T> inline
typename Array<T>::const_reverse_iterator
Array<T>::rend() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return const_reverse_iterator(begin());
#else
  return vec().rend();
#endif
}


template<typename T> inline
typename Array<T>::size_type
Array<T>::size() const
{
  return vec().size();
}


template<typename T> inline
typename Array<T>::size_type
Array<T>::max_size() const
{
  return std::numeric_limits<size_type>::max();
}


template<typename T> inline
void
Array<T>::resize(size_type new_size, const value_type& x)
{
  vec(true).resize(new_size,x);
}


template<typename T> inline
typename Array<T>::size_type
Array<T>::capacity() const
{
  return vec().capacity();
}


template<typename T> inline
bool Array<T>::empty() const
{
  return vec().empty();
}


template<typename T> inline
void Array<T>::reserve(size_type n)
{
  vec(true).reserve(n);
}


template<typename T> inline
typename Array<T>::reference
Array<T>::operator[](size_type i)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assertIndex(i);
#endif
  return vec()[i];
}


template<typename T> inline
typename Array<T>::const_reference
Array<T>::operator[](size_type i) const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assertIndex(i);
#endif
  return vec()[i];
}


template<typename T> inline
typename Array<T>::reference
Array<T>::at(size_type i)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assertIndex(i);
#endif
  return vec().at(i);
}


template<typename T> inline
typename Array<T>::const_reference
Array<T>::at(size_type i) const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assertIndex(i);
#endif
  return vec().at(i);
}


template<typename T> inline
typename Array<T>::reference
Array<T>::front()
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assertNotNull();
#endif
  return vec().front();
}


template<typename T> inline
typename Array<T>::const_reference
Array<T>::front() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assertNotNull();
#endif
  return vec().front();
}


template<typename T> inline
typename Array<T>::reference
Array<T>::back()
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assertNotNull();
#endif
  return vec().back();
}


template<typename T> inline
typename Array<T>::const_reference
Array<T>::back() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assertNotNull();
#endif
  return vec().back();
}


template<typename T> inline
void Array<T>::push_back(const value_type& x)
{
  vec(true).push_back(x);
}


template<typename T> inline
void Array<T>::pop_back()
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assertNotNull();
#endif
  vec(true).pop_back();
}


// 2009/11/13:: rabartl: After moving to a full RCPNode tracing and lookup
// model, I had to how modifying functions like insert(...) and erase(...)
// work which have active iterators controled by the client and yet need to
// allow the structure of the container change.  The way these troublesome
// functions work is that first the raw std::vector iterator is extracted.
// The function vec(true, true) then deletes the strong iterators but there is
// still a weak ArrayRCP object that is owned by the client which is being
// passed into this function.  The issue is that the design of ArrayRCP is
// such that the RCPNode object is not removed but instead remains in order to
// perform runtime checking.


template<typename T> inline
typename Array<T>::iterator
Array<T>::insert(iterator position, const value_type& x)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  // Assert a valid iterator and get vector iterator
  const typename std::vector<T>::iterator raw_poss = raw_position(position);
  const difference_type i = position - begin();
  vec(true, true).insert(raw_poss, x);
  return begin() + i;
#else
  return vec_.insert(position, x);
#endif
}


template<typename T> inline
void Array<T>::insert(iterator position, size_type n, const value_type& x)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  const typename std::vector<T>::iterator raw_poss = raw_position(position);
  vec(true, true).insert(raw_poss, n, x);
#else
  return vec_.insert(position, n, x);
#endif
}


template<typename T> template<typename InputIterator> inline
void Array<T>::insert(iterator position, InputIterator first, InputIterator last)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  const typename std::vector<T>::iterator raw_poss = raw_position(position);
  vec(true, true).insert(raw_poss, first, last);
#else
  return vec_.insert(position, first, last);
#endif
}


template<typename T> inline
typename Array<T>::iterator
Array<T>::erase(iterator position)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assertNotNull();
  // Assert a valid iterator and get vector iterator
  const typename std::vector<T>::iterator raw_poss = raw_position(position);
  const difference_type i = position - begin();
  vec(true, true).erase(raw_poss);
  return begin() + i;
#else
  return vec_.erase(position);
#endif
}


template<typename T> inline
typename Array<T>::iterator
Array<T>::erase(iterator first, iterator last)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  if (empty()) {
    TEUCHOS_ASSERT(first == begin());
    TEUCHOS_ASSERT(last == end());
    return end();
  }
  assertNotNull();
  // Assert a valid iterator and get vector iterator
  const typename std::vector<T>::iterator raw_first = raw_position(first);
  const typename std::vector<T>::iterator raw_last = raw_position(last);
  const difference_type i = first - begin();
  vec(true,true).erase(raw_first,raw_last);
  return begin() + i;
#else
  return vec_.erase(first,last);
#endif
}


template<typename T> inline
void Array<T>::swap(Array& x)
{
  vec(true).swap(x.vec());
}


template<typename T> inline
void Array<T>::clear()
{
  vec(true).clear();
}


// Non-standard functions


template<typename T> inline
Array<T>& Array<T>::append(const T& x)
{
  this->push_back(x);
  return *this;
}


template<typename T> inline
void Array<T>::remove(int i)
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  assertIndex(i);
#endif
  // Erase the i-th element of this array.
  this->erase( this->begin() + i );
}


template<typename T> inline
int Array<T>::length() const
{
  return this->size();
}


template<typename T> inline
std::string Array<T>::toString() const
{
  return (*this)().toString(); // Use ArrayView<T>::toString()
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
T* Array<T>::getRawPtr()
{
  return ( size() ? &(*this)[0] : 0 );
}


template<typename T> inline
const T* Array<T>::getRawPtr() const
{
  return ( size() ? &(*this)[0] : 0 );
}


// Conversions to and from std::vector


template<typename T> inline
Array<T>::Array( const std::vector<T> &v ) :
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  vec_(new std::vector<T>(v))
#else
  vec_(v)
#endif
{}


template<typename T> inline
std::vector<T> Array<T>::toVector() const
{
  if (!size())
    return std::vector<T>();
  std::vector<T> v(begin(),end());
  return v;
}


template<typename T> inline
Array<T>& Array<T>::operator=( const std::vector<T> &v )
{
  vec(true) = v;
  return *this;
}


// Views


template<typename T> inline
ArrayView<T> Array<T>::view( size_type offset, size_type size_in )
{
  if (size_in) {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    return ArrayView<T>(this->begin().persistingView(offset, size_in));
#else
    return arrayView( &vec()[offset], size_in );
#endif
  }
  return Teuchos::null;
}


template<typename T> inline
ArrayView<const T> Array<T>::view( size_type offset, size_type size_in ) const
{
  if (size_in) {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    return ArrayView<const T>(this->begin().persistingView(offset, size_in));
#else
    return arrayView( &vec()[offset], size_in );
#endif
  }
  return Teuchos::null;
  // NOTE: Above, we use a different implementation to call the const version
  // of begin() instead of the non-const version.  This sets up a different
  // ArrayRCP object that gets checked.
}


template<typename T> inline
ArrayView<T> Array<T>::operator()( size_type offset, size_type size_in )
{
  return view(offset, size_in);
}


template<typename T> inline
ArrayView<const T> Array<T>::operator()( size_type offset, size_type size_in ) const
{
  return view(offset, size_in);
}


template<typename T> inline
ArrayView<T> Array<T>::operator()()
{
  if (!size())
    return null;
  return this->view(0, size());
}


template<typename T> inline
ArrayView<const T> Array<T>::operator()() const
{
  if (!size())
    return null;
  return this->view(0, size());
}


template<typename T> inline
Array<T>::operator ArrayView<T>()
{
  return this->operator()();
}


template<typename T> inline
Array<T>::operator ArrayView<const T>() const
{
  return this->operator()();
}


// private


template<typename T>
std::vector<T>&
Array<T>::vec( bool isStructureBeingModified, bool activeIter )
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  (void)activeIter;
  if (isStructureBeingModified) {
    // Give up my ArrayRCPs used for iterator access since the array we be
    // getting modifed!  Any clients that have views through weak pointers
    // better not touch them!
    extern_arcp_ = null;
    extern_carcp_ = null;
  }
  return *vec_;
#else
  // get rid of "unused parameter" warnings
  (void)isStructureBeingModified;
  (void)activeIter;
  return vec_;
#endif
}


template<typename T> inline
const std::vector<T>&
Array<T>::vec() const
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  return *vec_;
#else
  return vec_;
#endif
}


template<typename T> inline
typename std::vector<T>::iterator
Array<T>::raw_position( iterator position )
{
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  const iterator first = this->begin();
  const iterator last = this->end();
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(first <= position && position <= last), DanglingReferenceError,
    "Error, this iterator is no longer valid for this Aray!"
    );
  // Note, above operator<=(...) functions will throw
  // IncompatibleIteratorsError if the iterators do not share the same
  // RCP_node object!
  return vec_->begin() + (position - this->begin());
#else
  return position;
#endif
}


template<typename T> inline
void Array<T>::assertIndex(int i) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !( 0 <= i && i < length() ), RangeError,
    "Array<T>::assertIndex(i): i="<<i<<" out of range [0, "<< length() << ")"
    );
}


template<typename T> inline
void Array<T>::assertNotNull() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !size(), NullReferenceError,
    typeName(*this)<<"::assertNotNull(): "
    "Error, the array has size zero!"
    );
}


} // namespace Teuchos


// Nonmember functions


template<typename T> inline
bool Teuchos::operator==( const Array<T> &a1, const Array<T> &a2 )
{ return (a1.vec() == a2.vec()); }


template<typename T> inline
bool Teuchos::operator!=( const Array<T> &a1, const Array<T> &a2 )
{ return (a1.vec() != a2.vec()); }


template<typename T> inline
void Teuchos::swap( Array<T> &a1, Array<T> &a2 )
{ a1.swap(a2); }


template<typename T> inline
bool Teuchos::operator<( const Array<T> &a1, const Array<T> &a2 )
{ return (a1.vec() < a2.vec()); }


template<typename T> inline
bool Teuchos::operator<=( const Array<T> &a1, const Array<T> &a2 )
{ return (a1.vec() <= a2.vec()); }


template<typename T> inline
bool Teuchos::operator>( const Array<T> &a1, const Array<T> &a2 )
{ return (a1.vec() > a2.vec()); }


template<typename T> inline
bool Teuchos::operator>=( const Array<T> &a1, const Array<T> &a2 )
{ return (a1.vec() >= a2.vec()); }


template<typename T> inline
std::ostream& Teuchos::operator<<(
  std::ostream& os, const Array<T>& array
  )
{
  return os << Teuchos::toString(array);
}


template<typename T> inline
int Teuchos::hashCode(const Array<T>& array)
{
  int rtn = hashCode(array.length());
  for (int i=0; i<array.length(); i++)
  {
    rtn += hashCode(array[i]);
  }
  return rtn;
}


template<typename T> inline
std::vector<T> Teuchos::createVector( const Array<T> &a )
{
  return a.toVector();
}


template<typename T> inline
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
  TEUCHOS_TEST_FOR_EXCEPTION(
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
  TEUCHOS_TEST_FOR_EXCEPT(c!='{'); // Should not throw!
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
    TEUCHOS_TEST_FOR_EXCEPTION(
      0 == entryStr.length(),
      InvalidArrayStringRepresentation,
      "Error, the std::string:\n"
      "----------\n"
      <<str<<
      "\n----------\n"
      "is not a valid array represntation because it has an empty array entry!"
      );
    // Remove the final '}' if this is the last entry and we did not
    // actually terminate the above getline(...) on ','
    bool found_end = false;
    if(entryStr[entryStr.length()-1]=='}') {
      entryStr = entryStr.substr(0,entryStr.length()-1);
      found_end = true;
      if( entryStr.length()==0 && a.size()==0 )
        return a; // This is the empty array "{}" (with any spaces in it!)
    }
    // Finally we can convert the entry and add it to the array!
    std::istringstream entryiss(entryStr);
    T entry;
    Teuchos::extractDataFromISS( entryiss, entry );
    // ToDo: We may need to define a traits class to allow us to specialized
    // how conversion from a std::string to a object is done!
    a.push_back(entry);
    // At the end of the loop body here, if we have reached the last '}'
    // then the input stream iss should be empty and iss.eof() should be
    // true, so the loop should terminate.  We put an std::exception test here
    // just in case something has gone wrong.
    TEUCHOS_TEST_FOR_EXCEPTION(
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


#endif // TEUCHOS_ARRAY_H
