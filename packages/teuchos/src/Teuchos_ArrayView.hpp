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

#ifndef TEUCHOS_ARRAY_VIEW_HPP
#define TEUCHOS_ARRAY_VIEW_HPP


#include "Teuchos_ENull.hpp"


namespace Teuchos {


template<class T> class ArrayRCP;


/** \brief Array view class.
 *
 * This class is designed to be used as a substitute for array arguments to
 * functions.  It aggregates a pointer to a contiguous array of data and the
 * size of that array.  In debug mode, it will perform runtime checks of all
 * usage.
 *
 * The <tt>ArrayView</tt> class has bind on construction semantics where once
 * an <tt>ArrayView</tt> object is bound to an array of data, it can not be
 * bound to a different view of data.  Therefore, a <tt>const ArrayView</tt>
 * object means that the data entries are <tt>const</tt> which is similar to
 * the behavior of a <tt>const std::vector</tt> object for instance.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup teuchos_mem_mng_grp
 */
template<class T>
class ArrayView {
public:

  /** \name std::vector typedefs */
  //@{

  /** \brief. */
  typedef Teuchos_Index Ordinal;

  /** \brief . */
  typedef T value_type;
  /** \brief . */
  typedef T* pointer;
  /** \brief . */
  typedef const T* const_pointer;
  /** \brief . */
  typedef T& reference;
  /** \brief . */
  typedef const T& const_reference;

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  /** \brief . */
  typedef ArrayRCP<T> iterator;
  /** \brief . */
  typedef ArrayRCP<const T> const_iterator;
#else
  /** \brief . */
  typedef pointer iterator;
  /** \brief . */
  typedef const_pointer const_iterator;
#endif

  /** \brief . */
  typedef std::reverse_iterator<iterator> reverse_iterator;
  /** \brief . */
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;


  /** \brief . */
  typedef typename std::vector<T>::size_type size_type;
  /** \brief . */
  typedef typename std::vector<T>::difference_type difference_type;
  /** \brief . */
  typedef typename std::vector<T>::allocator_type allocator_type;

  //@}

  //! @name Constructors/Initializers 
  //@{

	/** \brief Initialize <tt>ArrayView<T></tt> to NULL.
	 *
	 * This allows clients to write code like:
	 \code
	 ArrayView<int> p = null;
	 \endcode
	 * or
	 \code
	 ArrayView<int> p;
	 \endcode
	 * and construct to <tt>NULL</tt>
	 */
	ArrayView( ENull null_arg = null );

	/** \brief Initialize from another <tt>ArrayView<T></tt> object.
	 *
	 * After construction, <tt>this</tt> and <tt>array</tt> will reference the
	 * same array.
	 *
	 * This form of the copy constructor is required even though the
	 * below more general templated version is sufficient since some
	 * compilers will generate this function automatically which will
	 * give an incorrect implementation.
	 *
	 * Postconditions:<ul>
	 * <li>???
	 * </ul>
	 */
	ArrayView(const ArrayView<T>& array);

	/** \brief Destroy the array view object.
	 */
	~ArrayView();

	/** \brief Copy the data from one array view object to this array view
   * object.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->size() == array.size()</tt>
	 * </ul>
	 */
	ArrayView<T>& operator=(const ArrayView<T>& array);

  //@}

  //! @name General query functions 
  //@{

  /** \brief The total number of items in the managed array. */
  Ordinal size() const;

  //@}

  //! @name Element Access Functions 
  //@{

	/** \brief Random object access.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->get() != NULL</tt>
   * <li><tt>0 <= offset && offset < this->size()</tt>
	 * </ul>
   */
	T& operator[](Ordinal offset) const;

  //@}

  //! @name Views 
  //@{

	/** \brief Return a non-const view of a contiguous range of elements.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->get() != NULL</tt>
   * <li><tt>0 <= lowerOffset</tt>
   * <li><tt>lowerOffset + size < this->size()</tt>
	 * </ul>
	 *
	 * <b>Postconditions:</b><ul>
   * <li><tt>returnVal.size() == size</tt>
	 * </ul>
   */
	ArrayView<T> operator()( Ordinal lowerOffset, Ordinal size );

	/** \brief Return a const view of a contiguous range of elements.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->get() != NULL</tt>
   * <li><tt>0 <= lowerOffset</tt>
   * <li><tt>lowerOffset + size < this->size()</tt>
	 * </ul>
	 *
	 * <b>Postconditions:</b><ul>
   * <li><tt>returnVal.size() == size</tt>
	 * </ul>
   */
	const ArrayView<T> operator()( Ordinal lowerOffset, Ordinal size ) const;

  //@}

  //! @name Standard Container-Like Functions 
  //@{

  /** \brief Return an iterator to beginning of the array of data.
   *
   * If <tt>HAVE_TEUCHOS_ARRAY_BOUNDSCHECK</tt> is defined then the iterator
   * returned is an <tt>ArrayView<T></tt> object and all operations are
   * checked at runtime.  When <tt>HAVE_TEUCHOS_ARRAY_BOUNDSCHECK</tt> is not
   * defined, the a raw pointer <tt>T*</tt> is returned for fast execution.
   *
   * <b>Postconditions:</b><ul>
   * <li>[this->get()!=NULL</tt>] <tt>&*return == this->get()</tt>
   * <li>[<tt>this->get()==NULL</tt>] <tt>return == (null or NULL)</tt>
   * </ul>
   */
  const_iterator begin() const;

  /** \brief Return an iterator to past the end of the array of data.
   *
   * If <tt>HAVE_TEUCHOS_ARRAY_BOUNDSCHECK</tt> is defined then the iterator
   * returned is an <tt>ArrayView<T></tt> object and all operations are
   * checked at runtime.  When <tt>HAVE_TEUCHOS_ARRAY_BOUNDSCHECK</tt> is not
   * defined, the a raw pointer <tt>T*</tt> is returned for fast execution.
   *
   * <b>Postconditions:</b><ul>
   * <li>[<tt>this->get()!=NULL</tt>] <tt>&*end == this->get()+(this->upperOffset()+1)</tt>
   * <li>[<tt>this->get()==NULL</tt>] <tt>return == (null or NULL)</tt>
   * </ul>
   */
  const_iterator end() const;

  //@}

  //! @name Assertion Functions. 
  //@{

	/** \brief Throws <tt>std::logic_error</tt> if <tt>this->size()==0</tt>,
   * otherwise returns reference to <tt>*this</tt>.
   */
	const ArrayView<T>& assert_not_null() const;

	/** \brief Throws <tt>std::logic_error</tt> if <tt>lowerOffset < 0 ||
   * this->size() <= lowerOffset+size</tt>.
   */
	const ArrayView<T>& assert_in_range( Ordinal lowerOffset, Ordinal size ) const;

  //@}

private:

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  ArrayRCP<T> aptr_;
#else
	T *ptr_;
  int size_;
#endif

};


/** \brief Construct a array to non-const data .
 *
 * \relates ArrayView
 */
template<class T>
ArrayView<T> arrayView( T* p, typename ArrayView<T>::Ordinal size );


/** \brief Construct a array to const data .
 *
 * \relates ArrayView
 */
template<class T>
const ArrayView<T> arrayView( const T* p, typename ArrayView<T>::Ordinal size );


/** \brief Get a new <tt>std::vector<T></tt> object out of an
 * <tt>ArrayView<T></tt> object.
 *
 * \relates ArrayView
 */
template<class T>
std::vector<T> create_std_vector( ArrayView<T> &ptr );


/** \brief Get a new <tt>std::vector<T></tt> object out of an
 * <tt>ArrayView<T></tt> object.
 *
 * \relates ArrayView
 */
template<class T>
const std::vector<T> create_std_vector( const ArrayView<T> &ptr );


/** \brief Output stream inserter.
 *
 * The implementation of this function just prints pointer addresses and
 * therefore puts no restrictions on the data types involved.
 *
 * \relates ArrayView
 */
template<class T>
std::ostream& operator<<( std::ostream& out, const ArrayView<T>& p );


} // end namespace Teuchos


// ////////////////////////////////
// Implementations


// ToDo: Add implementations!


#endif	// TEUCHOS_ARRAY_VIEW_HPP
