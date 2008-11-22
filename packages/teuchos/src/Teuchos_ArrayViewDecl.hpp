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

#ifndef TEUCHOS_ARRAY_VIEW_DECL_HPP
#define TEUCHOS_ARRAY_VIEW_DECL_HPP


#include "Teuchos_RCPNode.hpp"
#include "Teuchos_ENull.hpp"
#include "Teuchos_NullIteratorTraits.hpp"
#include "Teuchos_ConstTypeTraits.hpp"


namespace Teuchos {


template<class T> class ArrayRCP;


/** \brief Array view class.
 *
 * This class is designed to be used as a substitute for array arguments to
 * functions.  It aggregates a pointer to a contiguous array of data and the
 * size of that array.  In debug mode, it will perform runtime checks of all
 * usage.
 *
 * The <tt>ArrayView</tt> class has the same shallow copy semantics of the #
 * <tt>Ptr</tt> class.  <tt>ArrayView</tt> is to <tt>ArrayRCP</tt> as
 * <tt>Ptr</tt> is to <tt>RCP</tt>.
 *
 * \section Teuchos_ArrayView_DesignDiscussion_sec Design Discussion
 *
 * This class is setup to allow derived classes that can only be allocated on
 * the stack.
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
#else
  /** \brief . */
  typedef pointer iterator;
#endif

  /** \brief . */
  typedef size_t size_type;
  /** \brief . */
  typedef ptrdiff_t difference_type;

  //@}

  //! @name Constructors/Destructors
  //@{

	/** \brief Initialize to NULL (implicitly or explicitly). */
	ArrayView( ENull null_arg = null );

	/** \brief Initialize view from raw memory.
	 *
   * \param p [in] Pointer to array of typed memory of size <tt>size</tt>.  If
   * <tt>p==0</tt>, then <tt>*this</tt> is a null view.  Note that the memory
   * pointed to by <tt>p</tt> can not go away until this view object is
   * destoryed!
   *
   * \param size [in] The size of the array that <tt>*this</tt> will represent
   * pointer to by <tt>p</tt>.  If <tt>p==0</tt> then <tt>size</tt> must be 0!
   *
	 * Preconditions:<ul>
	 * <li>[<tt>p!=0</tt>] <tt>size > 0</tt>
	 * <li>[<tt>p==0</tt>] <tt>size == 0</tt>
	 * </ul>
   *
	 * Postconditions:<ul>
	 * <li>???
	 * </ul>
	 */
	ArrayView( T* p, Ordinal size );

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

  /** \brief Non-const view of an std::vector<T> .*/
	ArrayView(
    std::vector<typename ConstTypeTraits<T>::NonConstType>& vec
    );

  /** \brief Const view of an std::vector<T> .*/
  ArrayView(
    const std::vector<typename ConstTypeTraits<T>::NonConstType>& vec
    );
  
  /** \brief Shallow copy assignment operator. */
	ArrayView<T>& operator=(const ArrayView<T>& array);

	/** \brief Destroy the array view object.
	 */
	~ArrayView();

  //@}

  //! @name General query functions 
  //@{

  /** \brief Returns true if the underlying pointer is null. */
  bool is_null() const;

  /** \brief The total number of items in the managed array. */
  Ordinal size() const;

  /** \brief Convert an ArrayView<T> to an <tt>std::string</tt> */
  std::string toString() const;

  //@}

  //! @name Element Access Functions 
  //@{

  /** \brief Return a raw pointer to beginning of array or NULL if unsized. */
  inline T* getRawPtr() const;

	/** \brief Random object access.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->get() != NULL</tt>
   * <li><tt>0 <= offset && offset < this->size()</tt>
	 * </ul>
   */
	T& operator[](Ordinal i) const;

  /** \brief Get the first element. */
  T& front() const;

  /** \brief Get the last element. */
  T& back() const;

  //@}

  //! @name Views 
  //@{

	/** \brief Return a view of a contiguous range of elements.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->get() != NULL</tt>
   * <li><tt>0 <= offset && offset + size <= this->size()</tt>
	 * </ul>
	 *
	 * <b>Postconditions:</b><ul>
   * <li><tt>returnVal.size() == size</tt>
	 * </ul>
   */
	ArrayView<T> view( Ordinal offset, Ordinal size ) const;

	/** \brief Return a view of a contiguous range of elements (calls view(offset,size)).
   */
	ArrayView<T> operator()( Ordinal offset, Ordinal size ) const;

	/** \brief Return a *this (just for compatibility with Array and ArrayPtr)
   */
	const ArrayView<T>& operator()() const;

  /** \brief Return an ArrayView<const T> of an ArrayView<T> object.
   *
   * WARNING!  If <tt>T</tt> is already const (e.g. <tt>const double</tt>)
   * then do not try to instantiate this function since it will not compile!
   */
  ArrayView<const T> getConst() const;

  /** \brief Impliict conversion from ArrayView<T> to ArrayView<const T>.
   *
   * WARNING!  If <tt>T</tt> is already const (e.g. <tt>const double</tt>)
   * then do not try to instantiate this function since it will not compile!
   */
  operator ArrayView<const T>() const;

  //@}

  /** \name Assignment */
  //@{

	/** \brief Copy the data from one array view object to this array view
   * object.
	 *
	 * <b>Preconditions:</b><ul>
	 * <li><tt>this->size() == array.size()</tt>
	 * </ul>
   *
   * WARNING!  If <tt>T</tt> is a const type (e.g. <tt>const double</tt>) then
   * do not try to instantiate this function since it will not compile!
   *
   * NOTE: This function is really like an operator=() function except that it
   * is declared const.  This is the correct behavior since a const ArrayView
   * simply means that we can not change what *this points to.  The type of
   * the template argument always determines if the underlyihng data is const
   * or not.
   */
	void assign(const ArrayView<const T>& array) const;

  //@}

  //! @name Standard Container-Like Functions 
  //@{

  /** \brief Return an iterator to beginning of the array of data.
   *
   * If <tt>HAVE_TEUCHOS_ARRAY_BOUNDSCHECK</tt> is defined then the iterator
   * returned is an <tt>ArrayRCP<T></tt> object and all operations are
   * checked at runtime.  When <tt>HAVE_TEUCHOS_ARRAY_BOUNDSCHECK</tt> is not
   * defined, the a raw pointer <tt>T*</tt> is returned for fast execution.
   *
   * <b>Postconditions:</b><ul>
   * <li>[this->get()!=NULL</tt>] <tt>&*return == this->get()</tt>
   * <li>[<tt>this->get()==NULL</tt>] <tt>return == (null or NULL)</tt>
   * </ul>
   */
  iterator begin() const;

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
  iterator end() const;

  //@}

  //! @name Assertion Functions. 
  //@{

	/** \brief Throws <tt>NullReferenceError</tt> if <tt>this->get()==NULL</tt>,
   * otherwise returns reference to <tt>*this</tt>.
   */
	const ArrayView<T>& assert_not_null() const;

	/** \brief Throws <tt>NullReferenceError</tt> if <tt>this->get()==NULL</tt>
   * or<tt>this->get()!=NULL</tt>, throws <tt>RangeError</tt> if <tt>(offset < 0 ||
   * this->size() < offset+size</tt>, otherwise returns reference to <tt>*this</tt>
   */
	const ArrayView<T>& assert_in_range( Ordinal offset, Ordinal size ) const;

  //@}

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  // I should make make this public but templated friends are
  // not very portable.
	explicit ArrayView( const ArrayRCP<T> &arcp );
#endif

private:

  // ///////////////////////
  // Private data members

	T *ptr_;
  int size_;
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  ArrayRCP<T> arcp_;
#endif

  void setUpIterators();

  // ///////////////////////
  // Private member functions

  void debug_assert_not_null() const
    {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
      assert_not_null();
#endif
    }

  void debug_assert_in_range( Ordinal offset, Ordinal size_in ) const
    {
      (void)offset; (void)size_in;
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
      assert_in_range(offset, size_in);
#endif
    }

  void debug_assert_valid_ptr() const
    {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
      arcp_.access_private_node().assert_valid_ptr(*this);
#endif
    }

public: // Bad bad bad

  // This is a very bad breach of encapsulation but it exists to avoid
  // problems with portability of tempalted friends
  T* access_private_ptr() const
    { return ptr_; }


};


/** \brief Construct a const or non-const view to const or non-const data.
 *
 * \relates ArrayView
 */
template<class T>
ArrayView<T> arrayView( T* p, typename ArrayView<T>::Ordinal size );


/** \brief Construct a non-const view of an std::vector.
 *
 * \relates ArrayView
 */
template<class T>
ArrayView<T> arrayViewFromVector( std::vector<T>& vec );


/** \brief Construct a const view of an std::vector.
 *
 * \relates ArrayView
 */
template<class T>
ArrayView<const T> arrayViewFromVector( const std::vector<T>& vec );


#ifndef __sun


// 2007/11/30: From some reason, the Sun C++ compile on sass9000 compains that
// a call to this function below is ambiguous.  However, if you just comment
// the function out, then the code on g++ (3.4.6 at least) will not compile.
// Therefore, I have no choice but to put in a hacked ifdef for the sun.


/** \brief Get a new <tt>std::vector<T></tt> object out of an
 * <tt>ArrayView<T></tt> object.
 *
 * Note that a copy of data is made!
 *
 * \relates ArrayView
 */
template<class T>
std::vector<T> createVector( const ArrayView<T> &av );


#endif // __sun


/** \brief Get a new <tt>std::vector<T></tt> object out of an
 * <tt>ArrayView<const T></tt> object.
 *
 * Note that a copy of data is made!
 *
 * \relates ArrayView
 */
template<class T>
std::vector<T> createVector( const ArrayView<const T> &av );


/** \brief Returns true if <tt>av.is_null()==true</tt>.
 *
 * \relates ArrayView
 */
template<class T>
bool is_null( const ArrayView<T> &av );


/** \brief Output stream inserter.
 *
 * The implementation of this function just prints pointer addresses and
 * therefore puts no restrictions on the data types involved.
 *
 * \relates ArrayView
 */
template<class T>
std::ostream& operator<<( std::ostream& out, const ArrayView<T>& av );


} // end namespace Teuchos


//
// Inline members
//


// ToDo: Fill in!


#endif	// TEUCHOS_ARRAY_VIEW_DECL_HPP
