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


#ifndef TEUCHOS_ARRAY_RCP_DECL_HPP
#define TEUCHOS_ARRAY_RCP_DECL_HPP


#include "Teuchos_RCP.hpp"
#include "Teuchos_Exceptions.hpp"
#include "Teuchos_ArrayViewDecl.hpp"


namespace Teuchos {

/** \brief Reference-counted smart pointer for managing arrays.
 * \tparam T The type of each element in the array.
 * \ingroup teuchos_mem_mng_grp
 *
 * \section Teuchos_ArrayRCP_Summary Summary
 *
 * ArrayRCP manages an array of objects of type T.  Like RCP, it uses
 * reference counting to decide when to deallocate the array.  This
 * lets you share references to the array, without worrying about who
 * is responsible for deallocating it.  We may thus call an ArrayRCP a
 * "shared array."
 *
 * \section Teuchos_ArrayRCP_Details Details
 *
 * Managing an array of objects is very different from managing a
 * pointer to an individual, possibly polymorphic, object.  For
 * example, while implicit conversions from derived to base types can
 * be useful when dealing with pointers to single objects, they often
 * cause problems when working with arrays of objects.  Therefore,
 * this class contains those capabilities of raw pointers that are
 * good when working with arrays of objects, but excludes those that
 * are bad, such as implicit conversions from derived to base types.
 * If you really do want a shared array with implicit conversions from
 * derived to base types, you may use an <tt>ArrayRCP<RCP<T> ></tt>.
 *
 * \section Teuchos_ArrayRCP_Bounds Optional bounds checking
 *
 * You may enable bounds checking and other safety checks for this
 * class by setting the <tt>Teuchos_ENABLE_DEBUG:BOOL=ON</tt> CMake
 * option when configuring your Trilinos build.  This option is off by
 * default.  It incurs a significant performance penalty and so is not
 * recommended for production builds.  Bounds checking requires that
 * you always create ArrayRCP instances with the correct range.  For
 * example, if you use one of the constructors that accepts a raw
 * pointer, you are responsible for supplying the correct number of
 * elements in the array.  Our bounds checking implementation does not
 * attempt to replace memory debugging tools such as the Memcheck tool
 * in <a href="http://en.wikipedia.org/wiki/Valgrind">Valgrind</a>.
 *
 * \section Teuchos_ArrayRCP_Req Requirements on the type T
 *
 * ArrayRCP imposes the following requirements on the type T of
 * elements in the array:
 * <ul>
 * <li> T must be default constructible.
 * <li> T must be copy constructible.
 * <li> TypeNameTraits must have a specialization for T.
 * </ul>
 */
template<class T>
class ArrayRCP {
public:

  //! @name Public types 
  //@{

  //! Integer index type used throughout ArrayRCP.
  typedef Teuchos_Ordinal Ordinal;

  //! Type representing the number of elements in an ArrayRCP or view thereof.
  typedef Ordinal size_type;

  //! Type representing the difference between two iterators of an ArrayRCP.
  typedef Ordinal difference_type;

  //! Category of ArrayRCP's iterator type.
  typedef std::random_access_iterator_tag iterator_category;

  //! Type of an ArrayRCP's iterator.
  typedef  T* iterator_type;

  //! Type of each array element.
  typedef  T value_type;

  //! Type of a (nonconstant) reference to an array element.
  typedef T& reference; 

  //! Type of a (constant) reference to an array element.
  typedef const T& const_reference; 

  //! Type of a (raw) (nonconstant) pointer to an array element.
  typedef T* pointer;

  //! Type of a (raw) (constant) pointer to an array element.
  typedef T* const_pointer;

  //! Type of each array element.
  typedef T  element_type;

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  //! Nonconstant iterator type used if bounds checking is enabled.
  typedef ArrayRCP<T> iterator;
  //! Constant iterator type used if bounds checking is enabled.
  typedef ArrayRCP<const T> const_iterator;
#else
  //! Nonconstant iterator type used if bounds checking is disabled.
  typedef T* iterator;
  //! Constant iterator type used if bounds checking is disabled.
  typedef const T* const_iterator;
#endif

  //@}

  //! @name Constructors/Destructors/Initializers 
  //@{

  /** \brief Initialize <tt>ArrayRCP<T></tt> to NULL.
   *
   * This lets users write code like:
   \code
   ArrayRCP<int> p = null;
   \endcode
   * or
   \code
   ArrayRCP<int> p;
   \endcode
   * Both lines of code above set the raw array pointer to <tt>NULL</tt>.
   */
  inline ArrayRCP( ENull null_arg = null );

  /** \brief Construct from a raw pointer and a valid range.
   *
   * \param p [in] Raw array pointer.
   * \param lowerOffset [in] Array index at which the range starts
   *   (zero if at the beginning of the range).
   * \param size [in] Number of array elements in the range.
   * \param has_ownership [in] True if the ArrayRCP is responsible for
   *   deallocating the raw array (using <tt>delete []</tt>) when the
   *   reference count goes to zero.  If false, the ArrayRCP does not
   *   deallocate the array.
   * \param rcpNodeLookup [in] Whether to perform RCPNode lookup.  The
   *   default value is fine for nearly all use cases.
   *
   * \post <tt>this->get() == p</tt>
   * \post <tt>this->lowerOffset() == lowerOffset</tt>
   * \post <tt>this->upperOffset() == size + lowerOffset - 1</tt>
   * \post <tt>this->has_ownership() == has_ownership</tt>
   *
   * \warning You should avoid manipulating raw pointers and use other
   *   methods to construct an ArrayRCP object instead!
   */
  inline ArrayRCP( T* p, size_type lowerOffset, size_type size,
    bool has_ownership, const ERCPNodeLookup rcpNodeLookup = RCP_ENABLE_NODE_LOOKUP );

  /** \brief Construct from a raw pointer, a valid range, and a deallocator.
   *
   * \param p [in] Raw array pointer.
   * \param lowerOffset [in] Array index at which the range starts
   *   (zero if at the beginning of the range).
   * \param size [in] Number of array elements in the range.
   * \param dealloc [in] Function (or object with an
   *   <tt>operator()(T*)</tt> method) responsible for deallocating
   *   the raw array when the reference count goes to zero.
   * \param has_ownership [in] True if the ArrayRCP is responsible for
   *   deallocating the raw array (using the given deallocator) when
   *   the reference count goes to zero.  If false, the ArrayRCP does
   *   not deallocate the array.
   *
   * \post <tt>this->get() == p</tt>
   * \post <tt>this->lowerOffset() == lowerOffset</tt>
   * \post <tt>this->upperOffset() == size + lowerOffset - 1</tt>
   * \post <tt>this->has_ownership() == has_ownership</tt>
   *
   * \warning You should avoid manipulating raw pointers and use other
   *   methods to construct an ArrayRCP object instead!
   */
  template<class Dealloc_T>
  inline ArrayRCP( T* p, size_type lowerOffset, size_type size, Dealloc_T dealloc,
    bool has_ownership );

  /** \brief Construct an array with the given number of elements.
   *
   * \param size [in] Number of elements in the array.
   * \param val [in] Value with which to fill all elements of the
   *   array.  
   *
   * This constructor fills the array as if with the following code:
     \code
     std::fill_n (begin (), n, val);
     \endcode
   *
   * \post <tt>this->lowerOffset() == 0</tt>
   * \post <tt>this->upperOffset() == size - 1</tt>
   * \post <tt>this->has_ownership() == true</tt>
   */
  inline explicit ArrayRCP( size_type size, const T& val = T() );

  /** \brief Initialize from another <tt>ArrayRCP<T></tt> object.
   *
   * After construction, <tt>this</tt> and <tt>r_ptr</tt> will
   * reference the same array.
   *
   * \post <tt>this->get() == r_ptr.get()</tt>
   * \post <tt>this->count() == r_ptr.count()</tt>
   * \post <tt>this->has_ownership() == r_ptr.has_ownership()</tt>
   * \post  If <tt>r_ptr.get() != NULL</tt> then <tt>r_ptr.count()</tt> is incremented by 1.
   *
   * \note To implementers: In compilers that conform to the C++
   *   standard, this copy constructor overload is unnecessary, since
   *   the more general templated version of the copy constructor
   *   below is sufficient.  However, not all compilers have handled
   *   this code correctly in the past.  This version ensures correct
   *   compilation with such compilers, without affecting compilers
   *   that correctly implement the C++ standard.
   */
  inline ArrayRCP(const ArrayRCP<T>& r_ptr);

  /** \brief Destructor, that decrements the reference count.
   *
   * If <tt>this->get() == NULL</tt> then the destructor does nothing.
   * Otherwise, it decrements the reference count of this and all
   * other references to the array.  If the final reference count is
   * zero, it also deallocates the array if owned (i.e., if
   * <tt>this->has_ownership()</tt> returns true).  Deallocation uses
   * the custom deallocator if one was supplied; otherwise it uses
   * <tt>delete []</tt>.
   */
  inline ~ArrayRCP();

  /** \brief Assignment operator: Makes <tt>*this</tt> reference the input array.
   *
   * If the input array is a reference to <tt>*this</tt> (that is, if
   * <tt>this->getRawPtr() == r_ptr.getRawPtr()</tt>), then this
   * method does nothing.  Otherwise, it does the following:
   * <ol>
   * <li> Decrements the reference count of <tt>*this</tt> (as
   *   if its destructor had been called)
   * <li> Makes <tt>*this</tt> a reference to the input array
   *   (thus incrementing its reference count)
   * </ol>
   * For example, after the following sample code is done, the array
   * to which x originally pointed on construction will have reference
   * count 2, and the array to which y originally pointed on
   * constructor will have reference count 0 (and will thus be
   * deallocated).
     \code
     ArrayRCP<T> x = arcp<T> (10);
     ArrayRCP<T> y = arcp<T> (42);
     x = y;     
     \endcode
   *
   * \post <tt>this->get() == r_ptr.get()</tt>
   * \post <tt>this->count() == r_ptr.count()</tt>
   * \post <tt>this->has_ownership() == r_ptr.has_ownership()</tt>
   * \post If <tt>r_ptr.get() != NULL</tt> then <tt>r_ptr.count()</tt> is incremented by 1
   */
  inline ArrayRCP<T>& operator=(const ArrayRCP<T>& r_ptr);

  //@}

  //! @name Object/Pointer Access Functions 
  //@{

  /** \brief Returns true if the underlying pointer is null. */
  inline bool is_null() const;

  /** \brief Pointer (<tt>-></tt>) access to members of underlying object for
   * current position.
   *
   * \pre <tt>this->get() != NULL</tt>
   * \pre <tt>this->lowerOffset() <= 0</tt>
   * \pre <tt>this->upperOffset() >= 0</tt>
   */
  inline T* operator->() const;

  /** \brief Dereference the underlying object for the current pointer
   * position.
   *
   * \pre <tt>this->get() != NULL</tt>
   * \pre <tt>this->lowerOffset() <= 0</tt>
   * \pre <tt>this->upperOffset() >= 0</tt>
   */
  inline T& operator*() const;

  /** \brief Get the raw C++ pointer to the underlying object.
   *
   * \pre [<tt>*this != null</tt>] <tt>this->lowerOffset() <= 0</tt>
   * \pre [<tt>*this != null</tt>] <tt>this->upperOffset() >= 0</tt>
   */
  inline T* get() const;

  /** \brief Get the raw C++ pointer to the underlying object.
   *
   * \pre [<tt>*this != null</tt>] <tt>this->lowerOffset() <= 0</tt>
   * \pre [<tt>*this != null</tt>] <tt>this->upperOffset() >= 0</tt>
   */
  inline T* getRawPtr() const;

  /** \brief Random object access.
   *
   * \pre <tt>this->get() != NULL</tt>
   * \pre <tt>this->lowerOffset() <= offset && offset <= this->upperOffset()</tt>
   */
  inline T& operator[](size_type offset) const;

  //@}

  //! @name Pointer Arithmetic Functions 
  //@{

  /** \brief Prefix increment of pointer (i.e. ++ptr).
   *
   * Does nothing if <tt>this->get() == NULL</tt>.
   *
   * \post [<tt>this->get()!=NULL</tt>] <tt>this->get()</tt> is incremented by <tt>1</tt>
   * \post [<tt>this->get()!=NULL</tt>] <tt>this->lowerOffset()</tt> is decremented by <tt>1</tt>
   * \post [<tt>this->get()!=NULL</tt>] <tt>this->upperOffset()</tt> is decremented by <tt>1</tt>
   */
  inline ArrayRCP<T>& operator++();

  /** \brief Postfix increment of pointer (i.e. ptr++).
   *
   * Does nothing if <tt>this->get() == NULL</tt>.
   *
   * \post <tt>this->get()</tt> is incremented by <tt>1</tt>
   * \post <tt>this->lowerOffset()</tt> is decremented by <tt>1</tt>
   * \post <tt>this->upperOffset()</tt> is decremented by <tt>1</tt>
   */
  inline ArrayRCP<T> operator++(int);

  /** \brief Prefix decrement of pointer (i.e. --ptr).
   *
   * Does nothing if <tt>this->get() == NULL</tt>.
   *
   * \post [<tt>this->get()!=NULL</tt>] <tt>this->get()</tt> is decremented by <tt>1</tt>
   * \post [<tt>this->get()!=NULL</tt>] <tt>this->lowerOffset()</tt> is incremented by <tt>1</tt>
   * \post [<tt>this->get()!=NULL</tt>] <tt>this->upperOffset()</tt> is incremented by <tt>1</tt>
   */
  inline ArrayRCP<T>& operator--();

  /** \brief Postfix decrement of pointer (i.e. ptr--).
   *
   * Does nothing if <tt>this->get() == NULL</tt>.
   *
   * \post <tt>this->get()</tt> is decremented by <tt>1</tt>
   * \post <tt>this->lowerOffset()</tt> is incremented by <tt>1</tt>
   * \post <tt>this->upperOffset()</tt> is incremented by <tt>1</tt>
   */
  inline ArrayRCP<T> operator--(int);

  /** \brief Pointer integer increment (i.e. ptr+=offset).
   *
   * Does nothing if <tt>this->get() == NULL</tt>.
   *
   * \post [<tt>this->get()!=NULL</tt>] <tt>this->get()</tt> is incremented by <tt>offset</tt>
   * \post [<tt>this->get()!=NULL</tt>] <tt>this->lowerOffset()</tt> is decremented by <tt>offset</tt>
   * \post [<tt>this->get()!=NULL</tt>] <tt>this->upperOffset()</tt> is decremented by <tt>offset</tt>
   */
  inline ArrayRCP<T>& operator+=(size_type offset);

  /** \brief Pointer integer increment (i.e. ptr-=offset).
   *
   * Does nothing if <tt>this->get() == NULL</tt>.
   *
   * \post [<tt>this->get()!=NULL</tt>] <tt>this->get()</tt> is decremented by <tt>offset</tt>
   * \post [<tt>this->get()!=NULL</tt>] <tt>this->lowerOffset()</tt> is incremented by <tt>offset</tt>
   * \post [<tt>this->get()!=NULL</tt>] <tt>this->upperOffset()</tt> is incremented by <tt>offset</tt>
   */
  inline ArrayRCP<T>& operator-=(size_type offset);

  /** \brief Pointer integer increment (i.e. ptr+offset).
   *
   * Returns a null pointer if <tt>this->get() == NULL</tt>.
   *
   * \post [<tt>this->get()!=NULL</tt>] <tt>return->get() == this->get() + offset</tt>
   * \post [<tt>this->get()!=NULL</tt>] <tt>return->lowerOffset() == this->lowerOffset() - offset</tt>
   * \post [<tt>this->get()!=NULL</tt>] <tt>return->upperOffset() == this->upperOffset() - offset</tt>
   *
   * Note that since implicit conversion of <tt>ArrayRCP<T></tt>
   * objects is not allowed that it does not help at all to make this function
   * into a non-member function.
   */
  inline ArrayRCP<T> operator+(size_type offset) const;

  /** \brief Pointer integer decrement (i.e. ptr-offset).
   *
   * Returns a null pointer if <tt>this->get() == NULL</tt>.
   *
   * \post [<tt>this->get()!=NULL</tt>] <tt>return->get() == this->get() - offset</tt>
   * \post [<tt>this->get()!=NULL</tt>] <tt>return->lowerOffset() == this->lowerOffset() + offset</tt>
   * \post [<tt>this->get()!=NULL</tt>] <tt>return->upperOffset() == this->upperOffset() + offset</tt>
   *
   * Note that since implicit conversion of <tt>ArrayRCP<T></tt>
   * objects is not allowed that it does not help at all to make this function
   * into a non-member function.
   */
  inline ArrayRCP<T> operator-(size_type offset) const;

  //@}

  //! @name Standard Container-Like Functions 
  //@{

  /** \brief Return an iterator to beginning of the array of data.
   *
   * If <tt>HAVE_TEUCHOS_ARRAY_BOUNDSCHECK</tt> is defined then the iterator
   * returned is an <tt>ArrayRCP<T></tt> object and all operations are
   * checked at runtime. When <tt>HAVE_TEUCHOS_ARRAY_BOUNDSCHECK</tt> is not
   * defined, the a raw pointer <tt>T*</tt> is returned for fast execution.
   *
   * \post [this->get()!=NULL</tt>] <tt>&*return == this->get()</tt>
   * \post [<tt>this->get()==NULL</tt>] <tt>return == (null or NULL)</tt>
   */
  inline iterator begin() const;

  /** \brief Return an iterator to past the end of the array of data.
   *
   * If <tt>HAVE_TEUCHOS_ARRAY_BOUNDSCHECK</tt> is defined then the iterator
   * returned is an <tt>ArrayRCP<T></tt> object and all operations are
   * checked at runtime. When <tt>HAVE_TEUCHOS_ARRAY_BOUNDSCHECK</tt> is not
   * defined, the a raw pointer <tt>T*</tt> is returned for fast execution.
   *
   * \post [<tt>this->get()!=NULL</tt>] <tt>&*end == this->get()+(this->upperOffset()+1)</tt>
   * \post [<tt>this->get()==NULL</tt>] <tt>return == (null or NULL)</tt>
   */
  inline iterator end() const;

  //@}

  //! @name ArrayRCP Views 
  //@{

  /** \brief Return object for only const access to data.
   *
   * This function should only compile successfully if the type <tt>T</tt> is
   * not already declared <tt>const</tt>!
   */
  inline ArrayRCP<const T> getConst() const;

  /** \brief Return a persisting view of a contiguous range of elements.
   *
   * \pre <tt>this->get() != NULL</tt>
   * \pre <tt>this->lowerOffset() <= lowerOffset</tt>
   * \pre <tt>lowerOffset + size - 1 <= this->upperOffset()</tt>
   *
   * \post <tt>return->get() == this->get() + lowerOffset</tt>
   * \post <tt>return->lowerOffset() == 0</tt>
   * \post <tt>return->upperOffset() == size-1</tt>
   *
   * \note A <tt>size==0</tt> view of even a null ArrayRCP is allowed.
   *   It returns a <tt>null</tt> view in that case.
   */
  inline ArrayRCP<T> persistingView( size_type lowerOffset, size_type size ) const;

  //@}

  //! @name Size and extent query functions 
  //@{

  /** \brief Return the lower offset to valid data. */
  inline size_type lowerOffset() const;

  /** \brief Return the upper offset to valid data. */
  inline size_type upperOffset() const;

  /** \brief The total number of items in the managed array
   * (i.e. <tt>upperOffset()-lowerOffset()+1</tt>).
   */
  inline size_type size() const;

  //@}

  //! @name ArrayView views 
  //@{

  /** \brief Return view of a contiguous range of elements.
   *
   * \pre <tt>this->get() != NULL</tt>
   * \pre <tt>this->lowerOffset() <= lowerOffset</tt>
   * \pre <tt>lowerOffset + size - 1 <= this->upperOffset()</tt>
   *
   * \post <tt>return->get() == this->get() + lowerOffset</tt>
   * \post <tt>return->lowerOffset() == 0</tt>
   * \post <tt>return->upperOffset() == size-1</tt>
   *
   * \note A <tt>size==0</tt> view of even a null ArrayRCP is allowed.
   *   It returns a <tt>null</tt> view in that case.
   */
  inline ArrayView<T> view( size_type lowerOffset, size_type size ) const;

  /** \brief Return a nonpersisting view of a contiguous range of elements.
   *
   * This is equivalent to calling <tt>view (offset, size)</tt>.
   */
  inline ArrayView<T> operator()( size_type lowerOffset, size_type size ) const;

  /** \brief Return a nonpersisting view of <tt>*this</tt>.
   *
   * \note This will return a null ArrayView if <tt>this->size() == 0</tt>.
   */
  inline ArrayView<T> operator()() const;

  //@}
  //! @name Implicit conversions
  //@{

  /** \brief Convert from ArrayRCP<T> to ArrayRCP<const T>. */
  inline operator ArrayRCP<const T>() const;

  //@}

  //! @name std::vector like and other misc functions
  //@{

  /** \brief Resize and assign n elements of val.
   *
   * \postconditions <tt>size() == n</tt>
   */
  inline void assign(size_type n, const T &val);

  /** \brief Resize and assign to iterator sequence [first, last)
   *
   * \postconditions <tt>size() == std::distance(first, last)</tt>
   *
   * This will not change the underlying pointer array if the size does not
   * change.
   */
  template<class Iter>
  inline void assign(Iter first, Iter last);

  /** \brief Deep copy the elements from one ArrayView object into this
   * object.
   *
   * This is equivalent to calling <tt>assign (av.begin (), av.end ())</tt>
   */
  inline void deepCopy(const ArrayView<const T>& av);

  //! Resize and append new elements if necessary.
  inline void resize(const size_type n, const T &val = T());

  /** \brief Resize to zero.
   *
   * \postconditions <tt>size()==0</tt>
   */
  inline void clear();

  //@}
  //! @name Reference counting
  //@{

  /** \brief Strength of the pointer.
   *
   * Return values:<ul>
   * <li><tt>RCP_STRONG</tt>: Underlying reference-counted object will be deleted
   *     when <tt>*this</tt> is destroyed if <tt>strong_count()==1</tt>. 
   * <li><tt>RCP_WEAK</tt>: Underlying reference-counted object will not be deleted
   *     when <tt>*this</tt> is destroyed if <tt>strong_count() > 0</tt>. 
   * <li><tt>RCP_STRENGTH_INVALID</tt>: <tt>*this</tt> is not strong or weak but
   *     is null.
   * </ul>
   */
  inline ERCPStrength strength() const;

  /** \brief Return whether the underlying object pointer is still valid.
   *
   * The underlying object will not be valid if the strong count has
   * gone to zero but the weak count has not.
   *
   * NOTE: Null is a valid object pointer.  If you want to know if there is a
   * non-null object and it is valid then <tt>!is_null() &&
   * is_valid_ptr()</tt> will be <tt>true</tt>.
   */
  inline bool is_valid_ptr() const;

  /** \brief Return the number of active <tt>RCP<></tt> objects that have a
   * "strong" reference to the underlying reference-counted object.
   *
   * \return If <tt>this->get() == NULL</tt> then this function returns 0.
   */
  inline int strong_count() const;

  /** \brief Return the number of active <tt>RCP<></tt> objects that have a
   * "weak" reference to the underlying reference-counted object.
   *
   * \return If <tt>this->get() == NULL</tt> then this function returns 0.
   */
  inline int weak_count() const;

  /** \brief Total count (strong_count() + weak_count()). */
  inline int total_count() const;

  /** \brief Give <tt>this</tt> and other <tt>ArrayRCP<></tt> objects
   * ownership of the underlying referenced array to delete it.
   *
   * See <tt>~ArrayRCP()</tt> above. This function does nothing if
   * <tt>this->get() == NULL</tt>.
   *
   * <b>Postconditions:</b><ul>
   * <li> If <tt>this->get() == NULL</tt> then
   * <ul>
   * <li><tt>this->has_ownership() == false</tt> (always!).
   * </ul>
   * <li> else
   * <ul>
   * <li><tt>this->has_ownership() == true</tt>
   * </ul>
   * </ul>
   */
  inline void set_has_ownership();

  /** \brief Returns true if <tt>this</tt> has ownership of object pointed to
   * by <tt>this->get()</tt> in order to deallocate it.
   *
   * See the above documentation for the destructor.
   *
   * \return If this->get() <tt>== NULL</tt> then this function always returns
   * <tt>false</tt>. Otherwise the value returned from this function depends
   * on which function was called most recently, if any;
   * <tt>set_has_ownership()</tt> (<tt>true</tt>) or <tt>release()</tt>
   * (<tt>false</tt>).
   */
  inline bool has_ownership() const;

  /** \brief Release the ownership of the underlying array.
   *
   * After this function is called then the client is responsible for deleting
   * the returned pointer no matter how many <tt>ref_count_ptr<T></tt> objects
   * have a reference to it. If <tt>this-></tt>get() <tt>== NULL</tt>, then
   * this call is meaningless.
   *
   * Note that this function does not have the exact same semantics as does
   * <tt>auto_ptr<T>::release()</tt>. In <tt>auto_ptr<T>::release()</tt>,
   * <tt>this</tt> is set to <tt>NULL</tt> while here in ArrayRCP<T>::
   * release() only an ownership flag is set and <tt>this</tt> still points to
   * the same array. It would be difficult to duplicate the behavior of
   * <tt>auto_ptr<T>::release()</tt> for this class.
   *
   * <b>Postconditions:</b><ul>
   * <li><tt>this->has_ownership() == false</tt>
   * </ul>
   *
   * \returns Returns the value of <tt>this->get()</tt>
   */
  inline T* release();

  /** \brief Create a new weak reference from another (strong) reference.
   *
   * A "weak" reference gives access to the array, without
   * incrementing its (strong) reference count.  This lets you have
   * access to the array, without affecting when it gets deallocated.
   *
   * \pre <tt>returnVal.is_valid_ptr()==true</tt>
   *
   * \post <tt>returnVal.get() == this->get()</tt>
   * \post <tt>returnVal.strong_count() == this->strong_count()</tt>
   * \post <tt>returnVal.weak_count() == this->weak_count()+1</tt>
   * \post <tt>returnVal.strength() == RCP_WEAK</tt>
   * \post <tt>returnVal.has_ownership() == this->has_ownership()</tt>
   */
  inline ArrayRCP<T> create_weak() const;

  /** \brief Create a new strong RCP object from another (weak) RCP object.
   *
   * A "weak" reference gives access to the array, without
   * incrementing its (strong) reference count.  This method lets you
   * "promote" a weak reference into a strong reference.  If the array
   * has been deallocated, the returned reference is null.
   *
   * \pre <tt>returnVal.is_valid_ptr()==true</tt>
   *
   * \post <tt>returnVal.get() == this->get()</tt>
   * \post <tt>returnVal.strong_count() == this->strong_count()+1</tt>
   * \post <tt>returnVal.weak_count() == this->weak_count()</tt>
   * \post <tt>returnVal.strength() == RCP_STRONG</tt>
   * \post <tt>returnVal.has_ownership() == this->has_ownership()</tt>
   */
  inline ArrayRCP<T> create_strong() const;

  /** \brief Returns true if the smart pointers share the same underlying reference-counted object.
   *
   * This method does more than just check if <tt>this->get() == r_ptr.get()</tt>.
   * It also checks to see if the underlying reference counting machinery is the
   * same.
   */
  template<class T2>
  inline bool shares_resource(const ArrayRCP<T2>& r_ptr) const;

  //@}
  //! @name Assertion Functions. 
  //@{

  /** \brief Throws <tt>NullReferenceError</tt> if <tt>this->get()==NULL</tt>,
   * otherwise returns reference to <tt>*this</tt>.
   */
  inline const ArrayRCP<T>& assert_not_null() const;

  /** \brief Throws <tt>NullReferenceError</tt> if <tt>this->get()==NULL</tt>
   * or<tt>this->get()!=NULL</tt>, throws <tt>RangeError</tt> if
   * <tt>(lowerOffset < this->lowerOffset() || this->upperOffset() <
   * upperOffset</tt>, otherwise returns reference to <tt>*this</tt>
   */
  inline const ArrayRCP<T>& assert_in_range( size_type lowerOffset, size_type size ) const;

  /** \brief If the object pointer is non-null, assert that it is still valid.
   *
   * If <tt>is_null()==false && strong_count()==0</tt>, this will throw
   * <tt>DanglingReferenceErorr</tt> with a great error message.
   *
   * If <tt>is_null()==true</tt>, then this will not throw any exception.
   *
   * In this context, null is a valid object.
   */
  inline const ArrayRCP<T>& assert_valid_ptr() const;

  //@}


  /** \name Deprecated */
  //@{

  /** \brief Returns <tt>strong_count()</tt> [deprecated]. */
  inline TEUCHOS_DEPRECATED int count() const;

  //@}

private:

  // //////////////////////////////////////////////////////////////
  // Private data members

  T *ptr_; // NULL if this pointer is null
  RCPNodeHandle node_; // NULL if this pointer is null
  size_type lowerOffset_; // 0 if this pointer is null
  size_type upperOffset_; // -1 if this pointer is null

  inline void debug_assert_not_null() const
    {
#ifdef TEUCHOS_REFCOUNTPTR_ASSERT_NONNULL
      assert_not_null();
#endif
    }

  inline void debug_assert_in_range( size_type lowerOffset_in,
    size_type size_in ) const
    {
      (void)lowerOffset_in; (void)size_in;
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
      assert_in_range(lowerOffset_in, size_in);
#endif
    }

  inline void debug_assert_valid_ptr() const
    {
#ifdef TEUCHOS_DEBUG
      assert_valid_ptr();
#endif
    }

public:


#ifndef DOXYGEN_COMPILE
  // These constructors should be private but I have not had good luck making
  // this portable (i.e. using friendship etc.) in the past
  // This is a very bad breach of encapsulation that is needed since MS VC++
  // 5.0 will not allow me to declare template functions as friends.
  ArrayRCP( T* p, size_type lowerOffset, size_type size,
    const RCPNodeHandle& node );
  T* access_private_ptr() const;
  RCPNodeHandle& nonconst_access_private_node();
  const RCPNodeHandle& access_private_node() const;
#endif

};  // end class ArrayRCP<...>


/** \brief Dummy specialization of ArrayRCP<void>.
 *
 * ArrayRCP<void> cannot be parsed because of reference and const_reference
 * typedefs resolving to "void &" and "const void &".  This full template
 * specialization ArrayRCP<void> neglects these. This will be mentioned in the
 * context of Kokkos::DefaultSparseOps<void>.  However, DefaultSparseOps<void>
 * is never instantiated, and until there is a need (and the semantics have
 * been decided), ArrayRCP<void> may not be instantiated either.
 */
template<>
class ArrayRCP<void> {
public:

  //! @name Public types 
  //@{

  /** \brief. */
  typedef Teuchos_Ordinal Ordinal;

  /** \brief . */
  typedef Ordinal size_type;
  /** \brief . */
  typedef Ordinal difference_type;
  /** \brief . */
  typedef std::random_access_iterator_tag iterator_category;
  /** \brief . */
  typedef  void* iterator_type;
  /** \brief . */
  typedef  void value_type;
  /** \brief . */
  // typedef T& reference;              // these are not valid
  /** \brief . */
  // typedef const T& const_reference;  // these are not valid
  /** \brief . */
  typedef void* pointer;
  /** \brief . */
  typedef void* const_pointer;
  /** \brief . */
  typedef void  element_type;

  /** \brief Default constructor, thows an exception.
   */
  inline ArrayRCP( );

  //@}

};  // end class ArrayRCP<void>

/** \brief Dummy specialization of ArrayRCP<const void>.
 *
 * See ArrayRCP<void> for details.
 */
template<>
class ArrayRCP<const void> {
public:

  //! @name Public types 
  //@{

  /** \brief. */
  typedef Teuchos_Ordinal Ordinal;

  /** \brief . */
  typedef Ordinal size_type;
  /** \brief . */
  typedef Ordinal difference_type;
  /** \brief . */
  typedef std::random_access_iterator_tag iterator_category;
  /** \brief . */
  typedef  const void* iterator_type;
  /** \brief . */
  typedef  const void value_type;
  /** \brief . */
  // typedef T& reference;              // these are not valid
  /** \brief . */
  // typedef const T& const_reference;  // these are not valid
  /** \brief . */
  typedef const void* pointer;
  /** \brief . */
  typedef const void* const_pointer;
  /** \brief . */
  typedef const void  element_type;

  /** \brief Default constructor, thows an exception.
   */
  inline ArrayRCP( );

  //@}

};  // end class ArrayRCP<void>

// 2008/09/22: rabartl: NOTE: I removed the TypeNameTraits<ArrayRCP<T> >
// specialization since I want to be able to print the type name of an
// ArrayRCP that does not have the type T fully defined!


/** \brief Traits specialization for ArrayRCP.
 *
 * \relates ArrayRCP
 */
template<typename T>
class NullIteratorTraits<ArrayRCP<T> > {
public:
  static ArrayRCP<T> getNull() { return null; }
};


/** \brief Wraps a preallocated array of data with the assumption to call the
 * array version of delete.
 *
 * \relates ArrayRCP
 */
template<class T>
ArrayRCP<T> arcp(
  T* p,
  typename ArrayRCP<T>::size_type lowerOffset,
  typename ArrayRCP<T>::size_type size,
  bool owns_mem = true
  );


/** \brief Wraps a preallocated array of data and uses a templated
 * deallocation strategy object to define deletion .
 *
 * \relates ArrayRCP
 */
template<class T, class Dealloc_T>
ArrayRCP<T> arcp(
  T* p,
  typename ArrayRCP<T>::size_type lowerOffset,
  typename ArrayRCP<T>::size_type size,
  Dealloc_T dealloc, bool owns_mem
  );

 
/** \brief Allocate a new array just given a dimension.
 *
 * <b>Warning!</b> The memory is allocated using <tt>new T[size]</tt> and is
 * *not* initialized (unless there is a default constructor for a user-defined
 * type).
 *
 * When called with 'size == 0' it returns a null ArrayRCP object.
 *
 * \relates ArrayRCP
 */
template<class T>
ArrayRCP<T> arcp( typename ArrayRCP<T>::size_type size );

 
/** \brief Allocate a new ArrayRCP object with a new RCPNode with memory
 * pointing to the initial node.
 *
 * The purpose of this function is to create a new "handle" to the array of
 * memory with its own seprate reference count.  The new ArrayRCP object will
 * have a new RCPNodeTmpl object that has a copy of the input ArrayRCP object
 * embedded in it.  This maintains the correct reference counting behaviors
 * but now gives a private count.  One would want to use arcpCloneNode(...) 
 * whenever it is important to keep a private reference count which is needed
 * for some types of use cases.
 *
 * \relates ArrayRCP
 */
template<class T>
ArrayRCP<T> arcpCloneNode( const ArrayRCP<T> &a );

 
/** \brief Allocate a new array by cloning data from an input array view.
 *
 * \relates ArrayRCP
 */
template<class T>
ArrayRCP<T> arcpClone( const ArrayView<const T> &v );


/** \brief Create an ArrayRCP with and also put in an embedded object.
 *
 * In this case the embedded object is destroyed (by setting to Embedded())
 * before the object at <tt>*p</tt> is destroyed.
 *
 * The embedded object can be extracted using <tt>getEmbeddedObj()</tt> and
 * <tt>getNonconstEmbeddedObject()</tt>.
 *
 * \relates ArrayRCP
 */
template<class T, class Embedded>
ArrayRCP<T>
arcpWithEmbeddedObjPreDestroy(
  T* p,
  typename ArrayRCP<T>::size_type lowerOffset,
  typename ArrayRCP<T>::size_type size,
  const Embedded &embedded,
  bool owns_mem = true
  );


/** \brief Create an ArrayRCP with and also put in an embedded object.
 *
 * In this case the embedded object is destroyed (by setting to Embedded())
 * after the object at <tt>*p</tt> is destroyed.
 *
 * The embedded object can be extracted using <tt>getEmbeddedObj()</tt> and
 * <tt>getNonconstEmbeddedObject()</tt>.
 *
 * \relates ArrayRCP
 */
template<class T, class Embedded>
ArrayRCP<T>
arcpWithEmbeddedObjPostDestroy(
  T* p,
  typename ArrayRCP<T>::size_type lowerOffset,
  typename ArrayRCP<T>::size_type size,
  const Embedded &embedded,
  bool owns_mem = true
  );


/** \brief Create an ArrayRCP with and also put in an embedded object.
 *
 * This function should be called when it is not important when the embedded
 * object is destroyed (by setting to Embedded()) with respect to when
 * <tt>*p</tt> is destroyed.
 *
 * The embedded object can be extracted using <tt>getEmbeddedObj()</tt> and
 * <tt>getNonconstEmbeddedObject()</tt>.
 *
 * \relates ArrayRCP
 */
template<class T, class Embedded>
ArrayRCP<T>
arcpWithEmbeddedObj(
  T* p,
  typename ArrayRCP<T>::size_type lowerOffset,
  typename ArrayRCP<T>::size_type size,
  const Embedded &embedded,
  bool owns_mem = true
  );


/** \brief Wrap an <tt>std::vector<T></tt> object as an
 * <tt>ArrayRCP<T></tt> object.
 *
 * \relates ArrayRCP
 */
template<class T>
ArrayRCP<T> arcp( const RCP<std::vector<T> > &v );


/** \brief Wrap a <tt>const std::vector<T></tt> object as an
 * <tt>ArrayRCP<const T></tt> object.
 *
 * \relates ArrayRCP
 */
template<class T>
ArrayRCP<const T> arcp( const RCP<const std::vector<T> > &v );


/** \brief Get an ArrayRCP object out of an ArrayView object.
 *
 * This conversion is required an proper in certain types of situations.  In a
 * debug build, a dangling reference will be detected with an exception being
 * thrown.
 *
 * \relates ArrayRCP
 */
template<class T>
ArrayRCP<T> arcpFromArrayView(const ArrayView<T> &av);


/** \brief Get an <tt>std::vector<T></tt> object out of an
 * <tt>ArrayRCP<T></tt> object that was created using the
 * <tt>arcp()</tt> above to wrap the std::vector in the first place..
 *
 * \relates ArrayRCP
 */
template<class T>
RCP<std::vector<T> > get_std_vector( const ArrayRCP<T> &ptr );


/** \brief Get a <tt>const std::vector<T></tt> object out of an
 * <tt>ArrayRCP<const T></tt> object that was created using the
 * <tt>arcp()</tt> above to wrap the std::vector in the first place.
 *
 * \relates ArrayRCP
 */
template<class T>
RCP<const std::vector<T> > get_std_vector( const ArrayRCP<const T> &ptr );


/** \brief Returns true if <tt>p.get()==NULL</tt>.
 *
 * \relates ArrayRCP
 */
template<class T>
bool is_null( const ArrayRCP<T> &p );


/** \brief Returns true if <tt>p.get()!=NULL</tt>.
 *
 * \relates ArrayRCP
 */
template<class T>
bool nonnull( const ArrayRCP<T> &p );


/** \brief Returns true if <tt>p.get()==NULL</tt>.
 *
 * \relates ArrayRCP
 */
template<class T>
bool operator==( const ArrayRCP<T> &p, ENull );


/** \brief Returns true if <tt>p.get()!=NULL</tt>.
 *
 * \relates ArrayRCP
 */
template<class T>
bool operator!=( const ArrayRCP<T> &p, ENull );


/** \brief .
 *
 * \relates ArrayRCP
 */
template<class T1, class T2>
bool operator==( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 );


/** \brief .
 *
 * \relates ArrayRCP
 */
template<class T1, class T2>
bool operator!=( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 );


/** \brief .
 *
 * \relates ArrayRCP
 */
template<class T1, class T2>
bool operator<( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 );


/** \brief .
 *
 * \relates ArrayRCP
 */
template<class T1, class T2>
bool operator<=( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 );


/** \brief .
 *
 * \relates ArrayRCP
 */
template<class T1, class T2>
bool operator>( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 );


/** \brief .
 *
 * \relates ArrayRCP
 */
template<class T1, class T2>
bool operator>=( const ArrayRCP<T1> &p1, const ArrayRCP<T2> &p2 );


/** \brief Returns difference of two ArrayRCP object</tt>.
 *
 * \relates ArrayRCP
 */
template<class T>
typename ArrayRCP<T>::difference_type
operator-( const ArrayRCP<T> &p1, const ArrayRCP<T> &p2 );


/** \brief Const cast of underlying <tt>ArrayRCP</tt> type from <tt>const T*</tt>
 * to <tt>T*</tt>.
 *
 * The function will compile only if (<tt>const_cast<T2*>(p1.get());</tt>)
 * compiles.
 *
 * \relates ArrayRCP
 */
template<class T2, class T1>
inline
ArrayRCP<T2> arcp_const_cast(const ArrayRCP<T1>& p1);


/** \brief Reinterpret cast of underlying <tt>ArrayRCP</tt> type from
 * <tt>T1*</tt> to <tt>T2*</tt>.
 *
 * The function will compile only if (<tt>reinterpret_cast<T2*>(p1.get());</tt>) compiles.
 *
 * <b>Warning!</b> Do not use this function unless you absolutely know what
 * you are doing. Doing a reinterpret cast is always a tricking thing and
 * must only be done by developers who are 100% comfortable with what they are
 * doing.
 *
 * \relates ArrayRCP
 */
template<class T2, class T1>
ArrayRCP<T2> arcp_reinterpret_cast(const ArrayRCP<T1>& p1);


/** \brief Reinterpret cast of underlying <tt>ArrayRCP</tt> type from
 * <tt>T1*</tt> to <tt>T2*</tt> where <tt>T2</tt> is a non-POD
 * (non-plain-old-data).
 *
 * The function will compile only if (<tt>reinterpret_cast<T2*>(p1.get());</tt>) compiles.
 *
 * This function is used to reinterpret cast an array of plain-old-data (POD)
 * (e.g. <tt>int</tt> or <tt>char</tt>) into an array of real objects.  The
 * constructors will be called on each of the memory locations with placement
 * new and the destructors will get called when the last RCP goes away.
 *
 * <b>Warning!</b> Do not use this function unless you absolutely know what
 * you are doing. Doing a reinterpret cast is always a tricking thing and
 * must only be done by developers who are 100% comfortable with what they are
 * doing.
 *
 * \relates ArrayRCP
 */
template<class T2, class T1>
ArrayRCP<T2> arcp_reinterpret_cast_nonpod(const ArrayRCP<T1>& p1, const T2& val=T2());


/** \brief Implicit case the underlying <tt>ArrayRCP</tt> type from
 * <tt>T1*</tt> to <tt>T2*</tt>.
 *
 * The function will compile only if (<tt>T2 *p = p1.get();</tt>) compiles.
 *
 * <b>Warning!</b> Do not use this function unless you absolutely know what you
 * are doing. While implicit casting of pointers to single objects is usually
 * 100% safe, implicit casting pointers to arrays of objects can be very
 * dangerous. One std::exception that is always safe is when you are implicit
 * casting an array of pointers to non-const objects to an array of const
 * pointers to const objects. For example, the following implicit conversion
 * from a array pointer objects <tt>aptr1</tt> of type
 * <tt>ArrayRCP<T*></tt> to 

 \code

 ArrayRCP<const T * const>
 aptr2 = arcp_implicit_cast<const T * const>(ptr1);

 \endcode

 * is always legal and safe to do.
 *
 * \relates ArrayRCP
 */
template<class T2, class T1>
inline
ArrayRCP<T2> arcp_implicit_cast(const ArrayRCP<T1>& p1);


/** \brief Set extra data associated with a <tt>ArrayRCP</tt> object.
 *
 * \param extra_data [in] Data object that will be set (copied)
 *
 * \param name [in] The name given to the extra data. The value of
 * <tt>name</tt> together with the data type <tt>T1</tt> of the extra data
 * must be unique from any other such data or the other data will be
 * overwritten.
 *
 * \param p [out] On output, will be updated with the input
 * <tt>extra_data</tt>
 *
 * \param destroy_when [in] Determines when <tt>extra_data</tt> will be
 * destroyed in relation to the underlying reference-counted object. If
 * <tt>destroy_when==PRE_DESTROY</tt> then <tt>extra_data</tt> will be deleted
 * before the underlying reference-counted object. If
 * <tt>destroy_when==POST_DESTROY</tt> (the default) then <tt>extra_data</tt>
 * will be deleted after the underlying reference-counted object.
 *
 * \param force_unique [in] Determines if this type and name pair must be
 * unique in which case if an object with this same type and name already
 * exists, then an std::exception will be thrown. The default is
 * <tt>true</tt> for safety.
 *
 * If there is a call to this function with the same type of extra
 * data <tt>T1</tt> and same arguments <tt>p</tt> and <tt>name</tt>
 * has already been made, then the current piece of extra data already
 * set will be overwritten with <tt>extra_data</tt>. However, if the
 * type of the extra data <tt>T1</tt> is different, then the extra
 * data can be added and not overwrite existing extra data. This
 * means that extra data is keyed on both the type and name. This
 * helps to minimize the chance that clients will unexpectedly
 * overwrite data by accident.
 *
 * When the last <tt>RefcountPtr</tt> object is removed and the
 * reference-count node is deleted, then objects are deleted in the following
 * order: (1) All of the extra data that where added with
 * <tt>destroy_when==PRE_DESTROY</tt> are first, (2) then the underlying
 * reference-counted object is deleted, and (3) the rest of the extra data
 * that was added with <tt>destroy_when==PRE_DESTROY</tt> is then deleted.
 * The order in which the objects are destroyed is not guaranteed. Therefore,
 * clients should be careful not to add extra data that has deletion
 * dependencies (instead consider using nested ArrayRCP objects as extra
 * data which will guarantee the order of deletion).
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>p->get() != NULL</tt> (throws <tt>NullReferenceError</tt>)
 * <li> If this function has already been called with the same template
 * type <tt>T1</tt> for <tt>extra_data</tt> and the same std::string <tt>name</tt>
 * and <tt>force_unique==true</tt>, then an <tt>std::invalid_argument</tt>
 * std::exception will be thrown.
 * </ul>
 *
 * Note, this function is made a non-member function to be consistent
 * with the non-member <tt>get_extra_data()</tt> functions.
 *
 * \relates ArrayRCP
 */
template<class T1, class T2>
void set_extra_data(
  const T1 &extra_data, const std::string& name,
  const Ptr<ArrayRCP<T2> > &p, EPrePostDestruction destroy_when = POST_DESTROY,
  bool force_unique = true );


/** \brief Get a non-const reference to extra data associated with a <tt>ArrayRCP</tt> object.
 *
 * \param p [in] Smart pointer object that extra data is being extracted from.
 *
 * \param name [in] Name of the extra data.
 *
 * \returns Returns a non-const reference to the extra_data object.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>p.get() != NULL</tt> (throws <tt>NullReferenceError</tt>)
 * <li><tt>name</tt> and <tt>T1</tt> must have been used in a previous
 * call to <tt>set_extra_data()</tt> (throws <tt>std::invalid_argument</tt>).
 * </ul>
 *
 * Note, this function must be a non-member function since the client
 * must manually select the first template argument.
 *
 * \relates ArrayRCP
 */
template<class T1, class T2>
T1& get_extra_data( ArrayRCP<T2>& p, const std::string& name );


/** \brief Get a const reference to extra data associated with a <tt>ArrayRCP</tt> object.
 *
 * \param p [in] Smart pointer object that extra data is being extracted from.
 *
 * \param name [in] Name of the extra data.
 *
 * \returns Returns a const reference to the extra_data object.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>p.get() != NULL</tt> (throws <tt>NullReferenceError</tt>)
 * <li><tt>name</tt> and <tt>T1</tt> must have been used in a previous
 * call to <tt>set_extra_data()</tt> (throws <tt>std::invalid_argument</tt>).
 * </ul>
 *
 * Note, this function must be a non-member function since the client
 * must manually select the first template argument.
 *
 * Also note that this const version is a false sense of security
 * since a client can always copy a const <tt>ArrayRCP</tt> object
 * into a non-const object and then use the non-const version to
 * change the data. However, its presence will help to avoid some
 * types of accidental changes to this extra data.
 *
 * \relates ArrayRCP
 */
template<class T1, class T2>
const T1& get_extra_data( const ArrayRCP<T2>& p, const std::string& name );


/** \brief Get a pointer to non-const extra data (if it exists) associated
 * with a <tt>ArrayRCP</tt> object.
 *
 * \param p [in] Smart pointer object that extra data is being extracted from.
 *
 * \param name [in] Name of the extra data.
 *
 * \returns Returns a non-const pointer to the extra_data object.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>p.get() != NULL</tt> (throws <tt>NullReferenceError</tt>)
 * </ul>
 *
 * <b>Postconditions:</b><ul>
 * <li> If <tt>name</tt> and <tt>T1</tt> have been used in a previous
 * call to <tt>set_extra_data()</tt> then <tt>return !=NULL</tt>
 * and otherwise <tt>return == NULL</tt>.
 * </ul>
 *
 * Note, this function must be a non-member function since the client
 * must manually select the first template argument.
 *
 * \relates ArrayRCP
 */
template<class T1, class T2>
T1* get_optional_extra_data( ArrayRCP<T2>& p, const std::string& name );


/** \brief Get a pointer to const extra data (if it exists) associated with a <tt>ArrayRCP</tt> object.
 *
 * \param p [in] Smart pointer object that extra data is being extracted from.
 *
 * \param name [in] Name of the extra data.
 *
 * \returns Returns a const pointer to the extra_data object if it exists.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>p.get() != NULL</tt> (throws <tt>NullReferenceError</tt>)
 * </ul>
 *
 * <b>Postconditions:</b><ul>
 * <li> If <tt>name</tt> and <tt>T1</tt> have been used in a previous
 * call to <tt>set_extra_data()</tt> then <tt>return !=NULL</tt>
 * and otherwise <tt>return == NULL</tt>.
 * </ul>
 *
 * Note, this function must be a non-member function since the client
 * must manually select the first template argument.
 *
 * Also note that this const version is a false sense of security
 * since a client can always copy a const <tt>ArrayRCP</tt> object
 * into a non-const object and then use the non-const version to
 * change the data. However, its presence will help to avoid some
 * types of accidental changes to this extra data.
 *
 * \relates ArrayRCP
 */
template<class T1, class T2>
const T1* get_optional_extra_data( const ArrayRCP<T2>& p, const std::string& name );


/** \brief Return a non-<tt>const</tt> reference to the underlying deallocator object.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>p.get() != NULL</tt> (throws <tt>NullReferenceError</tt>)
 * <li> The deallocator object type used to construct <tt>p</tt> is same as <tt>Dealloc_T</tt>
 * (throws <tt>NullReferenceError</tt>)
 * </ul>
 *
 * \relates ArrayRCP
 */
template<class Dealloc_T, class T>
Dealloc_T& get_nonconst_dealloc( const ArrayRCP<T>& p );


/** \brief Return a <tt>const</tt> reference to the underlying deallocator object.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>p.get() != NULL</tt> (throws <tt>NullReferenceError</tt>)
 * <li> The deallocator object type used to construct <tt>p</tt> is same as <tt>Dealloc_T</tt>
 * (throws <tt>NullReferenceError</tt>)
 * </ul>
 *
 * Note that the <tt>const</tt> version of this function provides only
 * a very ineffective attempt to avoid accidental changes to the
 * deallocation object. A client can always just create a new
 * non-<tt>const</tt> <tt>ArrayRCP<T></tt> object from any
 * <tt>const</tt> <tt>ArrayRCP<T></tt> object and then call the
 * non-<tt>const</tt> version of this function.
 *
 * \relates ArrayRCP
 */
template<class Dealloc_T, class T>
const Dealloc_T& get_dealloc( const ArrayRCP<T>& p );


/** \brief Return a pointer to the underlying non-<tt>const</tt> deallocator
 * object if it exists.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>p.get() != NULL</tt> (throws <tt>NullReferenceError</tt>)
 * </ul>
 *
 * <b>Postconditions:</b><ul>
 * <li> If the deallocator object type used to construct <tt>p</tt> is same as <tt>Dealloc_T</tt>
 * then <tt>return!=NULL</tt>, otherwise <tt>return==NULL</tt>
 * </ul>
 *
 * \relates ArrayRCP
 */
template<class Dealloc_T, class T>
const Dealloc_T* get_optional_dealloc( const ArrayRCP<T>& p );


/** \brief Return a pointer to the underlying <tt>const</tt> deallocator
 * object if it exists.
 *
 * <b>Preconditions:</b><ul>
 * <li><tt>p.get() != NULL</tt> (throws <tt>NullReferenceError</tt>)
 * </ul>
 *
 * <b>Postconditions:</b><ul>
 * <li> If the deallocator object type used to construct <tt>p</tt> is same as <tt>Dealloc_T</tt>
 * then <tt>return!=NULL</tt>, otherwise <tt>return==NULL</tt>
 * </ul>
 *
 * Note that the <tt>const</tt> version of this function provides only
 * a very ineffective attempt to avoid accidental changes to the
 * deallocation object. A client can always just create a new
 * non-<tt>const</tt> <tt>ArrayRCP<T></tt> object from any
 * <tt>const</tt> <tt>ArrayRCP<T></tt> object and then call the
 * non-<tt>const</tt> version of this function.
 *
 * \relates ArrayRCP
 */
template<class Dealloc_T, class T>
Dealloc_T* get_optional_nonconst_dealloc( const ArrayRCP<T>& p );


/** \brief Get a const reference to an embedded object that was set by calling
 * <tt>arcpWithEmbeddedObjPreDestroy()</tt>,
 * <tt>arcpWithEmbeddedObjPostDestory()</tt>, or <tt>arcpWithEmbeddedObj()</tt>.
 *
 * \relates ArrayRCP
 */
template<class TOrig, class Embedded, class T>
const Embedded& getEmbeddedObj( const ArrayRCP<T>& p );


/** \brief Get a const reference to an embedded object that was set by calling
 * <tt>arcpWithEmbeddedObjPreDestroy()</tt>,
 * <tt>arcpWithEmbeddedObjPostDestory()</tt>, or <tt>arcpWithEmbeddedObj()</tt>.
 *
 * \relates ArrayRCP
 */
template<class TOrig, class Embedded, class T>
Embedded& getNonconstEmbeddedObj( const ArrayRCP<T>& p );


/** \brief Output stream inserter.
 *
 * The implementation of this function just print pointer addresses and
 * therefore puts not restrictions on the data types involved.
 *
 * \relates ArrayRCP
 */
template<class T>
std::ostream& operator<<( std::ostream& out, const ArrayRCP<T>& p );


} // end namespace Teuchos


#endif  // TEUCHOS_ARRAY_RCP_DECL_HPP
