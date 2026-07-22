// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_ARRAY_VIEW_DECL_HPP
#define TEUCHOS_ARRAY_VIEW_DECL_HPP


#include "Teuchos_RCPNode.hpp"
#include "Teuchos_ENull.hpp"
#include "Teuchos_NullIteratorTraits.hpp"
#include <vector>

namespace Teuchos {

// Forward declaration; ArrayView uses ArrayRCP in debug mode.
template<class T> class ArrayRCP;

/** \brief Nonowning array view.
 * \tparam T The type of each element in the array.
 *
 * This class provides a nonowning view of a one-dimensional array
 * with zero or more entries.  It holds a pointer to the data, and the
 * number of entries in the view.  "Nonowning" means that it does not
 * manage the array's memory.  This means two things.  First,
 * ArrayView's destructor does not deallocate the array.  Second, if
 * the array's memory is deallocated while the ArrayView is in scope,
 * any further use of the ArrayView or its iterators will result in
 * undefined behavior.
 *
 * The <tt>ArrayView</tt> class has the same shallow copy semantics of
 * the <tt>Ptr</tt> class.  <tt>ArrayView</tt> is to <tt>ArrayRCP</tt>
 * as <tt>Ptr</tt> is to <tt>RCP</tt>.
 *
 * \section Teuchos_ArrayView_Bounds Optional bounds checking
 *
 * You may enable bounds checking and other safety checks for this
 * class by setting the <tt>Teuchos_ENABLE_DEBUG:BOOL=ON</tt> CMake
 * option when configuring your Trilinos build.  This option is off by
 * default.  It incurs a significant performance penalty and so is not
 * recommended for production builds.  Bounds checking requires that
 * you always create ArrayView instances with the correct range.  For
 * example, if you use one of the constructors that accepts a raw
 * pointer, you are responsible for supplying the correct number of
 * elements in the array.  Our bounds checking implementation does not
 * attempt to replace memory debugging tools such as the Memcheck tool
 * in <a href="http://en.wikipedia.org/wiki/Valgrind">Valgrind</a>.
 *
 * \section Teuchos_ArrayView_Req Requirements on the type T
 *
 * ArrayView imposes the following requirements on the type T of
 * elements in the array:
 * <ul>
 * <li> T must be default constructible.
 * <li> T must be copy constructible.
 * <li> TypeNameTraits must have a specialization for T.
 * </ul>
 *
 * \section Teuchos_ArrayView_DesignDiscussion_sec Design discussion
 *
 * This class has a partial specialization for <tt>const T</tt> that
 * omits the conversion operator <tt>operator ArrayView<const T>()
 * const</tt>, and the assign() method (which performs a deep copy).
 * The conversion operator does not make sense if T is already
 * <tt>const T'</tt> for some type <tt>T'</tt>, and the assign()
 * method does not make sense if the right-hand side of the assignment
 * is const.
 *
 * Partial specialization results in duplicated code, so Teuchos
 * developers should be careful to make modifications in both the
 * fully generic implementation and in the partial specialization.
 *
 * We considered avoiding most of the duplication by making
 * <tt>ArrayView<T></tt> and its partial specialization
 * <tt>ArrayView<const T></tt> inherit from a common base class, which
 * contains all the common code.  However, the circular dependency
 * between ArrayRCP and ArrayView would have complicated this
 * solution.  We chose instead the simple "partial specialization
 * without a common base class" solution, which does not interfere
 * with the ArrayRCP / ArrayView circular dependency.
 *
 * \ingroup teuchos_mem_mng_grp
 */
template<class T>
class ArrayView {
public:
  //! @name Public typedefs
  //@{

  //! Integer index type used throughout ArrayView.
  typedef Teuchos_Ordinal Ordinal;

  //! Type representing the number of elements in an ArrayRCP or view thereof.
  typedef Ordinal size_type;

  //! Type representing the difference between two size_type values.
  typedef Ordinal difference_type;

  //! Type of each array element.
  typedef T value_type;

  /// \brief Type of a pointer to an array element.
  ///
  /// It may be const or nonconst, depending on T.
  typedef T* pointer;

  //! Type of a const pointer to an array element.
  typedef const T* const_pointer;

  /// \brief Type of a reference to an array element.
  ///
  /// It may be const or nonconst, depending on T.
  typedef T& reference;

  //! Type of a const reference to an array element.
  typedef const T& const_reference;

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  //! Type of a nonconst iterator.
  typedef ArrayRCP<T> iterator;
  //! Type of a const iterator.
  typedef ArrayRCP<const T> const_iterator;
#else
  //! Type of a nonconst iterator.
  typedef pointer iterator;
  //! Type of a const iterator.
  typedef const_pointer const_iterator;
#endif

  //@}
  //! @name Constructors/Destructors
  //@{

  //! Constructor that initializes to NULL (implicitly or explicitly).
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
   */
  ArrayView (T* p, size_type size,
             const ERCPNodeLookup rcpNodeLookup = RCP_ENABLE_NODE_LOOKUP);

  /** \brief Initialize from another <tt>ArrayView<T></tt> object.
   *
   * After construction, <tt>this</tt> and <tt>array</tt> will reference the
   * same array.
   *
   * This form of the copy constructor is required even though the
   * below more general templated version is sufficient since some
   * compilers will generate this function automatically which will
   * give an incorrect implementation.
   */
  ArrayView (const ArrayView<T>& array);

  //! Create a nonconst view of an std::vector<T>.
  ArrayView (std::vector<typename std::remove_const_t<T>>& vec);

  //! Create a const view of an std::vector<T>.
  ArrayView (const std::vector<typename std::remove_const_t<T>>& vec);

  //! Shallow copy assignment operator.
  ArrayView<T>& operator= (const ArrayView<T>& array);

  //! Destructor.
  ~ArrayView();

  //@}
  //! @name General query functions
  //@{

  /** \brief Returns true if the underlying pointer is null. */
  bool is_null() const;

  /** \brief The total number of items in the managed array. */
  size_type size() const;

  /** \brief Convert an ArrayView<T> to an <tt>std::string</tt> */
  std::string toString() const;

  //@}
  //! @name Element Access Functions
  //@{

  /** \brief Return a raw pointer to beginning of array or NULL if unsized. */
  inline T* getRawPtr() const;

  /// \brief Return a raw pointer to beginning of array.
  ///
  /// Same semantics as \c getRawPtr (which see).
  inline T* data() const;

  /** \brief Random object access.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->get() != NULL</tt>
   * <li><tt>0 <= i && i < this->size()</tt>
   * </ul>
   */
  T& operator[](size_type i) const;

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
   *
   * NOTE: A <tt>size==0</tt> view of even a null ArrayView is allowed and
   * returns a <tt>null</tt> view.
   */
  ArrayView<T> view( size_type offset, size_type size ) const;

  /** \brief Return a view of a contiguous range of elements (calls
   * view(offset, size)).
   */
  ArrayView<T> operator()( size_type offset, size_type size ) const;

  //! \brief Return <tt>*this</tt> (just for compatibility with Array and ArrayPtr).
  const ArrayView<T>& operator() () const;

  //! Return a const view of a possibly nonconst view.
  ArrayView<const T> getConst() const;

  /** \brief Implicitly convert an ArrayView<T> to an ArrayView<const T>.
   *
   * \note This conversion operator does not exist if T is already a
   *   const type (that is, if T is <tt>const T'</tt> for some type
   *   <tt>T'</tt>).  In that case, the assignment operator and copy
   *   constructor achieve the same syntactic effect.
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
   * \note This method does not exist if T is already a const type
   *   (that is, if T is <tt>const T'</tt> for some type <tt>T'</tt>).
   *   This is because assignment to a const right-hand side does not
   *   make sense.
   *
   * \note This function does modify the right-hand side's data.
   *   However, it is declared const, because it does not change the
   *   right-hand side's pointer or the number of entries in the view.
   *   The pointer is const, even though the data (to which the
   *   pointer points) are not.
   */
  void assign (const ArrayView<const T>& array) const;

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
   * <li>[<tt>this->get()!=NULL</tt>] <tt>&*return == this->get()</tt>
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
  const ArrayView<T>& assert_in_range( size_type offset, size_type size ) const;

  //@}

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  // I should make these private but templated friends are not very portable.
  // Besides, if a client directly calls this it will not compile in an
  // optimized build.

  explicit ArrayView( const ArrayRCP<T> &arcp );

  explicit ArrayView( T* p, size_type size, const ArrayRCP<T> &arcp );

#endif

private:
  T *ptr_; //<! Pointer to the data
  size_type size_; //<! Number of entries in the view
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  ArrayRCP<T> arcp_;
#endif

  void setUpIterators(const ERCPNodeLookup rcpNodeLookup = RCP_ENABLE_NODE_LOOKUP);

  void debug_assert_not_null() const {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    assert_not_null();
#endif
  }

  void debug_assert_in_range( size_type offset, size_type size_in ) const {
    (void)offset; (void)size_in;
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    assert_in_range(offset, size_in);
#endif
  }

  void debug_assert_valid_ptr() const {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    arcp_.access_private_node().assert_valid_ptr(*this);
#endif
  }

public: // Bad bad bad

  // This is a very bad breach of encapsulation but it exists to avoid
  // problems with portability of tempalted friends
  T* access_private_ptr() const
    { return ptr_; }

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  ArrayRCP<T> access_private_arcp() const
    { return arcp_; }
#endif

};


/** \brief Partial specialization of ArrayView for const T.
 *
 * The main documentation for ArrayView explains why this class needs
 * a partial specialization for const types.
 *
 * \ingroup teuchos_mem_mng_grp
 */
template<class T>
class ArrayView<const T> {
public:
  typedef Teuchos_Ordinal Ordinal;
  typedef Ordinal size_type;
  typedef Ordinal difference_type;
  typedef const T value_type;
  typedef const T* pointer;
  typedef const T* const_pointer;
  typedef const T& reference;
  typedef const T& const_reference;

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  typedef ArrayRCP<const T> iterator;
  typedef ArrayRCP<const T> const_iterator;
#else
  typedef pointer iterator;
  typedef const_pointer const_iterator;
#endif

  ArrayView( ENull null_arg = null );

  ArrayView (const T* p, size_type size,
             const ERCPNodeLookup rcpNodeLookup = RCP_ENABLE_NODE_LOOKUP );

  ArrayView (const ArrayView<const T>& array);

  ArrayView (std::vector<typename std::remove_const_t<T>>& vec);

  ArrayView (const std::vector<typename std::remove_const_t<T>>& vec);

  ArrayView<const T>& operator= (const ArrayView<const T>& array);

  ~ArrayView();

  bool is_null() const;

  size_type size() const;

  std::string toString() const;

  inline const T* getRawPtr() const;

  inline const T* data() const;

  const T& operator[] (size_type i) const;

  const T& front() const;

  const T& back() const;

  ArrayView<const T> view( size_type offset, size_type size ) const;

  ArrayView<const T> operator()( size_type offset, size_type size ) const;

  const ArrayView<const T>& operator()() const;

  /** \brief Return a const view of *this.
   *
   * This object is already const (this is a specialization for const
   * T), so this method is trivial; it just returns *this.
   */
  ArrayView<const T> getConst() const;

  iterator begin() const;

  iterator end() const;

  const ArrayView<const T>& assert_not_null() const;

  const ArrayView<const T>& assert_in_range( size_type offset, size_type size ) const;

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK

  // I should make these private but templated friends are not very
  // portable.  Besides, if a client directly calls this it will not
  // compile in an optimized build.

  explicit ArrayView (const ArrayRCP<const T> &arcp );

  explicit ArrayView (const T* p, size_type size, const ArrayRCP<const T> &arcp );

#endif

private:
  const T* ptr_;
  size_type size_;
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  ArrayRCP<const T> arcp_;
#endif

  void setUpIterators(const ERCPNodeLookup rcpNodeLookup = RCP_ENABLE_NODE_LOOKUP);

  void debug_assert_not_null() const {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    assert_not_null();
#endif
  }

  void debug_assert_in_range( size_type offset, size_type size_in ) const {
    (void)offset; (void)size_in;
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    assert_in_range(offset, size_in);
#endif
  }

  void debug_assert_valid_ptr() const {
#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
    arcp_.access_private_node().assert_valid_ptr(*this);
#endif
  }

public: // Bad bad bad

  // This is a very bad breach of encapsulation but it exists to avoid
  // problems with portability of templated friends
  const T* access_private_ptr() const { return ptr_; }

#ifdef HAVE_TEUCHOS_ARRAY_BOUNDSCHECK
  ArrayRCP<const T> access_private_arcp() const { return arcp_; }
#endif
};



/** \brief Construct a const or non-const view to const or non-const data.
 *
 * \relates ArrayView
 */
template<class T>
ArrayView<T> arrayView( T* p, typename ArrayView<T>::size_type size );


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


/** \brief Returns true if <tt>av.get()!=NULL</tt>.
 *
 * \relates ArrayView
 */
template<class T>
bool nonnull( const ArrayView<T> &av );


/** \brief Output stream inserter.
 *
 * The implementation of this function just prints pointer addresses and
 * therefore puts no restrictions on the data types involved.
 *
 * \relates ArrayView
 */
template<class T>
std::ostream& operator<<( std::ostream& out, const ArrayView<T>& av );


/** \brief Const cast of underlying <tt>ArrayView</tt> type from <tt>const T*</tt>
 * to <tt>T*</tt>.
 *
 * The function will compile only if (<tt>const_cast<T2*>(p1.get());</tt>)
 * compiles.
 *
 * \relates ArrayView
 */
template<class T2, class T1>
ArrayView<T2> av_const_cast(const ArrayView<T1>& p1);


/** \brief Reinterpret cast of underlying <tt>ArrayView</tt> type from
 * <tt>T1*</tt> to <tt>T2*</tt>.
 *
 * The function will compile only if (<tt>reinterpret_cast<T2*>(p1.get());</tt>) compiles.
 *
 * <b>Warning!</b> Do not use this function unless you absolutely know what
 * you are doing. Doing a reinterpret cast is always a tricking thing and
 * must only be done by developers who are 100% comfortable with what they are
 * doing.
 *
 * \relates ArrayView
 */
template<class T2, class T1>
ArrayView<T2> av_reinterpret_cast(const ArrayView<T1>& p1);


} // end namespace Teuchos


//
// Inline members
//


// ToDo: Fill in!


#endif  // TEUCHOS_ARRAY_VIEW_DECL_HPP
