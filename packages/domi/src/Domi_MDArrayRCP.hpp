// @HEADER
// ***********************************************************************
//
//     Domi: Multi-dimensional Distributed Linear Algebra Services
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef DOMI_MDARRAYRCP_HPP
#define DOMI_MDARRAYRCP_HPP

// Standard includes
#include <cstdarg>

// Teuchos includes
#include "Teuchos_ArrayRCPDecl.hpp"

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_Utils.hpp"
#include "Domi_MDArrayView.hpp"

namespace Domi
{

// I put these non-member template functions here for the same reason
// that Ross did the same thing for the Teuchos::Array class.  See
// Teuchos_Array.hpp for details.
template< typename T > class MDArrayRCP;

/** \brief Equality operator.
 *
 * \relates MDArrayRCP
 */
template< typename T >
bool operator==(const MDArrayRCP< T > & a1,
                const MDArrayRCP< T > & a2);

/** \brief MDArray/MDArrayRCP equality operator.
 *
 * \relates MDArray
 * \relates MDArrayRCP
 */
template< typename T >
bool operator==(const MDArray< T > & a1,
                const MDArrayRCP< T > & a2);

/** \brief MDArrayRCP/MDArray equality operator.
 *
 * \relates MDArray
 * \relates MDArrayRCP
 */
template< typename T >
bool operator==(const MDArrayRCP< T > & a1,
                const MDArray< T > & a2);

/** \brief MDArrayView/MDArrayRCP equality operator.
 *
 * \relates MDArrayView
 * \relates MDArrayRCP
 */
template< typename T >
bool operator==(const MDArrayView< T > & a1,
                const MDArrayRCP< T > & a2);

/** \brief MDArrayRCP/MDArrayView equality operator.
 *
 * \relates MDArrayView
 * \relates MDArrayRCP
 */
template< typename T >
bool operator==(const MDArrayRCP< T > & a1,
                const MDArrayView< T > & a2);

/** \brief Inequality operator.
 *
 * \relates MDArrayRCP
 */
template< typename T >
bool operator!=(const MDArrayRCP< T > & a1,
                const MDArrayRCP< T > & a2);

/** \brief MDArray/MDArrayRCP inequality operator.
 *
 * \relates MDArray
 * \relates MDArrayRCP
 */
template< typename T >
bool operator!=(const MDArray< T > & a1,
                const MDArrayRCP< T > & a2);

/** \brief MDArrayView/MDArrayRCP inequality operator.
 *
 * \relates MDArrayView
 * \relates MDArrayRCP
 */
template< typename T >
bool operator!=(const MDArrayView< T > & a1,
                const MDArrayRCP< T > & a2);

/** \brief MDArrayRCP/MDArrayView inequality operator.
 *
 * \relates MDArrayView
 * \relates MDArrayRCP
 */
template< typename T >
bool operator!=(const MDArrayRCP< T > & a1,
                const MDArrayView< T > & a2);

/** \brief Memory-safe, reference-counted, templated,
 * multi-dimensional array class
 *
 * <tt>MDArrayRCP</tt> is a reference counted class similar to
 * <tt>Teuchos::ArrayRCP</tt>, except that it interprets its data
 * buffer as a multi-dimensional array instead of a one-dimensional
 * array.
 */

template< typename T >
class MDArrayRCP
{
public:

  /** \name Public types */
  //@{

  /** \brief Value type */
  typedef T value_type;

  /** \brief Pointer type */
  typedef T* pointer;

  /** \brief Const pointer type */
  typedef const T* const_pointer;

  /** \brief Reference type */
  typedef T& reference;

  /** \brief Const reference type */
  typedef const T& const_reference;

  //@}

  /** \name Constructors, Destructor, Initializers */
  //@{

  /** \brief Default constructor
   *
   * \param null_arg [in] Optional null pointer specifier
   */
  inline
  MDArrayRCP(Teuchos::ENull null_arg = Teuchos::null);

  /** \brief Constructor with <tt>Teuchos::ArrayView</tt> source,
   *  dimensions, and optional storage order flag.
   *
   * \param array [in] <tt>Teuchos::ArrayView</tt> of data buffer
   *
   * \param dims [in] An array that defines the lengths of each
   *        dimension.  The most convenient way to specify dimensions
   *        is with a Tuple returned by the non-member
   *        <tt>Teuchos::tuple<T>()</tt> function.
   *
   * \param layout [in] Specifies the order data elements are stored
   *        in memory (default DEFAULT_ORDER)
   *
   * The <tt>MDArrayRCP</tt> does not take ownership of the data
   * buffer.
   */
  inline
  MDArrayRCP(const Teuchos::ArrayView< T > & array,
             const Teuchos::ArrayView< dim_type > & dims,
             Layout layout = DEFAULT_ORDER);

  /** \brief Constructor with dimensions, default value and optional
   *  storage order flag.
   *
   * \param dims [in] An array that defines the lengths of each
   *        dimension.  The most convenient way to specify dimensions
   *        is with a Tuple returned by the non-member
   *        <tt>Teuchos::tuple<T>()</tt> function.
   *
   * \param val [in] Default array fill value
   *
   * \param layout [in] Specifies the order data elements are stored
   *        in memory
   *
   * This constructor allocates new memory and takes ownership of it.
   */
  inline explicit
  MDArrayRCP(const Teuchos::ArrayView< dim_type > & dims,
             const_reference val = T(),
             Layout layout = DEFAULT_ORDER);

  /** \brief Constructor with dimensions and storage order flag.
   *
   * \param dims [in] An array that defines the lengths of each
   *        dimension.  The most convenient way to specify dimensions
   *        is with a Tuple returned by the non-member
   *        <tt>Teuchos::tuple<T>()</tt> function.
   *
   * \param layout [in] Specifies the order data elements are stored
   *        in memory
   *
   * This constructor allocates new memory and takes ownership of it.
   */
  inline explicit
  MDArrayRCP(const Teuchos::ArrayView< dim_type > & dims,
             Layout layout);

  /** \brief Low-level view constructor
   *
   * \param dims [in] An array that defines the lengths of each
   *        dimension.  The most convenient way to specify dimensions
   *        is with a Tuple returned by the non-member
   *        <tt>Teuchos::tuple<T>()</tt> function.
   *
   * \param strides [in] An array that defines the data strides of
   *        each dimension.  The most convenient way to specify
   *        strides is with a Tuple returned by the non-member
   *        <tt>Teuchos::tuple<T>()</tt> function.
   *
   * \param data [in] a pointer to the data buffer used by the
   *        MDArrayRCP.  The MDArrayRCP will not take ownership of the
   *        data.
   *
   * \param layout [in] Specifies the order data elements are stored
   *        in memory
   *
   */
  inline explicit
  MDArrayRCP(const Teuchos::ArrayView< dim_type > & dims,
             const Teuchos::ArrayView< size_type > & strides,
             T * data,
             Layout layout = DEFAULT_ORDER);

  /** \brief Shallow copy constructor
   *
   * \param r_ptr [in] Source reference counted pointer
   */
  inline
  MDArrayRCP(const MDArrayRCP< T > & r_ptr);

  /** \brief Deep copy constructor
   *
   * \param source [in] Source MDArrayView
   */
  MDArrayRCP(const MDArrayView< T > & source);

  /** \brief Destructor
   */
  inline ~MDArrayRCP();

  /** \brief Assignment operator
   *
   * \param r_ptr [in] Source reference counted pointer
   */
  inline MDArrayRCP< T > & operator=(const MDArrayRCP< T > & r_ptr);

  //@}

  /** \name Attribute accessor methods */
  //@{

  /** \brief Return the number of dimensions
   */
  inline int numDims() const;

  /** \brief Return the array of dimensions
   */
  inline const Teuchos::Array< dim_type > & dimensions() const;

  /** \brief Return the dimension of the given axis
   *
   * \param axis [in] The axis being queried (0 for the first axis,
   *        1 for the second axis, and so forth)
   */
  inline dim_type dimension(int axis) const;

  /** \brief Return the total size of the <tt>MDArrayRCP</tt>
   */
  inline size_type size() const;

  /** \brief Return the indexing strides
   */
  inline const Teuchos::Array< size_type > & strides() const;

  /** \brief Return the underlying <tt>Teuchos::ArrayRCP</tt>
   */
  inline const Teuchos::ArrayRCP< T > & arrayRCP() const;

  /** \brief Return the storage order
   */
  inline const Layout layout() const;

  //@}

  /** \name Iterator classes and methods */
  //@{

  friend class MDIterator< MDArrayRCP< T > >;
  friend class MDIterator< MDArrayRCP< const T > >;
  friend class MDRevIterator< MDArrayRCP< T > >;
  friend class MDRevIterator< MDArrayRCP< const T > >;

  typedef MDIterator< MDArrayRCP< T > >          iterator;
  typedef MDIterator< MDArrayRCP< const T > >    const_iterator;
  typedef MDRevIterator< MDArrayRCP< T > >       reverse_iterator;
  typedef MDRevIterator< MDArrayRCP< const T > > const_reverse_iterator;

  /** \brief Return the beginning iterator
   */
  iterator begin();

  /** \brief Return the ending iterator
   */
  iterator end();

  /** \brief Return the beginning const_iterator
   */
  const_iterator begin() const;

  /** \brief Return the ending const_iterator
   */
  const_iterator end() const;

  /** \brief Return the beginning const_iterator
   */
  const_iterator cbegin() const;

  /** \brief Return the ending const_iterator
   */
  const_iterator cend() const;

  /** \brief Return the beginning reverse_iterator
   */
  reverse_iterator rbegin();

  /** \brief Return the ending reverse_iterator
   */
  reverse_iterator rend();

  /** \brief Return the beginning const_reverse_iterator
   */
  const_reverse_iterator crbegin() const;

  /** \brief Return the ending const_reverse_iterator
   */
  const_reverse_iterator crend() const;

  //@}

  /** \name Object/pointer access functions */
  //@{

  /** \brief Return true if the underlying pointer is null
   */
  inline bool is_null() const;

  /** \brief Pointer <tt>-></tt> access to members of underlying data buffer
   */
  inline pointer operator->() const;

  /** \brief Dereference the underlying data buffer
   */
  inline reference operator*();

  /** \brief Get the raw C++ pointer to the underlying data buffer
   */
  inline pointer get() const;

  //@}

  /** \name Conversions to MDArrayView */
  //@{

  /** \brief Perform an explicit conversion to a non-const
   *  <tt>MDArrayView<T></tt>
   */
  MDArrayView< T > mdArrayView();

  /** \brief Perform an explicit conversion to a const
   *  <tt>MDArrayView<T></tt>
   */
  const MDArrayView< T > mdArrayView() const;

  /** \brief Perform an explicit conversion to a non-const
   *  <tt>MDArrayView<const T></tt>
   */
  MDArrayView< const T > mdArrayViewConst();

  /** \brief Perform an explicit conversion to a const
   *  <tt>MDArrayView<const T></tt>
   */
  const MDArrayView< const T > mdArrayViewConst() const;

  /** \brief Perform an implicit conversion to a non-const
   *  <tt>MDArrayView</tt>
   */
  inline operator MDArrayView< T >() const;

  /** \brief Perform an implicit conversion to a const
   *  <tt>MDArrayView</tt>
   */
  inline operator MDArrayView< const T >() const;

  //@}

  /** \name Indexing operators that return <tt>MDArrayView</tt>s */
  //@{

  /** \brief Sub-array access operator.  The returned
   *  <tt>MDArrayView</tt> object will have one fewer dimensions than
   *  the calling <tt>MDArrayRCP</tt>.
   *
   * \param i [in] Index of the desired sub-array.  Note that to
   *        obtain expected behavior, you should always chain together
   *        <tt>n</tt> square bracket operators when referencing an
   *        <tt>n</tt>-dimensional <tt>MDArrayRCP</tt>.
   */
  MDArrayView< T > operator[](dim_type i);

  /** \brief Sub-array const access operator.  The returned
   *  <tt>MDArrayView</tt> object will have one fewer dimensions than
   *  the calling <tt>MDArrayRCP</tt>.
   *
   * \param i [in] Index of the desired sub-array.  Note that to
   *        obtain expected behavior, you should always chain together
   *        <tt>n</tt> square bracket operators when referencing an
   *        <tt>n</tt>-dimensional <tt>MDArrayRCP</tt>.
   */
  const MDArrayView< T > operator[](dim_type i) const;

  /** \brief Sub-array access operator.  The returned
   *  <tt>MDArrayView</tt> object will have the same number of
   *  dimensions as the calling <tt>MDArrayRCP</tt>.
   *
   * \param s [in] Slice representing the bounds of the desired
   *        sub-array.  Note that to obtain expected behavior, you
   *        should always chain together <tt>n</tt> square bracket
   *        operators when referencing an <tt>n</tt>-dimensional
   *        <tt>MDArrayRCP</tt>.
   */
  MDArrayView< T > operator[](Slice s);

  /** \brief Sub-array const access operator.  The returned
   *  <tt>MDArrayView</tt> object will have the same number of
   *  dimensions as the calling <tt>MDArrayRCP</tt>.
   *
   * \param s [in] Slice representing the bounds of the desired
   *        sub-array.  Note that to obtain expected behavior, you
   *        should always chain together <tt>n</tt> square bracket
   *        operators when referencing an <tt>n</tt>-dimensional
   *        <tt>MDArrayRCP</tt>.
   */
  const MDArrayView< T > operator[](Slice s) const;

  /** \brief Conversion to <tt>MDArrayView</tt>
   */
  MDArrayView< T > operator()();

  /** \brief Conversion to const <tt>MDArrayView</tt>
   */
  const MDArrayView< T > operator()() const;

  //@}

  /** \name Indexing operators that return a reference to a single array element */
  //@{

  /** \brief Non-const 1D element access operator
   *
   * \param i [in] 1D index.
   *
   * This operator should only be used with a 1D <tt>MDArrayRCP</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the
   * <tt>MDArrayRCP</tt> is not 1D, an exception will be thrown.
   */
  inline T & operator()(dim_type i);

  /** \brief Non-const 2D element access operator
   *
   * \param i [in] first 2D index.
   *
   * \param j [in] second 2D index.
   *
   * This operator should only be used with a 2D <tt>MDArrayRCP</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the
   * <tt>MDArrayRCP</tt> is not 2D, an exception will be thrown.
   */
  inline T & operator()(dim_type i, dim_type j);

  /** \brief Non-const 3D element access operator
   *
   * \param i [in] first 3D index.
   *
   * \param j [in] second 3D index.
   *
   * \param k [in] third 3D index.
   *
   * This operator should only be used with a 3D <tt>MDArrayRCP</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the
   * <tt>MDArrayRCP</tt> is not 3D, an exception will be thrown.
   */
  inline T & operator()(dim_type i, dim_type j, dim_type k);

  /** \brief Non-const 4D element access operator
   *
   * \param i [in] first 4D index.
   *
   * \param j [in] second 4D index.
   *
   * \param k [in] third 4D index.
   *
   * \param m [in] fourth 4D index.
   *
   * This operator should only be used with a 4D <tt>MDArrayRCP</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the
   * <tt>MDArrayRCP</tt> is not 4D, an exception will be thrown.
   */
  inline T & operator()(dim_type i, dim_type j, dim_type k, dim_type m);

  /** \brief Non-const 5D element access operator
   *
   * \param i [in] first 5D index.
   *
   * \param j [in] second 5D index.
   *
   * \param k [in] third 5D index.
   *
   * \param m [in] fourth 5D index.
   *
   * \param n [in] fifth 5D index.
   *
   * This operator should only be used with a 5D <tt>MDArrayRCP</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the
   * <tt>MDArrayRCP</tt> is not 5D, an exception will be thrown.
   */
  inline T & operator()(dim_type i, dim_type j, dim_type k, dim_type m,
                        dim_type n);

  /** \brief Non-const 6D and higher element access operator
   *
   * \param i [in] first index.
   *
   * \param j [in] second index.
   *
   * \param k [in] third index.
   *
   * \param m [in] fourth index.
   *
   * \param n [in] fifth index.
   *
   * \param p [in] sixth index.
   *
   * \param ... [in] seventh and higher indexes.
   *
   * This operator should only be used with a 6D and higher
   * <tt>MDArrayRCP</tt>s.  If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and
   * the <tt>MDArrayRCP</tt> is less than 6D, an exception will be
   * thrown.
   */
  inline T & operator()(dim_type i, dim_type j, dim_type k, dim_type m,
                        dim_type n, dim_type p, ...);

  /** \brief Const 1D element access operator
   *
   * \param i [in] 1D index.
   *
   * This operator should only be used with a 1D <tt>MDArrayRCP</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the
   * <tt>MDArrayRCP</tt> is not 1D, an exception will be thrown.
   */
  inline const T & operator()(dim_type i) const;

  /** \brief Const 2D element access operator
   *
   * \param i [in] first 2D index.
   *
   * \param j [in] second 2D index.
   *
   * This operator should only be used with a 2D <tt>MDArrayRCP</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the
   * <tt>MDArrayRCP</tt> is not 2D, an exception will be thrown.
   */
  inline const T & operator()(dim_type i, dim_type j) const;

  /** \brief Const 3D element access operator
   *
   * \param i [in] first 3D index.
   *
   * \param j [in] second 3D index.
   *
   * \param k [in] third 3D index.
   *
   * This operator should only be used with a 3D <tt>MDArrayRCP</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the
   * <tt>MDArrayRCP</tt> is not 3D, an exception will be thrown.
   */
  inline const T & operator()(dim_type i, dim_type j, dim_type k) const;

  /** \brief Const 4D element access operator
   *
   * \param i [in] first 4D index.
   *
   * \param j [in] second 4D index.
   *
   * \param k [in] third 4D index.
   *
   * \param m [in] fourth 4D index.
   *
   * This operator should only be used with a 4D <tt>MDArrayRCP</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the
   * <tt>MDArrayRCP</tt> is not 4D, an exception will be thrown.
   */
  inline const T & operator()(dim_type i, dim_type j, dim_type k,
                              dim_type m) const;

  /** \brief Const 5D element access operator
   *
   * \param i [in] first 5D index.
   *
   * \param j [in] second 5D index.
   *
   * \param k [in] third 5D index.
   *
   * \param m [in] fourth 5D index.
   *
   * \param n [in] fifth 5D index.
   *
   * This operator should only be used with a 5D <tt>MDArrayRCP</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the
   * <tt>MDArrayRCP</tt> is not 5D, an exception will be thrown.
   */
  inline const T & operator()(dim_type i, dim_type j, dim_type k,
                              dim_type m, dim_type n) const;

  /** \brief Const 6D and higher element access operator
   *
   * \param i [in] first index.
   *
   * \param j [in] second index.
   *
   * \param k [in] third index.
   *
   * \param m [in] fourth index.
   *
   * \param n [in] fifth index.
   *
   * \param p [in] sixth index.
   *
   * \param ... [in] seventh and higher indexes.
   *
   * This operator should only be used with a 6D and higher
   * <tt>MDArrayRCP</tt>s.  If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and
   * the <tt>MDArrayRCP</tt> is less than 6D, an exception will be
   * thrown.
   */
  inline const T & operator()(dim_type i, dim_type j, dim_type k,
                              dim_type m, dim_type n, dim_type p, ...) const;


  //@}

  /** \name Teuchos::Array-like and std::vector-like methods */
  //@{

  /** \brief Assign a value to all elements of the <tt>MDArrayRCP</tt>
   *
   * \param value [in] The value to be assigned
   */
  void assign(const_reference value);

  /** \brief Non-const single element access method with bounds
   * checking
   *
   * \param i, ... [in] Indexes representing the location of the
   *        single element of the <tt>MDArrayRCP</tt> to be accessed.
   *        Note that this method assumes that the user will provide
   *        the same number of arguments as the number of dimensions
   *        of the <tt>MDArrayRCP</tt>.
   */
  reference at(dim_type i, ...);

  /** \brief Const single element access method with bounds checking
   *
   * \param i, ... [in] Indexes representing the location of the
   *        single element of the <tt>MDArrayRCP</tt> to be accessed.
   *        Note that this method assumes that the user will provide
   *        the same number of arguments as the number of dimensions
   *        of the <tt>MDArrayRCP</tt>.
   */
  const_reference at(dim_type i, ...) const;

  /** \brief Return the capacity of the underlying <tt>Teuchos::ArrayRCP</tt>
   */
  inline size_type capacity() const;

  /** \brief Clear the <tt>MDArrayRCP</tt>
   */
  void clear();

  /** \brief Return whether the <tt>MDArrayRCP</tt> is empty
   */
  inline bool empty() const;

  /** \brief Return the maximum allowable size for the
   *  <tt>MDArrayRCP</tt>
   */
  inline size_type max_size() const;

  /** \brief Resize the <tt>MDArrayRCP</tt> based on the given
   * dimensions
   *
   * \param dims [in] An array that defines the new lengths of each
   *        dimension.  The most convenient way to specify dimensions
   *        is with a Tuple returned by the non-member
   *        <tt>Teuchos::tuple<T>()</tt> function.
   */
  void resize(const Teuchos::ArrayView< dim_type > & dims);

  /** \brief Return true if <tt>MDArrayRCP</tt> has been compiled with
   *  bounds checking on.
   */
  inline static bool hasBoundsChecking();

  /** \brief Convert the <tt>MDArrayRCP</tt> to a string
   *  representation
   */
  std::string toString() const;

  /** \brief Return a const raw pointer to the beginning of the
   *  <tt>MDArrayRCP</tt> or NULL if unsized.
   */
  inline const_pointer getRawPtr() const;

  /** \brief Return a raw pointer to the beginning of the
   *  <tt>MDArrayRCP</tt> or NULL if unsized.
   */
  inline pointer getRawPtr();

  //@}

  // These operators are declared as friends so that the compiler will
  // do automatic type conversion.

  /** \name Non-member operators and functions */
  //@{

  /** \brief Equality operator.
   */
  template< typename T2 >
  friend bool operator==(const MDArrayRCP< T2 > & a1,
                         const MDArrayRCP< T2 > & a2);

  /** \brief MDArray/MDArrayRCP equality operator.
   */
  template< typename T2 >
  friend bool operator==(const MDArray< T2 > & a1,
                         const MDArrayRCP< T2 > & a2);

  /** \brief MDArrayRCP/MDArray equality operator.
   */
  template< typename T2 >
  friend bool operator==(const MDArrayRCP< T2 > & a1,
                         const MDArray< T2 > & a2);

  /** \brief MDArrayRCP/MDArrayView equality operator.
   */
  template< typename T2 >
  friend bool operator==(const MDArrayRCP< T2 > & a1,
                         const MDArrayView< T2 > & a2);

  /** \brief MDArrayView/MDArrayRCP equality operator.
   */
  template< typename T2 >
  friend bool operator==(const MDArrayView< T2 > & a1,
                         const MDArrayRCP< T2 > & a2);

  /** \brief Inequality operator.
   */
  template< typename T2 >
  friend bool operator!=(const MDArrayRCP< T2 > & a1,
                         const MDArrayRCP< T2 > & a2);

  /** \brief MDArray/MDArrayRCP inequality operator.
   */
  template< typename T2 >
  friend bool operator!=(const MDArray< T2 > & a1,
                         const MDArrayRCP< T2 > & a2);

  /** \brief MDArrayRCP/MDArray inequality operator.
   */
  template< typename T2 >
  friend bool operator!=(const MDArrayRCP< T2 > & a1,
                         const MDArray< T2 > & a2);

  /** \brief MDArrayRCP/MDArrayView inequality operator.
   */
  template< typename T2 >
  friend bool operator!=(const MDArrayRCP< T2 > & a1,
                         const MDArrayView< T2 > & a2);

  /** \brief MDArrayView/MDArrayRCP inequality operator.
   */
  template< typename T2 >
  friend bool operator!=(const MDArrayView< T2 > & a1,
                         const MDArrayRCP< T2 > & a2);

  /** \brief Stream output operator
   */
  template< typename T2 >
  friend std::ostream & operator<<(std::ostream & os,
                                   const MDArrayRCP< T2 > & a);

  //@}

private:
  Teuchos::Array< dim_type >  _dimensions;
  Teuchos::Array< size_type > _strides;
  Teuchos::ArrayRCP< T >      _array;
  Layout                      _layout;
  pointer                     _ptr;

  // Used for array bounds checking
  void assertAxis(int axis) const;

  // Used for array bounds checking
  void assertIndex(dim_type i, int axis) const;
};

////////////////////////////////////////////////////////////////////////

/////////////////////
// Implementations //
/////////////////////

template< typename T >
MDArrayRCP< T >::MDArrayRCP(Teuchos::ENull null_arg) :
  _dimensions(Teuchos::tuple< dim_type >(0)),
  _strides(Teuchos::tuple< size_type >(1)), 
  _array(),
  _layout(DEFAULT_ORDER),
  _ptr()
{
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayRCP< T >::MDArrayRCP(const Teuchos::ArrayView< T > & array,
			    const Teuchos::ArrayView< dim_type > & dims,
			    Layout layout) :
  _dimensions(dims),
  _strides(computeStrides< size_type, dim_type >(dims, layout)),
  _array(array.getRawPtr(), 0, array.size(), false),
  _layout(layout),
  _ptr(_array.getRawPtr())
{
  TEUCHOS_TEST_FOR_EXCEPTION(array.size() < computeSize(dims),
			     RangeError,
			     "Teuchos::ArrayView size too small for "
                             "dimensions");
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayRCP< T >::MDArrayRCP(const Teuchos::ArrayView< dim_type > & dims,
			    const T & val,
			    Layout layout) :
  _dimensions(dims),
  _strides(computeStrides< size_type, dim_type >(dims, layout)),
  _array(computeSize(dims), val),
  _layout(layout),
  _ptr(_array.getRawPtr())
{
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayRCP< T >::MDArrayRCP(const Teuchos::ArrayView< dim_type > & dims,
			    Layout layout) :
  _dimensions(dims),
  _strides(computeStrides< size_type, dim_type >(dims, layout)),
  _array(computeSize(dims)),
  _layout(layout),
  _ptr(_array.getRawPtr())
{
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayRCP< T >::MDArrayRCP(const Teuchos::ArrayView< dim_type > & dims,
                            const Teuchos::ArrayView< size_type > & strides,
                            T * data,
                            Layout layout) :
  _dimensions(dims),
  _strides(strides),
  _array(data, 0, computeSize(dims, strides), false),
  _layout(layout),
  _ptr(data)
{
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayRCP< T >::MDArrayRCP(const MDArrayRCP< T > & r_ptr) :
  _dimensions(r_ptr._dimensions),
  _strides(r_ptr._strides),
  _array(r_ptr._array),
  _layout(r_ptr._layout),
  _ptr(_array.getRawPtr())
{
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayRCP< T >::MDArrayRCP(const MDArrayView< T > & source) :
  _dimensions(source.dimensions()),
  _strides(computeStrides< size_type, dim_type >(source.dimensions(),
                                                 source.layout())),
  _array(computeSize(source.dimensions())),
  _layout(source.layout()),
  _ptr(_array.getRawPtr())
{
  // Copy the values from the MDArrayView to the MDArrayRCP
  iterator thisit = begin();
  typename MDArrayView< T >::const_iterator srcit = source.cbegin();
  for ( ; srcit != source.cend(); ++thisit, ++srcit)
  {
    *thisit = *srcit;
  }
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayRCP< T >::~MDArrayRCP()
{
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayRCP< T > &
MDArrayRCP< T >::operator=(const MDArrayRCP< T > & r_ptr)
{
  _dimensions = r_ptr._dimensions;
  _strides    = r_ptr._strides;
  _array      = r_ptr._array;
  _layout     = r_ptr._layout;
  _ptr        = r_ptr._ptr;
  return *this;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
int
MDArrayRCP< T >::numDims() const
{
  return _dimensions.size();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const Teuchos::Array< dim_type > &
MDArrayRCP< T >::dimensions() const
{
  return _dimensions;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
dim_type
MDArrayRCP< T >::dimension(int axis) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  assertAxis(axis);
#endif
  return _dimensions[axis];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
size_type
MDArrayRCP< T >::size() const
{
  return _array.size();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const Teuchos::Array< size_type > &
MDArrayRCP< T >::strides() const
{
  return _strides;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const Teuchos::ArrayRCP< T > &
MDArrayRCP< T >::arrayRCP() const
{
  return _array;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const Layout
MDArrayRCP< T >::layout() const
{
  return _layout;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayRCP< T >::iterator
MDArrayRCP< T >::begin()
{
  return iterator(*this);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayRCP< T >::iterator
MDArrayRCP< T >::end()
{
  // Return the iterator corresponding to the last element
  return iterator(*this, true);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayRCP< T >::const_iterator
MDArrayRCP< T >::begin() const
{
  return const_iterator(*this);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayRCP< T >::const_iterator
MDArrayRCP< T >::end() const
{
  // Return the iterator corresponding to the last element
  return const_iterator(*this, true);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayRCP< T >::const_iterator
MDArrayRCP< T >::cbegin() const
{
  return const_iterator(*this);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayRCP< T >::const_iterator
MDArrayRCP< T >::cend() const
{
  // Return the iterator corresponding to the last element
  return const_iterator(*this, true);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayRCP< T >::reverse_iterator
MDArrayRCP< T >::rbegin()
{
  return reverse_iterator(*this);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayRCP< T >::reverse_iterator
MDArrayRCP< T >::rend()
{
  // Return the reverse_iterator corresponding to the last element
  return reverse_iterator(*this, true);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayRCP< T >::const_reverse_iterator
MDArrayRCP< T >::crbegin() const
{
  return const_reverse_iterator(*this);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayRCP< T >::const_reverse_iterator
MDArrayRCP< T >::crend() const
{
  // Return the reverse_iterator corresponding to the last element
  return const_reverse_iterator(*this, true);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool
MDArrayRCP< T >::is_null() const
{
  return _array.is_null();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
T *
MDArrayRCP< T >::operator->() const
{
  return _array.operator->();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
T &
MDArrayRCP< T >::operator*()
{
  return _array.operator*();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
T *
MDArrayRCP< T >::get() const
{
  return _array.get();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayView< T >
MDArrayRCP< T >::mdArrayView()
{
  return MDArrayView< T >(_array(), _dimensions, _layout);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const MDArrayView< T >
MDArrayRCP< T >::mdArrayView() const
{
  Teuchos::Array< dim_type > dims(_dimensions);
  return MDArrayView< T >(_array(), dims(), _layout);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayView< const T >
MDArrayRCP< T >::mdArrayViewConst()
{
  return MDArrayView< const T >(_array(), _dimensions, _layout);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const MDArrayView< const T >
MDArrayRCP< T >::mdArrayViewConst() const
{
  return MDArrayView< const T >(_array(), _dimensions, _layout);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayRCP< T >::operator MDArrayView< T >() const
{
  return mdArrayView();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayRCP< T >::operator MDArrayView< const T >() const
{
  return mdArrayViewConst();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayView< T >
MDArrayRCP< T >::operator[](dim_type i)
{
  // Note: array bounds checking, if active, will be performed by the
  // MDArrayView class
  return mdArrayView()[i];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const MDArrayView< T >
MDArrayRCP< T >::operator[](dim_type i) const
{
  // Note: array bounds checking, if active, will be performed by the
  // MDArrayView class
  return mdArrayView()[i];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const MDArrayView< T >
MDArrayRCP< T >::operator[](Slice s) const
{
  // Note: Slices produce safe indexes
  return mdArrayView()[s];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayView< T >
MDArrayRCP< T >::operator[](Slice s)
{
  // Note: Slices produce safe indexes
  return mdArrayView()[s];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayView< T >
MDArrayRCP< T >::operator()()
{
  return mdArrayView();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const MDArrayView< T >
MDArrayRCP< T >::operator()() const
{
  return mdArrayView();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
T &
MDArrayRCP< T >::operator()(dim_type i)
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() != 1), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with 1 index"
    );
  assertIndex(i, 0);
#endif
  return _ptr[i * _strides[0]];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
T &
MDArrayRCP< T >::operator()(dim_type i,
                            dim_type j)
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() != 2), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with 2 indexes"
    );
  assertIndex(i, 0);
  assertIndex(j, 1);
#endif
  return _ptr[i * _strides[0] + j * _strides[1]];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
T &
MDArrayRCP< T >::operator()(dim_type i,
                            dim_type j,
                            dim_type k)
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() != 3), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with 3 indexes"
    );
  assertIndex(i, 0);
  assertIndex(j, 1);
  assertIndex(k, 2);
#endif
  return _ptr[i * _strides[0] + j * _strides[1] + k * _strides[2]];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
T &
MDArrayRCP< T >::operator()(dim_type i,
                            dim_type j,
                            dim_type k,
                            dim_type m)
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() != 4), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with 4 indexes"
    );
  assertIndex(i, 0);
  assertIndex(j, 1);
  assertIndex(k, 2);
  assertIndex(m, 3);
#endif
  return _ptr[i * _strides[0] + j * _strides[1] + k * _strides[2] +
              m * _strides[3]];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
T &
MDArrayRCP< T >::operator()(dim_type i,
                            dim_type j,
                            dim_type k,
                            dim_type m,
                            dim_type n)
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() != 5), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with 5 indexes"
    );
  assertIndex(i, 0);
  assertIndex(j, 1);
  assertIndex(k, 2);
  assertIndex(m, 3);
  assertIndex(n, 4);
#endif
  return _ptr[i * _strides[0] + j * _strides[1] + k * _strides[2] +
              m * _strides[3] + n * _strides[4]];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
T &
MDArrayRCP< T >::operator()(dim_type i,
                            dim_type j,
                            dim_type k,
                            dim_type m,
                            dim_type n,
                            dim_type p,
                            ...)
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() < 6), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with too many indexes"
    );
  assertIndex(i, 0);
  assertIndex(j, 1);
  assertIndex(k, 2);
  assertIndex(m, 3);
  assertIndex(n, 4);
  assertIndex(p, 5);
#endif
  va_list indexes;
  size_type offset = i * _strides[0] + j * _strides[1] + k * _strides[2] +
                     m * _strides[3] + n * _strides[4] + p * _strides[5];
  va_start(indexes, p);
  for (int axis = 6; axis < _dimensions.size(); axis++)
  {
    dim_type q = va_arg(indexes, dim_type);
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
    assertIndex(q, axis);
#endif
    offset += q * _strides[axis];
  }
  va_end(indexes);
  return _ptr[offset];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const T &
MDArrayRCP< T >::operator()(dim_type i) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() != 1), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with 1 index"
    );
  assertIndex(i, 0);
#endif
  return _ptr[i * _strides[0]];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const T &
MDArrayRCP< T >::operator()(dim_type i,
                            dim_type j) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() != 2), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with 2 indexes"
    );
  assertIndex(i, 0);
  assertIndex(j, 1);
#endif
  return _ptr[i * _strides[0] + j * _strides[1]];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const T &
MDArrayRCP< T >::operator()(dim_type i,
                            dim_type j,
                            dim_type k) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() != 3), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with 3 indexes"
    );
  assertIndex(i, 0);
  assertIndex(j, 1);
  assertIndex(k, 2);
#endif
  return _ptr[i * _strides[0] + j * _strides[1] + k * _strides[2]];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const T &
MDArrayRCP< T >::operator()(dim_type i,
                            dim_type j,
                            dim_type k,
                            dim_type m) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() != 4), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with 4 indexes"
    );
  assertIndex(i, 0);
  assertIndex(j, 1);
  assertIndex(k, 2);
  assertIndex(m, 3);
#endif
  return _ptr[i * _strides[0] + j * _strides[1] + k * _strides[2] +
              m * _strides[3]];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const T &
MDArrayRCP< T >::operator()(dim_type i,
                            dim_type j,
                            dim_type k,
                            dim_type m,
                            dim_type n) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() != 5), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with 5 indexes"
    );
  assertIndex(i, 0);
  assertIndex(j, 1);
  assertIndex(k, 2);
  assertIndex(m, 3);
  assertIndex(n, 4);
#endif
  return _ptr[i * _strides[0] + j * _strides[1] + k * _strides[2] +
              m * _strides[3] + n * _strides[4]];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const T &
MDArrayRCP< T >::operator()(dim_type i,
                            dim_type j,
                            dim_type k,
                            dim_type m,
                            dim_type n,
                            dim_type p,
                            ...) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() < 6), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with too many indexes"
    );
  assertIndex(i, 0);
  assertIndex(j, 1);
  assertIndex(k, 2);
  assertIndex(m, 3);
  assertIndex(n, 4);
  assertIndex(p, 5);
#endif
  va_list indexes;
  size_type offset = i * _strides[0] + j * _strides[1] + k * _strides[2] +
                     m * _strides[3] + n * _strides[4] + p * _strides[5];
  va_start(indexes, p);
  for (int axis = 6; axis < _dimensions.size(); axis++)
  {
    dim_type q = va_arg(indexes, dim_type);
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
    assertIndex(q, axis);
#endif
    offset += q * _strides[axis];
  }
  va_end(indexes);
  return _ptr[offset];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
void
MDArrayRCP< T >::assign(const_reference value)
{
  for (iterator it = begin(); it != end(); ++it)
    *it = value;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
T &
MDArrayRCP< T >::at(dim_type i, ...)
{
  assertIndex(i, 0);
  va_list indexes;
  size_type offset = i * _strides[0];
  va_start(indexes, i);
  for (int axis = 1; axis < _dimensions.size(); axis++)
  {
    dim_type j = va_arg(indexes, dim_type);
    assertIndex(j, axis);
    offset += j * _strides[axis];
  }
  va_end(indexes);
  return _array[offset];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const T &
MDArrayRCP< T >::at(dim_type i, ...) const
{
  assertIndex(i, 0);
  va_list indexes;
  size_type offset = i * _strides[0];
  va_start(indexes, i);
  for (int axis = 1; axis < _dimensions.size(); axis++)
  {
    dim_type j = va_arg(indexes, dim_type);
    assertIndex(j, axis);
    offset += j * _strides[axis];
  }
  va_end(indexes);
  return _array[offset];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
size_type
MDArrayRCP< T >::capacity() const
{
  return _array.capacity();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
void
MDArrayRCP< T >::clear()
{
  _dimensions.resize(1);
  _dimensions[0] = 0;
  _strides.resize(1);
  _strides[0] = 1;
  _array.clear();
  _ptr = _array.getRawPtr();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool
MDArrayRCP< T >::empty() const
{
  return (_array.size() == 0);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
size_type
MDArrayRCP< T >::max_size() const
{
  return _array.max_size();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
void
MDArrayRCP< T >::resize(const Teuchos::ArrayView< dim_type > & dims)
{
  _dimensions.assign(dims.begin(), dims.end());
  _strides = computeStrides< size_type, dim_type >(dims, _layout);
  _array.resize(computeSize(dims));
  _ptr = _array.getRawPtr();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool
MDArrayRCP< T >::hasBoundsChecking()
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  return true;
#else
  return false;
#endif
}

////////////////////////////////////////////////////////////////////////

template< typename T >
std::string
MDArrayRCP< T >::toString() const
{
  return mdArrayView().toString();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const T *
MDArrayRCP< T >::getRawPtr() const
{
  return _array.getRawPtr();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
T *
MDArrayRCP< T >::getRawPtr()
{
  return _array.getRawPtr();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool operator==(const MDArrayRCP< T > & a1,
                const MDArrayRCP< T > & a2)
{
  return (a1() == a2());
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool operator==(const MDArray< T > & a1,
                const MDArrayRCP< T > & a2)
{
  return (a1() == a2());
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool operator==(const MDArrayRCP< T > & a1,
                const MDArray< T > & a2)
{
  return (a1() == a2());
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool operator==(const MDArrayView< T > & a1,
                const MDArrayRCP< T > & a2)
{
  return (a1 == a2());
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool operator==(const MDArrayRCP< T > & a1,
                const MDArrayView< T > & a2)
{
  return (a1() == a2);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool operator!=(const MDArrayRCP< T > & a1,
                const MDArrayRCP< T > & a2)
{
  return not (a1 == a2);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool operator!=(const MDArray< T > & a1, const MDArrayRCP< T > & a2)
{
  return (a1() != a2());
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool operator!=(const MDArrayRCP< T > & a1, const MDArray< T > & a2)
{
  return (a1() != a2());
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool operator!=(const MDArrayView< T > & a1,
                const MDArrayRCP< T > & a2)
{
  return (a1 != a2());
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool operator!=(const MDArrayRCP< T > & a1,
                const MDArrayView< T > & a2)
{
  return (a1() != a2);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
std::ostream & operator<<(std::ostream & os,
                          const MDArrayRCP< T > & a)
{
  os << a.toString();
  return os;
}

//////////////////////////
// Private implementations
//////////////////////////

template< typename T >
void
MDArrayRCP< T >::assertAxis(int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(0 <= axis && axis < _dimensions.size()),
    RangeError,
    "MDArrayRCP<T>::assertAxis(axis=" << axis << "): out of "
    << "range axis in [0, " << _dimensions.size() << ")"
    );
}

////////////////////////////////////////////////////////////////////////

template< typename T >
void
MDArrayRCP< T >::assertIndex(dim_type i, int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(0 <= i && i < _dimensions[axis]), RangeError,
    "MDArrayRCP<T>::assertIndex(i=" << i << ",axis=" << axis << "): out of "
    << "range i in [0, " << _dimensions[axis] << ")"
    );
}

}  // End namespace Domi

#endif  // DOMI_MDARRAYRCP_HPP
