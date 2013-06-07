// @HEADER
// ***********************************************************************
//
//            Domi: Multidimensional Datastructures Package
//                 Copyright (2013) Sandia Corporation
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef DOMI_MDARRAYRCP_HPP
#define DOMI_MDARRAYRCP_HPP

#include "Teuchos_ArrayRCPDecl.hpp"
#include "Domi_MDArray_Utils.hpp"
#include "Domi_MDArrayView.hpp"

#include <cstdarg>

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
bool operator==(const MDArrayRCP< T > & a1, const MDArrayRCP< T > & a2);

/** \brief Inequality operator.
 *
 * \relates MDArrayRCP
 */
template< typename T >
bool operator!=(const MDArrayRCP< T > & a1, const MDArrayRCP< T > & a2);

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

  /** \brief Ordinal */
  typedef Teuchos_Ordinal Ordinal;

  /** \brief Size type */
  typedef Ordinal size_type;

  /** \brief Value type */
  typedef T value_type;

  /** \brief Pointer type */
  typedef T* pointer_type;

  //@}

  /** \name Constructors, Destructor, Initializers */
  //@{

  /** \brief Default constructor
   *
   * \param null_arg [in] Optional null pointer specifier
   */
  inline MDArrayRCP(Teuchos::ENull null_arg = Teuchos::null);

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
   * \param storageOrder [in] Specifies the order data elements are
   *        stored in memory (default DEFAULT_ORDER)
   *
   * The <tt>MDArrayRCP</tt> does not take ownership of the data
   * buffer.
   */
  inline MDArrayRCP(const Teuchos::ArrayView< T > & array,
		    const Teuchos::ArrayView< size_type > & dims,
		    EStorageOrder storageOrder=DEFAULT_ORDER);

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
   * \param storageOrder [in] Specifies the order data elements are
   *        stored in memory
   *
   * This constructor allocates new memory and takes ownership of it.
   */
  inline explicit MDArrayRCP(const Teuchos::ArrayView< size_type > & dims,
			     const T & val=T(),
			     EStorageOrder storageOrder=DEFAULT_ORDER);

  /** \brief Constructor with dimensions and storage order flag.
   *
   * \param dims [in] An array that defines the lengths of each
   *        dimension.  The most convenient way to specify dimensions
   *        is with a Tuple returned by the non-member
   *        <tt>Teuchos::tuple<T>()</tt> function.
   *
   * \param storageOrder [in] Specifies the order data elements are
   *        stored in memory
   *
   * This constructor allocates new memory and takes ownership of it.
   */
  inline explicit MDArrayRCP(const Teuchos::ArrayView< size_type > & dims,
			     EStorageOrder storageOrder);

  /** \brief Shallow copy constructor
   *
   * \param r_ptr [in] Source reference counted pointer
   */
  inline MDArrayRCP(const MDArrayRCP< T > & r_ptr);

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
  inline int num_dims() const;

  /** \brief Return the array of dimensions
   */
  inline const Teuchos::Array< typename Teuchos::Array< T >::size_type > &
  dimensions() const;

  /** \brief Return the dimension of the given axis
   *
   * \param axis [in] The axis being queried (0 for the first axis,
   *        for the second axis, and so forth)
   */
  inline typename Teuchos::Array< T >::size_type dimension(int axis) const;

  /** \brief Return the total size of the <tt>MDArrayRCP</tt>
   */
  inline typename Teuchos::Array< T >::size_type size() const;

  /** \brief Return the indexing strides
   */
  inline const Teuchos::Array< typename Teuchos::Array< T >::size_type > &
  strides() const;

  /** \brief Return the underlying <tt>Teuchos::ArrayRCP</tt>
   */
  inline const Teuchos::ArrayRCP< T > & arrayRCP() const;

  /** \brief Return the storage order
   */
  inline const EStorageOrder storage_order() const;

  //@}

  /** \name Object/pointer access functions */
  //@{

  /** \brief Return true if the underlying pointer is null
   */
  inline bool is_null() const;

  /** \brief Pointer <tt>-></tt> access to members of underlying data buffer
   */
  inline T * operator->() const;

  /** \brief Dereference the underlying data buffer
   */
  inline T & operator*();

  /** \brief Get the raw C++ pointer to the underlying data buffer
   */
  inline T * get() const;

  //@}

  /** \name Conversions to MDArrayView */
  //@{

  /** \brief Perform an explicit conversion to a non-const
   *  <tt>MDArrayView</tt>
   */
  MDArrayView< T > mdArrayView();

  /** \brief Perform an explicit conversion to a const
   *  <tt>MDArrayView</tt>
   */
  MDArrayView< const T > mdArrayViewConst();

  /** \brief Perform an implicit conversion to a non-const
   *  <tt>MDArrayView</tt>
   */
  inline operator MDArrayView< T >();

  /** \brief Perform an implicit conversion to a const
   *  <tt>MDArrayView</tt>
   */
  inline operator MDArrayView< const T >();

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
  MDArrayView< T > operator[](size_type i);

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

  /** \brief Conversion to <tt>MDArrayView</tt>
   */
  MDArrayView< T > operator()();

  //@}

  /** \name Indexing operators that return a reference to a single array element */
  //@{

  /** \brief Non-const 1D element access operator
   *
   * \param i [in] 1D index.
   *
   * This operator should only be used with a 1D <tt>MDArrayRCP</tt>.
   * If Domi_ENABLE_ABC is true and the <tt>MDArrayRCP</tt> is not
   * 1D, an exception will be thrown.
   */
  inline T & operator()(size_type i);

  /** \brief Non-const 2D element access operator
   *
   * \param i [in] first 2D index.
   *
   * \param j [in] second 2D index.
   *
   * This operator should only be used with a 2D <tt>MDArrayRCP</tt>.
   * If Domi_ENABLE_ABC is true and the <tt>MDArrayRCP</tt> is not
   * 2D, an exception will be thrown.
   */
  inline T & operator()(size_type i, size_type j);

  /** \brief Non-const 3D element access operator
   *
   * \param i [in] first 3D index.
   *
   * \param j [in] second 3D index.
   *
   * \param k [in] third 3D index.
   *
   * This operator should only be used with a 3D <tt>MDArrayRCP</tt>.
   * If Domi_ENABLE_ABC is true and the <tt>MDArrayRCP</tt> is not
   * 3D, an exception will be thrown.
   */
  inline T & operator()(size_type i, size_type j, size_type k);

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
   * If Domi_ENABLE_ABC is true and the <tt>MDArrayRCP</tt> is not
   * 4D, an exception will be thrown.
   */
  inline T & operator()(size_type i, size_type j, size_type k, size_type m);

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
   * If Domi_ENABLE_ABC is true and the <tt>MDArrayRCP</tt> is not
   * 5D, an exception will be thrown.
   */
  inline T & operator()(size_type i, size_type j, size_type k, size_type m,
                        size_type n);

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
   * <tt>MDArrayRCP</tt>s.  If Domi_ENABLE_ABC is true and the
   * <tt>MDArrayRCP</tt> is less than 6D, an exception will be thrown.
   */
  inline T & operator()(size_type i, size_type j, size_type k, size_type m,
                        size_type n, size_type p, ...);

  /** \brief Const 1D element access operator
   *
   * \param i [in] 1D index.
   *
   * This operator should only be used with a 1D <tt>MDArrayRCP</tt>.
   * If Domi_ENABLE_ABC is true and the <tt>MDArrayRCP</tt> is not
   * 1D, an exception will be thrown.
   */
  inline const T & operator()(size_type i) const;

  /** \brief Const 2D element access operator
   *
   * \param i [in] first 2D index.
   *
   * \param j [in] second 2D index.
   *
   * This operator should only be used with a 2D <tt>MDArrayRCP</tt>.
   * If Domi_ENABLE_ABC is true and the <tt>MDArrayRCP</tt> is not
   * 2D, an exception will be thrown.
   */
  inline const T & operator()(size_type i, size_type j) const;

  /** \brief Const 3D element access operator
   *
   * \param i [in] first 3D index.
   *
   * \param j [in] second 3D index.
   *
   * \param k [in] third 3D index.
   *
   * This operator should only be used with a 3D <tt>MDArrayRCP</tt>.
   * If Domi_ENABLE_ABC is true and the <tt>MDArrayRCP</tt> is not
   * 3D, an exception will be thrown.
   */
  inline const T & operator()(size_type i, size_type j, size_type k) const;

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
   * If Domi_ENABLE_ABC is true and the <tt>MDArrayRCP</tt> is not
   * 4D, an exception will be thrown.
   */
  inline const T & operator()(size_type i, size_type j, size_type k,
                              size_type m) const;

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
   * If Domi_ENABLE_ABC is true and the <tt>MDArrayRCP</tt> is not
   * 5D, an exception will be thrown.
   */
  inline const T & operator()(size_type i, size_type j, size_type k,
                              size_type m, size_type n) const;

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
   * <tt>MDArrayRCP</tt>s.  If Domi_ENABLE_ABC is true and the
   * <tt>MDArrayRCP</tt> is less than 6D, an exception will be thrown.
   */
  inline const T & operator()(size_type i, size_type j, size_type k,
                              size_type m, size_type n, size_type p, ...) const;


  //@}

  /** \name Teuchos::Array-like and std::vector-like methods */
  //@{

  /** \brief Assign a value to all elements of the <tt>MDArrayRCP</tt>
   *
   * \param value [in] The value to be assigned
   */
  void assign(const T & value);

  /** \brief Non-const single element access method with bounds
   * checking
   *
   * \param i, ... [in] Indexes representing the location of the
   *        single element of the <tt>MDArrayRCP</tt> to be accessed.
   *        Note that this method assumes that the user will provide
   *        the same number of arguments as the number of dimensions
   *        of the <tt>MDArrayRCP</tt>.
   */
  T & at(size_type i, ...);

  /** \brief Const single element access method with bounds checking
   *
   * \param i, ... [in] Indexes representing the location of the
   *        single element of the <tt>MDArrayRCP</tt> to be accessed.
   *        Note that this method assumes that the user will provide
   *        the same number of arguments as the number of dimensions
   *        of the <tt>MDArrayRCP</tt>.
   */
  const T & at(size_type i, ...) const;

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
  void resize(const Teuchos::ArrayView< size_type > & dims);

  /** \brief Return true if <tt>MDArrayRCP</tt> has been compiled with
   *  bounds checking on. */
  inline static bool hasBoundsChecking();

  /** \brief Convert the <tt>MDArrayRCP</tt> to a string
   *  representation
   */
  std::string toString();

  /** \brief Return a raw pointer to the beginning of the
   *  <tt>MDArrayRCP</tt> or NULL if unsized. */
  inline T * getRawPtr();

  /** \brief Return a const raw pointer to the beginning of the
   *  <tt>MDArrayRCP</tt> or NULL if unsized. */
  inline const T * getRawPtr() const;

  //@}

private:
  Teuchos::Array< size_type > _dimensions;
  Teuchos::Array< size_type > _strides;
  Teuchos::ArrayRCP< T >      _array;
  EStorageOrder               _storage_order;
  T *                         _ptr;

  // Used for array bounds checking
  void assertIndex(size_type i, int axis) const;
};

//////////////////
// Implementations
//////////////////

template< typename T >
MDArrayRCP< T >::MDArrayRCP(Teuchos::ENull null_arg) :
  _dimensions(Teuchos::tuple< size_type >(0)),
  _strides(Teuchos::tuple< size_type >(1)), 
  _array(),
  _storage_order(DEFAULT_ORDER),
  _ptr()
{
}

template< typename T >
MDArrayRCP< T >::MDArrayRCP(const Teuchos::ArrayView< T > & array,
			    const Teuchos::ArrayView< size_type > & dims,
			    EStorageOrder storageOrder) :
  _dimensions(dims),
  _strides(computeStrides(dims, storageOrder)),
  _array(array.getRawPtr(), 0, array.size(), false),
  _storage_order(storageOrder),
  _ptr(_array.getRawPtr())
{
  TEUCHOS_TEST_FOR_EXCEPTION(array.size() < computeSize(dims),
			     Teuchos::RangeError,
			     "Teuchos::ArrayView size too small for "
                             "dimensions");
}

template< typename T >
MDArrayRCP< T >::MDArrayRCP(const Teuchos::ArrayView< size_type > & dims,
			    const T & val,
			    EStorageOrder storageOrder) :
  _dimensions(dims),
  _strides(computeStrides(dims, storageOrder)),
  _array(computeSize(dims), val),
  _storage_order(storageOrder),
  _ptr(_array.getRawPtr())
{
}

template< typename T >
MDArrayRCP< T >::MDArrayRCP(const Teuchos::ArrayView< size_type > & dims,
			    EStorageOrder storageOrder) :
  _dimensions(dims),
  _strides(computeStrides(dims, storageOrder)),
  _array(computeSize(dims)),
  _storage_order(storageOrder),
  _ptr(_array.getRawPtr())
{
}

template< typename T >
MDArrayRCP< T >::MDArrayRCP(const MDArrayRCP< T > & r_ptr) :
  _dimensions(r_ptr._dimensions),
  _strides(r_ptr._strides),
  _array(r_ptr._array),
  _storage_order(r_ptr._storage_order),
  _ptr(_array.getRawPtr())
{
}

template< typename T >
MDArrayRCP< T >::~MDArrayRCP()
{
}

template< typename T >
MDArrayRCP< T > &
MDArrayRCP< T >::operator=(const MDArrayRCP< T > & r_ptr)
{
  _dimensions    = r_ptr._dimensions;
  _strides       = r_ptr._strides;
  _array         = r_ptr._array;
  _storage_order = r_ptr._storage_order;
  _ptr           = r_ptr._ptr;
}

template< typename T >
int
MDArrayRCP< T >::num_dims() const
{
  return _dimensions.size();
}

template< typename T >
const Teuchos::Array< typename Teuchos::Array< T >::size_type > &
MDArrayRCP< T >::dimensions() const
{
  return _dimensions;
}

template< typename T >
typename Teuchos::Array< T >::size_type
MDArrayRCP< T >::dimension(int axis) const
{
  return _dimensions[axis];
}

template< typename T >
typename Teuchos::Array< T >::size_type
MDArrayRCP< T >::size() const
{
  return _array.size();
}

template< typename T >
const Teuchos::Array< typename Teuchos::Array< T >::size_type > &
MDArrayRCP< T >::strides() const
{
  return _strides;
}

template< typename T >
const Teuchos::ArrayRCP< T > &
MDArrayRCP< T >::arrayRCP() const
{
  return _array;
}

template< typename T >
const EStorageOrder
MDArrayRCP< T >::storage_order() const
{
  return _storage_order;
}

template< typename T >
bool
MDArrayRCP< T >::is_null() const
{
  return _array.is_null();
}

template< typename T >
T *
MDArrayRCP< T >::operator->() const
{
  return _array.operator->();
}

template< typename T >
T &
MDArrayRCP< T >::operator*()
{
  return _array.operator*();
}

template< typename T >
T *
MDArrayRCP< T >::get() const
{
  return _array.get();
}

template< typename T >
MDArrayView< T >
MDArrayRCP< T >::mdArrayView()
{
  return MDArrayView< T >(_array(), _dimensions, _storage_order);
}

template< typename T >
MDArrayView< const T >
MDArrayRCP< T >::mdArrayViewConst()
{
  return MDArrayView< const T >(_array(), _dimensions, _storage_order);
}

template< typename T >
MDArrayRCP< T >::operator MDArrayView< T >()
{
  return mdArrayView();
}

template< typename T >
MDArrayRCP< T >::operator MDArrayView< const T >()
{
  return mdArrayViewConst();
}

template< typename T >
MDArrayView< T >
MDArrayRCP< T >::operator[](size_type i)
{
  // Note: array bounds checking, if active, will be performed by the
  // MDArrayView class
  return mdArrayView()[i];
}

template< typename T >
MDArrayView< T >
MDArrayRCP< T >::operator[](Slice s)
{
  // Note: Slices produce safe indexes
  return mdArrayView()[s];
}

template< typename T >
MDArrayView< T >
MDArrayRCP< T >::operator()()
{
  return mdArrayView();
}

template< typename T >
T &
MDArrayRCP< T >::operator()(typename MDArrayRCP< T >::size_type i)
{
#ifdef Domi_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() != 1), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with 1 index"
    );
  assertIndex(i, 0);
#endif
  return _ptr[i * _strides[0]];
}

template< typename T >
T &
MDArrayRCP< T >::operator()(typename MDArrayRCP< T >::size_type i,
                            typename MDArrayRCP< T >::size_type j)
{
#ifdef Domi_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() != 2), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with 2 indexes"
    );
  assertIndex(i, 0);
  assertIndex(j, 1);
#endif
  return _ptr[i * _strides[0] + j * _strides[1]];
}

template< typename T >
T &
MDArrayRCP< T >::operator()(typename MDArrayRCP< T >::size_type i,
                            typename MDArrayRCP< T >::size_type j,
                            typename MDArrayRCP< T >::size_type k)
{
#ifdef Domi_ENABLE_ABC
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

template< typename T >
T &
MDArrayRCP< T >::operator()(typename MDArrayRCP< T >::size_type i,
                            typename MDArrayRCP< T >::size_type j,
                            typename MDArrayRCP< T >::size_type k,
                            typename MDArrayRCP< T >::size_type m)
{
#ifdef Domi_ENABLE_ABC
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

template< typename T >
T &
MDArrayRCP< T >::operator()(typename MDArrayRCP< T >::size_type i,
                            typename MDArrayRCP< T >::size_type j,
                            typename MDArrayRCP< T >::size_type k,
                            typename MDArrayRCP< T >::size_type m,
                            typename MDArrayRCP< T >::size_type n)
{
#ifdef Domi_ENABLE_ABC
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

template< typename T >
T &
MDArrayRCP< T >::operator()(typename MDArrayRCP< T >::size_type i,
                            typename MDArrayRCP< T >::size_type j,
                            typename MDArrayRCP< T >::size_type k,
                            typename MDArrayRCP< T >::size_type m,
                            typename MDArrayRCP< T >::size_type n,
                            typename MDArrayRCP< T >::size_type p,
                            ...)
{
#ifdef Domi_ENABLE_ABC
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
  for (size_type axis = 6; axis < _dimensions.size(); axis++)
  {
    size_type q = va_arg(indexes, size_type);
#ifdef Domi_ENABLE_ABC
    assertIndex(q, axis);
#endif
    offset += q * _strides[axis];
  }
  va_end(indexes);
  return _ptr[offset];
}

template< typename T >
const T &
MDArrayRCP< T >::operator()(typename MDArrayRCP< T >::size_type i) const
{
#ifdef Domi_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() != 1), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with 1 index"
    );
  assertIndex(i, 0);
#endif
  return _ptr[i * _strides[0]];
}

template< typename T >
const T &
MDArrayRCP< T >::operator()(typename MDArrayRCP< T >::size_type i,
                            typename MDArrayRCP< T >::size_type j) const
{
#ifdef Domi_ENABLE_ABC
  TEUCHOS_TEST_FOR_EXCEPTION(
    (_dimensions.size() != 2), RangeError,
    "Attempt to access " << _dimensions.size() << "D array with 2 indexes"
    );
  assertIndex(i, 0);
  assertIndex(j, 1);
#endif
  return _ptr[i * _strides[0] + j * _strides[1]];
}

template< typename T >
const T &
MDArrayRCP< T >::operator()(typename MDArrayRCP< T >::size_type i,
                            typename MDArrayRCP< T >::size_type j,
                            typename MDArrayRCP< T >::size_type k) const
{
#ifdef Domi_ENABLE_ABC
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

template< typename T >
const T &
MDArrayRCP< T >::operator()(typename MDArrayRCP< T >::size_type i,
                            typename MDArrayRCP< T >::size_type j,
                            typename MDArrayRCP< T >::size_type k,
                            typename MDArrayRCP< T >::size_type m) const
{
#ifdef Domi_ENABLE_ABC
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

template< typename T >
const T &
MDArrayRCP< T >::operator()(typename MDArrayRCP< T >::size_type i,
                            typename MDArrayRCP< T >::size_type j,
                            typename MDArrayRCP< T >::size_type k,
                            typename MDArrayRCP< T >::size_type m,
                            typename MDArrayRCP< T >::size_type n) const
{
#ifdef Domi_ENABLE_ABC
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

template< typename T >
const T &
MDArrayRCP< T >::operator()(typename MDArrayRCP< T >::size_type i,
                            typename MDArrayRCP< T >::size_type j,
                            typename MDArrayRCP< T >::size_type k,
                            typename MDArrayRCP< T >::size_type m,
                            typename MDArrayRCP< T >::size_type n,
                            typename MDArrayRCP< T >::size_type p,
                            ...) const
{
#ifdef Domi_ENABLE_ABC
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
  for (size_type axis = 6; axis < _dimensions.size(); axis++)
  {
    size_type q = va_arg(indexes, size_type);
#ifdef Domi_ENABLE_ABC
    assertIndex(q, axis);
#endif
    offset += q * _strides[axis];
  }
  va_end(indexes);
  return _ptr[offset];
}

template< typename T >
void
MDArrayRCP< T >::assign(const T & value)
{
  std::cout << "MDArray<T>::assign(const T &) will not be implemented until "
            << "iterators are implemented." << std::endl;
}

template< typename T >
T &
MDArrayRCP< T >::at(size_type i, ...)
{
  assertIndex(i, 0);
  va_list indexes;
  size_type offset = i * _strides[0];
  va_start(indexes, i);
  for (size_type axis = 1; axis < _dimensions.size(); axis++)
  {
    size_type j = va_arg(indexes, size_type);
    assertIndex(j, axis);
    offset += j * _strides[axis];
  }
  va_end(indexes);
  return _array[offset];
}

template< typename T >
const T &
MDArrayRCP< T >::at(size_type i, ...) const
{
  assertIndex(i, 0);
  va_list indexes;
  size_type offset = i * _strides[0];
  va_start(indexes, i);
  for (size_type axis = 1; axis < _dimensions.size(); axis++)
  {
    size_type j = va_arg(indexes, size_type);
    assertIndex(j, axis);
    offset += j * _strides[axis];
  }
  va_end(indexes);
  return _array[offset];
}

template< typename T >
typename MDArrayRCP< T >::size_type
MDArrayRCP< T >::capacity() const
{
  return _array.capacity();
}

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

template< typename T >
bool
MDArrayRCP< T >::empty() const
{
  return (_array.size() == 0);
}

template< typename T >
typename MDArrayRCP< T >::size_type
MDArrayRCP< T >::max_size() const
{
  return _array.max_size();
}

template< typename T >
void
MDArrayRCP< T >::resize(const Teuchos::ArrayView< size_type > & dims)
{
  _dimensions.assign(dims.begin(), dims.end());
  _strides = computeStrides(dims, _storage_order);
  _array.resize(computeSize(dims));
  _ptr = _array.getRawPtr();
}

template< typename T >
bool
MDArrayRCP< T >::hasBoundsChecking()
{
#ifdef Domi_ENABLE_ABC
  return true;
#else
  return false;
#endif
}

template< typename T >
std::string
MDArrayRCP< T >::toString()
{
  return mdArrayView().toString();
}

template< typename T >
T *
MDArrayRCP< T >::getRawPtr()
{
  return _array.getRawPtr();
}

template< typename T >
const T *
MDArrayRCP< T >::getRawPtr() const
{
  return _array.getRawPtr();
}

template< typename T >
bool operator==(const MDArrayRCP< T > & a1, const MDArrayRCP< T > & a2)
{
  if (a1._dimensions != a2._dimensions) return false;
  if (a1._storage_order == a2._storage_order)
    return (a1._array == a2._array);
  // This should be changed to compare individual elements when more
  // sophisticated iterators are available
  return false;
}

template< typename T >
bool operator!=(const MDArrayRCP< T > & a1, const MDArrayRCP< T > & a2)
{
  return not (a1 == a2);
}

template< typename T >
std::ostream & operator<<(std::ostream & os, const MDArrayRCP< T > & a)
{
  os << a.toString();
  return os;
}

//////////////////////////
// Private implementations
//////////////////////////

template< typename T >
void
MDArrayRCP< T >::assertIndex(size_type i, int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(0 <= i && i < _dimensions[axis]), Teuchos::RangeError,
    "MDArrayRCP<T>::assertIndex(i=" << i << ",axis=" << axis << "): out of "
    << "range i in [0, " << _dimensions[axis] << ")"
    );
}

}  // End namespace Domi

#endif  // DOMI_MDARRAYRCP_HPP
