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

#ifndef DOMI_MDARRAYVIEW_HPP
#define DOMI_MDARRAYVIEW_HPP

// Standar includes
#include <cstdarg>

// Teuchos includes
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_ConstTypeTraits.hpp"
#include "Teuchos_RCPNode.hpp"

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_Exceptions.hpp"
#include "Domi_Utils.hpp"
#include "Domi_Slice.hpp"
#include "Domi_MDIterator.hpp"
#include "Domi_MDRevIterator.hpp"

namespace Domi
{

// I put these non-member template functions here for the same reason
// that Ross did the same thing for the Teuchos::Array class.  See
// Teuchos_Array.hpp for details.
template< typename T > class MDArray;
template< typename T > class MDArrayView;
template< typename T > class MDArrayRCP;

/** \brief Equality operator.
 *
 * \relates MDArrayView
 */
template< typename T >
bool operator==(const MDArrayView< T > & a1,
                const MDArrayView< T > & a2);

/** \brief Inequality operator.
 *
 * \relates MDArrayView
 */
template< typename T >
bool operator!=(const MDArrayView< T > & a1,
                const MDArrayView< T > & a2);

/** \brief Non-member swap
 *
 * \relates MDArray
 */
template< typename T >
void swap(MDArrayView< T > & a1,
          MDArrayView< T > & a2);

/** \brief Memory-safe templated multi-dimensional array view class
 *
 * The <tt>MDArrayView</tt> class is to the <tt>MDArray</tt> class as
 * the <tt>Teuchos::ArrayView</tt> class is to the
 * <tt>Teuchos::Array</tt> class.  Its primary function is to provide
 * sub-views into both <tt>MDArray</tt> objects and other
 * <tt>MDArrayView</tt> objects.
 *
 * <tt>MDArrayView</tt> objects have a similar interface to the
 * <tt>MDArray</tt> class and much of the behavior is the same as is
 * described in the <tt>MDArray</tt> documentation.  The big difference
 * is that <tt>MDArrayView</tt> objects are constructed with
 * pre-existing data buffers.
 *
 * \ingroup domi_mem_mng_grp
 */
template< typename T >
class MDArrayView
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

  /** @name Constructors and Destructor */
  //@{

  /** \brief Default constructor
   *
   * \param null_arg [in] Enumerated constant denoting null
   *        construction
   *
   * Returns a view into an array of 1 dimension of length zero with
   * NULL pointer.
   */
  MDArrayView(Teuchos::ENull null_arg = Teuchos::null);

  /** \brief Constructor with a source <tt>Teuchos::ArrayView</tt>,
   *   dimensions, and optional storage order
   *
   * \param array [in] <tt>Teuchos::ArrayView</tt> of the data buffer
   *
   * \param dims [in] An array that defines the lengths of each
   *        dimension.  The most convenient way to specify dimensions
   *        is with a Tuple returned by the non-member
   *        <tt>Teuchos::tuple<T>()</tt> function.
   *
   * \param layout [in] An enumerated value specifying the internal
   *        storage order of the <tt>MDArrayView</tt>
   */
  MDArrayView(const Teuchos::ArrayView< T > & array,
	      const Teuchos::ArrayView< dim_type > & dims,
	      const Layout layout = DEFAULT_ORDER);

  /** \brief Constructor with a source <tt>Teuchos::ArrayView</tt>,
   *   dimensions, strides, and optional storage order
   *
   * \param array [in] <tt>Teuchos::ArrayView</tt> of the data buffer
   *
   * \param dims [in] An array that defines the lengths of each
   *        dimension.
   *
   * \param strides [in] An array that defines the strides between
   *        elements along each axis.
   *
   * \param layout [in] An enumerated value specifying the internal
   *        storage order of the <tt>MDArrayView</tt>
   *
   * Note that this constructor was introduced specifically to
   * implement the getConst() method, which converts an MDArrayView< T >
   * to an MDArrayView< const T >.  The presence of both strides and
   * layout in the input arguments makes this an easy constructor to
   * call incorrectly, and is not advised for general use.
   */
  MDArrayView(const Teuchos::ArrayView< T > & array,
	      const Teuchos::Array< dim_type > & dims,
              const Teuchos::Array< size_type > & strides,
	      const Layout layout = DEFAULT_ORDER);

  /** \brief Copy constructor
   *
   * \param array [in] The source <tt>MDArrayView</tt> to be copied
   */
  MDArrayView(const MDArrayView< T > & array);

  /* \brief Parent/single index sub-array view constructor
   *
   * \param parent [in] an MDArrayView, from which this MDArrayView
   *        will be derived
   *
   * \param axis [in] the axis to which this index ordinal applies
   *
   * \param index [in] the ordinal that defines this sub-array
   */
  MDArrayView(const MDArrayView< T > & parent,
              int axis,
              dim_type index);

  /* \brief Parent/single slice sub-array view constructor
   *
   * \param parent [in] an MDArrayView, from which this MDArrayView
   *        will be derived
   *
   * \param axis [in] the axis to which this index ordinal applies
   *
   * \param slice [in] the slice that defines this sub-array
   */
  MDArrayView(const MDArrayView< T > & parent,
              int axis,
              Slice slice);

  /** \brief Assignment operator
   *
   * \param array [in] The source <tt>MDArrayView</tt> to be copied
   */
  MDArrayView< T > & operator=(const MDArrayView< T > & array);

  /** \brief Destructor
   */
  ~MDArrayView();

  //@}

  /** \name Attribute accesor methods */
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

  /** \brief Return the total size of the <tt>MDArrayView</tt>
   */
  inline size_type size();

  /** \brief Return the indexing strides
   */
  inline const Teuchos::Array< size_type > & strides() const;

  /** \brief Return the underlying <tt>Teuchos::ArrayView</tt>
   */
  inline const Teuchos::ArrayView< T > & arrayView() const;

  /** \brief Return the storage order
   */
  inline Layout layout() const;

  /** \brief Return whether the MDArrayView is contiguous in memory
   */
  inline bool contiguous() const;

  //@}

  /** \name Iterator classes and methods */
  //@{

  friend class MDIterator< MDArrayView< T > >;
  friend class MDIterator< MDArrayView< const T > >;
  friend class MDRevIterator< MDArrayView< T > >;
  friend class MDRevIterator< MDArrayView< const T > >;

  typedef MDIterator< MDArrayView< T > >          iterator;
  typedef MDIterator< MDArrayView< const T > >    const_iterator;
  typedef MDRevIterator< MDArrayView< T > >       reverse_iterator;
  typedef MDRevIterator< MDArrayView< const T > > const_reverse_iterator;

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

  /** \name Indexing operators that return MDArrayViews */
  //@{

  /** \brief Return an MDArrayView< const T > of an MDArrayView< T > object
   * 
   * WARNING!  If <tt>T</tt> is already const (e.g. <tt>const double</tt>)
   * then do not try to instantiate this function since it will not compile!
   */
  MDArrayView< const T > getConst() const;

  /** \brief Sub-array access operator.  The returned
   *  <tt>MDArrayView</tt> object will have one fewer dimensions than
   *  the calling <tt>MDArrayView</tt>.
   *
   * \param i [in] Index of the desired sub-array.  Note that to
   *        obtain expected behavior, you should always chain together
   *        <tt>n</tt> square bracket operators when referencing an
   *        <tt>n</tt>-dimensional <tt>MDArrayView</tt>.
   */
  MDArrayView< T > operator[](dim_type i);

  /** \brief Sub-array const access operator.  The returned
   *  <tt>MDArrayView</tt> object will have one fewer dimensions than
   *  the calling <tt>MDArrayView</tt>.
   *
   * \param i [in] Index of the desired sub-array.  Note that to
   *        obtain expected behavior, you should always chain together
   *        <tt>n</tt> square bracket operators when referencing an
   *        <tt>n</tt>-dimensional <tt>MDArrayView</tt>.
   */
  const MDArrayView< T > operator[](dim_type i) const;

  /** \brief Sub-array access operator.  The returned
   *  <tt>MDArrayView</tt> object will have the same number of
   *  dimensions as the calling <tt>MDArrayView</tt>.
   *
   * \param s [in] Slice representing the bounds of the desired
   *        sub-array.  Note that to obtain expected behavior, you
   *        should always chain together <tt>n</tt> square bracket
   *        operators when referencing an <tt>n</tt>-dimensional
   *        <tt>MDArrayView</tt>.
   */
  MDArrayView< T > operator[](Slice s);

  /** \brief Sub-array const access operator.  The returned
   *  <tt>MDArrayView</tt> object will have the same number of
   *  dimensions as the calling <tt>MDArrayView</tt>.
   *
   * \param s [in] Slice representing the bounds of the desired
   *        sub-array.  Note that to obtain expected behavior, you
   *        should always chain together <tt>n</tt> square bracket
   *        operators when referencing an <tt>n</tt>-dimensional
   *        <tt>MDArrayView</tt>.
   */
  const MDArrayView< T > operator[](Slice s) const;

  //@}

  /** \name Indexing operators that return a reference to a single array element */
  //@{

  /** \brief Non-const 1D element access operator
   *
   * \param i [in] 1D index.
   *
   * This operator should only be used with a 1D <tt>MDArrayView</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the <tt>MDArrayView</tt> is not
   * 1D, an exception will be thrown.
   */
  inline T & operator()(dim_type i);

  /** \brief Non-const 2D element access operator
   *
   * \param i [in] first 2D index.
   *
   * \param j [in] second 2D index.
   *
   * This operator should only be used with a 2D <tt>MDArrayView</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the <tt>MDArrayView</tt> is not
   * 2D, an exception will be thrown.
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
   * This operator should only be used with a 3D <tt>MDArrayView</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the <tt>MDArrayView</tt> is not
   * 3D, an exception will be thrown.
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
   * This operator should only be used with a 4D <tt>MDArrayView</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the <tt>MDArrayView</tt> is not
   * 4D, an exception will be thrown.
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
   * This operator should only be used with a 5D <tt>MDArrayView</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the <tt>MDArrayView</tt> is not
   * 5D, an exception will be thrown.
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
   * <tt>MDArrayView</tt>s.  If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the
   * <tt>MDArrayView</tt> is less than 6D, an exception will be
   * thrown.
   */
  inline T & operator()(dim_type i, dim_type j, dim_type k, dim_type m,
                        dim_type n, dim_type p, ...);

  /** \brief Const 1D element access operator
   *
   * \param i [in] 1D index.
   *
   * This operator should only be used with a 1D <tt>MDArrayView</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the <tt>MDArrayView</tt> is not
   * 1D, an exception will be thrown.
   */
  inline const T & operator()(dim_type i) const;

  /** \brief Const 2D element access operator
   *
   * \param i [in] first 2D index.
   *
   * \param j [in] second 2D index.
   *
   * This operator should only be used with a 2D <tt>MDArrayView</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the <tt>MDArrayView</tt> is not
   * 2D, an exception will be thrown.
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
   * This operator should only be used with a 3D <tt>MDArrayView</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the <tt>MDArrayView</tt> is not
   * 3D, an exception will be thrown.
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
   * This operator should only be used with a 4D <tt>MDArrayView</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the <tt>MDArrayView</tt> is not
   * 4D, an exception will be thrown.
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
   * This operator should only be used with a 5D <tt>MDArrayView</tt>.
   * If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the <tt>MDArrayView</tt> is not
   * 5D, an exception will be thrown.
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
   * <tt>MDArrayView</tt>s.  If HAVE_DOMI_ARRAY_BOUNDSCHECK is true and the
   * <tt>MDArrayView</tt> is less than 6D, an exception will be
   * thrown.
   */
  inline const T & operator()(dim_type i, dim_type j, dim_type k,
                              dim_type m, dim_type n, dim_type p, ...) const;


  //@}

  /** \name Teuchos::Array-like and std::vector-like methods */
  //@{

  /** \brief Assign a value to all elements of the
   * <tt>MDArrayView</tt>
   *
   * \param value [in] The value to be assigned
   */
  void assign(const T & value);

  /** \brief Non-const single element access method with bounds
   * checking
   *
   * \param i, ... [in] Indexes representing the location of the
   *        single element of the <tt>MDArrayView</tt> to be accessed.
   *        Note that this method assumes that the user will provide
   *        the same number of arguments as the number of dimensions
   *        of the <tt>MDArrayView</tt>.
   */
  T & at(dim_type i, ...);

  /** \brief Const single element access method with bounds checking
   *
   * \param i, ... [in] Indexes representing the location of the
   *        single element of the <tt>MDArrayView</tt> to be accessed.
   *        Note that this method assumes that the user will provide
   *        the same number of arguments as the number of dimensions
   *        of the <tt>MDArrayView</tt>.
   */
  const T & at(dim_type i, ...) const;

  /** \brief Return true if <tt>MDArrayView</tt> has been compiled
   *  with bounds checking on.
   */
  inline static bool hasBoundsChecking();

  /** \brief Convert the <tt>MDArrayView</tt> to a string
   *  representation
   */
  std::string toString() const;

  /** \brief Return a const raw pointer to the beginning of the
   *  <tt>MDArrayView</tt> or NULL if unsized.
   */
  inline const T * getRawPtr() const;

  /** \brief Return a raw pointer to the beginning of the
   *  <tt>MDArrayView</tt> or NULL if unsized.
   */
  inline T * getRawPtr();

  // These operators are declared as friends so that the compiler will
  // do automatic type conversion.

  /** \name Non-member operators and functions */
  //@{

  /** \brief Equality operator.
   */
  template< typename T2 >
  friend bool operator==(const MDArrayView< T2 > & a1, const MDArrayView< T2 > & a2);

  /** \brief Inequality operator.
   */
  template< typename T2 >
  friend bool operator!=(const MDArrayView< T2 > & a1, const MDArrayView< T2 > & a2);

  /** \brief Swap function
   */
  template< typename T2 >
  friend void swap(MDArrayView< T2 > & a1, MDArrayView< T2 > & a2);

  /** \brief Stream output operator
   *
   * This operator calls the <tt>MDArrayView toString()</tt> method.
   */
  template< typename T2 >
  friend std::ostream & operator<<(std::ostream & os,
				   const MDArrayView< T2 > & a);

  //@}

private:

  Teuchos::Array< dim_type >  _dimensions;
  Teuchos::Array< size_type > _strides;
  Teuchos::ArrayView< T >     _array;
  Layout                      _layout;
  pointer                     _ptr;
  int                         _next_axis;

  // Method provided for aiding in array bounds checking.  It raises
  // an exception when the given axis is out of range of valid
  // dimensions.
  void assertAxis(int axis) const;

  // Method provided for aiding in array bounds checking.  It raises
  // an exception when the given index i is out of range along the
  // given axis.
  void assertIndex(dim_type i, int axis) const;

  // Private implementation of the toString() method that includes an
  // indent argument.  The publicly available version takes no
  // arguments and calls this one with initial indent value of zero.
  // This method can call itself recursively.
  std::string toString(int indent) const;

};  // class MDArrayView

/////////////////////
// Implementations //
/////////////////////

template< typename T >
MDArrayView< T >::
MDArrayView(Teuchos::ENull null_arg) :
  _dimensions(Teuchos::tuple< dim_type >(0)),
  _strides(Teuchos::tuple< size_type >(1)),
  _array(),
  _layout(DEFAULT_ORDER),
  _ptr(),
  _next_axis(0)
{
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayView< T >::MDArrayView(const Teuchos::ArrayView< T > & array,
			      const Teuchos::ArrayView< dim_type > & dims,
			      const Layout layout) :
  _dimensions(dims),
  _strides(computeStrides< size_type, dim_type >(dims, layout)),
  _array(array),
  _layout(layout),
  _ptr(_array.getRawPtr()),
  _next_axis(0)
{
  TEUCHOS_TEST_FOR_EXCEPTION(array.size() < computeSize(dims),
			     RangeError,
			     "Teuchos::ArrayView size too small for "
                             "dimensions");
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayView< T >::MDArrayView(const Teuchos::ArrayView< T > & array,
			      const Teuchos::Array< dim_type > & dims,
			      const Teuchos::Array< size_type > & strides,
			      const Layout layout) :
  _dimensions(dims),
  _strides(strides),
  _array(array),
  _layout(layout),
  _ptr(_array.getRawPtr()),
  _next_axis(0)
{
  const size_type required = computeSize< const size_type,
                                          const dim_type >(dims(),strides());
  TEUCHOS_TEST_FOR_EXCEPTION(
    array.size() < required,
    RangeError,
    "Teuchos::ArrayView size too small for "
    "dimensions and strides");
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayView< T >::MDArrayView(const MDArrayView< T > & array) :
  _dimensions(array._dimensions),
  _strides(array._strides),
  _array(array._array),
  _layout(array._layout),
  _ptr(_array.getRawPtr()),
  _next_axis(0)
{
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayView< T >::MDArrayView(const MDArrayView< T > & parent,
                              int axis,
                              dim_type index) :
  _dimensions(),
  _strides(),
  _array(),
  _layout(parent._layout),
  _ptr(),
  _next_axis(0)
{
  // Make sure axis and index are valid
  parent.assertAxis(axis);
  parent.assertIndex(index, axis);
  // Find the offset to the new MDArrayView
  size_type offset = index * parent._strides[axis];
  // Compute the dimensions of the new MDArrayView
  size_type n = parent._dimensions.size();
  // Compute the new dimensions and strides
  if (n == 1)
  {
    _dimensions.push_back(1);
    _strides.push_back(1);
  }
  else
  {
    for (int myAxis = 0; myAxis < n; myAxis++)
      if (myAxis != axis)
      {
	_dimensions.push_back(parent._dimensions[myAxis]);
	_strides.push_back(parent._strides[myAxis]);
      }
  }
  // Compute the new MDArrayView and pointer
  _array = parent._array.view(offset,
                              computeSize(_dimensions(),
                                          _strides()));
  _ptr = _array.getRawPtr();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayView< T >::MDArrayView(const MDArrayView< T > & parent,
                              int axis,
                              Slice slice) :
  _dimensions(parent._dimensions),
  _strides(parent._strides),
  _array(),
  _layout(parent._layout),
  _ptr(),
  _next_axis(0)
{
  // Make sure axis is valid
  parent.assertAxis(axis);
  // Note: the Slice.bounds() method produces safe indexes
  Slice bounds = slice.bounds(_dimensions[axis]);
  // Find the offset to the new MDArrayView
  size_type offset = bounds.start() * _strides[axis];
  // Compute the dimensions of the new MDArrayView
  _dimensions[axis] = (bounds.stop() - bounds.start()) / bounds.step();
  // Compute the strides of the new MDArrayView
  _strides[axis] *= bounds.step();
  // Compute the new MDArrayView and pointer
  _array = parent._array.view(offset,
                              computeSize(_dimensions(),
                                          _strides()));
  _ptr = _array.getRawPtr();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayView< T > &
MDArrayView< T >::operator=(const MDArrayView< T > & array)
{
  _dimensions = array._dimensions;
  _strides    = array._strides;
  _array      = array._array;
  _layout     = array._layout;
  _ptr        = array._ptr;
  _next_axis  = array._next_axis;
  return *this;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayView< T >::~MDArrayView()
{
}

////////////////////////////////////////////////////////////////////////

template< typename T >
int
MDArrayView< T >::numDims() const
{
  return _dimensions.size();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const Teuchos::Array< dim_type > &
MDArrayView< T >::dimensions() const
{
  return _dimensions;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
dim_type
MDArrayView< T >::dimension(int axis) const
{
#ifdef HAVE_DOMI_ARRAY_BOUNDSCHECK
  assertAxis(axis);
#endif
  return _dimensions[axis];
}

////////////////////////////////////////////////////////////////////////

template< typename T >
size_type
MDArrayView< T >::size()
{
  return computeSize(_dimensions(), _strides());
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const Teuchos::Array< size_type > &
MDArrayView< T >::strides() const
{
  return _strides;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const Teuchos::ArrayView< T > &
MDArrayView< T >::arrayView() const
{
  return _array;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
Layout
MDArrayView< T >::layout() const
{
  return _layout;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool
MDArrayView< T >::contiguous() const
{
  // Temporarily compute the strides this MDArrayView would have if
  // its memory were contiguous with no stride gaps
  Teuchos::Array< size_type > contig_strides =
    computeStrides< size_type, dim_type >(_dimensions, _layout);
  // If these strides are the same as the actual strides, then the
  // MDArrayView is contiguous
  return (contig_strides == _strides);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayView< T >::iterator
MDArrayView< T >::begin()
{
  return iterator(*this);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayView< T >::iterator
MDArrayView< T >::end()
{
  // Return the iterator corresponding to the last element
  return iterator(*this, true);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayView< T >::const_iterator
MDArrayView< T >::begin() const
{
  return const_iterator(getConst());
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayView< T >::const_iterator
MDArrayView< T >::end() const
{
  // Return the iterator corresponding to the last element
  return const_iterator(getConst(), true);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayView< T >::const_iterator
MDArrayView< T >::cbegin() const
{
  return const_iterator(getConst());
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayView< T >::const_iterator
MDArrayView< T >::cend() const
{
  // Return the iterator corresponding to the last element
  return const_iterator(getConst(), true);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayView< T >::reverse_iterator
MDArrayView< T >::rbegin()
{
  return reverse_iterator(*this);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayView< T >::reverse_iterator
MDArrayView< T >::rend()
{
  // Return the reverse_iterator corresponding to the last element
  return reverse_iterator(*this, true);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayView< T >::const_reverse_iterator
MDArrayView< T >::crbegin() const
{
  return const_reverse_iterator(*this);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
typename MDArrayView< T >::const_reverse_iterator
MDArrayView< T >::crend() const
{
  // Return the reverse_iterator corresponding to the last element
  return const_reverse_iterator(*this, true);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayView< const T >
MDArrayView< T >::getConst() const
{
  return MDArrayView< const T >(_array.getConst(),
                                _dimensions,
                                _strides,
                                _layout);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayView< T >
MDArrayView< T >::operator[](dim_type i)
{
  // Construct the new MDArrayView
  MDArrayView< T > result(*this, _next_axis, i);
  // Correct the next axis of the new MDArrayView
  if (_next_axis < result.numDims())
    result._next_axis = _next_axis;
  // Return the result
  return result;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const MDArrayView< T >
MDArrayView< T >::operator[](dim_type i) const
{
  // Construct the new MDArrayView
  MDArrayView< T > result(*this, _next_axis, i);
  // Correct the next axis of the new MDArrayView
  if (_next_axis < result.numDims())
    result._next_axis = _next_axis;
  // Return the result
  return result;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
MDArrayView< T >
MDArrayView< T >::operator[](Slice s)
{
  // Construct the new MDArrayView
  MDArrayView< T > result(*this, _next_axis, s);
  // Correct the next axis of the new MDArrayView
  result._next_axis = _next_axis + 1;
  if (result._next_axis >= _dimensions.size())
    result._next_axis = 0;
  // Return the result
  return result;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const MDArrayView< T >
MDArrayView< T >::operator[](Slice s) const
{
  // Construct the new MDArrayView
  MDArrayView< T > result(*this, _next_axis, s);
  // Correct the next axis of the new MDArrayView
  result._next_axis = _next_axis + 1;
  if (result._next_axis >= _dimensions.size())
    result._next_axis = 0;
  // Return the result
  return result;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
T &
MDArrayView< T >::operator()(dim_type i)
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
MDArrayView< T >::operator()(dim_type i,
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
MDArrayView< T >::operator()(dim_type i,
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
MDArrayView< T >::operator()(dim_type i,
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
MDArrayView< T >::operator()(dim_type i,
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
MDArrayView< T >::operator()(dim_type i,
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
MDArrayView< T >::operator()(dim_type i) const
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
MDArrayView< T >::operator()(dim_type i,
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
MDArrayView< T >::operator()(dim_type i,
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
MDArrayView< T >::operator()(dim_type i,
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
MDArrayView< T >::operator()(dim_type i,
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
MDArrayView< T >::operator()(dim_type i,
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
MDArrayView< T >::assign(const T & value)
{
  for (iterator it = begin(); it != end(); ++it)
    *it = value;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
T &
MDArrayView< T >::at(dim_type i, ...)
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
MDArrayView< T >::at(dim_type i, ...) const
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
bool
MDArrayView< T >::hasBoundsChecking()
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
MDArrayView< T >::toString() const
{
  // Call the private version of toString() with an indentation level
  // of zero
  return toString(0);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
std::string
MDArrayView< T >::toString(int indent) const
{
  std::stringstream ss;
  ss << "[";
  for (size_type i = 0; i < _dimensions[0]; i++)
  {
    // If 1D, output an element, else output a sub-array recursively
    if (numDims() == 1)
      ss << operator()(i);
    else
      ss << operator[](i).toString(indent+1);
    // If not the last item along this axis, output a comma
    if (i < _dimensions[0]-1)
    {
      ss << ",";
      // If 1D, follow the comma with a space, else a newline and
      // indentation spacing
      if (numDims() == 1)
	ss << " ";
      else
      {
	ss << std::endl << " ";
	for (int ii = 0; ii < indent; ii++) ss << " ";
      }
    }
  }
  ss << "]";
  return ss.str();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
const T *
MDArrayView< T >::getRawPtr() const
{
  return _array.getRawPtr();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
T *
MDArrayView< T >::getRawPtr()
{
  return _array.getRawPtr();
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool operator==(const MDArrayView< T > & a1, const MDArrayView< T > & a2)
{
  if (a1._dimensions != a2._dimensions) return false;
  if (a1._layout     != a2._layout    ) return false;
  typename MDArrayView< T >::const_iterator it1 = a1.begin();
  typename MDArrayView< T >::const_iterator it2 = a2.begin();
  for ( ; it1 != a1.end() && it2 != a2.end(); ++it1, ++it2)
  {
    if (*it1 != *it2) return false;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////

template< typename T >
bool operator!=(const MDArrayView< T > & a1, const MDArrayView< T > & a2)
{
  return not (a1 == a2);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
void swap(MDArrayView< T > & a1, MDArrayView< T > & a2)
{
  a1.swap(a2);
}

////////////////////////////////////////////////////////////////////////

template< typename T >
std::ostream & operator<<(std::ostream & os, const MDArrayView< T > & a)
{
  os << a.toString();
  return os;
}

/////////////////////////////
// Private implementations //
/////////////////////////////

template< typename T >
void
MDArrayView< T >::assertAxis(int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(0 <= axis && axis < _dimensions.size()),
    RangeError,
    "MDArrayView<T>::assertAxis(axis=" << axis << "): out of "
    << "range axis in [0, " << _dimensions.size() << ")"
    );
}

template< typename T >
void
MDArrayView< T >::assertIndex(dim_type i, int axis) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(0 <= i && i < _dimensions[axis]), RangeError,
    "MDArrayView<T>::assertIndex(i=" << i << ",axis=" << axis << "): out of "
    << "range i in [0, " << _dimensions[axis] << ")"
    );
}

}  // End namespace Domi

#endif  // DOMI_MDARRAYVIEW_HPP
