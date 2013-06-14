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

#ifndef DOMI_MDITERATOR_HPP
#define DOMI_MDITERATOR_HPP

namespace Domi
{

/** \brief Iterator class suitable for multi-dimensional arrays
 *
 * This iterator is specialized for multi-dimensional arrays.  It can
 * be used like a standard iterator, without knowledge of the number
 * of dimensions of the array.  It is designed to iterate only over
 * valid elements of the array.  It is possible, especially with
 * sliced views into a parent array, for the data buffer to have
 * stride gaps that do not belong to the array.  This iterator avoids
 * those stride gaps by keeping an internal record of the multi-
 * dimensional index, and only iterating over valid indices.
 *
 * To apply MDIterator to all three multi-dimensional array types
 * (MDArray, MDArrayView and MDArrayRCP), the class is templated on
 * parameter class MDARRAY, which is intended to be any one of the
 * three MDArray types.  The MDARRAY class is expected to support the
 * the typedefs Ordinal, size_type and value_type, the num_dims()
 * method and the _ptr, _strides, and _dimensions attributes.
 *
 * It is intended that the array class that will use the MDIterator
 * will declare the MDIterator to be a friend and to typedef the fully
 * qualified class.  For example, within the MDArray< T > class:
 *
 *   \code
 *   friend class MDIterator< MDArray< T > >;
 *   typedef MDIterator< MDArray< T > > iterator;
 *   \endcode
 *
 * and declare begin() and end() methods:
 *
 *   \code
 *   const iterator begin() const { return iterator(*this      ); }
 *   const iterator end()   const { return iterator(*this, true); }
 *   \endcode
 */
template< class MDARRAY >
class MDIterator
{
public:

  /** \name MDARRAY typedefs */
  //@{

  typedef typename MDARRAY::Ordinal    Ordinal;
  typedef typename MDARRAY::size_type  size_type;
  typedef typename MDARRAY::value_type value_type;

  //@}

  /** \name Constructors and Destructor */
  //@{

  /** \brief MDARRAY constructor
   *
   *  \param mdarray [in] The multi-dimensional array object on which
   *         will be iterated.
   *
   *  \param end_index [in] If true, set the internal index to the
   *         MDARRAY end() index.  If false, set the internal index to
   *         the MDARRAY begin() index.  Default false.
   *
   *  Produces an iterator with index corresponding to either the
   *  MDARRAY <tt>begin()</tt> or <tt>end()</tt> methods, depending on
   *  the value of the <tt>end_index</tt> argument.
   */
  MDIterator(const MDARRAY & mdarray,
             bool end_index = false) :
    _mdarray(mdarray),
    _index(mdarray.num_dims())
  {
    if (end_index)
      assign_end_index();
    else
      _index.assign(_mdarray.num_dims(), 0);
  }

  /** \brief Index constructor
   *
   *  \param mdarray [in] The multi-dimensional array object on which
   *         will be iterated.
   *
   *  \param index [in] A Teuchos::ArrayView that specifies where the
   *         internal index of the iterator should start.
   *
   *  Produces an iterator with index corresponding to the given
   *  index.
   */
  MDIterator(MDARRAY & mdarray,
             Teuchos::ArrayView< size_type > & index) :
    _mdarray(mdarray),
    _index(index)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      (_mdarray.num_dims() != _index.size()), Teuchos::RangeError,
      "Input array has " << _mdarray.num_dims() << " dimensions, while index "
      "has " << _index.size());
#ifdef DOMI_ENABLE_ABC
    for (_axis = 0; _axis < _index.size(); ++_axis)
      assert_index(_index[_axis], _axis);
#endif
  }

  /** \brief Copy constructor
   *
   * \param source [in] The source <tt>MDIterator</tt> to be copied
   */
  MDIterator(const MDIterator< MDARRAY > & source) :
    _mdarray(source._mdarray),
    _index(source._index)
  {
  }

  /** \brief Destructor
   */
  ~MDIterator() { }

  //@}

  /** \name Standard Operators */
  //@{

  /** \brief Assignment operator
   *
   *  \param source [in] Source iterator for assignment
   */
  MDIterator & operator=(const MDIterator< MDARRAY > & source)
  {
    _mdarray = source._mdarray;
    _index   = source._index;
  }

  /** \brief Equality operator
   *
   *  \param other [in] Iterator to be compared to
   */
  bool operator==(const MDIterator< MDARRAY > & other) const
  {
    // If underlying MDARRAYs are different, then return not equal
    if (_mdarray._ptr != other._mdarray._ptr) return false;
    // If any of the current index values differ, then return not equal 
    for (_axis = 0; _axis < _index.size(); _axis++)
      if (_index[_axis] != other._index[_axis]) return false;
    // Return equal
    return true;
  }

  /** \brief Inequality operator
   *
   *  \param other [in] Iterator to be compared to
   */
  bool operator!=(const MDIterator< MDARRAY > & other) const
  {
    return !(*this == other);
  }

  /** \brief Dereferencing operator */
  value_type & operator*()
  {
    Ordinal offset = 0;
    for (_axis=0; _axis < _index.size(); ++_axis)
      offset += _index[_axis] * _mdarray._strides[_axis];
    return _mdarray._ptr[offset];
  }

  /** \brief Prefix increment operator */
  MDIterator & operator++()
  {
    _axis = 0;
    _done = false;
    while (not _done)
    {
      _index[_axis]++;
      _done = (_index[_axis] < _mdarray._dimensions[_axis]);
      if (not _done)
      {
        _index[_axis] = 0;
        _axis++;
        if (_axis >= _index.size())
        {
          _done = true;
          assign_end_index();
        }
      }
    }
    return *this;
  }

  /** \brief Postfix increment operator */
  MDIterator operator++(int)
  {
    MDIterator result(*this);
    ++(*this);
    return result;
  }

  /** \brief Prefix decrement operator */
  MDIterator & operator--()
  {
    _axis = 0;
    _done = false;
    while (not _done)
    {
      _index[_axis]--;
      _done = (_index[_axis] >= 0);
      if (not _done)
      {
        _index[_axis] = _mdarray._dimensions[_axis];
        _axis++;
        if (_axis >= _index.size())
        {
          _done = true;
          assign_end_index();
        }
      }
    }
    return *this;
  }

  /** \brief Postfix decrement operator */
  MDIterator operator--(int)
  {
    MDIterator result(*this);
    --(*this);
    return result;
  }

  //@}

private:

  // A reference to the multi-dimensional array being iterated
  const MDARRAY & _mdarray;

  // The multi-dimensional index of the current iterate
  Teuchos::Array< size_type > _index;

  // A temporary value indicating the axis of the index currently
  // being incremented or decremented
  mutable size_type _axis;

  // A temporary value indicating whether an increment or decrement
  // operation is complete
  mutable bool _done;

  // We need an index that is recognized as the end index.  It must
  // not be a valid index for the MDARRAY.  Since there are a nearly
  // infinite number of indexes that could serve as the end index,
  // this method should always be used to assign the index of an end
  // iterator.
  void assign_end_index()
  {
    // We choose the end index to be equal to the MDARRAY dimensions,
    // where each index value is one greater than the largest valid
    // index for that axis.
    for (size_type axis = 0; axis < _index.size(); ++axis)
      _index[axis] = _mdarray._dimensions[axis];
  }

  // Assert that the given index is valid for the given axis
  void assert_index(size_type i, int axis) const
  {
  TEUCHOS_TEST_FOR_EXCEPTION(
    !(0 <= i && i < _mdarray._dimensions[axis]), Teuchos::RangeError,
    "MDIterator<MDARRAY>::assertIndex(i=" << i << ",axis=" << axis << "): out of "
    << "range i in [0, " << _mdarray._dimensions[axis] << ")"
    );
  }

};

}  // namespace Domi

#endif
