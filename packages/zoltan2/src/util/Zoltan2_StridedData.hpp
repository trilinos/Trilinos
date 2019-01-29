// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef _ZOLTAN2_STRIDEDDATA_HPP_
#define _ZOLTAN2_STRIDEDDATA_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Environment.hpp>
#include <typeinfo>

/*! \file Zoltan2_StridedData.hpp
 *  \brief This file defines the StridedData class.
 */

namespace Zoltan2{

/*!  \brief The StridedData class manages lists of weights or coordinates.
 *
 * A likely representation for multi-dimensional weights or coordinates is
 *  an array ordered by identifier by dimension, or vice versa. The purposes 
 *  of this class:
 *
 *   \li to make it easy for applications to supply these arrays to Zoltan2.
 *   \li to provide a [] operator for algorithm's access to strided arrays.
 *   \li to create contiguous arrays of weights or coordinates if needed by 
 *         the algorithm.
 *   \li to create a contiguous array of weights using the data type required
 *        by the algorithm.
 *
 * If we need to modify weights (say convert floating point weights to
 * integers for a third party library) methods to do this could be added here.
 */

template<typename lno_t, typename scalar_t>
class StridedData {
private:
  ArrayRCP<const scalar_t> vec_;
  int stride_;

public:

  /*! \brief Constructor
   *
   *  \param x  x[0] is the first element of the array.  The subsequent
   *               elements are at x[i * \c stride].
   *  \param stride   the stride of the elements in the strided array.
   */
  StridedData(ArrayRCP<const scalar_t> x, int stride) :  
    vec_(x), stride_(stride)
    { }

  /*! \brief Default constructor.  A zero-length strided array.
   */
  StridedData(): vec_(), stride_(0)
    { }

  /*! \brief Return the length of the strided array.
   *
   *  The length may be longer than the number of values because the
   *  stride may be greater than one.
   */
  lno_t size() const { return vec_.size(); }

  /*! \brief Access an element of the input array. 
   *
   *   \param idx  The logical index of the element in the strided array. 
   *
   *    For performance, this is inline and no error checking.
   */
  scalar_t operator[](lno_t idx) const { return vec_[idx*stride_]; }

  /*! \brief Create a contiguous array of the required type, 
                 perhaps for a TPL.
   *
   *  TODO: if are there particular conversions we would like to do
   *   for TPLs, we can add methods to do that.  Here we just
   *   essentially cast.  If the cast is not valid (like double to float)
   *   an exception is thrown.
   */

  template <typename T> void getInputArray(ArrayRCP<const T> &array) const;

  /*! \brief Get a reference counted pointer to the input.
      \param vec  on return is a reference counted pointer to the input.  
                       Element \c k is at <tt>vec[k*stride]</tt>.
      \param stride is describes the layout of the input in \c vec.
   */
  void getStridedList(ArrayRCP<const scalar_t> &vec, int &stride) const
  {
    vec = vec_;
    stride = stride_;
  }

  /*! \brief Get the raw input information.
      \param len on return is the length of storage at \c vec.  This will
                be some multiple of stride.
      \param vec  on return is a pointer to the input.
                  Element \c k is at <tt>vec[k*stride]</tt>.
      \param stride is describes the layout of the input in \c vec.
   */
  void getStridedList(size_t &len, const scalar_t *&vec, int &stride) const
  {
    len = vec_.size();
    if (len != 0) vec = vec_.getRawPtr();
    else          vec = NULL;
    stride = stride_;
  }

  /*! \brief Assignment operator
   */
  StridedData & operator= (const StridedData &sInput)
  {
    if (this != &sInput)
      sInput.getStridedList(vec_, stride_);

    return *this;
  }
};

// Helper function needed for T=scalar_t specialization
// Separate function needed // because, with T != scalar_t, 
// "array = vec_" // would not compile; 
// ArrayRCP does not overload "=" operator for different types.
// Separate helper function needed (outside StridedData class)
// because cannot specialize member function of templated class.
template<typename scalar_t, typename T>
static void getInputArrayHelper(
  ArrayRCP<const T> &target,
  const ArrayRCP<const scalar_t> &src)
{
  // Create a copy of desired type T
  // From logic in getInputArray, we know stride == 1.
  size_t n = src.size();
  T *tmp = new T [n];

  if (!tmp){
    std::cerr << "Error: " << __FILE__ << ", " << __LINE__<< std::endl;
    std::cerr << n << " objects" << std::endl;
    throw std::bad_alloc();
  }

  for (size_t i=0; i < n; i++){
    tmp[i] = static_cast<T>(src[i]);
  }
  target = arcp(tmp, 0, n);
}

// Specialization with T == scalar_t:  just copy ArrayRCP
template<typename scalar_t>
static void getInputArrayHelper(
  ArrayRCP<const scalar_t> &target,
  const ArrayRCP<const scalar_t> &src)
{
  target = src;
}

// Member function for getting unstrided view/copy of StridedData.
template<typename lno_t, typename scalar_t>
  template<typename T>
     void StridedData<lno_t, scalar_t>::getInputArray(
       ArrayRCP<const T> &array) const
{
  if (vec_.size() < 1){
    array = ArrayRCP<const T>();
  }
  else if (stride_ > 1) {
    // Create an unstrided copy
    size_t n = vec_.size() / stride_;
    T *tmp = new T [n];

    if (!tmp){
      std::cerr << "Error: " << __FILE__ << ", " << __LINE__<< std::endl;
      std::cerr << n << " objects" << std::endl;
      throw std::bad_alloc();
    }

    for (size_t i=0,j=0; i < n; i++,j+=stride_){
      tmp[i] = static_cast<T>(vec_[j]);
    }
    array = arcp(tmp, 0, n);
  }
  else { // stride == 1
    Zoltan2::getInputArrayHelper<scalar_t, T>(array, vec_);
  }
  return;
}

}  // namespace Zoltan2

#endif
