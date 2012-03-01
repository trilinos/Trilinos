#ifndef _ZOLTAN2_STRIDEDINPUT_HPP_
#define _ZOLTAN2_STRIDEDINPUT_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Environment.hpp>
#include <typeinfo>

/*! \file Zoltan2_StridedInput.hpp
 *  \brief This file defines the StridedInput class.
 */

namespace Zoltan2{

/*!  \brief The StridedInput class manages lists of weights or coordinates.
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
class StridedInput {
private:
  RCP<const Environment> env_;
  ArrayView<const scalar_t> vec_;
  int stride_;

public:

  /*! \brief Constructor
   *
   *  \param x  x[0] is the first element of the array.  The subsequent
   *               elements are at x[i * \c stride].
   *  \stride   the stride of the elements in the strided array.
   */
  StridedInput(ArrayView<const scalar_t> x, lno_t stride) :  
    vec_(x), stride_(stride) { }

  /*! \brief Default constructor
   */
  StridedInput(): vec_(), stride_(0) { }

  /*! \brief Access an element of the input array. 
   *
   *   \idx  The logical index of the element in the strided array. 
   *
   *    For performance, this is inline and no error checking.
   */
  scalar_t operator[](lno_t idx) const { return vec_[idx*stride_]; }

  /*! \brief Create a contiguous array of the required type, perhaps for a TPL.
   *
   *  TODO: if are there particular conversions we would like to do
   *   for TPLs, we can add methods to do that.  Here we just
   *   essentially cast.  If the cast is not valid (like double to float)
   *   an exception is thrown.
   */

  template <typename T> void getInputArray(ArrayRCP<const T> &array);

  /*! \brief The raw input information.
      \param len on return is the length of storage at \c vec.  This will
                be some multiple of stride.
      \param vec  on return is a pointer to the input.  Element k is at vec[k*stride].
      \param stride is describes the layout of the input in \c vec.
          
   */
  void getStridedList(size_t &len, const scalar_t *&vec, int &stride) const
  {
    len = vec_.size();
    vec = vec_.getRawPtr();
    stride = stride_;
  }

  /*! \brief Assignment operator
   */
  StridedInput & operator= (const StridedInput &sInput)
  {
    if (this != &sInput){
      size_t length;
      const scalar_t *vec;
      sInput.getStridedList(length, vec, stride_);
      vec_ = ArrayView<const scalar_t>(vec, length);
    }

    return *this;
  }
};

template<typename lno_t, typename scalar_t>
  template<typename T>
     void StridedInput<lno_t, scalar_t>::getInputArray(ArrayRCP<const T> &array)
{
  if (vec_.size() < 1){
    array = ArrayRCP<const T>();
  }
  else if (stride_==1 && typeid(T()) == typeid(scalar_t())){
    array = Teuchos::arcpFromArrayView<const T>(vec_);
  }
  else{
    Environment env;           // a default environment for error reporting
    size_t n = vec_.size() / stride_;
    T *tmp = new T [n];
    env.localMemoryAssertion(__FILE__, __LINE__, n, tmp);
    for (lno_t i=0,j=0; i < n; i++,j+=stride_){
      tmp[i] = Teuchos::as<T>(vec_[j]);
    }
    array = arcp(tmp, 0, n, true);
  }
  
  return;
}

}  // namespace Zoltan2

#endif
