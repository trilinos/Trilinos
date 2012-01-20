#ifndef _ZOLTAN2_STRIDEDINPUT_HPP_
#define _ZOLTAN2_STRIDEDINPUT_HPP_

#include <Zoltan2_Standards.hpp>
#include <typeinfo>

/*! \file Zoltan2_StridedInput.hpp
 *  \brief This file defines the StridedInput class.
 */

namespace Zoltan2{

/*!  *  \brief The StridedInput class manages lists of weights or coordinates.
 *
 * A likely representation for multi-dimensional weights or coordinates is
 *  an array ordered by identifier by dimension. The purposes of this class:
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
   *    x[0] is the first element of the array.  The subsequent
   *  elements are at x[i*stride].
   */
  StridedInput(RCP<const Environment> env, ArrayView<const scalar_t> x, 
    lno_t stride) :  env_(env), vec_(x), stride_(stride) 
  {
  }

  /*! \brief Constructor
   */
  StridedInput(): env_(rcp(new Environment)), vec_(), stride_(0) { }

  /*! \brief Access an element of the input array. 
   *
   *   For performance, no error checking.
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
  void getStridedList(size_t &len, const scalar_t *&vec, int &stride)
  {
    len = vec_.size();
    vec = vec_.getRawPtr();
    stride = stride_;
  }

#if 0
  /*! \brief Return the Environment object for the list.
      \return The Environment with which the list was contructed.
   */
  RCP<const Environment> getEnv() { return env_;}

  /*! \brief Assignment operator
   */
  StridedInput & operator= (const StridedInput &sInput)
  {
    if (this != &sInput){
      env_ = sInput.getEnv();
      size_t length;
      const scalar_t *vec;
      sInput.getStridedList(length, vec, stride_);
      vec_ = ArrayView(vec, length);
    }

    return *this;
  }
#endif


};

template<typename lno_t, typename scalar_t>
  template<typename T>
     void StridedInput<lno_t, scalar_t>::getInputArray(ArrayRCP<const T> &array)
{
  size_t n = vec_.size();

  if (n < 1){
    array = ArrayRCP<const T>();
  }
  else if (stride_==1 && typeid(T()) == typeid(scalar_t())){
    array = Teuchos::arcpFromArrayView<const T>(vec_);
  }
  else{
    T *tmp = new T [n];
    Z2_LOCAL_MEMORY_ASSERTION(*env_, n, tmp);
    for (lno_t i=0,j=0; i < n; i++,j+=stride_){
      tmp[i] = Teuchos::as<T>(vec_[j]);
    }
    array = arcp(tmp, 0, n, true);
  }
  
  return;
}

}  // namespace Zoltan2

#endif
