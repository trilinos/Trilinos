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

/*! \file Zoltan2_BasicVectorInput.hpp
    \brief Defines the BasicVectorInput class. 
*/

#ifndef _ZOLTAN2_BASICVECTORINPUT_HPP_
#define _ZOLTAN2_BASICVECTORINPUT_HPP_

#include <Zoltan2_VectorInput.hpp>
#include <Zoltan2_StridedData.hpp>

namespace Zoltan2 {

  /*!  \brief BasicVectorInput represents a vector (plus optional weights)
                supplied by the user as pointers to strided arrays.

    Input adapters provide access for Zoltan2 to the user's data.  The
    methods in the interface must be defined by users.  Many built-in
    adapters are already defined for common data structures, such as
    Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t is the data type for weights and vector element values.
    \li \c lno_t is the integral data type used by Zoltan2 for local indices and local counts.
    \li \c gno_t is the integral data type used by Zoltan2 to represent global indices and global counts.
    \li \c gid_t is the data type used by the application for global Ids.  If the application's global Id data type is a Teuchos Ordinal, then \c gid_t and \c gno_t are the same.  Otherwise, the application global Ids will be mapped to Teuchos Ordinals for use by Zoltan2 internally.  (Teuchos Ordinals are those data types for which traits are defined in Trilinos/packages/teuchos/src/Teuchos_OrdinalTraits.hpp.)
    \li \c node_t is a sub class of Kokkos::StandardNodeMemoryModel, which is used to optimize performance on many-core and multi-core architectures.  If you don't use Kokkos, you can ignore this data type.

    The template parameter (\c User) is a C++ class type which provides the
    actual data types with which the Zoltan2 library will be compiled, through
    a Traits mechanism.  \c User may be the
    actual class used by application to represent a vector, or it may be
    the empty helper class \c BasicUserTypes with which a Zoltan2 user
    can easily supply the data types for the library.


    The \c scalar_t type, representing use data such as matrix values, is
    used by Zoltan2 for weights, coordinates, part sizes and
    quality metrics.
    Some User types (like Tpetra::CrsMatrix) have an inherent scalar type,
    and some
    (like Tpetra::CrsGraph) do not.  For such objects, the scalar type is
    set by Zoltan2 to \c float.  If you wish to change it to double, set
    the second template parameter to \c double.


    BasicVectorInput may be a single vector or a set of corresponding vectors
    which have with the
    same global identifiers and the same distribution across processes.


 \todo Global identifiers should be optional.  If the user gives us
    gids in the input adapter, we will include them in the solution.
  \todo Is there any reason to specify coordinates for vector elements?

*/

template <typename User>
  class BasicVectorInput : public VectorInput<User> {

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef VectorInput<User>   base_adapter_t;
  typedef User user_t;

#endif

  /*! \brief Constructor for one vector with no weights.
   *
   *  \param numIds  the local length of the vector
   *  \param ids     pointer to the global ids of the local vector elements
   *  \param elements  pointer to the elements corresponding to the ids
   *  \param elementStride  the k'th element is at elements[k*elementStride]
   *         and elements is of lenth at least <tt>numIds * elementStride</tt>.
   *  
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this InputAdapter.
   */

  BasicVectorInput(lno_t numIds, const gid_t *ids, const scalar_t *elements, 
    int elementStride=1):
      env_(rcp(new Environment)), 
      numIds_(numIds), globalNumIds_(0), idList_(ids),
      numVectors_(1), elements_(), numWeights_(0), weights_()
  {
    vector<const scalar_t *> values;
    vector<int> strides;
    vector<const scalar_t *> emptyValues;
    vector<int> emptyStrides;

    values.push_back(elements);
    strides.push_back(elementStride);

    createBasicVector(values, strides, emptyValues, emptyStrides);
  }

  /*! \brief Constructor for one vector.
   *
   *  \param numIds  the local length of the vector
   *  \param ids     pointer to the global ids of the local vector elements
   *  \param elements  pointer to the elements corresponding to the ids
   *  \param elementStrides  the k'th element is at elements[k*elementStride]
   *                 and elements is of lenth numIds * elementStride
   *  \param weights  a list of pointers to arrays of weights.
   *      The number of weights per element is assumed to be
   *      \c weights.size().
   *  \param weightStrides  a list of strides for the \c weights.
   *     The weight for weight dimension \c n for \c ids[k] should be
   *     found at <tt>weights[n][weightStrides[n] * k]</tt>.
   *     If \c weightStrides.size() is zero, it is assumed all strides are one
   *  
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this InputAdapter.
   */

  BasicVectorInput(lno_t numIds, const gid_t *ids, 
    const scalar_t *elements, int elementStride,
    vector<const scalar_t *> &weights, vector<int> &weightStrides):
      env_(rcp(new Environment)), 
      numIds_(numIds), globalNumIds_(0), idList_(ids),
      numVectors_(1), elements_(), 
      numWeights_(weights.size()), weights_()
  {
    vector<const scalar_t *> values;
    vector<int> strides;

    values.push_back(elements);
    strides.push_back(elementStride);

    createBasicVector(values, strides, weights, weightStrides);
  }

  /*! \brief Constructor for a set of vectors.
   *
   *  \param numIds  the local length of the vectors
   *  \param ids     a pointer to the global ids of the local vector elements
   *  \param elements a list of pointers to the vector elements
   *          corresponding to the \c numIds ids.  The number of vectors
   *          assumed to be \c elements.size().
   *  \param elementStrides The strides for the \c elements list.
   *           The vector element for vector \c n for \c ids[k] should be
   *           found at <tt>elements[n][elementStrides[n] * k]</tt>.
   *           If \c elementStrides.size() is zero, it is assumed
   *           all strides are one.
   *  \param weights  a list of pointers to arrays of weights.
   *      The number of weights per vector element is assumed to be
   *      \c weights.size().
   *  \param weightStrides  a list of strides for the \c weights.
   *     The weight for weight dimension \c n for \c ids[k] should be
   *     found at <tt>weights[n][weightStrides[n] * k]</tt>.
   *     If \c weightStrides.size() is zero, it is assumed all strides are one.
   *  
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this InputAdapter.
   */

  BasicVectorInput(lno_t numIds, const gid_t *ids, 
    vector<const scalar_t *> &elements,  vector<int> &elementStrides,
    vector<const scalar_t *> &weights, vector<int> &weightStrides):
      env_(rcp(new Environment)), 
      numIds_(numIds), globalNumIds_(0), idList_(ids),
      numVectors_(elements.size()), elements_(),
      numWeights_(weights.size()), weights_()
  {
    createBasicVector(elements, elementStrides, weights, weightStrides);
  }

  /*! Destructor
   */
  ~BasicVectorInput() {};

  ////////////////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////////////////

  string inputAdapterName() const {return string("BasicVector");}

  size_t getLocalNumberOfObjects() const { return numIds_;}

  int getNumberOfWeightsPerObject() const { return numWeights_;}

  size_t getObjectWeights(int dim, const scalar_t *&wgt, int &stride) const
  {
    return getVectorWeights(dim, wgt, stride);
  }

  ////////////////////////////////////////////////////
  // The VectorInput interface.
  ////////////////////////////////////////////////////

  int getNumberOfVectors() const { return numVectors_;}

  int getNumberOfWeights() const { return numWeights_;}

  size_t getLocalLength() const { return numIds_; }

  size_t getGlobalLength() const { return globalNumIds_;}

  size_t getVector(const gid_t *&ids, 
     const scalar_t *&element, int &stride) const
  {
    return getVector(0, ids, element, stride);
  }

  size_t getVector(int i, const gid_t *&ids, 
     const scalar_t *&element, int &stride) const;

  size_t getVectorWeights(int dimension, 
     const scalar_t *&weights, int &stride) const;

private:

  // A default environment.  An Environment is an internal Zoltan2
  // class, so input adapters don't usually have one.  But we create
  // one here so we can use it for error handling.

  RCP<const Environment> env_;

  lno_t numIds_;
  gno_t globalNumIds_;

  const gid_t *idList_;

  int numVectors_;
  ArrayRCP<StridedData<lno_t, scalar_t> > elements_ ;

  int numWeights_;
  ArrayRCP<StridedData<lno_t, scalar_t> > weights_;

  void createBasicVector(
    vector<const scalar_t *> &elements,  vector<int> &elementStrides,
    vector<const scalar_t *> &weights, vector<int> &weightStrides);

};

////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////

template <typename User>
  size_t BasicVectorInput<User>::getVector(int i, const gid_t *&ids, 
    const scalar_t *&element, int &stride) const
{
  env_->localInputAssertion(__FILE__, __LINE__, "invalid vector number",
    i >= 0 && i < numVectors_, BASIC_ASSERTION);
  
  ids = idList_;

  size_t length;

  elements_[i].getStridedList(length, element, stride);

  return length;
}

template <typename User>
  size_t BasicVectorInput<User>::getVectorWeights(int dimension, 
    const scalar_t *&weights, int &stride) const
{
  env_->localInputAssertion(__FILE__, __LINE__,  "invalid weight dimension",
    dimension >= 0 && dimension < numWeights_, BASIC_ASSERTION);

  size_t length;

  weights_[dimension].getStridedList(length, weights, stride);

  return length;
}

template <typename User>
  void BasicVectorInput<User>::createBasicVector(
   vector<const scalar_t *> &elements,  vector<int> &elementStrides,
   vector<const scalar_t *> &weights, vector<int> &weightStrides)
{
  typedef StridedData<lno_t,scalar_t> input_t;

  gno_t tmp = numIds_;
  try{
    reduceAll<int, gno_t>(*(env_->comm_), Teuchos::REDUCE_SUM, 1, 
       &tmp, &globalNumIds_);
  }
  Z2_THROW_OUTSIDE_ERROR(*env_);

  if (numIds_){
    int stride = 1;
    elements_ = arcp(new input_t [numVectors_], 0, numVectors_, true);
    for (int v=0; v < numVectors_; v++){
      if (elementStrides.size())
        stride = elementStrides[v];
      ArrayRCP<const scalar_t> eltV(elements[v], 0, stride*numIds_, false); 
      elements_[v] = input_t(eltV, stride);
    }

    if (numWeights_){
      stride = 1;
      weights_ = arcp(new input_t [numWeights_], 0, numWeights_, true);
      for (int w=0; w < numWeights_; w++){
        if (weightStrides.size())
          stride = weightStrides[w];
        ArrayRCP<const scalar_t> wgtV(weights[w], 0, stride*numIds_, false); 
        weights_[w] = input_t(wgtV, stride);
      }
    }
  }
}
  
}  //namespace Zoltan2
  
#endif
