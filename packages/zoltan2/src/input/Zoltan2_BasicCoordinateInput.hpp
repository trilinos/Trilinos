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

/*! \file Zoltan2_BasicCoordinateAdapter.hpp
    \brief Defines the BasicCoordinateAdapter class.
*/

#ifndef _ZOLTAN2_BASICCOORDINATEADAPTER_HPP_
#define _ZOLTAN2_BASICCOORDINATEADAPTER_HPP_

#include <Zoltan2_CoordinateInput.hpp>
#include <Zoltan2_StridedData.hpp>
#include <vector>

namespace Zoltan2 {

  /*!  \brief BasicCoordinateAdapter represents geometric coordinates that are
                supplied by the user as pointers to strided arrays.

    Adapters provide access from Zoltan2 to the user's data.  The
    methods in the interface must be defined by users.  Many built-in
    adapters are already defined for common data structures, such as
    Tpetra and Epetra objects and C-language pointers to arrays.

    Data types:
    \li \c scalar_t is the data type for weights and coordinates
    \li \c lno_t is the integral data type used by Zoltan2 for local indices and local counts.
    \li \c gno_t is the integral data type used by Zoltan2 to represent global indices and global counts.
    \li \c gid_t is the data type used by the application for global Ids.  If the application's global Id data type is a Teuchos Ordinal, then \c gid_t and \c gno_t are the same.  Otherwise, the application global Ids will be mapped to Teuchos Ordinals for use by Zoltan2 internally.  (Teuchos Ordinals are those data types for which traits are defined in Trilinos/packages/teuchos/src/Teuchos_OrdinalTraits.hpp.)
    \li \c node_t is a sub class of Kokkos::StandardNodeMemoryModel, which is used to optimize performance on many-core and multi-core architectures.  If you don't use Kokkos, you can ignore this data type.

    The template parameter (\c User) is a C++ class type which provides the
    actual data types with which the Zoltan2 library will be compiled, through
    a Traits mechanism.  \c User may be the
    actual class used by application to represent coordinates, or it may be
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
   
*/

template <typename User>
  class BasicCoordinateAdapter : public CoordinateAdapter<User> {

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  typedef typename InputTraits<User>::scalar_t    scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef CoordinateAdapter<User>   base_adapter_t;
  typedef User user_t;

#endif

  /*! \brief Constructor for dimension 1, 2 or 3 and no weights.
   *
   * \param numIds The number of local coordinates.
   * \param ids    The global identifiers for the coordinates.
   * \param x      A pointer to the first dimension of the coordinates.
   * \param y      A pointer to the second dimension, if any.
   * \param z      A pointer to the third dimension, if any.
   * \param xStride  The stride for the \c x list.  The \x coordinate
   *          for point \c ids[n]  should be found at <tt>x[xStride * n]</tt>.
   * \param yStride  The stride for the \c y list.  The \y coordinate
   *          for point \c ids[n]  should be found at <tt>y[yStride * n]</tt>.
   * \param zStride  The stride for the \c z list.  The \z coordinate
   *          for point \c ids[n]  should be found at <tt>z[zStride * n]</tt>.
   *  
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this Adapter.
   */

  BasicCoordinateAdapter(lno_t numIds, const gid_t *ids,
    const scalar_t *x, const scalar_t *y, const scalar_t *z,
    int xStride=1, int yStride=1, int zStride=1);

  /*! \brief Constructor for arbitrary dimension with weights.
   *
   *  \param numIds   the local number of coordinates.
   *  \param ids     is a pointer to the coordinate global Ids. 
   *  \param values a list of pointers to the coordinate values
   *          corresponding to the \c numIds ids.  The coordinate
   *          dimension is taken to be \c values.size().
   *  \param valueStrides The strides for the \c values list.  
   *           The coordinate for dimension \c n for \c ids[k] should be
   *           found at <tt>values[n][valueStrides[n] * k]</tt>.
   *           If \c valueStrides.size() is zero, it is assumed 
   *           all strides are one.
   *  \param weights  a list of pointers to arrays of weights.  
   *      The number of weights per coordinate is assumed to be
   *      \c weights.size(). 
   *  \param weightStrides  a list of strides for the \c weights.
   *     The weight for weight dimension \c n for \c ids[k] should be
   *     found at <tt>weights[n][weightStrides[n] * k]</tt>.
   *     If \c weightStrides.size() is zero, it is assumed all strides are one.
   *  
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this Adapter.
   */

  BasicCoordinateAdapter(lno_t numIds, const gid_t *ids, 
    vector<const scalar_t *> &values,  vector<int> &valueStrides,
    vector<const scalar_t *> &weights, vector<int> &weightStrides);

  /*! Destructor
   */
  ~BasicCoordinateAdapter() {};

  ////////////////////////////////////////////////////////////////
  // The Adapter interface.
  ////////////////////////////////////////////////////////////////

  size_t getLocalNumberOfObjects() const { return numIds_;}

  int getNumberOfWeightsPerObject() const { return numWeights_;}

  size_t getObjectWeights(int idx, const scalar_t *&wgt, int &stride) const
  {
    return getCoordinateWeights(idx, wgt, stride);
  }

  ////////////////////////////////////////////////////
  // The CoordinateAdapter interface.
  ////////////////////////////////////////////////////

  int getCoordinateDimension() const { return dimension_;}

  int getNumberOfWeights() const { return numWeights_;}  

  size_t getLocalNumberOfCoordinates() const { return numIds_; }

  size_t getCoordinates(int dim, const gid_t *&gids, const scalar_t *&coords, 
    int &stride) const
  {
    if (dim < 0 || dim >= dimension_) {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid dimension " << dim << std::endl;
      throw std::runtime_error(emsg.str());
    }

    gids = idList_;
    
    size_t length;

    coords_[dim].getStridedList(length, coords, stride);

    return length;
  }

  size_t getCoordinateWeights(int idx, const scalar_t *&weights, 
    int &stride) const
  {
    if (idx < 0 || idx >= numWeights_) {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
    }
    
    size_t length;

    weights_[idx].getStridedList(length, weights, stride);

    return length;
  }

private:
  void initializeData(
    vector<const scalar_t *> &values,  vector<int> &valueStrides,
    vector<const scalar_t *> &weights, vector<int> &weightStrides);

  lno_t numIds_;
  const gid_t *idList_;

  int dimension_;
  ArrayRCP<StridedData<lno_t, scalar_t> > coords_;

  int numWeights_;
  ArrayRCP<StridedData<lno_t, scalar_t> > weights_;
};

/////////////////////////////////////////////////////////////////
// Definitions
/////////////////////////////////////////////////////////////////

template <typename User>
  BasicCoordinateAdapter<User>::BasicCoordinateAdapter( 
    lno_t numIds, const gid_t *ids,
    const scalar_t *x, const scalar_t *y, const scalar_t *z,
    int xStride, int yStride, int zStride):
      numIds_(numIds), idList_(ids), 
      dimension_(0), coords_(), 
      numWeights_(0), weights_()
{
  vector<const scalar_t *> values;
  vector<int> strides;
  vector<const scalar_t *> emptyValues;
  vector<int> emptyStrides;

  if (x){
    values.push_back(x);
    strides.push_back(xStride);
    dimension_++;
    if (y){
      values.push_back(y);
      strides.push_back(yStride);
      dimension_++;
      if (z){
        values.push_back(z);
        strides.push_back(zStride);
        dimension_++;
      }
    }
  }

  initializeData(values, strides, emptyValues, emptyStrides);
}

template <typename User>
  BasicCoordinateAdapter<User>::BasicCoordinateAdapter( 
    lno_t numIds, const gid_t *ids, 
    vector<const scalar_t *> &values,  vector<int> &valueStrides,
    vector<const scalar_t *> &weights, vector<int> &weightStrides):
      numIds_(numIds), idList_(ids), 
      dimension_(values.size()), coords_(),
      numWeights_(weights.size()), weights_()
{
  initializeData(values, valueStrides, weights, weightStrides);
}

template <typename User>
  void BasicCoordinateAdapter<User>::initializeData(
    vector<const scalar_t *> &values,  vector<int> &valueStrides,
    vector<const scalar_t *> &weights, vector<int> &weightStrides)
{
  typedef StridedData<lno_t,scalar_t> input_t;

  coords_ = arcp(new input_t [dimension_], 0, dimension_, true);

  if (numWeights_ > 0)
    weights_ = arcp(new input_t [numWeights_], 0, numWeights_, true);

  if (numIds_){
    int stride = 1;
    for (int x=0; x < dimension_; x++){
      if (valueStrides.size())
        stride = valueStrides[x];
      ArrayRCP<const scalar_t> coordV(values[x], 0, stride*numIds_, false); 
      coords_[x] = input_t(coordV, stride);
    }

    if (numWeights_){
      stride = 1;
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
