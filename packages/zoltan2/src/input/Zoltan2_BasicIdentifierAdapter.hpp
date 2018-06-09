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

/*! \file Zoltan2_BasicIdentifierAdapter.hpp
    \brief Defines the BasicIdentifierAdapter class.
*/

#ifndef _ZOLTAN2_BASICIDENTIFIERADAPTER_HPP_
#define _ZOLTAN2_BASICIDENTIFIERADAPTER_HPP_

#include <Zoltan2_IdentifierAdapter.hpp>
#include <Zoltan2_StridedData.hpp>

namespace Zoltan2 {

/*! \brief This class represents a collection of global Identifiers
 *           and their associated weights, if any.
 *
 *  The user supplies the identifiers and weights by way of pointers
 *    to arrays.  
 *
 *  The template parameter (\c User) is a C++ class type which provides the
 *  actual data types with which the Zoltan2 library will be compiled, through
 *  a Traits mechanism.  \c User may be the
 *  actual class used by application to represent coordinates, or it may be
 *  the empty helper class \c BasicUserTypes with which a Zoltan2 user
 *  can easily supply the data types for the library.
 *
 *  The \c scalar_t type, representing use data such as matrix values, is
 *  used by Zoltan2 for weights, coordinates, part sizes and
 *  quality metrics.
 *  Some User types (like Tpetra::CrsMatrix) have an inherent scalar type,
 *  and some
 *  (like Tpetra::CrsGraph) do not.  For such objects, the scalar type is
 *  set by Zoltan2 to \c float.  If you wish to change it to double, set
 *  the second template parameter to \c double.
 */

template <typename User>
  class BasicIdentifierAdapter: public IdentifierAdapter<User> {

public:
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::part_t   part_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef User user_t;

  /*! \brief Constructor
   *  \param numIds is the number of identifiers in the list
   *  \param ids should point to a list of numIds identifiers.
   *  \param weights  a list of pointers to arrays of weights.
   *      The number of weights per identifier is assumed to be
   *      \c weights.size().
   *  \param weightStrides  a list of strides for the \c weights.
   *     The weight for weight index \c n for \c ids[k] should be
   *     found at <tt>weights[n][weightStrides[n] * k]</tt>.
   *     If \c weightStrides.size() is zero, it is assumed all strides are one.
   *
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this Adapter.
   */
  BasicIdentifierAdapter(lno_t numIds, const gno_t *idPtr, 
    std::vector<const scalar_t *> &weights, std::vector<int> &weightStrides);

  /*! \brief Constructor
   *  \param numIds is the number of identifiers in the list
   *  \param ids should point to a list of numIds identifiers.
   *
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this Adapter.
   */
  BasicIdentifierAdapter(lno_t numIds, const gno_t *idPtr):
      numIds_(numIds), idList_(idPtr), weights_() {}

  ////////////////////////////////////////////////////////////////
  // The Adapter interface.
  ////////////////////////////////////////////////////////////////

  size_t getLocalNumIDs() const { return numIds_; }

  void getIDsView(const gno_t *&Ids) const { Ids = idList_; }

  int getNumWeightsPerID() const { return weights_.size(); }

  void getWeightsView(const scalar_t *&weights, int &stride, int idx) const {
    if (idx < 0 || idx >= weights_.size()) {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid weight index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
    }
    size_t length;
    weights_[idx].getStridedList(length, weights, stride);
  }

private:
  lno_t numIds_;
  const gno_t *idList_;
  ArrayRCP<StridedData<lno_t, scalar_t> > weights_;
};

////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////

template <typename User>
  BasicIdentifierAdapter<User>::BasicIdentifierAdapter(
    lno_t numIds, const gno_t *idPtr,
    std::vector<const scalar_t *> &weights, std::vector<int> &weightStrides):
      numIds_(numIds), idList_(idPtr), weights_()
{
  typedef StridedData<lno_t, scalar_t> input_t;
  size_t numWeights = weights.size();

  if (numWeights > 0){
    weights_ = arcp(new input_t [numWeights], 0, numWeights, true);

    if (numIds > 0){
      for (size_t i = 0; i < numWeights; i++){
        int stride = weightStrides.size() ? weightStrides[i] : 1;
        ArrayRCP<const scalar_t> wgtV(weights[i], 0, stride * numIds, false);
        weights_[i] = input_t(wgtV, stride);
      }
    }
  }
}
  
}  //namespace Zoltan2
  
#endif
