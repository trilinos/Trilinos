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

/*! \file Zoltan2_BasicVectorAdapter.hpp
    \brief Defines the BasicVectorAdapter class.
*/

#ifndef _ZOLTAN2_BASICVECTORADAPTER_HPP_
#define _ZOLTAN2_BASICVECTORADAPTER_HPP_

#include <Zoltan2_VectorAdapter.hpp>
#include <Zoltan2_StridedData.hpp>

namespace Zoltan2 {

/*!  \brief BasicVectorAdapter represents a vector (plus optional weights)
            supplied by the user as pointers to strided arrays.

    BasicVectorAdapter may be a single vector or multivector (set of 
    corresponding vectors with the same global identifiers and 
    distribution across processes).  A constructor specifically for use
    of BasicVectorAdapter to represent geometric coordinates is also provided.

    Data types:
    \li \c scalar_t is the data type for weights and vector entry values.
    \li \c lno_t is the integral data type used by Zoltan2 for local indices and local counts.
    \li \c gno_t is the integral data type used by Zoltan2 to represent global indices and global counts.
    \li \c zgid_t is the data type used by the application for global Ids.  If the application's global Id data type is a Teuchos Ordinal, then \c zgid_t and \c gno_t are the same.  Otherwise, the application global Ids will be mapped to Teuchos Ordinals for use by Zoltan2 internally.  (Teuchos Ordinals are those data types for which traits are defined in Trilinos/packages/teuchos/src/Teuchos_OrdinalTraits.hpp.)
    \li \c node_t is a sub class of KokkosClassic::StandardNodeMemoryModel, which is used to optimize performance on many-core and multi-core architectures.  If you don't use Kokkos, you can ignore this data type.

    The template parameter (\c User) is a C++ class type which provides the
    actual data types with which the Zoltan2 library will be compiled, through
    a Traits mechanism.  \c User may be the
    actual class used by application to represent a vector, or it may be
    the empty helper class \c BasicUserTypes with which a Zoltan2 user
    can easily supply the data types for the library.
*/

template <typename User>
  class BasicVectorAdapter : public VectorAdapter<User> {

public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t lno_t;
  typedef typename InputTraits<User>::gno_t gno_t;
  typedef typename InputTraits<User>::zgid_t zgid_t;
  typedef typename InputTraits<User>::part_t   part_t;
  typedef typename InputTraits<User>::node_t node_t;
  typedef VectorAdapter<User> base_adapter_t;
  typedef User user_t;

#endif

  /*! \brief Constructor for one vector with (optionally) one weight.
   *
   *  \param numIds  the local length of the vector
   *  \param ids     pointer to the global ids of the local vector entries
   *  \param entries  pointer to the entries corresponding to the ids
   *  \param entryStride  the k'th entry is at entries[k*entryStride]
   *         and entries is of length at least <tt>numIds * entryStride</tt>.
   *  \param usewgts flag indicating whether weights are provided on any process
   *  \param wgts    the number of weights per vector entry
   *  \param wgtStride  the weight for the k'th entry is at wgts[k*wgtStride]
   *         and wgts is of length at least <tt>numIds * wgtStride</tt>.
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this Adapter.
   */

  BasicVectorAdapter(lno_t numIds, const zgid_t *ids,
                     const scalar_t *entries, int entryStride=1,
                     bool usewgts=false,
                     const scalar_t *wgts=NULL, int wgtStride=1):
      numIds_(numIds), idList_(ids),
      numEntriesPerID_(1), entries_(),
      numWeights_(usewgts==true), weights_()
  {
    std::vector<const scalar_t *> values;
    std::vector<int> strides;
    std::vector<const scalar_t *> weightValues;
    std::vector<int> weightStrides;

    values.push_back(entries);
    strides.push_back(entryStride);
    if (usewgts) {
      weightValues.push_back(wgts);
      weightStrides.push_back(wgtStride);
    }

    createBasicVector(values, strides, weightValues, weightStrides);
  }

  /*! \brief Constructor for multivector (a set of vectors sharing the same
   *         global numbering and data distribution across processes).
   *
   *  \param numIds  the local length of each vector
   *  \param ids     a pointer to the global ids of the local vector entries
   *  \param entries a std::vector of pointers to the vector entries
   *          corresponding to the \c numIds ids.  The number of vectors
   *          assumed to be \c entries.size().
   *  \param entryStride The strides for \c entries.
   *           The vector entry for vector \c n for \c ids[k] should be
   *           found at <tt>entries[n][entryStride[n] * k]</tt>.
   *           If \c entryStride.size() is zero, it is assumed
   *           all strides are one.
   *  \param weights  a list of pointers to arrays of weights.
   *      The number of weights per vector entry is assumed to be
   *      \c weights.size().
   *  \param weightStrides  a list of strides for the \c weights.
   *     The weight for weight index \c n for \c ids[k] should be
   *     found at <tt>weights[n][weightStrides[n] * k]</tt>.
   *     If \c weightStrides.size() is zero, it is assumed all strides are one.
   *
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this Adapter.
   */

  BasicVectorAdapter(lno_t numIds, const zgid_t *ids,
    std::vector<const scalar_t *> &entries,  std::vector<int> &entryStride,
    std::vector<const scalar_t *> &weights, std::vector<int> &weightStrides):
      numIds_(numIds), idList_(ids),
      numEntriesPerID_(entries.size()), entries_(),
      numWeights_(weights.size()), weights_()
  {
    createBasicVector(entries, entryStride, weights, weightStrides);
  }

  /*! \brief A simple constructor for coordinate-based problems with
   *         dimension 1, 2 or 3 and (optionally) one weight per coordinate.
   *
   * \param numIds The number of local coordinates.
   * \param ids    The global identifiers for the coordinates.
   * \param x      A pointer to the x-dimension coordinates.
   * \param y      A pointer to the y-dimension coordinates, if any.
   * \param z      A pointer to the z-dimension coordinates, if any.
   * \param xStride  The stride for the \c x array.  The \x coordinate
   *          for point \c ids[n]  should be found at <tt>x[xStride * n]</tt>.
   *          Default = 1.
   * \param yStride  The stride for the \c y array.  The \y coordinate
   *          for point \c ids[n]  should be found at <tt>y[yStride * n]</tt>.
   *          Default = 1.
   * \param zStride  The stride for the \c z array.  The \z coordinate
   *          for point \c ids[n]  should be found at <tt>z[zStride * n]</tt>.
   *          Default = 1.
   * \param usewgts flag indicating whether weights are provided on any process.
   * \param wgts    the number of weights per vector entry
   * \param wgtStride the weight for the k'th coordinate is at wgts[k*wgtStride]
   *         and wgts is of length at least <tt>numIds * wgtStride</tt>.
   *
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this Adapter.
   */

  BasicVectorAdapter(lno_t numIds, const zgid_t *ids,
                     const scalar_t *x, const scalar_t *y,
                     const scalar_t *z,
                     int xStride=1, int yStride=1, int zStride=1,
                     bool usewgts=false, const scalar_t *wgts=NULL,
                     int wgtStride=1) :
      numIds_(numIds), idList_(ids), numEntriesPerID_(0), entries_(),
      numWeights_(usewgts==true), weights_()
  {
    std::vector<const scalar_t *> values, weightValues;
    std::vector<int> strides, weightStrides;

    if (x){
      values.push_back(x);
      strides.push_back(xStride);
      numEntriesPerID_++;
      if (y){
        values.push_back(y);
        strides.push_back(yStride);
        numEntriesPerID_++;
        if (z){
          values.push_back(z);
          strides.push_back(zStride);
          numEntriesPerID_++;
        }
      }
    }
    if (usewgts) {
      weightValues.push_back(wgts);
      weightStrides.push_back(wgtStride);
    }
    createBasicVector(values, strides, weightValues, weightStrides);
  }


  ~BasicVectorAdapter() {};

  ////////////////////////////////////////////////////////////////
  // The Adapter interface.
  ////////////////////////////////////////////////////////////////

  size_t getLocalNumIDs() const { return numIds_;}

  void getIDsView(const zgid_t *&ids) const {ids = idList_;}

  int getNumWeightsPerID() const { return numWeights_;}

  void getWeightsView(const scalar_t *&weights, int &stride, int idx) const
  {
    if (idx < 0 || idx >= numWeights_) {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid vector index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
    }
    size_t length;
    weights_[idx].getStridedList(length, weights, stride);
  }

  ////////////////////////////////////////////////////
  // The VectorAdapter interface.
  ////////////////////////////////////////////////////

  int getNumEntriesPerID() const { return numEntriesPerID_;}

  void getEntriesView(const scalar_t *&entries, int &stride, int idx = 0) const
  {
    if (idx < 0 || idx >= numEntriesPerID_) {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__
           << "  Invalid vector index " << idx << std::endl;
      throw std::runtime_error(emsg.str());
    }
    size_t length;
    entries_[idx].getStridedList(length, entries, stride);
  }

private:

  lno_t numIds_;
  const zgid_t *idList_;

  int numEntriesPerID_;
  ArrayRCP<StridedData<lno_t, scalar_t> > entries_ ;

  int numWeights_;
  ArrayRCP<StridedData<lno_t, scalar_t> > weights_;

  void createBasicVector(
    std::vector<const scalar_t *> &entries,  std::vector<int> &entryStride,
    std::vector<const scalar_t *> &weights, std::vector<int> &weightStrides)
  {
    typedef StridedData<lno_t,scalar_t> input_t;

    if (numIds_){
      int stride = 1;
      entries_ = arcp(new input_t[numEntriesPerID_], 0, numEntriesPerID_, true);
      for (int v=0; v < numEntriesPerID_; v++) {
        if (entryStride.size()) stride = entryStride[v];
        ArrayRCP<const scalar_t> eltV(entries[v], 0, stride*numIds_, false);
        entries_[v] = input_t(eltV, stride);
      }

      if (numWeights_) {
        stride = 1;
        weights_ = arcp(new input_t [numWeights_], 0, numWeights_, true);
        for (int w=0; w < numWeights_; w++){
          if (weightStrides.size()) stride = weightStrides[w];
          ArrayRCP<const scalar_t> wgtV(weights[w], 0, stride*numIds_, false);
          weights_[w] = input_t(wgtV, stride);
        }
      }
    }
  }
};

}  //namespace Zoltan2

#endif
