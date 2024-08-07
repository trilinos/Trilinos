// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_BasicVectorAdapter.hpp
    \brief Defines the BasicVectorAdapter class.
*/

#ifndef _ZOLTAN2_BASICVECTORADAPTER_HPP_
#define _ZOLTAN2_BASICVECTORADAPTER_HPP_

#include "Zoltan2_Adapter.hpp"
#include <Zoltan2_StridedData.hpp>
#include <Zoltan2_VectorAdapter.hpp>

namespace Zoltan2 {

/*!  \brief BasicVectorAdapter represents a vector (plus optional weights)
            supplied by the user as pointers to strided arrays.

    BasicVectorAdapter may be a single vector or multivector (set of
    corresponding vectors with the same global identifiers and
    distribution across processes).  A constructor specifically for use
    of BasicVectorAdapter to represent geometric coordinates is also provided.

    Data types:
    \li \c scalar_t is the data type for weights and vector entry values.
    \li \c lno_t is the integral data type used by Zoltan2 for local indices
           and local counts.
    \li \c gno_t is the data type used by the application for global Ids; must
           be a Teuchos Ordinal.  (Teuchos Ordinals are those data types for
          which traits are defined in
          Trilinos/packages/teuchos/src/Teuchos_OrdinalTraits.hpp.)
    \li \c node_t is a Kokkos CPU Node type.  If you don't use Kokkos, you can
          ignore this data type.

    The template parameter (\c User) is a C++ class type which provides the
    actual data types with which the Zoltan2 library will be compiled, through
    a Traits mechanism.  \c User may be the
    actual class used by application to represent a vector, or it may be
    the empty helper class \c BasicUserTypes with which a Zoltan2 user
    can easily supply the data types for the library.
*/

template <typename User> class BasicVectorAdapter : public VectorAdapter<User> {

public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS

  using scalar_t = typename InputTraits<User>::scalar_t;
  using lno_t = typename InputTraits<User>::lno_t;
  using gno_t = typename InputTraits<User>::gno_t;
  using offset_t = typename InputTraits<User>::offset_t;
  using part_t = typename InputTraits<User>::part_t;
  using node_t = typename InputTraits<User>::node_t;
  using user_t = User;
  using Base = BaseAdapter<User>;

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

  BasicVectorAdapter(lno_t numIds, const gno_t *ids, const scalar_t *entries,
                     int entryStride = 1, bool usewgts = false,
                     const scalar_t *wgts = NULL, int wgtStride = 1)
      : numIds_(numIds), idList_(ids), numEntriesPerID_(1),
        numWeights_(usewgts == true) {
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

  BasicVectorAdapter(lno_t numIds, const gno_t *ids,
                     std::vector<const scalar_t *> &entries,
                     std::vector<int> &entryStride,
                     std::vector<const scalar_t *> &weights,
                     std::vector<int> &weightStrides)
      : numIds_(numIds), idList_(ids), numEntriesPerID_(entries.size()),
        numWeights_(weights.size()) {
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

  BasicVectorAdapter(lno_t numIds, const gno_t *ids, const scalar_t *x,
                     const scalar_t *y, const scalar_t *z, int xStride = 1,
                     int yStride = 1, int zStride = 1, bool usewgts = false,
                     const scalar_t *wgts = NULL, int wgtStride = 1)
      : numIds_(numIds), idList_(ids), numEntriesPerID_(0),
        numWeights_(usewgts == true) {
    std::vector<const scalar_t *> values, weightValues;
    std::vector<int> strides, weightStrides;

    if (x) {
      values.push_back(x);
      strides.push_back(xStride);
      numEntriesPerID_++;
      if (y) {
        values.push_back(y);
        strides.push_back(yStride);
        numEntriesPerID_++;
        if (z) {
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

  ////////////////////////////////////////////////////////////////
  // The Adapter interface.
  ////////////////////////////////////////////////////////////////

  size_t getLocalNumIDs() const { return numIds_; }

  void getIDsView(const gno_t *&ids) const { ids = idList_; }

  void getIDsKokkosView(typename Base::ConstIdsDeviceView &ids) const {
    ids = this->kIds_;
  }

  void getIDsHostView(typename Base::ConstIdsHostView &ids) const override {
    auto hostIds = Kokkos::create_mirror_view(this->kIds_);
    Kokkos::deep_copy(hostIds, this->kIds_);
    ids = hostIds;
  }

  void getIDsDeviceView(typename Base::ConstIdsDeviceView &ids) const override {
    ids = this->kIds_;
  }

  int getNumWeightsPerID() const { return numWeights_; }

  virtual void
  getWeightsKokkos2dView(typename Base::WeightsDeviceView &wgt) const {
    wgt = kWeights_;
  }

  void getWeightsView(const scalar_t *&weights, int &stride, int idx) const {
    AssertCondition((idx >= 0) and (idx < numWeights_),
                    "Invalid weight index.");

    size_t length;
    weights_[idx].getStridedList(length, weights, stride);
  }

  void getWeightsHostView(typename Base::WeightsHostView1D &hostWgts,
                          int idx = 0) const override {
    AssertCondition((idx >= 0) and (idx < numWeights_),
                    "Invalid weight index.");

    auto weightsDevice =
        typename Base::WeightsDeviceView1D("weights", kWeights_.extent(0));
    getWeightsDeviceView(weightsDevice, idx);

    hostWgts = Kokkos::create_mirror_view(weightsDevice);
    Kokkos::deep_copy(hostWgts, weightsDevice);
  }

  void getWeightsHostView(typename Base::WeightsHostView &wgts) const override {
    auto hostWeights = Kokkos::create_mirror_view(kWeights_);
    Kokkos::deep_copy(hostWeights, kWeights_);
    wgts = hostWeights;
  }

  void getWeightsDeviceView(typename Base::WeightsDeviceView1D &deviceWgts,
                            int idx = 0) const override {
    AssertCondition((idx >= 0) and (idx < numWeights_),
                    "Invalid weight index.");

    const auto size = kWeights_.extent(0);
    deviceWgts = typename Base::WeightsDeviceView1D("weights", size);

    Kokkos::parallel_for(
        size, KOKKOS_CLASS_LAMBDA(const int id) {
          deviceWgts(id) = kWeights_(id, idx);
        });

    Kokkos::fence();
  }

  void
  getWeightsDeviceView(typename Base::WeightsDeviceView &wgts) const override {
    wgts = kWeights_;
  }

  ////////////////////////////////////////////////////
  // The VectorAdapter interface.
  ////////////////////////////////////////////////////

  int getNumEntriesPerID() const { return numEntriesPerID_; }

  void getEntriesView(const scalar_t *&entries, int &stride,
                      int idx = 0) const {
    if (idx < 0 || idx >= numEntriesPerID_) {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__ << "  Invalid vector index " << idx
           << std::endl;
      throw std::runtime_error(emsg.str());
    }
    size_t length;
    entries_[idx].getStridedList(length, entries, stride);
  }

  void getEntriesKokkosView(
      typename AdapterWithCoords<User>::CoordsDeviceView &entries) const {
    entries = kEntries_;
  }

  void getEntriesHostView(typename AdapterWithCoords<User>::CoordsHostView
                              &entries) const override {
    auto hostEntries = Kokkos::create_mirror_view(kEntries_);
    Kokkos::deep_copy(hostEntries, kEntries_);
    entries = hostEntries;
  }

  void getEntriesDeviceView(typename AdapterWithCoords<User>::CoordsDeviceView
                                &entries) const override {
    entries = kEntries_;
  }

private:
  lno_t numIds_;
  const gno_t *idList_;

  int numEntriesPerID_;
  int numWeights_;

  // Old API variable members
  ArrayRCP<StridedData<lno_t, scalar_t>> entries_;
  ArrayRCP<StridedData<lno_t, scalar_t>> weights_;

  // New API variable members
  typename Base::IdsDeviceView kIds_;
  typename AdapterWithCoords<User>::CoordsDeviceView kEntries_;
  typename Base::WeightsDeviceView kWeights_;

  void createBasicVector(std::vector<const scalar_t *> &entries,
                         std::vector<int> &entryStride,
                         std::vector<const scalar_t *> &weights,
                         std::vector<int> &weightStrides) {
    typedef StridedData<lno_t, scalar_t> input_t;

    if (numIds_) {
      // make kokkos ids
      kIds_ = typename Base::IdsDeviceView(
          Kokkos::ViewAllocateWithoutInitializing("ids"), numIds_);
      auto host_kIds_ = Kokkos::create_mirror_view(kIds_);
      for (int n = 0; n < numIds_; ++n) {
        host_kIds_(n) = idList_[n];
      }
      Kokkos::deep_copy(kIds_, host_kIds_);

      // make coordinates
      int stride = 1;
      entries_ = arcp(new input_t[numEntriesPerID_], 0, numEntriesPerID_, true);
      for (int v = 0; v < numEntriesPerID_; v++) {
        if (entryStride.size())
          stride = entryStride[v];
        ArrayRCP<const scalar_t> eltV(entries[v], 0, stride * numIds_, false);
        entries_[v] = input_t(eltV, stride);
      }

      // setup kokkos entries
      kEntries_ = typename AdapterWithCoords<User>::CoordsDeviceView(
          Kokkos::ViewAllocateWithoutInitializing("entries"), numIds_,
          numEntriesPerID_);

      size_t length;
      const scalar_t *entriesPtr;

      auto host_kokkos_entries = Kokkos::create_mirror_view(this->kEntries_);

      for (int idx = 0; idx < numEntriesPerID_; idx++) {
        entries_[idx].getStridedList(length, entriesPtr, stride);
        size_t fill_index = 0;
        for (int n = 0; n < numIds_; ++n) {
          host_kokkos_entries(fill_index++, idx) = entriesPtr[n * stride];
        }
      }
      Kokkos::deep_copy(this->kEntries_, host_kokkos_entries);
    }

    if (numWeights_) {
      int stride = 1;
      weights_ = arcp(new input_t[numWeights_], 0, numWeights_, true);
      for (int w = 0; w < numWeights_; w++) {
        if (weightStrides.size())
          stride = weightStrides[w];
        ArrayRCP<const scalar_t> wgtV(weights[w], 0, stride * numIds_, false);
        weights_[w] = input_t(wgtV, stride);
      }

      // set up final view with weights
      kWeights_ = typename Base::WeightsDeviceView(
          Kokkos::ViewAllocateWithoutInitializing("kokkos weights"), numIds_,
          numWeights_);

      // setup weights
      auto host_weight_temp_values =
          Kokkos::create_mirror_view(this->kWeights_);
      for (int idx = 0; idx < numWeights_; ++idx) {
        const scalar_t *weightsPtr;
        size_t length;
        weights_[idx].getStridedList(length, weightsPtr, stride);
        size_t fill_index = 0;
        for (size_t n = 0; n < length; n += stride) {
          host_weight_temp_values(fill_index++, idx) = weightsPtr[n];
        }
      }
      Kokkos::deep_copy(this->kWeights_, host_weight_temp_values);
    }
  }
};

} // namespace Zoltan2

#endif
