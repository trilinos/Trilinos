// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_BasicIdentifierAdapter.hpp
    \brief Defines the BasicIdentifierAdapter class.
*/

#ifndef _ZOLTAN2_BASICIDENTIFIERADAPTER_HPP_
#define _ZOLTAN2_BASICIDENTIFIERADAPTER_HPP_

#include <Kokkos_Core.hpp>
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
class BasicIdentifierAdapter : public IdentifierAdapter<User> {

public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  using scalar_t = typename InputTraits<User>::scalar_t;
  using lno_t = typename InputTraits<User>::lno_t;
  using gno_t = typename InputTraits<User>::gno_t;
  using part_t = typename InputTraits<User>::part_t;
  using node_t = typename InputTraits<User>::node_t;
  using device_t = typename node_t::device_type;
  using user_t = User;

  using Base = IdentifierAdapter<User>;
#endif

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
                         std::vector<const scalar_t *> &weights,
                         std::vector<int> &weightStrides);

  /*! \brief Constructor
   *  \param numIds is the number of identifiers in the list
   *  \param ids should point to a list of numIds identifiers.
   *
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this Adapter.
   */
  BasicIdentifierAdapter(lno_t numIds, const gno_t *idPtr)
      : localNumIDs_(numIds), idList_(idPtr), weights_() {}

  /*! \brief Constructor
   *  \param ids should point to a View of identifiers.
   *  \param weights  a list of pointers to arrays of weights.
   *      The number of weights per identifier is assumed to be
   *      \c weights.extent(1).
   *  \param weightStrides  a list of strides for the \c weights.
   *     The weight for weight index \c n for \c ids[k] should be
   *     found at <tt>weights[n][weightStrides[n] * k]</tt>.
   *     If \c weightStrides.size() is zero, it is assumed all strides are one.
   *
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this Adapter.
   */
  BasicIdentifierAdapter(typename Base::IdsDeviceView &ids,
                         typename Base::WeightsDeviceView &weights);

  ////////////////////////////////////////////////////////////////
  // The Adapter interface.
  ////////////////////////////////////////////////////////////////

  size_t getLocalNumIDs() const { return localNumIDs_; }

  void getIDsView(const gno_t *&ids) const { ids = idList_; }

  void getIDsKokkosView(typename Base::ConstIdsDeviceView &ids) const override {
    ids = idsView_;
  }

  void getIDsDeviceView(typename Base::ConstIdsDeviceView &ids) const {
    ids = idsView_;
  }

  void getIDsHostView(typename Base::ConstIdsHostView &ids) const {
    auto hostIds = Kokkos::create_mirror_view(idsView_);
    Kokkos::deep_copy(hostIds, idsView_);
    ids = hostIds;
  }

  int getNumWeightsPerID() const { return numWeightsPerID_; }

  void getWeightsView(const scalar_t *&wgt, int &stride, int idx = 0) const {
    if (idx < 0 || idx >= weights_.size()) {
      std::ostringstream emsg;
      emsg << __FILE__ << ":" << __LINE__ << "  Invalid weight index " << idx
           << std::endl;
      throw std::runtime_error(emsg.str());
    }
    size_t length;
    weights_[idx].getStridedList(length, wgt, stride);
  }

  void
  getWeightsKokkosView(typename Base::WeightsDeviceView &wgts) const override {
    wgts = weightsView_;
  }

  void getWeightsDeviceView(typename Base::WeightsDeviceView1D &deviceWgts,
                            int idx = 0) const override {
    AssertCondition((idx >= 0) and (idx < numWeightsPerID_),
                    "Invalid weight index.");

    const auto size = weightsView_.extent(0);
    deviceWgts = typename Base::WeightsDeviceView1D("deviceWgts", size);

    Kokkos::parallel_for(
        size, KOKKOS_CLASS_LAMBDA(const int id) {
          deviceWgts(id) = weightsView_(id, idx);
        });

    Kokkos::fence();
  }

  void
  getWeightsDeviceView(typename Base::WeightsDeviceView &wgts) const override {
    wgts = weightsView_;
  }

  void getWeightsHostView(typename Base::WeightsHostView1D &wgts,
                          int idx = 0) const override {
    AssertCondition((idx >= 0) and (idx < numWeightsPerID_),
                    "Invalid weight index.");

    auto weightsDevice =
        typename Base::WeightsDeviceView1D("weights", weightsView_.extent(0));
    getWeightsDeviceView(weightsDevice, idx);

    wgts = Kokkos::create_mirror_view(weightsDevice);
    Kokkos::deep_copy(wgts, weightsDevice);
  }

  void getWeightsHostView(typename Base::WeightsHostView &wgts) const override {
    auto hostWgts = Kokkos::create_mirror_view(weightsView_);
    Kokkos::deep_copy(hostWgts, weightsView_);
    wgts = hostWgts;
  }

private:
  lno_t localNumIDs_ = 0;
  const gno_t *idList_;
  ArrayRCP<StridedData<lno_t, scalar_t>> weights_;
  int numWeightsPerID_ = 0;

  typename Base::IdsDeviceView idsView_;
  typename Base::WeightsDeviceView weightsView_;
};

////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////

template <typename User>
BasicIdentifierAdapter<User>::BasicIdentifierAdapter(
    lno_t numIds, const gno_t *idPtr, std::vector<const scalar_t *> &weights,
    std::vector<int> &weightStrides)
    : localNumIDs_(numIds), idList_(idPtr), weights_() {
  typedef StridedData<lno_t, scalar_t> input_t;
  numWeightsPerID_ = weights.size();

  if (numWeightsPerID_ > 0) {
    weights_ = arcp(new input_t[numWeightsPerID_], 0, numWeightsPerID_, true);

    if (numIds > 0) {
      for (int i = 0; i < numWeightsPerID_; i++) {
        int stride = weightStrides.size() ? weightStrides[i] : 1;
        ArrayRCP<const scalar_t> wgtV(weights[i], 0, stride * numIds, false);
        weights_[i] = input_t(wgtV, stride);
      }
    }
  }
}

template <typename User>
BasicIdentifierAdapter<User>::BasicIdentifierAdapter(
    typename Base::IdsDeviceView &ids,
    typename Base::WeightsDeviceView &weights) {
  idsView_ = typename Base::IdsDeviceView("idsView_", ids.extent(0));
  Kokkos::deep_copy(idsView_, ids);

  weightsView_ = typename Base::WeightsDeviceView(
      "weightsView_", weights.extent(0), weights.extent(1));
  Kokkos::deep_copy(weightsView_, weights);
  localNumIDs_ = idsView_.extent(0);
  numWeightsPerID_ = weights.extent(1);
}

} // namespace Zoltan2

#endif
