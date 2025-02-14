// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_BasicKokkosIdentifierAdapter.hpp
    \brief Defines the BasicKokkosIdentifierAdapter class.
*/

#ifndef _ZOLTAN2_BASICKOKKOSIDENTIFIERADAPTER_HPP_
#define _ZOLTAN2_BASICKOKKOSIDENTIFIERADAPTER_HPP_

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
  class BasicKokkosIdentifierAdapter: public IdentifierAdapter<User> {

public:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
  using scalar_t = typename InputTraits<User>::scalar_t;
  using lno_t    = typename InputTraits<User>::lno_t;
  using gno_t    = typename InputTraits<User>::gno_t;
  using part_t   = typename InputTraits<User>::part_t;
  using node_t   = typename InputTraits<User>::node_t;
  using device_t = typename node_t::device_type;
  using user_t   = User;

  using Base = IdentifierAdapter<User>;
#endif

  /*! \brief Constructor
   *  \param ids should point to a View of identifiers.
   *  \param weights  a list of pointers to arrays of weights.
   *      The number of weights per identifier is assumed to be
   *      \c weights.extent(0).
   *  \param weightStrides  a list of strides for the \c weights.
   *     The weight for weight index \c n for \c ids[k] should be
   *     found at <tt>weights[n][weightStrides[n] * k]</tt>.
   *     If \c weightStrides.size() is zero, it is assumed all strides are one.
   *
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this Adapter.
   */
  BasicKokkosIdentifierAdapter(
    Kokkos::View<gno_t*, device_t> &ids,
    Kokkos::View<scalar_t**, device_t> &weights);

  ////////////////////////////////////////////////////////////////
  // The Adapter interface.
  ////////////////////////////////////////////////////////////////

  size_t getLocalNumIDs() const override {
    return localNumIds_;
  }

  void getIDsView(const gno_t *&ids) const override {
    auto h_ids = Kokkos::create_mirror_view(idsView_);
    Kokkos::deep_copy(h_ids, idsView_);
    ids = h_ids.data();
  }

  void getIDsKokkosView(Kokkos::View<const gno_t *, device_t> &ids) const override {
    ids = idsView_;
  }

  void getIDsDeviceView(
      typename Base::ConstIdsDeviceView &ids) const {
    ids = idsView_;
  }

  void getIDsHostView(
      typename Base::ConstIdsHostView &ids) const {
    auto hostIds = Kokkos::create_mirror_view(idsView_);
    Kokkos::deep_copy(hostIds, idsView_);
    ids = hostIds;
  }

  int getNumWeightsPerID() const override {
    return numWeightsPerID_;
  }

  void getWeightsView(const scalar_t *&wgt, int &stride,
                      int idx = 0) const override {
    auto h_wgts_2d = Kokkos::create_mirror_view(weightsView_);
    Kokkos::deep_copy(h_wgts_2d,weightsView_);

    wgt = Kokkos::subview(h_wgts_2d, Kokkos::ALL, idx).data();
    stride = 1;
  }

  void getWeightsKokkosView(Kokkos::View<scalar_t **, device_t> &wgts) const override {
    wgts = weightsView_;
  }

  void getWeightsDeviceView(typename Base::WeightsDeviceView &wgts) const override {
    wgts = weightsView_;
  }

  void getWeightsHostView(typename Base::WeightsHostView &wgts) const override {
    auto hostWgts = Kokkos::create_mirror_view(weightsView_);
    Kokkos::deep_copy(hostWgts, weightsView_);
    wgts = hostWgts;
  }

private:
  Kokkos::View<gno_t *, device_t> idsView_;
  Kokkos::View<scalar_t **, device_t> weightsView_;
  lno_t localNumIds_ = 0;
  int numWeightsPerID_ = 0;
};

////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////

template <typename User>
BasicKokkosIdentifierAdapter<User>::BasicKokkosIdentifierAdapter(
    Kokkos::View<gno_t *, device_t> &ids,
    Kokkos::View<scalar_t **, device_t> &weights) {
  localNumIds_ = ids.extent(0);
  numWeightsPerID_ = weights.extent(1);
  idsView_ = Kokkos::View<gno_t *, device_t>("idsView_", localNumIds_);
  Kokkos::deep_copy(idsView_, ids);

  weightsView_ = Kokkos::View<scalar_t **, device_t>("weightsView_",
                                                         weights.extent(0),
                                                         weights.extent(1));
  Kokkos::deep_copy(weightsView_, weights);
}

} //namespace Zoltan2

#endif
