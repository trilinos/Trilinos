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

using Kokkos::ALL;

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
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::part_t   part_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef typename node_t::device_type         device_t;
  typedef User user_t;

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
    return idsView_.extent(0);
  }

  void getIDsView(const gno_t *&ids) const override {
    auto kokkosIds = idsView_.view_host();
    ids = kokkosIds.data();
  }

  void getIDsKokkosView(Kokkos::View<const gno_t *, device_t> &ids) const override {
    ids = idsView_.view_device();
  }

  int getNumWeightsPerID() const override {
    return weightsView_.extent(1);
  }

  void getWeightsView(const scalar_t *&wgt, int &stride,
                      int idx = 0) const override
  {
    auto h_wgts_2d = weightsView_.view_host();

    wgt = Kokkos::subview(h_wgts_2d, Kokkos::ALL, idx).data();
    stride = 1;
  }

  void getWeightsKokkosView(Kokkos::View<scalar_t **, device_t> &wgts) const override {
    wgts = weightsView_.template view<device_t>();
  }

private:
  Kokkos::DualView<gno_t *, device_t> idsView_;
  Kokkos::DualView<scalar_t **, device_t> weightsView_;
};

////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////

template <typename User>
BasicKokkosIdentifierAdapter<User>::BasicKokkosIdentifierAdapter(
    Kokkos::View<gno_t *, device_t> &ids,
    Kokkos::View<scalar_t **, device_t> &weights)
{
  idsView_ = Kokkos::DualView<gno_t *, device_t>("idsView_", ids.extent(0));
  Kokkos::deep_copy(idsView_.h_view, ids);

  weightsView_ = Kokkos::DualView<scalar_t **, device_t>("weightsView_",
                                                         weights.extent(0),
                                                         weights.extent(1));
  Kokkos::deep_copy(weightsView_.h_view, weights);

  weightsView_.modify_host();
  weightsView_.sync_host();
  weightsView_.template sync<device_t>();

  idsView_.modify_host();
  idsView_.sync_host();
  idsView_.template sync<device_t>();
}

} //namespace Zoltan2

#endif
