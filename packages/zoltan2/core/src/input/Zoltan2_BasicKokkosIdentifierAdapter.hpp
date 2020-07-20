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
    Kokkos::View<gno_t*, typename node_t::device_type> &ids,
    Kokkos::View<scalar_t**, typename node_t::device_type> &weights);

  ////////////////////////////////////////////////////////////////
  // The Adapter interface.
  ////////////////////////////////////////////////////////////////

  size_t getLocalNumIDs() const {
    return idsView_.extent(0);
  }

  void getIDsKokkosView(Kokkos::View<const gno_t *,
    typename node_t::device_type> &ids) const {
    ids = idsView_;
  }

  int getNumWeightsPerID() const {
    return weightsView_.extent(1);
  }

  void getWeightsKokkosView(Kokkos::View<scalar_t **,
    typename node_t::device_type> &wgts) const {
    wgts = weightsView_;
  }

private:
  Kokkos::View<const gno_t *, typename node_t::device_type> idsView_;
  Kokkos::View<scalar_t **, typename node_t::device_type> weightsView_;
};

////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////

template <typename User>
  BasicKokkosIdentifierAdapter<User>::BasicKokkosIdentifierAdapter(
    Kokkos::View<gno_t*, typename node_t::device_type> &ids,
    Kokkos::View<scalar_t**, typename node_t::device_type> &weights):
      idsView_(ids), weightsView_(weights) {
}
  
} //namespace Zoltan2
  
#endif
