// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_GRAPH_TOOLS_HPP__
#define __TACHO_GRAPH_TOOLS_HPP__

/// \file Tacho_GraphTools.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Graph.hpp"
#include "Tacho_Util.hpp"

namespace Tacho {

class GraphTools {
public:
  typedef typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type host_device_type;
  typedef Kokkos::View<ordinal_type *, host_device_type> ordinal_type_array;

private:
  ordinal_type_array _perm, _peri;

  // status flag
  bool _is_ordered, _verbose;

public:
  GraphTools() = default;
  GraphTools(const GraphTools &b) = default;
  ~GraphTools() = default;

  ///
  /// construction of scotch graph
  ///
  GraphTools(const Graph &g) {
    _is_ordered = false;
    _verbose = false;

    // input
    const ordinal_type m = g.NumRows();

    // output
    _perm = ordinal_type_array(do_not_initialize_tag("PermutationArray"), m);
    _peri = ordinal_type_array(do_not_initialize_tag("InvPermutationArray"), m);
  }

  ///
  /// setup parameters
  ///

  void setVerbose(const bool verbose) { _verbose = verbose; }

  ///
  /// reorder
  ///
  void reorder(const ordinal_type verbose = 0) {
    // do nothing
    // later implement AMD
    const ordinal_type m = _perm.extent(0);
    for (ordinal_type i = 0; i < m; ++i) {
      _perm(i) = i;
      _peri(i) = i;
    }
  }

  ordinal_type_array PermVector() const { return _perm; }
  ordinal_type_array InvPermVector() const { return _peri; }

  std::ostream &showMe(std::ostream &os, const bool detail = false) const {
    std::streamsize prec = os.precision();
    os.precision(4);
    os << std::scientific;

    if (_is_ordered)
      os << " -- Native Ordering -- " << std::endl << "  PERM     PERI" << std::endl;
    else
      os << " -- Not Ordered -- " << std::endl;

    if (detail) {
      const ordinal_type w = 6, m = _perm.extent(0);
      for (ordinal_type i = 0; i < m; ++i)
        os << std::setw(w) << _perm[i] << "   " << std::setw(w) << _peri[i] << "   " << std::endl;
    }
    os.unsetf(std::ios::scientific);
    os.precision(prec);

    return os;
  }
};

} // namespace Tacho
#endif
