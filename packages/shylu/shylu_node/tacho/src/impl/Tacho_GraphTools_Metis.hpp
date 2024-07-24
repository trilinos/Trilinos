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
#ifndef __TACHO_GRAPH_TOOLS_METIS_HPP__
#define __TACHO_GRAPH_TOOLS_METIS_HPP__

/// \file Tacho_GraphTools_Scotch.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#if defined(TACHO_HAVE_METIS)
#include "Tacho_Graph.hpp"

#include "trilinos_amd.h"
#include "metis.h"

namespace Tacho {

class GraphTools_Metis {
public:
  typedef typename UseThisDevice<Kokkos::DefaultHostExecutionSpace>::type host_device_type;

  typedef Kokkos::View<idx_t *, host_device_type> idx_t_array;
  typedef Kokkos::View<ordinal_type *, host_device_type> ordinal_type_array;

private:
  // metis main data structure
  idx_t _nvts;
  idx_t_array _xadj, _adjncy, _vwgt;

  idx_t _options[METIS_NOPTIONS];

  // metis output
  idx_t_array _perm_t, _peri_t;
  ordinal_type_array _perm, _peri;

  // status flag
  bool _is_ordered, _verbose;

public:
  GraphTools_Metis();
  GraphTools_Metis(const GraphTools_Metis &b);

  ///
  /// construction of scotch graph
  ///
  GraphTools_Metis(const Graph &g);
  virtual ~GraphTools_Metis();

  ///
  /// setup metis parameters
  ///

  void setVerbose(const bool verbose);
  void setOption(const int id, const idx_t value);

  template <typename ordering_type>
  ordering_type amd_order(ordering_type n, const ordering_type *xadj,
                                           const ordering_type *adjncy,
                          ordering_type *perm,
                          double *control, double *info);

  ///
  /// reorder by metis
  ///

  void reorder(const ordinal_type verbose = 0);

  ordinal_type_array PermVector() const;
  ordinal_type_array InvPermVector() const;

  std::ostream &showMe(std::ostream &os, const bool detail = false) const;
};

} // namespace Tacho
#endif
#endif
