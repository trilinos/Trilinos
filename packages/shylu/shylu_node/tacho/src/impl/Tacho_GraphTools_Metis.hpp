// clang-format off
/* =====================================================================================
Copyright 2022 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

SCR#:2790.0

This file is part of Tacho. Tacho is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Kyungjoo Kim at <kyukim@sandia.gov,https://github.com/kyungjoo-kim>

Sandia National Laboratories, Albuquerque, NM, USA
===================================================================================== */
// clang-format on
#ifndef __TACHO_GRAPH_TOOLS_METIS_HPP__
#define __TACHO_GRAPH_TOOLS_METIS_HPP__

/// \file Tacho_GraphTools_Scotch.hpp
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

#if defined(TACHO_HAVE_METIS)
#include "Tacho_Graph.hpp"
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
