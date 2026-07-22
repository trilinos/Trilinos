// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_PointGenerator_hpp__
#define __Panzer_PointGenerator_hpp__

#include "Kokkos_DynRankView.hpp"

namespace panzer {

/** A class that uses run time polymorphism
  * to generate the reference cell points. This is a utility for
  * PointDescriptor.
  */
class PointGenerator {
public:
  //! Get the points for a particular topology
  virtual Kokkos::DynRankView<double> getPoints(const shards::CellTopology & topo) const = 0;

  //! Get the points for a particular topology
  virtual int numPoints(const shards::CellTopology & topo) const = 0;

  //! Check if the generator can generate points for the given topology
  virtual bool hasPoints(const shards::CellTopology & topo) const = 0;
};

}

#endif
