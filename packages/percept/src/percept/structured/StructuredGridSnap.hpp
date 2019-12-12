// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_StructuredGridSnap_hpp
#define percept_StructuredGridSnap_hpp

#include <percept/structured/BlockStructuredGrid.hpp>
#if !STK_PERCEPT_LITE
#  if defined(STK_BUILT_IN_SIERRA)
#    include <cgns/Iocgns_DatabaseIO.h>
#  else
#    include <Iocgns_DatabaseIO.h>
#  endif
#endif

class PGeom;
class PGeomAssocStructured;

namespace percept {

  class StructuredGridSnap {

  public:
    std::shared_ptr<BlockStructuredGrid> m_grid;
    int m_debug;

    StructuredGridSnap(std::shared_ptr<BlockStructuredGrid> input, int debug=0) :
      m_debug(debug)
    {
      m_grid = input;
      m_zero_distance_tol = 1e-6;
    }

    void snap_to_geometry(std::shared_ptr<PGeom> pgeom);

    void snap_to_geometry(std::shared_ptr<PGeom> pgeom,
                          PGeomAssocStructured &passoc);

  private:
    void snap_to_surface(std::shared_ptr<PGeom> pgeom,
                         std::shared_ptr<StructuredBlock> sgi,
                         const std::array<unsigned, 3> &range_min,
                         const std::array<unsigned, 3> &range_max,
                         const int surf_id);

    void snap_to_curve(std::shared_ptr<PGeom> pgeom,
                       std::shared_ptr<StructuredBlock> sgi,
                       const std::array<unsigned, 3> &range_min,
                       const std::array<unsigned, 3> &range_max,
                       const int curve_id);

  private:
    double m_zero_distance_tol; // distance tolerance for comparing two points

  };


}


#endif
