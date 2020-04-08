// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef adapt_sierra_element_CellTopology_hpp
#define adapt_sierra_element_CellTopology_hpp

//#include <SIERRA_code_types.h>

#include <Shards_CellTopology.hpp>

#include <adapt/sierra_element/percept_code_types.hpp>

  namespace percept { 
    namespace Elem {

      enum TopologyId {
        INVALID       = 0,
        NODE_0,

        EDGE_2,
        EDGE_3,

        FACE_TRI_3,
        FACE_TRI_4,
        FACE_TRI_6,
        FACE_QUAD_4,
        FACE_QUAD_5,
        FACE_QUAD_8,
        FACE_QUAD_9,

        PARTICLE_1,
        ROD_2,
        ROD_3,
        SHELL_LINE_2,
        SHELL_LINE_3,
        SHELL_TRI_3,
        //  SHELL_TRI_4,
        SHELL_TRI_6,
        SHELL_QUAD_4,
        SHELL_QUAD_9,

        SOLID_HEX_8,
        SOLID_HEX_20,
        SOLID_HEX_27,
        SOLID_TET_4,
        SOLID_TET_8,
        SOLID_TET_10,
        SOLID_WEDGE_6,
        SOLID_WEDGE_15,
        SOLID_PYRAMID_5,
        SOLID_PYRAMID_13,

        SHELL_QUAD_8,
        SOLID_WEDGE_18,
        SOLID_PYRAMID_14,

        FACE_PENT_5,
        FACE_HEX_6
      };

      using namespace shards;

      typedef shards::CellTopology CellTopology;

      template <class T>
      CellTopology getCellTopology() {
        return CellTopology(shards::getCellTopologyData<T>());
      }

      TopologyId getCellTopologyId(const CellTopology &cell_topology);

      CellTopology getBasicCellTopology(const char *name);

      CellTopology edgeCellTopology(const CellTopology &cell_topology, UInt ordinal);
      CellTopology faceCellTopology(const CellTopology &cell_topology, UInt ordinal);
      CellTopology sideCellTopology(const CellTopology &cell_topology, UInt ordinal);

      int getFaceEdge(const CellTopology &cell_topology, unsigned face, unsigned edge_of_face);

    } // namespace Elem
  } // namespace percept

#endif // adapt_sierra_element_CellTopology_hpp

