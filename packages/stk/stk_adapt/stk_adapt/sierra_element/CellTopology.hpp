/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef stk_adapt_sierra_element_CellTopology_hpp
#define stk_adapt_sierra_element_CellTopology_hpp

//#include <SIERRA_code_types.h>

#include <Shards_CellTopology.hpp>

#include <stk_adapt/sierra_element/stk_percept_code_types.hpp>

namespace stk { 
  namespace adapt { 
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

      CellTopology getBasicCellTopology(TopologyId id);
      CellTopology getBasicCellTopology(const char *name);

      CellTopology nodeCellTopology(const CellTopology &cell_topology, UInt ordinal);
      CellTopology edgeCellTopology(const CellTopology &cell_topology, UInt ordinal);
      CellTopology faceCellTopology(const CellTopology &cell_topology, UInt ordinal);
      CellTopology sideCellTopology(const CellTopology &cell_topology, UInt ordinal);

      bool isElement(const CellTopology &cell_topology, unsigned spatial_dimension);
      bool isSolidElement(const CellTopology &cell_topology, unsigned spatial_dimension);
      bool isShellElement(const CellTopology &cell_topology);
      bool isRodElement(const CellTopology &cell_topology);
      bool isParticleElement(const CellTopology &cell_topology);

      int getEdgeNode(const CellTopology &cell_topology, unsigned edge, unsigned node_of_edge);
      int getFaceNode(const CellTopology &cell_topology, unsigned face, unsigned node_of_face);
      int getSideNode(const CellTopology &cell_topology, unsigned side, unsigned node_of_side);
      int getFaceEdge(const CellTopology &cell_topology, unsigned face, unsigned edge_of_face);

      const unsigned *getNodesOfEdge(const CellTopology &cell_topology, unsigned edge);
      const unsigned *getNodesOfFace(const CellTopology &cell_topology, unsigned face);
      const unsigned *getNodesOfSide(const CellTopology &cell_topology, unsigned side);

      bool isCellTopologySubsetOf(const CellTopology &cell_topology, const CellTopology &richer);

      unsigned getParametricDimension(const CellTopology &cell_topology) ;

      int findReversePermutation(const CellTopologyData &top, int permutation_ord);

    } // namespace Elem
  } // namespace adapt
} // namespace stk

// namespace shards {

//   inline bool operator<(const CellTopology &c1, const CellTopology &c2) {
//     return c1.getTopology() < c2.getTopology();
//   }

//   inline bool operator==(const CellTopology &c1, const CellTopology &c2) {
//     return c1.getTopology() == c2.getTopology();
//   }

//   inline bool operator!=(const CellTopology &c1, const CellTopology &c2) {
//     return !(c1 == c2);
//   }

// } // namespace shards

#endif // stk_adapt_sierra_element_CellTopology_hpp

