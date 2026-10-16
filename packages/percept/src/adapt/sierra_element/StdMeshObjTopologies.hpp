// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef adapt_sierra_element_StdMeshObjTopologies_hpp
#define adapt_sierra_element_StdMeshObjTopologies_hpp

#include <Shards_BasicTopologies.hpp>

  namespace percept {
    enum { END_UINT_ARRAY = 99999 };

    namespace Elem {

      namespace StdMeshObjTopologies {

        struct RefinementTopologyExtraEntry
        {
          unsigned ordinal_of_node;               // ordinal of node in the total list of nodes - corresponds to the shards node ordinal
          unsigned rank_of_subcell;               // rank of the subcell this node is associated with                                   
          unsigned ordinal_of_subcell;            // ordinal of the subcell in the shards numbering (e.g. edge # 3)
          unsigned ordinal_of_node_on_subcell;    // ordinal of the node on the subcell (whcih node it is on a subcell that has multiple nodes)
          unsigned num_nodes_on_subcell;          // how many nodes exist on the subcell                                                       
          double parametric_coordinates[3];
        };

        typedef RefinementTopologyExtraEntry RefTopoX[];

        template<class Topo>
        class RefTopoX1
        {
        public:
          typedef RefinementTopologyExtraEntry RefTopoX2[];
        };

        template<class Topo>
        class RefinementTopologyExtra
        {
        public:
          static RefTopoX refinement_topology;
        };

        // Explicit instantiation declarations to suppress -Wundefined-var-template
        // Definitions are provided in StdMeshObjTopologies.cpp via GeneratedRefinementTable.hpp
        template<> RefTopoX RefinementTopologyExtra<shards::Line<2>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Beam<2>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::ShellLine<2>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::ShellLine<3>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Quadrilateral<4>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Triangle<3>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::ShellTriangle<3>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::ShellTriangle<6>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::ShellQuadrilateral<4>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::ShellQuadrilateral<8>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Tetrahedron<4>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Hexahedron<8>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Wedge<6>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Wedge<18>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Wedge<15>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Pyramid<5>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Pyramid<13>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Line<3>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Beam<3>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Triangle<6>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Quadrilateral<8>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Quadrilateral<9>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Hexahedron<27>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Hexahedron<20>>::refinement_topology;
        template<> RefTopoX RefinementTopologyExtra<shards::Tetrahedron<10>>::refinement_topology;

        void bootstrap();

      } // namespace StdMeshObjTopologies
    } // namespace Elem
  } // namespace percept

#endif // adapt_sierra_element_StdMeshObjTopologies_hpp
