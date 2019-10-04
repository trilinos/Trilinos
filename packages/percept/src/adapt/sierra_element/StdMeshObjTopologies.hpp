// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef adapt_sierra_element_StdMeshObjTopologies_hpp
#define adapt_sierra_element_StdMeshObjTopologies_hpp

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

        void bootstrap();

      } // namespace StdMeshObjTopologies
    } // namespace Elem
  } // namespace percept

#endif // adapt_sierra_element_StdMeshObjTopologies_hpp
