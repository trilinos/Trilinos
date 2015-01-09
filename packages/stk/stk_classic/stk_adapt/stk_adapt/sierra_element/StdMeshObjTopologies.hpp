/*--------------------------------------------------------------------*/
/*    Copyright 1997 - 2009 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef stk_adapt_sierra_element_StdMeshObjTopologies_hpp
#define stk_adapt_sierra_element_StdMeshObjTopologies_hpp

namespace stk_classic {
  namespace adapt {
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
  } // namespace adapt
} // namespace stk_classic

#endif // stk_adapt_sierra_element_StdMeshObjTopologies_hpp
