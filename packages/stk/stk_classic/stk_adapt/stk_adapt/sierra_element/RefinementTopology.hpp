/**-------------------------------------------------------------------*
 *    Copyright 1999 - 2009 Sandia Corporation.                       *
 *    Under the terms of Contract DE-AC04-94AL85000, there is a       *
 *    non-exclusive license for use of this work by or on behalf      *
 *    of the U.S. Government.  Export of this program may require     *
 *    a license from the United States Government.                    *
 *--------------------------------------------------------------------*/

#ifndef stk_adapt_sierra_element_RefinementTopology_hpp
#define stk_adapt_sierra_element_RefinementTopology_hpp

//#include <stk_percept_code_types.h>

#include <utility>

#include <stk_adapt/sierra_element/RefinementKey.hpp>
#include <stk_adapt/sierra_element/CellTopology.hpp>

namespace stk_classic { namespace adapt {
  namespace Elem {

    class MeshObjTopology;
    class RefinementKey;
    class MeshObjRefinementTopology;

    class RefinementTopology
    {
    public:
      RefinementTopology(MeshObjTopology *mesh_obj_topology,
                         UInt num_child, const MeshObjTopology * const *child_topology,
                         UInt num_child_nodes, const UInt * const *child_nodes,
                         const UInt num_edges, const UInt * const *edge_node,
                         const UInt num_faces, const UInt * const *face_node,
                         UInt num_orientations, const UInt * const *perm_node, const UInt * const *perm_edge, 
                         bool homogeneous_child);

      RefinementTopology(const CellTopology &mesh_obj_topology,
                         UInt num_child, const MeshObjTopology * const *child_topology,
                         UInt num_child_nodes, const UInt * const *child_nodes,
                         const UInt num_edges, const UInt * const *edge_node,
                         const UInt num_faces, const UInt * const *face_node,
                         UInt num_orientations, const UInt * const *perm_node, const UInt * const *perm_edge, 
                         bool homogeneous_child);

      ~RefinementTopology();

    private:
      RefinementTopology(const RefinementTopology &);
      RefinementTopology &operator=(const RefinementTopology &);

    public:
      /** Number of refinement child topologies. */
      UInt num_child() const {
        return m_numChild;
      }

      /** Query if the refined mesh object topologies are homogeneous. */
      bool homogeneous_child() const {
        return m_homogeneousChild;
      }

      const CellTopology *child_cell_topology() const {
        return m_childCellTopology;
      }
  
      CellTopology child_cell_topology(UInt ordinal) const {
        return (ordinal < m_numChild) ? m_childCellTopology[ordinal] : CellTopology();
      }

      /** Total number of unique nodes of the connected child objects.
       *  Nodes that are shared by child objects are only counted once.
       *  Every node of the parent is also one of the child nodes.
       */
      UInt num_child_nodes() const {
        return m_numChildNodes;
      }

      const UInt * const * child_nodes() const {
        return m_childNode;
      }

      /** Map ordinals of child object topology nodes to child nodes.
       *  Array dimension is [ child_topology(child)->num_nodes() ]
       *  @pre child < num_child()
       *  @pre node_of_child < child_topology(child)->num_nodes()
       *  @post 0 <= child_node(child)[i] < num_child_nodes
       */
      const UInt * child_node(UInt child) const {
        return (child < m_numChild) ? m_childNode[child] : NULL;
      }

      const UInt *edge_node(UInt edge) const {
        if (edge >= m_numEdges)
          return NULL;
    
        return m_edgeNode[edge];
      }
  
      const UInt *face_node(UInt face) const {
        if (face >= m_numFaces)
          return NULL;
    
        return m_faceNode[face];
      }
  
      const UInt *node_permutation(UInt orientation) const;

      /** Permutation vector from the expected orientation
       *  to the input permutation.
       */
      const UInt *edge_permutation(UInt orientation) const;

      bool query_refinement_topology(RefinementKey &object_key, MeshObjRefinementTopology & refTop) const;

      bool refine_rivara_tri(RefinementKey &, MeshObjRefinementTopology & refTop) const;
      bool refine_rivara_tet(RefinementKey &, MeshObjRefinementTopology & refTop) const;

      /*------------------------------------------------------------------*/
      /** Mapping of parent->face->child to parent->child->face
       * @pre  face_ordinal < num_faces()
       * @pre  face_child_ordinal < face_topology(face_ordinal)->num_child()
       * @return (child_ordinal, child_face_ordinal)
       */
      std::pair<UInt,UInt> child_face(const UInt face_ordinal ,
                                      const UInt face_child_ordinal) const;

      /** Mapping of parent->edge->child to parent->child->edge
       * @pre  edge_ordinal < getEdgeCount()
       * @pre  edge_child_ordinal < edge_topology(edge_ordinal)->num_child()
       * @return (child_ordinal, child_edge_ordinal)
       */
      std::pair<UInt,UInt> child_edge(const UInt edge_ordinal ,
                                      const UInt edge_child_ordinal) const;

    private:
      CellTopology                          m_cellTopology;
      UInt                                  m_numChild;
      CellTopology *                        m_childCellTopology;
      UInt                                  m_numChildNodes;
      const UInt * const *                  m_childNode;
      const UInt                            m_numEdges;  
      const UInt * const *                  m_edgeNode;             ///< Node ordinals map to edges
      ///< m_edgeNode[ m_numEdges ][ num-nodes-of-edge ]

      const UInt                            m_numFaces;
      const UInt * const *                  m_faceNode;             ///< Mapping of node ordinals to faces
      ///< m_faceNode[ m_numFaces ][ num-nodes-of-face ]
   
      UInt                                  m_numOrientations;
      const UInt * const *                  m_nodePermutation;
      const UInt * const *                  m_edgePermutation;

      bool                                  m_homogeneousChild;
    };


    /**
     * Class to hold refinement information for partial and heterogeneous
     */
    class MeshObjRefinementTopology
    {
    public:
      MeshObjRefinementTopology();
      ~MeshObjRefinementTopology();

      UInt num_child() const;

      UInt num_child_nodes() const ;

      CellTopology child_cell_topology(UInt child) const ;

      const UInt * child_node(UInt child) const ;

      bool homogeneous_refinement() const;

      bool full_refinement() const;
      /*------------------------------------------------------------------*/
      /** Mapping of parent->face->child to parent->child->face
       * @pre  face_ordinal < num_faces()
       * @pre  face_child_ordinal < face_topology(face_ordinal)->num_child()
       * @arg  objTop need to have objects regular topology as well
       * @return (child_ordinal, child_face_ordinal)
       */
      std::pair<UInt,UInt> child_face(const UInt face_ordinal ,
                                      const UInt face_child_ordinal,
                                      const Elem::CellTopology & objTop ,
                                      const RefinementKey &objDesiredKey) const ;

      /** Mapping of parent->edge->child to parent->child->edge
       * @pre  edge_ordinal < getEdgeCount()
       * @pre  edge_child_ordinal < edge_topology(edge_ordinal)->num_child()
       * @arg  objTop need to have objects regular topology as well
       * @return (child_ordinal, child_edge_ordinal)
       */
      std::pair<UInt,UInt> child_edge(const UInt edge_ordinal ,
                                      const UInt edge_child_ordinal,
                                      const Elem::CellTopology & objTop) const;


    public: // private:
      UInt                                  m_numChild;
      UInt                                  m_numChildNodes;
      const CellTopology *                  m_childCellTopology;
      UInt * *                              m_childNode;
      bool                                  m_homogeneous; /*whether refinement is self simmilar*/
      bool                                  m_fullRefinement; /*whether all edges are to be refined or not*/

    private:
      MeshObjRefinementTopology(const MeshObjRefinementTopology&);
      MeshObjRefinementTopology&operator=(const MeshObjRefinementTopology&);
    };


    const RefinementTopology *getRefinementTopology(const Elem::CellTopology &cell_topology);
    const UInt *getRefinementEdgeNode(const Elem::CellTopology &cell_topology, UInt edge);
    const UInt *getRefinementFaceNode(const Elem::CellTopology &cell_topology, UInt face);
    const UInt *getRefinementEdgePermutation(const Elem::CellTopology &cell_topology, UInt permutation_ordinal);

  } // namespace Elem
} // namespace adapt
} // namespace stk_classic

#endif // stk_adapt_sierra_element_RefinementTopology_hpp
