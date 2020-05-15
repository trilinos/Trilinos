// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef adapt_sierra_element_RefinementTopology_hpp
#define adapt_sierra_element_RefinementTopology_hpp

//#include <percept_code_types.h>

#include <utility>

#include <adapt/sierra_element/CellTopology.hpp>

namespace percept {
  namespace Elem {

    class MeshObjTopology;

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

      std::vector< std::pair<UInt,UInt> >
      get_children_on_ordinal(const UInt face_ordinal) const;

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


    const RefinementTopology *getRefinementTopology(const Elem::CellTopology &cell_topology);
    const UInt *getRefinementEdgeNode(const Elem::CellTopology &cell_topology, UInt edge);
    const UInt *getRefinementFaceNode(const Elem::CellTopology &cell_topology, UInt face);
    const UInt *getRefinementEdgePermutation(const Elem::CellTopology &cell_topology, UInt permutation_ordinal);

  } // namespace Elem
} // namespace percept

#endif // adapt_sierra_element_RefinementTopology_hpp
