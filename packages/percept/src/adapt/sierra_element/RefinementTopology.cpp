// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.




#include <algorithm>
#include <vector>
#include <map>

#include <adapt/sierra_element/percept_code_types.hpp>

//#include <stk_util/parallel/Exception.hpp>

#include <adapt/sierra_element/MeshObjTopology.hpp>
#include <adapt/sierra_element/RefinementTopology.hpp>
//#include <adapt/sierra_element/RefnementKey.hpp>
#include <adapt/sierra_element/CellTopology.hpp>

#include <percept/Util.hpp>

using namespace percept;

  namespace percept {
    namespace Elem {

    namespace {

      typedef std::map<CellTopology, const RefinementTopology *> CellTopologyRefinementMap;

      CellTopologyRefinementMap &
      get_cell_topology_refinement_map()
      {
        /* %TRACE[SPEC]% */  /* %TRACE% */
        static CellTopologyRefinementMap s_cellTopologyRefinementMap;

        return s_cellTopologyRefinementMap;
      }

    } // namespace <<unnamed>


    const RefinementTopology *
    getRefinementTopology(
                          const Elem::CellTopology &    cell_topology)
    {
      return get_cell_topology_refinement_map()[cell_topology];
    }


    const UInt *
    getRefinementEdgeNode(
                          const Elem::CellTopology &    cell_topology,
                          UInt                  edge)
    {
      return get_cell_topology_refinement_map()[cell_topology]->edge_node(edge);
    }


    const UInt *
    getRefinementFaceNode(
                          const Elem::CellTopology &    cell_topology,
                          UInt                  face)
    {
      return get_cell_topology_refinement_map()[cell_topology]->face_node(face);
    }


    const UInt *
    getRefinementEdgePermutation(
                                 const Elem::CellTopology &    cell_topology,
                                 UInt                  permutation_ordinal)
    {
      return get_cell_topology_refinement_map()[cell_topology] ? get_cell_topology_refinement_map()[cell_topology]->edge_permutation(permutation_ordinal) : 0;
    }


    RefinementTopology::RefinementTopology(
                                           MeshObjTopology *                     mesh_obj_topology,
                                           UInt                                  num_child,
                                           const MeshObjTopology * const *       child_topology,
                                           UInt                                  num_child_nodes,
                                           const UInt * const *                  child_nodes,
                                           const UInt                            num_edges, 
                                           const UInt * const *                  edge_node,
                                           const UInt                            num_faces, 
                                           const UInt * const *                  face_node,
                                           UInt                                  num_orientations,
                                           const UInt * const *                  perm_node,
                                           const UInt * const *                  perm_edge,
                                           bool                                  homogeneous_child)
      : m_cellTopology(mesh_obj_topology->getCellTopology()),
        m_numChild(num_child),
        m_childCellTopology(0),
        m_numChildNodes(num_child_nodes),
        m_childNode(child_nodes),
        m_numEdges(num_edges),
        m_edgeNode(edge_node),
        m_numFaces(num_faces),
        m_faceNode(face_node),
        m_numOrientations(num_orientations), 
        m_nodePermutation(perm_node),
        m_edgePermutation(perm_edge),
        m_homogeneousChild(homogeneous_child)
    {
      //  mesh_obj_topology->set_refinement_topology(this);

      m_childCellTopology = new CellTopology[num_child];
      for (size_t i = 0; i < num_child; ++i)
        m_childCellTopology[i] = child_topology[i]->getCellTopology();

      get_cell_topology_refinement_map()[mesh_obj_topology->getCellTopology()] = this;
    }

    RefinementTopology::~RefinementTopology()
    {
      delete[] m_childCellTopology;
    }


    const UInt *
    RefinementTopology::node_permutation(
                                         UInt                  permutation_ordinal) const
    {
      return permutation_ordinal < m_numOrientations ? m_nodePermutation[permutation_ordinal] : NULL;
    }


    const UInt *
    RefinementTopology::edge_permutation(
                                         UInt                  permutation_ordinal) const
    {
      return permutation_ordinal < m_numOrientations ? m_edgePermutation[permutation_ordinal] : NULL;
    }


    std::vector< std::pair<UInt,UInt> >
    RefinementTopology::get_children_on_ordinal(const UInt face_ordinal) const
    {
      std::vector< std::pair<UInt,UInt> > result;
      std::pair<UInt,UInt> entry(0,0);

      bool flag = (face_ordinal < m_cellTopology.getFaceCount());

      if (flag) {

        const CellTopology & faceTop = faceCellTopology(m_cellTopology, face_ordinal);

        const UInt * const faceN      = face_node(face_ordinal);

        for(unsigned face_child_ordinal=0; face_child_ordinal < getRefinementTopology(faceCellTopology(m_cellTopology, face_ordinal))->num_child(); ++face_child_ordinal)
        {
            const UInt * const faceChildN = Elem::getRefinementTopology(faceTop)->child_node(face_child_ordinal);

            const CellTopology faceChildTop = getRefinementTopology(faceTop)->child_cell_topology(face_child_ordinal);

            flag = false;
            UInt childIndex, faceIndex;

            for (childIndex = 0; childIndex < num_child(); ++childIndex) {

              const UInt * const           childN   =   child_node(    childIndex);
              const CellTopology childTop = child_cell_topology(childIndex);

              for (faceIndex = 0; faceIndex < childTop.getFaceCount(); ++faceIndex) {

                if (faceCellTopology(childTop, faceIndex) == faceChildTop) {

                  const UInt * const childFaceN = Elem::getRefinementTopology(childTop)->face_node(faceIndex);

                  UInt i = 0;

                  for (; i < faceChildTop.getNodeCount() &&
                         faceN[ faceChildN[i] ] == childN[ childFaceN[i] ]; ++i);

                  if (i == faceChildTop.getNodeCount()) {
                    VERIFY_TRUE(! flag);
                    flag = true;

                    entry.first  = childIndex;
                    entry.second = faceIndex;

                    result.push_back(entry);
                  }
                }
              }
            }
          }
      }

      if (! flag) {
        //         throw RuntimeError() << "CellTopology[" << m_cellTopology.getName() << "]::child_face("
        //                              << face_ordinal << "," << face_child_ordinal
        //                              << ") could not find a match." << std::endl << StackTrace;
        THROW( "CellTopology[" << m_cellTopology.getName() << "]::get_children_on_ordinal("
               << face_ordinal
               << ") could not find a match." );
      }

      return result;
    }


  } // namespace Elem
} // namespace percept 

