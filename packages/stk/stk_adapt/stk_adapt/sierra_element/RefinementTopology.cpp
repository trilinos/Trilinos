/*--------------------------------------------------------------------*/
/*    Copyright 1999 - 2009 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

/* Copyright 1999-2002 Sandia Corporation, Albuquerque, NM. */


#include <algorithm>
#include <vector>
#include <map>

#include <stk_adapt/sierra_element/stk_percept_code_types.hpp>

//#include <stk_util/parallel/Exception.hpp>

#include <stk_adapt/sierra_element/MeshObjTopology.hpp>
#include <stk_adapt/sierra_element/RefinementTopology.hpp>
//#include <stk_adapt/sierra_element/RefnementKey.hpp>
#include <stk_adapt/sierra_element/CellTopology.hpp>

#include <stk_percept/Util.hpp>

namespace stk { 

  using namespace percept;

  namespace adapt {
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


    /*--------------------------------------------------------------------*/
    /*----MeshObjRefinementTopology----------------------------------*/
    /*--------------------------------------------------------------------*/
    MeshObjRefinementTopology::MeshObjRefinementTopology()
      : m_numChild(0),
        m_numChildNodes(0),
        m_childCellTopology(0),
        m_childNode(0),
        m_homogeneous(true),
        m_fullRefinement(false)
    {}


    MeshObjRefinementTopology::~MeshObjRefinementTopology()
    {
      //Destructor
      //if homogeneous do not delete data as they refer to static variables
      //if not homogeneous release refinement pattern data

      if(!m_homogeneous){

        //Delete 2 dimensional array of childNode
        for(UInt i = 0; i<m_numChild; ++i){
          delete [] m_childNode[i];
        }
        delete [] m_childNode;

        //delete [] m_childTopology;
      }
    }


    UInt
    MeshObjRefinementTopology::num_child() const
    {
      return m_numChild;
    }

    UInt
    MeshObjRefinementTopology::num_child_nodes()  const
    {
      return m_numChildNodes;
    }


    CellTopology
    MeshObjRefinementTopology::child_cell_topology(
                                                   UInt                  child)  const
    {
      //  return (child < m_numChild) ? m_childTopology[child]->getCellTopology() : CellTopology();
      return (child < m_numChild) ? m_childCellTopology[child] : CellTopology();
    }


    const UInt *
    MeshObjRefinementTopology::child_node(
                                          UInt                  child)  const
    {
      return (child < m_numChild) ? m_childNode[child] : NULL;
    }

    bool
    MeshObjRefinementTopology::homogeneous_refinement()  const
    {
      return m_homogeneous;
    }


    bool
    MeshObjRefinementTopology::full_refinement()  const
    {
      return m_fullRefinement;
    }


    std::pair<UInt,UInt>
    MeshObjRefinementTopology::child_face(
                                          const UInt            face_ordinal,
                                          const UInt            face_child_ordinal,
                                          const CellTopology &  objTop,
                                          const RefinementKey & objDesiredKey)  const
    {
      std::pair<UInt,UInt> result(0,0);

      bool flag = (face_ordinal < objTop.getFaceCount()) &&
        (face_child_ordinal < getRefinementTopology(faceCellTopology(objTop, face_ordinal))->num_child());

      if (flag) {

        // Match nodes of
        //   face_topology(face_ordinal)->child_nodes(face_child_node)[*]
        // to
        //   child_topology(face_child.first)->face_nodes(face_child.second)[*]



        const CellTopology faceTop = Elem::faceCellTopology(objTop, face_ordinal);

        const UInt * const faceN      = getRefinementTopology(objTop)->face_node(face_ordinal);

        //const UInt * const faceChildN = faceTop.child_node(face_child_ordinal);

        //faceChildN needs to be replaced with faceDesiredTop.child_node
        //Need face desired key

        std::vector<UInt> sideEdgePattern;
        sideEdgePattern.clear();
        std::vector<UInt> objEdgePattern = objDesiredKey.ordered_cut_edges(objTop.getEdgeCount());
        for(UInt jEdge = 0; jEdge<objEdgePattern.size();++jEdge){
          for(UInt iEdge = 0; iEdge< faceTop.getEdgeCount();++iEdge){
            Int faceEdgeOrdinal = getFaceEdge(objTop, face_ordinal,iEdge);
            if(faceEdgeOrdinal == (Int) objEdgePattern[jEdge]) {
              sideEdgePattern.push_back(iEdge);
            }
          }
        }

        RefinementKey faceDesiredKey;
        faceDesiredKey.assign_value(sideEdgePattern);
        MeshObjRefinementTopology faceDesiredTop;

        if (! getRefinementTopology(faceTop)->query_refinement_topology(faceDesiredKey,faceDesiredTop)) {
          //throw RuntimeError() << "Failed query of face refinement topology."
          //                     << std::endl << StackTrace;
          THROW( "Failed query of face refinement topology.");
        }

        const UInt * const faceChildN = faceDesiredTop.child_node(face_child_ordinal);

        const CellTopology faceChildTop = getRefinementTopology(faceTop)->child_cell_topology(face_child_ordinal);

        flag = false;
        std::pair<UInt,UInt> cf;

        for (cf.first = 0; cf.first < num_child(); ++cf.first) {
          const UInt * const           childN   =   child_node(    cf.first);
          const CellTopology         childTop = child_cell_topology(cf.first);

          for (cf.second = 0; cf.second < childTop.getFaceCount(); ++cf.second) {
            if (Elem::faceCellTopology(childTop, cf.second) == faceChildTop) {

              const UInt * const childFaceN = getRefinementTopology(childTop)->face_node(cf.second);

              UInt i = 0;

              for (; i < faceChildTop.getNodeCount() &&
                     faceN[ faceChildN[i] ] == childN[ childFaceN[i] ]; ++i);

              if (i == faceChildTop.getNodeCount()) {
                VERIFY_TRUE_ON(! flag);
                flag = true;
                result.first  = cf.first;
                result.second = cf.second;
              }
            }
          }
        }
      }

      if (! flag)
        {
          //         throw RuntimeError() << "MeshObjRefinementTopology[" << objTop.getName() << "]::child_face("
          //                              << face_ordinal << "," << face_child_ordinal
          //                              << ") could not find a match." << std::endl << StackTrace;
          THROW("MeshObjRefinementTopology[" << objTop.getName() << "]::child_face("
                << face_ordinal << "," << face_child_ordinal
                << ") could not find a match." );

        }

      return result;
    }


    std::pair<UInt,UInt>
    MeshObjRefinementTopology::child_edge(
                                          const UInt edge_ordinal,
                                          const UInt edge_child_ordinal,
                                          const CellTopology & objTop)  const
    {
      std::pair<UInt,UInt> result(0,0);

      bool flag = (edge_ordinal < objTop.getEdgeCount()) &&
        (edge_child_ordinal < getRefinementTopology(edgeCellTopology(objTop, edge_ordinal))->num_child());

      if (flag) {

        // Match nodes of
        //   edge_topology(edge_ordinal)->child_nodes(edge_child_node)[*]
        // to
        //   child_topology(edge_child.first)->edge_nodes(edge_child.second)[*]

        const CellTopology edgeTop = edgeCellTopology(objTop, edge_ordinal);

        const UInt * const edgeN      = getRefinementTopology(objTop)->edge_node(edge_ordinal);
        const UInt * const edgeChildN = getRefinementTopology(edgeTop)->child_node(edge_child_ordinal);

        const CellTopology edgeChildTop = getRefinementTopology(edgeTop)->child_cell_topology(edge_child_ordinal);

        flag = false;
        std::pair<UInt,UInt> ce;
        for (ce.first = 0; ce.first < num_child(); ++ce.first) {

          const UInt * const           childN   = child_node(    ce.first);
          const CellTopology childTop = child_cell_topology(ce.first);

          for (ce.second = 0; ce.second < childTop.getEdgeCount(); ++ce.second) {

            if (Elem::edgeCellTopology(childTop, ce.second) == edgeChildTop) {

              const UInt * const childEdgeN = getRefinementTopology(childTop)->edge_node(ce.second);
              UInt i = 0;

              for (; i < edgeChildTop.getNodeCount() &&
                     edgeN[ edgeChildN[i] ] == childN[ childEdgeN[i] ]; ++i);


              if (i == edgeChildTop.getNodeCount()) {
                flag = true;
                result.first  = ce.first;
                result.second = ce.second;

              }
            }
          }
        }
      }

      if (! flag)
        {
          //         throw RuntimeError() << "MeshObjRefinementTopology[" << objTop.getName() << "]::child_edge("
          //                              << edge_ordinal << "," << edge_child_ordinal
          //                              << ") could not find a match." << std::endl << StackTrace;
          THROW("MeshObjRefinementTopology[" << objTop.getName() << "]::child_edge("
                << edge_ordinal << "," << edge_child_ordinal
                << ") could not find a match." );
        }

      return result;
    }


    std::pair<UInt,UInt>
    RefinementTopology::child_face(
                                   const UInt face_ordinal,
                                   const UInt face_child_ordinal) const
    {
      std::pair<UInt,UInt> result(0,0);

      bool flag = (face_ordinal < m_cellTopology.getFaceCount()) &&
        (face_child_ordinal < getRefinementTopology(faceCellTopology(m_cellTopology, face_ordinal))->num_child());

      if (flag) {

        // Match nodes of
        //   face_topology(face_ordinal)->child_nodes(face_child_node)[*]
        // to
        //   child_topology(face_child.first)->face_nodes(face_child.second)[*]

        const CellTopology & faceTop = faceCellTopology(m_cellTopology, face_ordinal);

        const UInt * const faceN      = face_node(face_ordinal);

        const UInt * const faceChildN = Elem::getRefinementTopology(faceTop)->child_node(face_child_ordinal);

        const CellTopology faceChildTop = getRefinementTopology(faceTop)->child_cell_topology(face_child_ordinal);

        flag = false;
        std::pair<UInt,UInt> cf;

        for (cf.first = 0; cf.first < num_child(); ++cf.first) {

          const UInt * const           childN   =   child_node(    cf.first);
          const CellTopology childTop = child_cell_topology(cf.first);

          for (cf.second = 0; cf.second < childTop.getFaceCount(); ++cf.second) {

            if (faceCellTopology(childTop, cf.second) == faceChildTop) {

              const UInt * const childFaceN = Elem::getRefinementTopology(childTop)->face_node(cf.second);

              UInt i = 0;

              for (; i < faceChildTop.getNodeCount() &&
                     faceN[ faceChildN[i] ] == childN[ childFaceN[i] ]; ++i);

              if (i == faceChildTop.getNodeCount()) {
                VERIFY_TRUE(! flag);
                flag = true;
                result.first  = cf.first;
                result.second = cf.second;
              }
            }
          }
        }
      }

      if (! flag) {
        //         throw RuntimeError() << "CellTopology[" << m_cellTopology.getName() << "]::child_face("
        //                              << face_ordinal << "," << face_child_ordinal
        //                              << ") could not find a match." << std::endl << StackTrace;
        THROW( "CellTopology[" << m_cellTopology.getName() << "]::child_face("
               << face_ordinal << "," << face_child_ordinal
               << ") could not find a match." );
      }

      return result;
    }


    std::pair<UInt,UInt>
    RefinementTopology::child_edge(
                                   const UInt    edge_ordinal,
                                   const UInt    edge_child_ordinal) const
    {
      std::pair<UInt,UInt> result(0,0);

      bool flag = (edge_ordinal < m_cellTopology.getEdgeCount()) &&
        (edge_child_ordinal < getRefinementTopology(edgeCellTopology(m_cellTopology, edge_ordinal))->num_child());

      if (flag) {

        // Match nodes of
        //   edge_topology(edge_ordinal)->child_nodes(edge_child_node)[*]
        // to
        //   child_topology(edge_child.first)->edge_nodes(edge_child.second)[*]

        const CellTopology & edgeTop = edgeCellTopology(m_cellTopology, edge_ordinal);

        const UInt * const edgeN      = edge_node(edge_ordinal);
        const UInt * const edgeChildN = getRefinementTopology(edgeTop)->child_node(edge_child_ordinal);

        const CellTopology edgeChildTop = getRefinementTopology(edgeTop)->child_cell_topology(edge_child_ordinal);

        flag = false;
        std::pair<UInt,UInt> ce;

        for (ce.first = 0; ce.first < num_child(); ++ce.first) {

          const UInt * const           childN   = child_node(    ce.first);
          const CellTopology childTop = child_cell_topology(ce.first);

          for (ce.second = 0; ce.second < childTop.getEdgeCount(); ++ce.second) {

            if (edgeCellTopology(childTop, ce.second) == edgeChildTop) {

              const UInt * const childEdgeN = getRefinementTopology(childTop)->edge_node(ce.second);

              UInt i = 0;

              for (; i < edgeChildTop.getNodeCount() &&
                     edgeN[ edgeChildN[i] ] == childN[ childEdgeN[i] ]; ++i);

              if (i == edgeChildTop.getNodeCount()) {
                VERIFY_TRUE(! flag);
                flag = true;
                result.first  = ce.first;
                result.second = ce.second;
              }
            }
          }
        }
      }

      if (! flag) {
        //         throw RuntimeError() << "CellTopology[" << m_cellTopology.getName() << "]::child_edge("
        //                              << edge_ordinal << "," << edge_child_ordinal
        //                              << ") could not find a match." << std::endl << StackTrace;
        THROW("CellTopology[" << m_cellTopology.getName() << "]::child_edge("
              << edge_ordinal << "," << edge_child_ordinal
              << ") could not find a match." );
      }

      return result;
    }


    bool
    RefinementTopology::query_refinement_topology(
                                                  RefinementKey &               objKey,
                                                  MeshObjRefinementTopology &   objRefTop) const
    {
      // check if partial refinement key defined
      if (0x0 != objKey.value()) {
        // check if tri or tet element
        if (this == getRefinementTopology(getCellTopology<shards::Triangle<3> >())) {
          return refine_rivara_tri(objKey, objRefTop);
        }
        else if (this == getRefinementTopology(getCellTopology<shards::Tetrahedron<4> >())) {
          return refine_rivara_tet(objKey, objRefTop);
        }
      }
      else {
        //Return normal homogeneous patterns due to default value of 0x0
        //or the fact that obj topology name is currently not supported for partial refinement
        //Need to return differnet pattens for all normal homogeneous objects
        //these refinement patterns will not be dynamically created but will
        //point to or be equal to the current static objects;
        objRefTop.m_homogeneous = true;
        objRefTop.m_fullRefinement = true;
        objRefTop.m_numChild = num_child();
        objRefTop.m_numChildNodes= num_child_nodes();
        objRefTop.m_childCellTopology = child_cell_topology();
        objRefTop.m_childNode = (UInt * *) child_nodes();
        return true;
      }
      return false;
    }


    // Performs logic for bisection of triangle based on refinement key
    // Bisections are performed in order of largest to smallest edge to be refined
    bool
    RefinementTopology::refine_rivara_tri(
                                          RefinementKey &               objKey,
                                          MeshObjRefinementTopology &   objRefTop) const
    {
      //partial refinement key defined so create pattern

      //refTop is not based on standard homogeneous patterns
      objRefTop.m_homogeneous = false;

      UInt objNumNodes = m_cellTopology.getNodeCount();

      //Vector of edges and order in which to be refined
      std::vector<UInt> edgePattern = objKey.ordered_cut_edges(m_cellTopology.getEdgeCount());

      UInt num_edges_cut = edgePattern.size();

      //number of child equals num_edges_cut +1
      objRefTop.m_numChild =num_edges_cut+1;

      //number of child nodes = m_cellTopology.getNodeCount() + num_edges_cut for 2D triangle
      objRefTop.m_numChildNodes = objNumNodes+num_edges_cut;

      //Determine childNode for different templates based on largest edge bisection
      UInt longEdge = edgePattern[0];

      //For first bisection pNode is new child node need to find the node
      //that is the opposite vertex of pNode
      UInt node0 = edge_node(longEdge)[0];
      UInt node1 = edge_node(longEdge)[1];
      UInt pNode = edge_node(longEdge)[2];
      UInt oppVertex = objNumNodes;
      for (UInt iNode = 0; iNode < objNumNodes; ++iNode) {
        if(iNode != node0 && iNode != node1) oppVertex = iNode;
      }
      VERIFY_TRUE(oppVertex != m_cellTopology.getNodeCount());

      //Define a vector of of vectors that hold child node maps
      //Store all possible child node maps to be genereated
      //First Bisection gives 2 child elements
      //Then each of those elements gives possibly 2 more child

      std::vector<std::vector<UInt> >  childNodeVec(6, std::vector<UInt>(3));

      //An array that tells wether a given possible child node map is active
      bool childActive[6]={false,false,false,false,false,false};

      //First Bisection

      childNodeVec[0][node0]=node0;
      childNodeVec[0][node1]=pNode;
      childNodeVec[0][oppVertex]=oppVertex;

      childNodeVec[1][node0]=pNode;
      childNodeVec[1][node1]=node1;
      childNodeVec[1][oppVertex]=oppVertex;

      childActive[0]=true;
      childActive[1]=true;

      UInt iChildFilled = 2;

      //Now determine further bisection if needed
      //If one of the elements formed in first bisection is refined than
      //They will no longer be active

      for (UInt iEdge=1;iEdge<edgePattern.size();++iEdge)
        {
          UInt edgeCut = edgePattern[iEdge];
          UInt edgeNode0 = edge_node(edgeCut)[0];
          UInt edgeNode1 = edge_node(edgeCut)[1];
          UInt qNode = edge_node(edgeCut)[2];

          //Determine which child element of first bisection is to be refined
          int iChild = -1;
          if(edgeNode0 == node0 || edgeNode1 == node0){
            iChild = 0;
          }
          if(edgeNode0 == node1 || edgeNode1 == node1){
            iChild = 1;
          }
          VERIFY_TRUE(iChild != -1);

          //Find childOppVertex
          UInt childOppVertex = num_child_nodes();
          UInt childOppVertexOrdinal = num_child_nodes();
          for(UInt iNode =0; iNode<objNumNodes;++iNode){
            if(childNodeVec[iChild][iNode] != edgeNode0 &&
               childNodeVec[iChild][iNode] != edgeNode1) {
              childOppVertexOrdinal = iNode;
              childOppVertex = childNodeVec[iChild][iNode];
            }
          }
          VERIFY_TRUE(childOppVertex != num_child_nodes());
          VERIFY_TRUE(childOppVertexOrdinal != num_child_nodes());

          //first child of recursive bisection
          childNodeVec[iChildFilled][edgeNode0]=edgeNode0;
          childNodeVec[iChildFilled][edgeNode1]=qNode;
          childNodeVec[iChildFilled][childOppVertexOrdinal]=childOppVertex;
          //second child of recursive bisection
          childNodeVec[iChildFilled+1][edgeNode0]=qNode;
          childNodeVec[iChildFilled+1][edgeNode1]=edgeNode1;
          childNodeVec[iChildFilled+1][childOppVertexOrdinal]=childOppVertex;

          childActive[iChild]=false;
          childActive[iChildFilled]=true;
          childActive[iChildFilled+1]=true;
          iChildFilled = iChildFilled + 2;
        }

      //Assign new childNode maps to MeshObjRefinementTopology
      objRefTop.m_childNode = new UInt *[objRefTop.m_numChild];
      UInt iChild=0;
      for(UInt i=0;i<6;++i){
        if(childActive[i]){
          objRefTop.m_childNode[iChild]=new UInt[4];
          for(UInt j=0; j <objNumNodes; ++j){
            objRefTop.m_childNode[iChild][j]=childNodeVec[i][j];

          }
          ++iChild;
        }
      }

      if(num_edges_cut == 3) {
        objRefTop.m_fullRefinement = true;
      }
      else{
        objRefTop.m_fullRefinement = false;
      }

      //Since partial refinement child topologies are same as parent
      //so use standard refinement pattern for this.  For heterogeneous elements
      //childTopology would need to be dynamically created
      objRefTop.m_childCellTopology = child_cell_topology();
      return true;
    }


    bool
    RefinementTopology::refine_rivara_tet(
                                          RefinementKey &               objKey,
                                          MeshObjRefinementTopology &   objRefTop) const
    {
      objRefTop.m_homogeneous = false;

      UInt objNumEdges = m_cellTopology.getEdgeCount();

      UInt objNumNodes = m_cellTopology.getNodeCount();

      std::vector<UInt> edgePattern = objKey.ordered_cut_edges(objNumEdges);

      objRefTop.m_numChild = 1;
      objRefTop.m_numChildNodes = objNumNodes;
      UInt longEdge = edgePattern[0];
      UInt node0 = edge_node(longEdge)[0];
      UInt node1 = edge_node(longEdge)[1];
      UInt pNode = edge_node(longEdge)[2];

      UInt oppEdge = objNumEdges;
      for(UInt i=0;i<objNumEdges;++i){
        if(node0 != edge_node(i)[0] && node0 != edge_node(i)[1] &&
           node1 != edge_node(i)[0] && node1 != edge_node(i)[1]) oppEdge = i;
      }
      VERIFY_TRUE(oppEdge != objNumEdges);

      UInt oppEdgeNode0 = edge_node(oppEdge)[0];
      UInt oppEdgeNode1 = edge_node(oppEdge)[1];

      std::vector<std::vector<UInt> >  childNodeVec(14, std::vector<UInt>(4));

      bool childActive[14]={false,false,false,false,false,false,
                            false,false,false,false,false,false,false,false};

      //First Bisection
      objRefTop.m_numChild = objRefTop.m_numChild+1;
      objRefTop.m_numChildNodes = objRefTop.m_numChildNodes +1;

      childNodeVec[0][node0]=node0;
      childNodeVec[0][node1]=pNode;
      childNodeVec[0][oppEdgeNode0]=oppEdgeNode0;
      childNodeVec[0][oppEdgeNode1]=oppEdgeNode1;

      childNodeVec[1][node0]=pNode;
      childNodeVec[1][node1]=node1;
      childNodeVec[1][oppEdgeNode0]=oppEdgeNode0;
      childNodeVec[1][oppEdgeNode1]=oppEdgeNode1;

      childActive[0]=true;
      childActive[1]=true;

      UInt iChildFilled = 2;

      //After first bisection iterate over edges to be refined
      for(UInt iEdge=1;iEdge<edgePattern.size();++iEdge){
        UInt edgeCut = edgePattern[iEdge];
        UInt edgeNode0 = edge_node(edgeCut)[0];
        UInt edgeNode1 = edge_node(edgeCut)[1];
        UInt qNode = edge_node(edgeCut)[2];

        objRefTop.m_numChildNodes = objRefTop.m_numChildNodes+1;

        //Find child that this edge exists on
        //All child elements that this edge  belong to need to be bisected
        for(UInt iChild = 0; iChild<childNodeVec.size();++iChild){
          if(childActive[iChild]){
            bool edgeNode0Found = false;
            bool edgeNode1Found = false;
            for(UInt iNode =0; iNode < objNumNodes; ++iNode){
              if(edgeNode0==childNodeVec[iChild][iNode]) edgeNode0Found = true;
              if(edgeNode1==childNodeVec[iChild][iNode]) edgeNode1Found = true;
            }
            if(edgeNode0Found && edgeNode1Found){

              //This child is to be bisected
              //Get m_edgeNode for this child tetrahedral
              UInt childEdge0[]={childNodeVec[iChild][0], childNodeVec[iChild][1]};
              UInt childEdge1[]={childNodeVec[iChild][1], childNodeVec[iChild][2]};
              UInt childEdge2[]={childNodeVec[iChild][2], childNodeVec[iChild][0]};
              UInt childEdge3[]={childNodeVec[iChild][0], childNodeVec[iChild][3]};
              UInt childEdge4[]={childNodeVec[iChild][1], childNodeVec[iChild][3]};
              UInt childEdge5[]={childNodeVec[iChild][2], childNodeVec[iChild][3]};
              UInt * newChildEdge[]={childEdge0,childEdge1,childEdge2,childEdge3,
                                     childEdge4,childEdge5};


              //Find oppEdgeNode0 and oppEdgeNode1 for this child Tet;
              UInt childOppEdge = objNumEdges;
              for(UInt i = 0; i <objNumEdges; ++i){
                if(edgeNode0 != newChildEdge[i][0] && edgeNode0 != newChildEdge[i][1] &&
                   edgeNode1 != newChildEdge[i][0] && edgeNode1 != newChildEdge[i][1]
                   ) childOppEdge = i;
              }

              UInt childOppEdgeNode0 = newChildEdge[childOppEdge][0];
              UInt childOppEdgeNode1 = newChildEdge[childOppEdge][1];

              objRefTop.m_numChild = objRefTop.m_numChild+1;

              //Two new child element node maps for this bisection
              //Deactivate element that these come frome
              childNodeVec[iChildFilled][edgeNode0]=edgeNode0;
              childNodeVec[iChildFilled][edgeNode1]=qNode;
              childNodeVec[iChildFilled][edge_node(childOppEdge)[0]]=childOppEdgeNode0;
              childNodeVec[iChildFilled][edge_node(childOppEdge)[1]]=childOppEdgeNode1;

              childNodeVec[iChildFilled+1][edgeNode0]=qNode;
              childNodeVec[iChildFilled+1][edgeNode1]=edgeNode1;
              childNodeVec[iChildFilled+1][edge_node(childOppEdge)[0]]=childOppEdgeNode0;
              childNodeVec[iChildFilled+1][edge_node(childOppEdge)[1]]=childOppEdgeNode1;


              childActive[iChild]=false;
              childActive[iChildFilled]=true;
              childActive[iChildFilled+1]=true;
              iChildFilled = iChildFilled + 2;
            }
          }
        }
      }

      objRefTop.m_childNode = new UInt *[objRefTop.m_numChild];
      UInt iChild=0;
      for(UInt i=0;i<14;++i){
        if(childActive[i]){
          objRefTop.m_childNode[iChild]=new UInt[4];
          for(UInt j=0; j <objNumNodes; ++j){
            objRefTop.m_childNode[iChild][j]=childNodeVec[i][j];

          }
          ++iChild;
        }
      }

      objRefTop.m_childCellTopology = child_cell_topology();

      if(objRefTop.m_numChild ==8){
        objRefTop.m_fullRefinement = true;
      }
      else{
        objRefTop.m_fullRefinement = false;
      }

      return true;
    }


  } // namespace Elem
} // namespace adapt 
} // namespace stk
