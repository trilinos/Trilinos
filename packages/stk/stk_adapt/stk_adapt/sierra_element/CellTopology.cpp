/*------------------------------------------------------------------------*/
/*               shards : Shared Discretization Tools                     */
/*                Copyright (2008, 2011) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/* Questions? Contact Pavel Bochev      (pbboche@sandia.gov)              */
/*                    H. Carter Edwards (hcedwar@sandia.gov)              */
/*                    Denis Ridzal      (dridzal@sandia.gov).             */
/*------------------------------------------------------------------------*/

//#define HAVE_SHARDS_DEBUG

#include <map>
#include <stdexcept>
#include <sstream>

#include <Shards_BasicTopologies.hpp>

//#include <stk_util/parallel/ExceptionReport.hpp>
//#include <stk_util/diag/StringUtil.hpp>

//#include <element/Elem_CellTopology.h>
#include <stk_adapt/sierra_element/CellTopology.hpp>

using namespace shards;

namespace stk { namespace adapt {
  namespace Elem {

    namespace {

      //srk typedef std::map<std::string, CellTopology, less_nocase<std::string> > CellTopologyNameMap;
      typedef std::map<std::string, CellTopology > CellTopologyNameMap;
      typedef std::map<TopologyId, CellTopology> CellTopologyIdMap;
      typedef std::map<CellTopology, TopologyId> IdCellTopologyMap;

      CellTopologyNameMap &
      get_cell_topology_name_map()
      {
        static CellTopologyNameMap        s_cellTopologyMap;

        if (s_cellTopologyMap.empty()) {
          s_cellTopologyMap[getCellTopology<Node>().getName()] = getCellTopology<Node>();
          s_cellTopologyMap[getCellTopology<Particle>().getName()] = getCellTopology<Particle>();
          s_cellTopologyMap[getCellTopology<Line<2> >().getName()] = getCellTopology<Line<2> >();
          s_cellTopologyMap[getCellTopology<Line<3> >().getName()] = getCellTopology<Line<3> >();
          s_cellTopologyMap[getCellTopology<ShellLine<2> >().getName()] = getCellTopology<ShellLine<2> >();
          s_cellTopologyMap[getCellTopology<ShellLine<3> >().getName()] = getCellTopology<ShellLine<3> >();
          s_cellTopologyMap[getCellTopology<Beam<2> >().getName()] = getCellTopology<Beam<2> >();
          s_cellTopologyMap[getCellTopology<Beam<3> >().getName()] = getCellTopology<Beam<3> >();
          s_cellTopologyMap[getCellTopology<Triangle<3> >().getName()] = getCellTopology<Triangle<3> >();
          s_cellTopologyMap[getCellTopology<Triangle<4> >().getName()] = getCellTopology<Triangle<4> >();
          s_cellTopologyMap[getCellTopology<Triangle<6> >().getName()] = getCellTopology<Triangle<6> >();
          s_cellTopologyMap[getCellTopology<ShellTriangle<3> >().getName()] = getCellTopology<ShellTriangle<3> >();
          //    s_cellTopologyMap[getCellTopology<ShellTriangle<4> >().getName()] = getCellTopology<ShellTriangle<4> >();
          s_cellTopologyMap[getCellTopology<ShellTriangle<6> >().getName()] = getCellTopology<ShellTriangle<6> >();
          s_cellTopologyMap[getCellTopology<Quadrilateral<4> >().getName()] = getCellTopology<Quadrilateral<4> >();
          s_cellTopologyMap[getCellTopology<Quadrilateral<8> >().getName()] = getCellTopology<Quadrilateral<8> >();
          s_cellTopologyMap[getCellTopology<Quadrilateral<9> >().getName()] = getCellTopology<Quadrilateral<9> >();
          s_cellTopologyMap[getCellTopology<ShellQuadrilateral<4> >().getName()] = getCellTopology<ShellQuadrilateral<4> >();
          s_cellTopologyMap[getCellTopology<ShellQuadrilateral<8> >().getName()] = getCellTopology<ShellQuadrilateral<8> >();
          s_cellTopologyMap[getCellTopology<ShellQuadrilateral<9> >().getName()] = getCellTopology<ShellQuadrilateral<9> >();
          s_cellTopologyMap[getCellTopology<Tetrahedron<4> >().getName()] = getCellTopology<Tetrahedron<4> >();
          s_cellTopologyMap[getCellTopology<Tetrahedron<8> >().getName()] = getCellTopology<Tetrahedron<8> >();
          s_cellTopologyMap[getCellTopology<Tetrahedron<10> >().getName()] = getCellTopology<Tetrahedron<10> >();
          s_cellTopologyMap[getCellTopology<Hexahedron<8> >().getName()] = getCellTopology<Hexahedron<8> >();
          s_cellTopologyMap[getCellTopology<Hexahedron<20> >().getName()] = getCellTopology<Hexahedron<20> >();
          s_cellTopologyMap[getCellTopology<Hexahedron<27> >().getName()] = getCellTopology<Hexahedron<27> >();
          s_cellTopologyMap[getCellTopology<Pyramid<5> >().getName()] = getCellTopology<Pyramid<5> >();
          s_cellTopologyMap[getCellTopology<Pyramid<13> >().getName()] = getCellTopology<Pyramid<13> >();
          s_cellTopologyMap[getCellTopology<Pyramid<14> >().getName()] = getCellTopology<Pyramid<14> >();
          s_cellTopologyMap[getCellTopology<Wedge<6> >().getName()] = getCellTopology<Wedge<6> >();
          s_cellTopologyMap[getCellTopology<Wedge<15> >().getName()] = getCellTopology<Wedge<15> >();
          s_cellTopologyMap[getCellTopology<Wedge<18> >().getName()] = getCellTopology<Wedge<18> >();
          s_cellTopologyMap[getCellTopology<Pentagon<5> >().getName()] = getCellTopology<Pentagon<5> >();
          s_cellTopologyMap[getCellTopology<Hexagon<6> >().getName()] = getCellTopology<Hexagon<6> >();
        }
  
        return s_cellTopologyMap;
      }


      CellTopologyIdMap &
      get_cell_topology_id_map()
      {
        static CellTopologyIdMap        s_cellTopologyMap;

        if (s_cellTopologyMap.empty()) {
          s_cellTopologyMap[NODE_0] = getCellTopology<Node>();
          s_cellTopologyMap[PARTICLE_1] = getCellTopology<Particle>();
          s_cellTopologyMap[EDGE_2] = getCellTopology<Line<2> >();
          s_cellTopologyMap[EDGE_3] = getCellTopology<Line<3> >();
          s_cellTopologyMap[SHELL_LINE_2] = getCellTopology<ShellLine<2> >();
          s_cellTopologyMap[SHELL_LINE_3] = getCellTopology<ShellLine<3> >();
          s_cellTopologyMap[ROD_2] = getCellTopology<Beam<2> >();
          s_cellTopologyMap[ROD_3] = getCellTopology<Beam<3> >();
          s_cellTopologyMap[FACE_TRI_3] = getCellTopology<Triangle<3> >();
          s_cellTopologyMap[FACE_TRI_4] = getCellTopology<Triangle<4> >();
          s_cellTopologyMap[FACE_TRI_6] = getCellTopology<Triangle<6> >();
          s_cellTopologyMap[SHELL_TRI_3] = getCellTopology<ShellTriangle<3> >();
          //    s_cellTopologyMap[SHELL_TRI_4] = getCellTopology<ShellTriangle<4> >();
          s_cellTopologyMap[SHELL_TRI_6] = getCellTopology<ShellTriangle<6> >();
          s_cellTopologyMap[FACE_QUAD_4] = getCellTopology<Quadrilateral<4> >();
          s_cellTopologyMap[FACE_QUAD_8] = getCellTopology<Quadrilateral<8> >();
          s_cellTopologyMap[FACE_QUAD_9] = getCellTopology<Quadrilateral<9> >();
          s_cellTopologyMap[SHELL_QUAD_4] = getCellTopology<ShellQuadrilateral<4> >();
          s_cellTopologyMap[SHELL_QUAD_8] = getCellTopology<ShellQuadrilateral<8> >();
          s_cellTopologyMap[SHELL_QUAD_9] = getCellTopology<ShellQuadrilateral<9> >();
          s_cellTopologyMap[SOLID_TET_4] = getCellTopology<Tetrahedron<4> >();
          s_cellTopologyMap[SOLID_TET_8] = getCellTopology<Tetrahedron<8> >();
          s_cellTopologyMap[SOLID_TET_10] = getCellTopology<Tetrahedron<10> >();
          s_cellTopologyMap[SOLID_HEX_8] = getCellTopology<Hexahedron<8> >();
          s_cellTopologyMap[SOLID_HEX_20] = getCellTopology<Hexahedron<20> >();
          s_cellTopologyMap[SOLID_HEX_27] = getCellTopology<Hexahedron<27> >();
          s_cellTopologyMap[SOLID_PYRAMID_5] = getCellTopology<Pyramid<5> >();
          s_cellTopologyMap[SOLID_PYRAMID_13] = getCellTopology<Pyramid<13> >();
          s_cellTopologyMap[SOLID_PYRAMID_14] = getCellTopology<Pyramid<14> >();
          s_cellTopologyMap[SOLID_WEDGE_6] = getCellTopology<Wedge<6> >();
          s_cellTopologyMap[SOLID_WEDGE_15] = getCellTopology<Wedge<15> >();
          s_cellTopologyMap[SOLID_WEDGE_18] = getCellTopology<Wedge<18> >();
          s_cellTopologyMap[FACE_PENT_5] = getCellTopology<Pentagon<5> >();
          s_cellTopologyMap[FACE_HEX_6] = getCellTopology<Hexagon<6> >();
        }
  
        return s_cellTopologyMap;
      }


      IdCellTopologyMap &
      get_id_cell_topology_map()
      {
        static IdCellTopologyMap        s_cellTopologyMap;

        if (s_cellTopologyMap.empty()) {
          s_cellTopologyMap[getCellTopology<Node>()] = NODE_0;
          s_cellTopologyMap[getCellTopology<Particle>()] = PARTICLE_1;
          s_cellTopologyMap[getCellTopology<Line<2> >()] = EDGE_2;
          s_cellTopologyMap[getCellTopology<Line<3> >()] = EDGE_3;
          s_cellTopologyMap[getCellTopology<ShellLine<2> >()] = SHELL_LINE_2;
          s_cellTopologyMap[getCellTopology<ShellLine<3> >()] = SHELL_LINE_3;
          s_cellTopologyMap[getCellTopology<Beam<2> >()] = ROD_2;
          s_cellTopologyMap[getCellTopology<Beam<3> >()] = ROD_3;
          s_cellTopologyMap[getCellTopology<Triangle<3> >()] = FACE_TRI_3;
          s_cellTopologyMap[getCellTopology<Triangle<4> >()] = FACE_TRI_4;
          s_cellTopologyMap[getCellTopology<Triangle<6> >()] = FACE_TRI_6;
          s_cellTopologyMap[getCellTopology<ShellTriangle<3> >()] = SHELL_TRI_3;
          //    s_cellTopologyMap[getCellTopology<ShellTriangle<4> >()] = SHELL_TRI_4;
          s_cellTopologyMap[getCellTopology<ShellTriangle<6> >()] = SHELL_TRI_6;
          s_cellTopologyMap[getCellTopology<Quadrilateral<4> >()] = FACE_QUAD_4;
          s_cellTopologyMap[getCellTopology<Quadrilateral<8> >()] = FACE_QUAD_8;
          s_cellTopologyMap[getCellTopology<Quadrilateral<9> >()] = FACE_QUAD_9;
          s_cellTopologyMap[getCellTopology<ShellQuadrilateral<4> >()] = SHELL_QUAD_4;
          s_cellTopologyMap[getCellTopology<ShellQuadrilateral<8> >()] = SHELL_QUAD_8;
          s_cellTopologyMap[getCellTopology<ShellQuadrilateral<9> >()] = SHELL_QUAD_9;
          s_cellTopologyMap[getCellTopology<Tetrahedron<4> >()] = SOLID_TET_4;
          s_cellTopologyMap[getCellTopology<Tetrahedron<8> >()] = SOLID_TET_8;
          s_cellTopologyMap[getCellTopology<Tetrahedron<10> >()] = SOLID_TET_10;
          s_cellTopologyMap[getCellTopology<Hexahedron<8> >()] = SOLID_HEX_8;
          s_cellTopologyMap[getCellTopology<Hexahedron<20> >()] = SOLID_HEX_20;
          s_cellTopologyMap[getCellTopology<Hexahedron<27> >()] = SOLID_HEX_27;
          s_cellTopologyMap[getCellTopology<Pyramid<5> >()] = SOLID_PYRAMID_5;
          s_cellTopologyMap[getCellTopology<Pyramid<13> >()] = SOLID_PYRAMID_13;
          s_cellTopologyMap[getCellTopology<Pyramid<14> >()] = SOLID_PYRAMID_14;
          s_cellTopologyMap[getCellTopology<Wedge<6> >()] = SOLID_WEDGE_6;
          s_cellTopologyMap[getCellTopology<Wedge<15> >()] = SOLID_WEDGE_15;
          s_cellTopologyMap[getCellTopology<Wedge<18> >()] = SOLID_WEDGE_18;
          s_cellTopologyMap[getCellTopology<Pentagon<5> >()] = FACE_PENT_5;
          s_cellTopologyMap[getCellTopology<Hexagon<6> >()] = FACE_HEX_6;
        }
  
        return s_cellTopologyMap;
      }

    } // namespace <unnamed>


    CellTopology
    getBasicCellTopology(
                         const char *          name)
    {
      CellTopologyNameMap &cell_topology_map = get_cell_topology_name_map();

      CellTopologyNameMap::const_iterator it = cell_topology_map.find(name);
  
      if (it == cell_topology_map.end())
        {
          //throw RuntimeError() << "Cell topology " << name << " is not defined";
          std::ostringstream msg;
          msg <<  "Cell topology " << name << " is not defined";
          throw std::runtime_error(msg.str());
        }

      return (*it).second;
    }


    CellTopology
    getBasicCellTopology(
                         TopologyId            id)
    {
      if (id == INVALID)
        return CellTopology();
  
      CellTopologyIdMap &cell_topology_map = get_cell_topology_id_map();

      CellTopologyIdMap::const_iterator it = cell_topology_map.find(id);
  
      if (it == cell_topology_map.end())
        {
          //throw RuntimeError() << "Cell topology " << id << " is not defined";
          std::ostringstream msg;
          msg << "Cell topology " << id << " is not defined";
          throw std::runtime_error(msg.str());
        }

      return (*it).second;
    }


    TopologyId
    getCellTopologyId(
                      const CellTopology &  cell_topology)
    {
      if (cell_topology.getCellTopologyData() == 0)
        return INVALID;
  
      IdCellTopologyMap &cell_topology_map = get_id_cell_topology_map();

      IdCellTopologyMap::const_iterator it = cell_topology_map.find(cell_topology);
  
      if (it == cell_topology_map.end())
        {
          //throw RuntimeError() << "Cell topology " << cell_topology.getName() << " is not defined";
          std::ostringstream msg;
          msg << "Cell topology " << cell_topology.getName() << " is not defined";
          throw std::runtime_error(msg.str());

        }

      return (*it).second;
    }


    Elem::CellTopology
    edgeCellTopology(
                     const Elem::CellTopology &    cell_topology,
                     UInt                          ordinal)
    {
      if (ordinal >= cell_topology.getEdgeCount())
        return NULL;
  
      return cell_topology.getCellTopologyData(1, ordinal);
    }
  

    Elem::CellTopology
    faceCellTopology(
                     const Elem::CellTopology &    cell_topology,
                     UInt                          ordinal)
    {
      if (ordinal >= cell_topology.getFaceCount())
        return NULL;
  
      return cell_topology.getCellTopologyData(2, ordinal);
    }


    /** Query 2D edge or 3D face topologies. */
    Elem::CellTopology
    sideCellTopology(
                     const Elem::CellTopology &    cell_topology,
                     UInt                          ordinal)
    {
      return cell_topology.getDimension() == 3 ? faceCellTopology(cell_topology, ordinal) : edgeCellTopology(cell_topology, ordinal);
    }


    bool isElement(const Elem::CellTopology &cell_topology, unsigned spatial_dimension) {
      bool is_element = isShellElement(cell_topology)
        || isRodElement(cell_topology) || isParticleElement(cell_topology)
        || cell_topology.getDimension() == spatial_dimension; // isSolidElement(cell_topology, spatial_dimension);

      return is_element;
    }


    bool isSolidElement(const Elem::CellTopology &cell_topology, unsigned spatial_dimension) {
      bool is_solid = cell_topology.getDimension() == spatial_dimension
        && !isShellElement(cell_topology) && !isRodElement(cell_topology) && !isParticleElement(cell_topology);

      return is_solid;
    }


    bool isShellElement(const Elem::CellTopology &cell_topology) {
      return cell_topology.getSideCount() == 2;
    }


    bool isRodElement(const Elem::CellTopology &cell_topology) {
      bool is_rod = cell_topology.getDimension() == 2
        && cell_topology.getSideCount() == 1;

      return is_rod;
    }


    bool isParticleElement(const Elem::CellTopology &cell_topology) {
      bool is_particle = cell_topology.getDimension() == 1
        && cell_topology.getNodeCount() == 1;

      return is_particle;
    }


    const unsigned *
    getNodesOfEdge(
                   const Elem::CellTopology &    cell_topology,
                   unsigned                      edge)
    {
      if (edge >= cell_topology.getEdgeCount())
        return NULL;
  
      // Get the topology to test the bounds of subcell_dim and subcell_ord.
      cell_topology.getCellTopologyData(1, edge);
      return cell_topology.getCellTopologyData()->subcell[1][edge].node;
    }


    const unsigned *
    getNodesOfFace(
                   const Elem::CellTopology &    cell_topology,
                   unsigned                      face)
    {
      if (face >= cell_topology.getFaceCount())
        return NULL;
  
      // Get the topology to test the bounds of subcell_dim and subcell_ord.
      cell_topology.getCellTopologyData(2, face);
      return cell_topology.getCellTopologyData()->subcell[2][face].node;
    }


    const unsigned *
    getNodesOfSide(
                   const Elem::CellTopology &    cell_topology,
                   unsigned                      side)
    {
      return cell_topology.getDimension() == 3 ? getNodesOfFace(cell_topology, side) : getNodesOfEdge(cell_topology, side);
    }


    int
    getEdgeNode(
                const Elem::CellTopology &    cell_topology,
                unsigned                      edge,
                unsigned                      node_of_edge)
    {
      return cell_topology.getNodeMap(1, edge, node_of_edge);
    }


    int
    getFaceNode(
                const Elem::CellTopology &    cell_topology,
                unsigned                      face,
                unsigned                      node_of_face)
    {
      return cell_topology.getNodeMap(2, face, node_of_face);
    }


    int
    getSideNode(
                const Elem::CellTopology &    cell_topology,
                unsigned                      side,
                unsigned                      node_of_side)
    {
      return cell_topology.getDimension() == 3 ? getFaceNode(cell_topology, side, node_of_side) : getEdgeNode(cell_topology, side, node_of_side);
    }


    int
    getFaceEdge(
                const Elem::CellTopology &    cell_topology,
                unsigned                      face,
                unsigned                      edge_of_face)
    {
      return ::mapCellFaceEdge(cell_topology.getBaseCellTopologyData(), face, edge_of_face);
    }


    Elem::CellTopology
    nodeCellTopology(
                     const Elem::CellTopology &    cell_topology,
                     UInt                          ordinal) 
    {
      return (ordinal < cell_topology.getNodeCount()) ? cell_topology.getCellTopologyData(0, ordinal) : NULL;
    }


    bool
    isCellTopologySubsetOf(
                           const Elem::CellTopology &    cell_topology,
                           const Elem::CellTopology &    richer) 
    {
      bool result = false;
      UInt i, j;

      if (0 < cell_topology.getVertexCount()) {
        result =
          cell_topology.getDimension()        == richer.getDimension() &&
          cell_topology.getVertexCount()      == richer.getVertexCount() &&
          cell_topology.getEdgeCount()        == richer.getEdgeCount() &&
          cell_topology.getFaceCount()        == richer.getFaceCount() &&
          cell_topology.getNodeCount()        <= richer.getNodeCount();

        for (i = 0; result && i < cell_topology.getEdgeCount(); ++i) {
          result = isCellTopologySubsetOf(edgeCellTopology(cell_topology, i), edgeCellTopology(richer, i));

          for (j = 0; result && j < edgeCellTopology(cell_topology, i).getNodeCount(); ++j)
            result = getEdgeNode(cell_topology, i, j) == getEdgeNode(richer, i, j);
        }

        for (i = 0; i < cell_topology.getFaceCount() && result; ++i) {
          result = isCellTopologySubsetOf(faceCellTopology(cell_topology, i), faceCellTopology(richer, i));

          for (j = 0; result && j < faceCellTopology(cell_topology, i).getNodeCount(); ++j)
            result = getFaceNode(cell_topology, i, j) == getFaceNode(richer, i, j);
        }
      }

      return result;
    }


    unsigned
    getParametricDimension(
                           const Elem::CellTopology &    cell_topology) 
    {
      unsigned parametric_dimension = cell_topology.getDimension();
      if (isShellElement(cell_topology) || isRodElement(cell_topology) || isParticleElement(cell_topology))
        --parametric_dimension;

      return parametric_dimension;
    }

    int
    findReversePermutation(
                           const CellTopologyData &      top,
                           int                           permutation_ord)
    {
      const int nv = top.vertex_count ;
      const int np = top.permutation_count ;  
      const unsigned * const perm_node = top.permutation[permutation_ord].node ;
      int p = 0 ;
      for ( ; p < np ; ++p ) {
        const unsigned * const reverse_perm_node = top.permutation[p].node ;
        int j = 0 ;
        for ( ; j < nv && reverse_perm_node[j] == perm_node[nv - 1 - j] ; ++j )
          ;
        if ( nv == j )
          break ;
      }
      return p ;
    }

  } // namespace Elem
} // namespace adapt
} // namespace stk
