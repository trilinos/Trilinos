/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

/*------------------------------------------------------------------------*/
/*               shards : Shared Discretization Tools                     */
/*         Copyright (2002-2008, 2010, 2011) National Technology &        */
/*                 Engineering Solutions of Sandia, LLC.                  */
/*                                                                        */
/* Under terms of Contract DE-NA0003525, there is a non-exclusive         */
/* license for use of this work by or on behalf of the U.S. Government.   */
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
#include <adapt/sierra_element/CellTopology.hpp>

using namespace shards;

namespace percept {
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

    int
    getFaceEdge(
                const Elem::CellTopology &    cell_topology,
                unsigned                      face,
                unsigned                      edge_of_face)
    {
      return ::mapCellFaceEdge(cell_topology.getBaseCellTopologyData(), face, edge_of_face);
    }

  } // namespace Elem
} // namespace percept

